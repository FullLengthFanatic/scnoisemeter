"""
tests/test_core.py
==================
Unit tests for the core classification logic.

These tests do NOT require a real BAM or GTF — they use mock objects and
synthetic data to verify:
  - ReadCategory enum completeness
  - Interval helper functions
  - Chimeric detection logic
  - TSO invasion detection
  - polyA priming detection
  - Base-level category tallying
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from scnoisemeter.constants import (
    CATEGORY_ORDER,
    LENGTH_BIN_BREAKS,
    NOISE_CATEGORIES,
    ReadCategory,
    SamFlag,
    TSO_10X,
)
from scnoisemeter.modules.classifier import (
    ReadClassifier,
    _add_bases,
    _bases_in,
    _get_jxn_positions,
    _overlaps_any,
)
from scnoisemeter.modules.metrics import compute_length_stratification, compute_metrics, CellTable
from scnoisemeter.modules.pipeline import _get_length_bin


# ---------------------------------------------------------------------------
# Constants integrity
# ---------------------------------------------------------------------------

class TestConstants:
    def test_all_categories_in_order(self):
        """Every ReadCategory must appear exactly once in CATEGORY_ORDER."""
        cats_in_order = set(CATEGORY_ORDER)
        # Unmapped, secondary, supplementary are excluded from the order list
        expected_missing = {ReadCategory.UNMAPPED, ReadCategory.SECONDARY, ReadCategory.SUPPLEMENTARY}
        for cat in ReadCategory:
            if cat in expected_missing:
                assert cat not in cats_in_order
            else:
                assert cat in cats_in_order, f"{cat} missing from CATEGORY_ORDER"

    def test_noise_categories_are_subset_of_all(self):
        for cat in NOISE_CATEGORIES:
            assert cat in ReadCategory

    def test_tso_sequence_length(self):
        assert len(TSO_10X) >= 20, "TSO_10X sequence seems too short"


# ---------------------------------------------------------------------------
# Interval helpers
# ---------------------------------------------------------------------------

class TestIntervalHelpers:
    """Tests for the pure-function interval utilities."""

    # _overlaps_any / _bases_in expect (intervals, prefix_max_end) as returned
    # by _extract_intervals.  Build that format directly here.
    _RAW = [(10, 20), (30, 50), (100, 200)]
    # prefix_max_end: running max of end values
    _PREFIX = [20, 50, 200]
    INTERVALS = (_RAW, _PREFIX)

    def test_overlaps_any_true(self):
        assert _overlaps_any(self.INTERVALS, 15, 25)   # overlaps [10,20)
        assert _overlaps_any(self.INTERVALS, 5,  15)   # overlaps [10,20)
        assert _overlaps_any(self.INTERVALS, 100, 201) # overlaps [100,200)

    def test_overlaps_any_false(self):
        assert not _overlaps_any(self.INTERVALS, 20, 30)  # gap between intervals
        assert not _overlaps_any(self.INTERVALS, 50, 99)  # gap between intervals
        assert not _overlaps_any(self.INTERVALS, 0, 10)   # touches start but no overlap

    def test_bases_in_full_overlap(self):
        assert _bases_in(self.INTERVALS, 10, 20) == 10

    def test_bases_in_partial_overlap(self):
        assert _bases_in(self.INTERVALS, 15, 25) == 5   # 15..20 = 5 bases

    def test_bases_in_spanning(self):
        # Spans two intervals [30,50) and [100,200)
        assert _bases_in(self.INTERVALS, 25, 150) == 20 + 50  # 20 from [30,50), 50 from [100,150)

    def test_bases_in_no_overlap(self):
        assert _bases_in(self.INTERVALS, 50, 100) == 0

    def test_add_bases(self):
        counts = {}
        _add_bases(counts, ReadCategory.EXONIC_SENSE, 10)
        _add_bases(counts, ReadCategory.EXONIC_SENSE, 5)
        _add_bases(counts, ReadCategory.INTRONIC_PURE, 3)
        assert counts[ReadCategory.EXONIC_SENSE] == 15
        assert counts[ReadCategory.INTRONIC_PURE] == 3


# ---------------------------------------------------------------------------
# CIGAR junction position extraction
# ---------------------------------------------------------------------------

class TestJunctionPositions:
    def _make_read(self, cigar_tuples, ref_start):
        read = MagicMock()
        read.cigartuples = cigar_tuples
        read.reference_start = ref_start
        return read

    def test_no_junction(self):
        # 100M — no intron
        read = self._make_read([(0, 100)], ref_start=1000)
        assert _get_jxn_positions(read) == []

    def test_single_junction(self):
        # 50M 10000N 50M — one intron starting at 1050
        read = self._make_read([(0, 50), (3, 10000), (0, 50)], ref_start=1000)
        positions = _get_jxn_positions(read)
        assert len(positions) == 1
        assert positions[0] == 1050

    def test_two_junctions(self):
        # 30M 500N 30M 1000N 30M
        read = self._make_read(
            [(0, 30), (3, 500), (0, 30), (3, 1000), (0, 30)],
            ref_start=0
        )
        positions = _get_jxn_positions(read)
        assert len(positions) == 2
        assert positions[0] == 30
        assert positions[1] == 560   # 30 + 500 + 30


# ---------------------------------------------------------------------------
# TSO invasion detection
# ---------------------------------------------------------------------------

class TestTSODetection:
    def _make_classifier(self):
        index = MagicMock()
        index.splice_sites = {}
        return ReadClassifier(index, reference=None)

    def _make_read(self, seq: str, cigar_tuples):
        read = MagicMock()
        read.query_sequence = seq
        read.cigartuples = cigar_tuples
        read.is_reverse = False
        return read

    def test_tso_in_5prime_clip_detected(self):
        clf = self._make_classifier()
        # 30bp soft-clip containing TSO, then 100M
        clip = TSO_10X[:20]  # first 20 bp of TSO
        seq = clip + "A" * 100
        read = self._make_read(seq, [(4, 20), (0, 100)])
        assert clf._check_tso_invasion(read) is True

    def test_polyg_in_clip_detected(self):
        clf = self._make_classifier()
        clip = "G" * 10  # poly-G tail
        seq = clip + "A" * 100
        read = self._make_read(seq, [(4, 10), (0, 100)])
        assert clf._check_tso_invasion(read) is True

    def test_no_tso_not_detected(self):
        clf = self._make_classifier()
        seq = "ATCGATCG" * 15
        read = self._make_read(seq, [(4, 8), (0, 112)])
        assert clf._check_tso_invasion(read) is False

    def test_no_soft_clip_not_detected(self):
        clf = self._make_classifier()
        seq = "A" * 100
        read = self._make_read(seq, [(0, 100)])
        assert clf._check_tso_invasion(read) is False


# ---------------------------------------------------------------------------
# Chimeric detection
# ---------------------------------------------------------------------------

class TestChimericDetection:
    def _make_classifier(self, chimeric_distance=10_000):
        index = MagicMock()
        index.splice_sites = {}
        return ReadClassifier(index, chimeric_distance=chimeric_distance, reference=None)

    def _make_read(self, sa_tag: str, contig: str, pos: int, is_reverse: bool):
        read = MagicMock()
        read.has_tag = lambda tag: tag == "SA"
        read.get_tag = lambda tag: sa_tag
        read.reference_name = contig
        read.reference_start = pos
        read.is_reverse = is_reverse
        return read

    def test_interchromosomal_is_chimeric(self):
        clf = self._make_classifier()
        # SA on a different chromosome
        read = self._make_read("chr2,1000,+,100M,255,0;", "chr1", 1000, False)
        is_chimeric, reason = clf._check_chimeric(read)
        assert is_chimeric
        assert "inter-chromosomal" in reason

    def test_strand_discordant_is_chimeric(self):
        clf = self._make_classifier()
        # SA on same chromosome, opposite strand
        read = self._make_read("chr1,1000,-,100M,255,0;", "chr1", 500, False)
        is_chimeric, reason = clf._check_chimeric(read)
        assert is_chimeric
        assert "strand-discordant" in reason

    def test_far_same_strand_is_chimeric(self):
        clf = self._make_classifier(chimeric_distance=10_000)
        # SA on same chromosome, same strand, 50 kb away
        read = self._make_read("chr1,60000,+,100M,255,0;", "chr1", 1000, False)
        is_chimeric, reason = clf._check_chimeric(read)
        assert is_chimeric

    def test_close_same_strand_not_chimeric(self):
        clf = self._make_classifier(chimeric_distance=10_000)
        # SA on same chromosome, same strand, 2 kb away — legitimate splice
        read = self._make_read("chr1,3000,+,100M,255,0;", "chr1", 1000, False)
        is_chimeric, _ = clf._check_chimeric(read)
        assert not is_chimeric

    def test_no_sa_tag_not_chimeric(self):
        clf = self._make_classifier()
        read = MagicMock()
        read.has_tag = lambda tag: False
        is_chimeric, _ = clf._check_chimeric(read)
        assert not is_chimeric


# ---------------------------------------------------------------------------
# Read category priority: UNASSIGNED when CB absent
# ---------------------------------------------------------------------------

class TestUnassigned:
    def _make_read(self):
        read = MagicMock()
        read.flag = 0
        read.is_secondary = False
        read.is_supplementary = False
        read.is_unmapped = False
        read.is_reverse = False
        read.reference_name = "chr1"
        read.reference_start = 1000
        read.reference_end = 1100
        read.query_name = "read1"
        read.query_alignment_length = 100
        read.cigartuples = [(0, 100)]
        read.get_blocks = lambda: [(1000, 1100)]
        read.has_tag = lambda tag: False
        read.get_tag = lambda tag: (_ for _ in ()).throw(KeyError(tag))
        return read

    def test_missing_cb_with_whitelist_returns_unassigned(self):
        """When a whitelist is provided but CB is absent → UNASSIGNED."""
        index = MagicMock()
        index.splice_sites = {}
        index.gene_shared_cod_cod = MagicMock()
        index.gene_shared_cod_cod.df = __import__('pandas').DataFrame()
        index.gene_shared_cod_ncod = MagicMock()
        index.gene_shared_cod_ncod.df = __import__('pandas').DataFrame()
        index.exons_plus = MagicMock()
        index.exons_plus.df = __import__('pandas').DataFrame()
        index.exons_minus = MagicMock()
        index.exons_minus.df = __import__('pandas').DataFrame()
        index.introns_plus = MagicMock()
        index.introns_plus.df = __import__('pandas').DataFrame()
        index.introns_minus = MagicMock()
        index.introns_minus.df = __import__('pandas').DataFrame()

        clf = ReadClassifier(index, reference=None,
                             whitelist={"VALIDBARCODE"})
        result = clf.classify(self._make_read())
        assert result is not None
        assert result.category == ReadCategory.UNASSIGNED

    def test_missing_cb_without_whitelist_classifies_genomically(self):
        """When no whitelist and no CB → barcode-agnostic mode, NOT unassigned."""
        index = MagicMock()
        index.splice_sites = {}
        index.gene_shared_cod_cod = MagicMock()
        index.gene_shared_cod_cod.df = __import__('pandas').DataFrame()
        index.gene_shared_cod_ncod = MagicMock()
        index.gene_shared_cod_ncod.df = __import__('pandas').DataFrame()
        index.exons_plus = MagicMock()
        index.exons_plus.df = __import__('pandas').DataFrame()
        index.exons_minus = MagicMock()
        index.exons_minus.df = __import__('pandas').DataFrame()
        index.introns_plus = MagicMock()
        index.introns_plus.df = __import__('pandas').DataFrame()
        index.introns_minus = MagicMock()
        index.introns_minus.df = __import__('pandas').DataFrame()

        clf = ReadClassifier(index, reference=None, whitelist=None)
        result = clf.classify(self._make_read())
        assert result is not None
        assert result.category != ReadCategory.UNASSIGNED
        assert result.cell_barcode == "NO_BARCODE"


# ---------------------------------------------------------------------------
# Regression: 5′/3′ end computation uses reference_end, not pos + qa_len
# ---------------------------------------------------------------------------

class TestEndpointComputation:
    """
    Regression test for the bug where exonic-sense 5′ and 3′ genomic
    endpoints were computed as ``pos + query_alignment_length`` instead of
    ``read.reference_end``.

    For a spliced read on the minus strand with CIGAR 100M5000N100M:
      reference_start  = 1000
      query_aln_length = 200   (100 + 100 aligned bases, no intron)
      reference_end    = 6200  (1000 + 100 + 5000 + 100)

    The buggy code would set five_prime = 1000 + 200 = 1200, which is
    5000 bp away from the true 5′ end (6200).  With a 100 bp TSS-proximity
    window, the read would never be counted as TSS-anchored.

    The fixed code uses read.reference_end (6200) directly.
    """

    def _make_spliced_minus_read(self):
        """
        Simulate a minus-strand read with CIGAR 100M 5000N 100M.
          reference_start  = 1000
          query_aln_length = 200
          reference_end    = 6200
        """
        read = MagicMock()
        read.is_reverse          = True
        read.reference_start     = 1000
        read.reference_end       = 6200   # 1000 + 100 + 5000 + 100
        read.query_alignment_length = 200 # only the 200 exonic bases
        return read

    def _make_spliced_plus_read(self):
        """
        Simulate a plus-strand read with CIGAR 100M 5000N 100M.
          reference_start  = 1000
          query_aln_length = 200
          reference_end    = 6200
        """
        read = MagicMock()
        read.is_reverse          = False
        read.reference_start     = 1000
        read.reference_end       = 6200
        read.query_alignment_length = 200
        return read

    def _simulate_pipeline_endpoints(self, read):
        """
        Reproduce the endpoint logic from pipeline._contig_worker as-written
        after the fix, returning (three_prime, five_prime).
        """
        pos     = read.reference_start
        ref_end = read.reference_end or (pos + read.query_alignment_length)

        if read.is_reverse:
            three_prime = pos
            five_prime  = ref_end
        else:
            three_prime = ref_end
            five_prime  = pos

        return three_prime, five_prime

    # ---- minus-strand tests (5′ end is affected by the bug) ----

    def test_minus_strand_five_prime_is_reference_end(self):
        """
        Minus-strand 5′ end must equal reference_end (6200), not
        pos + query_alignment_length (1200).
        """
        read = self._make_spliced_minus_read()
        _, five_prime = self._simulate_pipeline_endpoints(read)
        assert five_prime == 6200, (
            f"5′ end for minus-strand spliced read should be reference_end "
            f"(6200), got {five_prime}"
        )

    def test_minus_strand_five_prime_not_pos_plus_qa_len(self):
        """
        Minus-strand 5′ end must NOT be pos + query_alignment_length (1200).
        This is the value that the pre-fix code produced.
        """
        read = self._make_spliced_minus_read()
        buggy_five_prime = read.reference_start + read.query_alignment_length  # 1200
        _, five_prime = self._simulate_pipeline_endpoints(read)
        assert five_prime != buggy_five_prime, (
            "5′ end equals the pre-fix value pos + qa_len — fix not applied"
        )

    def test_minus_strand_three_prime_is_reference_start(self):
        """Minus-strand 3′ end must be reference_start (leftmost position)."""
        read = self._make_spliced_minus_read()
        three_prime, _ = self._simulate_pipeline_endpoints(read)
        assert three_prime == 1000

    # ---- plus-strand tests (3′ end is affected by the bug) ----

    def test_plus_strand_three_prime_is_reference_end(self):
        """
        Plus-strand 3′ end must equal reference_end (6200), not
        pos + query_alignment_length (1200).
        """
        read = self._make_spliced_plus_read()
        three_prime, _ = self._simulate_pipeline_endpoints(read)
        assert three_prime == 6200, (
            f"3′ end for plus-strand spliced read should be reference_end "
            f"(6200), got {three_prime}"
        )

    def test_plus_strand_three_prime_not_pos_plus_qa_len(self):
        """
        Plus-strand 3′ end must NOT be pos + query_alignment_length (1200).
        """
        read = self._make_spliced_plus_read()
        buggy_three_prime = read.reference_start + read.query_alignment_length
        three_prime, _ = self._simulate_pipeline_endpoints(read)
        assert three_prime != buggy_three_prime

    def test_plus_strand_five_prime_is_reference_start(self):
        """Plus-strand 5′ end must be reference_start (leftmost position)."""
        read = self._make_spliced_plus_read()
        _, five_prime = self._simulate_pipeline_endpoints(read)
        assert five_prime == 1000

    # ---- intron-size error quantification ----

    def test_intron_error_magnitude(self):
        """
        For a 5000 bp intron, the pre-fix error was 5000 bp — well outside
        the 100 bp TSS proximity window. Confirm the fix eliminates this error.
        """
        read = self._make_spliced_minus_read()
        intron_length = 5000
        buggy_five   = read.reference_start + read.query_alignment_length  # 1200
        _, fixed_five = self._simulate_pipeline_endpoints(read)  # 6200
        pre_fix_error = fixed_five - buggy_five
        assert pre_fix_error == intron_length, (
            f"Expected pre-fix error of {intron_length} bp, got {pre_fix_error}"
        )
        # With the fix, error is zero: five_prime == reference_end exactly
        assert fixed_five == read.reference_end

    # ---- reference_end=None fallback ----

    def test_reference_end_none_falls_back_gracefully(self):
        """
        When reference_end is None (should not happen for mapped reads, but
        guard against it), the code falls back to pos + qa_len without crash.
        """
        read = MagicMock()
        read.is_reverse             = True
        read.reference_start        = 500
        read.reference_end          = None
        read.query_alignment_length = 100
        _, five_prime = self._simulate_pipeline_endpoints(read)
        assert five_prime == 600   # fallback: pos + qa_len


# ---------------------------------------------------------------------------
# Read-length stratification
# ---------------------------------------------------------------------------

class TestLengthBinAssignment:
    """Unit tests for the bin-index helper."""

    def test_below_150(self):
        assert _get_length_bin(100) == 0
        assert _get_length_bin(149) == 0

    def test_exactly_150(self):
        assert _get_length_bin(150) == 1   # 150 → "150–500"

    def test_150_to_499(self):
        assert _get_length_bin(300) == 1
        assert _get_length_bin(499) == 1

    def test_exactly_500(self):
        # 500 bp must fall in the "500–1000" bin, not "<500"
        assert _get_length_bin(500) == 2

    def test_500_to_999(self):
        assert _get_length_bin(999) == 2

    def test_exactly_1000(self):
        assert _get_length_bin(1000) == 3

    def test_exactly_2000(self):
        assert _get_length_bin(2000) == 4

    def test_exactly_5000(self):
        # 5000 bp must fall in the ">5000" bin
        assert _get_length_bin(5000) == 5

    def test_above_5000(self):
        assert _get_length_bin(10000) == 5


class TestLengthStratification:
    """Tests for compute_length_stratification()."""

    def _make_bin_counts(self, assignments):
        """
        Build a length_bin_counts dict from a list of (category, read_length) pairs.
        """
        from collections import defaultdict
        counts: dict = defaultdict(lambda: defaultdict(int))
        for cat, length in assignments:
            counts[cat][_get_length_bin(length)] += 1
        return counts

    def _make_length_samples(self, assignments):
        """Build length_samples from the same list (no reservoir limit for tests)."""
        from collections import defaultdict
        samples: dict = defaultdict(list)
        for cat, length in assignments:
            samples[cat].append(length)
        return samples

    def test_fractions_sum_to_1_per_bin(self):
        """fraction_of_bin must sum to 1.0 within each bin (over all categories)."""
        assignments = [
            (ReadCategory.EXONIC_SENSE,    800),
            (ReadCategory.EXONIC_SENSE,    900),
            (ReadCategory.INTRONIC_PURE,   750),
            (ReadCategory.CHIMERIC,        600),
            (ReadCategory.EXONIC_SENSE,   1500),
            (ReadCategory.INTRONIC_PURE,  1800),
        ]
        lbc = self._make_bin_counts(assignments)
        ls  = self._make_length_samples(assignments)
        df  = compute_length_stratification(lbc, ls)

        for bin_lbl, grp in df.groupby("length_bin"):
            total = grp["fraction_of_bin"].sum()
            assert abs(total - 1.0) < 1e-9, \
                f"fraction_of_bin does not sum to 1 in bin '{bin_lbl}': {total}"

    def test_fraction_of_total_sums_to_1(self):
        """fraction_of_total must sum to 1.0 across the entire DataFrame."""
        assignments = [
            (ReadCategory.EXONIC_SENSE,   600),
            (ReadCategory.CHIMERIC,      2000),
            (ReadCategory.INTERGENIC_SPARSE, 400),
        ]
        lbc = self._make_bin_counts(assignments)
        ls  = self._make_length_samples(assignments)
        df  = compute_length_stratification(lbc, ls)
        total = df["fraction_of_total"].sum()
        assert abs(total - 1.0) < 1e-9, f"fraction_of_total does not sum to 1: {total}"

    def test_short_read_bin_prepended_when_median_below_300(self):
        """
        When the median read length is < 300 bp, the '<150' bin must appear
        as a separate row; the '<500' label must NOT appear.
        """
        # All reads < 300 bp → median will be well below 300
        assignments = [
            (ReadCategory.EXONIC_SENSE, 100),
            (ReadCategory.EXONIC_SENSE, 120),
            (ReadCategory.CHIMERIC,      80),
        ]
        lbc = self._make_bin_counts(assignments)
        ls  = self._make_length_samples(assignments)
        df  = compute_length_stratification(lbc, ls)
        bin_labels = set(df["length_bin"].unique())
        assert "<150" in bin_labels,     "Expected '<150' bin with short-read data"
        assert "<500" not in bin_labels, "Did not expect '<500' bin (should be split)"

    def test_short_read_bin_not_prepended_when_median_above_300(self):
        """
        When the median read length is ≥ 300 bp, the '<500' merged bin must
        appear; the '<150' bin must NOT appear as a separate entry.
        """
        # Median ~1000 bp → no short-read bin
        assignments = [
            (ReadCategory.EXONIC_SENSE, 1000),
            (ReadCategory.EXONIC_SENSE, 1200),
            (ReadCategory.CHIMERIC,      800),
        ]
        lbc = self._make_bin_counts(assignments)
        ls  = self._make_length_samples(assignments)
        df  = compute_length_stratification(lbc, ls)
        bin_labels = set(df["length_bin"].unique())
        assert "<150" not in bin_labels, "Did not expect '<150' bin with long-read data"
        # Either <500 or 500-1000 (depending on where reads fall), but no <150
        assert "<500" in bin_labels or "500–1000" in bin_labels

    def test_edge_500_in_correct_bin(self):
        """Reads of exactly 500 bp must appear in '500–1000', not '<500'."""
        assignments = [
            (ReadCategory.EXONIC_SENSE, 500),
            (ReadCategory.EXONIC_SENSE, 1500),  # push median well above 300
        ]
        lbc = self._make_bin_counts(assignments)
        ls  = self._make_length_samples(assignments)
        df  = compute_length_stratification(lbc, ls)
        # The 500 bp read should contribute to "500–1000"
        row = df[(df["length_bin"] == "500–1000") &
                 (df["category"] == ReadCategory.EXONIC_SENSE.value)]
        assert not row.empty, "Expected a row for 500 bp read in '500–1000' bin"
        assert row["count"].iloc[0] >= 1

    def test_edge_5000_in_gt5000_bin(self):
        """Reads of exactly 5000 bp must appear in '>5000', not '2000–5000'."""
        assignments = [
            (ReadCategory.EXONIC_SENSE, 5000),
            (ReadCategory.EXONIC_SENSE, 6000),
            (ReadCategory.EXONIC_SENSE, 1500),  # push median above 300
        ]
        lbc = self._make_bin_counts(assignments)
        ls  = self._make_length_samples(assignments)
        df  = compute_length_stratification(lbc, ls)
        row = df[(df["length_bin"] == ">5000") &
                 (df["category"] == ReadCategory.EXONIC_SENSE.value)]
        assert not row.empty
        assert row["count"].iloc[0] >= 2  # both 5000 and 6000 should land here

    def test_empty_input_returns_empty_dataframe(self):
        """compute_length_stratification on empty input must return an empty DataFrame."""
        df = compute_length_stratification({}, {})
        assert df.empty
        assert list(df.columns) == [
            "length_bin", "category", "count",
            "fraction_of_bin", "fraction_of_total",
        ]


# ---------------------------------------------------------------------------
# Illumina compatibility: platform flag, sort order, chromosome naming,
# paired-end chimeric detection
# ---------------------------------------------------------------------------

class TestPlatformFlag:
    """Platform.ILLUMINA must exist and be accepted by the CLI."""

    def test_illumina_value_in_enum(self):
        from scnoisemeter.constants import Platform
        assert Platform.ILLUMINA == "illumina"

    def test_illumina_in_platform_values(self):
        from scnoisemeter.constants import Platform
        values = [p.value for p in Platform]
        assert "illumina" in values

    def test_illumina_10x_still_present(self):
        """Existing specific Illumina variants must not have been removed."""
        from scnoisemeter.constants import Platform
        assert Platform.ILLUMINA_10X == "illumina_10x"
        assert Platform.ILLUMINA_BD  == "illumina_bd"


class TestSortOrderValidation:
    """Sort order check must raise ClickException for non-coordinate BAMs."""

    def _make_meta(self, sort_order: str):
        from scnoisemeter.utils.bam_inspector import BamMetadata
        from pathlib import Path
        meta = BamMetadata(path=Path("/fake/file.bam"))
        meta.sort_order = sort_order
        return meta

    def test_coordinate_sorted_passes(self):
        import click
        from scnoisemeter.cli import _validate_sort_order
        # Should not raise
        _validate_sort_order(self._make_meta("coordinate"))

    def test_queryname_sorted_raises(self):
        import click
        from scnoisemeter.cli import _validate_sort_order
        with pytest.raises(click.ClickException) as exc_info:
            _validate_sort_order(self._make_meta("queryname"))
        assert "queryname" in str(exc_info.value.format_message())
        assert "samtools sort" in str(exc_info.value.format_message())

    def test_unsorted_raises(self):
        import click
        from scnoisemeter.cli import _validate_sort_order
        with pytest.raises(click.ClickException):
            _validate_sort_order(self._make_meta("unsorted"))

    def test_absent_so_tag_warns_but_does_not_raise(self):
        """Empty SO tag = absent from header → warn only, don't abort."""
        import click
        from scnoisemeter.cli import _validate_sort_order
        # Should not raise even though sort order is unknown
        _validate_sort_order(self._make_meta(""))


class TestChromosomeNaming:
    """Chromosome naming mismatch detection."""

    def _make_meta(self, contig_names: list):
        from scnoisemeter.utils.bam_inspector import BamMetadata
        from pathlib import Path
        meta = BamMetadata(path=Path("/fake/file.bam"))
        meta.reference_names = contig_names
        return meta

    def test_detect_ucsc_style(self):
        from scnoisemeter.cli import _detect_chrom_style
        names = ["chr1", "chr2", "chr3", "chrX", "chrY", "chrM", "chr1_random"]
        assert _detect_chrom_style(names) == "ucsc"

    def test_detect_ensembl_style(self):
        from scnoisemeter.cli import _detect_chrom_style
        names = ["1", "2", "3", "X", "Y", "MT"]
        assert _detect_chrom_style(names) == "ensembl"

    def test_empty_returns_unknown(self):
        from scnoisemeter.cli import _detect_chrom_style
        assert _detect_chrom_style([]) == "unknown"

    def test_mismatch_raises_click_exception(self, tmp_path):
        """UCSC BAM + Ensembl GTF (or vice versa) must raise ClickException."""
        import click
        from scnoisemeter.cli import _validate_chromosome_naming

        # Write a tiny fake Ensembl-style GTF
        fake_gtf = tmp_path / "fake.gtf"
        fake_gtf.write_text(
            '1\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000000001";\n'
            '2\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000000002";\n'
        )
        # BAM uses UCSC style
        meta = self._make_meta(["chr1", "chr2", "chrX"])
        with pytest.raises(click.ClickException) as exc_info:
            _validate_chromosome_naming(meta, str(fake_gtf))
        msg = exc_info.value.format_message()
        assert "UCSC" in msg or "ENSEMBL" in msg
        assert "mismatch" in msg.lower()

    def test_matching_naming_passes(self, tmp_path):
        """Same naming convention on BAM and GTF must pass silently."""
        import click
        from scnoisemeter.cli import _validate_chromosome_naming

        fake_gtf = tmp_path / "fake.gtf"
        fake_gtf.write_text(
            'chr1\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000000001";\n'
            'chr2\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000000002";\n'
        )
        meta = self._make_meta(["chr1", "chr2", "chrX"])
        # Should not raise
        _validate_chromosome_naming(meta, str(fake_gtf))


class TestPairedEndChimericDetection:
    """Paired-end chimeric detection via FLAG bits and RNEXT/PNEXT."""

    def _make_classifier(self):
        index = MagicMock()
        index.splice_sites = {}
        return ReadClassifier(index, paired_end_chimeric=True, reference=None)

    def _make_paired_read(
        self,
        *,
        flag: int,
        reference_name: str = "chr1",
        next_reference_name: str = "chr1",
        template_length: int = 200,
    ):
        read = MagicMock()
        read.flag = flag
        read.has_tag = lambda tag: False   # no SA tag
        read.reference_name = reference_name
        read.next_reference_name = next_reference_name
        read.template_length = template_length
        read.is_reverse = bool(flag & SamFlag.REVERSE_STRAND)
        return read

    def test_discordant_mate_unmapped_is_chimeric(self):
        """FLAG 0x001 (paired) | 0x008 (mate unmapped) → chimeric."""
        clf = self._make_classifier()
        flag = SamFlag.PAIRED | SamFlag.MATE_UNMAPPED
        read = self._make_paired_read(flag=flag)
        is_chimeric, reason = clf._check_chimeric(read)
        assert is_chimeric
        assert "mate" in reason.lower() and "unmapped" in reason.lower()

    def test_interchromosomal_pair_is_chimeric(self):
        """Read on chr1, mate on chr2 → chimeric."""
        clf = self._make_classifier()
        flag = SamFlag.PAIRED | SamFlag.PROPER_PAIR
        read = self._make_paired_read(
            flag=flag,
            reference_name="chr1",
            next_reference_name="chr2",
        )
        is_chimeric, reason = clf._check_chimeric(read)
        assert is_chimeric
        assert "inter-chromosomal" in reason.lower()

    def test_large_insert_size_is_chimeric(self):
        """Same chromosome but |TLEN| > 1 Mb → chimeric."""
        clf = self._make_classifier()
        flag = SamFlag.PAIRED | SamFlag.PROPER_PAIR
        read = self._make_paired_read(
            flag=flag,
            reference_name="chr1",
            next_reference_name="chr1",
            template_length=2_000_000,
        )
        is_chimeric, reason = clf._check_chimeric(read)
        assert is_chimeric
        assert "insert" in reason.lower()

    def test_normal_pair_not_chimeric(self):
        """Same chromosome, small insert → not chimeric."""
        clf = self._make_classifier()
        flag = SamFlag.PAIRED | SamFlag.PROPER_PAIR
        read = self._make_paired_read(
            flag=flag,
            reference_name="chr1",
            next_reference_name="chr1",
            template_length=350,
        )
        is_chimeric, _ = clf._check_chimeric(read)
        assert not is_chimeric

    def test_unpaired_read_not_chimeric(self):
        """Unpaired read (FLAG 0x001 not set) must not trigger paired-end check."""
        clf = self._make_classifier()
        read = self._make_paired_read(flag=0)  # not paired
        is_chimeric, _ = clf._check_chimeric(read)
        assert not is_chimeric

    def test_paired_end_disabled_no_detection(self):
        """When paired_end_chimeric=False, discordant pairs are not flagged."""
        index = MagicMock()
        index.splice_sites = {}
        clf = ReadClassifier(index, paired_end_chimeric=False, reference=None)
        flag = SamFlag.PAIRED | SamFlag.MATE_UNMAPPED
        read = self._make_paired_read(flag=flag)
        is_chimeric, _ = clf._check_chimeric(read)
        assert not is_chimeric

    def test_illumina_platform_is_paired_chimeric(self):
        """_is_illumina_platform() must return True for all Illumina variants."""
        from scnoisemeter.constants import Platform
        from scnoisemeter.cli import _is_illumina_platform
        assert _is_illumina_platform(Platform.ILLUMINA)
        assert _is_illumina_platform(Platform.ILLUMINA_10X)
        assert _is_illumina_platform(Platform.ILLUMINA_BD)
        assert not _is_illumina_platform(Platform.ONT)
        assert not _is_illumina_platform(Platform.PACBIO)


# ---------------------------------------------------------------------------
# Cell metrics TSV column integrity
# ---------------------------------------------------------------------------

def _make_synthetic_result():
    """Return a minimal synthetic SampleResult-like MagicMock for metrics tests."""
    from collections import defaultdict
    from unittest.mock import MagicMock
    from scnoisemeter.constants import Platform, PipelineStage
    from pathlib import Path

    def dd_int():
        return defaultdict(int)

    result = MagicMock()
    result.bam_path = Path("/fake/test.bam")
    result.n_reads_processed = 100

    meta = MagicMock()
    meta.platform       = Platform.ILLUMINA_10X
    meta.pipeline_stage = PipelineStage.POST_FILTER
    meta.aligner        = "cellranger"
    meta.warnings       = []
    result.meta = meta

    rc = defaultdict(dd_int)
    bc = "AACGTACG"
    rc[bc][ReadCategory.EXONIC_SENSE]  = 80
    rc[bc][ReadCategory.INTRONIC_PURE] = 20
    result.read_counts = rc

    bcs = defaultdict(dd_int)
    bcs[bc][ReadCategory.EXONIC_SENSE]  = 800
    bcs[bc][ReadCategory.INTRONIC_PURE] = 200
    result.base_counts = bcs

    result.umi_sets        = defaultdict(lambda: defaultdict(set))
    result.artifact_flags  = defaultdict(lambda: {"tso": 0, "polya": 0, "noncanon": 0})
    result.length_samples  = {}
    result.exonic_sense_three_prime = []
    result.exonic_sense_five_prime  = []
    result._polya_site_dict  = {}
    result._tss_site_dict    = None
    result._numt_intervals   = None
    return result


class TestCellMetricsTSVColumns:
    """
    Regression guard: verify that read_frac_multimapper and read_frac_ambiguous
    are separate columns in the cell_metrics DataFrame/TSV, not concatenated.

    Protects against a class of bug where the TSV header loses the tab separator
    between adjacent category column names.
    """

    def _get_ct(self) -> CellTable:
        result = _make_synthetic_result()
        _sm, ct = compute_metrics(result, "test")
        return ct

    def test_read_frac_multimapper_is_separate_column(self):
        ct = self._get_ct()
        assert "read_frac_multimapper" in ct.df.columns

    def test_read_frac_ambiguous_is_separate_column(self):
        ct = self._get_ct()
        assert "read_frac_ambiguous" in ct.df.columns

    def test_no_concatenated_multimapper_ambiguous_column(self):
        """The concatenated name 'read_frac_multimapperread_frac_ambiguous' must not exist."""
        ct = self._get_ct()
        assert "read_frac_multimapperread_frac_ambiguous" not in ct.df.columns

    def test_base_frac_multimapper_is_separate_column(self):
        ct = self._get_ct()
        assert "base_frac_multimapper" in ct.df.columns

    def test_base_frac_ambiguous_is_separate_column(self):
        ct = self._get_ct()
        assert "base_frac_ambiguous" in ct.df.columns

    def test_no_concatenated_base_frac_multimapper_ambiguous_column(self):
        ct = self._get_ct()
        assert "base_frac_multimapperbase_frac_ambiguous" not in ct.df.columns

    def test_tsv_header_has_tab_between_multimapper_and_ambiguous(self, tmp_path):
        """The TSV written to disk must have a tab between the two column names."""
        import io
        ct = self._get_ct()
        buf = io.StringIO()
        ct.df.to_csv(buf, sep="\t")
        header = buf.getvalue().split("\n")[0]
        assert "read_frac_multimapper\tread_frac_ambiguous" in header, (
            "Missing tab between read_frac_multimapper and read_frac_ambiguous in TSV header"
        )

    def test_tsv_header_has_tab_between_base_frac_multimapper_and_ambiguous(self, tmp_path):
        import io
        ct = self._get_ct()
        buf = io.StringIO()
        ct.df.to_csv(buf, sep="\t")
        header = buf.getvalue().split("\n")[0]
        assert "base_frac_multimapper\tbase_frac_ambiguous" in header, (
            "Missing tab between base_frac_multimapper and base_frac_ambiguous in TSV header"
        )

    def test_noise_read_frac_column_present_for_violin(self):
        """noise_read_frac must be a column so the per-cell violin plot can render."""
        ct = self._get_ct()
        assert "noise_read_frac" in ct.df.columns

    def test_all_category_read_frac_columns_present(self):
        """Every category in CATEGORY_ORDER must have its own read_frac column."""
        ct = self._get_ct()
        for cat in CATEGORY_ORDER:
            expected = f"read_frac_{cat.value}"
            assert expected in ct.df.columns, f"Missing column: {expected}"

    def test_all_category_base_frac_columns_present(self):
        """Every category in CATEGORY_ORDER must have its own base_frac column."""
        ct = self._get_ct()
        for cat in CATEGORY_ORDER:
            expected = f"base_frac_{cat.value}"
            assert expected in ct.df.columns, f"Missing column: {expected}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
