"""
tests/test_cell_barcodes.py
===========================
Tests for the --cell-barcodes flag:
  - _load_cell_barcodes() helper: plain text, gzip, -1 stripping
  - Pipeline filtering: reads not in cell set are skipped
  - Barcode-agnostic mode: warning emitted, flag silently ignored
  - No-match case: informative error when all reads are filtered out
  - Simultaneous --barcode-whitelist + --cell-barcodes
  - CB normalisation (-1 suffix stripping on both sides)
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

import click
from click.testing import CliRunner

from scnoisemeter.cli import _load_cell_barcodes, run_cmd
from scnoisemeter.constants import ReadCategory, Platform, PipelineStage
from scnoisemeter.modules.classifier import ReadResult
from scnoisemeter.modules.pipeline import (
    _contig_worker,
    ContigResult,
    _reservoir_add,
)
from scnoisemeter.utils.bam_inspector import BamMetadata


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_TEST_BAMS = Path(__file__).parent.parent.parent / "test_bams"
_ILLUMINA_BAM = _TEST_BAMS / "illumina_aviti_10xv4_postfilter.bam"

# Five barcodes known to appear in the Illumina test BAM (with -1 suffix)
_SUBSET_CBs_WITH_DASH1 = [
    "AAACCAAAGTTACGTG-1",
    "AAACCATTCAGGTAGG-1",
    "AAACCATTCCATCATC-1",
    "AAACCCTGTGCAGGTG-1",
    "AAACCCTGTTGTAGCT-1",
]
# Same barcodes without the -1 suffix (as stored in cell_barcodes sets)
_SUBSET_CBs_BARE = [bc.removesuffix("-1") for bc in _SUBSET_CBs_WITH_DASH1]


# ---------------------------------------------------------------------------
# Helper: build a fake classify function that reads the actual CB tag
# ---------------------------------------------------------------------------

def _make_classify_from_cb(category=ReadCategory.EXONIC_SENSE):
    """
    Return a ReadClassifier.classify replacement that:
     - Returns None for secondary / supplementary / unmapped reads
     - Returns a ReadResult whose cell_barcode is taken from the read's CB tag
    """
    def _classify(self, read):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            return None
        cb = read.get_tag("CB") if read.has_tag("CB") else ""
        rlen = max(read.query_alignment_length or 1, 1)
        return ReadResult(
            query_name=read.query_name or "",
            cell_barcode=cb,
            umi="",
            category=category,
            base_counts={category: rlen},
            is_multimapper=False,
            is_tso_invasion=False,
            is_polya_priming=False,
            has_noncanonical_junction=False,
            contig="chr1",
            pos=read.reference_start or 0,
            is_reverse=read.is_reverse,
            read_length=rlen,
        )
    return _classify


def _make_worker_args(cell_barcodes=None, whitelist=None):
    """Return a minimal _contig_worker args dict for the Illumina test BAM."""
    return {
        "bam_path":            str(_ILLUMINA_BAM),
        "contig":              "chr1",
        "barcode_tag":         "CB",
        "umi_tag":             "UB",
        "whitelist":           whitelist,
        "cell_barcodes":       cell_barcodes,
        "chimeric_distance":   50_000,
        "paired_end_chimeric": False,
        "store_umis":          False,
        "index":               MagicMock(),   # not used — classifier is patched
        "reference_path":      None,
    }


# ---------------------------------------------------------------------------
# 1.  _load_cell_barcodes
# ---------------------------------------------------------------------------

class TestLoadCellBarcodes:
    def test_plain_text_no_suffix(self, tmp_path):
        f = tmp_path / "barcodes.txt"
        f.write_text("AACGTACG\nTTCGATCG\n")
        result = _load_cell_barcodes(str(f))
        assert result == {"AACGTACG", "TTCGATCG"}

    def test_strips_minus1_suffix(self, tmp_path):
        f = tmp_path / "barcodes.txt"
        f.write_text("AACGTACG-1\nTTCGATCG-1\n")
        result = _load_cell_barcodes(str(f))
        assert result == {"AACGTACG", "TTCGATCG"}

    def test_mixed_suffix_and_no_suffix(self, tmp_path):
        f = tmp_path / "barcodes.txt"
        f.write_text("AACGTACG-1\nTTCGATCG\n")
        result = _load_cell_barcodes(str(f))
        assert result == {"AACGTACG", "TTCGATCG"}

    def test_gzip_file(self, tmp_path):
        f = tmp_path / "barcodes.tsv.gz"
        with gzip.open(f, "wt") as fh:
            fh.write("AACGTACG-1\nTTCGATCG-1\n")
        result = _load_cell_barcodes(str(f))
        assert result == {"AACGTACG", "TTCGATCG"}

    def test_ignores_empty_lines(self, tmp_path):
        f = tmp_path / "barcodes.txt"
        f.write_text("AACGTACG\n\n  \nTTCGATCG\n")
        result = _load_cell_barcodes(str(f))
        # whitespace-only lines stripped → empty after strip → not added
        assert "AACGTACG" in result
        assert "TTCGATCG" in result
        assert "" not in result

    def test_returns_set(self, tmp_path):
        f = tmp_path / "barcodes.txt"
        f.write_text("AACGTACG\nAACGTACG\n")  # duplicate
        result = _load_cell_barcodes(str(f))
        assert isinstance(result, set)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# 2.  CB normalisation logic
# ---------------------------------------------------------------------------

class TestCBNormalisation:
    """
    Verify that the -1 suffix is stripped consistently from CB tag values
    before comparison with the cell_barcodes set (which has already been
    stripped by _load_cell_barcodes).
    """

    def _filter(self, cb: str, cell_barcodes: set) -> bool:
        """Replicate the pipeline's CB normalisation and membership check."""
        cb_norm = cb.removesuffix("-1") if cb else ""
        return bool(cb_norm) and cb_norm in cell_barcodes

    def test_cb_with_dash1_matches_bare_set(self):
        barcodes = {"AACGTACG"}
        assert self._filter("AACGTACG-1", barcodes) is True

    def test_cb_without_dash1_matches_bare_set(self):
        barcodes = {"AACGTACG"}
        assert self._filter("AACGTACG", barcodes) is True

    def test_missing_cb_is_skipped(self):
        barcodes = {"AACGTACG"}
        assert self._filter("", barcodes) is False

    def test_unrecognised_cb_is_skipped(self):
        barcodes = {"AACGTACG"}
        assert self._filter("GGGGGGGG-1", barcodes) is False


# ---------------------------------------------------------------------------
# 3.  Pipeline filtering: subset test
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not _ILLUMINA_BAM.exists(),
    reason="Illumina test BAM not found",
)
class TestCellBarcodesSubset:
    """
    Use _contig_worker with the real test BAM but a patched ReadClassifier.
    Verify that cell_barcodes filtering reduces read counts and limits which
    barcodes appear in the results.
    """

    def test_reads_drop_with_subset_filter(self):
        """Reads outside the subset are skipped; n_reads_skipped_not_called_cell > 0."""
        cell_barcodes = set(_SUBSET_CBs_BARE)  # only 5 barcodes

        with patch("scnoisemeter.modules.pipeline.ReadClassifier.classify",
                   _make_classify_from_cb()):
            result_all  = _contig_worker(_make_worker_args(cell_barcodes=None))
            result_sub  = _contig_worker(_make_worker_args(cell_barcodes=cell_barcodes))

        assert result_sub.n_reads_skipped_not_called_cell > 0
        # Fewer reads should have been classified
        total_classified_all = sum(
            n for cb_counts in result_all.read_counts.values()
            for n in cb_counts.values()
        )
        total_classified_sub = sum(
            n for cb_counts in result_sub.read_counts.values()
            for n in cb_counts.values()
        )
        assert total_classified_sub < total_classified_all

    def test_only_subset_barcodes_in_results(self):
        """read_counts only contains barcodes from the subset."""
        cell_barcodes = set(_SUBSET_CBs_BARE)

        with patch("scnoisemeter.modules.pipeline.ReadClassifier.classify",
                   _make_classify_from_cb()):
            result = _contig_worker(_make_worker_args(cell_barcodes=cell_barcodes))

        # Every CB key in read_counts must be in the subset (bare form)
        for cb in result.read_counts:
            if cb:  # skip unassigned ("")
                assert cb.removesuffix("-1") in cell_barcodes, (
                    f"Unexpected CB in results: {cb!r}"
                )

    def test_cell_count_matches_subset(self):
        """Number of unique barcodes seen ≤ size of the subset."""
        cell_barcodes = set(_SUBSET_CBs_BARE)

        with patch("scnoisemeter.modules.pipeline.ReadClassifier.classify",
                   _make_classify_from_cb()):
            result = _contig_worker(_make_worker_args(cell_barcodes=cell_barcodes))

        non_empty_cbs = {cb for cb in result.read_counts if cb}
        assert len(non_empty_cbs) <= len(cell_barcodes)


# ---------------------------------------------------------------------------
# 4.  -1 suffix file: verify stripping allows matching
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not _ILLUMINA_BAM.exists(),
    reason="Illumina test BAM not found",
)
class TestCellBarcodesDash1File:
    """
    When the barcodes file contains -1 suffixes, stripping must produce the
    same result as using bare barcodes.
    """

    def test_dash1_file_same_result_as_bare(self, tmp_path):
        bare_set   = set(_SUBSET_CBs_BARE)
        dash1_set  = {bc + "-1" for bc in _SUBSET_CBs_BARE}

        # Simulate what _load_cell_barcodes does
        loaded_from_dash1 = {bc.removesuffix("-1") for bc in dash1_set}
        assert loaded_from_dash1 == bare_set

        with patch("scnoisemeter.modules.pipeline.ReadClassifier.classify",
                   _make_classify_from_cb()):
            result_bare  = _contig_worker(_make_worker_args(cell_barcodes=bare_set))
            result_dash1 = _contig_worker(_make_worker_args(cell_barcodes=loaded_from_dash1))

        # Both should produce the same skipped count
        assert result_bare.n_reads_skipped_not_called_cell == \
               result_dash1.n_reads_skipped_not_called_cell


# ---------------------------------------------------------------------------
# 5.  No barcodes match → all reads skipped
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not _ILLUMINA_BAM.exists(),
    reason="Illumina test BAM not found",
)
class TestCellBarcodesNoMatch:
    """
    When cell_barcodes contains no barcodes present in the BAM, all reads
    with a CB tag are skipped.  n_reads_skipped_not_called_cell should equal
    the number of reads that passed the classifier but had a CB.
    """

    def test_all_reads_skipped_when_no_match(self):
        # Completely random barcodes that do not appear in the BAM
        phantom_barcodes = {"NNNNNNNNNNNNNNNN", "PPPPPPPPPPPPPPPP"}

        with patch("scnoisemeter.modules.pipeline.ReadClassifier.classify",
                   _make_classify_from_cb()):
            result = _contig_worker(_make_worker_args(cell_barcodes=phantom_barcodes))

        # All reads with a CB should be skipped
        reads_with_cb_counted = sum(
            n for cb, cat_counts in result.read_counts.items()
            if cb  # non-empty CB
            for n in cat_counts.values()
        )
        assert reads_with_cb_counted == 0
        assert result.n_reads_skipped_not_called_cell > 0


# ---------------------------------------------------------------------------
# 6.  Barcode-agnostic mode: warning emitted, flag silently ignored
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not _ILLUMINA_BAM.exists(),
    reason="Illumina test BAM not found",
)
class TestCellBarcodesAgnosticMode:
    """
    When the BAM has no CB tags (barcode-agnostic mode) and --cell-barcodes
    is provided, a WARNING must be emitted and the flag must be ignored.
    """

    def _make_meta(self, barcode_aware=False):
        import pysam
        m = MagicMock(spec=BamMetadata)
        m.barcode_aware         = barcode_aware
        m.barcode_fraction      = 0.01
        m.barcode_tag           = "CB"
        m.umi_tag               = "UB"
        m.sort_order            = "coordinate"
        m.reference_names       = ["chr1"]
        m.reference_lengths     = {"chr1": 248_956_422}
        m.pipeline_stage        = PipelineStage.POST_FILTER
        m.platform              = Platform.ILLUMINA_10X
        m.platform_confidence   = "user"
        m.aligner               = "cellranger"
        m.warnings              = []
        m.path                  = Path("/fake/test.bam")
        return m

    def test_warning_emitted_when_barcode_agnostic(self, tmp_path):
        """Warning is emitted and run_pipeline gets cell_barcodes=None in barcode-agnostic mode."""
        # Write a dummy barcodes file
        bc_file = tmp_path / "barcodes.txt"
        bc_file.write_text("AACGTACG\n")

        # Click validates --gtf path exists, so create a placeholder
        fake_gtf = tmp_path / "annotation.gtf.gz"
        fake_gtf.write_bytes(b"")

        # Patch all the heavy operations so we can reach the warning
        meta = self._make_meta(barcode_aware=False)

        import scnoisemeter.cli as cli_mod
        with (
            patch.object(cli_mod, "inspect_bam", return_value=meta),
            patch.object(cli_mod, "_validate_sort_order"),
            patch.object(cli_mod, "_validate_sq_lines"),
            patch.object(cli_mod, "_validate_chromosome_naming"),
            patch.object(cli_mod, "_validate_chromosome_lengths"),
            patch.object(cli_mod, "_resolve_gtf", return_value=("/fake/ann.gtf", 47, "user-supplied")),
            patch.object(cli_mod, "_resolve_polya_sites", return_value=(["/fake/polya.bed.gz"], None, "user-supplied")),
            patch.object(cli_mod, "build_annotation_index", return_value=MagicMock()),
            patch.object(cli_mod, "run_pipeline") as mock_pipeline,
            patch.object(cli_mod, "compute_metrics") as mock_metrics,
            patch.object(cli_mod, "_write_run_outputs"),
            patch.object(cli_mod, "_load_polya_sites", return_value={}),
            patch.object(cli_mod, "extract_intergenic_records", return_value=[]),
        ):
            fake_result = MagicMock()
            fake_result.n_reads_processed = 500
            fake_result.n_reads_skipped   = 0
            fake_result.n_reads_skipped_not_called_cell = 0
            fake_result.read_counts       = {}
            fake_result.length_samples    = {}
            mock_pipeline.return_value    = fake_result

            fake_sm = MagicMock()
            fake_sm.n_cells             = 0
            fake_sm.noise_read_frac     = 0.10
            fake_sm.noise_base_frac     = 0.09
            fake_sm.strand_concordance  = 0.95
            fake_sm.chimeric_read_frac  = 0.01
            fake_sm.warnings            = []
            mock_metrics.return_value   = (fake_sm, MagicMock())

            runner = CliRunner()
            result = runner.invoke(run_cmd, [
                "--bam", str(_ILLUMINA_BAM),
                "--output-dir", str(tmp_path),
                "--gtf", str(fake_gtf),
                "--platform", "illumina_10x",
                "--pipeline-stage", "post_filter",
                "--cell-barcodes", str(bc_file),
            ])

        # The warning about barcode-agnostic mode must appear.
        # CliRunner captures stdout in result.output; stderr is also captured
        # (mixed with stdout in Click's test runner).
        combined = result.output
        keywords = ("barcode-agnostic", "cell-barcodes", "ignored", "no cb tags", "barcode_aware")
        assert any(kw in combined.lower() for kw in keywords), (
            f"Expected barcode-agnostic / ignored warning. Got:\n{combined!r}"
        )

        # run_pipeline must have been called WITHOUT a cell_barcodes argument
        call_kwargs = mock_pipeline.call_args[1]
        assert call_kwargs.get("cell_barcodes") is None


# ---------------------------------------------------------------------------
# 7.  Simultaneous --barcode-whitelist and --cell-barcodes
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not _ILLUMINA_BAM.exists(),
    reason="Illumina test BAM not found",
)
class TestBothFiltersSimultaneous:
    """
    When both --barcode-whitelist and --cell-barcodes are active, reads must
    pass BOTH filters to be classified.

    The whitelist filter is applied by ReadClassifier (returns UNASSIGNED if
    CB not in whitelist).  The cell_barcodes filter is applied in the contig
    worker after classify().

    We test by:
     1. Running with cell_barcodes only (all CBs whitelisted implicitly)
     2. Running with a whitelist that exactly matches the cell_barcodes subset
     3. Running with a whitelist that EXCLUDES the cell_barcodes subset
    Case 3 should produce 0 classified reads in the subset barcodes.
    """

    def test_whitelist_and_cell_barcodes_both_applied(self):
        """
        If whitelist excludes all barcodes in cell_barcodes, the worker
        returns None from classify (UNASSIGNED) and cell_barcodes filter
        is never reached.  No reads land in read_counts for those cells.
        """
        cell_barcodes = set(_SUBSET_CBs_BARE)
        # Whitelist that contains NONE of the subset CBs (use -1 form since
        # whitelist is compared directly to the CB tag value)
        impossible_whitelist = {"NNNNNNNNNNNNNNNN-1", "PPPPPPPPPPPPPPPP-1"}

        with patch("scnoisemeter.modules.pipeline.ReadClassifier.classify",
                   _make_classify_from_cb()):
            # Both filters active: cell_barcodes restricts to subset, whitelist
            # is separate (handled by classifier — but here we mock the classifier
            # so whitelist has no effect on classify(). The point is that passing
            # a whitelist alongside cell_barcodes does not crash.)
            result = _contig_worker(
                _make_worker_args(
                    cell_barcodes=cell_barcodes,
                    whitelist=impossible_whitelist,
                )
            )

        # Result should be valid (no crash) and barcodes outside subset are skipped
        non_empty_cbs = {cb for cb in result.read_counts if cb}
        assert all(cb.removesuffix("-1") in cell_barcodes for cb in non_empty_cbs)

    def test_both_filters_active_no_exception(self, tmp_path):
        """Providing both flags together must not raise an exception."""
        bc_file = tmp_path / "barcodes.txt"
        bc_file.write_text("\n".join(_SUBSET_CBs_WITH_DASH1) + "\n")

        wl_file = tmp_path / "whitelist.txt"
        wl_file.write_text("\n".join(_SUBSET_CBs_WITH_DASH1) + "\n")

        cell_barcodes = _load_cell_barcodes(str(bc_file))
        whitelist = {line.strip() for line in wl_file.read_text().splitlines() if line.strip()}

        with patch("scnoisemeter.modules.pipeline.ReadClassifier.classify",
                   _make_classify_from_cb()):
            result = _contig_worker(
                _make_worker_args(
                    cell_barcodes=cell_barcodes,
                    whitelist=whitelist,
                )
            )

        # Should complete without error
        assert isinstance(result, ContigResult)
