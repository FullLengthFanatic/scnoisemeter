"""
tests/test_robustness_transfers.py
==================================
Screening items 8, 9, 10 and 11: checks that refuse to report a configuration
error as a result, plus two new pieces of evidence.

Each of the two guards is asserted in both directions -- it fires on the
failure it exists for, and stays silent on the legitimate case that looks
superficially similar.  A guard that only ever passes is worse than none,
because it reads as reassurance.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import click
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from scnoisemeter.cli import (  # noqa: E402
    _check_barcode_gate,
    _check_site_coordinate_range,
)
from scnoisemeter.constants import RT_SWITCH_REPEAT_MIN  # noqa: E402
from scnoisemeter.modules.classifier import ReadClassifier  # noqa: E402
from scnoisemeter.modules.intergenic_profiler import (  # noqa: E402
    IntergenicLocus,
    _assign_local_read_share,
)


class _Meta:
    def __init__(self, *, barcode_aware=True, barcode_fraction=0.97, tag="CB",
                 reference_lengths=None):
        self.barcode_aware = barcode_aware
        self.barcode_fraction = barcode_fraction
        self.barcode_tag = tag
        self.reference_lengths = reference_lengths or {}
        self.warnings: list = []


class _SM:
    def __init__(self, classified, unassigned, n_cells=0):
        self.n_reads_classified = classified
        self.n_reads_unassigned = unassigned
        self.n_cells = n_cells


# ---------------------------------------------------------------------------
# Item 8 -- barcode gate
# ---------------------------------------------------------------------------

class TestBarcodeGate:
    def test_stops_when_the_gate_rejects_almost_everything(self):
        with pytest.raises(click.ClickException, match="no usable cell barcode"):
            _check_barcode_gate(_SM(1_000_000, 940_000, 12), _Meta(), strict=True)

    def test_no_strict_downgrades_to_a_recorded_warning(self):
        meta = _Meta()
        _check_barcode_gate(_SM(1_000_000, 940_000, 12), meta, strict=False)
        assert len(meta.warnings) == 1
        # the observed numbers must be in the message, not just a threshold
        assert "940,000" in meta.warnings[0] and "12 cells" in meta.warnings[0]

    def test_healthy_sample_is_silent(self):
        meta = _Meta()
        _check_barcode_gate(_SM(1_000_000, 30_000, 5_000), meta, strict=True)
        assert meta.warnings == []

    def test_barcode_free_bam_is_not_a_fault(self):
        """A BAM with no CB tags is legitimately 100% unassigned."""
        meta = _Meta(barcode_aware=False)
        _check_barcode_gate(_SM(1_000_000, 1_000_000), meta, strict=True)
        assert meta.warnings == []

    def test_non_numeric_counters_skip_rather_than_raise(self):
        """This guard runs before every report is written.

        A mocked or absent counter must not be the thing that fails a run, so
        a value that cannot be compared to a threshold is skipped.
        """
        _check_barcode_gate(MagicMock(), MagicMock(), strict=True)
        _check_barcode_gate(_SM(None, None), _Meta(), strict=True)
        _check_barcode_gate(_SM(0, 0), _Meta(), strict=True)

    def test_boolean_counter_is_not_treated_as_one(self):
        _check_barcode_gate(_SM(True, True), _Meta(), strict=True)


# ---------------------------------------------------------------------------
# Item 9 -- assembly mismatch that contig names cannot reveal
# ---------------------------------------------------------------------------

class TestSiteCoordinateRange:
    REFS = {"chr1": 1_000_000}

    def test_in_range_atlas_passes_through_unchanged(self):
        sites = {("chr1", "+"): [100, 500_000, 999_999]}
        meta = _Meta(reference_lengths=self.REFS)
        assert _check_site_coordinate_range(
            sites, self.REFS, label="polyA site atlas", meta=meta
        ) is sites
        assert meta.warnings == []

    def test_wrong_assembly_is_discarded_not_reported_as_zero(self):
        """Every contig name matches; every coordinate is off the end.

        _check_site_contig_overlap passes this, which is the point: name
        agreement is not assembly agreement.
        """
        sites = {("chr1", "+"): [5_000_000, 6_000_000, 7_000_000]}
        meta = _Meta(reference_lengths=self.REFS)
        out = _check_site_coordinate_range(
            sites, self.REFS, label="polyA site atlas", meta=meta
        )
        assert out is None
        assert meta.warnings and "assembly mismatch" in meta.warnings[0]

    def test_a_few_out_of_range_positions_only_warn(self):
        sites = {("chr1", "+"): [i * 10 for i in range(1, 96)] + [2_000_000] * 5}
        meta = _Meta(reference_lengths=self.REFS)
        out = _check_site_coordinate_range(
            sites, self.REFS, label="polyA site atlas", meta=meta
        )
        assert out is sites
        assert meta.warnings and "under-count" in meta.warnings[0]

    def test_contigs_absent_from_the_bam_are_not_counted_against_the_atlas(self):
        sites = {("chrUn_xyz", "+"): [9_999_999]}
        meta = _Meta(reference_lengths=self.REFS)
        assert _check_site_coordinate_range(
            sites, self.REFS, label="x", meta=meta
        ) is sites

    def test_absent_inputs_pass_through(self):
        assert _check_site_coordinate_range(None, self.REFS, label="x") is None
        assert _check_site_coordinate_range({}, self.REFS, label="x") == {}
        sites = {("chr1", "+"): [1]}
        assert _check_site_coordinate_range(sites, None, label="x") is sites


# ---------------------------------------------------------------------------
# Item 10 -- RT template-switching direct repeat
# ---------------------------------------------------------------------------

class _Ref:
    def __init__(self, seq):
        self.seq = seq

    def fetch(self, _contig, start, end):
        return self.seq[start:end]


def _genome(donor_up: str, acceptor_up: str, donor=20, acceptor=40) -> str:
    seq = list("N" * 60)
    seq[donor - len(donor_up):donor] = list(donor_up)
    seq[acceptor - len(acceptor_up):acceptor] = list(acceptor_up)
    return "".join(seq)


class TestRTSwitchRepeat:
    @staticmethod
    def _classifier(seq):
        c = ReadClassifier.__new__(ReadClassifier)
        c.reference = _Ref(seq)
        return c

    def test_shared_repeat_is_measured(self):
        c = self._classifier(_genome("ACGTACGT", "ACGTACGT"))
        assert c._direct_repeat_len("chr1", 20, 40) == 8

    def test_partial_repeat_is_the_common_suffix(self):
        c = self._classifier(_genome("AAAACGTA", "TTTACGTA"))
        assert c._direct_repeat_len("chr1", 20, 40) == 5

    def test_threshold_separates_signature_from_coincidence(self):
        long_repeat = self._classifier(_genome("ACGTACGT", "ACGTACGT"))
        short = self._classifier(_genome("AAAAAAAG", "TTTTTTTG"))
        assert long_repeat._direct_repeat_len("chr1", 20, 40) >= RT_SWITCH_REPEAT_MIN
        assert short._direct_repeat_len("chr1", 20, 40) < RT_SWITCH_REPEAT_MIN

    def test_n_bases_do_not_count_as_a_repeat(self):
        """An unplaced or masked region matches itself trivially."""
        c = self._classifier("N" * 60)
        assert c._direct_repeat_len("chr1", 20, 40) == 0

    def test_no_reference_reports_nothing(self):
        c = ReadClassifier.__new__(ReadClassifier)
        c.reference = None
        assert c._has_rt_switch_repeat(MagicMock()) is False

    def test_boundaries_too_close_to_the_contig_start(self):
        c = self._classifier("ACGT" * 15)
        assert c._direct_repeat_len("chr1", 2, 40) == 0


# ---------------------------------------------------------------------------
# Item 11 -- local read share
# ---------------------------------------------------------------------------

def _locus(contig, start, n_reads):
    return IntergenicLocus(
        contig=contig, start=start, end=start + 500, strand="+",
        n_reads=n_reads, n_barcodes=1, has_splice_evidence=False,
        is_monoexonic=True, polya_run_downstream=False, near_polya_site=False,
        poisson_pvalue=1.0, poisson_pvalue_adj=1.0,
    )


class TestLocalReadShare:
    def test_isolated_window_holds_all_local_mass(self):
        loci = [_locus("chr1", 1_000_000, 100)]
        _assign_local_read_share(loci)
        assert loci[0].local_read_share == pytest.approx(1.0)

    def test_diffuse_neighbourhood_splits_the_share(self):
        loci = [_locus("chr1", 2_000_000 + i * 500, 100) for i in range(20)]
        _assign_local_read_share(loci)
        middle = loci[10]
        assert middle.local_read_share == pytest.approx(1 / 20)

    def test_contigs_are_independent(self):
        loci = [_locus("chr1", 500_000, 100), _locus("chr2", 500_000, 100)]
        _assign_local_read_share(loci)
        assert all(x.local_read_share == pytest.approx(1.0) for x in loci)

    def test_dominant_window_beside_a_trace_one(self):
        loci = [_locus("chr1", 1_000_000, 990), _locus("chr1", 1_001_000, 10)]
        _assign_local_read_share(loci)
        assert loci[0].local_read_share == pytest.approx(0.99)
        assert loci[1].local_read_share == pytest.approx(0.01)

    def test_empty_input(self):
        _assign_local_read_share([])
