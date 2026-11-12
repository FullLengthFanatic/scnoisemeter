"""
tests/test_intergenic_columns.py
================================
Screening items 5 and 6, shipped as reported columns with promotion unchanged.

Both change the same quantities as each other and as the outstanding
significance-model decision, so neither gates anything yet. These tests pin
that: the columns are populated, and the categories do not move because of
them.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from scnoisemeter.modules.intergenic_profiler import (  # noqa: E402
    IntergenicLocus,
    IntergenicReadRecord,
    _assign_strand_asymmetry,
    _score_locus,
)


def _locus(n_reads, strand_fraction):
    return IntergenicLocus(
        contig="chr1", start=0, end=500, strand="+", n_reads=n_reads,
        n_barcodes=1, has_splice_evidence=False, is_monoexonic=True,
        polya_run_downstream=False, near_polya_site=False,
        poisson_pvalue=1.0, poisson_pvalue_adj=1.0,
        strand_fraction=strand_fraction,
    )


def _record(barcode="CB1", umi="AAAA", pos=1000):
    return IntergenicReadRecord(
        contig="chr1", start=pos, end=pos + 100, strand="+",
        cell_barcode=barcode, has_junction=False, three_prime=pos + 100,
        umi=umi,
    )


class TestMoleculeCount:
    """Item 5: molecules counted and reported, never gated on."""

    def _score(self, records):
        return _score_locus(
            0, records, background_rate=1e-9, n_tests=1, min_barcodes=1,
            reference=None, polya_sites=None, repeat_intervals=None, alpha=0.01,
        )

    def test_pcr_duplicates_collapse_to_one_molecule(self):
        """30 reads of one molecule is a PCR artifact, not 30 units of support."""
        locus = self._score([_record(umi="AAAA") for _ in range(30)])
        assert locus.n_reads == 30
        assert locus.n_molecules == 1

    def test_distinct_umis_count_separately(self):
        recs = [_record(umi=u) for u in ("AAAA", "CCCC", "GGGG")]
        assert self._score(recs).n_molecules == 3

    def test_same_umi_in_different_cells_are_different_molecules(self):
        recs = [_record(barcode="CB1", umi="AAAA"), _record(barcode="CB2", umi="AAAA")]
        assert self._score(recs).n_molecules == 2

    def test_absent_umis_report_minus_one_not_zero(self):
        """No UMIs available is not the same claim as no molecules."""
        locus = self._score([_record(umi="") for _ in range(5)])
        assert locus.n_molecules == -1

    def test_molecule_count_does_not_change_the_category(self):
        """Promotion is unchanged in this release: the column is evidence only."""
        many = self._score([_record(umi=f"UMI{i}") for i in range(30)])
        one = self._score([_record(umi="AAAA") for _ in range(30)])
        assert many.n_molecules == 30 and one.n_molecules == 1
        assert many.category == one.category


class TestStrandAsymmetry:
    """Item 6: the test that replaced the decoy null, and why."""

    def test_same_fraction_at_greater_depth_is_stronger_evidence(self):
        """The existing strand_fraction threshold is depth-blind; this is not."""
        shallow, deep = _locus(4, 0.75), _locus(100, 0.75)
        _assign_strand_asymmetry([shallow, deep], stranded=True)
        assert shallow.strand_asymmetry_pvalue > 0.3
        assert deep.strand_asymmetry_pvalue < 1e-5
        assert deep.strand_asymmetry_pvalue < shallow.strand_asymmetry_pvalue

    @pytest.mark.parametrize("n,frac,expected", [
        (4, 0.75, 10 / 16),   # 3-1 or 1-3, or more lopsided
        (2, 1.0, 0.5),        # 2-0 or 0-2
        (1, 1.0, 1.0),        # one read is always one-sided
        (4, 0.5, 1.0),        # a tie
    ])
    def test_pvalue_is_two_sided(self, n, frac, expected):
        """The dominant strand is picked after counting, so both tails count."""
        locus = _locus(n, frac)
        _assign_strand_asymmetry([locus], stranded=True)
        assert locus.strand_asymmetry_pvalue == pytest.approx(expected)

    def test_balanced_window_is_not_significant_however_deep(self):
        big_balanced = _locus(200, 0.52)
        _assign_strand_asymmetry([big_balanced], stranded=True)
        assert big_balanced.strand_asymmetry_pvalue > 0.05

    def test_one_sided_window_is_significant(self):
        clean = _locus(100, 0.99)
        _assign_strand_asymmetry([clean], stranded=True)
        assert clean.strand_asymmetry_pvalue < 1e-20

    def test_unstranded_library_gets_no_number_rather_than_a_weak_one(self):
        """Reads divide 50/50 by design, so the test would measure the protocol."""
        locus = _locus(100, 0.99)
        _assign_strand_asymmetry([locus], stranded=False)
        assert math.isnan(locus.strand_asymmetry_pvalue)

    def test_default_is_absent_not_zero(self):
        """An uncomputed p-value must not read as an overwhelming one."""
        assert math.isnan(_locus(10, 0.9).strand_asymmetry_pvalue)

    @pytest.mark.parametrize("n,frac", [(0, 0.0), (1, 1.0)])
    def test_degenerate_windows_do_not_raise(self, n, frac):
        locus = _locus(n, frac)
        _assign_strand_asymmetry([locus], stranded=True)
        assert locus.strand_asymmetry_pvalue == locus.strand_asymmetry_pvalue or n == 0


class TestRejectedDecoyNull:
    """The wrong-strand decoy was implemented, measured and rejected.

    Within any window the minority-strand count cannot exceed the dominant one,
    so a decoy distribution built from minority counts sits below the target
    distribution even under pure noise. This test states the property that
    makes such an estimator anti-conservative, so nobody reintroduces it.
    """

    def test_minority_never_exceeds_dominant_so_a_decoy_is_biased(self):
        import random
        random.seed(0)
        for _ in range(2000):
            n = random.randint(1, 200)
            plus = sum(1 for _ in range(n) if random.random() < 0.5)
            dominant, minority = max(plus, n - plus), min(plus, n - plus)
            assert minority <= dominant
