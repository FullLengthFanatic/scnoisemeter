"""
tests/test_priming_evidence.py
==============================
The 0.9 internal-priming work, transferred from flightcollapse.

Each class pins one behaviour that a plausible future edit would break, and the
option defaults are asserted explicitly: every setting here changes a reported
number, so the release must be installable without moving one.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from scnoisemeter.cli import (  # noqa: E402
    _load_polya_sites,
    _polya_atlas_density,
    _warn_if_atlas_too_dense,
)
from scnoisemeter.constants import (  # noqa: E402
    POLYA_MIN_SUPPORT_DEFAULT,
    POLYA_TAIL_MIN_LENGTH,
    PRIMING_CRITERION_DEFAULT,
)
from scnoisemeter.modules.classifier import ReadClassifier  # noqa: E402


class _Read:
    """Minimal stand-in: the tail measurement reads only these three fields."""

    def __init__(self, cigartuples, query_sequence, is_reverse=False):
        self.cigartuples = cigartuples
        self.query_sequence = query_sequence
        self.is_reverse = is_reverse


@pytest.fixture
def classifier():
    return ReadClassifier.__new__(ReadClassifier)


# ---------------------------------------------------------------------------
# Untemplated 3' polyA tail from the read's own soft clip
# ---------------------------------------------------------------------------

class TestPolyATail:
    def test_soft_clipped_a_run_is_a_tail(self, classifier):
        read = _Read([(0, 50), (4, 12)], "G" * 50 + "A" * 12)
        clip, run, frac = classifier._polya_tail_stats(read)
        assert (clip, run) == (12, 12)
        assert frac == pytest.approx(1.0)
        assert classifier._has_polya_tail(read) is True

    def test_short_run_is_not_a_tail(self, classifier):
        read = _Read([(0, 50), (4, 4)], "G" * 50 + "A" * 4)
        assert classifier._polya_tail_stats(read)[1] < POLYA_TAIL_MIN_LENGTH
        assert classifier._has_polya_tail(read) is False

    def test_hard_clip_reports_no_tail_rather_than_wrong_bases(self, classifier):
        """Hard-clipped bases are absent from query_sequence.

        Indexing the aligned sequence as though it were the clip would read
        real alignment bases as a tail, so a hard-clipped BAM must report 0.
        """
        read = _Read([(0, 50), (5, 12)], "G" * 50)
        assert classifier._polya_tail_stats(read) == (0, 0, 0.0)
        assert classifier._has_polya_tail(read) is False

    def test_reverse_strand_tail_is_the_leading_clip(self, classifier):
        """On a minus-strand alignment the tail precedes the block as poly-T."""
        read = _Read([(4, 12), (0, 50)], "T" * 12 + "G" * 50, is_reverse=True)
        assert classifier._polya_tail_stats(read)[1] == 12
        assert classifier._has_polya_tail(read) is True

    def test_reverse_strand_ignores_the_trailing_clip(self, classifier):
        read = _Read([(0, 50), (4, 12)], "G" * 50 + "A" * 12, is_reverse=True)
        assert classifier._has_polya_tail(read) is False

    def test_adapter_after_the_tail_still_counts(self, classifier):
        """The run decides, not the clip's A-fraction.

        A genuine 12 bp tail inside a 60 bp clip that also holds adapter or TSO
        sequence has an A-fraction near 0.2; requiring both would silently
        demand a ~45 bp tail.
        """
        read = _Read([(0, 50), (4, 60)], "G" * 50 + "A" * 12 + "C" * 48)
        clip, run, frac = classifier._polya_tail_stats(read)
        assert (clip, run) == (60, 12)
        assert frac < 0.25
        assert classifier._has_polya_tail(read) is True


# ---------------------------------------------------------------------------
# Atlas confidence filtering and density
# ---------------------------------------------------------------------------

def _write_atlas(tmp_path: Path) -> Path:
    """Three PolyASite-shaped records: support 5.6, 44.4 and 88.9 percent."""
    rows = [
        ("chr1", "1000", "1001", "chr1:1001:+", "0.0001", "+", "5.6", "1", "20", "DI", "NaN"),
        ("chr1", "5000", "5001", "chr1:5001:+", "0.0500", "+", "44.4", "1", "80", "TE", "AAUAAA"),
        ("chr1", "9000", "9001", "chr1:9001:-", "0.2000", "-", "88.9", "1", "99", "TE", "AUUAAA"),
    ]
    path = tmp_path / "atlas.bed.gz"
    with gzip.open(path, "wt") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return path


class TestAtlasSupportFilter:
    def test_default_keeps_every_record(self, tmp_path):
        sites = _load_polya_sites(str(_write_atlas(tmp_path)), chrom_style="ucsc")
        assert sum(len(v) for v in sites.values()) == 3

    def test_threshold_drops_low_support_records(self, tmp_path):
        sites = _load_polya_sites(
            str(_write_atlas(tmp_path)), chrom_style="ucsc", min_support=25.0
        )
        assert sum(len(v) for v in sites.values()) == 2
        assert 1000 not in sites.get(("chr1", "+"), [])

    def test_plain_bed3_is_unaffected_by_the_threshold(self, tmp_path):
        path = tmp_path / "plain.bed"
        path.write_text("chr1\t100\t101\nchr1\t200\t201\n")
        sites = _load_polya_sites(str(path), chrom_style="ucsc", min_support=50.0)
        assert sum(len(v) for v in sites.values()) == 2

    def test_pas_hexamer_column_is_read(self, tmp_path):
        pas: dict = {}
        sites = _load_polya_sites(
            str(_write_atlas(tmp_path)), chrom_style="ucsc", pas_out=pas
        )
        assert sum(len(v) for v in sites.values()) == 3
        # only the two records naming a hexamer
        assert sum(len(v) for v in pas.values()) == 2
        assert 1000 not in pas.get(("chr1", "+"), [])

    def test_pas_subset_survives_an_earlier_load_without_it(self, tmp_path):
        """discover and run-plate workers load without pas_out and write the cache."""
        atlas = str(_write_atlas(tmp_path))
        _load_polya_sites(atlas, chrom_style="ucsc")
        pas: dict = {}
        _load_polya_sites(atlas, chrom_style="ucsc", pas_out=pas)
        assert sum(len(v) for v in pas.values()) == 2

    def test_default_is_off_so_upgrading_moves_no_number(self):
        assert POLYA_MIN_SUPPORT_DEFAULT == 0.0
        assert PRIMING_CRITERION_DEFAULT == "a_run"


class TestAtlasDensity:
    def test_density_is_measured_against_the_reference(self):
        # two sites 50 bp wide either side on a 10 kb contig -> ~2%
        sites = {("chr1", "+"): [1000, 5000]}
        frac = _polya_atlas_density(sites, {"chr1": 10_000}, proximity=50)
        assert frac == pytest.approx(200 / 10_000, abs=1e-6)

    def test_overlapping_windows_are_not_double_counted(self):
        sites = {("chr1", "+"): [1000, 1010, 1020]}
        frac = _polya_atlas_density(sites, {"chr1": 10_000}, proximity=50)
        assert frac == pytest.approx(120 / 10_000, abs=1e-6)

    def test_dense_atlas_records_a_warning(self):
        class _Meta:
            reference_lengths = {"chr1": 1000}
            warnings: list = []

        meta = _Meta()
        meta.warnings = []
        dense = {("chr1", "+"): list(range(0, 1000, 20))}
        _warn_if_atlas_too_dense(dense, meta, 0.0)
        assert meta.warnings and "polyA site" in meta.warnings[0]

    def test_sparse_atlas_is_silent(self):
        class _Meta:
            reference_lengths = {"chr1": 10_000_000}
            warnings: list = []

        meta = _Meta()
        meta.warnings = []
        _warn_if_atlas_too_dense({("chr1", "+"): [5000]}, meta, 25.0)
        assert meta.warnings == []


class TestPrimingCriterionOption:
    def test_unknown_criterion_is_rejected(self):
        from scnoisemeter.modules.annotation import AnnotationIndex
        with pytest.raises(ValueError, match="a_run|perc_a|either"):
            ReadClassifier(AnnotationIndex.__new__(AnnotationIndex),
                           priming_criterion="percent")
