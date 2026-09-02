"""
tests/test_regressions_v081.py
==============================
Regression guards for the v0.8.1 correctness and performance fixes.

Each test here pins behaviour that was previously wrong or that a plausible
future edit would silently break.  They are deliberately cheap: no external
annotation, no BAM.
"""

from __future__ import annotations

import random
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from scnoisemeter.cli import _check_site_contig_overlap, _harmonise_chrom_style
from scnoisemeter.modules.annotation import _extract_splice_junctions, _extract_splice_sites
from scnoisemeter.modules.intergenic_profiler import (
    IntergenicReadRecord,
    _near_polya_site,
    _repeat_overlap_fraction,
    _repeat_overlap_fraction_merged,
    merge_repeat_intervals,
)
from scnoisemeter.modules.metrics import SampleMetrics, to_multiqc_json

# ---------------------------------------------------------------------------
# Categorical groupby in the annotation build
# ---------------------------------------------------------------------------

def _categorical_exon_frame() -> pd.DataFrame:
    """Exon frame with Categorical Chromosome/Strand, as pr.read_gtf returns."""
    rows = [
        # chr1 + : two-exon transcript -> one intron 2000-3000
        ("chr1", 1000, 2000, "+", "TX1"),
        ("chr1", 3000, 4000, "+", "TX1"),
        # chr2 - : three-exon transcript -> introns 6000-7000 and 8000-9000
        ("chr2", 5000, 6000, "-", "TX2"),
        ("chr2", 7000, 8000, "-", "TX2"),
        ("chr2", 9000, 9500, "-", "TX2"),
        # single-exon transcript contributes no junction
        ("chr3", 100, 200, "+", "TX3"),
    ]
    df = pd.DataFrame(rows, columns=["Chromosome", "Start", "End", "Strand", "transcript_id"])
    df["Chromosome"] = df["Chromosome"].astype("category")
    df["Strand"] = df["Strand"].astype("category")
    return df


class TestSpliceJunctionExtraction:
    """`observed=False` on Categorical keys iterated the category product.

    The bug was pure cost, not a wrong answer, so the guard is that the junction
    set is exactly right *and* that adding unused categories cannot change it.
    """

    EXPECTED = {
        "chr1": {(2000, 3000, "+")},
        "chr2": {(6000, 7000, "-"), (8000, 9000, "-")},
    }

    def test_exact_junctions_with_categorical_keys(self):
        assert _extract_splice_junctions(_categorical_exon_frame()) == self.EXPECTED

    def test_unused_categories_do_not_change_the_result(self):
        df = _categorical_exon_frame()
        # Contigs present in the dtype but carrying no exon rows: this is the
        # normal state after subsetting a genome-wide GTF.
        df["Chromosome"] = df["Chromosome"].cat.add_categories(
            [f"chr{i}" for i in range(4, 23)]
        )
        assert _extract_splice_junctions(df) == self.EXPECTED

    def test_single_exon_transcript_yields_no_junction(self):
        assert "chr3" not in _extract_splice_junctions(_categorical_exon_frame())

    def test_splice_sites_cover_both_exon_boundaries(self):
        sites = _extract_splice_sites(_categorical_exon_frame())
        assert sites["chr1"] == {(1000, "+"), (2000, "+"), (3000, "+"), (4000, "+")}
        assert (9500, "-") in sites["chr2"]
        # keys are plain strings even though the column is Categorical
        assert all(isinstance(k, str) for k in sites)

    def test_missing_transcript_id_returns_empty(self):
        df = _categorical_exon_frame().drop(columns=["transcript_id"])
        assert _extract_splice_junctions(df) == {}


# ---------------------------------------------------------------------------
# Repeat overlap: merged fast path vs exhaustive reference
# ---------------------------------------------------------------------------

def _reference_overlap_fraction(ivls, records) -> float:
    """Exhaustive implementation the fast path must reproduce."""
    overlap = total = 0
    for record in records:
        start, end = int(record.start), int(record.end)
        total += max(end - start, 0)
        covered_end = start
        for ivl_start, ivl_end in sorted(ivls):
            if ivl_start >= end:
                break
            lo, hi = max(start, ivl_start), min(end, ivl_end)
            if hi > lo:
                lo = max(lo, covered_end)
                if hi > lo:
                    overlap += hi - lo
                    covered_end = hi
    return overlap / total if total else 0.0


def _rec(start: int, end: int) -> IntergenicReadRecord:
    return IntergenicReadRecord("chr1", start, end, "+", "CB", False, end)


class TestMergeRepeatIntervals:

    def test_merges_overlapping_and_adjacent(self):
        merged = merge_repeat_intervals({"chr1": [(150, 200), (100, 180), (200, 250)]})
        assert merged["chr1"] == [(100, 250)]

    def test_drops_empty_and_inverted_intervals(self):
        merged = merge_repeat_intervals({"chr1": [(100, 100), (300, 200), (400, 500)]})
        assert merged["chr1"] == [(400, 500)]

    def test_is_idempotent(self):
        once = merge_repeat_intervals({"chr1": [(10, 20), (15, 30), (100, 110)]})
        assert merge_repeat_intervals(once) == once

    def test_preserves_contig_keys(self):
        merged = merge_repeat_intervals({"chr1": [(1, 2)], "GL000009.2": [(3, 4)]})
        assert set(merged) == {"chr1", "GL000009.2"}


class TestRepeatOverlapFastPath:

    def test_matches_reference_on_randomised_inputs(self):
        rng = random.Random(20240917)
        for _ in range(500):
            pos = 0
            disjoint = []
            for _ in range(rng.randint(1, 30)):
                pos += rng.randint(1, 30)
                disjoint.append((pos, pos + rng.randint(1, 25)))
                pos += 25
            overlapping = [
                (s, s + rng.randint(1, 60))
                for s in (rng.randint(0, 400) for _ in range(rng.randint(1, 25)))
            ]
            records = [
                _rec(s, s + rng.randint(1, 120))
                for s in sorted(rng.randint(0, max(pos - 1, 1)) for _ in range(rng.randint(1, 5)))
            ]
            for ivls in (disjoint, overlapping):
                expected = _reference_overlap_fraction(ivls, records)
                assert _repeat_overlap_fraction({"chr1": ivls}, "chr1", records) == (
                    pytest.approx(expected)
                )

    def test_wrapper_accepts_unsorted_input(self):
        ivls = [(1000, 1100), (100, 200), (500, 600)]
        records = [_rec(150, 550)]
        assert _repeat_overlap_fraction({"chr1": ivls}, "chr1", records) == pytest.approx(
            _reference_overlap_fraction(ivls, records)
        )

    def test_read_far_along_the_contig_is_scored(self):
        """The old scan started at index 0; the fast path binary-searches."""
        ivls = [(i * 1000, i * 1000 + 500) for i in range(1, 5001)]
        merged = merge_repeat_intervals({"chr1": ivls})
        records = [_rec(4_000_000, 4_000_500)]
        assert _repeat_overlap_fraction_merged(merged, "chr1", records) == pytest.approx(1.0)

    def test_interval_spanning_the_read_start_is_not_missed(self):
        # The only interval starts before the read and ends after it, so the
        # binary search must step back one position.
        merged = merge_repeat_intervals({"chr1": [(0, 10_000)]})
        assert _repeat_overlap_fraction_merged(
            merged, "chr1", [_rec(5_000, 5_100)]
        ) == pytest.approx(1.0)

    def test_unknown_contig_and_empty_inputs(self):
        merged = merge_repeat_intervals({"chr1": [(0, 100)]})
        assert _repeat_overlap_fraction_merged(merged, "chr99", [_rec(0, 50)]) == 0.0
        assert _repeat_overlap_fraction_merged(merged, "chr1", []) == 0.0
        assert _repeat_overlap_fraction(None, "chr1", [_rec(0, 50)]) == 0.0


# ---------------------------------------------------------------------------
# Undefined metrics are absent, not zero
# ---------------------------------------------------------------------------

class TestUndefinedMetricsAreAbsent:
    """A metric with no denominator must not be rendered as a measurement of 0.

    The cohort reader already treats a blank TSV value as absent, so the only
    thing that was missing is SampleMetrics actually leaving them unset.
    """

    @staticmethod
    def _blank() -> SampleMetrics:
        return SampleMetrics(
            sample_name="s", bam_path="s.bam", platform="unknown",
            pipeline_stage="custom", aligner="",
        )

    @pytest.mark.parametrize("field", [
        "unmapped_read_frac",
        "low_mapq_read_frac",
        "per_cell_noise_median",
        "per_cell_noise_iqr",
    ])
    def test_default_is_none_not_zero(self, field):
        assert getattr(self._blank(), field) is None

    def test_multiqc_json_omits_unmeasured_fields(self):
        data = to_multiqc_json(self._blank())["data"]["s"]
        for field in ("unmapped_read_frac", "low_mapq_read_frac",
                      "per_cell_noise_median", "per_cell_noise_iqr"):
            assert field not in data

    def test_multiqc_json_includes_measured_fields(self):
        sm = self._blank()
        sm.unmapped_read_frac = 0.125
        sm.per_cell_noise_median = 0.0
        data = to_multiqc_json(sm)["data"]["s"]
        assert data["unmapped_read_frac"] == pytest.approx(0.125)
        # a genuine zero must survive; only *absent* is suppressed
        assert data["per_cell_noise_median"] == pytest.approx(0.0)
        assert "low_mapq_read_frac" not in data


# ---------------------------------------------------------------------------
# Site-atlas chromosome naming (v0.9.0-dev)
# ---------------------------------------------------------------------------

class TestHarmoniseChromStyle:
    """The old helper only stripped `chr`, so an Ensembl-named atlas against a
    UCSC-named BAM matched nothing while staying non-empty."""

    ENSEMBL_ATLAS = {("1", "+"): [100, 200], ("MT", "-"): [50], ("X", "+"): [7]}
    UCSC_ATLAS = {("chr1", "+"): [100, 200], ("chrM", "-"): [50]}

    def test_adds_prefix_for_ucsc_bam(self):
        out = _harmonise_chrom_style(self.ENSEMBL_ATLAS, "ucsc")
        assert set(out) == {("chr1", "+"), ("chrM", "-"), ("chrX", "+")}
        assert out[("chr1", "+")] == [100, 200]

    def test_strips_prefix_for_ensembl_bam(self):
        out = _harmonise_chrom_style(self.UCSC_ATLAS, "ensembl")
        # chrM -> MT, not M: Ensembl and GENCODE spell the mitochondrion MT, so
        # a plain prefix strip yields a key that never matches the BAM.
        assert set(out) == {("1", "+"), ("MT", "-")}

    def test_mito_contig_uses_each_style_own_spelling(self):
        assert ("chrM", "-") in _harmonise_chrom_style({("MT", "-"): [1]}, "ucsc")
        assert ("chrMT", "-") not in _harmonise_chrom_style({("MT", "-"): [1]}, "ucsc")
        assert ("MT", "-") in _harmonise_chrom_style({("chrM", "-"): [1]}, "ensembl")
        assert ("M", "-") not in _harmonise_chrom_style({("chrM", "-"): [1]}, "ensembl")

    @pytest.mark.parametrize("style,other,contigs", [
        ("ucsc", "ensembl", ["1", "MT", "X"]),
        ("ensembl", "ucsc", ["chr1", "chrM", "chrX"]),
    ])
    def test_round_trip_is_identity(self, style, other, contigs):
        """Harmonising to the other style and back must recover the input.

        This is the property that catches an asymmetric special case: without
        chrM -> MT, MT survives one hop and is lost on the return.
        """
        original = {(c, "+"): [1] for c in contigs}
        there = _harmonise_chrom_style(original, style)
        assert _harmonise_chrom_style(there, other) == original

    def test_mito_lookup_succeeds_in_both_directions(self):
        """The end-to-end consequence: the proximity test must fire on the mito
        contig whichever naming the BAM uses."""
        ucsc_bam = _harmonise_chrom_style({("MT", "-"): [500]}, "ucsc")
        assert _near_polya_site(ucsc_bam, "chrM", 505, proximity=50, strand="-")
        ensembl_bam = _harmonise_chrom_style({("chrM", "-"): [500]}, "ensembl")
        assert _near_polya_site(ensembl_bam, "MT", 505, proximity=50, strand="-")

    def test_noop_when_styles_already_agree(self):
        assert _harmonise_chrom_style(self.UCSC_ATLAS, "ucsc") == self.UCSC_ATLAS
        assert _harmonise_chrom_style(self.ENSEMBL_ATLAS, "ensembl") == self.ENSEMBL_ATLAS

    def test_noop_on_unknown_style_or_empty(self):
        assert _harmonise_chrom_style(self.ENSEMBL_ATLAS, "unknown") == self.ENSEMBL_ATLAS
        assert _harmonise_chrom_style({}, "ucsc") == {}

    def test_plain_string_keys_are_supported(self):
        # BED3 / legacy cache representation carries no strand in the key
        assert set(_harmonise_chrom_style({"1": [5], "2": [6]}, "ucsc")) == {"chr1", "chr2"}

    def test_lookup_succeeds_after_harmonisation(self):
        """The end-to-end point: the proximity test must actually fire."""
        atlas = _harmonise_chrom_style(self.ENSEMBL_ATLAS, "ucsc")
        assert _near_polya_site(atlas, "chr1", 105, proximity=50, strand="+")
        # and would not have, before harmonisation
        assert not _near_polya_site(self.ENSEMBL_ATLAS, "chr1", 105, proximity=50, strand="+")


class TestSiteContigOverlapGuard:

    def test_discards_atlas_sharing_no_contig(self):
        sites = {("1", "+"): [100]}
        assert _check_site_contig_overlap(
            sites, ["chr1", "chr2"], label="polyA site atlas"
        ) is None

    def test_keeps_atlas_with_any_shared_contig(self):
        sites = {("chr1", "+"): [100], ("weird_contig", "+"): [5]}
        assert _check_site_contig_overlap(
            sites, ["chr1", "chr2"], label="polyA site atlas"
        ) is sites

    def test_records_a_warning_on_the_metadata(self):
        class _Meta:
            warnings: list = []
        meta = _Meta()
        meta.warnings = []
        _check_site_contig_overlap(
            {("1", "+"): [100]}, ["chr1"], label="polyA site atlas", meta=meta
        )
        assert len(meta.warnings) == 1
        assert "polyA site atlas" in meta.warnings[0]
        # the message must name both namings so the user can act on it
        assert "'1'" in meta.warnings[0] and "'chr1'" in meta.warnings[0]

    def test_passthrough_for_absent_or_empty_inputs(self):
        assert _check_site_contig_overlap(None, ["chr1"], label="x") is None
        assert _check_site_contig_overlap({}, ["chr1"], label="x") == {}
        # no reference names to compare against: cannot judge, so keep
        sites = {("1", "+"): [1]}
        assert _check_site_contig_overlap(sites, [], label="x") is sites

    def test_near_polya_site_is_none_safe(self):
        assert _near_polya_site(None, "chr1", 100) is False


# ---------------------------------------------------------------------------
# Version strings must not drift apart
# ---------------------------------------------------------------------------

class TestVersionConsistency:
    """0.8.1 bumped pyproject and __init__ but not the four doc banners.

    Nothing linked them mechanically, so the docs shipped claiming 0.8.0 while
    the package reported 0.8.1. These assertions make the next bump fail loudly
    if any one of them is forgotten.
    """

    ROOT = Path(__file__).parent.parent

    # file -> pattern capturing the version in that file's banner
    DOC_BANNERS = {
        "README.md": r"Version \*\*(\d+\.\d+\.\d+)\*\* deliberately avoids",
        "docs/documentation.md": r"Version (\d+\.\d+\.\d+) operational reference\.",
        "docs/methods_annotated.md": r"companion for version (\d+\.\d+\.\d+)\*\*",
        "docs/whitepaper.md": r"Version (\d+\.\d+\.\d+)\. This document explains",
    }

    @staticmethod
    def _package_version() -> str:
        from scnoisemeter import __version__
        return __version__

    def test_pyproject_matches_dunder_version(self):
        text = (self.ROOT / "pyproject.toml").read_text()
        match = re.search(r'^version = "(\d+\.\d+\.\d+)"', text, re.MULTILINE)
        assert match, "no version field found in pyproject.toml"
        assert match.group(1) == self._package_version()

    @pytest.mark.parametrize("rel_path", sorted(DOC_BANNERS))
    def test_doc_banner_matches_package_version(self, rel_path):
        text = (self.ROOT / rel_path).read_text()
        match = re.search(self.DOC_BANNERS[rel_path], text)
        assert match, (
            f"{rel_path}: version banner not found. If the wording changed, "
            f"update DOC_BANNERS in this test rather than deleting the check."
        )
        assert match.group(1) == self._package_version(), (
            f"{rel_path} claims version {match.group(1)} but the package "
            f"reports {self._package_version()}"
        )

    def test_changelog_has_a_section_for_the_current_version(self):
        text = (self.ROOT / "CHANGELOG.md").read_text()
        assert f"## {self._package_version()}" in text, (
            f"CHANGELOG.md has no '## {self._package_version()}' section"
        )
