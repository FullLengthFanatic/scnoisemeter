"""
tests/test_discover.py
======================
Unit tests for the discover subcommand BAM inspector.

All tests use mock BAM headers — no real BAMs, no network.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch, PropertyMock

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from scnoisemeter.constants import Platform, PipelineStage, Chemistry
from scnoisemeter.utils.discover_inspector import (
    DiscoverBamInfo,
    _detect_chemistry,
    _detect_paired_end,
    inspect_bam_for_discover,
    format_discovery_table,
)
from scnoisemeter.utils.bam_inspector import BamMetadata


# ---------------------------------------------------------------------------
# Helpers: build mock BamMetadata with known platform/stage
# ---------------------------------------------------------------------------

def _make_meta(
    platform=Platform.UNKNOWN,
    pipeline_stage=PipelineStage.CUSTOM,
    sort_order="coordinate",
    aligner="",
    pipeline_hints=None,
    bam_path="/fake/file.bam",
):
    meta = BamMetadata(path=Path(bam_path))
    meta.platform = platform
    meta.pipeline_stage = pipeline_stage
    meta.sort_order = sort_order
    meta.aligner = aligner
    meta.pipeline_hints = pipeline_hints or []
    return meta


def _make_info(meta, has_index=True, chemistry="unknown", n_reads=2000):
    info = DiscoverBamInfo(meta=meta)
    info.has_index = has_index
    info.chemistry = chemistry
    info.n_reads_approx = n_reads
    return info


# ---------------------------------------------------------------------------
# Platform detection — via inspect_bam / bam_inspector
# (these test that the existing bam_inspector rules produce the right platform
#  — we mock inspect_bam to return known results)
# ---------------------------------------------------------------------------

class TestPlatformDetection:
    """
    Test that the existing platform detection in bam_inspector works correctly
    for the scenarios we care about in discover.  We mock inspect_bam to return
    pre-built BamMetadata objects with the platform already set.
    """

    def _run_discover(self, meta, has_index=True, tmp_path=None):
        """Helper: call inspect_bam_for_discover with a mocked inspect_bam."""
        bam_path = tmp_path / "test.bam" if tmp_path else Path("/fake/test.bam")
        bai_path = Path(str(bam_path) + ".bai")

        with (
            patch("scnoisemeter.utils.discover_inspector.inspect_bam", return_value=meta),
            patch.object(Path, "exists", return_value=has_index),
            patch("scnoisemeter.utils.discover_inspector._detect_chemistry"),
            patch("scnoisemeter.utils.discover_inspector._detect_paired_end"),
            patch("scnoisemeter.utils.discover_inspector._get_approx_read_count"),
        ):
            return inspect_bam_for_discover(bam_path)

    def test_pbmm2_pg_gives_pacbio(self, tmp_path):
        """pbmm2 in @PG → pacbio."""
        meta = _make_meta(
            platform=Platform.PACBIO,
            aligner="pbmm2",
            pipeline_hints=["pbmm2"],
        )
        info = self._run_discover(meta, tmp_path=tmp_path)
        assert info.platform == Platform.PACBIO

    def test_minimap2_epi2me_gives_ont(self, tmp_path):
        """minimap2 + epi2me in @CO → ont."""
        meta = _make_meta(
            platform=Platform.ONT,
            aligner="minimap2",
            pipeline_hints=["minimap2", "wf-single-cell"],
        )
        info = self._run_discover(meta, tmp_path=tmp_path)
        assert info.platform == Platform.ONT

    def test_minimap2_without_epi2me_gives_ont(self, tmp_path):
        """minimap2 without epi2me marker → ont (generic)."""
        meta = _make_meta(
            platform=Platform.ONT,
            aligner="minimap2",
            pipeline_hints=["minimap2"],
        )
        info = self._run_discover(meta, tmp_path=tmp_path)
        assert info.platform == Platform.ONT

    def test_star_with_starsolo_gives_illumina(self, tmp_path):
        """STAR + STARsolo in @PG → illumina_10x."""
        meta = _make_meta(
            platform=Platform.ILLUMINA_10X,
            aligner="STAR",
            pipeline_hints=["STAR", "STARsolo"],
        )
        info = self._run_discover(meta, tmp_path=tmp_path)
        assert info.platform in (Platform.ILLUMINA, Platform.ILLUMINA_10X)

    def test_star_alone_gives_smartseq(self, tmp_path):
        """Bare STAR (no STARsolo/cellranger) → smartseq."""
        meta = _make_meta(
            platform=Platform.SMARTSEQ,
            aligner="STAR",
            pipeline_hints=["STAR"],
        )
        info = self._run_discover(meta, tmp_path=tmp_path)
        assert info.platform == Platform.SMARTSEQ

    def test_cellranger_cl_path_gives_illumina(self, tmp_path):
        """cellranger in CL path → illumina."""
        meta = _make_meta(
            platform=Platform.ILLUMINA_10X,
            aligner="cellranger",
            pipeline_hints=["cellranger"],
        )
        info = self._run_discover(meta, tmp_path=tmp_path)
        assert info.platform in (Platform.ILLUMINA, Platform.ILLUMINA_10X)


# ---------------------------------------------------------------------------
# _refine_platform_from_pipeline_hints — unit tests for Smart-seq detection
# ---------------------------------------------------------------------------

class TestRefinePlatformHints:
    """Direct tests for the bam_inspector platform-refinement logic."""

    def _make_meta_unknown(self, hints):
        from scnoisemeter.utils.bam_inspector import BamMetadata, _refine_platform_from_pipeline_hints
        meta = BamMetadata(path=Path("/fake.bam"))
        meta.pipeline_hints = hints
        meta.platform = Platform.UNKNOWN
        _refine_platform_from_pipeline_hints(meta)
        return meta

    def test_star_alone_gives_smartseq(self):
        meta = self._make_meta_unknown(["STAR"])
        assert meta.platform == Platform.SMARTSEQ

    def test_star_with_starsolo_gives_illumina_10x(self):
        meta = self._make_meta_unknown(["STAR", "STARsolo"])
        assert meta.platform == Platform.ILLUMINA_10X

    def test_star_with_cellranger_gives_illumina_10x(self):
        meta = self._make_meta_unknown(["cellranger", "STAR"])
        assert meta.platform == Platform.ILLUMINA_10X

    def test_cellranger_alone_gives_illumina_10x(self):
        meta = self._make_meta_unknown(["cellranger"])
        assert meta.platform == Platform.ILLUMINA_10X

    def test_bd_signals_give_illumina_bd(self):
        meta = self._make_meta_unknown(["STARsolo", "rhapsody"])
        assert meta.platform == Platform.ILLUMINA_BD

    def test_wf_single_cell_gives_ont(self):
        meta = self._make_meta_unknown(["minimap2", "wf-single-cell"])
        assert meta.platform == Platform.ONT

    def test_pbmm2_gives_pacbio(self):
        meta = self._make_meta_unknown(["pbmm2"])
        assert meta.platform == Platform.PACBIO


# ---------------------------------------------------------------------------
# Pipeline stage detection
# ---------------------------------------------------------------------------

class TestPipelineStageDetection:
    def test_cb_present_high_fraction_gives_post_filter(self, tmp_path):
        """CB in >50% reads → post_filter (set by inspect_bam, checked here)."""
        meta = _make_meta(
            platform=Platform.ONT,
            pipeline_stage=PipelineStage.POST_FILTER,
        )
        meta.barcode_fraction = 0.92
        meta.barcode_aware = True

        with (
            patch("scnoisemeter.utils.discover_inspector.inspect_bam", return_value=meta),
            patch.object(Path, "exists", return_value=True),
            patch("scnoisemeter.utils.discover_inspector._detect_chemistry"),
            patch("scnoisemeter.utils.discover_inspector._detect_paired_end"),
            patch("scnoisemeter.utils.discover_inspector._get_approx_read_count"),
        ):
            info = inspect_bam_for_discover(tmp_path / "test.bam")

        assert info.pipeline_stage == PipelineStage.POST_FILTER

    def test_cb_absent_gives_pre_filter(self, tmp_path):
        """CB absent → pre_filter."""
        meta = _make_meta(
            platform=Platform.PACBIO,
            pipeline_stage=PipelineStage.PRE_FILTER,
        )
        meta.barcode_fraction = 0.0
        meta.barcode_aware = False

        with (
            patch("scnoisemeter.utils.discover_inspector.inspect_bam", return_value=meta),
            patch.object(Path, "exists", return_value=True),
            patch("scnoisemeter.utils.discover_inspector._detect_chemistry"),
            patch("scnoisemeter.utils.discover_inspector._detect_paired_end"),
            patch("scnoisemeter.utils.discover_inspector._get_approx_read_count"),
        ):
            info = inspect_bam_for_discover(tmp_path / "test.bam")

        assert info.pipeline_stage == PipelineStage.PRE_FILTER


# ---------------------------------------------------------------------------
# Blocking issue detection
# ---------------------------------------------------------------------------

class TestBlockingIssues:
    def test_queryname_sort_is_flagged(self, tmp_path):
        """queryname sort → run_issues contains sort warning."""
        meta = _make_meta(platform=Platform.ONT, sort_order="queryname")

        with (
            patch("scnoisemeter.utils.discover_inspector.inspect_bam", return_value=meta),
            patch.object(Path, "exists", return_value=True),
            patch("scnoisemeter.utils.discover_inspector._detect_chemistry"),
            patch("scnoisemeter.utils.discover_inspector._detect_paired_end"),
            patch("scnoisemeter.utils.discover_inspector._get_approx_read_count"),
        ):
            info = inspect_bam_for_discover(tmp_path / "test.bam")

        assert not info.can_run
        assert any("sort" in issue for issue in info.run_issues)

    def test_no_index_skipped_with_warning(self, tmp_path):
        """No .bai → run_issues contains index warning, inspect_bam not called."""
        bam_path = tmp_path / "test.bam"
        bam_path.write_bytes(b"fake")  # file exists

        with patch("scnoisemeter.utils.discover_inspector.inspect_bam") as mock_insp:
            info = inspect_bam_for_discover(bam_path)

        # inspect_bam should NOT have been called (no index → early return)
        mock_insp.assert_not_called()
        assert not info.can_run
        assert any("index" in issue for issue in info.run_issues)

    def test_unknown_platform_is_blocking(self, tmp_path):
        """Unknown platform → run_issues in --run-all mode."""
        meta = _make_meta(platform=Platform.UNKNOWN, sort_order="coordinate")

        with (
            patch("scnoisemeter.utils.discover_inspector.inspect_bam", return_value=meta),
            patch.object(Path, "exists", return_value=True),
            patch("scnoisemeter.utils.discover_inspector._detect_chemistry"),
            patch("scnoisemeter.utils.discover_inspector._detect_paired_end"),
            patch("scnoisemeter.utils.discover_inspector._get_approx_read_count"),
        ):
            info = inspect_bam_for_discover(tmp_path / "test.bam")

        assert not info.can_run
        assert any("platform" in issue for issue in info.run_issues)


# ---------------------------------------------------------------------------
# Chemistry detection
# ---------------------------------------------------------------------------

class TestChemistryDetection:
    def _make_info_for_chemistry(self, header: dict) -> DiscoverBamInfo:
        meta = _make_meta()
        info = _make_info(meta)
        with patch("pysam.AlignmentFile") as mock_af:
            mock_bam = MagicMock()
            mock_bam.header.to_dict.return_value = header
            mock_af.return_value.__enter__.return_value = mock_bam
            _detect_chemistry(Path("/fake/test.bam"), info)
        return info

    def test_10x_v3_in_co(self):
        header = {"CO": ["10x Chromium v3 library"]}
        info = self._make_info_for_chemistry(header)
        assert info.chemistry == Chemistry.TENX_V3.value

    def test_10x_v4_gem_x_in_co(self):
        header = {"CO": ["GEM-X library processing"]}
        info = self._make_info_for_chemistry(header)
        assert info.chemistry == Chemistry.TENX_V4.value

    def test_bd_rhapsody_in_co(self):
        header = {"CO": ["BD Rhapsody WTA library"]}
        info = self._make_info_for_chemistry(header)
        assert info.chemistry == Chemistry.BD_RHAPSODY_WTA.value

    def test_unknown_header_gives_unknown(self):
        header = {"CO": ["some random pipeline"]}
        info = self._make_info_for_chemistry(header)
        assert info.chemistry == "unknown"

    def test_10x_in_pg_cl(self):
        header = {
            "PG": [{"ID": "samtools", "CL": "/path/to/10x_cellranger/bin/samtools sort"}]
        }
        info = self._make_info_for_chemistry(header)
        assert info.chemistry in (Chemistry.TENX_V3.value, Chemistry.TENX_V4.value)


# ---------------------------------------------------------------------------
# Table formatting (smoke test)
# ---------------------------------------------------------------------------

class TestFormatTable:
    def test_empty_infos(self):
        result = format_discovery_table([])
        assert "no bam" in result.lower()

    def test_single_row(self):
        meta = _make_meta(
            platform=Platform.ONT,
            pipeline_stage=PipelineStage.POST_FILTER,
            sort_order="coordinate",
        )
        info = _make_info(meta, has_index=True, n_reads=4200)
        info.chemistry = "10x_v3"
        table = format_discovery_table([info])
        assert "ont" in table
        assert "post_filter" in table
        assert "4.2k" in table

    def test_unknown_platform_flagged(self):
        meta = _make_meta(platform=Platform.UNKNOWN, sort_order="coordinate")
        info = _make_info(meta)
        table = format_discovery_table([info])
        assert "⚠" in table

    def test_queryname_sort_flagged(self):
        meta = _make_meta(platform=Platform.ONT, sort_order="queryname")
        info = _make_info(meta)
        table = format_discovery_table([info])
        assert "⚠" in table
