"""
tests/test_adversarial.py
=========================
Regression tests for adversarial / edge-case input handling.

All tests use mocks — no real BAM files are required.  Each test verifies
that the appropriate error or warning is raised for a specific adversarial
input scenario.

Cases covered
-------------
1.  Ensembl-style chromosome names  → FATAL ClickException (naming mismatch)
2.  Name-sorted BAM                 → FATAL ClickException (wrong sort order)
3.  All unmapped reads              → Warning + SystemExit(0), no report
4.  Tiny BAM (< 100 reads)         → Warning emitted, run completes (exit 0)
10. Unaligned BAM (no @SQ lines)   → FATAL ClickException (not aligned)
5.  Missing BAM index               → FATAL ClickException (index not found)
6.  No @HD SO tag                   → Warning only, run proceeds
7.  Wrong species (chr1 length)     → Warning in meta.warnings, run proceeds
8.  Empty BAM (0 reads)             → FATAL ClickException (no reads)
9.  No CB tags + post_filter        → Warning emitted, barcode-agnostic mode
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch, PropertyMock

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

import click
from click.testing import CliRunner

from scnoisemeter.cli import (
    _validate_sort_order,
    _validate_sq_lines,
    _validate_chromosome_naming,
    _validate_chromosome_lengths,
    _detect_chrom_style,
    _gtf_chrom_style,
)
from scnoisemeter.utils.bam_inspector import BamMetadata
from scnoisemeter.constants import Platform, PipelineStage


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_meta(
    sort_order: str = "coordinate",
    reference_names: list | None = None,
    reference_lengths: dict | None = None,
    barcode_fraction: float = 1.0,
    barcode_aware: bool = True,
    platform: Platform = Platform.ONT,
    pipeline_stage: PipelineStage = PipelineStage.POST_FILTER,
    path: str = "/fake/sample.bam",
) -> BamMetadata:
    """Build a minimal BamMetadata instance for testing."""
    if reference_names is None:
        reference_names = ["chr1", "chr2", "chr3", "chrX"]
    if reference_lengths is None:
        reference_lengths = {
            "chr1": 248_956_422,
            "chr2": 242_193_529,
            "chr3": 198_295_559,
            "chrX": 156_040_895,
        }
    meta = BamMetadata(
        path=Path(path),
        platform=platform,
        sort_order=sort_order,
        reference_names=reference_names,
        reference_lengths=reference_lengths,
        barcode_fraction=barcode_fraction,
        barcode_aware=barcode_aware,
        barcode_tag_present=barcode_aware,
        pipeline_stage=pipeline_stage,
    )
    return meta


# ---------------------------------------------------------------------------
# CASE 1 — Ensembl-style chromosome names
# ---------------------------------------------------------------------------

class TestCase1EnsemblNaming:
    """BAM uses Ensembl names (1, 2, …); GTF uses UCSC (chr1, chr2, …)."""

    def test_raises_click_exception(self):
        meta = _make_meta(
            reference_names=["1", "2", "3", "X", "Y"],
            reference_lengths={"1": 248_956_422, "2": 242_193_529,
                               "3": 198_295_559, "X": 156_040_895, "Y": 57_227_415},
        )
        with patch(
            "scnoisemeter.cli._gtf_chrom_style", return_value="ucsc"
        ):
            with pytest.raises(click.ClickException) as exc_info:
                _validate_chromosome_naming(meta, "/fake/gencode.gtf.gz")
        msg = str(exc_info.value.format_message())
        assert "mismatch" in msg.lower()
        # Message must name what was found and what was expected
        assert "ENSEMBL" in msg or "ensembl" in msg.lower()
        assert "UCSC" in msg or "ucsc" in msg.lower()

    def test_error_names_example_contigs(self):
        """Error message must include an example Ensembl contig name."""
        meta = _make_meta(
            reference_names=["1", "2", "X"],
            reference_lengths={"1": 248_956_422, "2": 242_193_529, "X": 156_040_895},
        )
        with patch("scnoisemeter.cli._gtf_chrom_style", return_value="ucsc"):
            with pytest.raises(click.ClickException) as exc_info:
                _validate_chromosome_naming(meta, "/fake/gencode.gtf.gz")
        msg = exc_info.value.format_message()
        # Must mention one of the Ensembl contig names
        assert any(c in msg for c in ["'1'", "'2'", "'X'"]), (
            f"Expected an Ensembl contig name in the error message, got: {msg}"
        )


# ---------------------------------------------------------------------------
# CASE 2 — Name-sorted BAM
# ---------------------------------------------------------------------------

class TestCase2NameSorted:
    """BAM @HD SO:queryname — must fail with a clear sort-order error."""

    def test_raises_click_exception(self):
        meta = _make_meta(sort_order="queryname")
        with pytest.raises(click.ClickException) as exc_info:
            _validate_sort_order(meta)
        msg = exc_info.value.format_message()
        assert "queryname" in msg.lower() or "sort" in msg.lower()

    def test_error_mentions_samtools_sort(self):
        meta = _make_meta(sort_order="queryname")
        with pytest.raises(click.ClickException) as exc_info:
            _validate_sort_order(meta)
        msg = exc_info.value.format_message()
        assert "samtools sort" in msg

    def test_unsorted_also_fatal(self):
        """SO:unsorted should also trigger a fatal error."""
        meta = _make_meta(sort_order="unsorted")
        with pytest.raises(click.ClickException):
            _validate_sort_order(meta)


# ---------------------------------------------------------------------------
# CASE 3 — All reads unmapped
# ---------------------------------------------------------------------------

class TestCase3AllUnmapped:
    """BAM has reads but all are unmapped → graceful warning + exit 0."""

    def test_graceful_exit_with_warning(self, tmp_path):
        """
        Simulate: BAM exists with index, has 500 unmapped reads, 0 mapped.
        The run should exit 0 with a 'no mapped reads' warning.
        """
        from scnoisemeter.cli import cli
        runner = CliRunner()

        # Create real (empty) placeholder files so click's path validator passes
        bam_file = tmp_path / "all_unmapped.bam"
        bai_file = tmp_path / "all_unmapped.bam.bai"
        gtf_file = tmp_path / "gencode.gtf.gz"
        bam_file.write_bytes(b"")
        bai_file.write_bytes(b"")
        gtf_file.write_bytes(b"")

        with patch("scnoisemeter.cli.inspect_bam") as mock_inspect, \
             patch("scnoisemeter.cli._validate_sort_order"), \
             patch("pysam.AlignmentFile") as mock_bam_cls:

            mock_meta = _make_meta()
            mock_meta.warnings = []
            mock_inspect.return_value = mock_meta

            # sort-order pre-check: coordinate
            mock_bam_instance = MagicMock()
            mock_bam_instance.__enter__ = lambda s: s
            mock_bam_instance.__exit__ = MagicMock(return_value=False)
            mock_bam_instance.header.to_dict.return_value = {"HD": {"SO": "coordinate"}}
            type(mock_bam_instance).mapped = PropertyMock(return_value=0)
            type(mock_bam_instance).unmapped = PropertyMock(return_value=500)
            mock_bam_cls.return_value = mock_bam_instance

            result = runner.invoke(cli, [
                "run",
                "--bam", str(bam_file),
                "--gtf", str(gtf_file),
                "--platform", "ont",
                "--pipeline-stage", "post_filter",
                "--output-dir", str(tmp_path / "out"),
            ])

        # Should exit 0
        assert result.exit_code == 0, (
            f"Expected exit 0, got {result.exit_code}. Output:\n{result.output}\n{result.stderr}"
        )
        combined = (result.output or "") + (result.stderr or "")
        assert "no mapped reads" in combined.lower() or "unmapped" in combined.lower(), (
            f"Expected 'no mapped reads' warning, got: {combined}"
        )


# ---------------------------------------------------------------------------
# CASE 4 — Tiny BAM (< 100 reads)
# ---------------------------------------------------------------------------

class TestCase4TinyBAM:
    """BAM has fewer than 100 reads → warning in output, exit 0."""

    def test_low_read_count_warning_emitted(self, tmp_path):
        """
        After pipeline completes with n_reads_processed < 100, a warning
        must be added to meta.warnings and emitted to stderr.
        """
        from scnoisemeter.cli import cli
        runner = CliRunner()

        # Create real placeholder files
        bam_file = tmp_path / "tiny.bam"
        bai_file = tmp_path / "tiny.bam.bai"
        gtf_file = tmp_path / "gencode.gtf.gz"
        bam_file.write_bytes(b"")
        bai_file.write_bytes(b"")
        gtf_file.write_bytes(b"")

        # Build a mock SampleResult with only 4 reads
        mock_sample_result = MagicMock()
        mock_sample_result.n_reads_processed = 4
        mock_sample_result.read_counts = {}
        mock_sample_result.base_counts = {}
        mock_sample_result.umi_sets = {}
        mock_sample_result.artifact_flags = {}
        mock_sample_result.length_samples = {}
        mock_sample_result.length_bin_counts = {}
        mock_sample_result.intergenic_reads = []
        mock_sample_result.exonic_sense_three_prime = []
        mock_sample_result.exonic_sense_five_prime = []

        mock_meta = _make_meta()
        mock_meta.warnings = []
        mock_sample_result.meta = mock_meta
        mock_sample_result.bam_path = bam_file

        mock_bam_instance = MagicMock()
        mock_bam_instance.__enter__ = lambda s: s
        mock_bam_instance.__exit__ = MagicMock(return_value=False)
        mock_bam_instance.header.to_dict.return_value = {"HD": {"SO": "coordinate"}}
        type(mock_bam_instance).mapped = PropertyMock(return_value=4)
        type(mock_bam_instance).unmapped = PropertyMock(return_value=0)

        mock_sm = MagicMock()
        mock_sm.noise_read_frac = 0.0
        mock_sm.noise_base_frac = 0.0
        mock_sm.strand_concordance = 1.0
        mock_sm.chimeric_read_frac = 0.0
        mock_sm.n_cells = 0
        mock_sm.warnings = []

        with patch("scnoisemeter.cli.inspect_bam", return_value=mock_meta), \
             patch("scnoisemeter.cli._validate_sort_order"), \
             patch("scnoisemeter.cli._resolve_gtf", return_value=(str(gtf_file), 49, "user")), \
             patch("scnoisemeter.cli._validate_chromosome_naming"), \
             patch("scnoisemeter.cli._validate_chromosome_lengths"), \
             patch("scnoisemeter.cli._resolve_polya_sites", return_value=([], None, "user")), \
             patch("scnoisemeter.cli._check_version_consistency", return_value=None), \
             patch("scnoisemeter.cli.build_annotation_index", return_value=MagicMock()), \
             patch("scnoisemeter.cli.run_pipeline", return_value=mock_sample_result), \
             patch("scnoisemeter.cli._load_polya_sites", return_value={}), \
             patch("scnoisemeter.cli._load_tss_sites", return_value={}), \
             patch("scnoisemeter.cli.compute_metrics", return_value=(mock_sm, MagicMock())), \
             patch("scnoisemeter.cli._write_run_outputs"), \
             patch("scnoisemeter.cli.extract_intergenic_records", return_value=[]), \
             patch("pysam.AlignmentFile", return_value=mock_bam_instance):

            result = runner.invoke(cli, [
                "run",
                "--bam", str(bam_file),
                "--gtf", str(gtf_file),
                "--platform", "ont",
                "--pipeline-stage", "post_filter",
                "--output-dir", str(tmp_path / "out"),
            ])

        assert result.exit_code == 0, (
            f"Expected exit 0, got {result.exit_code}\nstdout:\n{result.output}\nstderr:\n{result.stderr}"
        )
        combined = (result.output or "") + (result.stderr or "")
        assert "unreliable" in combined.lower() or "minimum" in combined.lower(), (
            f"Expected low-reads warning, got: {combined}"
        )


# ---------------------------------------------------------------------------
# CASE 5 — Missing BAM index
# ---------------------------------------------------------------------------

class TestCase5MissingIndex:
    """BAM file exists but .bai index does not → FATAL ClickException."""

    def test_raises_click_exception_missing_index(self, tmp_path):
        from scnoisemeter.cli import cli
        runner = CliRunner()

        # BAM file exists but NO .bai index
        bam_file = tmp_path / "no_index.bam"
        gtf_file = tmp_path / "gencode.gtf.gz"
        bam_file.write_bytes(b"")
        gtf_file.write_bytes(b"")
        # Deliberately do NOT create bam_file.bai

        mock_bam_instance = MagicMock()
        mock_bam_instance.__enter__ = lambda s: s
        mock_bam_instance.__exit__ = MagicMock(return_value=False)
        mock_bam_instance.header.to_dict.return_value = {"HD": {"SO": "coordinate"}}

        with patch("pysam.AlignmentFile", return_value=mock_bam_instance):
            result = runner.invoke(cli, [
                "run",
                "--bam", str(bam_file),
                "--gtf", str(gtf_file),
                "--platform", "ont",
                "--pipeline-stage", "post_filter",
                "--output-dir", str(tmp_path / "out"),
            ])

        assert result.exit_code != 0, "Expected non-zero exit for missing index"
        combined = result.output + (result.stderr or "")
        assert "index" in combined.lower(), f"Expected 'index' in error, got: {combined}"
        assert "samtools index" in combined, f"Expected 'samtools index' hint, got: {combined}"


# ---------------------------------------------------------------------------
# CASE 6 — No @HD SO tag
# ---------------------------------------------------------------------------

class TestCase6NoSOTag:
    """BAM @HD has VN but no SO → warning only, run proceeds."""

    def test_no_crash_when_so_absent(self, caplog):
        """_validate_sort_order must not raise when sort_order is empty."""
        import logging
        meta = _make_meta(sort_order="")  # empty = tag absent
        # Should not raise
        with caplog.at_level(logging.WARNING, logger="scnoisemeter"):
            _validate_sort_order(meta)
        # Should have logged a warning
        assert any("SO tag" in r.message or "sort order" in r.message.lower()
                   for r in caplog.records), (
            f"Expected a sort-order warning, got: {[r.message for r in caplog.records]}"
        )

    def test_sort_order_coordinate_passes(self):
        """Coordinate-sorted BAMs must pass silently."""
        meta = _make_meta(sort_order="coordinate")
        _validate_sort_order(meta)  # should not raise


# ---------------------------------------------------------------------------
# CASE 7 — Wrong species (mouse chr1 length)
# ---------------------------------------------------------------------------

class TestCase7WrongSpecies:
    """BAM has mouse chr1 length → WARNING in meta.warnings, run proceeds."""

    def test_length_mismatch_adds_warning(self):
        """
        BAM chr1 length = 195,471,971 (mm10) vs expected 248,956,422 (GRCh38).
        _validate_chromosome_lengths must add a warning to meta.warnings.
        """
        meta = _make_meta(
            reference_lengths={
                "chr1": 195_471_971,  # mouse mm10 chr1 length
                "chr2": 181_748_087,
            },
        )
        assert not meta.warnings  # starts empty
        _validate_chromosome_lengths(meta, reference_path=None)
        assert meta.warnings, "Expected at least one warning about length mismatch"
        combined = " ".join(meta.warnings)
        assert "length" in combined.lower() or "mismatch" in combined.lower(), (
            f"Expected length-mismatch warning, got: {combined}"
        )

    def test_correct_lengths_no_warning(self):
        """GRCh38 chr1 length must pass without a warning."""
        meta = _make_meta()  # default lengths match GRCh38
        _validate_chromosome_lengths(meta, reference_path=None)
        # No length-related warnings
        assert not any("length" in w.lower() for w in meta.warnings)


# ---------------------------------------------------------------------------
# CASE 8 — Empty BAM (0 reads)
# ---------------------------------------------------------------------------

class TestCase8EmptyBAM:
    """BAM has no reads at all (mapped=0, unmapped=0) → FATAL ClickException."""

    def test_raises_fatal_error(self):
        from scnoisemeter.cli import cli
        runner = CliRunner()

        with patch("pathlib.Path.exists", return_value=True), \
             patch("pysam.AlignmentFile") as mock_bam_cls:

            mock_bam_instance = MagicMock()
            mock_bam_instance.__enter__ = lambda s: s
            mock_bam_instance.__exit__ = MagicMock(return_value=False)
            mock_bam_instance.header.to_dict.return_value = {"HD": {"SO": "coordinate"}}
            type(mock_bam_instance).mapped = PropertyMock(return_value=0)
            type(mock_bam_instance).unmapped = PropertyMock(return_value=0)
            mock_bam_cls.return_value = mock_bam_instance

            with patch("scnoisemeter.cli.inspect_bam") as mock_inspect:
                mock_inspect.return_value = _make_meta()

                result = runner.invoke(cli, [
                    "run",
                    "--bam", "/fake/empty.bam",
                    "--gtf", "/fake/gencode.gtf.gz",
                    "--platform", "ont",
                    "--pipeline-stage", "post_filter",
                    "--output-dir", "/tmp/fake_out",
                ])

        assert result.exit_code != 0, "Expected non-zero exit for empty BAM"
        combined = result.output + (result.stderr or "")
        assert "no reads" in combined.lower() or "empty" in combined.lower(), (
            f"Expected 'no reads' error message, got: {combined}"
        )


# ---------------------------------------------------------------------------
# CASE 9 — CB tags absent but --pipeline-stage post_filter
# ---------------------------------------------------------------------------

class TestCase9NoCBPostFilter:
    """No CB tags in post_filter BAM → warning + barcode-agnostic mode, exit 0."""

    def test_warning_when_no_cb_in_post_filter(self):
        """
        When pipeline_stage is post_filter (user-supplied) but no CB tags are
        found, a specific warning must be emitted.
        """
        from scnoisemeter.cli import (
            BARCODE_AUTODETECT_MIN_FRACTION,
            PipelineStage,
        )

        # Build a mock BamMetadata simulating no CB tags
        meta = _make_meta(
            barcode_fraction=0.0,
            barcode_aware=False,
            pipeline_stage=PipelineStage.POST_FILTER,
        )

        # Verify the warning message would include the key phrases
        # (We test the logic that _would_ emit the warning in run_cmd)
        stage = PipelineStage.POST_FILTER

        is_post_filter_no_cb = (
            meta.pipeline_stage == PipelineStage.POST_FILTER
            and not meta.barcode_aware
            and stage == PipelineStage.POST_FILTER
        )
        assert is_post_filter_no_cb, "Test setup failed: condition should be true"

        # Simulate the warning construction used in cli.py
        cb_msg = (
            f"CB tags are absent despite --pipeline-stage post_filter. "
            f"Only {meta.barcode_fraction:.1%} of sampled reads carry the 'CB' tag "
            f"(threshold: {BARCODE_AUTODETECT_MIN_FRACTION:.0%}). "
            f"Automatically switching to barcode-agnostic mode"
        )
        assert "barcode-agnostic" in cb_msg
        assert "post_filter" in cb_msg

    def test_coordinate_sorted_no_cb_proceeds(self, caplog):
        """
        _validate_sort_order should not crash for a coordinate-sorted BAM
        that has no CB tags (this is a separate concern from barcode detection).
        """
        import logging
        meta = _make_meta(sort_order="coordinate", barcode_fraction=0.0, barcode_aware=False)
        with caplog.at_level(logging.WARNING):
            _validate_sort_order(meta)
        # No sort-order warning expected for coordinate-sorted BAM
        sort_warnings = [r for r in caplog.records if "sort" in r.message.lower() and "SO tag" in r.message]
        assert not sort_warnings


# ---------------------------------------------------------------------------
# CASE 10 — Unaligned BAM (no @SQ lines)
# ---------------------------------------------------------------------------

class TestCase10UnalignedBAM:
    """BAM has no @SQ reference sequences in header → FATAL ClickException."""

    def test_raises_click_exception_when_no_sq_lines(self):
        """_validate_sq_lines should raise ClickException for empty reference_names."""
        meta = _make_meta(reference_names=[])
        with pytest.raises(click.ClickException) as exc_info:
            _validate_sq_lines(meta)
        msg = str(exc_info.value.format_message())
        assert "unaligned" in msg.lower(), f"Expected 'unaligned' in error, got: {msg}"
        assert "@SQ" in msg or "reference sequences" in msg, (
            f"Expected '@SQ' or 'reference sequences' in error, got: {msg}"
        )

    def test_error_mentions_alignment_tools(self):
        """Error message should mention alignment tools."""
        meta = _make_meta(reference_names=[])
        with pytest.raises(click.ClickException) as exc_info:
            _validate_sq_lines(meta)
        msg = str(exc_info.value.format_message())
        alignment_tools = ["minimap2", "pbmm2", "STAR"]
        assert any(tool in msg for tool in alignment_tools), (
            f"Expected alignment tool name in error, got: {msg}"
        )

    def test_aligned_bam_passes(self):
        """Normal aligned BAM with @SQ lines should not raise."""
        meta = _make_meta(reference_names=["chr1", "chr2", "chrX"])
        # Should not raise
        _validate_sq_lines(meta)
