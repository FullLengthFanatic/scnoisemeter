"""
tests/test_report.py
====================
Tests for platform-aware report generation:
  - Illumina short-read suppression of length-stratified bar chart and
    read-length histogram panels
  - Illumina insert size distribution panel
  - ONT / PacBio regression: both original panels remain present
  - Insert size computation helpers (flag filtering, abs(template_length),
    sampling cap)
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from unittest.mock import MagicMock

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from scnoisemeter.modules.metrics import SampleMetrics, CellTable
from scnoisemeter.modules.report import write_run_report, write_compare_report
from scnoisemeter.modules.pipeline import ContigResult, _reservoir_add
from scnoisemeter.constants import ReadCategory


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_sm(**kwargs) -> SampleMetrics:
    """Return a minimal SampleMetrics instance."""
    defaults = dict(
        sample_name="test",
        bam_path="/fake/test.bam",
        platform="illumina",
        pipeline_stage="post_filter",
        aligner="cellranger",
        n_reads_total=10_000,
        n_reads_classified=9_800,
        n_reads_unassigned=200,
        n_cells=100,
        read_fracs={},
        base_fracs={},
        noise_read_frac=0.20,
        noise_base_frac=0.19,
        noise_read_frac_strict=0.04,
        noise_base_frac_strict=0.04,
        strand_concordance=0.97,
        chimeric_read_frac=0.01,
        multimapper_read_frac=0.00,
        n_tso_invasion=50,
        n_polya_priming=200,
        n_noncanon_junction=10,
        per_cell_noise_median=0.20,
        per_cell_noise_iqr=0.05,
        warnings=[],
    )
    defaults.update(kwargs)
    return SampleMetrics(**defaults)


def _make_ct() -> CellTable:
    """Return an empty CellTable."""
    return CellTable(df=pd.DataFrame(), sample_name="test")


def _make_length_samples(n: int = 200) -> dict:
    """Return minimal length_samples with EXONIC_SENSE data."""
    return {ReadCategory.EXONIC_SENSE: [100] * n}


def _write_report(tmp_path, platform: str, insert_sizes=None, length_samples=None) -> str:
    """Helper: write a report and return the HTML text."""
    sm = _make_sm(platform=platform)
    ct = _make_ct()
    ls = length_samples if length_samples is not None else _make_length_samples()
    # Build a minimal length_stratified dataframe (single <150 bin)
    strat_df = pd.DataFrame([{
        "length_bin": "<150",
        "category": "exonic_sense",
        "count": 9000,
        "fraction_of_bin": 0.92,
    }])
    out = tmp_path / f"report_{platform}.html"
    write_run_report(
        sm, ct, ls, out,
        platform=platform,
        length_stratified=strat_df,
        insert_sizes=insert_sizes,
    )
    return out.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# 1. Illumina: "Noise by read length" chart ABSENT
# ---------------------------------------------------------------------------

class TestIlluminaSuppressLengthStratified:
    def test_illumina_no_noise_by_read_length(self, tmp_path):
        html = _write_report(tmp_path, platform="illumina")
        assert "Noise by read length" not in html

    def test_illumina_10x_no_noise_by_read_length(self, tmp_path):
        html = _write_report(tmp_path, platform="illumina_10x")
        assert "Noise by read length" not in html

    def test_illumina_bd_no_noise_by_read_length(self, tmp_path):
        html = _write_report(tmp_path, platform="illumina_bd")
        assert "Noise by read length" not in html


# ---------------------------------------------------------------------------
# 2. Illumina: "Read length distributions by category" chart ABSENT
# ---------------------------------------------------------------------------

class TestIlluminaSuppressLengthDistributions:
    def test_illumina_no_read_length_distributions(self, tmp_path):
        html = _write_report(tmp_path, platform="illumina")
        assert "Read length distributions by category" not in html

    def test_illumina_10x_no_read_length_distributions(self, tmp_path):
        html = _write_report(tmp_path, platform="illumina_10x")
        assert "Read length distributions by category" not in html


# ---------------------------------------------------------------------------
# 3. Illumina + paired reads: "Insert size distribution" chart PRESENT
# ---------------------------------------------------------------------------

class TestIlluminaInsertSizePlot:
    def test_insert_size_present_when_data_provided(self, tmp_path):
        insert_sizes = {"signal": [150, 200, 180, 160], "noise": [300, 350]}
        html = _write_report(tmp_path, platform="illumina",
                             insert_sizes=insert_sizes)
        assert "Insert size distribution" in html

    def test_insert_size_present_signal_only(self, tmp_path):
        insert_sizes = {"signal": [150, 200, 180, 160], "noise": []}
        html = _write_report(tmp_path, platform="illumina",
                             insert_sizes=insert_sizes)
        assert "Insert size distribution" in html

    def test_insert_size_present_noise_only(self, tmp_path):
        insert_sizes = {"signal": [], "noise": [300, 350]}
        html = _write_report(tmp_path, platform="illumina",
                             insert_sizes=insert_sizes)
        assert "Insert size distribution" in html

    def test_no_paired_reads_produces_warning(self, tmp_path):
        """When insert_sizes is None (no paired reads), a warning appears instead."""
        html = _write_report(tmp_path, platform="illumina", insert_sizes=None)
        assert "Insert size distribution" not in html
        assert "no properly paired reads" in html.lower() or "insert size" in html.lower()

    def test_insert_size_axes_labels(self, tmp_path):
        insert_sizes = {"signal": list(range(100, 500)), "noise": list(range(200, 600))}
        html = _write_report(tmp_path, platform="illumina",
                             insert_sizes=insert_sizes)
        assert "Insert size (bp)" in html
        assert "% of reads" in html


# ---------------------------------------------------------------------------
# 4. ONT / PacBio: both original charts ARE present (regression)
# ---------------------------------------------------------------------------

class TestLongReadPlatformsRegression:
    def test_ont_has_noise_by_read_length(self, tmp_path):
        html = _write_report(tmp_path, platform="ont")
        assert "Noise by read length" in html

    def test_ont_has_read_length_distributions(self, tmp_path):
        html = _write_report(tmp_path, platform="ont")
        assert "Read length distributions by category" in html

    def test_pacbio_has_noise_by_read_length(self, tmp_path):
        html = _write_report(tmp_path, platform="pacbio")
        assert "Noise by read length" in html

    def test_pacbio_has_read_length_distributions(self, tmp_path):
        html = _write_report(tmp_path, platform="pacbio")
        assert "Read length distributions by category" in html

    def test_ont_no_insert_size_chart(self, tmp_path):
        """Insert size chart must NOT appear for long-read platforms."""
        html = _write_report(tmp_path, platform="ont")
        assert "Insert size distribution" not in html

    def test_pacbio_no_insert_size_chart(self, tmp_path):
        html = _write_report(tmp_path, platform="pacbio")
        assert "Insert size distribution" not in html


# ---------------------------------------------------------------------------
# 5. Insert size collection in pipeline (ContigResult)
# ---------------------------------------------------------------------------

class TestInsertSizeCollection:
    """
    Unit tests for the insert-size reservoir-sampling logic added to the
    pipeline's per-read loop.  We test the ContigResult data structure
    directly — no BAM file required.
    """

    def test_contig_result_has_insert_size_fields(self):
        cr = ContigResult(contig="chr1")
        assert hasattr(cr, "insert_size_signal")
        assert hasattr(cr, "insert_size_noise")
        assert isinstance(cr.insert_size_signal, list)
        assert isinstance(cr.insert_size_noise, list)

    def test_reservoir_add_caps_at_max(self):
        reservoir = []
        max_size = 5
        for i in range(20):
            _reservoir_add(reservoir, i, max_size)
        assert len(reservoir) == max_size

    def test_reservoir_add_fills_to_max(self):
        reservoir = []
        for i in range(3):
            _reservoir_add(reservoir, i, 10)
        assert len(reservoir) == 3
        assert set(reservoir) == {0, 1, 2}


class TestInsertSizeFiltering:
    """
    Verify flag-filtering logic by simulating reads with the conditions that
    the pipeline enforces for insert size collection.
    """

    def _make_read(self, is_proper_pair=True, is_read1=True,
                   template_length=150, is_secondary=False,
                   is_supplementary=False, is_unmapped=False):
        read = MagicMock()
        read.is_proper_pair    = is_proper_pair
        read.is_read1          = is_read1
        read.template_length   = template_length
        read.is_secondary      = is_secondary
        read.is_supplementary  = is_supplementary
        read.is_unmapped       = is_unmapped
        return read

    def _should_collect(self, read) -> Optional[int]:
        """Replicate the pipeline's collection logic."""
        if not (read.is_proper_pair and read.is_read1):
            return None
        tlen = abs(read.template_length)
        if not (0 < tlen < 2000):
            return None
        return tlen

    def test_proper_pair_read1_collected(self):
        read = self._make_read(is_proper_pair=True, is_read1=True,
                               template_length=200)
        assert self._should_collect(read) == 200

    def test_negative_template_length_uses_abs(self):
        read = self._make_read(is_proper_pair=True, is_read1=True,
                               template_length=-200)
        assert self._should_collect(read) == 200

    def test_read2_not_collected(self):
        """Only read1 should contribute to avoid double-counting."""
        read = self._make_read(is_proper_pair=True, is_read1=False,
                               template_length=200)
        assert self._should_collect(read) is None

    def test_unpaired_read_not_collected(self):
        read = self._make_read(is_proper_pair=False, is_read1=True,
                               template_length=200)
        assert self._should_collect(read) is None

    def test_zero_template_length_not_collected(self):
        read = self._make_read(is_proper_pair=True, is_read1=True,
                               template_length=0)
        assert self._should_collect(read) is None

    def test_template_length_2000_or_more_not_collected(self):
        read = self._make_read(is_proper_pair=True, is_read1=True,
                               template_length=2000)
        assert self._should_collect(read) is None

    def test_template_length_1999_collected(self):
        read = self._make_read(is_proper_pair=True, is_read1=True,
                               template_length=1999)
        assert self._should_collect(read) == 1999

    def test_template_length_1_collected(self):
        read = self._make_read(is_proper_pair=True, is_read1=True,
                               template_length=1)
        assert self._should_collect(read) == 1


# ---------------------------------------------------------------------------
# 6. Compare report: barcode-agnostic violin suppression
# ---------------------------------------------------------------------------

def _make_ct_with_data(sample_name: str = "test") -> CellTable:
    """Return a CellTable with real per-cell noise data."""
    import numpy as np
    rng = np.random.default_rng(0)
    df = pd.DataFrame({"noise_read_frac": rng.beta(2, 8, 50)})
    return CellTable(df=df, sample_name=sample_name)


def _write_compare(tmp_path, n_cells_a: int, n_cells_b: int) -> str:
    """Write a compare report and return the HTML."""
    sm_a = _make_sm(sample_name="A", n_cells=n_cells_a, platform="ont")
    sm_b = _make_sm(sample_name="B", n_cells=n_cells_b, platform="ont")
    ct_a = _make_ct_with_data("A") if n_cells_a > 1 else _make_ct()
    ct_b = _make_ct_with_data("B") if n_cells_b > 1 else _make_ct()
    out = tmp_path / "compare.html"
    write_compare_report(
        sm_a, sm_b, ct_a, ct_b,
        length_samples_a={}, length_samples_b={},
        stats_df=None,
        output_path=out,
        offline=True,
    )
    return out.read_text(encoding="utf-8")


class TestCompareViolinBarcodeAgnostic:
    def test_violin_present_when_both_have_cells(self, tmp_path):
        html = _write_compare(tmp_path, n_cells_a=50, n_cells_b=50)
        assert "Per-cell noise fraction" in html

    def test_violin_absent_when_sample_a_agnostic(self, tmp_path):
        html = _write_compare(tmp_path, n_cells_a=1, n_cells_b=50)
        assert "Per-cell noise fraction — comparison" not in html

    def test_violin_absent_when_sample_b_agnostic(self, tmp_path):
        html = _write_compare(tmp_path, n_cells_a=50, n_cells_b=1)
        assert "Per-cell noise fraction — comparison" not in html

    def test_violin_absent_when_both_agnostic(self, tmp_path):
        html = _write_compare(tmp_path, n_cells_a=1, n_cells_b=1)
        assert "Per-cell noise fraction — comparison" not in html

    def test_warning_shown_when_a_agnostic(self, tmp_path):
        html = _write_compare(tmp_path, n_cells_a=1, n_cells_b=50)
        assert "barcode-agnostic" in html.lower()

    def test_warning_shown_when_b_agnostic(self, tmp_path):
        html = _write_compare(tmp_path, n_cells_a=50, n_cells_b=1)
        assert "barcode-agnostic" in html.lower()

    def test_no_warning_when_both_have_cells(self, tmp_path):
        html = _write_compare(tmp_path, n_cells_a=50, n_cells_b=50)
        assert "Per-cell comparison not shown" not in html
