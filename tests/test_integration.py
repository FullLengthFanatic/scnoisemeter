"""
tests/test_integration.py
=========================
End-to-end smoke tests using the tiny synthetic BAM and GTF from conftest.py.

These tests exercise paths that the mocked unit tests cannot reach:
  - pyranges annotation index building from a real GTF
  - pysam BAM reading in the chromosome workers
  - interaction between AnnotationIndex and ReadClassifier
  - metric computation on a real SampleResult

All fixtures are session-scoped so the annotation index and BAM are built
once per test session, not once per test.

See tests/conftest.py for the GTF layout and BAM read inventory.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pytest

from scnoisemeter.constants import ReadCategory
from scnoisemeter.modules.annotation import build_annotation_index
from scnoisemeter.modules.metrics import compute_metrics
from scnoisemeter.modules.pipeline import run_pipeline
from scnoisemeter.utils.bam_inspector import inspect_bam

from tests.conftest import BARCODE_CELL1, BARCODE_CELL2


# ---------------------------------------------------------------------------
# Session-scoped integration fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def annotation_index(tiny_gtf):
    """Build the annotation index from the tiny GTF (no cache, no side effects)."""
    return build_annotation_index(tiny_gtf, cache=False)


@pytest.fixture(scope="session")
def bam_meta(tiny_bam):
    """Inspect the tiny BAM and return its BamMetadata."""
    return inspect_bam(tiny_bam)


@pytest.fixture(scope="session")
def pipeline_result(tiny_bam, annotation_index, bam_meta):
    """
    Run the full pipeline on the tiny BAM.

    threads=1 keeps the test single-process.
    contigs=["chr1"] is passed explicitly because _select_contigs uses a
    1 Mb minimum — chr1 at GRCh38 length passes, but being explicit avoids
    any version-dependent filtering changes breaking the fixture.
    chrM is added automatically by run_pipeline via the mito-contig path.
    """
    return run_pipeline(
        tiny_bam,
        annotation_index,
        bam_meta,
        threads=1,
        store_umis=True,
        contigs=["chr1"],
    )


@pytest.fixture(scope="session")
def sample_metrics(pipeline_result):
    metrics, _ = compute_metrics(pipeline_result, "tiny_test")
    return metrics


@pytest.fixture(scope="session")
def cell_table(pipeline_result):
    _, ct = compute_metrics(pipeline_result, "tiny_test")
    return ct


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _total_reads_for_category(result, category: ReadCategory) -> int:
    return sum(
        cb_counts.get(category, 0)
        for cb_counts in result.read_counts.values()
    )


# ---------------------------------------------------------------------------
# Annotation index
# ---------------------------------------------------------------------------

class TestAnnotationIndex:
    def test_exons_plus_nonempty(self, annotation_index):
        assert len(annotation_index.exons_plus.df) > 0

    def test_introns_plus_nonempty(self, annotation_index):
        assert len(annotation_index.introns_plus.df) > 0

    def test_exons_minus_nonempty(self, annotation_index):
        """GENE_B is on the minus strand — must appear in exons_minus."""
        assert len(annotation_index.exons_minus.df) > 0

    def test_intergenic_nonempty(self, annotation_index):
        assert len(annotation_index.intergenic.df) > 0

    def test_gene_a_exon_coordinates(self, annotation_index):
        """GENE_A exon 1 (GTF 1000-2000) should be Start=999, End=2000 (0-based)."""
        df = annotation_index.exons_plus.df
        gene_a = df[df["gene_id"] == "ENSG00000000001"]
        assert not gene_a.empty
        starts = set(gene_a["Start"].tolist())
        assert 999 in starts, f"Expected exon1 Start=999, got {starts}"

    def test_intron_between_exon1_and_exon2(self, annotation_index):
        """GENE_A intron 1 should span 2000–2999 (0-based)."""
        df = annotation_index.introns_plus.df
        gene_a = df[df["gene_id"] == "ENSG00000000001"]
        assert not gene_a.empty
        intron_starts = set(gene_a["Start"].tolist())
        assert 2000 in intron_starts, (
            f"Expected intron1 Start=2000, got {intron_starts}"
        )

    def test_gene_b_minus_strand(self, annotation_index):
        """GENE_B should appear only in exons_minus, not exons_plus."""
        plus_ids  = set(annotation_index.exons_plus.df["gene_id"].tolist())
        minus_ids = set(annotation_index.exons_minus.df["gene_id"].tolist())
        assert "ENSG00000000002" in minus_ids
        assert "ENSG00000000002" not in plus_ids

    def test_no_mito_in_index(self, annotation_index):
        """Mitochondrial contigs must be excluded from the annotation index."""
        for attr in ("exons_plus", "exons_minus", "introns_plus", "introns_minus"):
            df = getattr(annotation_index, attr).df
            assert "chrM" not in df["Chromosome"].values, (
                f"chrM found in {attr} — mito filtering broken"
            )

    def test_intergenic_end_below_sentinel(self, annotation_index):
        """
        No intergenic interval should extend to the 2-billion-bp sentinel.
        The previous _manual_complement() added End=2_000_000_000 after the
        last gene on each chromosome, inflating the Poisson denominator ~24×.
        """
        df = annotation_index.intergenic.df
        assert (df["End"] < 2_000_000_000).all(), (
            "Intergenic interval(s) still reach the 2B sentinel — "
            "_manual_complement fix not applied"
        )


# ---------------------------------------------------------------------------
# BAM inspection
# ---------------------------------------------------------------------------

class TestBamInspection:
    def test_platform_detected_as_ont(self, bam_meta):
        """minimap2 @PG record should resolve to ONT platform."""
        from scnoisemeter.constants import Platform
        assert bam_meta.platform == Platform.ONT

    def test_barcode_aware(self, bam_meta):
        """All but possibly one read carry CB tags → barcode_aware=True."""
        assert bam_meta.barcode_aware is True

    def test_reference_names_include_chr1_and_chrm(self, bam_meta):
        assert "chr1" in bam_meta.reference_names
        assert "chrM" in bam_meta.reference_names

    def test_chr1_length_is_grch38(self, bam_meta):
        assert bam_meta.reference_lengths.get("chr1") == 248_956_422


# ---------------------------------------------------------------------------
# Pipeline execution
# ---------------------------------------------------------------------------

class TestPipelineExecution:
    def test_reads_processed_positive(self, pipeline_result):
        assert pipeline_result.n_reads_processed > 0

    def test_both_barcodes_present(self, pipeline_result):
        assert BARCODE_CELL1 in pipeline_result.read_counts
        assert BARCODE_CELL2 in pipeline_result.read_counts

    def test_exonic_sense_classified(self, pipeline_result):
        n = _total_reads_for_category(pipeline_result, ReadCategory.EXONIC_SENSE)
        assert n >= 1, "Expected at least one EXONIC_SENSE read"

    def test_exonic_antisense_classified(self, pipeline_result):
        n = _total_reads_for_category(pipeline_result, ReadCategory.EXONIC_ANTISENSE)
        assert n >= 1, "Expected at least one EXONIC_ANTISENSE read"

    def test_multimapper_classified(self, pipeline_result):
        n = _total_reads_for_category(pipeline_result, ReadCategory.MULTIMAPPER)
        assert n >= 1, "Expected at least one MULTIMAPPER read (NH=3 read in test BAM)"

    def test_chimeric_classified(self, pipeline_result):
        n = _total_reads_for_category(pipeline_result, ReadCategory.CHIMERIC)
        assert n >= 1, "Expected at least one CHIMERIC read (SA tag, 48 kbp gap)"

    def test_intronic_classified(self, pipeline_result):
        pure = _total_reads_for_category(pipeline_result, ReadCategory.INTRONIC_PURE)
        boundary = _total_reads_for_category(pipeline_result, ReadCategory.INTRONIC_BOUNDARY)
        assert pure + boundary >= 1, (
            "Expected at least one intronic read (pure or boundary)"
        )

    def test_intergenic_classified(self, pipeline_result):
        cats = (
            ReadCategory.INTERGENIC_SPARSE,
            ReadCategory.INTERGENIC_HOTSPOT,
            ReadCategory.INTERGENIC_NOVEL,
            ReadCategory.INTERGENIC_REPEAT,
        )
        n = sum(_total_reads_for_category(pipeline_result, c) for c in cats)
        assert n >= 1, "Expected at least one intergenic read"

    def test_mitochondrial_classified(self, pipeline_result):
        n = _total_reads_for_category(pipeline_result, ReadCategory.MITOCHONDRIAL)
        assert n >= 1, "Expected at least one MITOCHONDRIAL read"

    def test_no_spurious_supplementary_counted(self, pipeline_result):
        """SUPPLEMENTARY reads must not appear in read_counts."""
        n = _total_reads_for_category(pipeline_result, ReadCategory.SUPPLEMENTARY)
        assert n == 0

    def test_base_counts_positive(self, pipeline_result):
        total_bases = sum(
            sum(cat_bases.values())
            for cat_bases in pipeline_result.base_counts.values()
        )
        assert total_bases > 0

    def test_umi_sets_populated(self, pipeline_result):
        """UMI sets must have at least one entry per classified cell."""
        assert len(pipeline_result.umi_sets) > 0


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

class TestMetrics:
    def test_n_reads_classified_positive(self, sample_metrics):
        assert sample_metrics.n_reads_classified > 0

    def test_noise_frac_in_unit_interval(self, sample_metrics):
        assert 0.0 <= sample_metrics.noise_read_frac <= 1.0
        assert 0.0 <= sample_metrics.noise_base_frac <= 1.0

    def test_noise_frac_positive(self, sample_metrics):
        """We injected chimeric and intronic reads so noise must be > 0."""
        assert sample_metrics.noise_read_frac > 0.0

    def test_strict_noise_leq_conservative(self, sample_metrics):
        """Strict noise is a subset of conservative noise."""
        assert sample_metrics.noise_read_frac_strict <= sample_metrics.noise_read_frac

    def test_strand_concordance_in_unit_interval(self, sample_metrics):
        assert 0.0 <= sample_metrics.strand_concordance <= 1.0

    def test_chimeric_frac_positive(self, sample_metrics):
        assert sample_metrics.chimeric_read_frac > 0.0

    def test_multimapper_frac_positive(self, sample_metrics):
        assert sample_metrics.multimapper_read_frac > 0.0

    def test_read_fracs_sum_to_one(self, sample_metrics):
        total = sum(sample_metrics.read_fracs.values())
        assert abs(total - 1.0) < 1e-6, f"read_fracs sum to {total}, expected 1.0"

    def test_n_cells_positive(self, sample_metrics):
        assert sample_metrics.n_cells >= 1


class TestCellTable:
    def test_two_cells_in_table(self, cell_table):
        assert len(cell_table.df) == 2

    def test_noise_column_present(self, cell_table):
        assert "noise_read_frac" in cell_table.df.columns

    def test_all_category_columns_present(self, cell_table):
        from scnoisemeter.constants import CATEGORY_ORDER
        for cat in CATEGORY_ORDER:
            col = f"read_frac_{cat.value}"
            assert col in cell_table.df.columns, f"Missing column: {col}"

    def test_per_cell_fracs_sum_to_one(self, cell_table):
        from scnoisemeter.constants import CATEGORY_ORDER
        frac_cols = [f"read_frac_{c.value}" for c in CATEGORY_ORDER]
        for bc, row in cell_table.df.iterrows():
            total = row[frac_cols].sum()
            assert abs(total - 1.0) < 1e-6, (
                f"read fracs for barcode {bc} sum to {total}, expected 1.0"
            )
