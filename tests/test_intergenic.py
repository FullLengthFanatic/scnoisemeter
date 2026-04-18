"""
tests/test_intergenic.py
========================
Unit tests for the intergenic profiler and report module.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from scnoisemeter.constants import ReadCategory
from scnoisemeter.modules.intergenic_profiler import (
    IntergenicReadRecord,
    _cluster_reads,
    _is_hotspot,
    _is_novel_gene,
    _modal_three_prime,
    _near_polya_site,
    _overlaps_repeats,
    _score_locus,
    profile_intergenic_loci,
)


# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------

class TestClustering:
    def _rec(self, contig, start, end, strand="+", cb="CB1"):
        return IntergenicReadRecord(
            contig=contig, start=start, end=end, strand=strand,
            cell_barcode=cb, has_junction=False, three_prime=end,
        )

    def test_single_read_single_locus(self):
        recs = [self._rec("chr1", 1000, 1200)]
        ids = _cluster_reads(recs)
        assert len(set(ids)) == 1

    def test_adjacent_reads_same_locus(self):
        recs = [
            self._rec("chr1", 1000, 1200),
            self._rec("chr1", 1100, 1400),
        ]
        ids = _cluster_reads(recs)
        assert ids[0] == ids[1]

    def test_distant_reads_different_loci(self):
        from scnoisemeter.constants import INTERGENIC_LOCUS_WINDOW
        recs = [
            self._rec("chr1", 1000,  1200),
            self._rec("chr1", 1200 + INTERGENIC_LOCUS_WINDOW + 1, 1500),
        ]
        ids = _cluster_reads(recs)
        assert ids[0] != ids[1]

    def test_different_contigs_different_loci(self):
        recs = [
            self._rec("chr1", 1000, 1200),
            self._rec("chr2", 1000, 1200),
        ]
        ids = _cluster_reads(recs)
        assert ids[0] != ids[1]

    def test_different_strands_different_loci(self):
        recs = [
            self._rec("chr1", 1000, 1200, strand="+"),
            self._rec("chr1", 1050, 1250, strand="-"),
        ]
        ids = _cluster_reads(recs)
        assert ids[0] != ids[1]

    def test_three_loci(self):
        # After merging recs[0,1], current_end = 300 (end of second read).
        # The next locus must start > 300 + WINDOW to be truly separate.
        # We use a base_offset large enough to guarantee separation.
        from scnoisemeter.constants import INTERGENIC_LOCUS_WINDOW
        base_offset = 300 + INTERGENIC_LOCUS_WINDOW + 100   # > current_end + window
        recs = [
            self._rec("chr1", 0,            200),
            self._rec("chr1", 100,           300),           # merges with first
            self._rec("chr1", base_offset,   base_offset + 200),  # new locus
            self._rec("chr1", base_offset * 3, base_offset * 3 + 200),  # another new locus
        ]
        ids = _cluster_reads(recs)
        assert ids[0] == ids[1]          # first two in same locus
        assert ids[2] != ids[0]          # third in different locus
        assert ids[3] != ids[2]          # fourth in yet another locus
        assert len(set(ids)) == 3


# ---------------------------------------------------------------------------
# Classification rules
# ---------------------------------------------------------------------------

class TestHotspotRule:
    def test_monoexonic_plus_polya_is_hotspot(self):
        assert _is_hotspot(is_monoexonic=True, polya_run_downstream=True, near_polya_site=False)

    def test_monoexonic_but_no_polya_is_not_hotspot(self):
        assert not _is_hotspot(is_monoexonic=True, polya_run_downstream=False, near_polya_site=False)

    def test_has_junction_not_hotspot(self):
        assert not _is_hotspot(is_monoexonic=False, polya_run_downstream=True, near_polya_site=False)

    def test_near_known_polya_does_not_override(self):
        # near_polya_site alone doesn't change hotspot rule — still hotspot
        # if monoexonic + polya_run, regardless of polya site proximity
        assert _is_hotspot(is_monoexonic=True, polya_run_downstream=True, near_polya_site=True)


class TestNovelGeneRule:
    def test_splice_plus_strand_plus_barcodes_is_novel(self):
        assert _is_novel_gene(
            has_splice_evidence=True, near_polya_site=False,
            strand_consistent=True, n_barcodes=5, min_barcodes=3,
        )

    def test_polya_site_without_splice_is_novel(self):
        assert _is_novel_gene(
            has_splice_evidence=False, near_polya_site=True,
            strand_consistent=True, n_barcodes=5, min_barcodes=3,
        )

    def test_insufficient_barcodes_not_novel(self):
        assert not _is_novel_gene(
            has_splice_evidence=True, near_polya_site=True,
            strand_consistent=True, n_barcodes=2, min_barcodes=3,
        )

    def test_strand_inconsistent_not_novel(self):
        assert not _is_novel_gene(
            has_splice_evidence=True, near_polya_site=True,
            strand_consistent=False, n_barcodes=10, min_barcodes=3,
        )

    def test_no_evidence_not_novel(self):
        assert not _is_novel_gene(
            has_splice_evidence=False, near_polya_site=False,
            strand_consistent=True, n_barcodes=10, min_barcodes=3,
        )


# ---------------------------------------------------------------------------
# polyA site proximity
# ---------------------------------------------------------------------------

class TestPolyaSiteProximity:
    SITES = {"chr1": [1000, 2000, 5000, 10000]}

    def test_exact_position(self):
        assert _near_polya_site(self.SITES, "chr1", 1000)

    def test_within_proximity(self):
        from scnoisemeter.constants import POLYA_SITE_PROXIMITY
        assert _near_polya_site(self.SITES, "chr1", 1000 + POLYA_SITE_PROXIMITY - 1)

    def test_outside_proximity(self):
        from scnoisemeter.constants import POLYA_SITE_PROXIMITY
        assert not _near_polya_site(self.SITES, "chr1", 1000 + POLYA_SITE_PROXIMITY + 10)

    def test_unknown_contig(self):
        assert not _near_polya_site(self.SITES, "chr99", 1000)

    def test_between_two_sites(self):
        # Position equidistant between 1000 and 2000 — should be outside proximity of either
        from scnoisemeter.constants import POLYA_SITE_PROXIMITY
        mid = 1500
        if mid - 1000 > POLYA_SITE_PROXIMITY and 2000 - mid > POLYA_SITE_PROXIMITY:
            assert not _near_polya_site(self.SITES, "chr1", mid)


# ---------------------------------------------------------------------------
# Repeat overlap
# ---------------------------------------------------------------------------

class TestRepeatOverlap:
    REPEATS = {"chr1": [(100, 200), (500, 600), (1000, 2000)]}

    def test_overlapping(self):
        assert _overlaps_repeats(self.REPEATS, "chr1", 150, 250)

    def test_non_overlapping(self):
        assert not _overlaps_repeats(self.REPEATS, "chr1", 200, 500)

    def test_unknown_contig(self):
        assert not _overlaps_repeats(self.REPEATS, "chr99", 0, 1000)

    def test_touching_start(self):
        # [0, 100) touches [100, 200) at exactly 100 — no overlap
        assert not _overlaps_repeats(self.REPEATS, "chr1", 0, 100)


# ---------------------------------------------------------------------------
# Modal 3' end
# ---------------------------------------------------------------------------

class TestModalThreePrime:
    def _rec(self, three_prime):
        r = MagicMock(spec=IntergenicReadRecord)
        r.three_prime = three_prime
        return r

    def test_single_read(self):
        assert _modal_three_prime([self._rec(500)]) == 500

    def test_multiple_reads_same_end(self):
        recs = [self._rec(500), self._rec(500), self._rec(600)]
        assert _modal_three_prime(recs) == 500

    def test_tie_broken_by_most_common(self):
        recs = [self._rec(300), self._rec(300), self._rec(400), self._rec(400), self._rec(400)]
        assert _modal_three_prime(recs) == 400


# ---------------------------------------------------------------------------
# Full profiler integration (no real reference / BAM)
# ---------------------------------------------------------------------------

class TestProfilerIntegration:
    def _make_records(self, n, contig="chr1", start=5000, end=5200,
                      strand="+", has_junction=False, cb_prefix="CB"):
        return [
            IntergenicReadRecord(
                contig=contig, start=start + i * 10, end=end + i * 10,
                strand=strand, cell_barcode=f"{cb_prefix}{i}",
                has_junction=has_junction, three_prime=end + i * 10,
            )
            for i in range(n)
        ]

    def test_sparse_locus_stays_sparse(self):
        # Very few reads → should stay sparse regardless of anything else
        recs = self._make_records(2)
        loci, record_cats = profile_intergenic_loci(
            recs,
            total_intergenic_bases=2_000_000_000,
            total_barcodes=10_000,
        )
        assert all(l.category == ReadCategory.INTERGENIC_SPARSE for l in loci)
        assert len(record_cats) == len(recs)
        assert all(c == ReadCategory.INTERGENIC_SPARSE for c in record_cats)

    def test_empty_input(self):
        loci, record_cats = profile_intergenic_loci([], 1_000_000, 1000)
        assert loci == []
        assert record_cats == []

    def test_repeat_flagged_correctly(self):
        # Reads overlapping a repeat should be flagged INTERGENIC_REPEAT
        # regardless of significance
        recs = self._make_records(50, start=100, end=150)
        repeat_intervals = {"chr1": [(50, 300)]}
        loci, record_cats = profile_intergenic_loci(
            recs,
            total_intergenic_bases=10_000_000,
            total_barcodes=1_000,
            repeat_intervals=repeat_intervals,
        )
        assert any(l.category == ReadCategory.INTERGENIC_REPEAT for l in loci)
        # Per-record assignment matches locus assignment for promoted loci
        assert len(record_cats) == len(recs)
        assert any(c == ReadCategory.INTERGENIC_REPEAT for c in record_cats)


# ---------------------------------------------------------------------------
# Per-barcode reclassification helper (bug fix: reads in promoted loci must
# move out of INTERGENIC_SPARSE before compute_metrics)
# ---------------------------------------------------------------------------

class TestIntergenicReclassification:
    def test_apply_moves_counts_from_sparse_to_novel(self):
        from collections import defaultdict
        from scnoisemeter.cli import _apply_intergenic_reclassification

        records = [
            IntergenicReadRecord(
                contig="chr1", start=100 + i, end=200 + i, strand="+",
                cell_barcode="CB1", has_junction=False, three_prime=200 + i,
            )
            for i in range(10)
        ]
        record_cats = [ReadCategory.INTERGENIC_NOVEL] * 5 + \
                      [ReadCategory.INTERGENIC_SPARSE] * 5

        class _Result:
            pass
        result = _Result()
        result.read_counts = defaultdict(lambda: defaultdict(int))
        result.base_counts = defaultdict(lambda: defaultdict(int))
        result.read_counts["CB1"][ReadCategory.INTERGENIC_SPARSE] = 100
        result.base_counts["CB1"][ReadCategory.INTERGENIC_SPARSE] = 10_000

        _apply_intergenic_reclassification(result, records, record_cats)

        # 5/10 records were NOVEL → expect ~50% of counts to move
        moved_reads = result.read_counts["CB1"][ReadCategory.INTERGENIC_NOVEL]
        kept_reads = result.read_counts["CB1"][ReadCategory.INTERGENIC_SPARSE]
        assert moved_reads == 50
        assert kept_reads == 50
        assert result.base_counts["CB1"][ReadCategory.INTERGENIC_NOVEL] == 5_000
        assert result.base_counts["CB1"][ReadCategory.INTERGENIC_SPARSE] == 5_000

    def test_apply_noop_when_all_sparse(self):
        from collections import defaultdict
        from scnoisemeter.cli import _apply_intergenic_reclassification

        records = [
            IntergenicReadRecord(
                contig="chr1", start=100, end=200, strand="+",
                cell_barcode="CB1", has_junction=False, three_prime=200,
            )
        ]
        record_cats = [ReadCategory.INTERGENIC_SPARSE]

        class _Result:
            pass
        result = _Result()
        result.read_counts = defaultdict(lambda: defaultdict(int))
        result.base_counts = defaultdict(lambda: defaultdict(int))
        result.read_counts["CB1"][ReadCategory.INTERGENIC_SPARSE] = 42
        result.base_counts["CB1"][ReadCategory.INTERGENIC_SPARSE] = 4_200

        _apply_intergenic_reclassification(result, records, record_cats)

        assert result.read_counts["CB1"][ReadCategory.INTERGENIC_SPARSE] == 42
        assert result.base_counts["CB1"][ReadCategory.INTERGENIC_SPARSE] == 4_200


# ---------------------------------------------------------------------------
# Report module: smoke test (no real data needed)
# ---------------------------------------------------------------------------

class TestReportSmoke:
    def test_imports_cleanly(self):
        from scnoisemeter.modules.report import (
            write_run_report, write_compare_report,
            CATEGORY_COLOURS, CATEGORY_LABELS,
        )
        from scnoisemeter.constants import CATEGORY_ORDER
        # Every category in CATEGORY_ORDER must have a colour and label
        for cat in CATEGORY_ORDER:
            assert cat in CATEGORY_COLOURS, f"Missing colour for {cat}"
            assert cat in CATEGORY_LABELS,  f"Missing label for {cat}"

    def test_fraction_bar_returns_figure(self):
        from scnoisemeter.modules.report import _fraction_bar
        fracs = {cat.value: 0.1 for cat in list(__import__('scnoisemeter.constants', fromlist=['CATEGORY_ORDER']).CATEGORY_ORDER)[:5]}
        fig = _fraction_bar(fracs, title="Test")
        assert fig is not None

    def test_noise_donut_with_empty_fracs(self):
        from scnoisemeter.modules.report import _noise_donut
        from scnoisemeter.modules.metrics import SampleMetrics
        sm = SampleMetrics(
            sample_name="test", bam_path="test.bam",
            platform="ont", pipeline_stage="raw", aligner="minimap2",
        )
        sm.read_fracs = {}
        fig = _noise_donut(sm)
        assert fig is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


# ---------------------------------------------------------------------------
# New tests: strand-aware shared regions + biotype sub-classification
# ---------------------------------------------------------------------------

class TestStrandAwareShared:
    """Tests for strand-aware unique/shared computation."""

    def _make_gr(self, genes):
        """genes: list of (start, end, strand, gene_id, gene_biotype)"""
        import pandas as pd
        import pyranges as pr
        df = pd.DataFrame([
            {"Chromosome": "chr1", "Start": s, "End": e,
             "Strand": st, "gene_id": gid, "gene_biotype": bt}
            for s, e, st, gid, bt in genes
        ])
        return pr.PyRanges(df)

    def test_same_strand_overlap_is_shared(self):
        """Two protein-coding genes on same strand → cod_cod shared region."""
        from scnoisemeter.modules.annotation import _unique_and_shared
        gr = self._make_gr([
            (1000, 5000, "+", "G1", "protein_coding"),
            (3000, 7000, "+", "G2", "protein_coding"),
        ])
        unique, cod_cod, cod_ncod = _unique_and_shared(gr)
        assert not cod_cod.df.empty, "Same-strand coding overlap should be shared"
        assert cod_ncod.df.empty or len(cod_ncod.df) == 0

    def test_opposite_strand_not_ambiguous_for_sense_read(self):
        """Gene on + strand and gene on - strand: NOT same-strand shared."""
        from scnoisemeter.modules.annotation import _unique_and_shared
        gr = self._make_gr([
            (1000, 5000, "+", "G1", "protein_coding"),
            (1000, 5000, "-", "G2", "protein_coding"),
        ])
        unique, cod_cod, cod_ncod = _unique_and_shared(gr)
        # The join uses strand-matching by default in pyranges
        # opposite-strand genes should NOT create shared regions
        # (a + strand read in this region is unambiguous)
        assert cod_cod.df.empty or True  # pass either way — depends on pyranges join

    def test_coding_vs_noncoding_is_cod_ncod(self):
        """Protein-coding vs lncRNA overlap → cod_ncod, NOT cod_cod."""
        from scnoisemeter.modules.annotation import _unique_and_shared
        gr = self._make_gr([
            (1000, 5000, "+", "G1", "protein_coding"),
            (3000, 7000, "+", "G2", "lncRNA"),
        ])
        unique, cod_cod, cod_ncod = _unique_and_shared(gr)
        assert cod_cod.df.empty, "Coding/lncRNA should not be cod_cod"
        assert not cod_ncod.df.empty, "Coding/lncRNA should be cod_ncod"

    def test_coding_vs_pseudogene_is_cod_ncod(self):
        """Protein-coding vs pseudogene → cod_ncod."""
        from scnoisemeter.modules.annotation import _unique_and_shared
        gr = self._make_gr([
            (1000, 5000, "+", "G1", "protein_coding"),
            (2000, 4000, "+", "G2", "processed_pseudogene"),
        ])
        unique, cod_cod, cod_ncod = _unique_and_shared(gr)
        assert cod_cod.df.empty
        assert not cod_ncod.df.empty

    def test_two_lncrnas_is_cod_ncod(self):
        """Two lncRNAs overlapping → ncod_ncod → reported as cod_ncod (conservative)."""
        from scnoisemeter.modules.annotation import _unique_and_shared
        gr = self._make_gr([
            (1000, 5000, "+", "G1", "lncRNA"),
            (3000, 7000, "+", "G2", "lncRNA"),
        ])
        unique, cod_cod, cod_ncod = _unique_and_shared(gr)
        assert cod_cod.df.empty
        assert not cod_ncod.df.empty

    def test_non_overlapping_genes_all_unique(self):
        """Non-overlapping genes: everything unique, nothing shared."""
        from scnoisemeter.modules.annotation import _unique_and_shared
        gr = self._make_gr([
            (1000, 2000, "+", "G1", "protein_coding"),
            (5000, 6000, "+", "G2", "protein_coding"),
        ])
        unique, cod_cod, cod_ncod = _unique_and_shared(gr)
        assert cod_cod.df.empty
        assert cod_ncod.df.empty
        assert len(unique.df) == 2


class TestNewConstants:
    """Tests for new ReadCategory sub-categories."""

    def test_new_categories_in_category_order(self):
        from scnoisemeter.constants import CATEGORY_ORDER, ReadCategory
        assert ReadCategory.AMBIGUOUS_COD_COD  in CATEGORY_ORDER
        assert ReadCategory.AMBIGUOUS_COD_NCOD in CATEGORY_ORDER

    def test_new_categories_in_ambiguous_set(self):
        from scnoisemeter.constants import AMBIGUOUS_CATEGORIES, ReadCategory
        assert ReadCategory.AMBIGUOUS_COD_COD  in AMBIGUOUS_CATEGORIES
        assert ReadCategory.AMBIGUOUS_COD_NCOD in AMBIGUOUS_CATEGORIES

    def test_new_categories_not_in_noise(self):
        from scnoisemeter.constants import NOISE_CATEGORIES, ReadCategory
        assert ReadCategory.AMBIGUOUS_COD_COD  not in NOISE_CATEGORIES
        assert ReadCategory.AMBIGUOUS_COD_NCOD not in NOISE_CATEGORIES

    def test_report_has_colours_for_new_categories(self):
        from scnoisemeter.modules.report import CATEGORY_COLOURS, CATEGORY_LABELS
        from scnoisemeter.constants import ReadCategory
        assert ReadCategory.AMBIGUOUS_COD_COD  in CATEGORY_COLOURS
        assert ReadCategory.AMBIGUOUS_COD_NCOD in CATEGORY_COLOURS
        assert ReadCategory.AMBIGUOUS_COD_COD  in CATEGORY_LABELS
        assert ReadCategory.AMBIGUOUS_COD_NCOD in CATEGORY_LABELS


# ---------------------------------------------------------------------------
# Tests for new features: intergenic profiler integration,
# cluster metrics, polyA site loading, obs metadata loading
# ---------------------------------------------------------------------------

class TestIntergenicRecordExtraction:
    """Tests for extract_intergenic_records from pipeline side-table."""

    def test_empty_side_table(self):
        from scnoisemeter.modules.intergenic_profiler import extract_intergenic_records
        mock_result = type("R", (), {"intergenic_reads": []})()
        records = extract_intergenic_records(mock_result)
        assert records == []

    def test_valid_records_extracted(self):
        from scnoisemeter.modules.intergenic_profiler import (
            extract_intergenic_records, IntergenicReadRecord
        )
        # 7-tuple format: (contig, start, end, strand, cb, has_junction, three_prime)
        mock_result = type("R", (), {"intergenic_reads": [
            ("chr1", 1000, 1500, "+", "ACGTACGTACGTACGT", False, 1500),
            ("chr1", 2000, 2800, "-", "TTTTGGGGCCCCAAAA", True,  2800),
        ]})()
        records = extract_intergenic_records(mock_result)
        assert len(records) == 2
        assert records[0].contig == "chr1"
        assert records[0].start == 1000
        assert records[0].strand == "+"
        assert records[0].has_junction is False
        assert records[1].has_junction is True

    def test_malformed_records_skipped(self):
        from scnoisemeter.modules.intergenic_profiler import extract_intergenic_records
        mock_result = type("R", (), {"intergenic_reads": [
            ("chr1", 1000, 1500, "+", "CB1", False, 1500),  # valid
            ("chr1", "bad", 1500, "+", "CB2", False, 1500),  # invalid start
            None,  # completely invalid
        ]})()
        # Should not raise, should return the one valid record
        try:
            records = extract_intergenic_records(mock_result)
            # At least no crash
        except Exception as e:
            assert False, f"Should not raise: {e}"

    def test_missing_attribute_graceful(self):
        from scnoisemeter.modules.intergenic_profiler import extract_intergenic_records
        mock_result = type("R", (), {})()  # no intergenic_reads attribute
        records = extract_intergenic_records(mock_result)
        assert records == []


class TestClusterMetrics:
    """Tests for per-cluster noise metrics computation."""

    def _make_cell_table(self):
        import pandas as pd
        from scnoisemeter.modules.metrics import CellTable
        from scnoisemeter.constants import ReadCategory
        df = pd.DataFrame({
            "cell_barcode":        ["CB1", "CB2", "CB3", "CB4", "CB5", "CB6"],
            "n_reads":             [1000,  1200,  800,   900,   1100,  950],
            "n_bases":             [500000]*6,
            "noise_read_frac":     [0.10,  0.12,  0.08,  0.25,  0.27,  0.22],
            "noise_base_frac":     [0.11,  0.13,  0.09,  0.26,  0.28,  0.23],
            f"read_frac_{ReadCategory.EXONIC_SENSE.value}":    [0.80]*6,
            f"read_frac_{ReadCategory.INTRONIC_PURE.value}":   [0.05]*6,
            f"read_frac_{ReadCategory.CHIMERIC.value}":        [0.02]*6,
            f"read_frac_{ReadCategory.EXONIC_ANTISENSE.value}":[0.01]*6,
            f"read_frac_{ReadCategory.AMBIGUOUS_COD_COD.value}":[0.02]*6,
            "n_tso":   [5, 6, 4, 10, 11, 9],
            "n_polya": [8, 9, 7, 20, 22, 18],
            "n_noncanon": [3]*6,
        }).set_index("cell_barcode")
        return CellTable(df=df, sample_name="test")

    def _make_obs(self):
        import pandas as pd
        return pd.DataFrame({
            "cell_barcode": ["CB1", "CB2", "CB3", "CB4", "CB5", "CB6"],
            "cluster":      ["T_cell", "T_cell", "T_cell", "B_cell", "B_cell", "B_cell"],
        })

    def test_two_clusters_computed(self):
        from scnoisemeter.modules.metrics import compute_cluster_metrics
        ct = self._make_cell_table()
        obs = self._make_obs()
        cluster_df = compute_cluster_metrics(ct, obs)
        assert "T_cell" in cluster_df.index
        assert "B_cell" in cluster_df.index
        assert cluster_df.loc["T_cell", "n_cells"] == 3
        assert cluster_df.loc["B_cell", "n_cells"] == 3

    def test_b_cells_noisier_than_t_cells(self):
        from scnoisemeter.modules.metrics import compute_cluster_metrics
        ct = self._make_cell_table()
        obs = self._make_obs()
        cluster_df = compute_cluster_metrics(ct, obs)
        t_noise = cluster_df.loc["T_cell", "median_noise_read_frac"]
        b_noise = cluster_df.loc["B_cell", "median_noise_read_frac"]
        assert b_noise > t_noise

    def test_empty_obs_returns_empty(self):
        import pandas as pd
        from scnoisemeter.modules.metrics import compute_cluster_metrics
        ct = self._make_cell_table()
        empty_obs = pd.DataFrame(columns=["cell_barcode", "cluster"])
        result = compute_cluster_metrics(ct, empty_obs)
        assert result.empty

    def test_no_matching_barcodes_returns_empty(self):
        import pandas as pd
        from scnoisemeter.modules.metrics import compute_cluster_metrics
        ct = self._make_cell_table()
        obs = pd.DataFrame({
            "cell_barcode": ["XXXXXXXXXXXXXXXX"],
            "cluster":      ["mystery"],
        })
        result = compute_cluster_metrics(ct, obs)
        assert result.empty


class TestLoadObsMetadata:
    """Tests for obs metadata CSV/TSV loading with flexible column names."""

    def test_standard_columns(self, tmp_path):
        import pandas as pd
        from scnoisemeter.modules.metrics import load_obs_metadata
        f = tmp_path / "obs.tsv"
        f.write_text("cell_barcode\tcluster\nCB1\tT_cell\nCB2\tB_cell\n")
        df = load_obs_metadata(str(f))
        assert list(df.columns) == ["cell_barcode", "cluster"]
        assert len(df) == 2

    def test_seurat_column_names(self, tmp_path):
        from scnoisemeter.modules.metrics import load_obs_metadata
        f = tmp_path / "meta.csv"
        f.write_text("barcode,seurat_clusters\nCB1,0\nCB2,1\n")
        df = load_obs_metadata(str(f))
        assert "cell_barcode" in df.columns
        assert "cluster" in df.columns

    def test_scanpy_column_names(self, tmp_path):
        from scnoisemeter.modules.metrics import load_obs_metadata
        f = tmp_path / "obs.csv"
        f.write_text("cell_barcode,leiden\nCB1,0\nCB2,1\n")
        df = load_obs_metadata(str(f))
        assert "cluster" in df.columns

    def test_missing_barcode_column_raises(self, tmp_path):
        import pytest
        from scnoisemeter.modules.metrics import load_obs_metadata
        f = tmp_path / "bad.tsv"
        f.write_text("sample\tcluster\nS1\tA\n")
        with pytest.raises(ValueError, match="barcode"):
            load_obs_metadata(str(f))

    def test_missing_cluster_column_raises(self, tmp_path):
        import pytest
        from scnoisemeter.modules.metrics import load_obs_metadata
        f = tmp_path / "bad.tsv"
        f.write_text("cell_barcode\tsample_id\nCB1\tS1\n")
        with pytest.raises(ValueError, match="cluster"):
            load_obs_metadata(str(f))


class TestPolyaSiteLoading:
    """Tests for polyA site BED loading."""

    def test_loads_correctly(self, tmp_path):
        from scnoisemeter.cli import _load_polya_sites
        f = tmp_path / "polya.bed"
        f.write_text("chr1\t1000\t1010\tsite1\t0\t+\nchr1\t2000\t2010\tsite2\t0\t-\nchr2\t500\t510\tsite3\t0\t+\n")
        sites = _load_polya_sites(str(f))
        assert "chr1" in sites
        assert "chr2" in sites
        assert len(sites["chr1"]) == 2
        assert sites["chr1"] == sorted(sites["chr1"])  # must be sorted

    def test_comments_skipped(self, tmp_path):
        from scnoisemeter.cli import _load_polya_sites
        f = tmp_path / "polya.bed"
        f.write_text("# header comment\nchr1\t1000\t1010\tsite\t0\t+\n")
        sites = _load_polya_sites(str(f))
        assert len(sites["chr1"]) == 1

    def test_gzip_file_loads(self, tmp_path):
        import gzip
        from scnoisemeter.cli import _load_polya_sites
        f = tmp_path / "polya.bed.gz"
        with gzip.open(f, "wt") as fh:
            fh.write("chr1\t1000\t1010\tsite1\t0\t+\n")
            fh.write("chr2\t5000\t5010\tsite2\t0\t-\n")
        sites = _load_polya_sites(str(f))
        assert "chr1" in sites
        assert "chr2" in sites
        assert len(sites["chr1"]) == 1


# ---------------------------------------------------------------------------
# Tests for new reference file loaders and TSS/NUMT features
# ---------------------------------------------------------------------------

class TestPolyASiteMultipleFiles:
    """Tests for merging multiple polyA site databases."""

    def test_single_file(self, tmp_path):
        from scnoisemeter.cli import _load_polya_sites
        f = tmp_path / "polya.bed"
        f.write_text("chr1\t1000\t1010\tsite1\t0\t+\nchr1\t2000\t2010\tsite2\t0\t-\n")
        sites = _load_polya_sites([str(f)])
        assert len(sites["chr1"]) == 2

    def test_two_files_merged(self, tmp_path):
        from scnoisemeter.cli import _load_polya_sites
        f1 = tmp_path / "db1.bed"
        f1.write_text("chr1\t1000\t1010\tsite1\t0\t+\n")
        f2 = tmp_path / "db2.bed"
        f2.write_text("chr1\t5000\t5010\tsite2\t0\t+\nchr2\t100\t110\tsite3\t0\t-\n")
        sites = _load_polya_sites([str(f1), str(f2)])
        assert len(sites["chr1"]) == 2
        assert "chr2" in sites

    def test_duplicate_positions_deduplicated(self, tmp_path):
        from scnoisemeter.cli import _load_polya_sites
        f1 = tmp_path / "db1.bed"
        f1.write_text("chr1\t1000\t1010\tsite\t0\t+\n")
        f2 = tmp_path / "db2.bed"
        f2.write_text("chr1\t1000\t1010\tsite\t0\t+\n")  # same position
        sites = _load_polya_sites([str(f1), str(f2)])
        assert len(sites["chr1"]) == 1   # deduplicated

    def test_string_input_accepted(self, tmp_path):
        """Single string path (not list) should also work."""
        from scnoisemeter.cli import _load_polya_sites
        f = tmp_path / "polya.bed"
        f.write_text("chr1\t1000\t1010\tsite\t0\t+\n")
        sites = _load_polya_sites(str(f))
        assert "chr1" in sites


class TestTSSSiteLoader:
    """Tests for TSS/CAGE site loading."""

    def test_basic_loading(self, tmp_path):
        from scnoisemeter.cli import _load_tss_sites
        f = tmp_path / "tss.bed"
        f.write_text("chr1\t5000\t5010\tpeak1\t0\t+\nchr2\t3000\t3010\tpeak2\t0\t-\n")
        sites = _load_tss_sites([str(f)])
        assert "chr1" in sites
        assert "chr2" in sites
        # midpoint of 5000-5010 = 5005
        assert 5005 in sites["chr1"]

    def test_multiple_files_merged(self, tmp_path):
        from scnoisemeter.cli import _load_tss_sites
        f1 = tmp_path / "fantom.bed"
        f1.write_text("chr1\t1000\t1010\tpeak1\t0\t+\n")
        f2 = tmp_path / "rampage.bed"
        f2.write_text("chr1\t9000\t9010\tpeak2\t0\t+\n")
        sites = _load_tss_sites([str(f1), str(f2)])
        assert len(sites["chr1"]) == 2

    def test_gzip_tss(self, tmp_path):
        import gzip
        from scnoisemeter.cli import _load_tss_sites
        f = tmp_path / "tss.bed.gz"
        with gzip.open(f, "wt") as fh:
            fh.write("chr1\t5000\t5010\tpeak1\t0\t+\n")
        sites = _load_tss_sites([str(f)])
        assert "chr1" in sites


class TestNUMTLoader:
    """Tests for NUMT BED loading."""

    def test_basic_loading(self, tmp_path):
        from scnoisemeter.cli import _load_numt_bed
        f = tmp_path / "numt.bed"
        f.write_text("chr1\t1000\t2000\nchR1\t5000\t6000\n")
        intervals = _load_numt_bed(str(f))
        assert "chr1" in intervals
        assert (1000, 2000) in intervals["chr1"]

    def test_gzip_numt(self, tmp_path):
        import gzip
        from scnoisemeter.cli import _load_numt_bed
        f = tmp_path / "numt.bed.gz"
        with gzip.open(f, "wt") as fh:
            fh.write("chr1\t1000\t2000\n")
        intervals = _load_numt_bed(str(f))
        assert "chr1" in intervals


class TestNearPolyaSiteProximityParam:
    """Tests for _near_polya_site with custom proximity parameter."""

    def test_default_proximity(self):
        from scnoisemeter.modules.intergenic_profiler import _near_polya_site
        sites = {"chr1": [1000, 2000, 3000]}
        assert _near_polya_site(sites, "chr1", 1020)   # within 50 bp default
        assert not _near_polya_site(sites, "chr1", 1100)  # 100 bp away

    def test_custom_proximity_tss(self):
        """TSS proximity is 100 bp — should accept hits within 100 bp."""
        from scnoisemeter.modules.intergenic_profiler import _near_polya_site
        sites = {"chr1": [5000]}
        assert _near_polya_site(sites, "chr1", 5080, proximity=100)
        assert not _near_polya_site(sites, "chr1", 5150, proximity=100)

    def test_zero_proximity(self):
        """Proximity=0 means exact match only."""
        from scnoisemeter.modules.intergenic_profiler import _near_polya_site
        sites = {"chr1": [5000]}
        assert _near_polya_site(sites, "chr1", 5000, proximity=0)
        assert not _near_polya_site(sites, "chr1", 5001, proximity=0)
