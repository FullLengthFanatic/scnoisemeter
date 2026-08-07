"""
Tests for the `cohort` reader and its comparability guards.

The rules under test exist because the author's own result directories break
them: half predate provenance recording, half were classified in
barcode-agnostic mode, and two were produced by different tool versions whose
TSO metric differs by ~700x for the same BAM.
"""

import json

import pandas as pd
import pytest

from scnoisemeter.modules.cohort import (
    build_composition_table,
    build_summary_table,
    load_cohort,
    ranking_metric,
)


def _write_results_dir(tmp_path, name, metrics, *, run_info=None, cells=None):
    """Create a minimal result directory of the shape `run` writes."""
    d = tmp_path / name
    d.mkdir()
    with open(d / f"{name}.read_metrics.tsv", "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in metrics.items():
            fh.write(f"{k}\t{v}\n")
    if run_info is not None:
        (d / f"{name}.run_info.json").write_text(json.dumps(run_info))
    if cells is not None:
        cells.to_csv(d / f"{name}.cell_metrics.tsv", sep="\t")
    return d


V072 = {
    "n_reads_classified": 1000,
    "n_cells": 2,
    "broad_noncanonical_read_frac": 0.20,
    "artifact_candidate_read_frac": 0.03,
    "read_frac_exonic_sense": 0.70,
    "read_frac_intronic_pure": 0.15,
    "read_frac_chimeric": 0.03,
    "n_tso_invasion": 10,
    "n_tso_concatemer": 0,
}

# Pre-v0.7 directories carry the old aggregate names and no *_strict row.
V061 = {
    "n_reads_classified": 2000,
    "n_cells": 1,
    "noise_read_frac": 0.10,
    "noise_base_frac": 0.11,
    "read_frac_exonic_sense": 0.85,
    "read_frac_intronic_pure": 0.08,
    "n_tso_invasion": 4,
}


class TestLegacyDirectories:
    def test_old_aggregate_names_are_mapped_forward(self, tmp_path):
        _write_results_dir(tmp_path, "old", V061)
        _write_results_dir(tmp_path, "new", V072)
        cohort = load_cohort([tmp_path / "old", tmp_path / "new"])
        old = next(s for s in cohort.samples if s.sample_id == "old")
        assert old.get("broad_noncanonical_read_frac") == 0.10
        assert old.get("broad_noncanonical_base_frac") == 0.11

    def test_metric_absent_from_the_file_stays_none(self, tmp_path):
        """Never zero: a metric that did not exist is not a measurement of none."""
        _write_results_dir(tmp_path, "old", V061)
        _write_results_dir(tmp_path, "new", V072)
        cohort = load_cohort([tmp_path / "old", tmp_path / "new"])
        old = next(s for s in cohort.samples if s.sample_id == "old")
        assert old.get("artifact_candidate_read_frac") is None
        assert old.flag_rate("n_tso_concatemer") is None

    def test_blank_value_is_absent_not_zero(self, tmp_path):
        metrics = dict(V072, n_numt_intervals_loaded="")
        _write_results_dir(tmp_path, "a", metrics)
        _write_results_dir(tmp_path, "b", V072)
        cohort = load_cohort([tmp_path / "a", tmp_path / "b"])
        assert cohort.samples[0].get("n_numt_intervals_loaded") is None

    def test_ranking_falls_back_when_preferred_metric_is_missing(self, tmp_path):
        _write_results_dir(tmp_path, "old", V061)
        _write_results_dir(tmp_path, "new", V072)
        cohort = load_cohort([tmp_path / "old", tmp_path / "new"])
        assert ranking_metric(cohort) == "broad_noncanonical_read_frac"

    def test_ranking_uses_artifact_candidates_when_all_report_it(self, tmp_path):
        _write_results_dir(tmp_path, "a", V072)
        _write_results_dir(tmp_path, "b", dict(V072, artifact_candidate_read_frac=0.05))
        cohort = load_cohort([tmp_path / "a", tmp_path / "b"])
        assert ranking_metric(cohort) == "artifact_candidate_read_frac"


class TestComparabilityWarnings:
    def test_version_mismatch_warns(self, tmp_path):
        _write_results_dir(tmp_path, "a", V072,
                           run_info={"scnoisemeter_version": "0.7.2", "platform": "pacbio"})
        _write_results_dir(tmp_path, "b", V072,
                           run_info={"scnoisemeter_version": "0.6.1", "platform": "pacbio"})
        cohort = load_cohort([tmp_path / "a", tmp_path / "b"])
        assert any("different scnoisemeter versions" in w for w in cohort.warnings)

    def test_same_version_does_not_warn(self, tmp_path):
        info = {"scnoisemeter_version": "0.8.0", "platform": "pacbio"}
        _write_results_dir(tmp_path, "a", V072, run_info=info)
        _write_results_dir(tmp_path, "b", V072, run_info=info)
        cohort = load_cohort([tmp_path / "a", tmp_path / "b"])
        assert not any("different scnoisemeter versions" in w for w in cohort.warnings)

    def test_gtf_mismatch_warns(self, tmp_path):
        _write_results_dir(tmp_path, "a", V072,
                           run_info={"scnoisemeter_version": "0.8.0", "gtf_version": 45})
        _write_results_dir(tmp_path, "b", V072,
                           run_info={"scnoisemeter_version": "0.8.0", "gtf_version": 42})
        cohort = load_cohort([tmp_path / "a", tmp_path / "b"])
        assert any("different GENCODE releases" in w for w in cohort.warnings)

    def test_mixed_strandedness_warns_and_is_flagged(self, tmp_path):
        _write_results_dir(tmp_path, "tenx", V072,
                           run_info={"scnoisemeter_version": "0.8.0", "unstranded": False})
        _write_results_dir(tmp_path, "smartseq", V072,
                           run_info={"scnoisemeter_version": "0.8.0", "unstranded": True})
        cohort = load_cohort([tmp_path / "tenx", tmp_path / "smartseq"])
        assert cohort.has_mixed_strandedness
        assert any("stranded and unstranded" in w for w in cohort.warnings)

    def test_missing_provenance_warns(self, tmp_path):
        _write_results_dir(tmp_path, "a", V072)
        _write_results_dir(tmp_path, "b", V072)
        cohort = load_cohort([tmp_path / "a", tmp_path / "b"])
        assert any("No run_info.json" in w for w in cohort.warnings)

    def test_barcode_agnostic_samples_are_named(self, tmp_path):
        cells = pd.DataFrame(
            {"broad_noncanonical_read_frac": [0.1, 0.2, 0.3]},
            index=["c1", "c2", "c3"],
        )
        _write_results_dir(tmp_path, "withcells", V072, cells=cells)
        _write_results_dir(tmp_path, "nocells", V072)
        cohort = load_cohort([tmp_path / "withcells", tmp_path / "nocells"])
        assert any("No per-cell data for" in w and "nocells" in w for w in cohort.warnings)
        assert {s.sample_id: s.has_cells for s in cohort.samples} == {
            "withcells": True, "nocells": False,
        }


class TestTables:
    def test_all_zero_categories_are_dropped(self, tmp_path):
        _write_results_dir(tmp_path, "a", V072)
        _write_results_dir(tmp_path, "b", V072)
        cohort = load_cohort([tmp_path / "a", tmp_path / "b"])
        metrics = cohort.category_metrics("read_frac_")
        assert "read_frac_exonic_sense" in metrics
        assert "read_frac_multimapper" not in metrics   # zero in every sample

    def test_summary_and_composition_have_one_row_per_sample(self, tmp_path):
        _write_results_dir(tmp_path, "a", V072)
        _write_results_dir(tmp_path, "b", V061)
        cohort = load_cohort([tmp_path / "a", tmp_path / "b"])
        assert len(build_summary_table(cohort)) == 2
        assert len(build_composition_table(cohort)) == 2

    def test_flag_rate_divides_by_classified_reads(self, tmp_path):
        _write_results_dir(tmp_path, "a", V072)
        _write_results_dir(tmp_path, "b", V072)
        cohort = load_cohort([tmp_path / "a", tmp_path / "b"])
        assert cohort.samples[0].flag_rate("n_tso_invasion") == pytest.approx(10 / 1000)

    def test_directory_without_metrics_is_skipped_with_a_warning(self, tmp_path):
        _write_results_dir(tmp_path, "a", V072)
        _write_results_dir(tmp_path, "b", V072)
        empty = tmp_path / "empty"
        empty.mkdir()
        cohort = load_cohort([tmp_path / "a", tmp_path / "b", empty])
        assert len(cohort) == 2
        assert any("no *.read_metrics.tsv" in w for w in cohort.warnings)


class TestSampleSheet:
    def test_sheet_controls_label_group_and_order(self, tmp_path):
        _write_results_dir(tmp_path, "a", V072)
        _write_results_dir(tmp_path, "b", dict(V072, artifact_candidate_read_frac=0.001))
        sheet = pd.DataFrame(
            [{"sample": "a", "label": "Kit A", "group": "long-read", "order": 1},
             {"sample": "b", "label": "Kit B", "group": "short-read", "order": 2}]
        ).set_index("sample")
        cohort = load_cohort([tmp_path / "a", tmp_path / "b"], sample_sheet=sheet)
        assert [s.label for s in cohort.samples] == ["Kit A", "Kit B"]
        assert [s.group for s in cohort.samples] == ["long-read", "short-read"]

    def test_without_sheet_samples_are_ordered_cleanest_first(self, tmp_path):
        _write_results_dir(tmp_path, "dirty", dict(V072, artifact_candidate_read_frac=0.30))
        _write_results_dir(tmp_path, "clean", dict(V072, artifact_candidate_read_frac=0.01))
        cohort = load_cohort([tmp_path / "dirty", tmp_path / "clean"])
        assert [s.sample_id for s in cohort.samples] == ["clean", "dirty"]


class TestFigures:
    """Figures must render without raising on the awkward real-world shapes."""

    def test_figures_render_for_a_mixed_cohort(self, tmp_path):
        from scnoisemeter.modules.report_figures import (
            _cohort_artifact_dots, _cohort_composition_bars,
            _cohort_deviation_heatmap, _cohort_percell_boxes,
        )
        cells = pd.DataFrame(
            {"broad_noncanonical_read_frac": [0.1, 0.2, 0.3, 0.4]},
            index=["c1", "c2", "c3", "c4"],
        )
        for i in range(4):
            _write_results_dir(
                tmp_path, f"s{i}",
                dict(V072, artifact_candidate_read_frac=0.01 * (i + 1)),
                run_info={"scnoisemeter_version": "0.8.0", "unstranded": i == 3},
                cells=cells if i < 2 else None,
            )
        cohort = load_cohort([tmp_path / f"s{i}" for i in range(4)])
        for fn in (_cohort_composition_bars, _cohort_deviation_heatmap,
                   _cohort_artifact_dots, _cohort_percell_boxes):
            assert fn(cohort).data, f"{fn.__name__} produced no traces"

    def test_near_zero_row_is_not_coloured_as_an_outlier(self, tmp_path):
        """A 0.0003 percentage-point spread must not paint the strongest cell."""
        from scnoisemeter.modules.report_figures import _cohort_deviation_heatmap
        for i in range(4):
            _write_results_dir(
                tmp_path, f"s{i}",
                dict(V072, read_frac_intergenic_hotspot=3e-6 if i == 0 else 0.0),
            )
        cohort = load_cohort([tmp_path / f"s{i}" for i in range(4)])
        fig = _cohort_deviation_heatmap(cohort)
        y = list(fig.data[0].y)
        row = y.index("Intergenic hotspot")
        assert set(fig.data[0].z[row]) == {0.0}
