"""
tests/test_reference_dependent_blanks.py
========================================
Four checks read the genome sequence: internal priming, non-canonical
junctions, RT-switch repeats and the intergenic hotspot call. Without
`--reference` none of them runs, and a 0 would say "looked, found none".
"""

from __future__ import annotations

import json

import pandas as pd
from click.testing import CliRunner

from scnoisemeter.cli import cli
from scnoisemeter.modules.metrics import SampleMetrics
from scnoisemeter.modules.report_figures import _artifact_flags

NEEDS_GENOME = ("n_polya_priming", "n_noncanon_junction", "n_rt_switch_repeat")
HOTSPOT = ("read_frac_intergenic_hotspot", "base_frac_intergenic_hotspot")


def _annotation_args(tiny_gtf, tmp_path):
    polya = tmp_path / "polya.bed"
    polya.write_text("chr1\t5990\t6010\tsite\t1\t+\n")
    return ["--gtf", str(tiny_gtf), "--polya-sites", str(polya), "--tss-db", "none",
            "--library-strand", "stranded", "--threads", "1", "--no-cache"]


def _run(tiny_bam, tiny_gtf, tmp_path, *extra):
    out = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "run", "--bam", str(tiny_bam), *_annotation_args(tiny_gtf, tmp_path),
        "--platform", "ont", "--output-dir", str(out), *extra,
    ])
    assert result.exit_code == 0, result.output
    return out


def _metrics(out):
    rows = (out / "tiny.read_metrics.tsv").read_text().splitlines()
    return dict(line.split("\t", 1) for line in rows)


def test_read_metrics_leaves_the_genome_checks_blank_without_reference(tiny_bam, tiny_gtf, tmp_path):
    metrics = _metrics(_run(tiny_bam, tiny_gtf, tmp_path))
    for key in NEEDS_GENOME:
        assert metrics[key] == "", key
    for key in HOTSPOT:
        assert key not in metrics, key
    assert metrics["n_tso_invasion"] == "0"      # needs no genome, so it was measured


def test_read_metrics_reports_the_genome_checks_with_reference(tiny_bam, tiny_gtf, tiny_fasta, tmp_path):
    metrics = _metrics(_run(tiny_bam, tiny_gtf, tmp_path, "--reference", str(tiny_fasta)))
    for key in NEEDS_GENOME:
        assert metrics[key].isdigit(), key
    for key in HOTSPOT:
        assert float(metrics[key]) >= 0, key


def test_multiqc_json_omits_the_genome_checks_without_reference(tiny_bam, tiny_gtf, tmp_path):
    out = _run(tiny_bam, tiny_gtf, tmp_path)
    data = json.loads((out / "tiny.multiqc.json").read_text())["data"]["tiny"]
    for key in ("polya_priming_frac", "rt_switch_repeat_frac", "read_frac_intergenic_hotspot"):
        assert key not in data, key
    assert "tso_invasion_frac" in data


def test_cell_metrics_leave_the_genome_checks_blank_without_reference(tiny_bam, tiny_gtf, tmp_path):
    out = _run(tiny_bam, tiny_gtf, tmp_path)
    cells = pd.read_csv(out / "tiny.cell_metrics.tsv", sep="\t", index_col=0)
    for column in ("n_polya", "n_noncanon", *HOTSPOT):
        assert cells[column].isna().all(), column
    assert cells["n_tso"].notna().all()


def test_compare_does_not_report_a_hotspot_share_without_reference(tiny_bam, tiny_gtf, tmp_path):
    out = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "compare", "--bam-a", str(tiny_bam), "--bam-b", str(tiny_bam),
        *_annotation_args(tiny_gtf, tmp_path), "--platform", "ont", "--output-dir", str(out),
    ])
    assert result.exit_code == 0, result.output
    stats = pd.read_csv(out / "comparison.stats.tsv", sep="\t", index_col="category")
    assert stats.loc["intergenic_hotspot", ["frac_sample_A", "frac_sample_B", "delta"]].isna().all()
    assert stats.loc["exonic_sense", ["frac_sample_A", "frac_sample_B", "delta"]].notna().all()


def test_report_names_the_flags_that_were_not_measured():
    sm = SampleMetrics(sample_name="s", bam_path="s.bam", platform="ont",
                       pipeline_stage="pre_filter", aligner="minimap2")
    sm.n_reads_classified = 1000
    sm.n_tso_invasion = 10
    sm.n_polya_priming = None
    sm.n_noncanon_junction = None
    bar = _artifact_flags(sm).data[0]
    text = dict(zip(bar.y, bar.text))
    assert text["TSO invasion"] == "1.00%"
    assert text["Internal polyA priming"] == "not measured (no --reference)"
    assert text["Non-canonical junction"] == "not measured (no --reference)"


def test_a_reference_that_cannot_be_opened_counts_as_no_reference(tiny_bam, tiny_gtf, tmp_path):
    broken = tmp_path / "broken.fa.gz"
    broken.write_text("not a bgzipped FASTA\n")
    metrics = _metrics(_run(tiny_bam, tiny_gtf, tmp_path, "--reference", str(broken)))
    for key in NEEDS_GENOME:
        assert metrics[key] == "", key
