"""
tests/test_cli_run.py
=====================
`scnoisemeter run` through the real command line, on the tiny fixtures.

The unit tests build result objects by hand. Nothing else ran the command a
user types, so a pre-write check that rejected every real barcode-aware result
set shipped unnoticed.
"""

from __future__ import annotations

import json

from click.testing import CliRunner

from scnoisemeter.cli import cli


def _run(tiny_bam, tiny_gtf, tmp_path, *extra):
    polya = tmp_path / "polya.bed"
    polya.write_text("chr1\t5990\t6010\tsite\t1\t+\n")
    out = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "run", "--bam", str(tiny_bam), "--gtf", str(tiny_gtf),
        "--polya-sites", str(polya), "--tss-db", "none",
        "--platform", "ont", "--library-strand", "stranded",
        "--threads", "1", "--no-cache", "--output-dir", str(out), *extra,
    ])
    return result, out


def _metrics(out):
    rows = (out / "tiny.read_metrics.tsv").read_text().splitlines()
    return dict(line.split("\t", 1) for line in rows)


def test_run_completes_and_writes_the_documented_files(tiny_bam, tiny_gtf, tmp_path):
    result, out = _run(tiny_bam, tiny_gtf, tmp_path)
    assert result.exit_code == 0, result.output
    for suffix in (".read_metrics.tsv", ".cell_metrics.tsv", ".intergenic_loci.tsv",
                   ".multiqc.json", ".run_info.json", ".qc.json", ".report.html"):
        assert (out / f"tiny{suffix}").exists(), suffix
    qc = json.loads((out / "tiny.qc.json").read_text())
    assert qc["n_failed_error"] == 0


def test_read_metrics_carries_every_flag_count(tiny_bam, tiny_gtf, tmp_path):
    result, out = _run(tiny_bam, tiny_gtf, tmp_path)
    assert result.exit_code == 0, result.output
    metrics = _metrics(out)
    for key in ("n_tso_invasion", "n_polya_priming", "n_polya_tail",
                "n_noncanon_junction", "n_rt_switch_repeat", "n_tso_concatemer"):
        assert key in metrics, key


def test_run_info_records_the_settings_that_change_a_number(tiny_bam, tiny_gtf, tmp_path):
    result, out = _run(tiny_bam, tiny_gtf, tmp_path,
                       "--polya-min-support", "25", "--priming-criterion", "perc_a")
    assert result.exit_code == 0, result.output
    info = json.loads((out / "tiny.run_info.json").read_text())
    assert info["priming"] == {"polya_min_support": 25.0, "criterion": "perc_a",
                               "anchor_veto_internal_priming": False}
    assert info["intergenic_reclassification_sampled"] is False


def test_cell_barcodes_matching_nothing_stops_the_run(tiny_bam, tiny_gtf, tmp_path):
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("TTTTTTTTTTTTTTTT-1\n")
    result, out = _run(tiny_bam, tiny_gtf, tmp_path, "--cell-barcodes", str(barcodes))
    assert result.exit_code != 0
    assert "0 reads remain" in result.output
    assert not (out / "tiny.read_metrics.tsv").exists()
