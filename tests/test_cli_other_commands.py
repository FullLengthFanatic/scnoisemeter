"""
tests/test_cli_other_commands.py
================================
`compare`, `discover` and `run-plate` through the real command line, on the
tiny fixtures.

tests/test_cli_run.py does the same for `run`. Until these existed the other
three commands were exercised only through mocks, and none of them ran the
pre-write checks that `run` has.
"""

from __future__ import annotations

import json

import pysam
from click.testing import CliRunner

from scnoisemeter.cli import cli


def _copy_bam(src, dst, *, program=None, drop_tags=()):
    """Copy a sorted BAM, optionally renaming its @PG program or dropping tags."""
    with pysam.AlignmentFile(str(src)) as fin:
        header = fin.header.to_dict()
        if program:
            header["PG"] = [{"ID": program, "PN": program, "VN": "1.0"}]
        out_header = pysam.AlignmentHeader.from_dict(header)
        with pysam.AlignmentFile(str(dst), "wb", header=out_header) as fout:
            for read in fin:
                fields = read.to_dict()
                fields["tags"] = [t for t in fields["tags"] if t[:2] not in drop_tags]
                fout.write(pysam.AlignedSegment.from_dict(fields, out_header))
    pysam.index(str(dst))
    return dst


def _annotation_args(tiny_gtf, tmp_path):
    polya = tmp_path / "polya.bed"
    polya.write_text("chr1\t5990\t6010\tsite\t1\t+\n")
    return ["--gtf", str(tiny_gtf), "--polya-sites", str(polya), "--tss-db", "none"]


def _qc(path):
    return json.loads(path.read_text())


def _metrics(path):
    return dict(line.split("\t", 1) for line in path.read_text().splitlines())


# ---------------------------------------------------------------------------
# compare
# ---------------------------------------------------------------------------

def _compare(tiny_bam, tiny_gtf, tmp_path, *extra):
    out = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "compare", "--bam-a", str(tiny_bam), "--bam-b", str(tiny_bam),
        *_annotation_args(tiny_gtf, tmp_path),
        "--platform", "ont", "--library-strand", "stranded",
        "--threads", "1", "--no-cache", "--output-dir", str(out), *extra,
    ])
    return result, out


def test_compare_writes_the_documented_files_for_a_nested_pair(tiny_bam, tiny_gtf, tmp_path):
    result, out = _compare(tiny_bam, tiny_gtf, tmp_path)
    assert result.exit_code == 0, result.output
    for name in ("comparison.metrics.tsv", "comparison.stats.tsv", "comparison.report.html",
                 "comparison.retention.tsv", "comparison.transitions.tsv",
                 "comparison.matching.tsv"):
        assert (out / name).exists(), name


def test_compare_checks_both_result_sets_before_writing(tiny_bam, tiny_gtf, tmp_path):
    result, out = _compare(tiny_bam, tiny_gtf, tmp_path)
    assert result.exit_code == 0, result.output
    for label in ("sample_A", "sample_B"):
        assert _qc(out / f"{label}.qc.json")["n_failed_error"] == 0


def test_compare_stops_when_the_whitelist_matches_no_barcode(tiny_bam, tiny_gtf, tmp_path):
    whitelist = tmp_path / "whitelist.txt"
    whitelist.write_text("TTTTTTTTTTTTTTTT\n")
    result, out = _compare(tiny_bam, tiny_gtf, tmp_path, "--barcode-whitelist", str(whitelist))
    assert result.exit_code != 0
    assert "carry no usable cell barcode" in result.output
    assert not (out / "comparison.metrics.tsv").exists()


# ---------------------------------------------------------------------------
# discover
# ---------------------------------------------------------------------------

def _discover(tiny_bam, tiny_gtf, tiny_fasta, tmp_path):
    bam_dir = tmp_path / "bams"
    bam_dir.mkdir()
    # The header has to name the platform: discover skips a BAM it cannot place.
    _copy_bam(tiny_bam, bam_dir / "cells.bam", program="pbmm2")
    out = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "discover", "--bam-dir", str(bam_dir), "--reference", str(tiny_fasta),
        *_annotation_args(tiny_gtf, tmp_path),
        "--run-all", "--threads", "1", "--output-dir", str(out),
    ])
    return result, out / "cells"


def test_discover_run_all_writes_a_result_set_per_bam(tiny_bam, tiny_gtf, tiny_fasta, tmp_path):
    result, out = _discover(tiny_bam, tiny_gtf, tiny_fasta, tmp_path)
    assert result.exit_code == 0, result.output
    for suffix in (".read_metrics.tsv", ".cell_metrics.tsv", ".multiqc.json",
                   ".run_info.json", ".report.html"):
        assert (out / f"cells{suffix}").exists(), suffix
    assert _metrics(out / "cells.read_metrics.tsv")["n_cells"] == "2"


def test_discover_checks_the_result_set_before_writing(tiny_bam, tiny_gtf, tiny_fasta, tmp_path):
    result, out = _discover(tiny_bam, tiny_gtf, tiny_fasta, tmp_path)
    assert result.exit_code == 0, result.output
    assert _qc(out / "cells.qc.json")["n_failed_error"] == 0


# ---------------------------------------------------------------------------
# run-plate
# ---------------------------------------------------------------------------

def _run_plate(tiny_bam, tiny_gtf, tmp_path, *extra, out_name="out"):
    plate_dir = tmp_path / "plate"
    if not plate_dir.exists():
        for well in ("A1", "A2"):
            well_dir = plate_dir / f"P1_{well}"
            well_dir.mkdir(parents=True)
            # Plate wells carry no cell barcode: the well is the cell.
            _copy_bam(tiny_bam, well_dir / "well.bam", drop_tags=("CB", "UB"))
    out = tmp_path / out_name
    result = CliRunner().invoke(cli, [
        "run-plate", "--plate-dir", str(plate_dir),
        *_annotation_args(tiny_gtf, tmp_path),
        "--platform", "smartseq", "--threads", "1", "--output-dir", str(out), *extra,
    ])
    return result, out / "P1"


def test_run_plate_writes_plate_and_per_well_tables(tiny_bam, tiny_gtf, tmp_path):
    result, out = _run_plate(tiny_bam, tiny_gtf, tmp_path)
    assert result.exit_code == 0, result.output
    for suffix in (".read_metrics.tsv", ".cell_metrics.tsv", ".per_well_metrics.tsv",
                   ".run_info.json", ".report.html"):
        assert (out / f"P1{suffix}").exists(), suffix
    wells = (out / "P1.per_well_metrics.tsv").read_text().splitlines()
    assert len(wells) == 3                      # header and two wells
    assert _metrics(out / "P1.read_metrics.tsv")["n_cells"] == "2"


def test_run_plate_checks_the_plate_result_before_writing(tiny_bam, tiny_gtf, tmp_path):
    result, out = _run_plate(tiny_bam, tiny_gtf, tmp_path)
    assert result.exit_code == 0, result.output
    assert _qc(out / "P1.qc.json")["n_failed_error"] == 0


def test_run_plate_parallel_wells_reports_what_serial_does(tiny_bam, tiny_gtf, tmp_path):
    serial, out_serial = _run_plate(tiny_bam, tiny_gtf, tmp_path, out_name="serial")
    parallel, out_parallel = _run_plate(tiny_bam, tiny_gtf, tmp_path,
                                        "--parallel-wells", "2", out_name="parallel")
    assert serial.exit_code == 0, serial.output
    assert parallel.exit_code == 0, parallel.output
    assert (_metrics(out_parallel / "P1.read_metrics.tsv")
            == _metrics(out_serial / "P1.read_metrics.tsv"))
