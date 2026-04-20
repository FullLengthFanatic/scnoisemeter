"""
cli.py
======
Command-line interface for scNoiseMeter.

Subcommands
-----------
  run      Classify reads in a single BAM and produce QC metrics.
  compare  Run on two BAMs (e.g. pre- and post-filter) and produce a
           benchmarking comparison report.

Usage examples
--------------
  # Single-sample QC
  scnoisemeter run \\
    --bam sample.bam \\
    --gtf gencode.v45.gtf.gz \\
    --barcode-whitelist whitelist_10x_v3.txt \\
    --chemistry 10x_v3 \\
    --platform ont \\
    --threads 8 \\
    --output-dir results/ \\
    --sample-name my_sample

  # Pre/post filter benchmarking
  scnoisemeter compare \\
    --bam-a raw.bam \\
    --bam-b filtered.bam \\
    --gtf gencode.v45.gtf.gz \\
    --barcode-whitelist whitelist_10x_v3.txt \\
    --threads 8 \\
    --output-dir results/
"""

from __future__ import annotations

import hashlib
import json
import logging
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path
from typing import Optional

import click

from scnoisemeter import __version__
from scnoisemeter.constants import (
    BARCODE_AUTODETECT_MIN_FRACTION,
    Chemistry,
    DEFAULT_CHIMERIC_DISTANCE,
    DEFAULT_THREADS,
    Platform,
    PipelineStage,
    ReadCategory,
)
from scnoisemeter.modules.annotation import build_annotation_index
from scnoisemeter.modules.metrics import (
    compute_metrics, to_multiqc_json,
    compute_cluster_metrics, load_obs_metadata,
    compute_length_stratification,
)
from scnoisemeter.modules.intergenic_profiler import (
    profile_intergenic_loci, extract_intergenic_records,
    compute_intergenic_bases,
)
from scnoisemeter.modules.report import write_run_report, write_compare_report
from scnoisemeter.modules.pipeline import run_pipeline
from scnoisemeter.utils.bam_inspector import inspect_bam, BamMetadata
from scnoisemeter.utils.annotation_fetcher import (
    fetch_latest_gencode_gtf,
    fetch_gencode_gtf_version,
    fetch_latest_polyasite_atlas,
    fetch_polyadb4_atlas,
    fetch_fantom5_cage_peaks,
    fetch_10x_whitelist,
    extract_gencode_version_from_filename,
    extract_polyasite_version_from_filename,
)

logger = logging.getLogger("scnoisemeter")


# ---------------------------------------------------------------------------
# Shared options (reused across subcommands via a decorator)
# ---------------------------------------------------------------------------

def _shared_options(func):
    """Decorator that attaches all options common to both subcommands."""
    decorators = [
        click.option("--gtf",               required=False, default=None, type=click.Path(exists=True),
                     help="GENCODE GTF (plain or .gz).  If omitted, the latest GENCODE human annotation "
                          "is downloaded automatically to ~/.cache/scnoisemeter/. "
                          "Takes precedence over --gtf-version."),
        click.option("--gtf-version",       default=None,   type=int,
                     help="GENCODE release number to auto-download (e.g. 42). "
                          "Ignored when --gtf is supplied. "
                          "Use --gtf-version 42 to match PolyASite 3.0 and avoid version mismatch warnings."),
        click.option("--barcode-whitelist",  default=None,   type=click.Path(exists=True), help="File of valid corrected barcodes, one per line."),
        click.option("--whitelist-db",
                     default="none",
                     show_default=True,
                     type=click.Choice(["10x_v3", "10x_v4", "none"], case_sensitive=False),
                     help="10x barcode whitelist to auto-download when --barcode-whitelist is not supplied. "
                          "10x_v3: 3M Chromium v3 barcodes (2018); "
                          "10x_v4: 3M Chromium v4 / 3' GEX barcodes (2023); "
                          "none: no whitelist filter (default)."),
        click.option("--barcode-tag",        default="CB",   show_default=True, help="BAM tag for corrected cell barcode."),
        click.option("--umi-tag",            default="UB",   show_default=True, help="BAM tag for corrected UMI."),
        click.option("--chemistry",
                     type=click.Choice([c.value for c in Chemistry], case_sensitive=False),
                     default=Chemistry.TENX_V3.value, show_default=True,
                     help="Library chemistry (sets barcode length expectation)."),
        click.option("--platform",
                     type=click.Choice(["auto"] + [p.value for p in Platform], case_sensitive=False),
                     default="auto", show_default=True,
                     help="Sequencing platform.  'auto' detects from BAM header."),
        click.option("--pipeline-stage",
                     type=click.Choice(["auto"] + [s.value for s in PipelineStage], case_sensitive=False),
                     default="auto", show_default=True,
                     help="Processing stage of the BAM.  'auto' detects from BAM header."),
        click.option("--chimeric-distance",  default=DEFAULT_CHIMERIC_DISTANCE, show_default=True,
                     type=int, help="Max same-strand intra-chromosomal SA distance (bp) below which a split is NOT chimeric."),
        click.option("--repeats",            default=None,   type=click.Path(exists=True), help="RepeatMasker BED for hg38 (optional)."),
        click.option("--reference",          default=None,   type=click.Path(exists=True), help="Reference FASTA (.fa/.fa.gz + .fai) for polyA and junction checks."),
        click.option("--threads",            default=DEFAULT_THREADS, show_default=True, type=int, help="Parallel worker processes."),
        click.option("--no-umi-dedup",       is_flag=True,   help="Skip UMI set tracking (reduces memory for very large datasets)."),
        click.option("--no-cache",           is_flag=True,   help="Do not read or write the annotation cache."),
        click.option("--exclude-biotypes",   default=None,   multiple=True,
                     help="Gene biotypes to exclude from annotation (repeatable)."),
        click.option("--output-dir",         required=True,  type=click.Path(), help="Directory for output files (created if absent)."),
        click.option("--obs-metadata",       default=None,   type=click.Path(exists=True),
                     help="Per-cell metadata TSV with 'cell_barcode' and 'cluster' columns "
                          "(e.g. Seurat meta.data or Scanpy obs). Enables per-cluster noise profiles."),
        click.option("--polya-sites",        multiple=True,  type=click.Path(exists=True),
                     help="PolyA site BED file(s) for the 3\'-anchored full-length fraction. "
                          "Accepts plain .bed or .bed.gz. Repeat the flag to merge multiple "
                          "databases. When supplied, --polya-db is ignored."),
        click.option("--polya-db",
                     default="polyasite3",
                     show_default=True,
                     type=click.Choice(["polyasite3", "polyadb4", "both"], case_sensitive=False),
                     help="polyA database to auto-download when --polya-sites is not supplied. "
                          "polyasite3: PolyASite 3.0 (scRNA-seq atlas, GENCODE v42, ~569 k sites); "
                          "polyadb4: PolyA_DB v4 main collection (bulk + long-read validated, "
                          "hg38, ~281 k sites); "
                          "both: merge both databases for maximum site coverage."),
        click.option("--tss-sites",          multiple=True,  type=click.Path(exists=True),
                     help="CAGE peak / TSS BED file(s) for the 5\'-anchored full-length metric. "
                          "Accepts plain .bed or .bed.gz. Repeat for multiple databases. "
                          "When supplied, --tss-db is ignored."),
        click.option("--tss-db",
                     default="fantom5",
                     show_default=True,
                     type=click.Choice(["fantom5", "none"], case_sensitive=False),
                     help="TSS database to auto-download when --tss-sites is not supplied. "
                          "fantom5: FANTOM5 phase 1+2 CAGE peaks, hg38 (~180 k peaks); "
                          "none: disable TSS metric entirely."),
        click.option("--numt-bed",           default=None,   type=click.Path(exists=True),
                     help="NUMT BED file (nuclear mitochondrial DNA segments, hg38 coordinates). "
                          "When provided, mitochondrial reads are cross-checked against NUMT "
                          "intervals and flagged separately in the output."),
        click.option("--offline",            is_flag=True,
                     help="Use only cached annotation files; never make network calls. "
                          "Raises an error if the cache is empty. "
                          "Ignored when --gtf / --polya-sites are supplied explicitly."),
        click.option("--verbose", "-v",      is_flag=True,   help="Enable debug logging."),
    ]
    for dec in reversed(decorators):
        func = dec(func)
    return func


# ---------------------------------------------------------------------------
# CLI root group
# ---------------------------------------------------------------------------

@click.group()
@click.version_option(__version__, prog_name="scnoisemeter")
def cli():
    """scNoiseMeter — quantify technical noise in scRNA-seq BAM files."""
    pass


# ---------------------------------------------------------------------------
# `run` subcommand
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Annotation resolution helpers (GTF + polyA atlas auto-download)
# ---------------------------------------------------------------------------

def _resolve_gtf(
    gtf_arg: Optional[str],
    gtf_version: Optional[int] = None,
    offline: bool = False,
):
    """
    Resolve the GTF path.

    Priority: --gtf > --gtf-version > auto-download latest.

    Returns (path_str, version_int_or_None, source_str).
    """
    if gtf_arg is not None:
        version = extract_gencode_version_from_filename(gtf_arg)
        return gtf_arg, version, "user-supplied"

    if gtf_version is not None:
        path, version = fetch_gencode_gtf_version(gtf_version, offline=offline)
        return str(path), version, "auto-downloaded"

    path, version = fetch_latest_gencode_gtf(offline=offline)
    return str(path), version, "auto-downloaded"


def _resolve_polya_sites(
    polya_sites_arg,
    polya_db: str = "polyasite3",
    hint_gencode_version: Optional[int] = None,
    offline: bool = False,
):
    """
    Resolve polyA site file paths.

    If *polya_sites_arg* is non-empty, those paths are used directly and
    *polya_db* is ignored.  Otherwise, the database(s) selected by *polya_db*
    are auto-downloaded (or read from cache).

    *polya_db* choices:
      polyasite3  PolyASite 3.0 (scRNA-seq, GENCODE v42)
      polyadb4    PolyA_DB v4 main collection (bulk + long-read, hg38)
      both        merge both databases

    Returns (paths_list, version_int_or_None, source_str).
    Version is the GENCODE version embedded in the atlas filename (PolyASite
    convention); PolyA_DB v4 returns None since it uses NCBI/RefSeq.
    """
    if not polya_sites_arg:
        polya_db = polya_db.lower()
        if polya_db == "polyadb4":
            path = fetch_polyadb4_atlas(offline=offline)
            return [str(path)], None, "auto-downloaded"
        elif polya_db == "both":
            pa_path, pa_version = fetch_latest_polyasite_atlas(
                hint_max_gencode_version=hint_gencode_version,
                offline=offline,
            )
            db4_path = fetch_polyadb4_atlas(offline=offline)
            return [str(pa_path), str(db4_path)], pa_version, "auto-downloaded"
        else:  # polyasite3 (default)
            path, version = fetch_latest_polyasite_atlas(
                hint_max_gencode_version=hint_gencode_version,
                offline=offline,
            )
            return [str(path)], version, "auto-downloaded"

    paths = list(polya_sites_arg)
    # Detect version from first recognisable filename (PolyASite 3.0 naming convention)
    version: Optional[int] = None
    for p in paths:
        version = extract_polyasite_version_from_filename(p)
        if version is not None:
            break
    return paths, version, "user-supplied"


def _resolve_tss_sites(
    tss_sites_arg,
    tss_db: str = "fantom5",
    offline: bool = False,
) -> list:
    """
    Resolve TSS site file paths.

    If *tss_sites_arg* is non-empty, those paths are used directly and
    *tss_db* is ignored.  Otherwise, *tss_db* controls auto-download:
      fantom5  FANTOM5 phase 1+2 CAGE peaks (hg38)
      none     Return empty list; TSS metric is disabled.

    Returns a list of resolved path strings (may be empty).
    """
    if tss_sites_arg:
        return list(tss_sites_arg)
    if tss_db.lower() == "none":
        return []
    # fantom5 (default)
    path = fetch_fantom5_cage_peaks(offline=offline)
    return [str(path)]


def _resolve_whitelist(
    barcode_whitelist_arg: Optional[str],
    whitelist_db: str = "none",
    offline: bool = False,
) -> Optional[set]:
    """
    Resolve the barcode whitelist.

    Priority: --barcode-whitelist > --whitelist-db > None.
    """
    if barcode_whitelist_arg is not None:
        return _load_whitelist(barcode_whitelist_arg)
    if whitelist_db.lower() == "none":
        return None
    path = fetch_10x_whitelist(whitelist_db.lower(), offline=offline)
    return _load_whitelist(str(path))


def _check_version_consistency(gtf_version: Optional[int], polya_version: Optional[int]) -> Optional[str]:
    """
    Emit a warning when the GENCODE versions embedded in the GTF and polyA atlas
    filenames differ by more than 2 major releases.

    Returns the warning string if emitted, otherwise None.
    """
    if gtf_version is None or polya_version is None:
        return None
    diff = abs(gtf_version - polya_version)
    if diff > 5:
        msg = (
            f"Warning: GTF is GENCODE v{gtf_version} but polyA atlas was built on "
            f"GENCODE v{polya_version}. "
            f"Genes with 3\u2032 UTRs annotated between "
            f"v{min(gtf_version, polya_version)} and v{max(gtf_version, polya_version)} "
            f"may have reduced polyA anchoring scores. "
            f"To resolve: pass --gtf-version {polya_version} to pin the GTF to match "
            f"the atlas, or --polya-db polyadb4 to use a database without GENCODE versioning."
        )
        click.echo(msg, err=True)
        return msg
    return None


# ---------------------------------------------------------------------------
# BAM validation helpers (chromosome naming + sort order)
# ---------------------------------------------------------------------------

def _detect_chrom_style(names: list) -> str:
    """
    Return 'ucsc' if the majority of primary contigs start with 'chr',
    'ensembl' if they do not, or 'unknown' if the list is empty or mixed.

    Only primary autosomes + sex chromosomes are considered (1–22, X, Y, MT)
    to avoid confusion from alt/patch contigs.
    """
    primary = [n for n in names if n in {
        "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
        "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
        "chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM",
        "1","2","3","4","5","6","7","8","9","10","11","12","13","14",
        "15","16","17","18","19","20","21","22","X","Y","MT",
    }]
    if not primary:
        return "unknown"
    n_ucsc = sum(1 for n in primary if n.startswith("chr"))
    n_ens  = sum(1 for n in primary if not n.startswith("chr"))
    if n_ucsc > n_ens:
        return "ucsc"
    if n_ens > n_ucsc:
        return "ensembl"
    return "unknown"


def _gtf_chrom_style(gtf_path: str) -> str:
    """
    Peek at the first 200 non-comment, non-empty lines of *gtf_path*
    to determine whether it uses UCSC ('chr1') or Ensembl ('1') naming.
    """
    import gzip as _gz
    opener = _gz.open if gtf_path.endswith(".gz") else open
    n_ucsc = n_ens = 0
    with opener(gtf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            chrom = line.split("\t", 1)[0]
            if chrom.startswith("chr"):
                n_ucsc += 1
            else:
                n_ens += 1
            if n_ucsc + n_ens >= 200:
                break
    if n_ucsc > n_ens:
        return "ucsc"
    if n_ens > n_ucsc:
        return "ensembl"
    return "unknown"


def _validate_chromosome_naming(meta: BamMetadata, gtf_path: str) -> None:
    """
    Compare chromosome naming styles between the BAM header and the GTF.
    Aborts with a clear message if they differ.
    """
    bam_style = _detect_chrom_style(meta.reference_names)
    gtf_style = _gtf_chrom_style(gtf_path)

    if bam_style == "unknown" or gtf_style == "unknown":
        logger.warning(
            "Could not determine chromosome naming style for BAM (%s) or GTF (%s) — "
            "proceeding, but annotation matching may fail.",
            bam_style, gtf_style,
        )
        return

    if bam_style != gtf_style:
        bam_example = next(
            (n for n in meta.reference_names if n.startswith("chr") or n in {"1","2","X","Y","MT"}),
            meta.reference_names[0] if meta.reference_names else "?",
        )
        gtf_example = "chr1" if gtf_style == "ucsc" else "1"
        raise click.ClickException(
            f"Chromosome naming mismatch between BAM and GTF.\n"
            f"  BAM uses {bam_style.upper()} style  (e.g. '{bam_example}')\n"
            f"  GTF uses {gtf_style.upper()} style  (e.g. '{gtf_example}')\n"
            f"\n"
            f"To fix: either supply a GTF that matches the BAM's naming convention,\n"
            f"or remap the BAM chromosome names with 'samtools reheader'.\n"
            f"Do NOT rely on scNoiseMeter to silently remap — all annotations would\n"
            f"fail to match and all reads would appear intergenic."
        )


def _validate_chromosome_lengths(meta: BamMetadata, reference_path: Optional[str] = None) -> None:
    """
    Cross-check primary chromosome lengths in the BAM header against the reference
    genome to catch wrong-species or wrong-assembly alignments.

    Strategy (in order of availability):
    1. If *reference_path* is given and a .fai index exists, parse expected lengths from it.
    2. Otherwise fall back to hardcoded GRCh38 / hg38 expected lengths for key chromosomes.

    Issues a WARNING (non-fatal) when a significant mismatch is detected.
    """
    # GRCh38 / hg38 primary chromosome lengths (authoritative)
    _GRCH38_LENGTHS = {
        "chr1":  248_956_422,  "1":  248_956_422,
        "chr2":  242_193_529,  "2":  242_193_529,
        "chr3":  198_295_559,  "3":  198_295_559,
        "chr4":  190_214_555,  "4":  190_214_555,
        "chr5":  181_538_259,  "5":  181_538_259,
        "chrX":  156_040_895,  "X":  156_040_895,
    }

    expected_lengths: dict = {}

    if reference_path:
        fai_path = reference_path + ".fai"
        try:
            with open(fai_path) as fh:
                for line in fh:
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        expected_lengths[parts[0]] = int(parts[1])
        except (OSError, ValueError):
            pass  # Fall through to hardcoded lengths

    if not expected_lengths:
        expected_lengths = _GRCH38_LENGTHS

    # Check a few primary chromosomes
    mismatches = []
    for chrom, expected_len in expected_lengths.items():
        if chrom not in meta.reference_lengths:
            continue
        bam_len = meta.reference_lengths[chrom]
        # Flag if lengths differ by more than 1 Mb (clearly different assembly)
        if abs(bam_len - expected_len) > 1_000_000:
            mismatches.append((chrom, bam_len, expected_len))
        if len(mismatches) >= 3:
            break  # Enough evidence

    if mismatches:
        examples = "; ".join(
            f"{chrom}: BAM={bam_len:,} bp, expected={exp_len:,} bp"
            for chrom, bam_len, exp_len in mismatches[:3]
        )
        msg = (
            f"Chromosome length mismatch between BAM and reference genome. "
            f"{examples}. "
            f"This may indicate a wrong-species or wrong-assembly alignment. "
            f"Verify that the BAM was aligned to GRCh38/hg38."
        )
        logger.warning(msg)
        meta.warnings.append(msg)


def _validate_sq_lines(meta: BamMetadata) -> None:
    """
    Abort if the BAM has no @SQ reference sequences (i.e. is unaligned).
    """
    if not meta.reference_names:
        raise click.ClickException(
            "BAM appears to be unaligned: no @SQ reference sequences found "
            "in header. scNoiseMeter requires a coordinate-sorted, aligned BAM. "
            "Please align your reads first using minimap2, pbmm2, or STAR."
        )


def _validate_sort_order(meta: BamMetadata) -> None:
    """
    Abort if the BAM is not coordinate-sorted.
    """
    so = meta.sort_order.lower()
    if so and so != "coordinate":
        raise click.ClickException(
            f"BAM sort order is '{meta.sort_order}', but scNoiseMeter requires "
            f"coordinate-sorted BAMs.\n"
            f"\n"
            f"Please sort your BAM first:\n"
            f"  samtools sort -o sorted.bam {meta.path.name}\n"
            f"  samtools index sorted.bam"
        )
    if not so:
        logger.warning(
            "BAM @HD SO tag is absent — cannot confirm sort order. "
            "Proceeding, but the BAM must be coordinate-sorted."
        )


def _is_illumina_platform(platform) -> bool:
    """Return True for any Illumina-family platform (includes Smart-seq — paired-end)."""
    illumina_vals = {Platform.ILLUMINA, Platform.ILLUMINA_10X, Platform.ILLUMINA_BD, Platform.SMARTSEQ}
    return platform in illumina_vals


def _apply_intergenic_reclassification(result, intergenic_records, record_categories):
    """
    Move per-barcode read/base counts from INTERGENIC_SPARSE into the promoted
    sub-categories (HOTSPOT/NOVEL/REPEAT) according to the profiler's per-record
    classification.

    ``result.intergenic_reads`` is a reservoir sample (cap 500_000) of the
    classifier's SPARSE output, so per-record reclassification is only
    exhaustive when the full sample fits in the reservoir.  For larger samples
    the reservoir is a uniform sub-sample; the per-barcode proportions
    computed from it are unbiased estimators of the true proportions, and we
    scale each barcode's SPARSE counts accordingly.

    Mutates ``result.read_counts`` and ``result.base_counts`` in place.
    """
    if not intergenic_records or not record_categories:
        return

    # Per-barcode counts of sampled records by promoted category
    per_cb_by_cat: dict[str, dict] = {}
    per_cb_total: dict[str, int] = {}
    for rec, new_cat in zip(intergenic_records, record_categories):
        cb = rec.cell_barcode
        per_cb_total[cb] = per_cb_total.get(cb, 0) + 1
        if new_cat != ReadCategory.INTERGENIC_SPARSE:
            per_cb_by_cat.setdefault(cb, {})
            per_cb_by_cat[cb][new_cat] = per_cb_by_cat[cb].get(new_cat, 0) + 1

    for cb, cat_counts in per_cb_by_cat.items():
        sparse_reads = result.read_counts[cb].get(ReadCategory.INTERGENIC_SPARSE, 0)
        sparse_bases = result.base_counts[cb].get(ReadCategory.INTERGENIC_SPARSE, 0)
        denom = per_cb_total[cb]
        if denom == 0 or sparse_reads == 0:
            continue
        moved_reads = 0
        moved_bases = 0
        for new_cat, n in cat_counts.items():
            frac = n / denom
            r_move = min(int(round(sparse_reads * frac)), sparse_reads - moved_reads)
            b_move = min(int(round(sparse_bases * frac)), sparse_bases - moved_bases)
            if r_move <= 0 and b_move <= 0:
                continue
            result.read_counts[cb][new_cat] += r_move
            result.base_counts[cb][new_cat] += b_move
            moved_reads += r_move
            moved_bases += b_move
        result.read_counts[cb][ReadCategory.INTERGENIC_SPARSE] -= moved_reads
        result.base_counts[cb][ReadCategory.INTERGENIC_SPARSE] -= moved_bases


@cli.command("run")
@click.option("--bam",           required=True, type=click.Path(exists=True), help="Input BAM (must have .bai index).")
@click.option("--sample-name",   default=None, help="Sample label used in output files and reports.  Defaults to BAM filename stem.")
@click.option("--cell-barcodes", default=None, type=click.Path(exists=True),
              help="Path to called-cell barcodes file (one per line, plain text or .gz). "
                   "Reads whose CB tag is not in this list are skipped entirely and do not "
                   "contribute to any metric.  Compatible with Cell Ranger "
                   "filtered_feature_bc_matrix/barcodes.tsv.gz.  Distinct from "
                   "--barcode-whitelist (which defines valid barcode sequences for the chemistry).")
@_shared_options
def run_cmd(
    bam, sample_name, cell_barcodes,
    gtf, gtf_version, barcode_whitelist, whitelist_db, barcode_tag, umi_tag,
    chemistry, platform, pipeline_stage, chimeric_distance,
    repeats, reference, threads, no_umi_dedup, no_cache,
    exclude_biotypes, output_dir, obs_metadata, polya_sites, polya_db,
    tss_sites, tss_db, numt_bed, offline, verbose,
):
    """
    Classify reads in a single BAM and produce QC metrics.

    Output files written to --output-dir:
      <sample>.read_metrics.tsv       Sample-wide scalar metrics
      <sample>.cell_metrics.tsv       Per-cell metrics table
      <sample>.multiqc.json           MultiQC custom content
      <sample>.length_distributions/  Per-category read-length histograms
    """
    _setup_logging(verbose)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    bam_path = Path(bam)
    sample_name = sample_name or bam_path.stem

    # Quick sort-order check from BAM header (does not require an index)
    # Done first so that name-sorted BAMs get a clear "wrong sort order" error
    # rather than a generic "no index" error.
    import pysam as _pysam
    with _pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as _bam_hdr:
        _hd = _bam_hdr.header.to_dict().get("HD", {})
        _so_early = _hd.get("SO", "").lower()
    if _so_early and _so_early != "coordinate":
        raise click.ClickException(
            f"BAM sort order is '{_so_early}', but scNoiseMeter requires "
            f"coordinate-sorted BAMs.\n"
            f"\n"
            f"Please sort your BAM first:\n"
            f"  samtools sort -o sorted.bam {bam_path.name}\n"
            f"  samtools index sorted.bam"
        )

    # Pre-check: BAM index must exist (pysam requires it for random access)
    bai_path = Path(str(bam_path) + ".bai")
    bai_path2 = bam_path.with_suffix(".bai")
    if not bai_path.exists() and not bai_path2.exists():
        raise click.ClickException(
            f"BAM index not found for '{bam_path.name}'.\n"
            f"\n"
            f"Please create an index first:\n"
            f"  samtools index {bam_path.name}\n"
            f"\n"
            f"This produces '{bam_path.name}.bai', which is required for\n"
            f"chromosome-level parallel processing."
        )

    # Load whitelist (explicit path takes priority over auto-download)
    whitelist = _resolve_whitelist(barcode_whitelist, whitelist_db, offline)

    # Resolve platform / stage overrides
    plat  = None if platform == "auto"       else Platform(platform)
    stage = None if pipeline_stage == "auto" else PipelineStage(pipeline_stage)

    # Inspect BAM
    meta = inspect_bam(
        bam_path,
        barcode_tag=barcode_tag,
        umi_tag=umi_tag,
        platform=plat,
        pipeline_stage=stage,
    )

    # Validate sort order (fatal if not coordinate-sorted; also warns if SO absent)
    _validate_sort_order(meta)

    # Validate that the BAM is aligned (has @SQ reference sequences)
    _validate_sq_lines(meta)

    # Warn if post_filter stage is declared but no corrected barcodes are found.
    # Smart-seq is excluded: one BAM per cell, no CB tags is expected.
    if (
        meta.pipeline_stage == PipelineStage.POST_FILTER
        and not meta.barcode_aware
        and stage == PipelineStage.POST_FILTER  # user explicitly requested post_filter
        and meta.platform != Platform.SMARTSEQ
    ):
        _cb_msg = (
            f"CB tags are absent despite --pipeline-stage post_filter. "
            f"Only {meta.barcode_fraction:.1%} of sampled reads carry the '{barcode_tag}' tag "
            f"(threshold: {BARCODE_AUTODETECT_MIN_FRACTION:.0%}). "
            f"Automatically switching to barcode-agnostic mode — all reads will be aggregated "
            f"as a single sample. Per-cell metrics will not be available. "
            f"If this is unexpected, verify that your BAM was produced after barcode correction."
        )
        logger.warning(_cb_msg)
        click.echo(f"\nWarning: {_cb_msg}", err=True)

    # Load and validate --cell-barcodes
    cell_barcode_set: Optional[set] = None
    if cell_barcodes is not None:
        if not meta.barcode_aware:
            _cb_ignored_msg = (
                "--cell-barcodes was supplied but this BAM is running in barcode-agnostic mode "
                f"(no CB tags found; only {meta.barcode_fraction:.1%} of sampled reads carry "
                f"the '{barcode_tag}' tag).  The --cell-barcodes filter will be ignored."
            )
            logger.warning(_cb_ignored_msg)
            click.echo(f"\nWarning: {_cb_ignored_msg}", err=True)
        else:
            cell_barcode_set = _load_cell_barcodes(cell_barcodes)

    # Pre-flight BAM read count check (fast — reads from .bai index counters)
    import pysam as _pysam
    with _pysam.AlignmentFile(str(bam_path), "rb") as _bam_check:
        total_mapped   = _bam_check.mapped
        total_unmapped = _bam_check.unmapped
    total_reads = total_mapped + total_unmapped

    if total_reads == 0:
        raise click.ClickException(
            f"BAM file contains no reads.\n"
            f"\n"
            f"'{bam_path.name}' appears to be completely empty. "
            f"Please verify the input file."
        )

    if total_mapped == 0:
        click.echo(
            f"\nWarning: No mapped reads found in BAM '{bam_path.name}'.\n"
            f"All {total_unmapped:,} reads are unmapped — there is nothing to classify.\n"
            f"Exiting without producing a report.",
            err=True,
        )
        raise SystemExit(0)

    # Resolve GTF: use cache or auto-download if --gtf was not supplied
    gtf, gtf_version, gtf_source = _resolve_gtf(gtf, gtf_version=gtf_version, offline=offline)

    # Validate chromosome naming vs GTF (fatal if mismatch)
    _validate_chromosome_naming(meta, gtf)

    # Validate chromosome lengths vs reference FASTA (warn if wrong species/assembly)
    _validate_chromosome_lengths(meta, reference_path=reference)

    # Resolve polyA sites: use cache or auto-download if --polya-sites omitted
    polya_paths, polya_version, polya_source = _resolve_polya_sites(
        polya_sites,
        polya_db=polya_db,
        hint_gencode_version=gtf_version if gtf_source == "auto-downloaded" else None,
        offline=offline,
    )

    # Version consistency check (warn if GTF and atlas differ by > 5 major versions)
    version_warning = _check_version_consistency(gtf_version, polya_version)

    # Build annotation index
    index = build_annotation_index(
        gtf,
        repeats_path=repeats,
        exclude_biotypes=list(exclude_biotypes),
        cache=not no_cache,
    )

    # Enable paired-end chimeric detection for Illumina platforms
    use_paired_chimeric = _is_illumina_platform(meta.platform)
    if use_paired_chimeric:
        logger.info("Illumina platform detected — enabling paired-end chimeric detection.")

    # Run pipeline
    result = run_pipeline(
        bam_path, index, meta,
        whitelist=whitelist,
        cell_barcodes=cell_barcode_set,
        chimeric_distance=chimeric_distance,
        paired_end_chimeric=use_paired_chimeric,
        threads=threads,
        store_umis=not no_umi_dedup,
        reference_path=reference,
    )

    # Check that --cell-barcodes didn't filter out every read
    if cell_barcode_set is not None and result.n_reads_skipped_not_called_cell > 0:
        reads_remaining = result.n_reads_processed
        if reads_remaining == 0:
            raise click.ClickException(
                f"0 reads remain after cell barcode filtering — "
                f"check that the barcode file matches this BAM.\n"
                f"\n"
                f"All {result.n_reads_skipped_not_called_cell:,} reads were skipped because "
                f"their CB tag did not appear in '{Path(cell_barcodes).name}'.\n"
                f"Verify that the barcodes file was produced from the same sample/run."
            )

    # Warn if total read count is very low (statistically unreliable)
    MIN_RELIABLE_READS = 100
    if result.n_reads_processed < MIN_RELIABLE_READS:
        _low_count_msg = (
            f"Only {result.n_reads_processed} reads were classified "
            f"(minimum recommended: {MIN_RELIABLE_READS}). "
            f"Results are statistically unreliable and should not be interpreted."
        )
        logger.warning(_low_count_msg)
        click.echo(f"\nWarning: {_low_count_msg}", err=True)
        # Also store in meta.warnings so it surfaces in the HTML report
        meta.warnings.append(_low_count_msg)

    # Detect BAM chromosome naming style for BED normalization
    _bam_chrom_style = _detect_chrom_style(meta.reference_names)

    # Attach polyA site dict (merged from resolved polya-site files)
    result._polya_site_dict = _load_polya_sites(polya_paths, chrom_style=_bam_chrom_style)

    # Attach TSS site dict (resolved: auto-download or user-supplied)
    tss_paths = _resolve_tss_sites(tss_sites, tss_db=tss_db, offline=offline)
    result._tss_site_dict = _load_tss_sites(tss_paths, chrom_style=_bam_chrom_style) if tss_paths else None

    # Attach NUMT intervals
    if numt_bed:
        result._numt_intervals = _load_numt_bed(numt_bed)
    else:
        result._numt_intervals = None

    # --- Intergenic profiler (must run BEFORE compute_metrics so promoted
    # loci contribute to the correct noise/ambiguous categories) ---
    intergenic_records = extract_intergenic_records(result)
    intergenic_loci    = []
    if intergenic_records:
        total_ig_bases  = compute_intergenic_bases(index)
        total_barcodes  = len(result.read_counts)

        polya_site_dict = result._polya_site_dict

        repeat_dict = None
        if index.repeats is not None and not index.repeats.df.empty:
            repeat_dict = {}
            for _, row in index.repeats.df.iterrows():
                c = row["Chromosome"]
                repeat_dict.setdefault(c, []).append(
                    (int(row["Start"]), int(row["End"]))
                )

        ref_handle = None
        if reference:
            try:
                import pysam as _pysam
                ref_handle = _pysam.FastaFile(reference)
            except Exception as _exc:
                logger.warning("Could not open reference FASTA for polyA context: %s", _exc)
        try:
            ig_loci, record_categories = profile_intergenic_loci(
                intergenic_records,
                total_intergenic_bases=total_ig_bases,
                total_barcodes=total_barcodes,
                reference=ref_handle,
                polya_sites=polya_site_dict,
                repeat_intervals=repeat_dict,
            )
        finally:
            if ref_handle is not None:
                ref_handle.close()
        intergenic_loci = ig_loci
        _apply_intergenic_reclassification(result, intergenic_records, record_categories)
        logger.info(
            "Intergenic profiling complete: %d loci characterised.", len(ig_loci)
        )

    # Compute metrics (after intergenic reclassification)
    platform_str = meta.platform.value
    sm, ct = compute_metrics(
        result, sample_name,
        platform=platform_str,
        unstranded=(meta.platform == Platform.SMARTSEQ),
    )

    # Attach annotation provenance for the HTML report
    sm._gtf_info = {"version": gtf_version, "source": gtf_source, "path": gtf}
    sm._polya_info = {
        "version": polya_version,
        "source": polya_source,
        "path": polya_paths[0] if polya_paths else None,
        "db": polya_db,
    }
    sm._tss_info = {
        "db": tss_db,
        "source": "user-supplied" if tss_sites else "auto-downloaded",
        "path": tss_paths[0] if tss_paths else None,
    } if tss_paths else None
    if version_warning:
        sm.warnings.append(version_warning)

    # Attach cell barcode info for the HTML report
    if cell_barcode_set is not None:
        sm._cell_barcodes_info = {
            "path": cell_barcodes,
            "n_barcodes": len(cell_barcode_set),
            "n_reads_skipped": result.n_reads_skipped_not_called_cell,
        }
    else:
        sm._cell_barcodes_info = None

    # --- Per-cluster metrics ---
    cluster_df = None
    if obs_metadata:
        try:
            obs_df = load_obs_metadata(obs_metadata)
            cluster_df = compute_cluster_metrics(ct, obs_df)
            logger.info("Per-cluster metrics computed for %d clusters.", len(cluster_df))
        except Exception as exc:
            logger.warning("Could not compute cluster metrics: %s", exc)

    # Write outputs
    _write_run_outputs(sm, ct, result, output_dir, sample_name,
                       cluster_df=cluster_df, intergenic_loci=intergenic_loci,
                       platform=meta.platform)

    click.echo(f"\n✓ scNoiseMeter run complete.  Results in: {output_dir}")
    click.echo(f"  Noise fraction (reads) : {sm.noise_read_frac:.2%}")
    click.echo(f"  Noise fraction (bases) : {sm.noise_base_frac:.2%}")
    click.echo(f"  Strand concordance     : {sm.strand_concordance:.2%}")
    click.echo(f"  Chimeric rate          : {sm.chimeric_read_frac:.2%}")
    click.echo(f"  Cells detected         : {sm.n_cells:,}")
    if cell_barcode_set is not None:
        click.echo(f"  Cell barcodes file     : {cell_barcodes} ({len(cell_barcode_set):,} barcodes)")
        click.echo(f"  Reads skipped (not called cell) : {result.n_reads_skipped_not_called_cell:,}")


# ---------------------------------------------------------------------------
# `compare` subcommand
# ---------------------------------------------------------------------------

@cli.command("compare")
@click.option("--bam-a",        required=True, type=click.Path(exists=True), help="BAM A (e.g. pre-filter / raw).")
@click.option("--bam-b",        required=True, type=click.Path(exists=True), help="BAM B (e.g. post-filter).")
@click.option("--label-a",      default="sample_A", show_default=True, help="Label for BAM A in reports.")
@click.option("--label-b",      default="sample_B", show_default=True, help="Label for BAM B in reports.")
@_shared_options
def compare_cmd(
    bam_a, bam_b, label_a, label_b,
    gtf, gtf_version, barcode_whitelist, whitelist_db, barcode_tag, umi_tag,
    chemistry, platform, pipeline_stage, chimeric_distance,
    repeats, reference, threads, no_umi_dedup, no_cache,
    exclude_biotypes, output_dir, obs_metadata, polya_sites, polya_db,
    tss_sites, tss_db, numt_bed, offline, verbose,
):
    """
    Compare noise profiles between two BAMs (e.g. pre- vs post-filter).

    Runs the full classification pipeline on both BAMs, computes paired
    proportion tests for each noise category, and writes a comparison report.

    Output files written to --output-dir:
      comparison.metrics.tsv          Side-by-side scalar metrics
      comparison.stats.tsv            Paired proportion test results
      comparison.report.html          Interactive HTML report
    """
    _setup_logging(verbose)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    whitelist = _resolve_whitelist(barcode_whitelist, whitelist_db, offline)
    plat  = None if platform == "auto"       else Platform(platform)
    stage = None if pipeline_stage == "auto" else PipelineStage(pipeline_stage)

    # Resolve GTF (use cache or auto-download if not supplied)
    gtf, gtf_version, gtf_source = _resolve_gtf(gtf, gtf_version=gtf_version, offline=offline)

    # Resolve polyA sites (use cache or auto-download if not supplied)
    polya_paths, polya_version, polya_source = _resolve_polya_sites(
        polya_sites,
        polya_db=polya_db,
        hint_gencode_version=gtf_version if gtf_source == "auto-downloaded" else None,
        offline=offline,
    )

    # Resolve TSS sites
    tss_paths = _resolve_tss_sites(tss_sites, tss_db=tss_db, offline=offline)

    # Version consistency check
    version_warning = _check_version_consistency(gtf_version, polya_version)

    # Build annotation index once (shared between both runs)
    index = build_annotation_index(
        gtf,
        repeats_path=repeats,
        exclude_biotypes=list(exclude_biotypes),
        cache=not no_cache,
    )

    results = {}
    for bam_path, label in [(bam_a, label_a), (bam_b, label_b)]:
        click.echo(f"\nProcessing {label} ({bam_path}) …")
        meta = inspect_bam(Path(bam_path), barcode_tag=barcode_tag, umi_tag=umi_tag,
                           platform=plat, pipeline_stage=stage)
        result = run_pipeline(
            bam_path, index, meta,
            whitelist=whitelist,
            chimeric_distance=chimeric_distance,
            threads=threads,
            store_umis=not no_umi_dedup,
            reference_path=reference,
        )
        _bam_cs = _detect_chrom_style(meta.reference_names)
        result._polya_site_dict = _load_polya_sites(polya_paths, chrom_style=_bam_cs)
        result._tss_site_dict = _load_tss_sites(tss_paths, chrom_style=_bam_cs) if tss_paths else None
        sm, ct = compute_metrics(
            result, label,
            platform=meta.platform.value,
            unstranded=(meta.platform == Platform.SMARTSEQ),
        )
        sm._gtf_info = {"version": gtf_version, "source": gtf_source, "path": gtf}
        sm._polya_info = {
            "version": polya_version,
            "source": polya_source,
            "path": polya_paths[0] if polya_paths else None,
            "db": polya_db,
        }
        sm._tss_info = {
            "db": tss_db,
            "source": "user-supplied" if tss_sites else "auto-downloaded",
            "path": tss_paths[0] if tss_paths else None,
        } if tss_paths else None
        if version_warning:
            sm.warnings.append(version_warning)
        results[label] = (sm, ct, result)

    # Write comparison outputs
    _write_compare_outputs(results, output_dir, label_a, label_b)

    click.echo(f"\n✓ Comparison complete.  Results in: {output_dir}")
    sm_a_obj, _, _ = results[label_a]
    sm_b_obj, _, _ = results[label_b]
    click.echo(f"  {label_a} noise fraction (reads): {sm_a_obj.noise_read_frac:.2%}")
    click.echo(f"  {label_b} noise fraction (reads): {sm_b_obj.noise_read_frac:.2%}")


# ---------------------------------------------------------------------------
# `discover` subcommand
# ---------------------------------------------------------------------------

@cli.command("discover")
@click.option("--bam-dir",     required=True,  type=click.Path(exists=True, file_okay=False),
              help="Directory to scan for .bam files.")
@click.option("--reference",   required=True,  type=click.Path(exists=True),
              help="Reference FASTA (.fa/.fa.gz + .fai) for polyA and junction checks.")
@click.option("--tss-sites",   multiple=True,  type=click.Path(exists=True),
              help="CAGE peak / TSS BED file(s) for the 5\u2032-anchored metric.")
@click.option("--threads",     default=DEFAULT_THREADS, show_default=True, type=int,
              help="Parallel worker processes.")
@click.option("--output-dir",  required=True,  type=click.Path(),
              help="Root output directory; each BAM gets its own subdirectory.")
@click.option("--gtf",         default=None,   type=click.Path(exists=True),
              help="GENCODE GTF (plain or .gz). Auto-downloaded if omitted. Takes precedence over --gtf-version.")
@click.option("--gtf-version", default=None,   type=int,
              help="GENCODE release number to auto-download (e.g. 42). Ignored when --gtf is supplied.")
@click.option("--polya-sites", multiple=True,  type=click.Path(exists=True),
              help="PolyA site BED file(s). Auto-downloaded if omitted; see --polya-db.")
@click.option("--polya-db",
              default="polyasite3", show_default=True,
              type=click.Choice(["polyasite3", "polyadb4", "both"], case_sensitive=False),
              help="polyA database to auto-download when --polya-sites is not supplied. "
                   "polyasite3 / polyadb4 / both (merge).")
@click.option("--tss-sites",   multiple=True,  type=click.Path(exists=True),
              help="CAGE peak / TSS BED file(s). Auto-downloaded if omitted; see --tss-db.")
@click.option("--tss-db",
              default="fantom5", show_default=True,
              type=click.Choice(["fantom5", "none"], case_sensitive=False),
              help="TSS database to auto-download when --tss-sites is not supplied. "
                   "fantom5 / none (disable TSS metric).")
@click.option("--run-all",     is_flag=True,
              help="Non-interactive: auto-run all BAMs that can be fully inferred. "
                   "BAMs with blocking issues are skipped with a warning.")
@click.option("--offline",     is_flag=True,
              help="Use only cached annotation files; never make network calls. "
                   "Raises an error if the cache is empty. "
                   "Ignored when --gtf / --polya-sites are supplied explicitly.")
@click.option("--verbose", "-v", is_flag=True, help="Enable debug logging.")
def discover_cmd(
    bam_dir, reference, tss_sites, tss_db, threads, output_dir,
    gtf, gtf_version, polya_sites, polya_db, run_all, offline, verbose,
):
    """
    Discover BAMs in a directory, infer their parameters, and run scNoiseMeter
    on the selected (or all) files.

    In interactive mode (default), a summary table is shown and the user
    selects which BAMs to process.  In --run-all mode, all BAMs with fully
    inferable parameters are run automatically.
    """
    from scnoisemeter.utils.discover_inspector import (
        DiscoverBamInfo, inspect_bam_for_discover,
        format_discovery_table, _collect_selected_indices,
        _prompt_platform, _prompt_chemistry, _normalise_platform,
    )
    from scnoisemeter.modules.pipeline import run_pipeline
    from scnoisemeter.modules.metrics import compute_metrics

    _setup_logging(verbose)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    bam_dir = Path(bam_dir)

    # Resolve shared annotation references upfront (so we only download once)
    click.echo("\nResolving annotation references …")
    gtf_path, gtf_version, gtf_source = _resolve_gtf(gtf, gtf_version=gtf_version, offline=offline)
    polya_paths, polya_version, polya_source = _resolve_polya_sites(
        polya_sites,
        polya_db=polya_db,
        hint_gencode_version=gtf_version if gtf_source == "auto-downloaded" else None,
        offline=offline,
    )
    tss_paths = _resolve_tss_sites(tss_sites, tss_db=tss_db, offline=offline)
    version_warning = _check_version_consistency(gtf_version, polya_version)

    # ------------------------------------------------------------------ #
    # STEP 1 — Discover and inspect all BAMs
    # ------------------------------------------------------------------ #
    bam_files = sorted(bam_dir.glob("*.bam"))
    if not bam_files:
        click.echo(f"No .bam files found in {bam_dir}", err=True)
        return

    click.echo(f"\nInspecting {len(bam_files)} BAM file(s) in {bam_dir} …\n")
    infos: list[DiscoverBamInfo] = []
    for bam_path in bam_files:
        click.echo(f"  {bam_path.name} …", nl=False)
        info = inspect_bam_for_discover(bam_path)
        infos.append(info)
        status = "ok" if info.can_run else "⚠ " + "; ".join(info.run_issues[:1])
        click.echo(f" {status}")

    # ------------------------------------------------------------------ #
    # STEP 2 — Summary table and selection
    # ------------------------------------------------------------------ #
    click.echo("\n" + format_discovery_table(infos))

    if run_all:
        # Non-interactive: take all that have no blocking issues
        selected_indices = [i for i, info in enumerate(infos) if info.can_run]
        skipped = [infos[i] for i in range(len(infos)) if i not in selected_indices]
        for s in skipped:
            click.echo(
                f"  Skipping {s.bam_path.name}: {'; '.join(s.run_issues)}", err=True
            )
    else:
        # Interactive selection
        selected_indices = _collect_selected_indices(infos)
        if not selected_indices:
            click.echo("No BAMs selected — exiting.")
            return

        # Resolve unknown platform/chemistry interactively
        for idx in list(selected_indices):
            info = infos[idx]
            if not info.has_index or (
                info.sort_order and info.sort_order.lower() not in ("coordinate", "")
            ):
                click.echo(
                    f"  Skipping {info.bam_path.name}: {'; '.join(info.run_issues)}",
                    err=True,
                )
                selected_indices.remove(idx)
                continue

            if info.platform == Platform.UNKNOWN:
                platform_val = _prompt_platform(info)
                if platform_val is None:
                    selected_indices.remove(idx)
                    continue
                info.meta.platform = Platform(platform_val)
                info.meta.platform_confidence = "user"
                # Remove the blocking "platform unknown" issue now that it's resolved
                info.run_issues = [r for r in info.run_issues if "platform" not in r]

    if not selected_indices:
        click.echo("No BAMs to run.")
        return

    # ------------------------------------------------------------------ #
    # STEP 3 — Run each selected BAM
    # ------------------------------------------------------------------ #
    results_summary = []

    for idx in selected_indices:
        info = infos[idx]
        stem = info.bam_path.stem
        bam_output_dir = output_dir / stem
        bam_output_dir.mkdir(parents=True, exist_ok=True)

        # Resolve chemistry and platform for this BAM
        # Normalise illumina_10x/illumina_bd → illumina for --platform flag
        platform_val = _normalise_platform(info.platform) if info.platform != Platform.UNKNOWN else "auto"
        stage_val = info.pipeline_stage.value

        # Construct and print the exact command
        cmd_parts = [
            "scnoisemeter", "run",
            "--bam", str(info.bam_path),
            "--reference", reference,
            "--gtf", gtf_path,
        ]
        for p in polya_paths:
            cmd_parts += ["--polya-sites", p]
        for ts in tss_paths:
            cmd_parts += ["--tss-sites", ts]
        cmd_parts += [
            "--platform", platform_val,
            "--pipeline-stage", stage_val,
            "--threads", str(threads),
            "--output-dir", str(bam_output_dir),
            "--sample-name", stem,
        ]
        if verbose:
            cmd_parts.append("--verbose")

        click.echo(f"\n{'='*60}")
        click.echo(f"Running: {stem}")
        click.echo("  " + " ".join(cmd_parts))
        click.echo("=" * 60)

        # Run inline (reuse the resolved annotation + pipeline directly)
        try:
            _run_single_bam_for_discover(
                info=info,
                gtf_path=gtf_path,
                gtf_version=gtf_version,
                gtf_source=gtf_source,
                polya_paths=polya_paths,
                polya_version=polya_version,
                polya_source=polya_source,
                polya_db=polya_db,
                version_warning=version_warning,
                tss_paths=tss_paths,
                tss_db=tss_db,
                reference=reference,
                threads=threads,
                output_dir=bam_output_dir,
                sample_name=stem,
                verbose=verbose,
            )
            results_summary.append((stem, bam_output_dir, "success", None))
            click.echo(f"  ✓ Completed: {bam_output_dir}")
        except Exception as exc:
            logger.error("Run failed for %s: %s", stem, exc, exc_info=verbose)
            results_summary.append((stem, bam_output_dir, "failed", str(exc)))
            click.echo(f"  ✗ Failed: {exc}", err=True)

    # ------------------------------------------------------------------ #
    # STEP 4 — Summary table
    # ------------------------------------------------------------------ #
    click.echo(f"\n{'='*60}")
    click.echo("DISCOVER SUMMARY")
    click.echo("=" * 60)
    _print_discover_summary(results_summary, output_dir)


def _run_single_bam_for_discover(
    *,
    info: "DiscoverBamInfo",
    gtf_path: str,
    gtf_version,
    gtf_source: str,
    polya_paths: list,
    polya_version,
    polya_source: str,
    polya_db: str,
    version_warning: Optional[str],
    tss_paths: list,
    tss_db: str,
    reference: str,
    threads: int,
    output_dir: Path,
    sample_name: str,
    verbose: bool,
) -> None:
    """Run the full scNoiseMeter pipeline on one BAM (used by discover_cmd)."""
    from scnoisemeter.modules.annotation import build_annotation_index
    from scnoisemeter.modules.intergenic_profiler import (
        profile_intergenic_loci, extract_intergenic_records, compute_intergenic_bases,
    )
    from scnoisemeter.modules.metrics import compute_metrics, compute_length_stratification
    from scnoisemeter.modules.pipeline import run_pipeline

    meta = info.meta
    bam_path = info.bam_path

    # Validate sort order and chromosome naming
    _validate_sort_order(meta)
    _validate_chromosome_naming(meta, gtf_path)

    # Build annotation index
    index = build_annotation_index(
        gtf_path,
        exclude_biotypes=[],
        cache=True,
    )

    use_paired_chimeric = _is_illumina_platform(meta.platform)

    result = run_pipeline(
        bam_path, index, meta,
        whitelist=None,
        chimeric_distance=DEFAULT_CHIMERIC_DISTANCE,
        paired_end_chimeric=use_paired_chimeric,
        threads=threads,
        store_umis=True,
        reference_path=reference,
    )

    _bam_cs = _detect_chrom_style(meta.reference_names)
    result._polya_site_dict = _load_polya_sites(polya_paths, chrom_style=_bam_cs)

    result._tss_site_dict = _load_tss_sites(tss_paths, chrom_style=_bam_cs) if tss_paths else None

    result._numt_intervals = None

    # Intergenic profiler (runs before compute_metrics so promoted loci
    # contribute to the correct categories)
    intergenic_records = extract_intergenic_records(result)
    intergenic_loci = []
    if intergenic_records:
        total_ig_bases = compute_intergenic_bases(index)
        ref_handle = None
        if reference:
            try:
                import pysam as _pysam
                ref_handle = _pysam.FastaFile(reference)
            except Exception as _exc:
                logger.warning("Could not open reference FASTA for polyA context: %s", _exc)
        try:
            ig_loci, record_categories = profile_intergenic_loci(
                intergenic_records,
                total_intergenic_bases=total_ig_bases,
                total_barcodes=len(result.read_counts),
                reference=ref_handle,
                polya_sites=result._polya_site_dict,
            )
        finally:
            if ref_handle is not None:
                ref_handle.close()
        intergenic_loci = ig_loci
        _apply_intergenic_reclassification(result, intergenic_records, record_categories)

    sm, ct = compute_metrics(
        result, sample_name,
        platform=meta.platform.value,
        unstranded=(meta.platform == Platform.SMARTSEQ),
    )
    sm._gtf_info = {"version": gtf_version, "source": gtf_source, "path": gtf_path}
    sm._polya_info = {
        "version": polya_version,
        "source": polya_source,
        "path": polya_paths[0] if polya_paths else None,
        "db": polya_db,
    }
    sm._tss_info = {
        "db": tss_db,
        "source": "auto-downloaded",
        "path": tss_paths[0] if tss_paths else None,
    } if tss_paths else None
    if version_warning:
        sm.warnings.append(version_warning)

    strat_df = compute_length_stratification(result.length_bin_counts, result.length_samples)
    _write_run_outputs(
        sm, ct, result, output_dir, sample_name,
        cluster_df=None,
        intergenic_loci=intergenic_loci,
        platform=meta.platform,
    )

    click.echo(f"  Noise conservative : {sm.noise_read_frac:.2%}")
    click.echo(f"  Chimeric rate      : {sm.chimeric_read_frac:.2%}")
    click.echo(f"  Strand concordance : {sm.strand_concordance:.2%}")


def _print_discover_summary(
    results: list,
    output_dir: Path,
) -> None:
    """Print a results summary table from discover runs."""
    if not results:
        click.echo("  (no runs completed)")
        return

    headers = ("Sample", "Status", "Output dir")
    rows = []
    for stem, bam_out, status, err in results:
        rows.append((stem, status, str(bam_out)))

        # Try to read metrics from TSV
        metrics_tsv = bam_out / f"{stem}.read_metrics.tsv"
        if metrics_tsv.exists():
            metrics = {}
            try:
                with open(metrics_tsv) as fh:
                    for line in fh:
                        k, _, v = line.strip().partition("\t")
                        metrics[k] = v
            except Exception:
                pass
            if metrics:
                click.echo(
                    f"\n  {stem}:\n"
                    f"    noise_conservative  = {float(metrics.get('noise_read_frac', 0)):.2%}\n"
                    f"    noise_strict        = {float(metrics.get('noise_read_frac_strict', metrics.get('noise_read_frac', 0))):.2%}\n"
                    f"    chimeric_rate       = {float(metrics.get('chimeric_read_frac', 0)):.2%}\n"
                    f"    strand_concordance  = {float(metrics.get('strand_concordance', 0)):.2%}"
                )
        elif err:
            click.echo(f"\n  {stem}: FAILED — {err}")


# ---------------------------------------------------------------------------
# `run-plate` subcommand
# ---------------------------------------------------------------------------

import re as _re

_PLATE_WELL_RE = _re.compile(r"^(?P<plate>.+)_(?P<well>[A-Pa-p]\d{1,2})$")

# ---------------------------------------------------------------------------
# Parallel well processing — module-level so ProcessPoolExecutor can pickle
# ---------------------------------------------------------------------------

_worker_state: dict = {}


def _plate_worker_init(gtf, repeats_path, exclude_biotypes,
                       polya_paths, tss_paths, chrom_style):
    """Initialiser for each worker process: loads shared read-only data once."""
    import logging as _log
    _log.disable(_log.WARNING)  # suppress per-well INFO/WARNING spam in workers

    from scnoisemeter.modules.annotation import build_annotation_index
    _worker_state["index"] = build_annotation_index(
        gtf,
        repeats_path=repeats_path,
        exclude_biotypes=list(exclude_biotypes),
        cache=True,
    )
    _worker_state["polya"] = (
        _load_polya_sites(polya_paths, chrom_style=chrom_style) if polya_paths else {}
    )
    _worker_state["tss"] = (
        _load_tss_sites(tss_paths, chrom_style=chrom_style) if tss_paths else None
    )


def _plate_well_task(task: dict) -> dict:
    """Worker task: inspect BAM, run pipeline, compute per-well metrics."""
    from scnoisemeter.constants import Platform
    from scnoisemeter.modules.metrics import compute_metrics
    from scnoisemeter.modules.pipeline import run_pipeline
    from scnoisemeter.utils.bam_inspector import inspect_bam

    bam_path = Path(task["bam_path"])
    well_id  = task["well_id"]
    plate_id = task["plate_id"]

    try:
        meta = inspect_bam(
            bam_path,
            barcode_tag=task["barcode_tag"],
            umi_tag=task["umi_tag"],
            platform=task["platform_override"],
            pipeline_stage=task["stage_override"],
        )
        _validate_sort_order(meta)

        well_result = run_pipeline(
            bam_path, _worker_state["index"], meta,
            whitelist=task["whitelist"],
            chimeric_distance=task["chimeric_distance"],
            paired_end_chimeric=_is_illumina_platform(meta.platform),
            threads=task["threads"],
            store_umis=task["store_umis"],
            reference_path=task["reference"],
        )

        well_result._polya_site_dict = _worker_state["polya"]
        well_result._tss_site_dict   = _worker_state["tss"]
        well_result._numt_intervals  = None

        well_sm, _ = compute_metrics(
            well_result, f"{plate_id}_{well_id}",
            platform=meta.platform.value,
            unstranded=(meta.platform == Platform.SMARTSEQ),
        )

        # Clear large dicts before pickling result back to the main process
        well_result._polya_site_dict = None
        well_result._tss_site_dict   = None

        return {
            "ok": True,
            "well_id":     well_id,
            "well_result": well_result,
            "well_sm":     well_sm,
            "platform":    meta.platform,
        }
    except Exception as exc:
        return {"ok": False, "well_id": well_id, "error": str(exc)}


def _discover_plate_wells(plate_dir: Path) -> dict:
    """
    Scan *plate_dir* for subdirectories whose names match PlateID_WellID.

    Returns
    -------
    dict mapping plate_id → list of (well_id, bam_path) tuples, sorted by
    well_id.  Only subdirectories that contain exactly one BAM file are
    included.  Subdirectories that do not match the pattern are reported as
    warnings but do not cause an error.
    """
    plates: dict = {}
    skipped: list = []

    for subdir in sorted(plate_dir.iterdir()):
        if not subdir.is_dir():
            continue
        m = _PLATE_WELL_RE.match(subdir.name)
        if not m:
            skipped.append(subdir.name)
            continue

        plate_id = m.group("plate")
        well_id  = m.group("well").upper()

        # Find the BAM — accept any *.bam inside the subdirectory
        bams = sorted(subdir.glob("*.bam"))
        if len(bams) == 0:
            logger.warning("No BAM found in %s — skipping well.", subdir)
            continue
        if len(bams) > 1:
            logger.warning(
                "Multiple BAMs in %s — using the first: %s", subdir, bams[0].name
            )
        bam_path = bams[0]

        plates.setdefault(plate_id, []).append((well_id, bam_path))

    if skipped:
        click.echo(
            f"\nNote: {len(skipped)} subdirectories in the plate directory do not match "
            f"the expected PlateID_WellID format (e.g. 881_A1) and were skipped:\n"
            f"  {', '.join(skipped[:10])}"
            + (" …" if len(skipped) > 10 else ""),
            err=True,
        )

    for plate_id in plates:
        plates[plate_id].sort(key=lambda t: t[0])

    return plates


@cli.command("run-plate")
@click.option(
    "--plate-dir", "plate_dir", required=True, type=click.Path(exists=True, file_okay=False),
    help=(
        "Directory containing per-well subdirectories named PlateID_WellID "
        "(e.g. 881_A1, 881_B3, 882_H12). "
        "Each subdirectory must contain exactly one BAM file with a .bai index. "
        "Wells not matching this pattern are skipped with a warning."
    ),
)
@click.option(
    "--sample-sheet", "sample_sheet", default=None, type=click.Path(exists=True),
    help=(
        "CSV file with per-well barcode metadata. "
        "Column detection is automatic: any header containing 'i7' or 'i5' is used as "
        "the barcode column; the key column can be 'Sample_Name' (PlateID_WellID values) "
        "or separate 'WellID'/'PlateID' columns. Headerless files are accepted — columns "
        "are assigned as Sample_Name, i7_sequence, i5_reverse_complement by position. "
        "Use --sequencer to control i7 RC orientation. "
        "Extra columns are carried through to the per-well metrics TSV."
    ),
)
@click.option(
    "--sequencer",
    type=click.Choice(["nextseq", "novaseq", "novaseqx", "aviti", "none"], case_sensitive=False),
    default="none",
    show_default=True,
    help=(
        "Sequencer platform, used to determine i7 orientation in the sample sheet. "
        "NextSeq 500/550/1000/2000 and NovaSeq X report the i7 as its reverse complement "
        "in FASTQ files. NovaSeq 6000 and Element AVITI report i7 in the designed orientation. "
        "When nextseq or novaseqx is passed, the i7 sequence from --sample-sheet is "
        "reverse-complemented before display. Has no effect without --sample-sheet."
    ),
)
@click.option("--sample-name", default=None,
              help="Override for the plate label in output filenames. Defaults to PlateID.")
@click.option(
    "--plate-id", "plate_ids", multiple=True,
    help=(
        "Process only this plate ID (repeatable). "
        "When omitted all discovered plates are processed. "
        "Example: --plate-id 881 --plate-id 882"
    ),
)
@click.option(
    "--parallel-wells", "parallel_wells", default=1, show_default=True, type=int,
    help=(
        "Number of wells to process in parallel using separate worker processes. "
        "Each worker loads the annotation index from cache and runs independently. "
        "Per-worker thread count is --threads // --parallel-wells (minimum 1). "
        "Recommended: set to the number of available CPU cores divided by 4."
    ),
)
@_shared_options
def run_plate_cmd(
    plate_dir, sample_sheet, sequencer, sample_name, plate_ids, parallel_wells,
    gtf, gtf_version, barcode_whitelist, whitelist_db, barcode_tag, umi_tag,
    chemistry, platform, pipeline_stage, chimeric_distance,
    repeats, reference, threads, no_umi_dedup, no_cache,
    exclude_biotypes, output_dir, obs_metadata, polya_sites, polya_db,
    tss_sites, tss_db, numt_bed, offline, verbose,
):
    """
    Classify reads from a plate of Smart-seq wells and produce per-plate metrics.

    Expects a directory of per-well subdirectories named PlateID_WellID
    (e.g. 881_A1). Reads from all wells on the same plate are pooled and
    reported together. Per-well statistics are written when --sample-sheet
    is supplied.

    Output files written to --output-dir (one set per plate):
      <PlateID>.read_metrics.tsv        Plate-wide aggregate metrics
      <PlateID>.cell_metrics.tsv        Placeholder (one row = whole plate)
      <PlateID>.per_well_metrics.tsv    Per-well summary table
      <PlateID>.report.html             Interactive HTML report
    """
    from scnoisemeter.modules.annotation import build_annotation_index
    from scnoisemeter.modules.intergenic_profiler import (
        profile_intergenic_loci, extract_intergenic_records, compute_intergenic_bases,
    )
    from scnoisemeter.modules.metrics import compute_metrics, compute_length_stratification
    from scnoisemeter.modules.pipeline import run_pipeline, merge_sample_results
    from scnoisemeter.utils.sample_sheet import parse_sample_sheet, lookup_well

    _setup_logging(verbose)

    plate_dir_path = Path(plate_dir)
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Discover wells
    # ------------------------------------------------------------------ #
    plates = _discover_plate_wells(plate_dir_path)
    if not plates:
        raise click.ClickException(
            f"No PlateID_WellID subdirectories found in '{plate_dir}'.\n"
            f"\n"
            f"Expected folder names like '881_A1', '881_B3', '882_H12'.\n"
            f"Each folder must contain one BAM file (.bam) and its index (.bam.bai)."
        )

    if plate_ids:
        unknown = sorted(set(plate_ids) - set(plates))
        if unknown:
            click.echo(f"Warning: plate ID(s) not found and will be skipped: {unknown}", err=True)
        plates = {k: v for k, v in plates.items() if k in plate_ids}
        if not plates:
            raise click.ClickException(
                f"None of the requested plate IDs were found: {list(plate_ids)}"
            )

    total_wells = sum(len(v) for v in plates.values())
    click.echo(
        f"\nDiscovered {total_wells} wells across {len(plates)} plate(s): "
        + ", ".join(sorted(plates.keys()))
    )

    # ------------------------------------------------------------------ #
    # Load optional sample sheet
    # ------------------------------------------------------------------ #
    well_metadata: dict = {}
    if sample_sheet:
        well_metadata = parse_sample_sheet(sample_sheet, sequencer=sequencer)
        click.echo(f"Sample sheet loaded: {len(well_metadata)} wells.")

        # Cross-check BAMs vs sample sheet — warn on mismatches, never fail.
        bams_without_sheet = [
            f"{pid}_{wid}"
            for pid, wlist in plates.items()
            for wid, _ in wlist
            if lookup_well(well_metadata, pid, wid) is None
        ]
        bam_wells_only = {wid.upper() for wlist in plates.values() for wid, _ in wlist}
        bam_combined   = {
            f"{pid.upper()}/{wid.upper()}"
            for pid, wlist in plates.items()
            for wid, _ in wlist
        }
        sheets_without_bam = []
        for key in well_metadata:
            if "/" in key:
                found = key in bam_combined
            else:
                found = key in bam_wells_only
            if not found:
                sheets_without_bam.append(key)

        if bams_without_sheet:
            click.echo(
                f"Warning: {len(bams_without_sheet)} BAM well(s) have no matching sample "
                f"sheet entry — barcode metadata will be absent for these wells.",
                err=True,
            )
        if sheets_without_bam:
            click.echo(
                f"Warning: {len(sheets_without_bam)} sample sheet entry/entries have no "
                f"corresponding BAM and will be ignored.",
                err=True,
            )

    # ------------------------------------------------------------------ #
    # Resolve annotation, polyA, and TSS (shared across all plates/wells)
    # ------------------------------------------------------------------ #

    # Need at least one BAM to check chromosome naming before loading GTF
    first_bam = next(iter(next(iter(plates.values()))))[1]

    # Whitelist
    whitelist = _resolve_whitelist(barcode_whitelist, whitelist_db, offline)

    # Resolve platform/stage overrides
    plat  = None if platform == "auto"       else Platform(platform)
    stage = None if pipeline_stage == "auto" else PipelineStage(pipeline_stage)

    # Resolve GTF
    gtf, gtf_version, gtf_source = _resolve_gtf(gtf, gtf_version=gtf_version, offline=offline)

    # Quick chromosome naming check against the first BAM
    _first_meta = inspect_bam(first_bam, barcode_tag=barcode_tag, umi_tag=umi_tag,
                               platform=plat, pipeline_stage=stage)
    _validate_chromosome_naming(_first_meta, gtf)

    # Resolve polyA + TSS
    polya_paths, polya_version, polya_source = _resolve_polya_sites(
        polya_sites, polya_db=polya_db,
        hint_gencode_version=gtf_version if gtf_source == "auto-downloaded" else None,
        offline=offline,
    )
    version_warning = _check_version_consistency(gtf_version, polya_version)
    tss_paths = _resolve_tss_sites(tss_sites, tss_db=tss_db, offline=offline)

    # Build annotation index ONCE for all wells
    index = build_annotation_index(
        gtf,
        repeats_path=repeats,
        exclude_biotypes=list(exclude_biotypes),
        cache=not no_cache,
    )

    # Pre-load polyA and TSS site dicts ONCE — reused for every well.
    # Loading from the compressed BED takes ~35 s; doing it per-well would
    # multiply that cost by the number of wells (thousands for a full plate run).
    _shared_bam_cs    = _detect_chrom_style(_first_meta.reference_names)
    _shared_polya     = _load_polya_sites(polya_paths, chrom_style=_shared_bam_cs)
    _shared_tss       = _load_tss_sites(tss_paths, chrom_style=_shared_bam_cs) if tss_paths else None

    # ------------------------------------------------------------------ #
    # Process plates
    # ------------------------------------------------------------------ #
    for plate_id, well_list in sorted(plates.items()):
        plate_label = sample_name or plate_id
        click.echo(f"\n{'='*60}")
        click.echo(f"Plate: {plate_id}  ({len(well_list)} wells)")
        click.echo("=" * 60)

        plate_output_dir = output_dir_path / plate_id
        plate_output_dir.mkdir(parents=True, exist_ok=True)

        plate_result: "SampleResult | None" = None
        per_well_rows: list = []
        n_failed = 0

        # Wells without a BAI index are skipped before dispatching
        indexed_wells = []
        for well_id, bam_path in well_list:
            bai_path = Path(str(bam_path) + ".bai")
            bai_path2 = bam_path.with_suffix(".bai")
            if not bai_path.exists() and not bai_path2.exists():
                click.echo(f"  [{well_id}] skipping — no .bai index for {bam_path.name}", err=True)
                n_failed += 1
            else:
                indexed_wells.append((well_id, bam_path))

        if not indexed_wells:
            click.echo(
                f"  Plate {plate_id} has no indexed wells — skipping "
                f"(index all BAMs with `samtools index` before running run-plate).",
                err=True,
            )
            continue

        if parallel_wells > 1 and indexed_wells:
            # ---------------------------------------------------------------- #
            # Parallel path — one worker process per well (up to parallel_wells)
            # ---------------------------------------------------------------- #
            _threads_per_well = max(1, threads // parallel_wells)
            click.echo(
                f"  Running {len(indexed_wells)} wells with {parallel_wells} parallel workers "
                f"({_threads_per_well} thread(s) each) …"
            )
            _tasks = [
                {
                    "bam_path":         str(bam_path),
                    "well_id":          well_id,
                    "plate_id":         plate_id,
                    "barcode_tag":      barcode_tag,
                    "umi_tag":          umi_tag,
                    "platform_override": plat,
                    "stage_override":   stage,
                    "whitelist":        whitelist,
                    "chimeric_distance": chimeric_distance,
                    "threads":          _threads_per_well,
                    "store_umis":       not no_umi_dedup,
                    "reference":        reference,
                }
                for well_id, bam_path in indexed_wells
            ]
            _well_results: dict = {}
            with ProcessPoolExecutor(
                max_workers=parallel_wells,
                initializer=_plate_worker_init,
                initargs=(
                    gtf, repeats, list(exclude_biotypes),
                    polya_paths, tss_paths, _shared_bam_cs,
                ),
            ) as executor:
                futures = {executor.submit(_plate_well_task, t): t["well_id"] for t in _tasks}
                try:
                    for future in as_completed(futures):
                        res = future.result()
                        wid = res["well_id"]
                        if res["ok"]:
                            _well_results[wid] = res
                            click.echo(
                                f"  [{wid}] {res['well_sm'].n_reads_total:,} reads, "
                                f"noise={res['well_sm'].noise_read_frac:.1%}"
                            )
                        else:
                            click.echo(f"  [{wid}] FAILED: {res['error']}", err=True)
                            n_failed += 1
                except BrokenProcessPool as pool_exc:
                    raise click.ClickException(
                        f"Worker process died during parallel well processing "
                        f"(completed {len(_well_results)}/{len(indexed_wells)} wells before failure). "
                        f"This usually means a worker was killed by the OS (OOM) or crashed. "
                        f"Try reducing --parallel-wells to lower peak memory usage. "
                        f"Error: {pool_exc}"
                    ) from pool_exc
                except Exception as pool_exc:
                    raise click.ClickException(
                        f"Worker pool failed during parallel well processing "
                        f"(completed {len(_well_results)}/{len(indexed_wells)} wells before failure). "
                        f"Error: {pool_exc}"
                    ) from pool_exc

            # Merge in original well order; reassign polya/tss after merge
            for well_id, bam_path in indexed_wells:
                if well_id not in _well_results:
                    continue
                res = _well_results[well_id]
                well_sm     = res["well_sm"]
                well_result = res["well_result"]
                well_meta_from_sheet = lookup_well(well_metadata, plate_id, well_id)
                row: dict = {
                    "plate_id":              plate_id,
                    "well_id":               well_id,
                    "n_reads_total":         well_sm.n_reads_total,
                    "n_reads_classified":    well_sm.n_reads_classified,
                    "noise_read_frac":       f"{well_sm.noise_read_frac:.4f}",
                    "noise_base_frac":       f"{well_sm.noise_base_frac:.4f}",
                    "strand_concordance":    f"{well_sm.strand_concordance:.4f}",
                    "chimeric_read_frac":    f"{well_sm.chimeric_read_frac:.4f}",
                    "multimapper_read_frac": f"{well_sm.multimapper_read_frac:.4f}",
                    "full_length_read_frac": (
                        f"{well_sm.full_length_read_frac:.4f}"
                        if well_sm.full_length_read_frac is not None else ""
                    ),
                    "n_tso_invasion":        well_sm.n_tso_invasion,
                    "n_polya_priming":       well_sm.n_polya_priming,
                    "bam_path":              str(bam_path),
                }
                if well_meta_from_sheet:
                    row.update({
                        "i5_name": well_meta_from_sheet.get("i5_name", ""),
                        "i7_name": well_meta_from_sheet.get("i7_name", ""),
                        "i5_seq":  well_meta_from_sheet.get("i5_seq",  ""),
                        "i7_seq":  well_meta_from_sheet.get("i7_seq",  ""),
                    })
                    for k, v in well_meta_from_sheet.items():
                        if k not in {"well_id", "plate_id", "i5_name", "i7_name",
                                     "i5_seq", "i7_seq"} and k not in row:
                            row[k] = v
                per_well_rows.append(row)
                if plate_result is None:
                    plate_result = well_result
                else:
                    merge_sample_results(plate_result, well_result)

            if plate_result is not None:
                plate_result._polya_site_dict = _shared_polya
                plate_result._tss_site_dict   = _shared_tss

        else:
            # ---------------------------------------------------------------- #
            # Serial path
            # ---------------------------------------------------------------- #
            for well_id, bam_path in indexed_wells:
                try:
                    meta = inspect_bam(
                        bam_path, barcode_tag=barcode_tag, umi_tag=umi_tag,
                        platform=plat, pipeline_stage=stage,
                    )
                    _validate_sort_order(meta)

                    well_result = run_pipeline(
                        bam_path, index, meta,
                        whitelist=whitelist,
                        chimeric_distance=chimeric_distance,
                        paired_end_chimeric=_is_illumina_platform(meta.platform),
                        threads=threads,
                        store_umis=not no_umi_dedup,
                        reference_path=reference,
                    )

                    well_result._polya_site_dict = _shared_polya
                    well_result._tss_site_dict   = _shared_tss
                    well_result._numt_intervals  = None

                    well_sm, _ = compute_metrics(
                        well_result, f"{plate_id}_{well_id}",
                        platform=meta.platform.value,
                        unstranded=(meta.platform == Platform.SMARTSEQ),
                    )

                    well_meta_from_sheet = lookup_well(well_metadata, plate_id, well_id)
                    row = {
                        "plate_id":              plate_id,
                        "well_id":               well_id,
                        "n_reads_total":         well_sm.n_reads_total,
                        "n_reads_classified":    well_sm.n_reads_classified,
                        "noise_read_frac":       f"{well_sm.noise_read_frac:.4f}",
                        "noise_base_frac":       f"{well_sm.noise_base_frac:.4f}",
                        "strand_concordance":    f"{well_sm.strand_concordance:.4f}",
                        "chimeric_read_frac":    f"{well_sm.chimeric_read_frac:.4f}",
                        "multimapper_read_frac": f"{well_sm.multimapper_read_frac:.4f}",
                        "full_length_read_frac": (
                            f"{well_sm.full_length_read_frac:.4f}"
                            if well_sm.full_length_read_frac is not None else ""
                        ),
                        "n_tso_invasion":        well_sm.n_tso_invasion,
                        "n_polya_priming":       well_sm.n_polya_priming,
                        "bam_path":              str(bam_path),
                    }
                    if well_meta_from_sheet:
                        row.update({
                            "i5_name": well_meta_from_sheet.get("i5_name", ""),
                            "i7_name": well_meta_from_sheet.get("i7_name", ""),
                            "i5_seq":  well_meta_from_sheet.get("i5_seq",  ""),
                            "i7_seq":  well_meta_from_sheet.get("i7_seq",  ""),
                        })
                        for k, v in well_meta_from_sheet.items():
                            if k not in {"well_id", "plate_id", "i5_name", "i7_name",
                                         "i5_seq", "i7_seq"} and k not in row:
                                row[k] = v

                    per_well_rows.append(row)
                    click.echo(
                        f"  [{well_id}] {well_sm.n_reads_total:,} reads, "
                        f"noise={well_sm.noise_read_frac:.1%}"
                    )

                    if plate_result is None:
                        plate_result = well_result
                    else:
                        merge_sample_results(plate_result, well_result)

                except Exception as exc:
                    logger.error("Well %s failed: %s", well_id, exc, exc_info=verbose)
                    click.echo(f"  [{well_id}] FAILED: {exc}", err=True)
                    n_failed += 1

        if plate_result is None:
            click.echo(f"  All wells failed for plate {plate_id} — skipping.", err=True)
            continue

        # ------------------------------------------------------------------ #
        # Plate-level metrics + outputs
        # ------------------------------------------------------------------ #
        plate_meta = plate_result.meta  # use last well's meta (platform is the same)
        plate_meta.warnings = []        # clear per-well warnings; plate report stays clean

        # Intergenic profiler at plate level (runs before compute_metrics)
        intergenic_records = extract_intergenic_records(plate_result)
        intergenic_loci = []
        if intergenic_records:
            total_ig_bases = compute_intergenic_bases(index)
            ref_handle = None
            if reference:
                try:
                    import pysam as _pysam
                    ref_handle = _pysam.FastaFile(reference)
                except Exception as _exc:
                    logger.warning("Could not open reference FASTA for polyA context: %s", _exc)
            try:
                ig_loci, record_categories = profile_intergenic_loci(
                    intergenic_records,
                    total_intergenic_bases=total_ig_bases,
                    total_barcodes=len(plate_result.read_counts),
                    reference=ref_handle,
                    polya_sites=getattr(plate_result, "_polya_site_dict", None) or {},
                )
            finally:
                if ref_handle is not None:
                    ref_handle.close()
            intergenic_loci = ig_loci
            _apply_intergenic_reclassification(plate_result, intergenic_records, record_categories)

        plate_sm, plate_ct = compute_metrics(
            plate_result, plate_label,
            platform=plate_meta.platform.value,
            unstranded=(plate_meta.platform == Platform.SMARTSEQ),
        )
        plate_sm._gtf_info   = {"version": gtf_version, "source": gtf_source, "path": gtf}
        plate_sm._polya_info = {
            "version": polya_version, "source": polya_source,
            "path": polya_paths[0] if polya_paths else None, "db": polya_db,
        }
        plate_sm._tss_info = {
            "db": tss_db, "source": "auto-downloaded",
            "path": tss_paths[0] if tss_paths else None,
        } if tss_paths else None
        plate_sm._cell_barcodes_info = None
        if version_warning:
            plate_sm.warnings.append(version_warning)

        _write_run_outputs(
            plate_sm, plate_ct, plate_result, plate_output_dir, plate_label,
            cluster_df=None, intergenic_loci=intergenic_loci,
            platform=plate_meta.platform,
        )

        # Per-well metrics TSV
        if per_well_rows:
            import pandas as _pd
            pw_path = plate_output_dir / f"{plate_label}.per_well_metrics.tsv"
            _pd.DataFrame(per_well_rows).to_csv(pw_path, sep="\t", index=False)
            logger.info("Wrote %s", pw_path)

        n_wells_ok = len(well_list) - n_failed
        click.echo(
            f"\n  Plate {plate_id} complete: {n_wells_ok}/{len(well_list)} wells processed."
        )
        click.echo(f"  Noise fraction (reads) : {plate_sm.noise_read_frac:.2%}")
        click.echo(f"  Strand concordance     : {plate_sm.strand_concordance:.2%}")
        click.echo(f"  Chimeric rate          : {plate_sm.chimeric_read_frac:.2%}")
        click.echo(f"  Results in             : {plate_output_dir}")


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _open_bed(path: str):
    """Open a BED file, handling both plain and gzip-compressed files."""
    import gzip
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def _strip_chr_if_needed(sites: dict, chrom_style: str) -> dict:
    """
    If the BAM uses Ensembl-style chromosome names (no chr prefix) but the BED
    uses UCSC-style (chr prefix), strip the chr prefix so lookups succeed.
    No-op when styles already match or when style is unknown.
    """
    if chrom_style != "ensembl":
        return sites
    # Check if the keys actually have the chr prefix before stripping
    has_chr = any(k.startswith("chr") for k in sites)
    if not has_chr:
        return sites
    stripped = {}
    for k, v in sites.items():
        stripped[k.removeprefix("chr")] = v
    logger.debug("Stripped 'chr' prefix from BED chromosome names to match Ensembl BAM naming.")
    return stripped


def _site_cache_path(paths, chrom_style: str, prefix: str) -> Optional[Path]:
    """
    Return a pickle cache path for a loaded BED site dict, keyed on source
    file mtimes + sizes + a content hash of the first 64 KB + chrom_style.
    The head-bytes hash guards against in-place edits that preserve mtime
    (rare: ``touch -d``, some editor save modes, NFS quirks).
    Returns None when caching is not applicable (no paths, or stat fails).
    """
    if isinstance(paths, str):
        paths = [paths]
    if not paths:
        return None
    path_objs = [Path(p) for p in paths]
    try:
        parts = []
        for p in path_objs:
            st = p.stat()
            with open(p, "rb") as fh:
                head = fh.read(65536)
            head_hash = hashlib.md5(head).hexdigest()[:10]
            parts.append(f"{p}:{st.st_mtime_ns}:{st.st_size}:{head_hash}")
        parts.append(chrom_style)
        h = hashlib.md5("|".join(parts).encode()).hexdigest()[:10]
    except OSError:
        return None
    cache_dir = Path.home() / ".cache" / "scnoisemeter"
    return cache_dir / f".scnoisemeter_{prefix}_{h}.pkl.gz"


def _load_polya_sites(paths, chrom_style: str = "ucsc") -> dict:
    """
    Load one or more polyA site BED files into a merged dict of
    contig → sorted list of midpoint positions.

    Accepts a single path string or a list/tuple of paths. Multiple
    databases (e.g. PolyASite 2.0 + APADB) are merged transparently —
    duplicate positions are collapsed after sorting.

    *chrom_style* is the style detected from the BAM header ('ucsc' or
    'ensembl').  When 'ensembl', the chr prefix is stripped from BED keys
    so lookups against Ensembl-aligned BAM reads succeed.

    Supported databases:
      - PolyASite 2.0  (atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz)
      - APADB 3.0      (APADB_3.0.hg38.bed.gz)
      - Any BED3+ file with chrom / start / end in the first three columns
    """
    import gzip as _gz, pickle as _pkl
    _cache = _site_cache_path(paths, chrom_style, "polya")
    if _cache and _cache.exists():
        try:
            with _gz.open(_cache, "rb") as fh:
                sites = _pkl.load(fh)
            logger.info("Loaded polyA site dict from cache: %s", _cache.name)
            return sites
        except Exception:
            _cache.unlink(missing_ok=True)

    if isinstance(paths, str):
        paths = [paths]
    sites: dict = {}
    total = 0
    for path in paths:
        with _open_bed(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                try:
                    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                except ValueError:
                    continue
                mid = (start + end) // 2
                sites.setdefault(chrom, []).append(mid)
                total += 1
        logger.info("Loaded polyA sites from %s", path)
    for chrom in sites:
        sites[chrom] = sorted(set(sites[chrom]))
    sites = _strip_chr_if_needed(sites, chrom_style)
    logger.info(
        "Total polyA site positions: %d across %d contigs",
        sum(len(v) for v in sites.values()), len(sites),
    )

    if _cache:
        try:
            _cache.parent.mkdir(parents=True, exist_ok=True)
            with _gz.open(_cache, "wb") as fh:
                _pkl.dump(sites, fh, protocol=5)
        except Exception:
            pass

    return sites


def _load_tss_sites(paths, chrom_style: str = "ucsc") -> dict:
    """
    Load one or more CAGE peak / TSS BED files for the 5'-anchored
    full-length metric.

    Accepts a single path or list of paths. Format: BED3+ (chrom, start, end).
    Midpoint of each peak is used as the reference TSS position.

    *chrom_style* is the style detected from the BAM header ('ucsc' or
    'ensembl').  When 'ensembl', the chr prefix is stripped from BED keys
    so lookups against Ensembl-aligned BAM reads succeed.

    Compatible databases:
      - FANTOM5 CAGE peaks (hg38.cage_peak_phase1and2combined_ann.bed.gz)
      - ENCODE RAMPAGE peaks
      - Any BED3+ file
    """
    import gzip as _gz, pickle as _pkl
    _cache = _site_cache_path(paths, chrom_style, "tss")
    if _cache and _cache.exists():
        try:
            with _gz.open(_cache, "rb") as fh:
                sites = _pkl.load(fh)
            logger.info("Loaded TSS site dict from cache: %s", _cache.name)
            return sites
        except Exception:
            _cache.unlink(missing_ok=True)

    if isinstance(paths, str):
        paths = [paths]
    sites: dict = {}
    for path in paths:
        with _open_bed(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                try:
                    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                except ValueError:
                    continue
                mid = (start + end) // 2
                sites.setdefault(chrom, []).append(mid)
        logger.info("Loaded TSS sites from %s", path)
    for chrom in sites:
        sites[chrom] = sorted(set(sites[chrom]))
    sites = _strip_chr_if_needed(sites, chrom_style)
    logger.info(
        "Total TSS positions: %d across %d contigs",
        sum(len(v) for v in sites.values()), len(sites),
    )

    if _cache:
        try:
            _cache.parent.mkdir(parents=True, exist_ok=True)
            with _gz.open(_cache, "wb") as fh:
                _pkl.dump(sites, fh, protocol=5)
        except Exception:
            pass

    return sites


def _load_repeats_bed(path: str) -> dict:
    """
    Load a RepeatMasker BED file into a dict of
    contig → list of (start, end) tuples.

    Compatible formats:
      - UCSC RepeatMasker track (rmsk.txt.gz or rmsk.bed.gz)
        Expected columns: chrom, start, end, [name, score, strand, ...]
      - Any BED3+ file
    """
    intervals: dict = {}
    with _open_bed(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            try:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            except ValueError:
                continue
            intervals.setdefault(chrom, []).append((start, end))
    n = sum(len(v) for v in intervals.values())
    logger.info("Loaded %d RepeatMasker intervals from %s", n, path)
    return intervals


def _load_numt_bed(path: str) -> dict:
    """
    Load a NUMT (Nuclear Mitochondrial DNA segment) BED file into a dict of
    contig → list of (start, end) tuples.

    NUMTs are nuclear genome segments that are copies of mitochondrial DNA.
    Reads from genuine mitochondrial transcripts may map to NUMT loci instead
    of (or in addition to) chrM, and NUMT-derived reads may be erroneously
    classified as mitochondrial. When a NUMT BED is provided, reads classified
    as MITOCHONDRIAL are cross-checked: if the read overlaps a NUMT interval
    in the nuclear genome, it is flagged as potentially NUMT-derived.

    Compatible sources:
      - UCSC Table Browser: hg38 > Variation & Repeats > NUMTs
        (numts.txt.gz converted to BED3)
      - Any BED3+ file with nuclear genome NUMT coordinates
    """
    intervals: dict = {}
    with _open_bed(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            try:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            except ValueError:
                continue
            intervals.setdefault(chrom, []).append((start, end))
    n = sum(len(v) for v in intervals.values())
    logger.info("Loaded %d NUMT intervals from %s", n, path)
    return intervals


def _write_run_outputs(sm, ct, result, output_dir: Path, sample_name: str,
                       cluster_df=None, intergenic_loci=None, platform=None) -> None:
    import json

    # Sample-wide metrics TSV
    metrics_path = output_dir / f"{sample_name}.read_metrics.tsv"
    rows = [("metric", "value")]
    for attr in ["n_reads_total", "n_reads_classified", "n_reads_unassigned",
                 "n_cells", "noise_read_frac", "noise_base_frac",
                 "strand_concordance", "chimeric_read_frac",
                 "multimapper_read_frac", "per_cell_noise_median",
                 "per_cell_noise_iqr", "n_tso_invasion", "n_polya_priming",
                 "n_noncanon_junction"]:
        rows.append((attr, getattr(sm, attr, "")))
    for cat_name, frac in sm.read_fracs.items():
        rows.append((f"read_frac_{cat_name}", f"{frac:.6f}"))
    for cat_name, frac in sm.base_fracs.items():
        rows.append((f"base_frac_{cat_name}", f"{frac:.6f}"))
    if sm.full_length_read_frac is not None:
        rows.append(("full_length_read_frac", f"{sm.full_length_read_frac:.6f}"))

    with open(metrics_path, "w") as fh:
        for k, v in rows:
            fh.write(f"{k}\t{v}\n")
    logger.info("Wrote %s", metrics_path)

    # Per-cell metrics TSV
    if not ct.df.empty:
        cell_path = output_dir / f"{sample_name}.cell_metrics.tsv"
        ct.df.to_csv(cell_path, sep="\t")
        logger.info("Wrote %s", cell_path)

    # MultiQC JSON
    mq_path = output_dir / f"{sample_name}.multiqc.json"
    with open(mq_path, "w") as fh:
        json.dump(to_multiqc_json(sm), fh, indent=2)
    logger.info("Wrote %s", mq_path)

    # Length distribution TSVs
    _write_length_distributions(result.length_samples, output_dir, sample_name)

    # Per-cluster metrics TSV
    if cluster_df is not None and not cluster_df.empty:
        cluster_path = output_dir / f"{sample_name}.cluster_metrics.tsv"
        cluster_df.to_csv(cluster_path, sep="\t")
        logger.info("Wrote %s", cluster_path)

    # Intergenic loci TSV
    if intergenic_loci:
        import pandas as pd
        loci_path = output_dir / f"{sample_name}.intergenic_loci.tsv"
        loci_rows = [
            {
                "contig":               l.contig,
                "start":                l.start,
                "end":                  l.end,
                "strand":               l.strand,
                "n_reads":              l.n_reads,
                "n_barcodes":           l.n_barcodes,
                "has_splice_evidence":  l.has_splice_evidence,
                "is_monoexonic":        l.is_monoexonic,
                "polya_run_downstream": l.polya_run_downstream,
                "near_polya_site":      l.near_polya_site,
                "poisson_pvalue_adj":   l.poisson_pvalue_adj,
                "category":             l.category.value,
            }
            for l in intergenic_loci
        ]
        pd.DataFrame(loci_rows).to_csv(loci_path, sep="\t", index=False)
        logger.info("Wrote %s", loci_path)

    # Length-stratified TSV
    strat_df = compute_length_stratification(
        result.length_bin_counts,
        result.length_samples,
    )
    _illumina_platform_vals = {"illumina", "illumina_10x", "illumina_bd", "smartseq"}
    is_illumina = (
        platform is not None
        and platform.value in _illumina_platform_vals
    )
    if not strat_df.empty:
        strat_path = output_dir / f"{sample_name}_length_stratified.tsv"
        if is_illumina:
            with open(strat_path, "w") as fh:
                fh.write(
                    "# NOTE: Single-bin output expected for Illumina short-read data "
                    "(all reads <150 bp). See insert_size field in cell_metrics for "
                    "per-molecule size information.\n"
                )
                strat_df.to_csv(fh, sep="\t", index=False)
        else:
            strat_df.to_csv(strat_path, sep="\t", index=False)
        logger.info("Wrote %s", strat_path)

    # Insert size distribution (Illumina only)
    insert_sizes = None
    if is_illumina:
        signal = getattr(result, "insert_size_signal", [])
        noise  = getattr(result, "insert_size_noise",  [])
        if signal or noise:
            insert_sizes = {"signal": signal, "noise": noise}

    # HTML report
    report_path = output_dir / f"{sample_name}.report.html"
    write_run_report(sm, ct, result.length_samples, report_path,
                     cluster_df=cluster_df,
                     intergenic_loci=intergenic_loci or [],
                     length_stratified=strat_df if not strat_df.empty else None,
                     platform=platform.value if platform is not None else None,
                     insert_sizes=insert_sizes)


def _write_length_distributions(length_samples: dict, output_dir: Path, sample_name: str) -> None:
    ld_dir = output_dir / f"{sample_name}.length_distributions"
    ld_dir.mkdir(exist_ok=True)
    for cat, lengths in length_samples.items():
        if not lengths:
            continue
        cat_name = cat.value if hasattr(cat, "value") else str(cat)
        path = ld_dir / f"{cat_name}.lengths.tsv"
        with open(path, "w") as fh:
            fh.write("read_length\n")
            for L in lengths:
                fh.write(f"{L}\n")


def _write_compare_outputs(results: dict, output_dir: Path, label_a: str, label_b: str) -> None:
    """Write side-by-side metrics TSV and statistical comparison."""
    import pandas as pd
    from scipy.stats import chi2_contingency

    sm_a, ct_a, _ = results[label_a]
    sm_b, ct_b, _ = results[label_b]

    # Side-by-side scalar metrics
    metrics_a = {k: v for k, v in vars(sm_a).items() if isinstance(v, (int, float))}
    metrics_b = {k: v for k, v in vars(sm_b).items() if isinstance(v, (int, float))}
    all_keys = sorted(set(metrics_a) | set(metrics_b))

    cmp_path = output_dir / "comparison.metrics.tsv"
    with open(cmp_path, "w") as fh:
        fh.write(f"metric\t{label_a}\t{label_b}\tdelta\n")
        for k in all_keys:
            va = metrics_a.get(k, float("nan"))
            vb = metrics_b.get(k, float("nan"))
            try:
                delta = vb - va
            except TypeError:
                delta = ""
            fh.write(f"{k}\t{va}\t{vb}\t{delta}\n")

    # Paired McNemar-style proportion tests per category
    # Note: BAM-B is a subset of BAM-A (post-filter ⊆ pre-filter).
    # We use a chi-squared test on the contingency table of read counts.
    stats_rows = []
    for cat in ReadCategory:
        n_a = sm_a.n_reads_classified * sm_a.read_fracs.get(cat.value, 0)
        n_b = sm_b.n_reads_classified * sm_b.read_fracs.get(cat.value, 0)
        total_a = sm_a.n_reads_classified
        total_b = sm_b.n_reads_classified
        if total_a == 0 or total_b == 0:
            continue
        contingency = [[int(n_a), int(total_a - n_a)],
                       [int(n_b), int(total_b - n_b)]]
        try:
            chi2, p, _, _ = chi2_contingency(contingency)
        except Exception:
            chi2, p = float("nan"), float("nan")
        stats_rows.append({
            "category": cat.value,
            f"frac_{label_a}": sm_a.read_fracs.get(cat.value, 0),
            f"frac_{label_b}": sm_b.read_fracs.get(cat.value, 0),
            "delta": sm_b.read_fracs.get(cat.value, 0) - sm_a.read_fracs.get(cat.value, 0),
            "chi2": chi2,
            "p_value": p,
        })

    if stats_rows:
        import pandas as pd
        stats_df = pd.DataFrame(stats_rows)
        # Bonferroni correction
        n_tests = stats_df["p_value"].notna().sum()
        stats_df["p_adjusted"] = (stats_df["p_value"] * n_tests).clip(upper=1.0)
        stats_path = output_dir / "comparison.stats.tsv"
        stats_df.to_csv(stats_path, sep="\t", index=False)
        logger.info("Wrote %s", stats_path)

        # HTML comparison report
        sm_a, ct_a, result_a = results[label_a]
        sm_b, ct_b, result_b = results[label_b]
        report_path = output_dir / "comparison.report.html"
        write_compare_report(
            sm_a, sm_b, ct_a, ct_b,
            result_a.length_samples, result_b.length_samples,
            stats_df,
            report_path,
        )
        logger.info("Wrote %s", report_path)


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _load_whitelist(path: Optional[str]) -> Optional[set]:
    if path is None:
        return None
    import gzip
    logger.info("Loading barcode whitelist: %s", path)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        whitelist = {line.strip() for line in fh if line.strip()}
    logger.info("  %d barcodes loaded.", len(whitelist))
    return whitelist


def _load_cell_barcodes(path: str) -> set:
    """
    Load a called-cell barcode file (plain text or gzip, one barcode per line).

    Trailing -1 suffixes are stripped from every entry so that Cell Ranger
    barcodes.tsv.gz files (which carry the -1) are normalised to bare
    sequences.  The CB tag values in the BAM are stripped the same way
    inside the pipeline worker before comparison.
    """
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    barcodes: set = set()
    with opener(path, "rt") as fh:
        for line in fh:
            bc = line.strip().removesuffix("-1")
            if bc:
                barcodes.add(bc)
    logger.info("Loaded %d called-cell barcodes from %s", len(barcodes), path)
    return barcodes


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stderr,
    )


if __name__ == "__main__":
    cli()
