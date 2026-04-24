"""
metrics.py
==========
Converts the raw per-cell count dictionaries from pipeline.py into clean,
interpretable metrics.

Two output levels are produced:

  SampleMetrics  –  sample-wide scalar metrics (QC mode / MultiQC output)
  CellMetrics    –  per-cell DataFrame (per-cell QC and benchmarking)

All fractions are computed twice:
  - read_frac_*  : fraction of classified reads in each category
  - base_frac_*  : fraction of classified bases in each category

"Classified" excludes UNMAPPED, SECONDARY, SUPPLEMENTARY from the denominator.
UNASSIGNED reads are included in the denominator (they are real reads, just
unattributable to a cell).

UMI complexity ratio
--------------------
  umi_complexity_<category> = unique_UMIs / total_reads  (per cell per category)
A ratio close to 1.0 means every read represents a distinct molecule
(consistent with genuine transcription or rare events).
A ratio << 1.0 means many reads share UMIs (expected for high-coverage
genuine transcripts; suspicious if seen in noise categories).

Full-length read fraction (long-read specific)
----------------------------------------------
Estimated from the fraction of exonic-sense reads whose length exceeds a
platform-specific threshold (default: 500 bp for ONT, 1000 bp for PacBio).
A proxy for RT completeness.

Strand concordance ratio
------------------------
  strand_concordance = exonic_sense / (exonic_sense + exonic_antisense)
Values < 0.95 suggest strand-switching artifacts or a non-stranded library.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from scnoisemeter.constants import (
    AMBIGUOUS_CATEGORIES,
    CATEGORY_ORDER,
    LENGTH_BIN_BREAKS,
    LENGTH_BIN_LABELS_LONG,
    LENGTH_BIN_LABELS_SHORT,
    LENGTH_SHORT_READ_THRESHOLD,
    NOISE_CATEGORIES,
    NOISE_CATEGORIES_STRICT,
    NOISE_CATEGORIES_STRICT_UNSTRANDED,
    NOISE_CATEGORIES_UNSTRANDED,
    ReadCategory,
)

# Ensure new sub-categories are excluded from noise (they go in ambiguous bucket)
_AMBIGUOUS_SUB = {ReadCategory.AMBIGUOUS_COD_COD, ReadCategory.AMBIGUOUS_COD_NCOD}
from scnoisemeter.modules.pipeline import SampleResult

logger = logging.getLogger(__name__)

# Minimum reads per cell to include in per-cell table (avoids noise from
# extremely low-coverage barcodes)
MIN_READS_PER_CELL = 10

# Fallback full-length read thresholds by platform (bp)
# Used only when no polyA site database is provided via --polya-sites.
# When a polyA site database IS provided, a read is considered full-length
# if its 3' end falls within POLYA_SITE_PROXIMITY bp of a known polyA site —
# a much more biologically meaningful criterion.
FULL_LENGTH_THRESHOLD = {
    "ont":      500,
    "pacbio":  1000,
    "default":  500,
}


# ---------------------------------------------------------------------------
# Output dataclasses
# ---------------------------------------------------------------------------

@dataclass
class SampleMetrics:
    """
    Sample-wide scalar QC metrics.

    All fraction fields are floats in [0, 1].
    All count fields are integers.
    """
    sample_name:           str
    bam_path:              str
    platform:              str
    pipeline_stage:        str
    aligner:               str
    n_reads_total:         int = 0
    n_reads_classified:    int = 0
    n_reads_unassigned:    int = 0
    n_cells:               int = 0

    # Per-category read fractions (sample-wide)
    read_fracs:            dict = field(default_factory=dict)

    # Per-category base fractions (sample-wide)
    base_fracs:            dict = field(default_factory=dict)

    # Aggregate noise fractions
    noise_read_frac:       float = 0.0
    noise_base_frac:       float = 0.0

    # Strict noise: only unambiguous RT/PCR artifacts (excludes INTRONIC_PURE
    # and INTRONIC_BOUNDARY which may represent genuine pre-mRNA capture).
    # Conservative noise (noise_read_frac above) includes these and is an
    # upper bound. Strict noise is a lower bound. True noise lies between them.
    noise_read_frac_strict: float = 0.0
    noise_base_frac_strict: float = 0.0

    # Strand concordance
    strand_concordance:    float = 0.0

    # Chimeric rate
    chimeric_read_frac:    float = 0.0

    # Multimapper rate
    multimapper_read_frac: float = 0.0

    # True when the protocol is non-stranded (e.g. Smart-seq2 / FLASH-seq).
    # When True, noise_read_frac excludes EXONIC_ANTISENSE — those reads are
    # genuine cDNA signal from the reverse strand, not artifacts.
    is_unstranded: bool = False

    # Soft-clip fraction (reported separately — see bam_inspector)
    # UMI complexity per category (mean across cells)
    umi_complexity:        dict = field(default_factory=dict)

    # Full-length read fraction (long-read specific)
    full_length_read_frac: Optional[float] = None

    # Artifact flags (sample-wide counts)
    n_tso_invasion:        int = 0
    n_polya_priming:       int = 0
    n_noncanon_junction:   int = 0

    # Optional metrics — None when the required reference file was not provided
    tss_anchored_frac:     Optional[float] = None   # 5'-anchored at TSS (CAGE)
    numt_read_frac:        Optional[float] = None   # mito reads flagged as NUMT

    # Per-cell summary statistics (median ± IQR of noise fraction)
    per_cell_noise_median: float = 0.0
    per_cell_noise_iqr:    float = 0.0

    # Warnings from the BAM inspector
    warnings:              list = field(default_factory=list)


@dataclass
class CellTable:
    """
    Per-cell metrics as a pandas DataFrame.

    Index: cell_barcode (str)
    Columns include: n_reads, n_bases, read_frac_*, base_frac_*,
    umi_complexity_*, noise_read_frac, noise_base_frac,
    n_tso, n_polya, n_noncanon.
    """
    df: pd.DataFrame
    sample_name: str


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_metrics(
    result: SampleResult,
    sample_name: str,
    *,
    platform: str = "default",
    min_reads_per_cell: int = MIN_READS_PER_CELL,
    unstranded: bool = False,
) -> tuple[SampleMetrics, CellTable]:
    """
    Compute sample-wide and per-cell metrics from a :class:`SampleResult`.

    Returns
    -------
    (SampleMetrics, CellTable)
    """
    meta = result.meta

    sm = SampleMetrics(
        sample_name=sample_name,
        bam_path=str(result.bam_path),
        platform=meta.platform.value,
        pipeline_stage=meta.pipeline_stage.value,
        aligner=meta.aligner,
        warnings=list(meta.warnings),
        is_unstranded=unstranded,
    )

    # ------------------------------------------------------------------
    # 1. Build per-cell DataFrame
    # ------------------------------------------------------------------
    rows = []
    for cb, cat_read_counts in result.read_counts.items():
        total_reads = sum(cat_read_counts.values())
        if total_reads < min_reads_per_cell:
            continue

        total_bases = sum(result.base_counts[cb].values())
        row = {"cell_barcode": cb, "n_reads": total_reads, "n_bases": total_bases}

        # Read fractions
        for cat in CATEGORY_ORDER:
            row[f"read_frac_{cat.value}"] = cat_read_counts.get(cat, 0) / total_reads if total_reads else 0.0

        # Base fractions
        for cat in CATEGORY_ORDER:
            row[f"base_frac_{cat.value}"] = (
                result.base_counts[cb].get(cat, 0) / total_bases if total_bases else 0.0
            )

        # UMI complexity per category
        for cat in CATEGORY_ORDER:
            n_reads_cat = cat_read_counts.get(cat, 0)
            n_umis_cat  = len(result.umi_sets[cb].get(cat, set()))
            row[f"umi_complexity_{cat.value}"] = (
                n_umis_cat / n_reads_cat if n_reads_cat > 0 else float("nan")
            )

        # Aggregate noise — use unstranded set for non-stranded protocols
        _noise_cats = NOISE_CATEGORIES_UNSTRANDED if unstranded else NOISE_CATEGORIES
        noise_reads = sum(cat_read_counts.get(cat, 0) for cat in _noise_cats)
        row["noise_read_frac"] = noise_reads / total_reads if total_reads else 0.0

        noise_bases = sum(result.base_counts[cb].get(cat, 0) for cat in _noise_cats)
        row["noise_base_frac"] = noise_bases / total_bases if total_bases else 0.0

        # Artifact flags
        flags = result.artifact_flags.get(cb, {})
        row["n_tso"]      = flags.get("tso",      0)
        row["n_polya"]    = flags.get("polya",     0)
        row["n_noncanon"] = flags.get("noncanon",  0)

        rows.append(row)

    cell_df = pd.DataFrame(rows).set_index("cell_barcode") if rows else pd.DataFrame()

    # ------------------------------------------------------------------
    # 2. Sample-wide aggregation
    # ------------------------------------------------------------------
    sm.n_cells = len(cell_df)

    # Aggregate read / base counts across all barcodes
    total_reads_all: dict[ReadCategory, int] = {}
    total_bases_all: dict[ReadCategory, int] = {}

    for cb, cat_counts in result.read_counts.items():
        for cat, n in cat_counts.items():
            total_reads_all[cat] = total_reads_all.get(cat, 0) + n

    for cb, cat_counts in result.base_counts.items():
        for cat, n in cat_counts.items():
            total_bases_all[cat] = total_bases_all.get(cat, 0) + n

    sm.n_reads_total      = result.n_reads_processed
    sm.n_reads_classified = sum(total_reads_all.values())
    sm.n_reads_unassigned = total_reads_all.get(ReadCategory.UNASSIGNED, 0)

    denom_reads = sm.n_reads_classified or 1
    denom_bases = sum(total_bases_all.values()) or 1

    sm.read_fracs = {
        cat.value: total_reads_all.get(cat, 0) / denom_reads
        for cat in CATEGORY_ORDER
    }
    sm.base_fracs = {
        cat.value: total_bases_all.get(cat, 0) / denom_bases
        for cat in CATEGORY_ORDER
    }

    # Aggregate noise fractions — use unstranded set for non-stranded protocols
    _noise_cats = NOISE_CATEGORIES_UNSTRANDED if unstranded else NOISE_CATEGORIES
    sm.noise_read_frac = sum(
        sm.read_fracs.get(cat.value, 0.0) for cat in _noise_cats
    )
    sm.noise_base_frac = sum(
        sm.base_fracs.get(cat.value, 0.0) for cat in _noise_cats
    )
    _noise_cats_strict = NOISE_CATEGORIES_STRICT_UNSTRANDED if unstranded else NOISE_CATEGORIES_STRICT
    sm.noise_read_frac_strict = sum(
        sm.read_fracs.get(cat.value, 0.0) for cat in _noise_cats_strict
    )
    sm.noise_base_frac_strict = sum(
        sm.base_fracs.get(cat.value, 0.0) for cat in _noise_cats_strict
    )

    # Strand concordance
    es = total_reads_all.get(ReadCategory.EXONIC_SENSE, 0)
    ea = total_reads_all.get(ReadCategory.EXONIC_ANTISENSE, 0)
    sm.strand_concordance = es / (es + ea) if (es + ea) > 0 else float("nan")

    # Chimeric / multimapper rates
    sm.chimeric_read_frac    = sm.read_fracs.get(ReadCategory.CHIMERIC.value, 0.0)
    sm.multimapper_read_frac = sm.read_fracs.get(ReadCategory.MULTIMAPPER.value, 0.0)

    # Artifact flag totals
    for cb, flags in result.artifact_flags.items():
        sm.n_tso_invasion      += flags.get("tso",      0)
        sm.n_polya_priming     += flags.get("polya",    0)
        sm.n_noncanon_junction += flags.get("noncanon", 0)

    # UMI complexity (mean across cells, per category)
    if not cell_df.empty:
        for cat in CATEGORY_ORDER:
            col = f"umi_complexity_{cat.value}"
            if col in cell_df.columns:
                sm.umi_complexity[cat.value] = float(
                    cell_df[col].replace([float("inf"), float("-inf")], float("nan")).mean(skipna=True)
                )

    # Per-cell noise summary stats
    if not cell_df.empty and "noise_read_frac" in cell_df.columns:
        vals = cell_df["noise_read_frac"].dropna().values
        if len(vals) > 0:
            sm.per_cell_noise_median = float(np.median(vals))
            sm.per_cell_noise_iqr    = float(np.percentile(vals, 75) - np.percentile(vals, 25))

    # Full-length read fraction
    # Preferred method: fraction of exonic-sense reads whose 3′ end falls
    # within POLYA_SITE_PROXIMITY bp of an annotated polyA site.
    # This is a direct measure of whether the molecule was captured from
    # a genuine polyadenylation site, regardless of read length.
    # Fallback: fraction of exonic-sense reads above a platform-specific
    # minimum length threshold (used when no polyA site database is available).
    polya_site_dict = getattr(result, "_polya_site_dict", None)
    es_three_prime = getattr(result, "exonic_sense_three_prime", [])

    if polya_site_dict and es_three_prime:
        # polyA-site-anchored full-length fraction
        from scnoisemeter.modules.intergenic_profiler import _near_polya_site
        n_near = sum(
            1 for contig, pos in es_three_prime
            if _near_polya_site(polya_site_dict, contig, pos)
        )
        sm.full_length_read_frac = float(n_near / len(es_three_prime))
        sm._polya_sites_used = True
    else:
        # Length-based fallback
        fl_threshold = FULL_LENGTH_THRESHOLD.get(
            platform, FULL_LENGTH_THRESHOLD["default"]
        )
        es_lengths = result.length_samples.get(ReadCategory.EXONIC_SENSE, [])
        if es_lengths:
            sm.full_length_read_frac = float(
                sum(1 for L in es_lengths if L >= fl_threshold) / len(es_lengths)
            )
        sm._polya_sites_used = False

    # --- 5'-anchored at TSS fraction (requires --tss-sites) ---
    # Fraction of exonic-sense reads whose 5' mapping end falls within
    # TSS_SITE_PROXIMITY bp of an annotated transcription start site.
    # Combined with the 3'-anchored metric this provides a true full-length
    # estimate: both ends anchored at known transcript boundaries.
    tss_site_dict   = getattr(result, "_tss_site_dict", None)
    es_five_prime   = getattr(result, "exonic_sense_five_prime", [])
    if tss_site_dict and es_five_prime:
        from scnoisemeter.modules.intergenic_profiler import _near_polya_site
        from scnoisemeter.constants import TSS_SITE_PROXIMITY
        n_tss = sum(
            1 for contig, pos in es_five_prime
            if _near_polya_site(tss_site_dict, contig, pos,
                                proximity=TSS_SITE_PROXIMITY)
        )
        sm.tss_anchored_frac = float(n_tss / len(es_five_prime))
    else:
        sm.tss_anchored_frac = None

    # --- NUMT interval summary (requires --numt-bed) ---
    # Reports the number of loaded NUMT intervals as metadata.
    # Full per-read NUMT disambiguation requires dual-alignment and is
    # deferred to a future version; for now we flag the presence of NUMT
    # intervals in the report so users are aware.
    numt_intervals = getattr(result, "_numt_intervals", None)
    if numt_intervals:
        sm.numt_read_frac = float(
            sum(len(v) for v in numt_intervals.values())
        )   # repurposed as interval count until dual-alignment is available
    else:
        sm.numt_read_frac = None

    ct = CellTable(df=cell_df, sample_name=sample_name)
    return sm, ct


# ---------------------------------------------------------------------------
# Read-length stratification
# ---------------------------------------------------------------------------

def compute_length_stratification(
    length_bin_counts: dict,
    length_samples: dict,
) -> pd.DataFrame:
    """
    Aggregate per-bin, per-category read counts into a tidy DataFrame.

    Parameters
    ----------
    length_bin_counts : dict
        From SampleResult.length_bin_counts: category → {bin_idx → count}.
        Bin indices 0–5 correspond to LENGTH_BIN_BREAKS breakpoints.
    length_samples : dict
        From SampleResult.length_samples: category → reservoir-sampled lengths.
        Used only to estimate the overall median read length to decide whether
        to activate the short-read <150 bp bin.

    Returns
    -------
    pd.DataFrame with columns:
        length_bin        – human-readable bin label
        category          – ReadCategory value string
        count             – exact read count in this (bin, category) cell
        fraction_of_bin   – count / total reads in this bin
        fraction_of_total – count / total reads across all bins
    """
    import bisect

    # Estimate median from reservoir samples
    all_lengths = [L for lengths in length_samples.values() for L in lengths]
    median_length = float(np.median(all_lengths)) if all_lengths else 1000.0

    use_short_bin = median_length < LENGTH_SHORT_READ_THRESHOLD

    if use_short_bin:
        bin_labels = LENGTH_BIN_LABELS_LONG
        # bin index maps directly to label index
        label_for_idx = {i: bin_labels[i] for i in range(len(bin_labels))}
    else:
        bin_labels = LENGTH_BIN_LABELS_SHORT
        # merge idx 0 (<150) and idx 1 (150–500) both into "<500"
        label_for_idx = {
            0: "<500",
            1: "<500",
            2: "500–1000",
            3: "1000–2000",
            4: "2000–5000",
            5: ">5000",
        }

    # Accumulate (bin_label, category) → count
    tally: dict = {}
    for cat, bin_counts in length_bin_counts.items():
        for bin_idx, n in bin_counts.items():
            label = label_for_idx.get(bin_idx, ">5000")
            key = (label, cat)
            tally[key] = tally.get(key, 0) + n

    if not tally:
        return pd.DataFrame(
            columns=["length_bin", "category", "count",
                     "fraction_of_bin", "fraction_of_total"]
        )

    total_all = sum(tally.values())

    rows = []
    for bin_label in bin_labels:
        bin_total = sum(n for (lbl, _), n in tally.items() if lbl == bin_label)
        for cat in CATEGORY_ORDER:
            n = tally.get((bin_label, cat), 0)
            rows.append({
                "length_bin":        bin_label,
                "category":          cat.value,
                "count":             n,
                "fraction_of_bin":   n / bin_total   if bin_total   > 0 else 0.0,
                "fraction_of_total": n / total_all   if total_all   > 0 else 0.0,
            })

    df = pd.DataFrame(rows)
    # Drop rows with no reads to keep the TSV compact
    df = df[df["count"] > 0].reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# MultiQC-compatible JSON output
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Per-cluster noise metrics
# ---------------------------------------------------------------------------

def compute_cluster_metrics(
    ct: CellTable,
    obs_metadata: "pd.DataFrame",
) -> "pd.DataFrame":
    """
    Compute per-cluster noise summary statistics.

    Parameters
    ----------
    ct : CellTable
        Per-cell metrics table from compute_metrics().
    obs_metadata : pd.DataFrame
        Two-column DataFrame with columns 'cell_barcode' and 'cluster'.
        Typically loaded from a user-supplied --obs-metadata TSV.

    Returns
    -------
    cluster_df : pd.DataFrame
        One row per cluster, columns include:
        n_cells, median_noise_read_frac, iqr_noise_read_frac,
        median_exonic_sense_frac, median_intronic_pure_frac,
        median_chimeric_frac, median_tso_rate, median_polya_rate.
    """
    if ct.df.empty or obs_metadata.empty:
        return pd.DataFrame()

    # Merge cell metrics with cluster labels
    merged = ct.df.reset_index().merge(
        obs_metadata.rename(columns={
            "cell_barcode": "cell_barcode",
            "cluster": "_cluster",
        }),
        on="cell_barcode",
        how="inner",
    )

    if merged.empty or "_cluster" not in merged.columns:
        logger.warning(
            "No cell barcodes in obs_metadata matched cell metrics table. "
            "Check that barcode formats match (e.g. no trailing -1 suffix)."
        )
        return pd.DataFrame()

    rows = []
    for cluster, grp in merged.groupby("_cluster", observed=False):
        row = {"cluster": cluster, "n_cells": len(grp)}

        for col, out_name in [
            ("noise_read_frac",        "median_noise_read_frac"),
            ("noise_base_frac",        "median_noise_base_frac"),
            (f"read_frac_{ReadCategory.EXONIC_SENSE.value}",    "median_exonic_sense_frac"),
            (f"read_frac_{ReadCategory.INTRONIC_PURE.value}",   "median_intronic_pure_frac"),
            (f"read_frac_{ReadCategory.CHIMERIC.value}",        "median_chimeric_frac"),
            (f"read_frac_{ReadCategory.EXONIC_ANTISENSE.value}","median_exonic_antisense_frac"),
            (f"read_frac_{ReadCategory.AMBIGUOUS_COD_COD.value}","median_ambiguous_cod_cod_frac"),
        ]:
            if col in grp.columns:
                vals = grp[col].dropna().values
                row[out_name] = float(np.median(vals)) if len(vals) else float("nan")
                iqr_col = out_name.replace("median_", "iqr_")
                row[iqr_col] = float(
                    np.percentile(vals, 75) - np.percentile(vals, 25)
                ) if len(vals) > 1 else 0.0

        for flag_col, out_name in [
            ("n_tso",     "median_tso_rate"),
            ("n_polya",   "median_polya_rate"),
            ("n_noncanon","median_noncanon_rate"),
        ]:
            if flag_col in grp.columns and "n_reads" in grp.columns:
                rates = grp[flag_col] / grp["n_reads"].replace(0, float("nan"))
                row[out_name] = float(rates.median())

        rows.append(row)

    cluster_df = pd.DataFrame(rows).set_index("cluster")
    return cluster_df


def load_obs_metadata(path: str) -> pd.DataFrame:
    """
    Load cell barcode → cluster mapping from a TSV or CSV file.

    The file must have at minimum two columns named 'cell_barcode' and
    'cluster'.  Additional columns are ignored.  This format is compatible
    with the output of Seurat's `write.csv(meta.data)` and Scanpy's
    `adata.obs.to_csv()` after selecting the relevant columns.
    """
    import os
    sep = "\t" if path.endswith(".tsv") or path.endswith(".txt") else ","
    df = pd.read_csv(path, sep=sep, dtype=str)
    df.columns = [c.strip() for c in df.columns]

    # Accept common column name variants
    rename_map = {}
    for col in df.columns:
        if col.lower() in ("barcode", "cell_barcode", "cellbarcode", "cb"):
            rename_map[col] = "cell_barcode"
        elif col.lower() in ("cluster", "seurat_clusters", "leiden",
                              "louvain", "cell_type", "celltype"):
            rename_map[col] = "cluster"

    df = df.rename(columns=rename_map)

    if "cell_barcode" not in df.columns:
        raise ValueError(
            f"Could not find a barcode column in {path}. "
            "Expected one of: cell_barcode, barcode, CB."
        )
    if "cluster" not in df.columns:
        raise ValueError(
            f"Could not find a cluster column in {path}. "
            "Expected one of: cluster, seurat_clusters, leiden, cell_type."
        )

    logger.info(
        "Loaded obs metadata: %d cells, %d clusters from %s",
        len(df), df["cluster"].nunique(), os.path.basename(path),
    )
    return df[["cell_barcode", "cluster"]]


def to_multiqc_json(sm: SampleMetrics) -> dict:
    """
    Format SampleMetrics as a MultiQC custom content JSON dict.

    The returned dict should be written as:
      {"id": "scnoisemeter", "data": {sample_name: {...}}}
    """
    data = {
        "noise_read_frac":          round(sm.noise_read_frac,          4),
        "noise_base_frac":          round(sm.noise_base_frac,          4),
        "noise_read_frac_strict":   round(sm.noise_read_frac_strict,   4),
        "noise_base_frac_strict":   round(sm.noise_base_frac_strict,   4),
        "strand_concordance":       round(sm.strand_concordance,       4),
        "chimeric_read_frac":       round(sm.chimeric_read_frac,       4),
        "multimapper_read_frac":    round(sm.multimapper_read_frac,    4),
        "n_cells":                  sm.n_cells,
        "n_reads_classified":       sm.n_reads_classified,
        "per_cell_noise_median":    round(sm.per_cell_noise_median,    4),
        "per_cell_noise_iqr":       round(sm.per_cell_noise_iqr,       4),
    }
    # Add per-category read fractions
    for cat_name, frac in sm.read_fracs.items():
        data[f"read_frac_{cat_name}"] = round(frac, 4)

    if sm.full_length_read_frac is not None:
        data["full_length_read_frac"] = round(sm.full_length_read_frac, 4)

    return {"id": "scnoisemeter", "data": {sm.sample_name: data}}
