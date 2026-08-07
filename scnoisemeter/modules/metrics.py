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

UMI sequence-diversity ratio
----------------------------
  umi_sequence_diversity_<category> = unique UMI strings / reads

This is not a molecule-complexity estimate: UMIs are not grouped by gene or
genomic molecule and sequencing-error correction is not attempted.  The old
``umi_complexity_*`` columns remain as compatibility aliases.

End-anchoring fractions
-----------------------
Reported from the same exonic-sense read against strand-aware TSS and polyA
atlases.  No read-length-only fallback is presented as full-length evidence.

Strand concordance ratio
------------------------
  strand_concordance = exonic_sense / (exonic_sense + exonic_antisense)
Values < 0.95 suggest strand-switching artifacts or a non-stranded library.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from scnoisemeter.modules.pipeline import SampleResult

from scnoisemeter.constants import (
    CATEGORY_ORDER,
    LENGTH_BIN_LABELS_LONG,
    LENGTH_BIN_LABELS_SHORT,
    LENGTH_SHORT_READ_THRESHOLD,
    ARTIFACT_CANDIDATE_CATEGORIES,
    NOISE_CATEGORIES,
    NOISE_CATEGORIES_UNSTRANDED,
    ReadCategory,
)

# Ensure new sub-categories are excluded from noise (they go in ambiguous bucket)
_AMBIGUOUS_SUB = {ReadCategory.AMBIGUOUS_COD_COD, ReadCategory.AMBIGUOUS_COD_NCOD}

logger = logging.getLogger(__name__)

# Minimum reads per cell to include in per-cell table (avoids noise from
# extremely low-coverage barcodes)
MIN_READS_PER_CELL = 10


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
    n_reads_mapped:        int = 0
    n_reads_unmapped:      int = 0
    n_primary_mapped:      int = 0
    n_secondary:           int = 0
    n_supplementary:       int = 0
    n_qcfail:              int = 0
    n_duplicate:           int = 0
    n_cells:               int = 0

    # Precise aliases: BAM index counters are alignment-record counts, and
    # classification is per mapped primary alignment (mates count separately).
    n_records_total:       int = 0
    n_records_mapped:      int = 0
    n_records_unmapped:    int = 0
    n_alignments_classified: int = 0

    # Per-category read fractions (sample-wide)
    read_fracs:            dict = field(default_factory=dict)

    # Per-category base fractions (sample-wide)
    base_fracs:            dict = field(default_factory=dict)

    # Descriptive aggregate composition.  broad_noncanonical_* is everything
    # outside canonical exonic-sense signal; artifact_candidate_* is the subset
    # carrying positive artifact evidence.  Neither is a causal noise estimate.
    broad_noncanonical_read_frac: float = 0.0
    broad_noncanonical_base_frac: float = 0.0
    artifact_candidate_read_frac: float = 0.0
    artifact_candidate_base_frac: float = 0.0

    # Strand concordance
    strand_concordance:    float = 0.0

    # Chimeric rate
    chimeric_read_frac:    float = 0.0

    # Multimapper rate
    multimapper_read_frac: float = 0.0
    unmapped_read_frac:    float = 0.0
    low_mapq_read_frac:    float = 0.0
    mean_mapq:             Optional[float] = None
    mean_edit_distance:    Optional[float] = None
    softclip_base_frac:    Optional[float] = None

    # True when the protocol is non-stranded (e.g. Smart-seq2 / FLASH-seq).
    # When True, broad_noncanonical_* excludes EXONIC_ANTISENSE — those reads are
    # genuine cDNA signal from the reverse strand, not artifacts.
    is_unstranded: bool = False

    # Soft-clip fraction (reported separately — see bam_inspector)
    # UMI sequence diversity per category (mean across cells). The historical
    # ``umi_complexity`` name is retained as a compatibility alias.
    umi_sequence_diversity: dict = field(default_factory=dict)
    umi_complexity:        dict = field(default_factory=dict)

    # End anchoring. ``full_length_read_frac`` is a deprecated alias for
    # both_ends_anchored_frac and is only populated when both atlases exist.
    full_length_read_frac: Optional[float] = None
    three_prime_anchored_frac: Optional[float] = None
    five_prime_anchored_frac: Optional[float] = None
    both_ends_anchored_frac: Optional[float] = None

    # Artifact flags (sample-wide counts)
    n_tso_invasion:        int = 0
    n_polya_priming:       int = 0
    n_noncanon_junction:   int = 0
    n_tso_concatemer:      int = 0
    n_discordant_pair:     int = 0
    n_low_mapq:            int = 0

    # Optional metrics — None when the required reference file was not provided
    tss_anchored_frac:     Optional[float] = None   # deprecated endpoint alias
    n_numt_intervals_loaded: Optional[int] = None

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
    umi_complexity_*, broad_noncanonical_read_frac, broad_noncanonical_base_frac,
    n_tso, n_polya, n_noncanon, n_tso_concat.
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
        # A BAM without CB tags has sample-level metrics but no cells. Plate
        # mode relabels this sentinel to the well ID before aggregation.
        if cb in {"", "NO_BARCODE"}:
            continue
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

        # UMI sequence diversity per category.  This is not gene-aware
        # deduplication and must not be interpreted as molecule complexity.
        for cat in CATEGORY_ORDER:
            n_reads_cat = cat_read_counts.get(cat, 0)
            n_umis_cat  = len(result.umi_sets[cb].get(cat, set()))
            invalid_umi = cat in getattr(result, "_invalid_umi_categories", set())
            diversity = (
                float("nan")
                if invalid_umi or not getattr(result, "umi_tracking_enabled", True)
                else n_umis_cat / n_reads_cat if n_reads_cat > 0
                else float("nan")
            )
            row[f"umi_sequence_diversity_{cat.value}"] = diversity
            row[f"umi_complexity_{cat.value}"] = diversity

        # Aggregate composition — the unstranded set drops EXONIC_ANTISENSE,
        # which is expected signal rather than an artifact in those protocols.
        _broad_cats = NOISE_CATEGORIES_UNSTRANDED if unstranded else NOISE_CATEGORIES
        broad_reads = sum(cat_read_counts.get(cat, 0) for cat in _broad_cats)
        broad_bases = sum(result.base_counts[cb].get(cat, 0) for cat in _broad_cats)
        row["broad_noncanonical_read_frac"] = broad_reads / total_reads if total_reads else 0.0
        row["broad_noncanonical_base_frac"] = broad_bases / total_bases if total_bases else 0.0

        # ARTIFACT_CANDIDATE_CATEGORIES excludes EXONIC_ANTISENSE already, so it
        # needs no stranded/unstranded split.
        candidate_reads = sum(cat_read_counts.get(cat, 0) for cat in ARTIFACT_CANDIDATE_CATEGORIES)
        candidate_bases = sum(result.base_counts[cb].get(cat, 0) for cat in ARTIFACT_CANDIDATE_CATEGORIES)
        row["artifact_candidate_read_frac"] = candidate_reads / total_reads if total_reads else 0.0
        row["artifact_candidate_base_frac"] = candidate_bases / total_bases if total_bases else 0.0

        # Artifact flags
        flags = result.artifact_flags.get(cb, {})
        row["n_tso"]        = flags.get("tso",        0)
        row["n_polya"]      = flags.get("polya",      0)
        row["n_noncanon"]   = flags.get("noncanon",   0)
        row["n_tso_concat"] = flags.get("tso_concat", 0)
        row["n_discordant_pair"] = flags.get("discordant_pair", 0)
        row["n_low_mapq"] = flags.get("low_mapq", 0)
        row["low_mapq_read_frac"] = (
            row["n_low_mapq"] / total_reads if total_reads else 0.0
        )

        alignment = getattr(result, "alignment_sums", {}).get(cb, {})
        row["mean_mapq"] = (
            alignment.get("mapq_sum", 0) / total_reads
            if total_reads else float("nan")
        )
        n_nm = alignment.get("nm_observed", 0)
        row["mean_edit_distance"] = (
            alignment.get("nm_sum", 0) / n_nm if n_nm else float("nan")
        )
        query_bases = alignment.get("query_bases", 0)
        row["softclip_base_frac"] = (
            alignment.get("softclip_bases", 0) / query_bases
            if query_bases else float("nan")
        )

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

    sm.n_reads_total      = result.n_reads_total or (
        result.n_reads_mapped_index + result.n_reads_unmapped
    )
    sm.n_reads_mapped     = result.n_reads_mapped_index
    sm.n_reads_unmapped   = result.n_reads_unmapped
    sm.n_primary_mapped   = result.n_primary_mapped
    sm.n_secondary        = result.n_secondary
    sm.n_supplementary    = result.n_supplementary
    sm.n_qcfail           = result.n_qcfail
    sm.n_duplicate        = result.n_duplicate
    sm.n_reads_classified = sum(total_reads_all.values())
    sm.n_reads_unassigned = total_reads_all.get(ReadCategory.UNASSIGNED, 0)
    sm.n_records_total = sm.n_reads_total
    sm.n_records_mapped = sm.n_reads_mapped
    sm.n_records_unmapped = sm.n_reads_unmapped
    sm.n_alignments_classified = sm.n_reads_classified

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

    # Aggregate composition — the unstranded set drops EXONIC_ANTISENSE, which
    # is expected signal rather than an artifact in those protocols.
    _broad_cats = NOISE_CATEGORIES_UNSTRANDED if unstranded else NOISE_CATEGORIES
    sm.broad_noncanonical_read_frac = sum(
        sm.read_fracs.get(cat.value, 0.0) for cat in _broad_cats
    )
    sm.broad_noncanonical_base_frac = sum(
        sm.base_fracs.get(cat.value, 0.0) for cat in _broad_cats
    )
    # ARTIFACT_CANDIDATE_CATEGORIES excludes EXONIC_ANTISENSE already, so it
    # needs no stranded/unstranded split.
    sm.artifact_candidate_read_frac = sum(
        sm.read_fracs.get(cat.value, 0.0) for cat in ARTIFACT_CANDIDATE_CATEGORIES
    )
    sm.artifact_candidate_base_frac = sum(
        sm.base_fracs.get(cat.value, 0.0) for cat in ARTIFACT_CANDIDATE_CATEGORIES
    )

    # Strand concordance
    es = total_reads_all.get(ReadCategory.EXONIC_SENSE, 0)
    ea = total_reads_all.get(ReadCategory.EXONIC_ANTISENSE, 0)
    sm.strand_concordance = es / (es + ea) if (es + ea) > 0 else float("nan")

    # Chimeric / multimapper rates
    sm.chimeric_read_frac    = sm.read_fracs.get(ReadCategory.CHIMERIC.value, 0.0)
    sm.multimapper_read_frac = sm.read_fracs.get(ReadCategory.MULTIMAPPER.value, 0.0)
    sm.unmapped_read_frac = (
        sm.n_reads_unmapped / sm.n_reads_total if sm.n_reads_total else 0.0
    )

    # Artifact flag totals
    for cb, flags in result.artifact_flags.items():
        sm.n_tso_invasion      += flags.get("tso",        0)
        sm.n_polya_priming     += flags.get("polya",      0)
        sm.n_noncanon_junction += flags.get("noncanon",   0)
        sm.n_tso_concatemer    += flags.get("tso_concat", 0)
        sm.n_discordant_pair   += flags.get("discordant_pair", 0)
        sm.n_low_mapq          += flags.get("low_mapq", 0)

    alignment_totals = {
        key: sum(
            values.get(key, 0)
            for values in getattr(result, "alignment_sums", {}).values()
        )
        for key in ("mapq_sum", "nm_sum", "nm_observed", "softclip_bases", "query_bases")
    }
    if sm.n_reads_classified:
        sm.mean_mapq = alignment_totals["mapq_sum"] / sm.n_reads_classified
        sm.low_mapq_read_frac = sm.n_low_mapq / sm.n_reads_classified
    if alignment_totals["nm_observed"]:
        sm.mean_edit_distance = (
            alignment_totals["nm_sum"] / alignment_totals["nm_observed"]
        )
    if alignment_totals["query_bases"]:
        sm.softclip_base_frac = (
            alignment_totals["softclip_bases"] / alignment_totals["query_bases"]
        )

    # UMI sequence diversity (mean across cells, per category)
    if not cell_df.empty:
        for cat in CATEGORY_ORDER:
            col = f"umi_sequence_diversity_{cat.value}"
            if col in cell_df.columns:
                value = float(
                    cell_df[col].replace([float("inf"), float("-inf")], float("nan")).mean(skipna=True)
                )
                sm.umi_sequence_diversity[cat.value] = value
                sm.umi_complexity[cat.value] = value

    # Per-cell noise summary stats
    if not cell_df.empty and "broad_noncanonical_read_frac" in cell_df.columns:
        vals = cell_df["broad_noncanonical_read_frac"].dropna().values
        if len(vals) > 0:
            sm.per_cell_noise_median = float(np.median(vals))
            sm.per_cell_noise_iqr    = float(np.percentile(vals, 75) - np.percentile(vals, 25))

    # Strand-aware endpoint anchoring on the *same sampled read*.  A 3'-anchor
    # alone is not full-length and no length-only fallback is used.
    polya_site_dict = getattr(result, "_polya_site_dict", None)
    tss_site_dict   = getattr(result, "_tss_site_dict", None)
    endpoints = getattr(result, "exonic_sense_endpoints", [])
    if endpoints and not unstranded:
        from scnoisemeter.modules.intergenic_profiler import _near_polya_site
        from scnoisemeter.constants import TSS_SITE_PROXIMITY
        three_hits = []
        five_hits = []
        for contig, strand, five_pos, three_pos, _read_key in endpoints:
            three_hits.append(
                bool(polya_site_dict) and _near_polya_site(
                    polya_site_dict, contig, three_pos, strand=strand
                )
            )
            five_hits.append(
                bool(tss_site_dict) and _near_polya_site(
                    tss_site_dict, contig, five_pos,
                    proximity=TSS_SITE_PROXIMITY, strand=strand,
                )
            )
        if polya_site_dict:
            sm.three_prime_anchored_frac = float(np.mean(three_hits))
        if tss_site_dict:
            sm.five_prime_anchored_frac = float(np.mean(five_hits))
            sm.tss_anchored_frac = sm.five_prime_anchored_frac
        if polya_site_dict and tss_site_dict:
            sm.both_ends_anchored_frac = float(np.mean([
                a and b for a, b in zip(three_hits, five_hits)
            ]))
            sm.full_length_read_frac = sm.both_ends_anchored_frac
    sm._polya_sites_used = bool(polya_site_dict)

    # NUMT BED provenance only. Do not expose an interval count as a read
    # fraction until per-read dual-alignment evidence is implemented.
    numt_intervals = getattr(result, "_numt_intervals", None)
    if numt_intervals:
        sm.n_numt_intervals_loaded = sum(len(v) for v in numt_intervals.values())

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
        n_cells, median_broad_noncanonical_read_frac,
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
            ("broad_noncanonical_read_frac", "median_broad_noncanonical_read_frac"),
            ("artifact_candidate_read_frac", "median_artifact_candidate_read_frac"),
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
            ("n_tso",       "median_tso_rate"),
            ("n_polya",     "median_polya_rate"),
            ("n_noncanon",  "median_noncanon_rate"),
            ("n_tso_concat","median_tso_concatemer_rate"),
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
        "broad_noncanonical_read_frac": round(sm.broad_noncanonical_read_frac, 4),
        "artifact_candidate_read_frac": round(sm.artifact_candidate_read_frac, 4),
        "broad_noncanonical_base_frac": round(sm.broad_noncanonical_base_frac, 4),
        "artifact_candidate_base_frac": round(sm.artifact_candidate_base_frac, 4),
        "strand_concordance":       round(sm.strand_concordance,       4),
        "chimeric_read_frac":       round(sm.chimeric_read_frac,       4),
        "multimapper_read_frac":    round(sm.multimapper_read_frac,    4),
        "unmapped_read_frac":       round(sm.unmapped_read_frac,       4),
        "low_mapq_read_frac":       round(sm.low_mapq_read_frac,       4),
        "n_cells":                  sm.n_cells,
        "n_reads_classified":       sm.n_reads_classified,
        "n_records_total":          sm.n_records_total,
        "n_alignments_classified":  sm.n_alignments_classified,
        "per_cell_noise_median":    round(sm.per_cell_noise_median,    4),
        "per_cell_noise_iqr":       round(sm.per_cell_noise_iqr,       4),
        "tso_invasion_frac":        round(sm.n_tso_invasion    / (sm.n_reads_classified or 1), 4),
        "polya_priming_frac":       round(sm.n_polya_priming   / (sm.n_reads_classified or 1), 4),
        "tso_concatemer_frac":      round(sm.n_tso_concatemer  / (sm.n_reads_classified or 1), 4),
        "discordant_pair_frac":     round(sm.n_discordant_pair / (sm.n_reads_classified or 1), 4),
    }
    # Add per-category read fractions
    for cat_name, frac in sm.read_fracs.items():
        data[f"read_frac_{cat_name}"] = round(frac, 4)

    if sm.full_length_read_frac is not None:
        data["full_length_read_frac"] = round(sm.full_length_read_frac, 4)
    if sm.three_prime_anchored_frac is not None:
        data["three_prime_anchored_frac"] = round(sm.three_prime_anchored_frac, 4)
    if sm.five_prime_anchored_frac is not None:
        data["five_prime_anchored_frac"] = round(sm.five_prime_anchored_frac, 4)
    if sm.both_ends_anchored_frac is not None:
        data["both_ends_anchored_frac"] = round(sm.both_ends_anchored_frac, 4)
    if sm.mean_mapq is not None:
        data["mean_mapq"] = round(sm.mean_mapq, 4)
    if sm.mean_edit_distance is not None:
        data["mean_edit_distance"] = round(sm.mean_edit_distance, 4)
    if sm.softclip_base_frac is not None:
        data["softclip_base_frac"] = round(sm.softclip_base_frac, 4)

    return {"id": "scnoisemeter", "data": {sm.sample_name: data}}
