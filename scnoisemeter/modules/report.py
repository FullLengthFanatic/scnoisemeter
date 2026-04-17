"""
report.py
=========
Generates interactive HTML reports using Plotly.

Single-sample report (run mode)
--------------------------------
  - Noise profile donut chart (read fractions by category)
  - Noise profile bar chart (base fractions by category)
  - Per-cell noise distribution (violin + strip plot)
  - Read length distributions per category (overlapping histograms)
  - Artifact flag summary (TSO, polyA, non-canonical junctions)
  - Sample metadata table
  - Warnings panel

Comparison report (compare mode)
---------------------------------
  - Side-by-side noise profile bars (A vs B)
  - Delta plot showing change per category with significance indicators
  - Per-cell noise distribution shift (paired violin)
  - Read length distribution overlays
  - Statistics table from comparison.stats.tsv

All plots are written into a single self-contained HTML file.
Plotly JS is loaded from CDN — requires internet access to render.
For fully offline reports, set offline=True (embeds the JS, ~3 MB).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

import pandas as pd
import plotly.io as pio

from scnoisemeter.constants import (
    AMBIGUOUS_CATEGORIES,
    NOISE_CATEGORIES,
)
from scnoisemeter.modules.metrics import CellTable, SampleMetrics
from scnoisemeter.modules.report_figures import (
    CATEGORY_COLOURS,
    CATEGORY_CRITERIA,
    CATEGORY_LABELS,
    _artifact_flags,
    _base_fraction_bars,
    _category_legend,
    _cluster_heatmap,
    _cluster_noise_plot,
    _comparison_bars,
    _comparison_lengths,
    _comparison_violin,
    _delta_plot,
    _fraction_bar,
    _insert_size_distribution,
    _intergenic_loci_plots,
    _length_distributions,
    _length_stratified_chart,
    _noise_bars,
    _noise_donut,
    _per_cell_violin,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def write_run_report(
    sm: SampleMetrics,
    ct: CellTable,
    length_samples: dict,
    output_path: Path,
    *,
    cluster_df: "Optional[pd.DataFrame]" = None,
    intergenic_loci: "Optional[list]" = None,
    length_stratified: "Optional[pd.DataFrame]" = None,
    platform: Optional[str] = None,
    insert_sizes: "Optional[dict]" = None,
    offline: bool = False,
) -> None:
    """
    Write a single-sample HTML report.

    Parameters
    ----------
    sm : SampleMetrics
    ct : CellTable
    length_samples : dict[ReadCategory, list[int]]
        Reservoir-sampled read lengths per category.
    output_path : Path
        Destination HTML file.
    cluster_df : pd.DataFrame, optional
        Per-cluster noise metrics from compute_cluster_metrics().
        When provided, adds cluster decomposition plots to the report.
    length_stratified : pd.DataFrame, optional
        Output of compute_length_stratification().  When provided, a stacked
        horizontal bar chart ("Noise by read length") is added after the
        classification bar charts.
    offline : bool
        If True, embed Plotly JS in the HTML (no internet required).
    """
    _illumina_platforms = {"illumina", "illumina_10x", "illumina_bd"}
    is_illumina = bool(platform and platform.lower() in _illumina_platforms)

    figures = []

    figures.append(_noise_donut(sm))
    figures.append(_category_legend())
    figures.append(_noise_bars(sm))
    figures.append(_base_fraction_bars(sm))

    # Length-stratified bar chart and read-length histograms are suppressed for
    # Illumina short-read data (all reads fall into the <150 bp bin, making
    # both panels uninformative).  Insert size distribution is shown instead.
    if not is_illumina:
        if length_stratified is not None and not length_stratified.empty:
            figures.append(_length_stratified_chart(length_stratified))

    if not ct.df.empty:
        figures.append(_per_cell_violin(ct))

    if not is_illumina and length_samples:
        figures.append(_length_distributions(length_samples))

    figures.append(_artifact_flags(sm))

    # Insert size distribution for Illumina paired-end data
    report_warnings = list(sm.warnings)
    if is_illumina:
        has_signal = bool(insert_sizes and insert_sizes.get("signal"))
        has_noise  = bool(insert_sizes and insert_sizes.get("noise"))
        if has_signal or has_noise:
            figures.append(_insert_size_distribution(insert_sizes))
        else:
            report_warnings.append(
                "Illumina BAM: no properly paired reads found in the sample — "
                "insert size distribution plot was not generated."
            )

    # Intergenic loci plots (only when profiler ran and found significant loci)
    if intergenic_loci:
        for fig_ig in _intergenic_loci_plots(intergenic_loci):
            if fig_ig.data:
                figures.append(fig_ig)

    # Per-cluster plots (only when obs metadata was provided)
    if cluster_df is not None and not cluster_df.empty:
        figures.append(_cluster_noise_plot(cluster_df))
        figures.append(_cluster_heatmap(cluster_df))

    # For Illumina short-read data, TSS and polyA anchoring metrics reflect
    # read-end proximity only (reads rarely span full transcripts at 50–150 bp).
    if is_illumina:
        report_warnings.append(
            "Illumina short-read data: TSS-anchored and polyA-anchored "
            "full-length metrics reflect read-end proximity only.  Because "
            "reads are 50–150 bp, they rarely span an entire transcript.  "
            "These values should NOT be compared directly with long-read "
            "TSS/polyA fractions; they indicate strand-correct 5\u2032 or 3\u2032 "
            "proximity, not full-length capture efficiency."
        )

    html = _assemble_html(
        figures=figures,
        title=f"scNoiseMeter — {sm.sample_name}",
        metadata_table=_metadata_table(sm),
        warnings=report_warnings,
        offline=offline,
    )

    output_path.write_text(html, encoding="utf-8")
    logger.info("Wrote report: %s", output_path)


def write_compare_report(
    sm_a: SampleMetrics,
    sm_b: SampleMetrics,
    ct_a: CellTable,
    ct_b: CellTable,
    length_samples_a: dict,
    length_samples_b: dict,
    stats_df: Optional[pd.DataFrame],
    output_path: Path,
    *,
    offline: bool = False,
) -> None:
    """Write a two-sample comparison HTML report."""
    figures = []

    figures.append(_comparison_bars(sm_a, sm_b))
    figures.append(_delta_plot(sm_a, sm_b, stats_df))

    warnings = list(set(sm_a.warnings + sm_b.warnings))

    # Only render per-cell comparison violin when both samples have real per-cell
    # data (n_cells > 1).  n_cells == 1 is the barcode-agnostic sentinel.
    if sm_a.n_cells > 1 and sm_b.n_cells > 1:
        if not ct_a.df.empty and not ct_b.df.empty:
            figures.append(_comparison_violin(ct_a, ct_b))
    else:
        warnings.append(
            "Per-cell comparison not shown: one or both samples are barcode-agnostic."
        )

    if length_samples_a and length_samples_b:
        figures.append(_comparison_lengths(length_samples_a, length_samples_b,
                                            sm_a.sample_name, sm_b.sample_name))

    html = _assemble_html(
        figures=figures,
        title=f"scNoiseMeter — {sm_a.sample_name} vs {sm_b.sample_name}",
        metadata_table=_comparison_metadata_table(sm_a, sm_b),
        warnings=warnings,
        offline=offline,
    )

    output_path.write_text(html, encoding="utf-8")
    logger.info("Wrote comparison report: %s", output_path)


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

def _metadata_table(sm: SampleMetrics) -> str:
    rows = [
        ("Sample name",     sm.sample_name),
        ("BAM path",        sm.bam_path),
        ("Platform",        sm.platform),
        ("Pipeline stage",  sm.pipeline_stage),
        ("Aligner",         sm.aligner or "unknown"),
        ("Total reads",     f"{sm.n_reads_total:,}"),
        ("Classified reads", f"{sm.n_reads_classified:,}"),
        ("Cells detected",
         f"{sm.n_cells:,}" if sm.n_cells > 1
         else "N/A — barcode-agnostic mode (no CB tags in BAM)"),
        ("Noise — conservative (reads)",
         f"{sm.noise_read_frac:.2%}  "
         "(includes intronic pure/boundary — upper bound; may contain genuine pre-mRNA)"),
        ("Noise — strict (reads)",
         f"{sm.noise_read_frac_strict:.2%}  "
         "(unambiguous RT/PCR artifacts only — lower bound)"),
        ("Noise (bases)",   f"{sm.noise_base_frac:.2%}"),
        ("Strand concordance", f"{sm.strand_concordance:.2%}"),
        ("Chimeric rate",   f"{sm.chimeric_read_frac:.2%}"),
    ]
    if sm.full_length_read_frac is not None:
        if getattr(sm, "_polya_sites_used", False):
            label = "3′-end at annotated polyA site"
            tooltip = "(fraction of exonic-sense reads ending within 50 bp of a known polyA site)"
        else:
            label = "Full-length read fraction"
            tooltip = "(length proxy — provide --polya-sites for the polyA-anchored metric)"
        rows.append((label, f"{sm.full_length_read_frac:.2%} {tooltip}"))

    if getattr(sm, "tss_anchored_frac", None) is not None:
        rows.append((
            "5′-end at annotated TSS",
            f"{sm.tss_anchored_frac:.2%}  "
            "(fraction of exonic-sense reads starting within 100 bp of a CAGE peak)",
        ))

    if getattr(sm, "numt_read_frac", None) is not None:
        rows.append((
            "NUMT intervals loaded",
            f"{int(sm.numt_read_frac):,}  "
            "(nuclear mitochondrial segments; full NUMT disambiguation requires dual-alignment)",
        ))

    # Cell barcodes row (present only when --cell-barcodes was used)
    cell_bc_info = getattr(sm, "_cell_barcodes_info", None)
    if cell_bc_info:
        import os as _os
        rows.append((
            "Cell barcodes",
            f"{_os.path.basename(cell_bc_info['path'])} · {cell_bc_info['n_barcodes']:,} barcodes loaded",
        ))

    # Annotation provenance rows — always shown regardless of resolution path
    import os
    gtf_info   = getattr(sm, "_gtf_info",   None)
    polya_info = getattr(sm, "_polya_info",  None)

    gtf_version = gtf_info.get("version") if gtf_info else None
    gtf_source  = gtf_info.get("source", "user-supplied") if gtf_info else None
    gtf_path    = gtf_info.get("path") if gtf_info else None

    if gtf_version:
        rows.append(("GTF annotation", f"GENCODE v{gtf_version} ({gtf_source})"))
    elif gtf_path:
        rows.append(("GTF annotation",
                     f"{os.path.basename(gtf_path)} ({gtf_source})"))
    else:
        rows.append(("GTF annotation", "not specified"))

    polya_version = polya_info.get("version") if polya_info else None
    polya_source  = polya_info.get("source", "user-supplied") if polya_info else None
    polya_path    = polya_info.get("path") if polya_info else None
    polya_db      = polya_info.get("db", "polyasite3") if polya_info else None

    mismatch = ""
    if (
        gtf_version and polya_version
        and abs(gtf_version - polya_version) > 5
    ):
        mismatch = " &#9888; version mismatch"

    if polya_db == "both" and polya_version:
        rows.append(("PolyA atlas",
                     f"PolyASite 3.0 GENCODE v{polya_version} + PolyA_DB v4 ({polya_source}{mismatch})"))
    elif polya_db == "both":
        rows.append(("PolyA atlas", f"PolyASite 3.0 + PolyA_DB v4 ({polya_source})"))
    elif polya_db == "polyadb4":
        rows.append(("PolyA atlas", f"PolyA_DB v4 ({polya_source})"))
    elif polya_version:
        rows.append(("PolyA atlas",
                     f"PolyASite 3.0 GENCODE v{polya_version} ({polya_source}{mismatch})"))
    elif polya_path:
        rows.append(("PolyA atlas",
                     f"{os.path.basename(polya_path)} ({polya_source}{mismatch})"))
    else:
        rows.append(("PolyA atlas", "not specified"))

    tss_info = getattr(sm, "_tss_info", None)
    if tss_info:
        tss_db_val  = tss_info.get("db", "fantom5")
        tss_source  = tss_info.get("source", "user-supplied")
        tss_path    = tss_info.get("path")
        if tss_db_val == "fantom5":
            rows.append(("TSS / CAGE atlas", f"FANTOM5 CAGE peaks ({tss_source})"))
        elif tss_path:
            rows.append(("TSS / CAGE atlas",
                         f"{os.path.basename(tss_path)} ({tss_source})"))
        else:
            rows.append(("TSS / CAGE atlas", f"{tss_db_val} ({tss_source})"))

    inner = "".join(
        f"<tr><td><b>{k}</b></td><td>{v}</td></tr>" for k, v in rows
    )
    return f"<table class='meta-table'>{inner}</table>"


def _comparison_metadata_table(sm_a: SampleMetrics, sm_b: SampleMetrics) -> str:
    keys = [
        ("Platform",        "platform"),
        ("Pipeline stage",  "pipeline_stage"),
        ("Total reads",     "n_reads_total"),
        ("Cells detected",  "n_cells"),
        ("Noise (reads)",   "noise_read_frac"),
        ("Noise (bases)",   "noise_base_frac"),
        ("Strand concordance", "strand_concordance"),
        ("Chimeric rate",   "chimeric_read_frac"),
    ]
    header = f"<tr><th>Metric</th><th>{sm_a.sample_name}</th><th>{sm_b.sample_name}</th></tr>"
    rows = []
    for label, attr in keys:
        va = getattr(sm_a, attr, "")
        vb = getattr(sm_b, attr, "")
        if isinstance(va, float):
            va = f"{va:.2%}" if attr.endswith("frac") or attr == "strand_concordance" else f"{va:.4f}"
        if isinstance(vb, float):
            vb = f"{vb:.2%}" if attr.endswith("frac") or attr == "strand_concordance" else f"{vb:.4f}"
        if isinstance(va, int):
            va = f"{va:,}"
        if isinstance(vb, int):
            vb = f"{vb:,}"
        rows.append(f"<tr><td><b>{label}</b></td><td>{va}</td><td>{vb}</td></tr>")
    return f"<table class='meta-table'>{header}{''.join(rows)}</table>"


def _assemble_html(
    figures: list,
    title: str,
    metadata_table: str,
    warnings: list,
    offline: bool,
) -> str:
    include_plotlyjs = "cdn"
    if offline:
        include_plotlyjs = True

    fig_divs = []
    for fig in figures:
        div = pio.to_html(fig, full_html=False, include_plotlyjs=False)
        fig_divs.append(f"<div class='plot-card'>{div}</div>")

    warning_html = ""
    if warnings:
        items = "".join(f"<li>{w}</li>" for w in warnings)
        warning_html = f"""
        <div class="warnings-panel">
          <h3>⚠ Warnings</h3>
          <ul>{items}</ul>
        </div>"""

    plotly_script = (
        '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'
        if not offline else ""
    )
    if offline:
        import plotly
        plotly_js = (Path(plotly.__file__).parent / "package_data" / "plotly.min.js").read_text()
        plotly_script = f"<script>{plotly_js}</script>"

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{title}</title>
  {plotly_script}
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
            background: #f5f6fa; color: #2d3436; }}
    header {{ background: #2d3436; color: white; padding: 24px 32px; }}
    header h1 {{ font-size: 1.5rem; font-weight: 600; }}
    header p  {{ font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }}
    .container {{ max-width: 1400px; margin: 0 auto; padding: 32px 24px; }}
    .section-title {{ font-size: 1.1rem; font-weight: 600; color: #636e72;
                      text-transform: uppercase; letter-spacing: 0.05em;
                      margin: 32px 0 12px; }}
    .plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(560px, 1fr));
                  gap: 20px; }}
    .plot-card {{ background: white; border-radius: 8px;
                  box-shadow: 0 1px 4px rgba(0,0,0,.08);
                  padding: 16px; overflow: hidden; }}
    .meta-table {{ width: 100%; border-collapse: collapse; font-size: 0.9rem; }}
    .meta-table td, .meta-table th {{ padding: 6px 12px; border-bottom: 1px solid #eee; }}
    .meta-table tr:hover td {{ background: #f9f9f9; }}
    .warnings-panel {{ background: #fff3cd; border: 1px solid #ffc107;
                       border-radius: 6px; padding: 16px 20px; margin-bottom: 24px; }}
    .warnings-panel h3 {{ margin-bottom: 8px; color: #856404; }}
    .warnings-panel ul {{ margin-left: 20px; font-size: 0.88rem; color: #533f03; }}
    footer {{ text-align: center; padding: 24px; font-size: 0.8rem; color: #b2bec3; }}
  </style>
</head>
<body>
  <header>
    <h1>scNoiseMeter</h1>
    <p>{title}</p>
  </header>
  <div class="container">
    {warning_html}
    <div class="section-title">Sample metadata</div>
    <div class="plot-card" style="margin-bottom:20px">{metadata_table}</div>
    <div class="section-title">Classification results</div>
    <div class="plot-grid">
      {''.join(fig_divs)}
    </div>
  </div>
  <footer>
    Generated by scNoiseMeter &nbsp;·&nbsp; github.com/scnoisemeter
  </footer>
</body>
</html>"""
