"""
report_figures.py
=================
Plotly figure builders and category metadata for the scNoiseMeter HTML reports.

Separated from report.py so each module stays under 500 lines.
Public API entry points (write_run_report, write_compare_report) live in report.py.
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from scnoisemeter.constants import (
    ADAPTIVE_PVALUE_THRESHOLD,
    CATEGORY_ORDER,
    LENGTH_BIN_LABELS_LONG,
    LENGTH_BIN_LABELS_SHORT,
    ReadCategory,
)
from scnoisemeter.modules.metrics import CellTable, SampleMetrics

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Colour palette — one colour per category, consistent across all plots
# ---------------------------------------------------------------------------

CATEGORY_COLOURS = {
    ReadCategory.EXONIC_SENSE:        "#2ecc71",   # green  — signal
    ReadCategory.EXONIC_ANTISENSE:    "#e74c3c",   # red    — opposite orientation
    ReadCategory.INTRONIC_JXNSPAN:    "#f39c12",   # amber  — ambiguous
    ReadCategory.INTRONIC_PURE:       "#e67e22",   # orange — intronic
    ReadCategory.INTRONIC_BOUNDARY:   "#d35400",   # dark orange
    ReadCategory.INTERGENIC_SPARSE:   "#95a5a6",   # grey   — not enriched
    ReadCategory.INTERGENIC_REPEAT:   "#7f8c8d",   # dark grey
    ReadCategory.INTERGENIC_HOTSPOT:  "#c0392b",   # dark red — artifact
    ReadCategory.INTERGENIC_NOVEL:    "#8e44ad",   # purple — candidate biology
    ReadCategory.INTERGENIC_ENRICHED: "#7f8c8d",   # grey — unresolved enrichment
    ReadCategory.CHIMERIC:            "#2980b9",   # blue   — chimeric
    ReadCategory.MITOCHONDRIAL:       "#1abc9c",   # teal
    ReadCategory.MULTIMAPPER:         "#bdc3c7",   # light grey
    ReadCategory.AMBIGUOUS:           "#ecf0f1",   # very light grey
    ReadCategory.AMBIGUOUS_COD_NCOD:  "#d2b4de",   # light purple — coding/noncoding overlap
    ReadCategory.AMBIGUOUS_COD_COD:   "#9b59b6",   # purple — coding/coding overlap
    ReadCategory.UNASSIGNED:          "#f0f0f0",   # near-white
}

# Human-readable labels for plot axes and legends
CATEGORY_LABELS = {
    ReadCategory.EXONIC_SENSE:        "Exonic sense",
    ReadCategory.EXONIC_ANTISENSE:    "Exonic antisense",
    ReadCategory.INTRONIC_JXNSPAN:    "Intronic junction-spanning",
    ReadCategory.INTRONIC_PURE:       "Intronic pure",
    ReadCategory.INTRONIC_BOUNDARY:   "Intronic boundary",
    ReadCategory.INTERGENIC_SPARSE:   "Intergenic sparse",
    ReadCategory.INTERGENIC_REPEAT:   "Intergenic repeat",
    ReadCategory.INTERGENIC_HOTSPOT:  "Intergenic hotspot",
    ReadCategory.INTERGENIC_NOVEL:    "Intergenic novel",
    ReadCategory.INTERGENIC_ENRICHED: "Intergenic enriched, unresolved",
    ReadCategory.CHIMERIC:            "Chimeric",
    ReadCategory.MITOCHONDRIAL:       "Mitochondrial",
    ReadCategory.MULTIMAPPER:         "Multi-mapper",
    ReadCategory.AMBIGUOUS:           "Ambiguous",
    ReadCategory.AMBIGUOUS_COD_NCOD:  "Ambiguous coding/noncoding",
    ReadCategory.AMBIGUOUS_COD_COD:   "Ambiguous coding/coding",
    ReadCategory.UNASSIGNED:          "Unassigned",
}

# One-line criterion for each category — shown as hover text in plots
# and in the legend table in the HTML report.
CATEGORY_CRITERIA: dict[ReadCategory, str] = {
    ReadCategory.UNMAPPED:            "Did not align to the reference genome.",
    ReadCategory.SECONDARY:           "SAM flag 0x100 — duplicate multi-mapper record, skipped.",
    ReadCategory.SUPPLEMENTARY:       "SAM flag 0x800 — split alignment partner, passed to chimeric detector only.",
    ReadCategory.MULTIMAPPER:         "Explicit multi-hit evidence: NH tag > 1, or non-empty XA when NH is absent. MAPQ 0 alone is reported separately as low MAPQ.",
    ReadCategory.CHIMERIC:            "SA evidence is inter-chromosomal, strand-discordant, or incompatible in query/genomic order; extreme paired insert size is a fallback.",
    ReadCategory.MITOCHONDRIAL:       "Maps to the mitochondrial chromosome (chrM / MT).",
    ReadCategory.EXONIC_SENSE:        "≥1 aligned base overlaps an annotated exon on the same strand as the read.",
    ReadCategory.EXONIC_ANTISENSE:    "≥1 aligned base overlaps an annotated exon on the opposite strand; interpretation depends on protocol strandedness and biology.",
    ReadCategory.INTRONIC_JXNSPAN:    "Intronic bases on an alignment with a CIGAR N; exact junction and motif evidence are checked independently.",
    ReadCategory.INTRONIC_BOUNDARY:   "Spans an exon–intron boundary (sense) with no CIGAR N at the junction — candidate incomplete reverse transcription.",
    ReadCategory.INTRONIC_PURE:       "Intronic alignment without CIGAR N or an exon–intron boundary; can represent pre-mRNA or technical capture.",
    ReadCategory.INTERGENIC_REPEAT:   "Outside gene bodies with ≥50% sampled aligned span overlapping supplied repeat intervals.",
    ReadCategory.INTERGENIC_HOTSPOT:  "Enriched fixed window, monoexonic, with strand-aware genomic A/T-run context and no nearby matched polyA site — internal-priming candidate.",
    ReadCategory.INTERGENIC_NOVEL:    "Intergenic, above significance threshold, strand-consistent, multi-barcode, AND splice or polyA evidence — candidate unannotated transcript.",
    ReadCategory.INTERGENIC_ENRICHED: "Fixed intergenic window with significant enrichment but insufficient positive evidence for either internal priming or a transcript candidate.",
    ReadCategory.INTERGENIC_SPARSE:   "Intergenic window not promoted by the global-rate enrichment/support screen.",
    ReadCategory.AMBIGUOUS:           "Maps to a region where multiple genes overlap and the sub-type (cod/cod or cod/ncod) cannot be determined.",
    ReadCategory.AMBIGUOUS_COD_COD:   "Maps to a same-strand region shared by two protein-coding gene bodies.",
    ReadCategory.AMBIGUOUS_COD_NCOD:  "Maps to a same-strand region shared by coding and non-coding gene bodies.",
    ReadCategory.UNASSIGNED:          "CB tag absent when a whitelist was provided, or CB tag not in the whitelist.",
}


# ---------------------------------------------------------------------------
# Figure builders
# ---------------------------------------------------------------------------

def _noise_donut(sm: SampleMetrics) -> go.Figure:
    """
    Summary overview: key scalar metrics plus non-exonic category bars.

    When exonic sense > 80%, a standard pie collapses other categories into an
    illegible sliver.  Instead we show:
      - Key composition numbers as annotations
      - A horizontal bar chart of ONLY the non-exonic-sense, non-unassigned
        categories, giving each category its own full-width bar.
    """
    exclude = {
        ReadCategory.EXONIC_SENSE, ReadCategory.UNASSIGNED,
        ReadCategory.SECONDARY, ReadCategory.UNMAPPED,
        ReadCategory.SUPPLEMENTARY,
    }
    cats   = [c for c in CATEGORY_ORDER
              if c not in exclude and sm.read_fracs.get(c.value, 0) > 0.0005]
    values = [sm.read_fracs[c.value] for c in cats]
    labels = [CATEGORY_LABELS[c] for c in cats]
    colors = [CATEGORY_COLOURS[c] for c in cats]

    es_frac        = sm.read_fracs.get(ReadCategory.EXONIC_SENSE.value, 0)
    noise_frac = getattr(sm, "broad_noncanonical_read_frac", sm.noise_read_frac) or 0
    noise_frac_strict = getattr(sm, "artifact_candidate_read_frac", sm.noise_read_frac_strict)
    strand_conc    = sm.strand_concordance
    max_val = max(values) if values else 0.01

    fig = go.Figure()

    # Horizontal bars for non-exonic/ambiguous categories
    criteria = [CATEGORY_CRITERIA.get(c, "") for c in cats]
    hover = [
        f"<b>{lbl}</b><br>{v:.2%}<br><i>{crit}</i>"
        for lbl, v, crit in zip(labels, values, criteria)
    ]
    fig.add_trace(go.Bar(
        x=values,
        y=labels,
        orientation="h",
        marker_color=colors,
        text=None,           # no text labels — values visible on hover and in the bar chart below
        hovertext=hover,
        hoverinfo="text",
        width=0.6,
    ))

    fig.update_layout(
        title=dict(
            text=(
                f"Read classification overview — "
                f"Exonic sense: <b>{es_frac:.1%}</b>  |  "
                f"Broad non-canonical: <b>{noise_frac:.1%}</b>  |  "
                f"Artifact candidates: <b>{noise_frac_strict:.1%}</b>  |  "
                f"Strand concordance: <b>{strand_conc:.1%}</b>"
                if strand_conc else
                f"Read classification overview — "
                f"Exonic sense: <b>{es_frac:.1%}</b>  |  "
                f"Broad non-canonical: <b>{noise_frac:.1%}</b>  |  "
                f"Artifact candidates: <b>{noise_frac_strict:.1%}</b>"
            ),
            x=0.5, font=dict(size=12),
        ),
        xaxis=dict(
            title="Read fraction",
            tickformat=".1%",
            range=[0, max_val * 1.05],
        ),
        yaxis=dict(autorange="reversed", automargin=True),
        showlegend=False,
        margin=dict(t=70, b=40, l=220, r=20),
        height=max(300, len(cats) * 38 + 100),
    )
    return fig


def _noise_bars(sm: SampleMetrics) -> go.Figure:
    """All-categories read fraction bar (context view)."""
    return _fraction_bar(sm.read_fracs, title="Read fraction — all categories")


def _base_fraction_bars(sm: SampleMetrics) -> go.Figure:
    """All-categories base fraction bar (context view)."""
    return _fraction_bar(sm.base_fracs, title="Base fraction — all categories")


def _fraction_bar(fracs: dict, title: str) -> go.Figure:
    """Simple horizontal bar chart, one bar per category with a value."""
    cats   = [c for c in CATEGORY_ORDER if fracs.get(c.value, 0) > 0.0001]
    values = [fracs[c.value] for c in cats]
    labels = [CATEGORY_LABELS[c] for c in cats]
    colors = [CATEGORY_COLOURS[c] for c in cats]

    if not values:
        return go.Figure()

    max_val = max(values)
    max_label_len = max(len(label) for label in labels)
    left_margin = max(220, max_label_len * 8)

    criteria = [CATEGORY_CRITERIA.get(c, "") for c in cats]
    hover = [
        f"<b>{lbl}</b><br>{v:.2%}<br><i>{crit}</i>"
        for lbl, v, crit in zip(labels, values, criteria)
    ]
    fig = go.Figure(go.Bar(
        x=values,
        y=labels,
        orientation="h",
        marker_color=colors,
        text=[f"{v:.2%}" for v in values],
        textposition=["inside" if v > max_val * 0.55 else "outside" for v in values],
        hovertext=hover,
        hoverinfo="text",
    ))
    fig.update_layout(
        title=dict(text=title, x=0.5),
        xaxis=dict(
            title="Fraction",
            tickformat=".1%",
            range=[0, min(max_val * 1.22, 1.0)],
        ),
        yaxis=dict(autorange="reversed", automargin=True),
        margin=dict(t=60, b=40, l=left_margin, r=80),
        height=max(300, len(cats) * 36 + 100),
        uniformtext=dict(minsize=9, mode="hide"),
    )
    return fig


def _length_stratified_chart(strat_df: "pd.DataFrame") -> go.Figure:
    """
    Stacked horizontal bar chart — one row per read-length bin.

    Each bar is coloured by category (using the standard CATEGORY_COLOURS
    palette) and shows the within-bin fraction of reads.  Total read count
    for each bin is annotated to the right of the bar.
    """
    if strat_df.empty:
        return go.Figure()

    # Determine which bin labels are present and their display order
    present_bins = set(strat_df["length_bin"].unique())
    bin_order = [b for b in LENGTH_BIN_LABELS_LONG if b in present_bins]
    # Fall back to SHORT labels if none of the LONG labels matched
    if not bin_order:
        bin_order = [b for b in LENGTH_BIN_LABELS_SHORT if b in present_bins]
    # Any remaining bins not in either standard list
    bin_order += [b for b in strat_df["length_bin"].unique() if b not in bin_order]

    fig = go.Figure()

    for cat in CATEGORY_ORDER:
        cat_str = cat.value
        cat_df  = strat_df[strat_df["category"] == cat_str]
        if cat_df.empty:
            continue

        label = CATEGORY_LABELS.get(cat, cat_str)
        color = CATEGORY_COLOURS.get(cat, "#cccccc")

        x_vals, y_vals, hover_texts = [], [], []
        for bin_lbl in bin_order:
            row = cat_df[cat_df["length_bin"] == bin_lbl]
            frac = float(row["fraction_of_bin"].iloc[0]) if not row.empty else 0.0
            n    = int(row["count"].iloc[0])             if not row.empty else 0
            x_vals.append(frac)
            y_vals.append(bin_lbl)
            hover_texts.append(
                f"<b>{label}</b><br>Bin: {bin_lbl}<br>"
                f"Count: {n:,}<br>Fraction of bin: {frac:.2%}"
            )

        fig.add_trace(go.Bar(
            x=x_vals,
            y=y_vals,
            name=label,
            orientation="h",
            marker_color=color,
            hovertext=hover_texts,
            hoverinfo="text",
        ))

    # Annotate total read count per bin on the right
    for bin_lbl in bin_order:
        bin_total = int(strat_df[strat_df["length_bin"] == bin_lbl]["count"].sum())
        fig.add_annotation(
            x=1.01,
            y=bin_lbl,
            text=f"{bin_total:,}",
            xref="paper",
            yref="y",
            showarrow=False,
            xanchor="left",
            font=dict(size=11),
        )

    # Layout: title sits above the legend, legend sits above the chart.
    # title at y=1.20 (yanchor="bottom") → title bottom is above legend top.
    # legend at y=1.05 (yanchor="bottom") → legend bottom clears the chart area.
    # margin t=170 provides ~170px for both title and legend in the top margin,
    # giving ≥ 20px gap between title and legend, and ≥ 15px between legend
    # and chart top — well above the 8 px / 12 px minimums required.
    fig.update_layout(
        barmode="stack",
        title=dict(text="Broad non-canonical composition by read length", x=0.5, y=0.98, yanchor="top"),
        xaxis=dict(title="Fraction of bin", tickformat=".0%", range=[0, 1.0]),
        yaxis=dict(autorange="reversed", automargin=True),
        legend=dict(
            orientation="h",
            yanchor="bottom", y=1.08,
            xanchor="right",  x=1.0,
            font=dict(size=10),
        ),
        margin=dict(t=160, b=60, l=120, r=100),
        height=max(360, len(bin_order) * 62 + 200),
    )
    return fig


def _per_cell_violin(ct: CellTable) -> go.Figure:
    """Violin + strip plot of the legacy broad-composition column."""
    df = ct.df
    if "noise_read_frac" not in df.columns:
        return go.Figure()

    vals = df["noise_read_frac"].dropna().values

    # With >5 000 cells, rendering individual outlier points causes Plotly to
    # silently drop the violin body.  Disable points for large datasets.
    show_points = "outliers" if len(vals) <= 5_000 else False

    fig = go.Figure()
    fig.add_trace(go.Violin(
        y=vals,
        name="broad non-canonical",
        box_visible=True,
        meanline_visible=True,
        fillcolor="#e74c3c",
        opacity=0.6,
        line_color="#c0392b",
        points=show_points,
        spanmode="hard",
    ))
    fig.update_layout(
        title=dict(text="Per-cell broad non-canonical composition", x=0.5),
        yaxis=dict(title="Read fraction", tickformat=".1%"),
        showlegend=False,
        height=380,
        margin=dict(t=60, b=40, l=70, r=20),
    )
    return fig


def _length_distributions(length_samples: dict) -> go.Figure:
    """
    Read length distributions split into two panels:
      Left  — signal categories (exonic sense, intronic junction-spanning,
               ambiguous) where length reflects transcript biology
      Right — other categories (chimeric, intronic, intergenic, antisense)

    Splitting avoids the problem of exonic-sense (many reads, broad range)
    completely dominating the plot and making all other categories invisible.
    Each panel uses its own y-axis scale.
    """

    SIGNAL_CATS = {
        ReadCategory.EXONIC_SENSE,
        ReadCategory.INTRONIC_JXNSPAN,
        ReadCategory.AMBIGUOUS_COD_COD,
        ReadCategory.AMBIGUOUS_COD_NCOD,
        ReadCategory.MITOCHONDRIAL,
        ReadCategory.MULTIMAPPER,
    }
    NOISE_CATS = {
        ReadCategory.CHIMERIC,
        ReadCategory.EXONIC_ANTISENSE,
        ReadCategory.INTRONIC_PURE,
        ReadCategory.INTRONIC_BOUNDARY,
        ReadCategory.INTERGENIC_SPARSE,
        ReadCategory.INTERGENIC_REPEAT,
        ReadCategory.INTERGENIC_HOTSPOT,
        ReadCategory.INTERGENIC_NOVEL,
        ReadCategory.INTERGENIC_ENRICHED,
    }

    # Only include categories with actual data
    signal_cats = [c for c in CATEGORY_ORDER
                   if c in SIGNAL_CATS and length_samples.get(c)]
    noise_cats  = [c for c in CATEGORY_ORDER
                   if c in NOISE_CATS  and length_samples.get(c)]

    if not signal_cats and not noise_cats:
        return go.Figure()

    # If only one group has data, use a single panel
    if not noise_cats:
        fig = go.Figure()
        for cat in signal_cats:
            fig.add_trace(go.Histogram(
                x=length_samples[cat], name=CATEGORY_LABELS[cat],
                marker_color=CATEGORY_COLOURS[cat],
                opacity=0.7, nbinsx=50, histnorm="percent",
            ))
        fig.update_layout(
            barmode="overlay",
            title=dict(text="Read length distributions by category", x=0.5),
            xaxis=dict(title="Read length (bp)", range=[0, 5000]),
            yaxis=dict(title="% of reads in category"),
            legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.2),
            height=420, margin=dict(t=60, b=120, l=70, r=20),
        )
        return fig

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["Signal categories", "Other categories"],
        horizontal_spacing=0.10,
    )

    for cat in signal_cats:
        fig.add_trace(go.Histogram(
            x=length_samples[cat],
            name=CATEGORY_LABELS[cat],
            marker_color=CATEGORY_COLOURS[cat],
            opacity=0.75, nbinsx=50, histnorm="percent",
            legendgroup="signal",
        ), row=1, col=1)

    for cat in noise_cats:
        fig.add_trace(go.Histogram(
            x=length_samples[cat],
            name=CATEGORY_LABELS[cat],
            marker_color=CATEGORY_COLOURS[cat],
            opacity=0.75, nbinsx=50, histnorm="percent",
            legendgroup="noise",
        ), row=1, col=2)

    fig.update_xaxes(title_text="Read length (bp)", range=[0, 5000], row=1, col=1)
    fig.update_xaxes(title_text="Read length (bp)", range=[0, 5000], row=1, col=2)
    fig.update_yaxes(title_text="% of reads in category", row=1, col=1)
    fig.update_yaxes(title_text="% of reads in category", row=1, col=2)

    fig.update_layout(
        barmode="overlay",
        title=dict(text="Read length distributions by category", x=0.5),
        legend=dict(orientation="h", x=0.5, xanchor="center",
                    y=-0.22, font=dict(size=10)),
        height=460,
        margin=dict(t=70, b=130, l=70, r=20),
    )
    return fig


def _insert_size_distribution(insert_sizes: dict) -> go.Figure:
    """
    Overlapping histograms of fragment insert size for Illumina paired-end data.

    Signal trace (green): reads classified as EXONIC_SENSE.
    Noise trace (red):    all other classified reads.

    Only properly paired reads (SAM flag 0x2) with 0 < |template_length| < 2000
    are included.  Sampled from at most 500,000 read pairs.
    """
    signal_sizes = insert_sizes.get("signal", [])
    noise_sizes  = insert_sizes.get("noise",  [])
    n_sampled    = len(signal_sizes) + len(noise_sizes)

    fig = go.Figure()

    if signal_sizes:
        fig.add_trace(go.Histogram(
            x=signal_sizes,
            name="Signal (exonic sense)",
            marker_color="#2ecc71",
            opacity=0.6,
            nbinsx=100,
            histnorm="percent",
            legendgroup="signal",
        ))

    if noise_sizes:
        fig.add_trace(go.Histogram(
            x=noise_sizes,
            name="Other classified categories",
            marker_color="#e74c3c",
            opacity=0.6,
            nbinsx=100,
            histnorm="percent",
            legendgroup="noise",
        ))

    fig.update_layout(
        barmode="overlay",
        title=dict(
            text=(
                "Insert size distribution (paired-end)"
                f"<br><sup>Sampled from \u2264500,000 read pairs (n\u2009=\u2009{n_sampled:,})</sup>"
            ),
            x=0.5,
        ),
        xaxis=dict(title="Insert size (bp)", range=[0, 1000]),
        yaxis=dict(title="% of reads"),
        legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.15),
        height=420,
        margin=dict(t=80, b=100, l=70, r=20),
    )
    return fig


def _artifact_flags(sm: SampleMetrics) -> go.Figure:
    """Horizontal bar chart of artifact flag rates — labels never cut off."""
    labels = ["TSO invasion", "Internal polyA priming", "Non-canonical junction", "TSO concatemer"]
    values = [sm.n_tso_invasion, sm.n_polya_priming, sm.n_noncanon_junction,
              getattr(sm, "n_tso_concatemer", 0)]
    colors = ["#e74c3c", "#e67e22", "#3498db", "#8e44ad"]
    descriptions = [
        "Reads with TSO sequence (or its reverse complement) in a soft-clip",
        "Reads ending at an A-rich (forward) / T-rich (reverse) reference stretch",
        "Splice junctions not matching GT-AG / GC-AG / AT-AC rule",
        "Reads containing more than one TSO sequence or its reverse complement",
    ]

    denom = max(sm.n_reads_classified, 1)
    fracs = [v / denom for v in values]
    max_frac = max(fracs) if fracs else 0.01

    hover = [
        f"<b>{lbl}</b><br>{f:.3%} of classified reads<br><i>{desc}</i>"
        for lbl, f, desc in zip(labels, fracs, descriptions)
    ]

    fig = go.Figure(go.Bar(
        x=fracs,
        y=labels,
        orientation="h",
        marker_color=colors,
        text=[f"{f:.3%}" for f in fracs],
        textposition=["inside" if f > max_frac * 0.5 else "outside" for f in fracs],
        hovertext=hover,
        hoverinfo="text",
    ))
    fig.update_layout(
        title=dict(text="Artifact flag rates (fraction of classified reads)", x=0.5),
        xaxis=dict(
            title="Fraction",
            tickformat=".2%",
            range=[0, max_frac * 1.05],  # tight — text inside or on bars
        ),
        yaxis=dict(autorange="reversed", automargin=True),
        height=280,
        margin=dict(t=60, b=40, l=200, r=60),
        uniformtext=dict(minsize=9, mode="hide"),
    )
    return fig


# ---------------------------------------------------------------------------
# Category legend table
# ---------------------------------------------------------------------------

def _category_legend() -> go.Figure:
    """
    Static table: category name (coloured) + one-line classification criterion.
    The first column (HTML swatch) is removed — Plotly table cells render the
    raw HTML string rather than interpreting it, so swatches appeared as
    garbled text.  Instead the category name cell is coloured directly.
    """
    show_cats = [
        c for c in CATEGORY_ORDER
        if c not in {ReadCategory.SECONDARY, ReadCategory.SUPPLEMENTARY,
                     ReadCategory.UNMAPPED}
    ]

    names    = [CATEGORY_LABELS.get(c, c.value) for c in show_cats]
    criteria = [CATEGORY_CRITERIA.get(c, "")     for c in show_cats]
    row_bg   = ["#f9f9f9" if i % 2 == 0 else "white" for i in range(len(show_cats))]
    name_colours = [CATEGORY_COLOURS.get(c, "#333") for c in show_cats]

    fig = go.Figure(go.Table(
        columnwidth=[200, 560],
        header=dict(
            values=["<b>Category</b>", "<b>Classification criterion</b>"],
            fill_color="#2c3e50",
            font=dict(color="white", size=11),
            align="left",
            height=30,
        ),
        cells=dict(
            values=[names, criteria],
            fill_color=[row_bg, row_bg],
            font=dict(size=10, color=[name_colours, ["#333"] * len(show_cats)]),
            align="left",
            height=26,
        ),
    ))
    fig.update_layout(
        title=dict(
            text="Category definitions — hover over any bar for details",
            x=0.5, font=dict(size=12, color="#555"),
        ),
        margin=dict(t=50, b=10, l=10, r=10),
        height=max(300, len(show_cats) * 28 + 80),
    )
    return fig


# ---------------------------------------------------------------------------
# Intergenic loci plot
# ---------------------------------------------------------------------------

def _intergenic_loci_plots(intergenic_loci: list) -> "list[go.Figure]":
    """
    Return two figures for intergenic windows:

    1. Horizontal bar chart — top 20 loci by read count (immediately readable).
    2. Scatter plot — adj. p-value (x) vs read count (y, log scale) for all
       loci, with y-jitter so overlapping points spread apart.  Point size is
       proportional to n_barcodes. A dashed line marks the configured adjusted-p threshold.
    """
    if not intergenic_loci:
        return []

    import math
    import random

    from scnoisemeter.constants import ReadCategory

    cat_colours = {
        ReadCategory.INTERGENIC_HOTSPOT: "#c0392b",
        ReadCategory.INTERGENIC_NOVEL:   "#8e44ad",
        ReadCategory.INTERGENIC_ENRICHED:"#7f8c8d",
        ReadCategory.INTERGENIC_REPEAT:  "#7f8c8d",
        ReadCategory.INTERGENIC_SPARSE:  "#bdc3c7",
    }
    cat_labels = {
        ReadCategory.INTERGENIC_HOTSPOT: "Internal priming hotspot",
        ReadCategory.INTERGENIC_NOVEL:   "Candidate unannotated transcript",
        ReadCategory.INTERGENIC_ENRICHED:"Enriched, unresolved",
        ReadCategory.INTERGENIC_REPEAT:  "Repeat-derived",
        ReadCategory.INTERGENIC_SPARSE:  "Sparse (below threshold)",
    }

    # ── Figure 1: top-20 bar chart ────────────────────────────────────────────
    top20 = sorted(intergenic_loci, key=lambda locus: locus.n_reads, reverse=True)[:20]
    bar_labels = [
        f"{locus.contig}:{locus.start:,}-{locus.end:,}" for locus in top20
    ]
    bar_values = [locus.n_reads for locus in top20]
    bar_colours = [
        cat_colours.get(locus.category, "#95a5a6") for locus in top20
    ]
    bar_hover = [
        f"{locus.contig}:{locus.start:,}-{locus.end:,}<br>"
        f"Reads: {locus.n_reads}<br>Barcodes: {locus.n_barcodes}<br>"
        f"Category: {cat_labels.get(locus.category, locus.category.value)}<br>"
        f"Adj. p: {locus.poisson_pvalue_adj:.2e}"
        for locus in top20
    ]

    fig_bar = go.Figure(go.Bar(
        x=bar_values,
        y=bar_labels,
        orientation="h",
        marker_color=bar_colours,
        hovertext=bar_hover,
        hoverinfo="text",
    ))
    fig_bar.update_layout(
        title=dict(text="Top 20 intergenic loci by read count", x=0.5),
        xaxis=dict(title="Read count"),
        yaxis=dict(autorange="reversed", automargin=True),
        height=max(300, len(top20) * 28 + 120),
        margin=dict(t=60, b=60, l=220, r=40),
        showlegend=False,
    )

    # ── Figure 2: scatter — adj. p-value vs read count, jittered ─────────────
    # Compute jitter range: ±5% of the log10 range of read counts
    all_reads = [
        locus.n_reads for locus in intergenic_loci if locus.n_reads > 0
    ]
    if all_reads:
        log_min = math.log10(max(1, min(all_reads)))
        log_max = math.log10(max(all_reads))
        jitter_scale = max(0.05, (log_max - log_min) * 0.05)
    else:
        jitter_scale = 0.1

    rng = random.Random(42)  # deterministic jitter

    by_cat: dict = {}
    for locus in intergenic_loci:
        by_cat.setdefault(locus.category, []).append(locus)

    fig_sc = go.Figure()
    for cat, loci in by_cat.items():
        x_vals, y_vals, sizes, hover = [], [], [], []
        for locus in loci:
            pval = max(locus.poisson_pvalue_adj, 1e-300)  # avoid log(0)
            jitter = rng.uniform(-jitter_scale, jitter_scale)
            log_reads = math.log10(max(1, locus.n_reads))
            x_vals.append(pval)
            y_vals.append(10 ** (log_reads + jitter))
            sizes.append(max(6, min(24, locus.n_barcodes * 2 + 4)))
            hover.append(
                f"{locus.contig}:{locus.start:,}-{locus.end:,}<br>"
                f"Reads: {locus.n_reads}<br>Barcodes: {locus.n_barcodes}<br>"
                f"Adj. p: {locus.poisson_pvalue_adj:.2e}<br>"
                f"polyA downstream: {locus.polya_run_downstream}<br>"
                f"Near polyA site: {locus.near_polya_site}"
            )
        fig_sc.add_trace(go.Scatter(
            x=x_vals,
            y=y_vals,
            mode="markers",
            name=cat_labels.get(cat, cat.value),
            marker=dict(
                color=cat_colours.get(cat, "#95a5a6"),
                size=sizes,
                opacity=0.72,
                line=dict(width=0.5, color="white"),
            ),
            hovertext=hover,
            hoverinfo="text",
        ))

    # Configured adjusted-p threshold.
    fig_sc.add_shape(
        type="line",
        x0=ADAPTIVE_PVALUE_THRESHOLD, x1=ADAPTIVE_PVALUE_THRESHOLD,
        y0=0, y1=1,
        xref="x", yref="paper",
        line=dict(color="#e74c3c", dash="dash", width=1.5),
    )
    fig_sc.add_annotation(
        x=ADAPTIVE_PVALUE_THRESHOLD, y=0.98,
        xref="x", yref="paper",
        text=f"adjusted p = {ADAPTIVE_PVALUE_THRESHOLD:g}",
        showarrow=False,
        font=dict(size=10, color="#e74c3c"),
        xanchor="left",
    )

    fig_sc.update_layout(
        title=dict(
            text="Intergenic windows — adjusted p-value vs read count"
                 "<br><sup>Point size ∝ distinct barcodes · y-axis jittered for readability</sup>",
            x=0.5,
        ),
        xaxis=dict(title="Adjusted p-value", type="log"),
        yaxis=dict(title="Read count (log scale)", type="log"),
        legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.15),
        height=450,
        margin=dict(t=90, b=110, l=70, r=30),
    )

    return [fig_bar, fig_sc]


# ---------------------------------------------------------------------------
# Per-cluster plots
# ---------------------------------------------------------------------------

def _cluster_noise_plot(cluster_df: "pd.DataFrame") -> go.Figure:
    """
    Grouped bar chart showing selected median category fractions per cluster.
    Clusters are sorted by the legacy broad-composition column.
    """
    if cluster_df.empty:
        return go.Figure()

    # Sort clusters by broad non-canonical fraction
    sort_col = (
        "median_broad_noncanonical_read_frac"
        if "median_broad_noncanonical_read_frac" in cluster_df.columns
        else "median_noise_read_frac"
    )
    if sort_col in cluster_df.columns:
        cluster_df = cluster_df.sort_values(sort_col, ascending=False)

    clusters = [str(c) for c in cluster_df.index]

    fig = go.Figure()
    components = [
        ("median_exonic_antisense_frac",  "Exonic antisense",   CATEGORY_COLOURS[ReadCategory.EXONIC_ANTISENSE]),
        ("median_intronic_pure_frac",     "Intronic pure",      CATEGORY_COLOURS[ReadCategory.INTRONIC_PURE]),
        ("median_chimeric_frac",          "Chimeric",           CATEGORY_COLOURS[ReadCategory.CHIMERIC]),
        ("median_ambiguous_cod_cod_frac", "Ambiguous cod/cod",  CATEGORY_COLOURS[ReadCategory.AMBIGUOUS_COD_COD]),
    ]

    for col, label, colour in components:
        if col not in cluster_df.columns:
            continue
        fig.add_trace(go.Bar(
            name=label,
            x=clusters,
            y=cluster_df[col].fillna(0).values,
            marker_color=colour,
            opacity=0.85,
        ))

    fig.update_layout(
        barmode="stack",
        title=dict(text="Per-cluster category decomposition (median read fraction)", x=0.5),
        xaxis=dict(title="Cluster", tickangle=-35),
        yaxis=dict(title="Median broad non-canonical fraction", tickformat=".1%"),
        legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.3),
        height=460,
        margin=dict(t=60, b=140, l=70, r=20),
    )
    return fig


def _cluster_heatmap(cluster_df: "pd.DataFrame") -> go.Figure:
    """
    Heatmap of available category and aggregate fractions across clusters.
    Each cell shows the median fraction for that cluster × category.
    """
    if cluster_df.empty:
        return go.Figure()

    metric_cols = [c for c in cluster_df.columns if c.startswith("median_") and "iqr" not in c]
    if not metric_cols:
        return go.Figure()

    clusters = [str(c) for c in cluster_df.index]
    short_names = [c.replace("median_", "").replace("_frac", "").replace("_", " ") for c in metric_cols]

    z = cluster_df[metric_cols].fillna(0).values.T.tolist()

    fig = go.Figure(go.Heatmap(
        z=z,
        x=clusters,
        y=short_names,
        colorscale="RdYlGn_r",
        text=[[f"{v:.1%}" for v in row] for row in z],
        texttemplate="%{text}",
        textfont=dict(size=9),
        colorbar=dict(title="Fraction", tickformat=".0%"),
    ))
    fig.update_layout(
        title=dict(text="Per-cluster metric heatmap", x=0.5),
        xaxis=dict(title="Cluster", tickangle=-35),
        yaxis=dict(autorange="reversed"),
        height=max(300, len(metric_cols) * 30 + 120),
        margin=dict(t=60, b=100, l=180, r=60),
    )
    return fig


# ---------------------------------------------------------------------------
# Comparison figures
# ---------------------------------------------------------------------------

def _comparison_bars(sm_a: SampleMetrics, sm_b: SampleMetrics) -> go.Figure:
    """
    Horizontal grouped bar chart: read fractions for A vs B per category.
    Horizontal orientation gives category names proper space and avoids
    the diagonal-label cutoff problem of vertical bars.
    Uses split panels (same logic as _fraction_bar) when one sample's
    exonic sense dominates.
    """
    cats = [c for c in CATEGORY_ORDER
            if sm_a.read_fracs.get(c.value, 0) > 0.0005
            or sm_b.read_fracs.get(c.value, 0) > 0.0005]
    labels  = [CATEGORY_LABELS[c] for c in cats]
    vals_a  = [sm_a.read_fracs.get(c.value, 0) for c in cats]
    vals_b  = [sm_b.read_fracs.get(c.value, 0) for c in cats]
    max_val = max(vals_a + vals_b) if (vals_a or vals_b) else 1
    max_label = max((len(label) for label in labels), default=20)
    left_margin = max(260, max_label * 8)

    fig = go.Figure()
    fig.add_trace(go.Bar(
        name=sm_a.sample_name,
        y=labels, x=vals_a,
        orientation="h",
        marker_color="#3498db",
        opacity=0.85,
        text=[f"{v:.2%}" for v in vals_a],
        textposition=["inside" if v > max_val * 0.55 else "outside" for v in vals_a],
        hovertemplate=f"{sm_a.sample_name}: %{{x:.4f}}<extra></extra>",
    ))
    fig.add_trace(go.Bar(
        name=sm_b.sample_name,
        y=labels, x=vals_b,
        orientation="h",
        marker_color="#e74c3c",
        opacity=0.85,
        text=[f"{v:.2%}" for v in vals_b],
        textposition=["inside" if v > max_val * 0.55 else "outside" for v in vals_b],
        hovertemplate=f"{sm_b.sample_name}: %{{x:.4f}}<extra></extra>",
    ))
    fig.update_layout(
        barmode="group",
        title=dict(text="Read fraction by category — comparison", x=0.5),
        xaxis=dict(title="Fraction", tickformat=".1%",
                   range=[0, min(max_val * 1.18, 1.0)]),
        yaxis=dict(autorange="reversed", automargin=True),
        legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.12),
        height=max(360, len(cats) * 44 + 140),
        margin=dict(t=60, b=80, l=left_margin, r=60),
        uniformtext=dict(minsize=9, mode="hide"),
    )
    return fig


def _delta_plot(
    sm_a: SampleMetrics,
    sm_b: SampleMetrics,
    stats_df: Optional[pd.DataFrame],
) -> go.Figure:
    """
    Horizontal descriptive delta plot: B minus A per category.
    Horizontal layout gives long category names the space they need.
    Red bars = B higher than A; green = B lower than A.
    """
    cats = [c for c in CATEGORY_ORDER
            if sm_a.read_fracs.get(c.value, 0) > 0.0005
            or sm_b.read_fracs.get(c.value, 0) > 0.0005]
    labels = [CATEGORY_LABELS[c] for c in cats]
    deltas = [sm_b.read_fracs.get(c.value, 0) - sm_a.read_fracs.get(c.value, 0)
              for c in cats]
    colors = ["#e74c3c" if d > 0 else "#2ecc71" for d in deltas]

    texts = [f"{d:+.2%}" for d in deltas]
    max_label = max((len(label) for label in labels), default=20)
    left_margin = max(260, max_label * 8)
    abs_max = max((abs(d) for d in deltas), default=0.01)

    fig = go.Figure(go.Bar(
        y=labels,
        x=deltas,
        orientation="h",
        marker_color=colors,
        text=texts,
        textposition=["inside" if abs(d) > abs_max * 0.55 else "outside" for d in deltas],
        hovertemplate="%{y}<br>Δ = %{x:.4f}<extra></extra>",
    ))
    fig.add_vline(x=0, line_dash="dash", line_color="black", line_width=1)
    fig.update_layout(
        title=dict(
            text=(
                f"Δ read fraction ({sm_b.sample_name} − {sm_a.sample_name})<br>"
                "<sup>Descriptive composition effect; paired-cell bootstrap "
                "intervals are written to comparison.stats.tsv</sup>"
            ),
            x=0.5,
        ),
        xaxis=dict(
            title="Δ fraction",
            tickformat="+.1%",
            zeroline=True,
            range=[-(abs_max * 1.35), abs_max * 1.35],
        ),
        yaxis=dict(autorange="reversed", automargin=True),
        height=max(360, len(cats) * 38 + 140),
        margin=dict(t=90, b=40, l=left_margin, r=60),
        uniformtext=dict(minsize=9, mode="hide"),
    )
    return fig


def _comparison_violin(ct_a: CellTable, ct_b: CellTable) -> go.Figure:
    """Side-by-side violin of per-cell noise fractions."""
    fig = go.Figure()
    for ct, colour in [(ct_a, "#3498db"), (ct_b, "#e74c3c")]:
        if "noise_read_frac" not in ct.df.columns:
            continue
        vals = ct.df["noise_read_frac"].dropna().values
        fig.add_trace(go.Violin(
            y=vals,
            name=ct.sample_name,
            box_visible=True,
            meanline_visible=True,
            fillcolor=colour,
            opacity=0.6,
            line_color=colour,
            points="outliers",
        ))
    fig.update_layout(
        title=dict(text="Per-cell broad non-canonical composition — comparison", x=0.5),
        yaxis=dict(title="Read fraction", tickformat=".1%"),
        height=400,
        margin=dict(t=60, b=40, l=70, r=20),
    )
    return fig


def _comparison_lengths(
    ls_a: dict, ls_b: dict, label_a: str, label_b: str,
) -> go.Figure:
    """Overlay of exonic-sense read length distributions for A vs B."""
    fig = go.Figure()
    for ls, label, colour in [(ls_a, label_a, "#3498db"), (ls_b, label_b, "#e74c3c")]:
        lengths = ls.get(ReadCategory.EXONIC_SENSE, [])
        if not lengths:
            continue
        fig.add_trace(go.Histogram(
            x=lengths,
            name=f"{label} (exonic-sense)",
            marker_color=colour,
            opacity=0.6,
            nbinsx=60,
            histnorm="percent",
        ))
    fig.update_layout(
        barmode="overlay",
        title=dict(text="Exonic-sense read length distribution — comparison", x=0.5),
        xaxis=dict(title="Read length (bp)", range=[0, 5000]),
        yaxis=dict(title="Percent"),
        height=380,
        margin=dict(t=60, b=60, l=70, r=20),
    )
    return fig
