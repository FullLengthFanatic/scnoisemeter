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
    INTERGENIC_SCATTER_MAX_SPARSE,
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

# Chosen with an OKLab/CVD validator rather than by eye, then checked against
# real sample compositions.  The set passes, on the CATEGORY_ORDER adjacency
# used by every stacked and ordered chart here: OKLCH lightness band, chroma
# floor (nothing reads as gray), colour-vision-deficiency separation
# (worst adjacent pair ΔE 8.1 protan), and the normal-vision floor
# (worst adjacent ΔE 17.6).  Several slots sit below 3:1 contrast on the white
# card, which is why the bar charts carry visible value labels and the report
# ships the category-definition table: identity is never colour alone.
#
# The previous palette failed three of those checks: intergenic_repeat and
# intergenic_enriched were the same hex, ambiguous/unassigned were invisible
# against the card, and intronic_pure/intronic_jxnspan were ΔE 7.8 apart in
# normal vision while sitting adjacent in every stack.
CATEGORY_COLOURS = {
    ReadCategory.EXONIC_SENSE:        "#78b95f",   # green — canonical signal, kept calm
    ReadCategory.EXONIC_ANTISENSE:    "#c50037",   # crimson — opposite orientation
    ReadCategory.INTRONIC_JXNSPAN:    "#e4a067",   # tan
    ReadCategory.INTRONIC_PURE:       "#746b00",   # dark olive
    ReadCategory.INTRONIC_BOUNDARY:   "#c67261",   # clay
    ReadCategory.INTERGENIC_SPARSE:   "#635dfd",   # indigo — not enriched
    ReadCategory.INTERGENIC_REPEAT:   "#00c8d7",   # cyan
    ReadCategory.INTERGENIC_HOTSPOT:  "#f33179",   # magenta — internal-priming candidate
    ReadCategory.INTERGENIC_NOVEL:    "#8b68ad",   # muted purple — candidate biology
    ReadCategory.INTERGENIC_ENRICHED: "#c2ae3e",   # gold — enriched, unresolved
    ReadCategory.CHIMERIC:            "#393eca",   # blue
    ReadCategory.MITOCHONDRIAL:       "#33a182",   # sea green
    ReadCategory.MULTIMAPPER:         "#89401c",   # brown
    ReadCategory.AMBIGUOUS:           "#e384fc",   # orchid
    ReadCategory.AMBIGUOUS_COD_NCOD:  "#8316a4",   # violet — coding/noncoding overlap
    ReadCategory.AMBIGUOUS_COD_COD:   "#9dabf7",   # periwinkle — coding/coding overlap
    ReadCategory.UNASSIGNED:          "#1d81c0",   # steel blue
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
    ReadCategory.MULTIMAPPER:         "Explicit NH tag > 1. MAPQ and producer-specific X* tags are not treated as generic multi-hit counts; low MAPQ is reported separately.",
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
    noise_frac = sm.broad_noncanonical_read_frac or 0
    noise_frac_strict = sm.artifact_candidate_read_frac
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
                # Headline on its own line, numbers as a short subtitle.  All
                # four metrics on the title line overflowed the 560 px report
                # grid cell and the last one was clipped mid-word.  The full
                # set, including strand concordance, is in the metadata table.
                "Read classification overview"
                "<br><sup>"
                f"exonic sense <b>{es_frac:.1%}</b> · "
                f"non-canonical <b>{noise_frac:.1%}</b> · "
                f"artifact <b>{noise_frac_strict:.1%}</b>"
                "</sup>"
            ),
            # Centre on the card, not on the plot area: the 220 px left margin
            # would otherwise push a centred title off the right edge.
            x=0.5, xanchor="center", xref="container", font=dict(size=13),
        ),
        xaxis=dict(
            title="Read fraction",
            tickformat=".1%",
            range=[0, max_val * 1.05],
        ),
        yaxis=dict(autorange="reversed", automargin=True),
        showlegend=False,
        margin=dict(t=92, b=40, l=220, r=20),   # room for the two-line title
        height=max(300, len(cats) * 38 + 120),
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

    # Determine which bin labels are present and their display order.  The two
    # label sets share their upper bins, so "did any LONG label match" is not a
    # usable test: against SHORT data it matched four of five bins and left
    # "<500" to fall through to the tail, putting the shortest reads last.
    # Pick whichever set covers more of the bins actually present.
    present_bins = set(strat_df["length_bin"].unique())
    labels = max((LENGTH_BIN_LABELS_LONG, LENGTH_BIN_LABELS_SHORT),
                 key=lambda ls: len(present_bins.intersection(ls)))
    bin_order = [b for b in labels if b in present_bins]
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

    # The legend sits below the plot, not above it.  Squeezing 13 categories
    # into the top margin worked at full width but wrapped to seven rows in a
    # 666 px report card and grew straight through the title.
    n_legend_rows = max(1, -(-len(fig.data) // 2))
    fig.update_layout(
        barmode="stack",
        title=dict(text="Broad non-canonical composition by read length",
                   x=0.5, xanchor="center", xref="container"),
        xaxis=dict(title="Fraction of bin", tickformat=".0%", range=[0, 1.0]),
        yaxis=dict(autorange="reversed", automargin=True),
        legend=dict(
            orientation="h",
            yanchor="top", y=-0.18,
            xanchor="center", x=0.5,
            font=dict(size=10),
        ),
        margin=dict(t=60, b=60 + n_legend_rows * 22, l=120, r=100),
        height=max(360, len(bin_order) * 62 + 140 + n_legend_rows * 22),
    )
    return fig


def _per_cell_violin(ct: CellTable) -> go.Figure:
    """Violin + strip plot of the legacy broad-composition column."""
    df = ct.df
    if "broad_noncanonical_read_frac" not in df.columns:
        return go.Figure()

    vals = df["broad_noncanonical_read_frac"].dropna().values

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


def _length_axis_max(*sample_groups) -> float:
    """
    Upper bound for a read-length axis, taken from the data.

    The previous fixed 5,000 bp ceiling was wrong in both directions: it clipped
    the tail (BD45 reaches 7,371 bp) while also wasting a third of the width on
    empty axis, since the 99.5th percentile of every real long-read sample here
    falls between 2,300 and 4,600 bp. Cutting at that percentile keeps a handful
    of outliers from flattening the distribution.
    """
    values = [v for group in sample_groups for series in group for v in series]
    if not values:
        return 5000.0
    values.sort()
    idx = min(len(values) - 1, int(0.995 * len(values)))
    return max(200.0, values[idx] * 1.05)


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
            xaxis=dict(title="Read length (bp)",
                       range=[0, _length_axis_max([length_samples[c] for c in signal_cats])]),
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

    _x_max = _length_axis_max(
        [length_samples[c] for c in signal_cats],
        [length_samples[c] for c in noise_cats],
    )
    fig.update_xaxes(title_text="Read length (bp)", range=[0, _x_max], row=1, col=1)
    fig.update_xaxes(title_text="Read length (bp)", range=[0, _x_max], row=1, col=2)
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
    # Each bar is already named on the axis, so hue carries no information here.
    # Four arbitrary colours only invited the reader to look for a link to the
    # read categories, which these flags are independent of.
    colors = ["#4a5568"] * len(labels)
    descriptions = [
        "Reads with TSO sequence (or its reverse complement) in a soft-clip",
        "Reads ending at an A-rich (forward) / T-rich (reverse) reference stretch",
        "Splice junctions not matching GT-AG / GC-AG / AT-AC rule",
        "Reads containing more than one TSO sequence or its reverse complement",
    ]

    denom = max(sm.n_reads_classified, 1)
    fracs = [v / denom for v in values]
    max_frac = max(fracs) if fracs else 0.01

    def _pct(f: float) -> str:
        """Enough decimals for this value, not for the largest bar.

        Keying precision to the axis maximum rendered a measured 196 reads as
        "0.00%" next to a 2.42% bar, which reads as "none found".
        """
        if f == 0:
            return "0%"
        dec = 2 if f >= 0.01 else 3 if f >= 0.001 else 4
        return f"{f:.{dec}%}"

    # Axis ticks follow the largest bar; the labels above carry their own scale.
    _dec = 2 if max_frac >= 0.01 else 3 if max_frac >= 0.001 else 4
    hover = [
        f"<b>{lbl}</b><br>{_pct(f)} of classified reads<br><i>{desc}</i>"
        for lbl, f, desc in zip(labels, fracs, descriptions)
    ]

    fig = go.Figure(go.Bar(
        x=fracs,
        y=labels,
        orientation="h",
        marker_color=colors,
        text=[_pct(f) for f in fracs],
        textposition=["inside" if f > max_frac * 0.5 else "outside" for f in fracs],
        hovertext=hover,
        hoverinfo="text",
    ))
    fig.update_layout(
        title=dict(text="Artifact flag rates (fraction of classified reads)",
                   x=0.5, xanchor="center", xref="container"),
        xaxis=dict(
            title="Fraction",
            # Flag rates run from ~1e-5 to ~1e-2 across real samples; a fixed
            # two-decimal percent collapsed the whole axis to "0.00%" repeated.
            tickformat=f".{_dec}%",
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

    # Same hues as everywhere else in the report; a second local dict here is
    # how intergenic_repeat and intergenic_enriched came to share a colour.
    cat_colours = CATEGORY_COLOURS
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

    rng = random.Random(42)  # deterministic jitter and sampling

    by_cat: dict = {}
    for locus in intergenic_loci:
        by_cat.setdefault(locus.category, []).append(locus)

    # Windows that were never promoted are, by definition, indistinguishable
    # background: tens of thousands of identical markers carry no information
    # and dominate the file size.  Sample them, keep every promoted window.
    n_sparse_total = len(by_cat.get(ReadCategory.INTERGENIC_SPARSE, []))
    n_sparse_shown = min(n_sparse_total, INTERGENIC_SCATTER_MAX_SPARSE)
    if n_sparse_total > n_sparse_shown:
        by_cat[ReadCategory.INTERGENIC_SPARSE] = rng.sample(
            by_cat[ReadCategory.INTERGENIC_SPARSE], n_sparse_shown
        )

    fig_sc = go.Figure()
    for cat, loci in by_cat.items():
        label = cat_labels.get(cat, cat.value)
        x_vals, y_vals, sizes, custom = [], [], [], []
        for locus in loci:
            jitter = rng.uniform(-jitter_scale, jitter_scale)
            log_reads = math.log10(max(1, locus.n_reads))
            x_vals.append(max(locus.poisson_pvalue_adj, 1e-300))  # avoid log(0)
            y_vals.append(10 ** (log_reads + jitter))
            sizes.append(max(6, min(24, locus.n_barcodes * 2 + 4)))
            # customdata carries the raw values; the hover string is built once
            # by Plotly from the template below.  Formatting a string per point
            # instead put tens of MB of duplicated markup in the HTML.
            custom.append((
                locus.contig, locus.start, locus.end, locus.n_reads,
                locus.n_barcodes, locus.poisson_pvalue_adj,
                locus.polya_run_downstream, locus.near_polya_site,
            ))
        fig_sc.add_trace(go.Scatter(
            x=x_vals,
            y=y_vals,
            mode="markers",
            name=label,
            marker=dict(
                color=cat_colours.get(cat, "#95a5a6"),
                size=sizes,
                opacity=0.72,
                line=dict(width=0.5, color="white"),
            ),
            customdata=custom,
            hovertemplate=(
                "%{customdata[0]}:%{customdata[1]:,}-%{customdata[2]:,}<br>"
                "Reads: %{customdata[3]:,}<br>Barcodes: %{customdata[4]:,}<br>"
                "Adj. p: %{customdata[5]:.2e}<br>"
                "polyA downstream: %{customdata[6]}<br>"
                "Near polyA site: %{customdata[7]}"
                f"<extra>{label}</extra>"
            ),
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

    _sampled_note = (
        f" · {n_sparse_shown:,} of {n_sparse_total:,} below-threshold windows sampled"
        if n_sparse_total > n_sparse_shown else ""
    )
    fig_sc.update_layout(
        title=dict(
            text="Intergenic windows — adjusted p-value vs read count"
                 "<br><sup>Point size ∝ distinct barcodes · y-axis jittered for readability"
                 f"{_sampled_note}</sup>",
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
    sort_col = "median_broad_noncanonical_read_frac"
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
        if "broad_noncanonical_read_frac" not in ct.df.columns:
            continue
        vals = ct.df["broad_noncanonical_read_frac"].dropna().values
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
        xaxis=dict(title="Read length (bp)", range=[0, _length_axis_max(
            [ls.get(ReadCategory.EXONIC_SENSE, []) for ls in (ls_a, ls_b)])]),
        yaxis=dict(title="Percent"),
        height=380,
        margin=dict(t=60, b=60, l=70, r=20),
    )
    return fig


# ---------------------------------------------------------------------------
# Cohort figures (N independent samples)
# ---------------------------------------------------------------------------

# Groups get one colour each; categories keep CATEGORY_COLOURS everywhere.
_GROUP_COLOURS = ["#393eca", "#c50037", "#33a182", "#c2ae3e",
                  "#8316a4", "#00c8d7", "#89401c", "#f33179"]


def _group_colour_map(cohort) -> dict:
    groups = [g for g in dict.fromkeys(s.group for s in cohort.samples) if g]
    return {g: _GROUP_COLOURS[i % len(_GROUP_COLOURS)] for i, g in enumerate(groups)}


def _cohort_composition_bars(cohort) -> go.Figure:
    """
    One stacked bar per sample, exonic sense excluded.

    Exonic sense runs 64-83% of reads, so including it would compress every
    difference between methods into the tail of the bar.  It is printed as text
    on the right instead, which keeps the magnitude available without spending
    the axis on it.  Fractions stay absolute (not renormalised), so bar length
    is directly the total non-canonical composition and the chart is both the
    ranking and the breakdown.
    """
    if not cohort.samples:
        return go.Figure()

    labels = [s.label for s in cohort.samples]
    shown = [c for c in CATEGORY_ORDER
             if c is not ReadCategory.EXONIC_SENSE
             and any((s.get(f"read_frac_{c.value}") or 0) > 0 for s in cohort.samples)]

    fig = go.Figure()
    for cat in shown:
        values = [s.get(f"read_frac_{cat.value}") or 0.0 for s in cohort.samples]
        fig.add_trace(go.Bar(
            x=values, y=labels, name=CATEGORY_LABELS[cat], orientation="h",
            marker_color=CATEGORY_COLOURS[cat],
            marker_line=dict(color="white", width=1),   # 2px surface gap between segments
            customdata=[[CATEGORY_LABELS[cat]]] * len(values),
            hovertemplate="%{y}<br>%{customdata[0]}: %{x:.2%}<extra></extra>",
        ))

    totals = [sum(s.get(f"read_frac_{c.value}") or 0.0 for c in shown)
              for s in cohort.samples]
    for label, total, s in zip(labels, totals, cohort.samples):
        es = s.get("read_frac_exonic_sense")
        fig.add_annotation(
            x=total, y=label, xref="x", yref="y",
            text=f"  exonic sense {es:.1%}" if es is not None else "  exonic sense n/a",
            showarrow=False, xanchor="left",
            font=dict(size=10, color="#636e72"),
        )

    x_max = max(totals) if totals else 1.0
    fig.update_layout(
        barmode="stack",
        title=dict(text="Alignment composition excluding exonic sense", x=0.5),
        xaxis=dict(title="Fraction of classified reads", tickformat=".0%",
                   range=[0, x_max * 1.32]),
        yaxis=dict(autorange="reversed", automargin=True),
        legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.18,
                    font=dict(size=10), traceorder="normal"),
        height=max(360, len(labels) * 34 + 220),
        margin=dict(t=60, b=120, l=190, r=40),
    )
    return fig


def _cohort_deviation_heatmap(cohort) -> go.Figure:
    """
    Samples x categories, coloured by deviation from the cohort median.

    Colouring by the raw fraction would only restate the composition figure.
    Colouring by deviation answers the question the composition figure cannot:
    which category is unusual for this method.  Cell text keeps the raw value
    so the absolute number is never lost.

    Below four samples a cohort median is not meaningful, so the colour falls
    back to the raw fraction and the title says so.
    """
    metrics = cohort.category_metrics("read_frac_")
    if not cohort.samples or not metrics:
        return go.Figure()

    labels = [s.label for s in cohort.samples]
    raw = [[s.get(m) for s in cohort.samples] for m in metrics]

    use_deviation = len(cohort.samples) >= 4
    # A row where every sample sits near zero has a spread of ~0, and scaling by
    # it would paint a 0.0003 percentage-point difference as the strongest
    # signal in the figure.  Such rows stay neutral; the text still shows the
    # values.  The same applies when too few samples reported the category for a
    # median to mean anything.
    MIN_ROW_MAX = 0.001      # 0.1% of reads
    MIN_ROW_SAMPLES = 3
    z, text = [], []
    for row in raw:
        present = [v for v in row if v is not None]
        median = sorted(present)[len(present) // 2] if present else 0.0
        spread = (max(present) - min(present)) if present else 0.0
        informative = (
            use_deviation and spread > 0
            and len(present) >= MIN_ROW_SAMPLES
            and max(present) >= MIN_ROW_MAX
        )
        z.append([
            None if v is None else ((v - median) / spread if informative else 0.0)
            for v in row
        ])
        text.append(["n/a" if v is None else f"{v:.2%}" for v in row])

    unstranded = {s.label for s in cohort.samples if s.is_unstranded}
    y_labels = [
        CATEGORY_LABELS.get(next(c for c in CATEGORY_ORDER
                                 if m == f"read_frac_{c.value}"), m)
        for m in metrics
    ]

    fig = go.Figure(go.Heatmap(
        z=z, x=labels, y=y_labels,
        colorscale="RdBu_r", zmid=0, zmin=-1, zmax=1,
        text=text, texttemplate="%{text}", textfont=dict(size=9),
        hovertemplate="%{x}<br>%{y}<br>%{text}<extra></extra>",
        colorbar=dict(title="vs cohort<br>median",
                      tickvals=[-1, 0, 1], ticktext=["lowest", "median", "highest"]),
    ))

    # Mark the cells whose meaning changes with library strandedness.
    if cohort.has_mixed_strandedness:
        for xi, label in enumerate(labels):
            if label not in unstranded:
                continue
            for yi, m in enumerate(metrics):
                if m in {"read_frac_exonic_antisense"}:
                    fig.add_shape(type="rect", x0=xi - 0.5, x1=xi + 0.5,
                                  y0=yi - 0.5, y1=yi + 0.5,
                                  line=dict(color="#2d3436", width=2, dash="dot"),
                                  fillcolor="rgba(0,0,0,0)")

    subtitle = ("Colour is deviation from the cohort median; cell text is the raw fraction"
                if use_deviation else
                "Fewer than 4 samples: no cohort median, so colour is not scaled")
    subtitle += "; rows that are near-zero everywhere are left neutral"
    if cohort.has_mixed_strandedness:
        subtitle += " · dotted cells are unstranded and not comparable here"

    fig.update_layout(
        title=dict(text=f"Category signature by sample<br><sup>{subtitle}</sup>", x=0.5),
        xaxis=dict(tickangle=-30, automargin=True),
        yaxis=dict(autorange="reversed", automargin=True),
        height=max(340, len(metrics) * 30 + 190),
        margin=dict(t=100, b=90, l=210, r=70),
    )
    return fig


def _cohort_artifact_dots(cohort) -> go.Figure:
    """
    One row per artifact flag, one dot per sample, log x-axis.

    These rates span roughly four orders of magnitude across real samples, so a
    linear axis collapses all but the largest onto zero.  Dots rather than bars
    because bar length measured from an arbitrary log baseline means nothing.
    The vertical tick is the cohort median for that flag.

    Every flag keeps its row even when no sample reports it, so a missing
    measurement is visible as an empty row rather than as an absent concept.
    """
    from scnoisemeter.modules.cohort import ARTIFACT_FLAG_COUNTS

    if not cohort.samples:
        return go.Figure()

    group_colours = _group_colour_map(cohort)
    has_groups = any(s.group for s in cohort.samples)
    n_samples = len(cohort.samples)

    row_labels, rates_by_row = [], []
    for key, flag_label in ARTIFACT_FLAG_COUNTS.items():
        rates = [(s, s.flag_rate(key)) for s in cohort.samples]
        n_missing = sum(1 for _, r in rates if r is None)
        n_zero = sum(1 for _, r in rates if r == 0)
        note = ""
        if n_missing:
            note = f"<br><sup>not reported by {n_missing}/{n_samples}</sup>"
        elif n_zero == n_samples:
            note = f"<br><sup>zero in all {n_samples}</sup>"
        row_labels.append(flag_label + note)
        rates_by_row.append([(s, r) for s, r in rates if r is not None and r > 0])

    fig = go.Figure()
    seen_groups: set = set()
    for label, present in zip(row_labels, rates_by_row):
        if not present:
            # Register the category so a flag nothing reported still gets a row.
            # Dropping it would read as "this artifact does not exist".
            fig.add_trace(go.Scatter(x=[None], y=[label], mode="markers",
                                     showlegend=False, hoverinfo="skip"))
        if present:
            values = sorted(r for _, r in present)
            median = values[len(values) // 2]
            fig.add_trace(go.Scatter(
                x=[median], y=[label], mode="markers",
                marker=dict(symbol="line-ns-open", size=26, color="#2d3436",
                            line=dict(width=2)),
                showlegend=False,
                hovertemplate=f"cohort median {median:.2e}<extra></extra>",
            ))
        for s, rate in present:
            group = s.group or "sample"
            fig.add_trace(go.Scatter(
                x=[rate], y=[label], mode="markers",
                name=group, legendgroup=group,
                showlegend=has_groups and group not in seen_groups,
                marker=dict(size=12, opacity=0.85,
                            color=group_colours.get(s.group, "#393eca"),
                            line=dict(width=1.5, color="white")),
                customdata=[[s.label]],
                hovertemplate="%{customdata[0]}<br>%{x:.3e}<extra></extra>",
            ))
            seen_groups.add(group)

    if not any(rates_by_row):
        return go.Figure()

    fig.update_layout(
        title=dict(
            text="Artifact flag rates by sample"
                 "<br><sup>Log axis: these rates span orders of magnitude · "
                 "vertical tick is the cohort median</sup>",
            x=0.5,
        ),
        xaxis=dict(title="Fraction of classified reads", type="log",
                   showgrid=True, dtick=1, tickformat=".0e",
                   minor=dict(showgrid=False)),
        yaxis=dict(automargin=True, categoryorder="array",
                   categoryarray=row_labels[::-1]),   # first flag at the top
        legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.3),
        height=max(300, len(row_labels) * 52 + 190),
        margin=dict(t=95, b=110 if has_groups else 80, l=210, r=50),
    )
    return fig


def _cohort_percell_boxes(cohort) -> go.Figure:
    """
    Per-cell composition spread, for the samples that have real per-cell data.

    A tight box means the composition is uniform across cells and therefore
    protocol-level; a long tail means a subpopulation drives it.  Samples run in
    barcode-agnostic mode have no per-cell values at all and are named in the
    subtitle rather than silently dropped.
    """
    col = "broad_noncanonical_read_frac"
    with_cells = [s for s in cohort.samples
                  if s.has_cells and col in (s.cell_df.columns if s.cell_df is not None else [])]
    if not with_cells:
        return go.Figure()

    import numpy as np

    group_colours = _group_colour_map(cohort)
    fig = go.Figure()
    for s in reversed(with_cells):        # keep composition-figure order top-down
        values = s.cell_df[col].dropna().to_numpy(dtype=float)
        if not len(values):
            continue
        # Summarise here rather than shipping every cell to the browser: one
        # sample in a real cohort has 176k cells, and the raw arrays alone put
        # megabytes of numbers into the HTML for a five-number summary.
        q1, med, q3 = (float(v) for v in np.percentile(values, [25, 50, 75]))
        iqr = q3 - q1
        lo = float(values[values >= q1 - 1.5 * iqr].min()) if len(values) else q1
        hi = float(values[values <= q3 + 1.5 * iqr].max()) if len(values) else q3
        fig.add_trace(go.Box(
            name=s.label, orientation="h", y=[s.label],
            q1=[q1], median=[med], q3=[q3], lowerfence=[lo], upperfence=[hi],
            mean=[float(values.mean())],
            marker_color=group_colours.get(s.group, "#393eca"),
            line=dict(width=1.5),
            hovertext=[f"{s.label}<br>{len(values):,} cells<br>"
                       f"median {med:.2%} · IQR {q1:.2%}–{q3:.2%}"],
            hoverinfo="text",
        ))

    excluded = [s.label for s in cohort.samples if s not in with_cells]
    subtitle = f"{len(with_cells)} of {len(cohort.samples)} samples have per-cell data"
    if excluded:
        subtitle += " · no cells in: " + ", ".join(excluded)

    fig.update_layout(
        title=dict(
            text=f"Per-cell non-canonical composition<br><sup>{subtitle}</sup>", x=0.5),
        xaxis=dict(title="Fraction of the cell's reads", tickformat=".0%"),
        yaxis=dict(automargin=True),
        showlegend=False,
        height=max(320, len(with_cells) * 42 + 170),
        margin=dict(t=95, b=60, l=190, r=40),
    )
    return fig
