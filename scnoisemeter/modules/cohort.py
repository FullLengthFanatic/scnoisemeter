"""
cohort.py
=========
Reads finished `run` / `run-plate` output directories and assembles them into a
cross-sample comparison.

Why this is separate from `compare`
-----------------------------------
`compare` answers "what did this filtering step do to these reads": it matches
reads by name between two BAMs and reports retention and category transitions.
That is only defined for a nested pair.

`cohort` answers a different question — "where are the artifacts in each of
these methods" — across samples that share no reads and usually no barcodes.
It therefore reads the metrics that `run` already wrote rather than the BAMs,
which also means a comparison across a dozen samples costs seconds instead of
the many hours a re-classification would take.

Three rules the real data forces
--------------------------------
1. A metric that is missing is reported as absent, never as zero.  Metrics have
   been added over time (``n_tso_concatemer`` does not exist before v0.7), and
   rendering an absent metric as 0 reads as a measurement of "none".
2. Results produced by different scnoisemeter versions are flagged loudly.
   Methodology changes move metrics by orders of magnitude, so pooling versions
   silently manufactures differences that are not biological.
3. Stranded and unstranded protocols are not comparable in every category.
   ``exonic_antisense`` and ``strand_concordance`` mean different things in a
   non-stranded library, so those cells are marked rather than quietly ranked.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from scnoisemeter.constants import CATEGORY_ORDER

logger = logging.getLogger(__name__)


# Pre-v0.7 result directories carry the old aggregate names.  Reading them is
# one-way: nothing writes these any more.
_LEGACY_METRIC_ALIASES = {
    "noise_read_frac":        "broad_noncanonical_read_frac",
    "noise_base_frac":        "broad_noncanonical_base_frac",
    "noise_read_frac_strict": "artifact_candidate_read_frac",
    "noise_base_frac_strict": "artifact_candidate_base_frac",
}

# Artifact flags are stored as counts; the cohort compares them as rates.
ARTIFACT_FLAG_COUNTS = {
    "n_tso_invasion":      "TSO invasion",
    "n_polya_priming":     "Internal polyA priming",
    "n_noncanon_junction": "Non-canonical junction",
    "n_tso_concatemer":    "TSO concatemer",
}

# Categories whose meaning depends on library strandedness.
STRAND_DEPENDENT_METRICS = {"read_frac_exonic_antisense",
                            "base_frac_exonic_antisense",
                            "strand_concordance"}


@dataclass
class CohortSample:
    """One finished result directory, parsed."""
    sample_id:   str
    label:       str
    results_dir: Path
    group:       Optional[str] = None
    order:       Optional[int] = None
    metrics:     dict = field(default_factory=dict)   # name -> float, absent keys omitted
    run_info:    Optional[dict] = None
    cell_df:     Optional[pd.DataFrame] = None

    @property
    def platform(self) -> Optional[str]:
        return (self.run_info or {}).get("platform")

    @property
    def version(self) -> Optional[str]:
        return (self.run_info or {}).get("scnoisemeter_version")

    @property
    def is_unstranded(self) -> bool:
        return bool((self.run_info or {}).get("unstranded"))

    @property
    def has_cells(self) -> bool:
        """True when per-cell values are real rather than the barcode-agnostic sentinel."""
        return self.cell_df is not None and len(self.cell_df) > 1

    def get(self, metric: str) -> Optional[float]:
        """Value of `metric`, or None when this sample never reported it."""
        return self.metrics.get(metric)

    def flag_rate(self, count_key: str) -> Optional[float]:
        """Artifact-flag count as a fraction of classified reads, or None if absent."""
        n = self.metrics.get(count_key)
        denom = self.metrics.get("n_reads_classified")
        if n is None or not denom:
            return None
        return n / denom


@dataclass
class Cohort:
    samples: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.samples)

    @property
    def labels(self) -> list:
        return [s.label for s in self.samples]

    @property
    def has_mixed_strandedness(self) -> bool:
        return len({s.is_unstranded for s in self.samples}) > 1

    def category_metrics(self, prefix: str = "read_frac_") -> list:
        """Category metric names present in at least one sample with a non-zero value.

        Four of the seventeen categories are zero in every real sample seen so
        far; carrying them into a heatmap wastes rows on nothing.
        """
        keep = []
        for cat in CATEGORY_ORDER:
            name = f"{prefix}{cat.value}"
            if any((s.get(name) or 0) > 0 for s in self.samples):
                keep.append(name)
        return keep


def _read_metrics_tsv(path: Path) -> dict:
    """Parse a two-column `metric<TAB>value` file into floats, skipping blanks."""
    metrics: dict = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] == "metric":
                continue
            key, raw = parts[0], parts[1].strip()
            if raw == "":
                continue          # written but not measured — stays absent
            try:
                metrics[key] = float(raw)
            except ValueError:
                continue          # non-numeric provenance rows, if any
    for old, new in _LEGACY_METRIC_ALIASES.items():
        if new not in metrics and old in metrics:
            metrics[new] = metrics[old]
    return metrics


def _single(results_dir: Path, pattern: str) -> Optional[Path]:
    hits = sorted(results_dir.glob(pattern))
    return hits[0] if hits else None


def load_cohort(results_dirs, sample_sheet: "Optional[pd.DataFrame]" = None) -> Cohort:
    """
    Load result directories into a :class:`Cohort`.

    Parameters
    ----------
    results_dirs : iterable of path-like
        Output directories previously written by `run` or `run-plate`.
    sample_sheet : DataFrame, optional
        Indexed by sample id, with optional `label`, `group` and `order` columns.
    """
    cohort = Cohort()
    seen: dict = {}

    for raw_dir in results_dirs:
        results_dir = Path(raw_dir)
        metrics_path = _single(results_dir, "*.read_metrics.tsv")
        if metrics_path is None:
            cohort.warnings.append(
                f"{results_dir}: no *.read_metrics.tsv found — directory skipped."
            )
            continue

        sample_id = metrics_path.name[: -len(".read_metrics.tsv")]
        if sample_id in seen:
            sample_id = f"{sample_id} ({results_dir.name})"

        run_info = None
        info_path = _single(results_dir, "*.run_info.json")
        if info_path is not None:
            try:
                run_info = json.loads(info_path.read_text())
            except (OSError, ValueError) as exc:
                cohort.warnings.append(f"{results_dir}: unreadable run_info.json ({exc}).")

        cell_df = None
        cells_path = _single(results_dir, "*.cell_metrics.tsv")
        if cells_path is not None:
            try:
                cell_df = pd.read_csv(cells_path, sep="\t", index_col=0)
                for old, new in _LEGACY_METRIC_ALIASES.items():
                    if new not in cell_df.columns and old in cell_df.columns:
                        cell_df[new] = cell_df[old]
            except (OSError, ValueError) as exc:
                cohort.warnings.append(f"{results_dir}: unreadable cell_metrics.tsv ({exc}).")

        sample = CohortSample(
            sample_id=sample_id,
            label=sample_id,
            results_dir=results_dir,
            metrics=_read_metrics_tsv(metrics_path),
            run_info=run_info,
            cell_df=cell_df,
        )

        if sample_sheet is not None and sample_id in sample_sheet.index:
            row = sample_sheet.loc[sample_id]
            sample.label = str(row.get("label") or sample_id)
            sample.group = str(row["group"]) if row.get("group") is not None else None
            try:
                sample.order = int(row["order"])
            except (KeyError, TypeError, ValueError):
                sample.order = None

        seen[sample_id] = sample
        cohort.samples.append(sample)

    if not cohort.samples:
        return cohort

    _order_samples(cohort)
    _check_comparability(cohort)
    return cohort


def ranking_metric(cohort: Cohort) -> str:
    """
    Metric to rank samples by.

    Prefer the artifact-candidate fraction, but fall back to the broad
    composition when some samples predate it: ranking a cohort by a metric most
    of its members never reported puts them all in an arbitrary tail.
    """
    preferred = "artifact_candidate_read_frac"
    if all(s.get(preferred) is not None for s in cohort.samples):
        return preferred
    return "broad_noncanonical_read_frac"


def _order_samples(cohort: Cohort) -> None:
    """Sample-sheet order when given, otherwise cleanest sample first."""
    if any(s.order is not None for s in cohort.samples):
        cohort.samples.sort(key=lambda s: (s.order is None, s.order or 0, s.label))
        return
    metric = ranking_metric(cohort)
    cohort.samples.sort(
        key=lambda s: (s.get(metric) is None, s.get(metric) or 0.0)
    )


def _check_comparability(cohort: Cohort) -> None:
    """Warn about the ways a cross-sample comparison can be misleading."""
    missing_info = [s.label for s in cohort.samples if s.run_info is None]
    if missing_info:
        cohort.warnings.append(
            "No run_info.json in: " + ", ".join(missing_info) +
            ". These directories predate provenance recording, so tool version, "
            "platform and annotation cannot be checked for comparability."
        )

    versions = {s.version for s in cohort.samples if s.version}
    if len(versions) > 1:
        cohort.warnings.append(
            "Samples were produced by different scnoisemeter versions ("
            + ", ".join(sorted(versions)) +
            "). Metric definitions changed between releases, so differences "
            "between samples may reflect the tool rather than the library. "
            "Re-run the older samples before drawing conclusions."
        )

    gtf_versions = {
        (s.run_info or {}).get("gtf_version")
        for s in cohort.samples if (s.run_info or {}).get("gtf_version")
    }
    if len(gtf_versions) > 1:
        cohort.warnings.append(
            "Samples were annotated with different GENCODE releases ("
            + ", ".join(str(v) for v in sorted(gtf_versions, key=str)) +
            "). Categories that depend on annotation completeness, in particular "
            "the intergenic ones, are not directly comparable across them."
        )

    if cohort.has_mixed_strandedness:
        unstranded = [s.label for s in cohort.samples if s.is_unstranded]
        cohort.warnings.append(
            "The cohort mixes stranded and unstranded protocols (unstranded: "
            + ", ".join(unstranded) +
            "). Exonic antisense and strand concordance are not comparable "
            "across that boundary and are marked in the figures."
        )

    no_cells = [s.label for s in cohort.samples if not s.has_cells]
    if no_cells:
        cohort.warnings.append(
            "No per-cell data for: " + ", ".join(no_cells) +
            ". These were classified in barcode-agnostic mode, so they are "
            "absent from the per-cell figure but present everywhere else."
        )


def build_summary_table(cohort: Cohort) -> pd.DataFrame:
    """One row per sample: identity, provenance and the headline numbers."""
    rows = []
    for s in cohort.samples:
        row = {
            "sample":      s.sample_id,
            "label":       s.label,
            "group":       s.group,
            "platform":    s.platform,
            "version":     s.version,
            "unstranded":  s.is_unstranded,
            "n_reads_classified": s.get("n_reads_classified"),
            "n_cells":     s.get("n_cells"),
            "broad_noncanonical_read_frac": s.get("broad_noncanonical_read_frac"),
            "artifact_candidate_read_frac": s.get("artifact_candidate_read_frac"),
            "strand_concordance": s.get("strand_concordance"),
            "chimeric_read_frac": s.get("chimeric_read_frac"),
        }
        for key, label in ARTIFACT_FLAG_COUNTS.items():
            row[f"{key}_rate"] = s.flag_rate(key)
        rows.append(row)
    return pd.DataFrame(rows)


def build_composition_table(cohort: Cohort) -> pd.DataFrame:
    """Samples x categories, read and base fractions, absent values left as NaN."""
    rows = []
    for s in cohort.samples:
        row = {"sample": s.sample_id, "label": s.label, "group": s.group}
        for cat in CATEGORY_ORDER:
            row[f"read_frac_{cat.value}"] = s.get(f"read_frac_{cat.value}")
            row[f"base_frac_{cat.value}"] = s.get(f"base_frac_{cat.value}")
        rows.append(row)
    return pd.DataFrame(rows)
