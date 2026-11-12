"""
selfcheck.py
============
Runtime invariant checks over a finished result set.

The test suite proves the code behaves on constructed inputs.  These checks ask
a different question, on real data, every run: does *this* result set satisfy
the properties the tool claims about its own output?

The distinction matters because the strongest claims scNoiseMeter makes are
structural rather than statistical -- categories are mutually exclusive, bases
are partitioned without double-counting, one composition is a subset of another
-- and each of those is a subtraction away from being verifiable.  A silent
violation would not look like a crash; it would look like a slightly wrong
number in a table someone later cites.

Severity
--------
``ERROR``   A property that cannot legitimately fail.  Its violation means a
            reported number is wrong, so the run stops unless ``--no-strict``.
``WARNING`` A property that a documented, legitimate configuration can break.
            Reservoir-estimated intergenic reclassification is the main one:
            promoted totals are scaled from a uniform sub-sample, so exact
            conservation is not expected to survive rounding.

Every check records the observed values, not just a verdict, so ``qc.json``
answers "by how much" without a second run.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from scnoisemeter.constants import (
    ARTIFACT_CANDIDATE_CATEGORIES,
    BROAD_NONCANONICAL_CATEGORIES,
    CATEGORY_ORDER,
    ReadCategory,
)

logger = logging.getLogger(__name__)

ERROR = "error"
WARNING = "warning"

# Composition fractions are rounded to 4 decimals on the way out, and the
# intergenic reclassification scales integer counts, so exact equality is the
# wrong test even when nothing is broken.
FRACTION_SUM_TOLERANCE = 1e-3
COUNT_RELATIVE_TOLERANCE = 1e-6


@dataclass
class Invariant:
    """One checked property of a result set."""
    name: str
    ok: bool
    severity: str
    statement: str
    observed: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "ok": self.ok,
            "severity": self.severity,
            "statement": self.statement,
            "observed": self.observed,
        }


def _num(value):
    """*value* as a number, or None when it is not one.

    NaN is how an unmeasured metric is carried (a category with no reads has
    no UMI diversity, a sample with no exonic reads has no strand concordance),
    so it is absent, not out of range.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    if value != value:
        return None
    return value


def run_invariants(sm, ct=None, result=None) -> list:
    """Check every structural property of a finished result set.

    Returns the full list, passing and failing alike: a check that ran and
    passed is evidence, and a report that lists only failures cannot be
    distinguished from a report whose checks never ran.

    A check whose inputs are absent is skipped rather than failed -- an
    optional metric that was never computed is not a violated invariant.
    """
    checks: list = []

    def add(name, ok, severity, statement, **observed):
        checks.append(Invariant(name=name, ok=bool(ok), severity=severity,
                                statement=statement, observed=observed))

    sampled = bool(getattr(result, "_intergenic_reclassification_sampled", False))
    # Reservoir-estimated promotion scales counts, so conservation is expected
    # to hold only up to rounding; saying so is more useful than a hard failure.
    conservation_severity = WARNING if sampled else ERROR

    # -- compositions sum to one ---------------------------------------------
    for label, fracs in (("read", getattr(sm, "read_fracs", None)),
                         ("base", getattr(sm, "base_fracs", None))):
        if not fracs:
            continue
        total = sum(v for v in fracs.values() if _num(v) is not None)
        add(f"{label}_fractions_sum_to_one",
            abs(total - 1.0) <= FRACTION_SUM_TOLERANCE,
            conservation_severity,
            f"The {label} composition over all categories sums to 1.",
            observed_sum=round(total, 6),
            tolerance=FRACTION_SUM_TOLERANCE,
            reservoir_estimated=sampled)

    # -- every fraction is a fraction ----------------------------------------
    out_of_range = {}
    for attr in (
        "broad_noncanonical_read_frac", "broad_noncanonical_base_frac",
        "artifact_candidate_read_frac", "artifact_candidate_base_frac",
        "strand_concordance", "chimeric_read_frac", "multimapper_read_frac",
        "unmapped_read_frac", "low_mapq_read_frac", "softclip_base_frac",
        "three_prime_anchored_frac", "three_prime_at_pas_site_frac",
        "five_prime_anchored_frac", "both_ends_anchored_frac",
        "full_length_read_frac", "per_cell_noise_median",
    ):
        value = _num(getattr(sm, attr, None))
        if value is not None and not (0.0 <= value <= 1.0):
            out_of_range[attr] = value
    for label, fracs in (("read_fracs", getattr(sm, "read_fracs", None) or {}),
                         ("base_fracs", getattr(sm, "base_fracs", None) or {})):
        for cat, value in fracs.items():
            v = _num(value)
            if v is not None and not (0.0 <= v <= 1.0):
                out_of_range[f"{label}[{cat}]"] = v
    add("fractions_within_unit_interval", not out_of_range, ERROR,
        "Every reported fraction lies in [0, 1].",
        out_of_range=out_of_range)

    # -- the narrow composition is a subset of the broad one -----------------
    # This is a claim about the taxonomy, not about any sample: the
    # artifact-candidate categories are a strict subset of the broad
    # non-canonical ones, so no sample can invert the two fractions.
    add("artifact_candidate_categories_subset_of_broad",
        ARTIFACT_CANDIDATE_CATEGORIES <= BROAD_NONCANONICAL_CATEGORIES,
        ERROR,
        "The artifact-candidate categories are a subset of the broad "
        "non-canonical categories.",
        artifact_candidate=sorted(c.value for c in ARTIFACT_CANDIDATE_CATEGORIES),
        broad_only=sorted(
            c.value for c in BROAD_NONCANONICAL_CATEGORIES - ARTIFACT_CANDIDATE_CATEGORIES
        ))
    for label in ("read", "base"):
        narrow = _num(getattr(sm, f"artifact_candidate_{label}_frac", None))
        broad = _num(getattr(sm, f"broad_noncanonical_{label}_frac", None))
        if narrow is None or broad is None:
            continue
        add(f"artifact_candidate_{label}_frac_not_above_broad",
            narrow <= broad + FRACTION_SUM_TOLERANCE, ERROR,
            f"artifact_candidate_{label}_frac does not exceed "
            f"broad_noncanonical_{label}_frac, because its categories are a "
            f"subset of them.",
            narrow=narrow, broad=broad)

    # -- endpoint anchoring is a conjunction ---------------------------------
    both = _num(getattr(sm, "both_ends_anchored_frac", None))
    three = _num(getattr(sm, "three_prime_anchored_frac", None))
    five = _num(getattr(sm, "five_prime_anchored_frac", None))
    if both is not None and three is not None and five is not None:
        add("both_ends_not_above_either_end",
            both <= min(three, five) + FRACTION_SUM_TOLERANCE, ERROR,
            "both_ends_anchored_frac does not exceed either single-end "
            "fraction, because it is measured on the same reads as their "
            "conjunction.",
            both_ends=both, three_prime=three, five_prime=five)
    pas = _num(getattr(sm, "three_prime_at_pas_site_frac", None))
    if pas is not None and three is not None:
        add("pas_anchoring_not_above_all_site_anchoring",
            pas <= three + FRACTION_SUM_TOLERANCE, WARNING,
            "three_prime_at_pas_site_frac does not exceed "
            "three_prime_anchored_frac, because hexamer-annotated sites are a "
            "subset of all sites. Expected to hold unless "
            "--anchor-veto-internal-priming is in use, which vetoes the "
            "unrestricted fraction but not the PAS subset.",
            pas=pas, all_sites=three)

    # -- UMI diversity is a ratio of unique strings to reads -----------------
    bad_umi = {
        cat: v for cat, v in (getattr(sm, "umi_sequence_diversity", None) or {}).items()
        if _num(v) is not None and not (0.0 <= v <= 1.0)
    }
    add("umi_diversity_within_unit_interval", not bad_umi, ERROR,
        "Every UMI sequence-diversity ratio lies in [0, 1]: unique UMI "
        "strings cannot outnumber the reads carrying them.",
        out_of_range=bad_umi)

    # -- read accounting -----------------------------------------------------
    classified = _num(getattr(sm, "n_reads_classified", None))
    if classified is not None and getattr(sm, "read_fracs", None):
        # Flag counts are counts of classified reads, so none can exceed the
        # denominator every one of them is reported over.
        over = {}
        for attr in ("n_tso_invasion", "n_polya_priming", "n_polya_tail",
                     "n_rt_switch_repeat", "n_noncanon_junction",
                     "n_tso_concatemer", "n_discordant_pair", "n_low_mapq"):
            v = _num(getattr(sm, attr, None))
            if v is not None and v > classified:
                over[attr] = v
        add("flag_counts_not_above_classified_reads", not over, ERROR,
            "No artifact-flag count exceeds n_reads_classified, the "
            "denominator each of them is reported over.",
            n_reads_classified=classified, exceeding=over)

    unassigned = _num(getattr(sm, "n_reads_unassigned", None))
    if classified is not None and unassigned is not None:
        add("unassigned_not_above_classified", unassigned <= classified, ERROR,
            "UNASSIGNED reads are part of the classified denominator, so they "
            "cannot exceed it.",
            n_reads_unassigned=unassigned, n_reads_classified=classified)

    # -- the cell table agrees with the scalar cell count --------------------
    df = getattr(ct, "df", None)
    n_cells = _num(getattr(sm, "n_cells", None))
    if df is not None and n_cells is not None and hasattr(df, "__len__"):
        add("n_cells_matches_cell_table", len(df) == n_cells, ERROR,
            "n_cells equals the number of rows in the per-cell table.",
            n_cells=n_cells, cell_table_rows=len(df))

    # -- strand concordance is the ratio it is documented to be --------------
    concordance = _num(getattr(sm, "strand_concordance", None))
    fracs = getattr(sm, "read_fracs", None) or {}
    sense = _num(fracs.get(ReadCategory.EXONIC_SENSE.value))
    anti = _num(fracs.get(ReadCategory.EXONIC_ANTISENSE.value))
    if concordance is not None and sense is not None and anti is not None:
        denom = sense + anti
        expected = (sense / denom) if denom > 0 else 0.0
        add("strand_concordance_matches_its_definition",
            abs(concordance - expected) <= 1e-2, WARNING,
            "strand_concordance equals exonic_sense / (exonic_sense + "
            "exonic_antisense), recomputed from the reported composition. "
            "Rounded composition fractions make this approximate.",
            reported=concordance, recomputed=round(expected, 6))

    # -- the category list itself --------------------------------------------
    add("category_order_has_no_duplicates",
        len(CATEGORY_ORDER) == len(set(CATEGORY_ORDER)), ERROR,
        "CATEGORY_ORDER lists each category once, so no read or base is "
        "counted twice when a composition is assembled.",
        n_categories=len(CATEGORY_ORDER))

    return checks


def failures(checks: list, severity: str = ERROR) -> list:
    return [c for c in checks if not c.ok and c.severity == severity]


def write_qc_report(checks: list, output_dir, sample_name: str) -> Optional[Path]:
    """Write ``<sample>.qc.json``: every check, its verdict and its numbers."""
    if not checks:
        return None
    path = Path(output_dir) / f"{sample_name}.qc.json"
    payload = {
        "sample_name": sample_name,
        "n_checks": len(checks),
        "n_failed_error": len(failures(checks, ERROR)),
        "n_failed_warning": len(failures(checks, WARNING)),
        "checks": [c.as_dict() for c in checks],
    }
    try:
        path.write_text(json.dumps(payload, indent=2, default=str))
    except OSError as exc:
        logger.warning("Could not write %s: %s", path.name, exc)
        return None
    return path


def summarise(checks: list) -> str:
    """One line naming what failed, or that nothing did."""
    errors = failures(checks, ERROR)
    warns = failures(checks, WARNING)
    if not errors and not warns:
        return f"All {len(checks)} result invariants hold."
    parts = []
    if errors:
        parts.append(f"{len(errors)} violated: " + ", ".join(c.name for c in errors))
    if warns:
        parts.append(f"{len(warns)} to check: " + ", ".join(c.name for c in warns))
    return "; ".join(parts)
