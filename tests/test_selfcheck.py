"""
tests/test_selfcheck.py
=======================
Screening item 7: runtime invariants over a finished result set.

These tests are about the checker, not the tool: for each invariant, a healthy
result set must pass and a specifically corrupted one must fail. An invariant
that cannot fail is not a check, and one that fires on valid data is worse than
none, so both directions are asserted.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from scnoisemeter.constants import CATEGORY_ORDER, ReadCategory  # noqa: E402
from scnoisemeter.modules.metrics import SampleMetrics  # noqa: E402
from scnoisemeter.modules.selfcheck import (  # noqa: E402
    ERROR,
    WARNING,
    failures,
    run_invariants,
    summarise,
    write_qc_report,
)


class _CT:
    def __init__(self, n_rows):
        self.df = pd.DataFrame(index=range(n_rows))


class _Result:
    def __init__(self, sampled=False):
        self._intergenic_reclassification_sampled = sampled


@pytest.fixture
def healthy():
    """A result set that satisfies every invariant."""
    sm = SampleMetrics(sample_name="s", bam_path="b", platform="p",
                       pipeline_stage="x", aligner="a")
    sm.read_fracs = {c.value: 0.0 for c in CATEGORY_ORDER}
    sm.base_fracs = {c.value: 0.0 for c in CATEGORY_ORDER}
    sm.read_fracs[ReadCategory.EXONIC_SENSE.value] = 0.8
    sm.read_fracs[ReadCategory.EXONIC_ANTISENSE.value] = 0.2
    sm.base_fracs[ReadCategory.EXONIC_SENSE.value] = 1.0
    sm.strand_concordance = 0.8
    sm.n_reads_classified = 1000
    sm.n_reads_unassigned = 10
    sm.n_cells = 3
    sm.broad_noncanonical_read_frac = 0.20
    sm.artifact_candidate_read_frac = 0.05
    sm.three_prime_anchored_frac = 0.90
    sm.five_prime_anchored_frac = 0.80
    sm.both_ends_anchored_frac = 0.75
    sm.umi_sequence_diversity = {"exonic_sense": 0.4}
    return sm, _CT(3)


def _names(checks, severity=ERROR):
    return {c.name for c in failures(checks, severity)}


class TestHealthyResultSet:
    def test_nothing_fires_on_a_valid_result(self, healthy):
        sm, ct = healthy
        checks = run_invariants(sm, ct, _Result())
        assert checks, "no invariants ran"
        assert failures(checks, ERROR) == []
        assert failures(checks, WARNING) == []

    def test_passing_checks_are_reported_not_omitted(self, healthy):
        """A report listing only failures cannot be told from one that never ran."""
        sm, ct = healthy
        checks = run_invariants(sm, ct, _Result())
        assert all(c.ok for c in checks)
        assert "invariants hold" in summarise(checks)


class TestEachInvariantCanFail:
    """One corruption per invariant, each asserted to fire by name."""

    def test_narrow_composition_above_broad(self, healthy):
        sm, ct = healthy
        sm.artifact_candidate_read_frac = 0.9
        assert "artifact_candidate_read_frac_not_above_broad" in _names(
            run_invariants(sm, ct, _Result()))

    def test_fraction_outside_unit_interval(self, healthy):
        sm, ct = healthy
        sm.strand_concordance = 1.4
        assert "fractions_within_unit_interval" in _names(
            run_invariants(sm, ct, _Result()))

    def test_per_category_fraction_outside_unit_interval(self, healthy):
        sm, ct = healthy
        sm.read_fracs[ReadCategory.INTRONIC_PURE.value] = -0.1
        assert "fractions_within_unit_interval" in _names(
            run_invariants(sm, ct, _Result()))

    def test_both_ends_above_a_single_end(self, healthy):
        sm, ct = healthy
        sm.both_ends_anchored_frac = 0.95
        assert "both_ends_not_above_either_end" in _names(
            run_invariants(sm, ct, _Result()))

    def test_flag_count_above_its_denominator(self, healthy):
        sm, ct = healthy
        sm.n_polya_tail = 5000
        assert "flag_counts_not_above_classified_reads" in _names(
            run_invariants(sm, ct, _Result()))

    def test_unassigned_above_classified(self, healthy):
        sm, ct = healthy
        sm.n_reads_unassigned = 5000
        assert "unassigned_not_above_classified" in _names(
            run_invariants(sm, ct, _Result()))

    def test_cell_count_disagrees_with_the_table(self, healthy):
        sm, ct = healthy
        sm.n_cells = 99
        assert "n_cells_matches_cell_table" in _names(
            run_invariants(sm, ct, _Result()))

    def test_umi_diversity_above_one(self, healthy):
        sm, ct = healthy
        sm.umi_sequence_diversity = {"exonic_sense": 1.7}
        assert "umi_diversity_within_unit_interval" in _names(
            run_invariants(sm, ct, _Result()))

    def test_composition_does_not_sum_to_one(self, healthy):
        sm, ct = healthy
        sm.read_fracs[ReadCategory.EXONIC_SENSE.value] = 0.4
        assert "read_fractions_sum_to_one" in _names(
            run_invariants(sm, ct, _Result()))


class TestSeverity:
    def test_reservoir_estimation_downgrades_conservation_to_a_warning(self, healthy):
        """Promoted totals are scaled from a uniform sub-sample.

        Exact conservation is then not expected to survive rounding, so a
        near-miss is something to note, not a wrong number.
        """
        sm, ct = healthy
        sm.read_fracs[ReadCategory.EXONIC_SENSE.value] = 0.4
        checks = run_invariants(sm, ct, _Result(sampled=True))
        assert "read_fractions_sum_to_one" in _names(checks, WARNING)
        assert "read_fractions_sum_to_one" not in _names(checks, ERROR)

    def test_without_reservoir_estimation_it_is_an_error(self, healthy):
        sm, ct = healthy
        sm.read_fracs[ReadCategory.EXONIC_SENSE.value] = 0.4
        checks = run_invariants(sm, ct, _Result(sampled=False))
        assert "read_fractions_sum_to_one" in _names(checks, ERROR)


class TestAbsentInputs:
    def test_absent_optional_metrics_are_skipped_not_failed(self):
        """An optional metric that was never computed is not a violation."""
        sm = SampleMetrics(sample_name="s", bam_path="b", platform="p",
                           pipeline_stage="x", aligner="a")
        checks = run_invariants(sm, None, None)
        assert failures(checks, ERROR) == []
        names = {c.name for c in checks}
        assert "both_ends_not_above_either_end" not in names
        assert "n_cells_matches_cell_table" not in names
        # the taxonomy check needs no sample data and must still run
        assert "artifact_candidate_categories_subset_of_broad" in names


class TestQCReport:
    def test_report_records_every_check_and_its_numbers(self, healthy, tmp_path):
        sm, ct = healthy
        sm.n_cells = 99
        checks = run_invariants(sm, ct, _Result())
        path = write_qc_report(checks, tmp_path, "sample1")
        assert path is not None and path.name == "sample1.qc.json"
        payload = json.loads(path.read_text())
        assert payload["n_checks"] == len(checks)
        assert payload["n_failed_error"] >= 1
        failed = [c for c in payload["checks"] if not c["ok"]]
        observed = failed[0]["observed"]
        # the observed values are the point: a verdict alone needs a second run
        assert observed and any(v is not None for v in observed.values())

    def test_no_checks_writes_no_file(self, tmp_path):
        assert write_qc_report([], tmp_path, "s") is None
