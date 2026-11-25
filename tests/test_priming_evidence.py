"""
tests/test_priming_evidence.py
==============================
The 0.9 internal-priming work, transferred from flightcollapse.

Each class pins one behaviour that a plausible future edit would break, and the
option defaults are asserted explicitly: every setting here changes a reported
number, so a default moves only in a release whose CHANGELOG names the numbers
it changes.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pysam
import pytest
from click.testing import CliRunner

sys.path.insert(0, str(Path(__file__).parent.parent))

from test_priming_site import (  # noqa: E402
    CONTIG_LEN,
    L_VECTORS,
    S_VECTORS,
    build,
    contig,
    segment,
)

from scnoisemeter.cli import (  # noqa: E402
    _build_per_well_metrics_row,
    _load_polya_sites,
    _polya_atlas_density,
    _warn_if_atlas_too_dense,
    cli,
)
from scnoisemeter.constants import (  # noqa: E402
    POLYA_MIN_SUPPORT_DEFAULT,
    POLYA_TAIL_MIN_LENGTH,
    PRIMING_CRITERION_DEFAULT,
    ReadCategory,
)
from scnoisemeter.modules.classifier import ReadClassifier  # noqa: E402
from scnoisemeter.modules.metrics import SampleMetrics  # noqa: E402
from scnoisemeter.modules.priming_site import observed_site  # noqa: E402


class _Read:
    """Minimal stand-in: the tail measurement reads only these three fields."""

    def __init__(self, cigartuples, query_sequence, is_reverse=False):
        self.cigartuples = cigartuples
        self.query_sequence = query_sequence
        self.is_reverse = is_reverse


@pytest.fixture
def classifier():
    return ReadClassifier.__new__(ReadClassifier)


# ---------------------------------------------------------------------------
# Clipped 3' polyA tail from the read's own soft clip
# ---------------------------------------------------------------------------

class TestPolyATail:
    def test_soft_clipped_a_run_is_a_tail(self, classifier):
        read = _Read([(0, 50), (4, 12)], "G" * 50 + "A" * 12)
        clip, run, frac = classifier._polya_tail_stats(read)
        assert (clip, run) == (12, 12)
        assert frac == pytest.approx(1.0)
        assert classifier._has_polya_tail(read) is True

    def test_short_run_is_not_a_tail(self, classifier):
        read = _Read([(0, 50), (4, 4)], "G" * 50 + "A" * 4)
        assert classifier._polya_tail_stats(read)[1] < POLYA_TAIL_MIN_LENGTH
        assert classifier._has_polya_tail(read) is False

    def test_hard_clip_reports_no_tail_rather_than_wrong_bases(self, classifier):
        """Hard-clipped bases are absent from query_sequence.

        Indexing the aligned sequence as though it were the clip would read
        real alignment bases as a tail, so a hard-clipped BAM must report 0.
        """
        read = _Read([(0, 50), (5, 12)], "G" * 50)
        assert classifier._polya_tail_stats(read) == (0, 0, 0.0)
        assert classifier._has_polya_tail(read) is False

    def test_reverse_strand_tail_is_the_leading_clip(self, classifier):
        """On a minus-strand alignment the tail precedes the block as poly-T."""
        read = _Read([(4, 12), (0, 50)], "T" * 12 + "G" * 50, is_reverse=True)
        assert classifier._polya_tail_stats(read)[1] == 12
        assert classifier._has_polya_tail(read) is True

    def test_reverse_strand_ignores_the_trailing_clip(self, classifier):
        read = _Read([(0, 50), (4, 12)], "G" * 50 + "A" * 12, is_reverse=True)
        assert classifier._has_polya_tail(read) is False

    def test_adapter_after_the_tail_still_counts(self, classifier):
        """The run decides, not the clip's A-fraction.

        A genuine 12 bp tail inside a 60 bp clip that also holds adapter or TSO
        sequence has an A-fraction near 0.2; requiring both would silently
        demand a ~45 bp tail.
        """
        read = _Read([(0, 50), (4, 60)], "G" * 50 + "A" * 12 + "C" * 48)
        clip, run, frac = classifier._polya_tail_stats(read)
        assert (clip, run) == (60, 12)
        assert frac < 0.25
        assert classifier._has_polya_tail(read) is True


# ---------------------------------------------------------------------------
# Atlas confidence filtering and density
# ---------------------------------------------------------------------------

def _write_atlas(tmp_path: Path) -> Path:
    """Three PolyASite-shaped records: support 5.6, 44.4 and 88.9 percent."""
    rows = [
        ("chr1", "1000", "1001", "chr1:1001:+", "0.0001", "+", "5.6", "1", "20", "DI", "NaN"),
        ("chr1", "5000", "5001", "chr1:5001:+", "0.0500", "+", "44.4", "1", "80", "TE", "AAUAAA"),
        ("chr1", "9000", "9001", "chr1:9001:-", "0.2000", "-", "88.9", "1", "99", "TE", "AUUAAA"),
    ]
    path = tmp_path / "atlas.bed.gz"
    with gzip.open(path, "wt") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return path


class TestAtlasSupportFilter:
    def test_unfiltered_keeps_every_record(self, tmp_path):
        sites = _load_polya_sites(
            str(_write_atlas(tmp_path)), chrom_style="ucsc", min_support=0.0
        )
        assert sum(len(v) for v in sites.values()) == 3

    def test_default_drops_rarely_supported_records(self, tmp_path):
        sites = _load_polya_sites(str(_write_atlas(tmp_path)), chrom_style="ucsc")
        assert sum(len(v) for v in sites.values()) == 2
        assert 1000 not in sites.get(("chr1", "+"), [])

    def test_threshold_drops_low_support_records(self, tmp_path):
        sites = _load_polya_sites(
            str(_write_atlas(tmp_path)), chrom_style="ucsc", min_support=25.0
        )
        assert sum(len(v) for v in sites.values()) == 2
        assert 1000 not in sites.get(("chr1", "+"), [])

    def test_plain_bed3_is_unaffected_by_the_threshold(self, tmp_path):
        path = tmp_path / "plain.bed"
        path.write_text("chr1\t100\t101\nchr1\t200\t201\n")
        sites = _load_polya_sites(str(path), chrom_style="ucsc", min_support=50.0)
        assert sum(len(v) for v in sites.values()) == 2

    def test_pas_hexamer_column_is_read(self, tmp_path):
        pas: dict = {}
        sites = _load_polya_sites(
            str(_write_atlas(tmp_path)), chrom_style="ucsc", min_support=0.0, pas_out=pas
        )
        assert sum(len(v) for v in sites.values()) == 3
        # only the two records naming a hexamer
        assert sum(len(v) for v in pas.values()) == 2
        assert 1000 not in pas.get(("chr1", "+"), [])

    def test_pas_subset_survives_an_earlier_load_without_it(self, tmp_path):
        """discover and run-plate workers load without pas_out and write the cache."""
        atlas = str(_write_atlas(tmp_path))
        _load_polya_sites(atlas, chrom_style="ucsc")
        pas: dict = {}
        _load_polya_sites(atlas, chrom_style="ucsc", pas_out=pas)
        assert sum(len(v) for v in pas.values()) == 2

    def test_option_defaults(self):
        assert POLYA_MIN_SUPPORT_DEFAULT == 25.0
        assert PRIMING_CRITERION_DEFAULT == "count"


class TestAtlasDensity:
    def test_density_is_measured_against_the_reference(self):
        # two sites 50 bp wide either side on a 10 kb contig -> ~2%
        sites = {("chr1", "+"): [1000, 5000]}
        frac = _polya_atlas_density(sites, {"chr1": 10_000}, proximity=50)
        assert frac == pytest.approx(200 / 10_000, abs=1e-6)

    def test_overlapping_windows_are_not_double_counted(self):
        sites = {("chr1", "+"): [1000, 1010, 1020]}
        frac = _polya_atlas_density(sites, {"chr1": 10_000}, proximity=50)
        assert frac == pytest.approx(120 / 10_000, abs=1e-6)

    def test_dense_atlas_records_a_warning(self):
        class _Meta:
            reference_lengths = {"chr1": 1000}
            warnings: list = []

        meta = _Meta()
        meta.warnings = []
        dense = {("chr1", "+"): list(range(0, 1000, 20))}
        _warn_if_atlas_too_dense(dense, meta, 0.0)
        assert meta.warnings and "polyA site" in meta.warnings[0]

    def test_sparse_atlas_is_silent(self):
        class _Meta:
            reference_lengths = {"chr1": 10_000_000}
            warnings: list = []

        meta = _Meta()
        meta.warnings = []
        _warn_if_atlas_too_dense({("chr1", "+"): [5000]}, meta, 25.0)
        assert meta.warnings == []


class TestPrimingCriterionOption:
    def test_unknown_criterion_is_rejected(self):
        from scnoisemeter.modules.annotation import AnnotationIndex
        index = AnnotationIndex.__new__(AnnotationIndex)
        with pytest.raises(ValueError, match=r"a_run\|perc_a\|either\|count"):
            ReadClassifier(index, priming_criterion="percent")
        assert ReadClassifier(index, priming_criterion="count").priming_criterion == "count"

    def test_unknown_site_rule_is_rejected(self):
        with pytest.raises(ValueError, match=r"clip\|a_run"):
            ReadClassifier(MagicMock(splice_sites={}), priming_site_rule="auto")

    def test_run_pipeline_without_a_site_rule_takes_the_platform_rule(self, tmp_path, monkeypatch):
        # A Python caller that names no rule gets what --priming-site-rule auto gives
        from scnoisemeter.constants import Platform
        from scnoisemeter.modules import pipeline
        from scnoisemeter.utils.bam_inspector import BamMetadata
        seen = []

        def worker(args):
            seen.append(args["priming_site_rule"])
            return pipeline.ContigResult(contig=args["contig"])

        monkeypatch.setattr(pipeline, "_contig_worker", worker)
        meta = BamMetadata(path=tmp_path / "10x.bam")
        meta.platform = Platform.ILLUMINA_10X
        pipeline.run_pipeline(tmp_path / "10x.bam", MagicMock(), meta, threads=1, contigs=["chr1"])
        assert seen == ["clip"]


# ---------------------------------------------------------------------------
# Per-read priming label and the count criterion, on the Appendix A vectors
# ---------------------------------------------------------------------------

_COMP = str.maketrans("ACGT", "TGCA")


@pytest.fixture
def classifier_for(tmp_path):
    """``(name, minus, criterion, rule)`` -> a classifier on the vector's contig, and its read.

    The FASTA holds chrS as the appendix plants it, or its reverse complement
    for the mirrored (minus) read.
    """
    opened = []

    def make(name, minus, criterion="count", rule="a_run"):
        v = {**L_VECTORS, **S_VECTORS}[name]
        seq = contig(v).translate(_COMP)[::-1] if minus else contig(v)
        fasta = tmp_path / f"{name}_{'minus' if minus else 'plus'}.fa"
        fasta.write_text(f">chrS\n{seq}\n")
        pysam.faidx(str(fasta))
        opened.append(pysam.FastaFile(str(fasta)))
        clf = ReadClassifier(MagicMock(splice_sites={}), reference=opened[-1],
                             priming_criterion=criterion, priming_site_rule=rule)
        return clf, build(v, minus)[0]

    yield make
    for fa in opened:
        fa.close()


def blank_result(clf, read):
    return clf._make_result(read, "", "", ReadCategory.EXONIC_SENSE, {}, False)


@pytest.mark.parametrize("minus", [False, True])
@pytest.mark.parametrize("name,rule,label", [
    ("L1", "a_run", "nontemplated"), ("L2", "a_run", "templated"),
    ("L4_6", "a_run", "ambiguous"), ("L7", "a_run", "templated"), ("L9", "a_run", "templated"),
    ("S1", "clip", "nontemplated"), ("S2", "clip", "templated"),
    ("L6", "a_run", "nontemplated"), ("S3", "clip", "nontemplated"),
    ("S7", "clip", "nontemplated"),                 # Review Focus 2
])
def test_label_on_observed_site(name, rule, label, minus, classifier_for):
    clf, read = classifier_for(name, minus, criterion="count", rule=rule)
    assert clf._priming_label(read) == (True, label)


@pytest.mark.parametrize("name,rule", [
    ("L3", "a_run"), ("L5", "a_run"), ("L10", "a_run"), ("S4", "clip"), ("S6", "clip"),
])
def test_no_site_no_label(name, rule, classifier_for):
    clf, read = classifier_for(name, False, criterion="count", rule=rule)
    assert clf._priming_label(read) == (False, None)


def test_trimmed_read_gets_the_same_site_and_label_as_the_untrimmed_one(classifier_for):
    # L2 keeps its poly(A); S2 is the same molecule after Cell Ranger trimming
    kept_clf, kept = classifier_for("L2", False, criterion="count", rule="a_run")
    trim_clf, trimmed = classifier_for("S2", False, criterion="count", rule="clip")
    assert observed_site(kept, "a_run") == observed_site(trimmed, "clip") == 1000
    assert kept_clf._priming_label(kept) == trim_clf._priming_label(trimmed) == (True, "templated")


def test_count_flag_is_the_templated_label(classifier_for):
    clf, read = classifier_for("L2", False, "count", "a_run")
    res = clf._decorate_independent_flags(read, blank_result(clf, read))
    assert res.is_polya_priming and res.priming_site_observed and res.priming_label == "templated"


@pytest.mark.parametrize("minus", [False, True])
@pytest.mark.parametrize("name,rule,label", [
    ("L1", "a_run", "nontemplated"), ("L4_6", "a_run", "ambiguous"), ("S1", "clip", "nontemplated"),
])
def test_count_flag_leaves_other_labels_unflagged(name, rule, label, minus, classifier_for):
    clf, read = classifier_for(name, minus, "count", rule)
    res = clf._decorate_independent_flags(read, blank_result(clf, read))
    assert (res.priming_site_observed, res.priming_label, res.is_polya_priming) == (True, label, False)


def test_old_criteria_are_untouched(classifier_for):
    clf, read = classifier_for("L2", False, "a_run", "a_run")
    res = clf._decorate_independent_flags(read, blank_result(clf, read))
    assert res.is_polya_priming is False            # alignment end 1010, no A-run downstream: the anchor conflict of the brief
    assert res.priming_label == "templated"         # the label is still reported


def test_no_reference_means_not_observable():        # Review Focus 5
    clf = ReadClassifier(MagicMock(splice_sites={}), reference=None, priming_criterion="count")
    assert clf._priming_label(build(L_VECTORS["L2"], False)[0]) == (False, None)


@pytest.mark.parametrize("minus", [False, True])
def test_site_near_contig_end_has_no_label(minus, classifier_for):   # Review Focus 3
    # site at 1990 of 2000; the mirror puts it at 10, a window below the contig start
    clf, read = classifier_for("L1_at_contig_end", minus, "count", "a_run")
    assert clf._priming_label(read) == (False, None)


def test_read_aligned_entirely_to_adenines_is_labelled(classifier_for):   # Review Focus 3
    # L9's genome: 15 A's at 990.  The read is those 15 A's only, so its site
    # is its first aligned base and the window holds them all.
    clf, _ = classifier_for("L9", False, "count", "a_run")
    assert clf._priming_label(segment(990, "15M", "A" * 15)) == (True, "templated")


@pytest.mark.parametrize("rule", ["a_run", "clip"])
def test_record_without_cigar_has_no_label(rule, classifier_for):
    # A mapped record with CIGAR '*' was classified before the label existed;
    # the label is computed on every read, so it must not raise on one.
    clf, _ = classifier_for("L1", False, "a_run", rule)
    read = segment(900, "100M30S", "C" * 100 + "A" * 30)
    read.cigartuples = None
    assert clf._priming_label(read) == (False, None)


# ---------------------------------------------------------------------------
# Label counts and shares over the reads with an observed site, through `run`
# ---------------------------------------------------------------------------

LABEL_COUNTS = ("n_priming_site_observed", "n_priming_templated",
                "n_priming_nontemplated", "n_priming_ambiguous")
LABEL_FRACS = ("priming_templated_frac", "priming_nontemplated_frac",
               "priming_ambiguous_frac")

# L1 non-templated, L2 templated, L4_6 ambiguous, L5 hard-clipped (no site)
ONE_CONTIG = {"chrV": ["L1", "L2", "L4_6", "L5"]}
TWO_CONTIGS = {"chrV1": ["L1", "L2"], "chrV2": ["L4_6", "L5"]}


def _vector_library(tmp_path: Path, layout: dict) -> tuple[Path, Path]:
    """BAM and FASTA of *layout*, ``{contig: [vector id, ...]}``.

    Each vector takes a CONTIG_LEN block of its contig, planted as the
    appendix plants chrS, and its plus-strand read moves with the block.
    """
    tmp_path.mkdir(parents=True, exist_ok=True)
    fasta = tmp_path / "vectors.fa"
    fasta.write_text("".join(
        f">{name}\n{''.join(contig(L_VECTORS[v]) for v in vectors)}\n"
        for name, vectors in layout.items()
    ))
    pysam.faidx(str(fasta))
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": CONTIG_LEN * len(vectors)}
               for name, vectors in layout.items()],
    })
    bam = tmp_path / "vectors.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for name, vectors in layout.items():
            for block, v in enumerate(vectors):
                fields = build(L_VECTORS[v], False)[0].to_dict()
                fields.update(name=v, ref_name=name,
                              ref_pos=str(int(fields["ref_pos"]) + block * CONTIG_LEN))
                out.write(pysam.AlignedSegment.from_dict(fields, header))
    pysam.index(str(bam))
    return bam, fasta


def _run_vectors(tiny_gtf, tmp_path: Path, layout: dict, *, reference: bool,
                 threads: int = 1) -> tuple[dict, dict]:
    """``run --platform ont`` (site rule a_run) on the layout: read_metrics rows, MultiQC data."""
    bam, fasta = _vector_library(tmp_path, layout)
    polya = tmp_path / "polya.bed"
    polya.write_text("chr1\t5990\t6010\tsite\t1\t+\n")
    out = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "run", "--bam", str(bam), "--gtf", str(tiny_gtf),
        "--polya-sites", str(polya), "--tss-db", "none",
        "--platform", "ont", "--library-strand", "stranded",
        "--threads", str(threads), "--no-cache", "--output-dir", str(out),
        *(["--reference", str(fasta)] if reference else []),
    ])
    assert result.exit_code == 0, result.output
    rows = (out / "vectors.read_metrics.tsv").read_text().splitlines()
    multiqc = json.loads((out / "vectors.multiqc.json").read_text())
    return dict(line.split("\t", 1) for line in rows), multiqc["data"]["vectors"]


def test_label_counts_are_over_reads_with_an_observed_site(tiny_gtf, tmp_path):
    metrics, multiqc = _run_vectors(tiny_gtf, tmp_path, ONE_CONTIG, reference=True)
    assert metrics["n_reads_classified"] == "4"           # L5 is classified, it has no site
    assert [metrics[k] for k in LABEL_COUNTS] == ["3", "1", "1", "1"]
    assert [metrics[k] for k in LABEL_FRACS] == ["0.333333"] * 3
    assert multiqc["n_priming_site_observed"] == 3
    assert [multiqc[k] for k in LABEL_FRACS] == [0.33333333] * 3
    assert metrics["n_polya_priming"] == metrics["n_priming_templated"] == "1"
    assert multiqc["polya_priming_frac"] == 0.25


def test_label_counts_are_not_observable_without_reference(tiny_gtf, tmp_path):   # Review Focus 5
    metrics, multiqc = _run_vectors(tiny_gtf, tmp_path, ONE_CONTIG, reference=False)
    assert [metrics[k] for k in LABEL_COUNTS] == [""] * 4
    assert [metrics[k] for k in LABEL_FRACS] == ["not observable"] * 3
    assert [k for k in multiqc if "priming" in k] == []


def test_no_read_with_a_site_leaves_the_shares_not_observable(tiny_gtf, tmp_path):
    # L3 ends right before its A10, L5 is hard-clipped: neither has a site
    metrics, multiqc = _run_vectors(tiny_gtf, tmp_path, {"chrV": ["L3", "L5"]}, reference=True)
    assert [metrics[k] for k in LABEL_COUNTS] == ["0"] * 4
    assert [metrics[k] for k in LABEL_FRACS] == ["not observable"] * 3
    assert multiqc["n_priming_site_observed"] == 0
    assert [multiqc[k] for k in LABEL_FRACS] == [None] * 3


@pytest.mark.parametrize("threads", [1, 2])
def test_label_counts_merge_across_contigs(threads, tiny_gtf, tmp_path):
    metrics, multiqc = _run_vectors(tiny_gtf, tmp_path, TWO_CONTIGS, reference=True,
                                    threads=threads)
    assert [metrics[k] for k in LABEL_COUNTS] == ["3", "1", "1", "1"]
    assert [multiqc[k] for k in LABEL_FRACS] == [0.33333333] * 3


@pytest.mark.parametrize("value", [None, 3])
def test_per_well_row_carries_the_label_counts(value):
    sm = SampleMetrics(sample_name="w", bam_path="w.bam", platform="smartseq",
                       pipeline_stage="pre_filter", aligner="STAR")
    for key in LABEL_COUNTS:
        setattr(sm, key, value)
    row = _build_per_well_metrics_row(sm, plate_id="P1", well_id="A1",
                                      bam_path="w.bam", well_meta_from_sheet=None)
    assert [row[k] for k in LABEL_COUNTS] == [value] * 4


# ---------------------------------------------------------------------------
# --anchor-veto-internal-priming reads the sixth endpoint field
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("veto,anchored", [(False, 1.0), (True, 0.75)])
def test_anchor_veto_drops_only_the_flagged_read(veto, anchored, tmp_path):
    from scnoisemeter.modules.metrics import compute_metrics
    from scnoisemeter.modules.pipeline import SampleResult
    from scnoisemeter.utils.bam_inspector import BamMetadata
    result = SampleResult(tmp_path / "s.bam", BamMetadata(path=tmp_path / "s.bam"))
    # four reads end at one poly(A) site; the sixth field is the per-read flag
    result.exonic_sense_endpoints.extend([
        ("chr1", "+", 1000, 5000, "templated", True),
        ("chr1", "+", 1100, 5000, "nontemplated", False),
        ("chr1", "+", 1200, 5000, "no_site", False),
        ("chr1", "+", 1300, 5000, "five_fields"),       # written before 0.9
    ])
    result._polya_site_dict = {("chr1", "+"): [5000]}
    result._anchor_veto_internal_priming = veto
    sm, _ = compute_metrics(result, "s")
    assert sm.three_prime_anchored_frac == anchored
