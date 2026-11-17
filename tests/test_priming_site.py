"""
tests/test_priming_site.py
==========================
The priming-site rules and the label of flight-polyafilter-bench, checked on
the shared synthetic vectors (Appendix A of the harmonisation plan; the ids
are the same in tecap and in the benchmark check).

Contig chrS is 2000 nt of C except the planted bases.  Reads are
transcript-sense on the plus strand.  The minus-strand variant is the mirror:
reverse-complemented read, reversed CIGAR, flag 16,
``pos' = 2000 - reference_end`` and ``site' = 2000 - site``.  Sequences use
the notation of the appendix: "C100 A5 G" is 100 C's, 5 A's and one G.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import NamedTuple, Optional

import pysam
import pytest

from scnoisemeter.constants import (
    PRIMING_CLIP_MIN,
    PRIMING_CLIP_PURITY,
    PRIMING_DEFINITION_SOURCE,
    PRIMING_LABEL_WINDOW,
    PRIMING_NONTEMPLATED_MAX,
    PRIMING_RUN_MIN,
    PRIMING_RUN_TOLERANCE,
    PRIMING_SITE_RULES,
    PRIMING_TEMPLATED_MIN,
    Platform,
)
from scnoisemeter.modules.priming_site import (
    observed_site,
    priming_label,
    site_rule_for,
    terminal_a_run,
)

CONTIG_LEN = 2000
_HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chrS", "LN": CONTIG_LEN}],
})
_COMP = str.maketrans("ACGT", "TGCA")


class LVector(NamedTuple):
    """A row of the long-read table: site is None where the appendix says none."""
    planted: dict
    pos: int
    cigar: str
    seq: str
    site: Optional[int]
    aligned_a: Optional[int] = None
    clip_a: Optional[int] = None
    a_window: Optional[int] = None
    label: Optional[str] = None


class SVector(NamedTuple):
    """A row of the short-read table; pa is the tag value, None for no tag."""
    planted: dict
    pos: int
    cigar: str
    seq: str
    pa: Optional[int]
    site: Optional[int]
    a_window: Optional[int] = None
    label: Optional[str] = None


L_VECTORS = {
    "L1": LVector({}, 900, "100M63S", "C100 A30 G33", 1000, 0, 30, 0, "nontemplated"),
    "L2": LVector({1000: "A10"}, 900, "110M63S", "C100 A10 A30 G33", 1000, 10, 30, 10, "templated"),
    "L3": LVector({1000: "A10"}, 900, "100M", "C100", None),
    **{f"L4_{k}": LVector({1002: f"A{k}"}, 900, "100M63S", "C100 A30 G33", 1000, 0, 30, k, "ambiguous")
       for k in (5, 6, 7)},
    "L5": LVector({}, 900, "100M63H", "C100", None),
    "L6": LVector({}, 900, "100M63S", "C100 A5 G A24 G33", 1000, 0, 29, 0, "nontemplated"),
    "L7": LVector({990: "AAAAGAAAAA"}, 900, "100M53S", "C90 AAAAGAAAAA A20 G33",
                  990, 9, 20, 9, "templated"),
    "L8": LVector({}, 900, "100M63S", "C100 A5 G A3 G A21 G32", None),
    "L9": LVector({990: "A15"}, 900, "105M", "C90 A15", 990, 15, 0, 15, "templated"),
    "L10": LVector({}, 900, "100M42S", "C100 A9 G33", None),
    "L11": LVector({}, 900, "100M64S", "C100 G A30 G33", None),
    # run found, but the window holds 10 nt: no label, not a read with a site
    "L1_at_contig_end": LVector({}, 1890, "100M10S", "C100 A10", 1990, 0, 10),
}

S_VECTORS = {
    "S1": SVector({}, 900, "100M30S", "C100 A30", 30, 1000, 0, "nontemplated"),
    "S2": SVector({1000: "A10"}, 900, "100M40S", "C100 A40", 40, 1000, 10, "templated"),
    "S3": SVector({}, 900, "100M30S", "C100 A30", None, 1000, 0, "nontemplated"),
    "S4": SVector({}, 900, "100M30S", "C100 A20 G10", None, None),
    "S5": SVector({}, 900, "100M19S", "C100 A19", None, None),
    "S6": SVector({}, 900, "100M", "C100", None, None),
    "S7": SVector({}, 900, "100M34S", "C100 AAAAGTG A27", 34, 1000, 0, "nontemplated"),
}


def expand(notation: str) -> str:
    """Bases of an appendix sequence: a letter with a count repeats, anything else is literal."""
    out = []
    for tok in notation.split():
        m = re.fullmatch(r"([ACGT])(\d+)", tok)
        out.append(m[1] * int(m[2]) if m else tok)
    return "".join(out)


def contig(v) -> str:
    """The plus-strand chrS of a vector: all C, planted bases written over."""
    bases = list("C" * CONTIG_LEN)
    for start, notation in v.planted.items():
        planted = expand(notation)
        bases[start:start + len(planted)] = planted
    return "".join(bases)


def segment(pos: int, cigar: str, seq: Optional[str], *, reverse: bool = False,
            pa: Optional[int] = None) -> pysam.AlignedSegment:
    seg = pysam.AlignedSegment(_HEADER)
    seg.query_name = "read"
    seg.query_sequence = seq
    seg.flag = 16 if reverse else 0
    seg.reference_id = 0
    seg.reference_start = pos
    seg.mapping_quality = 60
    seg.cigarstring = cigar
    if pa is not None:
        seg.set_tag("pa", pa)
    return seg


def build(v, minus: bool) -> tuple[pysam.AlignedSegment, Optional[int]]:
    """The vector's read and its expected site; ``minus`` gives the mirrored read."""
    pa = getattr(v, "pa", None)
    read = segment(v.pos, v.cigar, expand(v.seq), pa=pa)
    if not minus:
        return read, v.site
    mirror = segment(CONTIG_LEN - read.reference_end,
                     "".join(reversed(re.findall(r"\d+\D", v.cigar))),
                     read.query_sequence.translate(_COMP)[::-1], reverse=True, pa=pa)
    return mirror, None if v.site is None else CONTIG_LEN - v.site


@pytest.mark.parametrize("name", L_VECTORS)
@pytest.mark.parametrize("minus", [False, True])
def test_a_run_rule(name, minus):
    v = L_VECTORS[name]
    read, want = build(v, minus)
    assert observed_site(read, "a_run") == (want if v.site is not None else None)


@pytest.mark.parametrize("name", [n for n, v in L_VECTORS.items() if v.site is not None])
@pytest.mark.parametrize("minus", [False, True])
def test_run_parts(name, minus):
    v = L_VECTORS[name]
    read, want = build(v, minus)
    assert terminal_a_run(read) == (want, v.aligned_a, v.clip_a)


@pytest.mark.parametrize("name", S_VECTORS)
@pytest.mark.parametrize("minus", [False, True])
def test_clip_rule(name, minus):
    v = S_VECTORS[name]
    read, want = build(v, minus)
    assert observed_site(read, "clip") == (want if v.site is not None else None)


def test_cell_ranger_trim_with_mismatches_keeps_the_alignment_end():      # Review Focus 2
    read, _ = build(S_VECTORS["S7"], False)
    assert observed_site(read, "clip") == 1000 and observed_site(read, "a_run") is None


def test_record_without_sequence_has_no_site():                           # Review Focus 5
    read = segment(pos=900, cigar="100M", seq=None)
    assert observed_site(read, "clip") is None and observed_site(read, "a_run") is None


def test_run_covering_the_whole_alignment():                              # Review Focus 3
    read = segment(pos=990, cigar="15M", seq="A" * 15)
    assert observed_site(read, "a_run") == 990


def test_unknown_rule_is_rejected():
    read, _ = build(L_VECTORS["L1"], False)
    with pytest.raises(ValueError):
        observed_site(read, "run_start")


def test_labels():
    assert [priming_label(k) for k in (0, 4, 5, 7, 8, 20)] == [
        "nontemplated", "nontemplated", "ambiguous", "ambiguous", "templated", "templated"]


@pytest.mark.parametrize("name", [n for n, v in {**L_VECTORS, **S_VECTORS}.items()
                                  if v.site is not None])
def test_label_column_follows_from_the_planted_genome(name):
    v = {**L_VECTORS, **S_VECTORS}[name]
    window = contig(v)[v.site:v.site + PRIMING_LABEL_WINDOW]
    if len(window) < PRIMING_LABEL_WINDOW:
        assert (v.a_window, v.label) == (None, None)
    else:
        assert (window.count("A"), priming_label(window.count("A"))) == (v.a_window, v.label)


def test_constants_equal_the_benchmark_definition():
    d = json.loads((Path(__file__).parent / "data" / "definition.json").read_text())
    short, long_ = d["site"]["short_read"], d["site"]["long_read"]
    assert (PRIMING_LABEL_WINDOW, PRIMING_TEMPLATED_MIN, PRIMING_NONTEMPLATED_MAX) == (
        d["label"]["window"], d["label"]["templated_min"], d["label"]["nontemplated_max"])
    assert (PRIMING_CLIP_MIN, PRIMING_CLIP_PURITY) == (short["clip_min"], short["clip_purity_min"])
    assert (PRIMING_RUN_MIN, PRIMING_RUN_TOLERANCE) == (long_["run_min"], long_["run_tolerance"])
    assert PRIMING_DEFINITION_SOURCE == f'{d["source"]} {d["version"]}'


_RULE_BY_PLATFORM = {
    Platform.ONT: "a_run",
    Platform.PACBIO: "a_run",
    Platform.ILLUMINA: "clip",
    Platform.ILLUMINA_10X: "clip",
    Platform.ILLUMINA_BD: "clip",
    Platform.SMARTSEQ: "clip",
    Platform.UNKNOWN: "a_run",
}


@pytest.mark.parametrize("platform", list(Platform))
def test_site_rule_for(platform):
    assert site_rule_for(platform) == _RULE_BY_PLATFORM[platform]
    assert site_rule_for(platform) in PRIMING_SITE_RULES
