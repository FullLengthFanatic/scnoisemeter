"""
priming_site.py
===============
The priming site and label of flight-polyafilter-bench, ported so that a read
gets its site exactly where the benchmark puts it.

The site is a 0-based boundary at the read's 3' end in transcript sense.  On a
forward read it is the first reference base after the last cDNA base
(``reference_end`` style); on a reverse read the last cDNA base is the base at
the site (``reference_start`` style).  The read strand decides orientation, as
in the rest of the classifier.

Two rules give an observed site; ``site_rule_for`` picks one per platform.

  clip    Short reads.  The alignment end, when the 3' soft clip holds the
          poly(A) or the read carries a ``pa`` tag (``clipped_polya_length``,
          port of ``scripts/01_bd67_preprocess.py:polya_softclip_length``).
  a_run   Long reads.  Where the terminal A run starts, aligned A's included
          (``terminal_a_run``, port of ``sense_3p_clip``, ``inward_pairs``
          and ``terminal_a_run`` in ``pafbench/bd67.py``).

A read without an observed site gets no label.  The label of a read with a
site comes from the A count of the 20 nt downstream of it
(``priming_label``).  This module does not read the genome.
"""

from __future__ import annotations

from typing import Iterator, Optional

import pysam

from scnoisemeter.constants import (
    PRIMING_CLIP_MIN,
    PRIMING_CLIP_PURITY,
    PRIMING_NONTEMPLATED_MAX,
    PRIMING_RUN_MIN,
    PRIMING_RUN_TOLERANCE,
    PRIMING_SITE_RULES,
    PRIMING_TEMPLATED_MIN,
    Platform,
)

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_ALIGNED, _QUERY, _REF = (0, 7, 8), (1, 4), (2, 3)
_CLIP_PLATFORMS = frozenset({
    Platform.ILLUMINA, Platform.ILLUMINA_10X, Platform.ILLUMINA_BD, Platform.SMARTSEQ,
})


def _alignment_end(read: pysam.AlignedSegment) -> int:
    return read.reference_start if read.is_reverse else read.reference_end


def _sense_3p_clip(read: pysam.AlignedSegment) -> str:
    """Soft-clipped bases after the aligned 3' end, in transcript sense."""
    cig, seq = read.cigartuples, read.query_sequence
    if read.is_reverse:
        n = cig[0][1] if cig[0][0] == 4 else 0
        return seq[:n].translate(_COMP)[::-1]
    n = cig[-1][1] if cig[-1][0] == 4 else 0
    return seq[len(seq) - n:] if n else ""


def _inward_pairs(read: pysam.AlignedSegment) -> Iterator[tuple[int, int]]:
    """(query, ref) pairs of the aligned bases from the transcript 3' end inward."""
    cig = read.cigartuples
    if read.is_reverse:
        q, r, step = 0, read.reference_start, 1
    else:
        q, r, step, cig = read.query_length - 1, read.reference_end - 1, -1, reversed(cig)
    for op, n in cig:
        if op in _ALIGNED:
            for k in range(n):
                yield q + step * k, r + step * k
            q, r = q + step * n, r + step * n
        elif op in _QUERY:
            q += step * n
        elif op in _REF:
            r += step * n


def terminal_a_run(read: pysam.AlignedSegment,
                   tolerance: int = PRIMING_RUN_TOLERANCE) -> tuple[int, int, int]:
    """``(site, aligned_a, clip_a)`` of the terminal A run at the read's 3' end.

    The run grows over the 3' soft clip, then inward over the aligned bases,
    with up to ``tolerance`` interior non-A bases from one shared budget; a
    non-A is taken only between A's.  ``site`` is the boundary where the run
    starts, or the alignment end when the run has no aligned base;
    ``aligned_a`` and ``clip_a`` are the A's in its aligned and clipped parts.
    """
    seq = read.query_sequence.translate(_COMP) if read.is_reverse else read.query_sequence
    clip, pairs, left = _sense_3p_clip(read), _inward_pairs(read), []

    def aligned(t):
        """t-th aligned sense base from the 3' end as (base, ref), None past the 5' end."""
        while len(left) <= t:
            p = next(pairs, None)
            if p is None:
                return None
            left.append((seq[p[0]], p[1]))
        return left[t]

    budget, j, n = tolerance, 0, len(clip)
    while j < n:
        if clip[j] != "A":
            if not (budget and j + 1 < n and clip[j + 1] == "A" and (j or aligned(0)[0] == "A")):
                break
            budget -= 1
        j += 1
    t = aligned_a = 0
    while (cur := aligned(t)) is not None:
        if cur[0] == "A":
            aligned_a += 1
        else:
            nxt = aligned(t + 1)
            if not (budget and nxt is not None and nxt[0] == "A" and (t or j)):
                break
            budget -= 1
        t += 1
    if t == 0:
        site = _alignment_end(read)
    else:
        site = left[t - 1][1] + 1 if read.is_reverse else left[t - 1][1]
    return site, aligned_a, clip[:j].count("A")


def clipped_polya_length(read: pysam.AlignedSegment) -> int:
    """A's of the poly(A) clipped at the read's 3' end, 0 below the threshold.

    The ``pa`` tag when present (Cell Ranger's count of the A's it trimmed).
    Otherwise the 3' soft clip, of at least PRIMING_CLIP_MIN bases, must hold
    at least PRIMING_CLIP_MIN A's making up PRIMING_CLIP_PURITY of it; on a
    reverse read that clip is the leading one and its A's read as T's.
    """
    if read.has_tag("pa"):
        return read.get_tag("pa")
    cig = read.cigartuples
    if not cig:
        return 0
    seq = read.query_sequence
    if read.is_reverse:
        if cig[0][0] != 4:
            return 0
        n = cig[0][1]
        if n < PRIMING_CLIP_MIN:
            return 0
        t_count = seq[:n].count("T")
        return t_count if t_count >= PRIMING_CLIP_MIN and t_count >= PRIMING_CLIP_PURITY * n else 0
    if cig[-1][0] != 4:
        return 0
    n = cig[-1][1]
    if n < PRIMING_CLIP_MIN:
        return 0
    a_count = seq[-n:].count("A")
    return a_count if a_count >= PRIMING_CLIP_MIN and a_count >= PRIMING_CLIP_PURITY * n else 0


def observed_site(read: pysam.AlignedSegment, rule: str) -> Optional[int]:
    """The read's priming site under ``rule``, or None when it has none.

    ``clip``: the alignment end when ``clipped_polya_length`` reaches
    PRIMING_CLIP_MIN.  ``a_run``: the run start when the terminal A run holds
    PRIMING_RUN_MIN A's, aligned and clipped together.  A record without a
    sequence has no site.
    """
    if rule not in PRIMING_SITE_RULES:
        raise ValueError(f"unknown priming site rule {rule!r}; expected one of {PRIMING_SITE_RULES}")
    if not read.query_sequence:
        return None
    if rule == "clip":
        return _alignment_end(read) if clipped_polya_length(read) >= PRIMING_CLIP_MIN else None
    site, aligned_a, clip_a = terminal_a_run(read)
    return site if aligned_a + clip_a >= PRIMING_RUN_MIN else None


def priming_label(a_count: int) -> str:
    """Label from the A count of the 20 nt downstream of the site."""
    if a_count >= PRIMING_TEMPLATED_MIN:
        return "templated"
    if a_count <= PRIMING_NONTEMPLATED_MAX:
        return "nontemplated"
    return "ambiguous"


def site_rule_for(platform: Platform) -> str:
    """``clip`` for short-read platforms, ``a_run`` for long reads and UNKNOWN."""
    return "clip" if platform in _CLIP_PLATFORMS else "a_run"
