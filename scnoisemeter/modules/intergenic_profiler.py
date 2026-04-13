"""
intergenic_profiler.py
======================
Post-classification second pass over intergenic reads.

The read classifier marks all intergenic reads as INTERGENIC_SPARSE by
default.  This module re-examines them in aggregate and promotes loci that
pass significance thresholds into more informative sub-categories:

  INTERGENIC_SPARSE   → stays sparse  (below adaptive threshold)
  INTERGENIC_REPEAT   → overlaps RepeatMasker annotation
  INTERGENIC_HOTSPOT  → passes threshold but shows internal-priming signature
  INTERGENIC_NOVEL    → passes threshold, strand-consistent, polyA-site-proximal

The Poisson background model
-----------------------------
We model the expected read depth at any intergenic position as:

  lambda = (total intergenic reads) / (total intergenic base-pairs)

For a window of size W bp, the expected count is lambda * W.
A locus is "significant" if its observed count has Poisson p < alpha
(Bonferroni-corrected for the number of windows tested).

This is conceptually identical to MACS2 peak calling and is well-validated
for this type of sparse signal-above-background problem.

Hotspot classification heuristics (rule-based, no ML)
------------------------------------------------------
A locus is called INTERGENIC_HOTSPOT (internal priming artifact) if ALL of:
  1. Monoexonic: all reads lack CIGAR N (no splice junctions)
  2. A-rich 3′ context: ≥ POLYA_RUN_MIN_LENGTH As within POLYA_CONTEXT_WINDOW
     bp downstream of the modal read 3′ end in the reference
  3. NOT near an annotated polyA site (> POLYA_SITE_PROXIMITY bp from any
     entry in the PolyASite 2.0 / APADB database, if provided)

A locus is called INTERGENIC_NOVEL if ALL of:
  1. Passes the Poisson significance threshold
  2. Strong strand consistency (≥ 80% reads on one strand)
  3. ≥ MIN_NOVEL_DISTINCT_BARCODES distinct cell barcodes
  4. At least one read has a CIGAR N (splice evidence) OR is within
     POLYA_SITE_PROXIMITY of an annotated polyA site

Everything else that passes the threshold but satisfies neither rule set
is reported as INTERGENIC_HOTSPOT (conservative default — flag for review).

Locus definition
----------------
Reads are grouped into loci by merging overlapping read alignments on
the same strand within INTERGENIC_LOCUS_WINDOW bp of each other.
This is a simple single-linkage clustering and is fast enough to run
on the intergenic read subset.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import poisson

from scnoisemeter.constants import (
    ADAPTIVE_MIN_BARCODES_ABSOLUTE,
    ADAPTIVE_MIN_BARCODES_FRACTION,
    ADAPTIVE_MIN_READS,
    ADAPTIVE_PVALUE_THRESHOLD,
    INTERGENIC_LOCUS_WINDOW,
    POLYA_CONTEXT_WINDOW,
    POLYA_RUN_MIN_LENGTH,
    POLYA_SITE_PROXIMITY,
    ReadCategory,
)

logger = logging.getLogger(__name__)

# Minimum fraction of reads on the dominant strand to call strand-consistency
STRAND_CONSISTENCY_MIN = 0.80

# Minimum distinct barcodes for a novel-gene candidate
MIN_NOVEL_DISTINCT_BARCODES = 3


# ---------------------------------------------------------------------------
# Per-locus result
# ---------------------------------------------------------------------------

@dataclass
class IntergenicLocus:
    """
    A cluster of intergenic reads that has been examined and classified.

    Attributes
    ----------
    contig, start, end      Genomic coordinates of the locus (0-based, half-open)
    strand                  Dominant strand ('+', '-', or '.' if unclear)
    n_reads                 Total reads mapping to the locus
    n_barcodes              Distinct cell barcodes
    has_splice_evidence     Any read at this locus has a CIGAR N operation
    is_monoexonic           All reads at this locus lack CIGAR N
    polya_run_downstream    A-rich stretch found downstream of modal 3′ end
    near_polya_site         Within POLYA_SITE_PROXIMITY of an annotated polyA site
    poisson_pvalue          Raw Poisson p-value
    poisson_pvalue_adj      Bonferroni-adjusted p-value
    category                Final sub-category assigned
    """
    contig:              str
    start:               int
    end:                 int
    strand:              str
    n_reads:             int
    n_barcodes:          int
    has_splice_evidence: bool
    is_monoexonic:       bool
    polya_run_downstream: bool
    near_polya_site:     bool
    poisson_pvalue:      float
    poisson_pvalue_adj:  float
    category:            ReadCategory = ReadCategory.INTERGENIC_SPARSE


# ---------------------------------------------------------------------------
# Input record (one per intergenic read, extracted from SampleResult)
# ---------------------------------------------------------------------------

@dataclass
class IntergenicReadRecord:
    """Minimal per-read record needed for intergenic profiling."""
    contig:        str
    start:         int      # reference_start (0-based)
    end:           int      # reference_end   (0-based exclusive)
    strand:        str      # '+' or '-'
    cell_barcode:  str
    has_junction:  bool     # any CIGAR N present
    three_prime:   int      # reference_end (used for polyA context check)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def profile_intergenic_loci(
    records: list[IntergenicReadRecord],
    total_intergenic_bases: int,
    total_barcodes: int,
    *,
    reference=None,           # pysam.FastaFile, optional
    polya_sites: Optional[dict] = None,  # contig → sorted list of positions
    repeat_intervals: Optional[dict] = None,  # contig → list of (start, end)
    alpha: float = ADAPTIVE_PVALUE_THRESHOLD,
) -> tuple[list[IntergenicLocus], dict]:
    """
    Profile intergenic reads and classify their loci.

    Parameters
    ----------
    records:
        All intergenic read records from the sample (extracted from
        SampleResult by the caller).
    total_intergenic_bases:
        Total intergenic base-pairs in the genome (from AnnotationIndex).
        Used as the Poisson background denominator.
    total_barcodes:
        Total number of distinct barcodes in the sample.
        Used to scale the minimum-barcodes threshold.
    reference:
        pysam.FastaFile for polyA context lookup.
    polya_sites:
        Dict mapping contig name → sorted list of known polyA site positions.
        If None, the near_polya_site flag is always False.
    repeat_intervals:
        Dict mapping contig name → sorted list of (start, end) tuples from
        RepeatMasker.  Used to flag INTERGENIC_REPEAT loci.
    alpha:
        Bonferroni-corrected significance threshold.

    Returns
    -------
    (loci, reclassification_map)

    loci:
        List of IntergenicLocus objects, one per significant locus.
    reclassification_map:
        Dict mapping (contig, start, end) → ReadCategory for reads that
        should be re-labelled from INTERGENIC_SPARSE to a richer category.
        The pipeline uses this to update the cell-level count tables.
    """
    if not records:
        return [], {}

    # ------------------------------------------------------------------
    # 1. Group reads into loci by single-linkage clustering per contig/strand
    # ------------------------------------------------------------------
    logger.info("  Clustering %d intergenic reads into loci …", len(records))
    locus_assignments = _cluster_reads(records)   # record_idx → locus_id
    locus_groups: dict[int, list[int]] = {}
    for idx, locus_id in enumerate(locus_assignments):
        locus_groups.setdefault(locus_id, []).append(idx)

    n_loci = len(locus_groups)
    logger.info("  %d intergenic loci identified.", n_loci)

    # ------------------------------------------------------------------
    # 2. Compute Poisson background rate
    # ------------------------------------------------------------------
    total_intergenic_reads = len(records)
    if total_intergenic_bases <= 0:
        logger.warning("total_intergenic_bases is 0 — Poisson background cannot be computed.")
        background_rate = 0.0
    else:
        background_rate = total_intergenic_reads / total_intergenic_bases

    # ------------------------------------------------------------------
    # 3. Minimum barcode threshold (adaptive)
    # ------------------------------------------------------------------
    min_barcodes = max(
        ADAPTIVE_MIN_BARCODES_ABSOLUTE,
        int(total_barcodes * ADAPTIVE_MIN_BARCODES_FRACTION),
    )

    # ------------------------------------------------------------------
    # 4. Score each locus
    # ------------------------------------------------------------------
    loci: list[IntergenicLocus] = []
    n_tests = n_loci  # Bonferroni denominator

    for locus_id, indices in locus_groups.items():
        locus_records = [records[i] for i in indices]
        locus = _score_locus(
            locus_id, locus_records,
            background_rate=background_rate,
            n_tests=n_tests,
            min_barcodes=min_barcodes,
            reference=reference,
            polya_sites=polya_sites,
            repeat_intervals=repeat_intervals,
            alpha=alpha,
        )
        loci.append(locus)

    # ------------------------------------------------------------------
    # 5. Build reclassification map for loci that are not INTERGENIC_SPARSE
    # ------------------------------------------------------------------
    reclassification_map: dict[tuple, ReadCategory] = {}
    for locus in loci:
        if locus.category != ReadCategory.INTERGENIC_SPARSE:
            key = (locus.contig, locus.start, locus.end)
            reclassification_map[key] = locus.category

    n_novel   = sum(1 for l in loci if l.category == ReadCategory.INTERGENIC_NOVEL)
    n_hotspot = sum(1 for l in loci if l.category == ReadCategory.INTERGENIC_HOTSPOT)
    n_repeat  = sum(1 for l in loci if l.category == ReadCategory.INTERGENIC_REPEAT)
    logger.info(
        "  Intergenic loci: %d novel candidates, %d hotspots, %d repeat-derived, %d sparse",
        n_novel, n_hotspot, n_repeat, n_loci - n_novel - n_hotspot - n_repeat,
    )

    return loci, reclassification_map


# ---------------------------------------------------------------------------
# Locus scoring
# ---------------------------------------------------------------------------

def _score_locus(
    locus_id: int,
    records: list[IntergenicReadRecord],
    *,
    background_rate: float,
    n_tests: int,
    min_barcodes: int,
    reference,
    polya_sites: Optional[dict],
    repeat_intervals: Optional[dict],
    alpha: float,
) -> IntergenicLocus:

    contig  = records[0].contig
    start   = min(r.start for r in records)
    end     = max(r.end   for r in records)
    n_reads = len(records)

    # Strand consistency
    n_plus  = sum(1 for r in records if r.strand == "+")
    n_minus = n_reads - n_plus
    dominant_strand = "+" if n_plus >= n_minus else "-"
    strand_frac = max(n_plus, n_minus) / n_reads if n_reads > 0 else 0.0
    strand_consistent = strand_frac >= STRAND_CONSISTENCY_MIN

    # Distinct barcodes
    n_barcodes = len({r.cell_barcode for r in records})

    # Splice evidence
    has_junction   = any(r.has_junction for r in records)
    is_monoexonic  = not has_junction

    # Poisson significance
    locus_width  = max(end - start, 1)
    expected     = background_rate * locus_width
    raw_pvalue   = 1.0 - poisson.cdf(n_reads - 1, expected) if expected > 0 else 1.0
    adj_pvalue   = min(raw_pvalue * n_tests, 1.0)   # Bonferroni
    significant  = (
        adj_pvalue < alpha
        and n_reads >= ADAPTIVE_MIN_READS
        and n_barcodes >= min_barcodes
    )

    # PolyA context check (requires reference)
    polya_run_downstream = False
    if reference is not None:
        modal_end = _modal_three_prime(records)
        polya_run_downstream = _check_polya_context(reference, contig, modal_end)

    # Annotated polyA site proximity
    near_polya = False
    if polya_sites is not None:
        modal_end = _modal_three_prime(records)
        near_polya = _near_polya_site(polya_sites, contig, modal_end)

    # RepeatMasker overlap
    overlaps_repeat = False
    if repeat_intervals is not None:
        overlaps_repeat = _overlaps_repeats(repeat_intervals, contig, start, end)

    # ------------------------------------------------------------------
    # Classification decision tree
    # ------------------------------------------------------------------
    if overlaps_repeat:
        category = ReadCategory.INTERGENIC_REPEAT

    elif not significant:
        category = ReadCategory.INTERGENIC_SPARSE

    elif _is_hotspot(is_monoexonic, polya_run_downstream, near_polya):
        category = ReadCategory.INTERGENIC_HOTSPOT

    elif _is_novel_gene(has_junction, near_polya, strand_consistent, n_barcodes, min_barcodes):
        category = ReadCategory.INTERGENIC_NOVEL

    else:
        # Significant but ambiguous — conservative default: flag as hotspot
        # (better to over-flag artifacts than over-claim novel genes)
        category = ReadCategory.INTERGENIC_HOTSPOT

    return IntergenicLocus(
        contig=contig,
        start=start,
        end=end,
        strand=dominant_strand,
        n_reads=n_reads,
        n_barcodes=n_barcodes,
        has_splice_evidence=has_junction,
        is_monoexonic=is_monoexonic,
        polya_run_downstream=polya_run_downstream,
        near_polya_site=near_polya,
        poisson_pvalue=raw_pvalue,
        poisson_pvalue_adj=adj_pvalue,
        category=category,
    )


# ---------------------------------------------------------------------------
# Classification rule functions (pure, testable)
# ---------------------------------------------------------------------------

def _is_hotspot(
    is_monoexonic: bool,
    polya_run_downstream: bool,
    near_polya_site: bool,
) -> bool:
    """
    Rule-based hotspot (internal priming artifact) classifier.

    Requires monoexonic reads AND A-rich 3′ context.
    The absence of a nearby annotated polyA site increases confidence
    but is not strictly required (the polyA site DB may be incomplete).
    """
    return is_monoexonic and polya_run_downstream


def _is_novel_gene(
    has_splice_evidence: bool,
    near_polya_site: bool,
    strand_consistent: bool,
    n_barcodes: int,
    min_barcodes: int,
) -> bool:
    """
    Rule-based novel gene candidate classifier.

    Requires: strand consistency AND minimum barcode support AND
    at least one of: splice evidence OR annotated polyA site proximity.
    """
    return (
        strand_consistent
        and n_barcodes >= min_barcodes
        and (has_splice_evidence or near_polya_site)
    )


# ---------------------------------------------------------------------------
# Read clustering (single-linkage, by contig + strand + proximity)
# ---------------------------------------------------------------------------

def _cluster_reads(records: list[IntergenicReadRecord]) -> list[int]:
    """
    Assign each read to a locus ID using single-linkage clustering.

    Two reads are in the same locus if they are on the same contig + strand
    and their genomic intervals overlap or are within INTERGENIC_LOCUS_WINDOW bp.

    Returns a list of locus IDs (integers), one per record, in the same order
    as the input records list.

    Algorithm: sort by (contig, strand, start), then sweep with a running
    maximum end position to detect adjacency.  O(n log n).
    """
    if not records:
        return []

    # Sort order: contig, strand, start
    order = sorted(
        range(len(records)),
        key=lambda i: (records[i].contig, records[i].strand, records[i].start),
    )

    locus_ids = [0] * len(records)
    current_locus = 0
    current_end = -1
    current_contig = ""
    current_strand = ""

    for idx in order:
        r = records[idx]
        if (
            r.contig != current_contig
            or r.strand != current_strand
            or r.start > current_end + INTERGENIC_LOCUS_WINDOW
        ):
            current_locus += 1
            current_end = r.end
            current_contig = r.contig
            current_strand = r.strand
        else:
            current_end = max(current_end, r.end)
        locus_ids[idx] = current_locus

    return locus_ids


# ---------------------------------------------------------------------------
# Helper: modal 3′ end position
# ---------------------------------------------------------------------------

def _modal_three_prime(records: list[IntergenicReadRecord]) -> int:
    """Return the most common read 3′ end position in the locus."""
    from collections import Counter
    counts = Counter(r.three_prime for r in records)
    return counts.most_common(1)[0][0]


# ---------------------------------------------------------------------------
# Helper: polyA context check
# ---------------------------------------------------------------------------

def _check_polya_context(reference, contig: str, position: int) -> bool:
    """
    Check for an A-run of >= POLYA_RUN_MIN_LENGTH within
    POLYA_CONTEXT_WINDOW bp downstream of *position*.
    """
    import re
    try:
        context = reference.fetch(contig, position, position + POLYA_CONTEXT_WINDOW).upper()
        return bool(re.search(f"A{{{POLYA_RUN_MIN_LENGTH},}}", context))
    except (ValueError, KeyError):
        return False


# ---------------------------------------------------------------------------
# Helper: annotated polyA site proximity
# ---------------------------------------------------------------------------

def _near_polya_site(
    polya_sites: dict,
    contig: str,
    position: int,
    proximity: int = None,
) -> bool:
    """
    Return True if *position* is within *proximity* bp of any annotated site
    on *contig*.  Defaults to POLYA_SITE_PROXIMITY (50 bp).

    Used for both polyA-site proximity (3' end anchoring) and TSS proximity
    (5' end anchoring) — the latter passes TSS_SITE_PROXIMITY (100 bp).

    polya_sites[contig] must be a sorted list of integer positions.
    Uses binary search (bisect) for O(log n) lookup.
    """
    import bisect
    if proximity is None:
        proximity = POLYA_SITE_PROXIMITY
    sites = polya_sites.get(contig)
    if not sites:
        return False

    idx = bisect.bisect_left(sites, position)
    for candidate_idx in [idx - 1, idx]:
        if 0 <= candidate_idx < len(sites):
            if abs(sites[candidate_idx] - position) <= proximity:
                return True
    return False


# ---------------------------------------------------------------------------
# Helper: RepeatMasker overlap
# ---------------------------------------------------------------------------

def _overlaps_repeats(
    repeat_intervals: dict,
    contig: str,
    start: int,
    end: int,
) -> bool:
    """Return True if [start, end) overlaps any repeat interval on contig."""
    ivls = repeat_intervals.get(contig, [])
    for ivl_start, ivl_end in ivls:
        if ivl_start >= end:
            break
        if ivl_end > start:
            return True
    return False


# ---------------------------------------------------------------------------
# Utility: compute total intergenic base-pairs from AnnotationIndex
# ---------------------------------------------------------------------------

def compute_intergenic_bases(index) -> int:
    """
    Sum the total number of intergenic base-pairs from an AnnotationIndex.
    Used as the Poisson background denominator.
    """
    if index.intergenic is None or index.intergenic.df.empty:
        return 0
    df = index.intergenic.df
    return int((df["End"] - df["Start"]).sum())


# ---------------------------------------------------------------------------
# Utility: extract IntergenicReadRecords from a SampleResult
# ---------------------------------------------------------------------------

def extract_intergenic_records(sample_result) -> list[IntergenicReadRecord]:
    """
    Extract per-read intergenic records from a SampleResult.

    Converts the lightweight side-table stored in SampleResult.intergenic_reads
    (a list of 7-tuples written by the pipeline worker) into typed
    IntergenicReadRecord objects for the profiler.

    Each tuple is: (contig, start, end, strand, cb, has_junction, three_prime)
    """
    records = []
    for rec in getattr(sample_result, "intergenic_reads", []):
        try:
            contig, start, end, strand, cb, has_junction, three_prime = rec
            records.append(IntergenicReadRecord(
                contig=contig,
                start=int(start),
                end=int(end),
                strand=strand,
                cell_barcode=str(cb),
                has_junction=bool(has_junction),
                three_prime=int(three_prime),
            ))
        except (ValueError, TypeError) as exc:
            logger.debug("Skipping malformed intergenic record: %s", exc)
    logger.info(
        "Extracted %d intergenic read records for profiling.",
        len(records)
    )
    return records
