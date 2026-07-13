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
  INTERGENIC_ENRICHED → passes threshold without more specific evidence

The Poisson background model
-----------------------------
We model the expected read depth at any intergenic position as:

  lambda = (total intergenic reads) / (total intergenic base-pairs)

For a window of size W bp, the expected count is lambda * W.
A locus is "significant" if its observed count has Poisson p < alpha
(Bonferroni-corrected for the number of windows tested).

This is a simple global-rate enrichment screen, not a replacement for a peak
caller: it has no local-background or mappability model.

Hotspot classification heuristics (rule-based, no ML)
------------------------------------------------------
A window is called INTERGENIC_HOTSPOT (internal-priming candidate) if ALL of:
  1. Monoexonic: all reads lack CIGAR N (no splice junctions)
  2. Strand-aware A/T-rich 3′ context within POLYA_CONTEXT_WINDOW bp of the
     modal read 3′ end in the reference
  3. NOT within POLYA_SITE_PROXIMITY bp of an annotated polyA site (otherwise
     routed to INTERGENIC_NOVEL when strand-consistent with enough barcodes)

A window is called INTERGENIC_NOVEL if ALL of:
  1. Passes the Poisson significance threshold
  2. Strong strand consistency (≥ 80% reads on one strand)
  3. ≥ MIN_NOVEL_DISTINCT_BARCODES distinct cell barcodes
  4. At least one read has a CIGAR N (splice evidence) OR is within
     POLYA_SITE_PROXIMITY of an annotated polyA site

Everything else that passes the threshold but satisfies neither rule set is
reported as INTERGENIC_ENRICHED, which deliberately makes no causal claim.

Locus definition
----------------
Reads are assigned by their strand-correct 3' coordinate to pre-defined,
non-overlapping INTERGENIC_LOCUS_WINDOW-bp genomic windows. Strand consistency
is measured after assignment and is therefore independent evidence.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
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
    INTERGENIC_REPEAT_MIN_FRACTION,
    POLYA_CONTEXT_WINDOW,
    POLYA_RUN_MIN_LENGTH,
    POLYA_SITE_PROXIMITY,
    ReadCategory,
)

logger = logging.getLogger(__name__)

# Minimum fraction of reads on the dominant strand to call strand-consistency
STRAND_CONSISTENCY_MIN = 0.80

# Minimum distinct barcodes for an unannotated-transcript candidate
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
    strand_fraction:     float = 0.0
    repeat_overlap_fraction: float = 0.0
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
    n_bases:       int = 0
    umi:           str = ""
    read_length:   int = 0
    read_key:      str = ""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def profile_intergenic_loci(
    records: list[IntergenicReadRecord],
    total_intergenic_bases: int | dict,
    total_barcodes: int,
    *,
    reference=None,           # pysam.FastaFile, optional
    polya_sites: Optional[dict] = None,  # contig → sorted list of positions
    repeat_intervals: Optional[dict] = None,  # contig → list of (start, end)
    alpha: float = ADAPTIVE_PVALUE_THRESHOLD,
    n_test_windows: Optional[int] = None,
) -> tuple[list[IntergenicLocus], list]:
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
    n_test_windows:
        Pre-defined genome-window correction denominator. CLI workflows pass
        every non-mitochondrial BAM window, a conservative upper bound.

    Returns
    -------
    (loci, record_categories)

    loci:
        List of IntergenicLocus objects, one per locus.
    record_categories:
        List of ReadCategory, aligned with ``records`` (same length and order).
        Records in promoted loci carry the new category
        (INTERGENIC_NOVEL / HOTSPOT / REPEAT); everything else stays
        INTERGENIC_SPARSE.  The pipeline uses this to move per-cell counts
        out of INTERGENIC_SPARSE before metrics are computed.
    """
    if not records:
        return [], []

    # ------------------------------------------------------------------
    # 1. Group reads into fixed, pre-defined 3'-end windows. Strand is not part
    #    of the key: it is measured inside each window as independent evidence.
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
    total_bases_value = (
        sum(total_intergenic_bases.values())
        if isinstance(total_intergenic_bases, dict)
        else int(total_intergenic_bases)
    )
    if total_bases_value <= 0:
        logger.warning("total_intergenic_bases is 0 — Poisson background cannot be computed.")
        background_rate = 0.0
    else:
        background_rate = total_intergenic_reads / total_bases_value

    # ------------------------------------------------------------------
    # 3. Minimum barcode threshold (adaptive)
    # ------------------------------------------------------------------
    min_barcodes = max(
        ADAPTIVE_MIN_BARCODES_ABSOLUTE,
        int(total_barcodes * ADAPTIVE_MIN_BARCODES_FRACTION),
    )

    # ------------------------------------------------------------------
    # 4. Score each locus and assign a category per record
    # ------------------------------------------------------------------
    loci: list[IntergenicLocus] = []
    record_categories: list[ReadCategory] = [ReadCategory.INTERGENIC_SPARSE] * len(records)
    # Correct for all fixed intergenic windows, including empty windows. Using
    # only observed, data-derived loci made the original p-values selective.
    if n_test_windows is not None:
        n_tests = max(1, int(n_test_windows))
    elif isinstance(total_intergenic_bases, dict):
        n_tests = max(1, sum(
            int(np.ceil(max(0, bases) / INTERGENIC_LOCUS_WINDOW))
            for bases in total_intergenic_bases.values()
        ))
    else:
        n_tests = max(1, int(np.ceil(total_bases_value / INTERGENIC_LOCUS_WINDOW)))

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
        if locus.category != ReadCategory.INTERGENIC_SPARSE:
            for i in indices:
                record_categories[i] = locus.category

    n_novel = sum(
        1 for locus in loci if locus.category == ReadCategory.INTERGENIC_NOVEL
    )
    n_hotspot = sum(
        1 for locus in loci if locus.category == ReadCategory.INTERGENIC_HOTSPOT
    )
    n_repeat = sum(
        1 for locus in loci if locus.category == ReadCategory.INTERGENIC_REPEAT
    )
    n_enriched = sum(
        1 for locus in loci if locus.category == ReadCategory.INTERGENIC_ENRICHED
    )
    logger.info(
        "  Intergenic windows: %d transcript candidates, %d hotspots, "
        "%d repeat-derived, %d enriched-unresolved, %d sparse",
        n_novel, n_hotspot, n_repeat, n_enriched,
        n_loci - n_novel - n_hotspot - n_repeat - n_enriched,
    )

    return loci, record_categories


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
    modal_end = _modal_three_prime(records)
    start = (modal_end // INTERGENIC_LOCUS_WINDOW) * INTERGENIC_LOCUS_WINDOW
    end = start + INTERGENIC_LOCUS_WINDOW
    n_reads = len(records)

    # Strand consistency
    n_plus  = sum(1 for r in records if r.strand == "+")
    n_minus = n_reads - n_plus
    dominant_strand = "+" if n_plus >= n_minus else "-"
    strand_frac = max(n_plus, n_minus) / n_reads if n_reads > 0 else 0.0
    strand_consistent = strand_frac >= STRAND_CONSISTENCY_MIN

    # Distinct barcodes
    n_barcodes = len({
        r.cell_barcode for r in records
        if r.cell_barcode not in {"", "NO_BARCODE"}
    })

    # Splice evidence
    has_junction   = any(r.has_junction for r in records)
    is_monoexonic  = not has_junction

    # Poisson significance
    expected     = background_rate * INTERGENIC_LOCUS_WINDOW
    raw_pvalue   = float(poisson.sf(n_reads - 1, expected)) if expected > 0 else 1.0
    adj_pvalue   = min(raw_pvalue * n_tests, 1.0)   # Bonferroni
    significant  = (
        adj_pvalue < alpha
        and n_reads >= ADAPTIVE_MIN_READS
    )

    # PolyA context check (requires reference)
    polya_run_downstream = False
    if reference is not None:
        polya_run_downstream = _check_polya_context(
            reference, contig, modal_end, strand=dominant_strand
        )

    # Annotated polyA site proximity
    near_polya = False
    if polya_sites is not None:
        near_polya = _near_polya_site(
            polya_sites, contig, modal_end, strand=dominant_strand
        )

    # RepeatMasker overlap
    repeat_fraction = 0.0
    if repeat_intervals is not None:
        repeat_fraction = _repeat_overlap_fraction(
            repeat_intervals, contig, records
        )
    overlaps_repeat = repeat_fraction >= INTERGENIC_REPEAT_MIN_FRACTION

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
        category = ReadCategory.INTERGENIC_ENRICHED

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
        strand_fraction=strand_frac,
        repeat_overlap_fraction=repeat_fraction,
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

    Requires monoexonic reads, an A-rich 3′ context, AND that the locus is
    NOT within POLYA_SITE_PROXIMITY bp of an annotated polyA site (such a
    locus with strand support is routed to INTERGENIC_NOVEL instead).
    """
    return is_monoexonic and polya_run_downstream and not near_polya_site


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
# Fixed-window assignment
# ---------------------------------------------------------------------------

def _cluster_reads(records: list[IntergenicReadRecord]) -> list[int]:
    """Assign reads to fixed genomic 3'-end windows.

    Fixed windows are defined before observing read density, avoiding the
    selective p-values produced by testing single-linkage clusters discovered
    from the same reads.  Strand is deliberately excluded from the key so the
    dominant-strand fraction remains independent evidence.
    """
    if not records:
        return []

    keys = [
        (r.contig, int(r.three_prime) // INTERGENIC_LOCUS_WINDOW)
        for r in records
    ]
    key_to_id = {
        key: idx + 1 for idx, key in enumerate(sorted(set(keys)))
    }
    return [key_to_id[key] for key in keys]


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

def _check_polya_context(reference, contig: str, position: int, strand: str = "+") -> bool:
    """
    Check for an A-run of >= POLYA_RUN_MIN_LENGTH within POLYA_CONTEXT_WINDOW bp
    of the locus transcript 3′ end (*position*), strand-aware.

    For a + strand locus the A-run sits downstream (higher coordinate) on the
    + strand.  For a - strand locus the transcript 3′ end is the lower
    coordinate and the same A-run reads as a T-run on the + strand just upstream
    (lower coordinate).  Checking only the + strand downstream under-detected
    minus-strand internal priming.
    """
    import re
    try:
        if strand == "-":
            lo = max(0, position - POLYA_CONTEXT_WINDOW)
            if lo >= position:
                return False
            context = reference.fetch(contig, lo, position).upper()
            return bool(re.search(f"T{{{POLYA_RUN_MIN_LENGTH},}}", context))
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
    strand: Optional[str] = None,
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
    keys = [contig]  # legacy BED3/cache representation
    if strand in {"+", "-"}:
        keys = [(contig, strand), (contig, "."), contig]
    elif strand is None:
        keys = [(contig, "+"), (contig, "-"), (contig, "."), contig]

    for key in keys:
        sites = polya_sites.get(key)
        if not sites:
            continue
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


def _repeat_overlap_fraction(
    repeat_intervals: dict,
    contig: str,
    records: list[IntergenicReadRecord],
) -> float:
    """Fraction of sampled aligned record span overlapping repeat union."""
    ivls = repeat_intervals.get(contig, [])
    if not ivls or not records:
        return 0.0
    overlap = 0
    total = 0
    for record in records:
        start, end = int(record.start), int(record.end)
        width = max(end - start, 0)
        total += width
        covered_end = start
        for ivl_start, ivl_end in sorted(ivls):
            if ivl_start >= end:
                break
            lo, hi = max(start, ivl_start), min(end, ivl_end)
            if hi > lo:
                lo = max(lo, covered_end)
                if hi > lo:
                    overlap += hi - lo
                    covered_end = hi
    return overlap / total if total else 0.0


# ---------------------------------------------------------------------------
# Utility: compute total intergenic base-pairs from AnnotationIndex
# ---------------------------------------------------------------------------

def compute_intergenic_bases(
    index,
    *,
    by_contig: bool = False,
    contig_lengths: Optional[dict[str, int]] = None,
) -> int | dict:
    """
    Sum intergenic base-pairs for the Poisson background denominator.

    When BAM contig lengths are supplied, the denominator is the exact
    complement of the union of gene bodies from coordinate 0 to each contig
    end. Without them, the annotation's between-gene gaps are used as a
    compatibility fallback.
    """
    if contig_lengths:
        from scnoisemeter.constants import MITO_CONTIG_NAMES

        frames = []
        for attr in ("gene_bodies_plus", "gene_bodies_minus"):
            obj = getattr(index, attr, None)
            if obj is not None and not obj.df.empty:
                frames.append(obj.df[["Chromosome", "Start", "End"]])
        genes = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
            columns=["Chromosome", "Start", "End"]
        )
        genes_by_contig = {
            str(chrom): sorted(
                (int(row.Start), int(row.End))
                for row in grp.itertuples(index=False)
            )
            for chrom, grp in genes.groupby("Chromosome", observed=False)
        }
        result = {}
        for contig, raw_length in contig_lengths.items():
            if contig in MITO_CONTIG_NAMES:
                continue
            length = max(0, int(raw_length))
            covered = 0
            current_start = current_end = None
            for start, end in genes_by_contig.get(str(contig), []):
                start, end = max(0, start), min(length, end)
                if end <= start:
                    continue
                if current_start is None:
                    current_start, current_end = start, end
                elif start > current_end:
                    covered += current_end - current_start
                    current_start, current_end = start, end
                else:
                    current_end = max(current_end, end)
            if current_start is not None:
                covered += current_end - current_start
            result[str(contig)] = max(0, length - covered)
        return result if by_contig else sum(result.values())

    if index.intergenic is None or index.intergenic.df.empty:
        return {} if by_contig else 0
    df = index.intergenic.df
    if by_contig:
        widths = (df["End"] - df["Start"]).groupby(
            df["Chromosome"], observed=False
        ).sum()
        return {str(chrom): int(width) for chrom, width in widths.items()}
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

    Current tuples contain:
      (contig, start, end, strand, cb, has_junction, three_prime,
       n_bases, umi, read_length, read_key)

    The legacy seven-field format is accepted for backwards compatibility.
    """
    records = []
    for rec in getattr(sample_result, "intergenic_reads", []):
        try:
            if len(rec) == 7:
                contig, start, end, strand, cb, has_junction, three_prime = rec
                n_bases, umi, read_length, read_key = (
                    int(end) - int(start), "", int(end) - int(start), ""
                )
            else:
                (
                    contig, start, end, strand, cb, has_junction, three_prime,
                    n_bases, umi, read_length, read_key,
                ) = rec
            records.append(IntergenicReadRecord(
                contig=contig,
                start=int(start),
                end=int(end),
                strand=strand,
                cell_barcode=str(cb),
                has_junction=bool(has_junction),
                three_prime=int(three_prime),
                n_bases=int(n_bases),
                umi=str(umi),
                read_length=int(read_length),
                read_key=str(read_key),
            ))
        except (ValueError, TypeError) as exc:
            logger.debug("Skipping malformed intergenic record: %s", exc)
    logger.info(
        "Extracted %d intergenic read records for profiling.",
        len(records)
    )
    return records
