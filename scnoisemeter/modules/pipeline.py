"""
pipeline.py
===========
Chromosome-level parallel execution of the read classifier.

Architecture
------------
One worker process per chromosome (or chromosome chunk for very large contigs).
Each worker:
  1. Opens the BAM by region using pysam.fetch(contig)
  2. Instantiates a ReadClassifier (annotation index is shared read-only)
  3. Streams reads, classifies each, and accumulates:
       - Per-cell read-level counts  (cell_barcode → category → n_reads)
       - Per-cell base-level counts  (cell_barcode → category → n_bases)
       - Per-cell UMI sets           (cell_barcode → category → set(umi))
       - Per-cell artifact flags     (cell_barcode → {tso, polya, noncano})
       - Sample-wide length lists    (category → list of read lengths, sampled)
  4. Returns a ContigResult to the main process, which merges all contigs.

Memory model
------------
Per-cell counts are plain dicts-of-dicts, never storing individual reads.
This keeps memory proportional to n_cells × n_categories, not n_reads.

UMI sets are stored per cell per category to compute UMI complexity ratios.
For large experiments (100k cells) this is ~100k × 17 categories × avg_set_size.
At 10 UMIs/cell/category (generous), that's 17M small Python sets = ~1–2 GB.
An optional --no-umi-dedup flag discards UMI sets to reduce memory.

Read-length sampling
--------------------
Storing all read lengths (500M reads × 8 bytes = 4 GB) is infeasible.
We reservoir-sample up to MAX_LENGTH_SAMPLE read lengths per category.
"""

from __future__ import annotations

import bisect
import logging
import random
from collections import defaultdict
from dataclasses import dataclass, field
from multiprocessing import Pool
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Picklable defaultdict factory functions
# (lambdas cannot be pickled by multiprocessing on Python 3.10)
# ---------------------------------------------------------------------------

def _dd_int():
    return defaultdict(int)

def _dd_set():
    return defaultdict(set)

def _dd_flags():
    return {"tso": 0, "polya": 0, "noncanon": 0}

def _dd_list():
    return []

import pysam

from scnoisemeter.constants import (
    BamTag,
    DEFAULT_CHIMERIC_DISTANCE,
    DEFAULT_THREADS,
    LENGTH_BIN_BREAKS,
    MITO_CONTIG_NAMES,
    ReadCategory,
)
from scnoisemeter.modules.annotation import AnnotationIndex
from scnoisemeter.modules.classifier import ReadClassifier
from scnoisemeter.utils.bam_inspector import BamMetadata

logger = logging.getLogger(__name__)

MAX_LENGTH_SAMPLE = 5_000   # per category, per contig worker


def _get_length_bin(read_length: int) -> int:
    """
    Return bin index (0–5) for *read_length*.

    Index mapping (matches LENGTH_BIN_BREAKS = [150, 500, 1000, 2000, 5000]):
      0 → <150 bp
      1 → 150–500 bp
      2 → 500–1000 bp
      3 → 1000–2000 bp
      4 → 2000–5000 bp
      5 → ≥5000 bp
    """
    return bisect.bisect_right(LENGTH_BIN_BREAKS, read_length)


# ---------------------------------------------------------------------------
# Per-contig result
# ---------------------------------------------------------------------------

@dataclass
class ContigResult:
    """
    Aggregated counts produced by one contig worker.

    All dicts keyed by cell_barcode (str).  Empty string = unassigned.
    """
    contig: str

    # read_counts[cb][category] = n_reads
    read_counts:  dict = field(default_factory=lambda: defaultdict(_dd_int))

    # base_counts[cb][category] = n_bases
    base_counts:  dict = field(default_factory=lambda: defaultdict(_dd_int))

    # umi_sets[cb][category] = set of UMI strings
    umi_sets:     dict = field(default_factory=lambda: defaultdict(_dd_set))

    # artifact_flags[cb] = {"tso": int, "polya": int, "noncanon": int}
    artifact_flags: dict = field(default_factory=lambda: defaultdict(_dd_flags))

    # length_samples[category] = list of read lengths (reservoir sample)
    length_samples: dict = field(default_factory=lambda: defaultdict(_dd_list))

    # length_bin_counts[category][bin_idx] = exact read count per length bin.
    # bin_idx 0–5 per _get_length_bin(); no sampling — every read is counted.
    length_bin_counts: dict = field(default_factory=lambda: defaultdict(_dd_int))

    # insert_size_signal / insert_size_noise: reservoir-sampled insert sizes
    # (abs(template_length)) for properly-paired Illumina reads, split by
    # signal (EXONIC_SENSE) vs noise (everything else).  Populated only when
    # paired_end_chimeric=True; empty list otherwise.
    insert_size_signal: list = field(default_factory=list)
    insert_size_noise:  list = field(default_factory=list)

    # intergenic_reads: lightweight side-table for intergenic profiler
    # Each entry: (contig, start, end, strand, cb, has_junction, three_prime)
    intergenic_reads: list = field(default_factory=list)

    # exonic_sense_three_prime: reservoir-sampled 3′ end positions of
    # exonic-sense reads, used for polyA-site-anchored full-length fraction
    exonic_sense_three_prime: list = field(default_factory=list)

    # exonic_sense_five_prime: reservoir-sampled 5′ start positions of
    # exonic-sense reads, used for TSS/CAGE-anchored full-length fraction
    exonic_sense_five_prime: list = field(default_factory=list)

    # Total reads processed (for progress / logging)
    n_reads_processed: int = 0
    n_reads_skipped:   int = 0
    n_reads_skipped_not_called_cell: int = 0


# ---------------------------------------------------------------------------
# Merged sample-wide result
# ---------------------------------------------------------------------------

@dataclass
class SampleResult:
    """
    Merged result across all contigs for one BAM file.

    Identical structure to ContigResult but covers all chromosomes.
    Also stores metadata from the BAM inspector.
    """
    bam_path:  Path
    meta:      BamMetadata

    read_counts:    dict = field(default_factory=lambda: defaultdict(_dd_int))
    base_counts:    dict = field(default_factory=lambda: defaultdict(_dd_int))
    umi_sets:       dict = field(default_factory=lambda: defaultdict(_dd_set))
    artifact_flags: dict = field(default_factory=lambda: defaultdict(_dd_flags))
    length_samples: dict = field(default_factory=lambda: defaultdict(_dd_list))

    # length_bin_counts[category][bin_idx] = exact read count per length bin
    length_bin_counts: dict = field(default_factory=lambda: defaultdict(_dd_int))

    insert_size_signal: list = field(default_factory=list)
    insert_size_noise:  list = field(default_factory=list)

    intergenic_reads: list = field(default_factory=list)
    exonic_sense_three_prime: list = field(default_factory=list)
    exonic_sense_five_prime:  list = field(default_factory=list)

    n_reads_total:     int = 0
    n_reads_processed: int = 0
    n_reads_skipped:   int = 0
    n_reads_skipped_not_called_cell: int = 0


# ---------------------------------------------------------------------------
# Worker function (must be module-level for pickling)
# ---------------------------------------------------------------------------

def _contig_worker(args: dict) -> ContigResult:
    """
    Process all primary alignments on a single contig.

    args keys:
      bam_path, contig, barcode_tag, umi_tag, whitelist,
      chimeric_distance, store_umis, index, reference_path
    """
    bam_path            = args["bam_path"]
    contig              = args["contig"]
    barcode_tag         = args["barcode_tag"]
    umi_tag             = args["umi_tag"]
    whitelist           = args["whitelist"]
    chimeric_dist       = args["chimeric_distance"]
    paired_end_chimeric = args.get("paired_end_chimeric", False)
    store_umis          = args["store_umis"]
    index               = args["index"]
    reference_path      = args.get("reference_path")
    cell_barcodes       = args.get("cell_barcodes")  # set or None

    reference = None
    if reference_path:
        try:
            reference = pysam.FastaFile(reference_path)
        except Exception as exc:
            logger.warning("Could not open reference FASTA (%s): %s — polyA/junction checks disabled.", reference_path, exc)

    classifier = ReadClassifier(
        index,
        barcode_tag=barcode_tag,
        umi_tag=umi_tag,
        whitelist=whitelist,
        chimeric_distance=chimeric_dist,
        paired_end_chimeric=paired_end_chimeric,
        reference=reference,
    )

    result = ContigResult(contig=contig)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # Handle mitochondrial contig name variants
        fetch_contig = contig
        try:
            for read in bam.fetch(fetch_contig):
                result.n_reads_processed += 1

                res = classifier.classify(read)
                if res is None:
                    result.n_reads_skipped += 1
                    continue

                cb  = res.cell_barcode
                cat = res.category

                # --cell-barcodes filter: skip reads not from called cells.
                # Normalise by stripping the trailing -1 suffix that Cell Ranger
                # appends to both CB tags and barcodes.tsv entries.
                if cell_barcodes is not None:
                    cb_norm = cb.removesuffix("-1") if cb else ""
                    if not cb_norm or cb_norm not in cell_barcodes:
                        result.n_reads_skipped_not_called_cell += 1
                        continue

                # Exact length-bin count (no sampling — one increment per read)
                result.length_bin_counts[cat][_get_length_bin(res.read_length)] += 1

                # Read-level counts
                result.read_counts[cb][cat] += 1

                # Base-level counts
                for bc_cat, n in res.base_counts.items():
                    result.base_counts[cb][bc_cat] += n

                # UMI sets
                if store_umis and res.umi:
                    result.umi_sets[cb][cat].add(res.umi)

                # Artifact flags
                if res.is_tso_invasion:
                    result.artifact_flags[cb]["tso"] += 1
                if res.is_polya_priming:
                    result.artifact_flags[cb]["polya"] += 1
                if res.has_noncanonical_junction:
                    result.artifact_flags[cb]["noncanon"] += 1

                # Read-length reservoir sample
                _reservoir_add(
                    result.length_samples[cat],
                    res.read_length,
                    MAX_LENGTH_SAMPLE,
                )

                # Insert size reservoir sample (Illumina paired-end only).
                # Collect from read1 only to avoid double-counting pairs.
                if paired_end_chimeric and read.is_proper_pair and read.is_read1:
                    tlen = abs(read.template_length)
                    if 0 < tlen < 2000:
                        if cat == ReadCategory.EXONIC_SENSE:
                            _reservoir_add(result.insert_size_signal, tlen, 500_000)
                        else:
                            _reservoir_add(result.insert_size_noise, tlen, 500_000)

                # Store intergenic read position for the profiler
                if cat in (ReadCategory.INTERGENIC_SPARSE,
                           ReadCategory.INTERGENIC_REPEAT):
                    import pysam as _pysam
                    result.intergenic_reads.append((
                        res.contig,
                        res.pos,
                        res.pos + res.read_length,
                        "-" if res.is_reverse else "+",
                        cb,
                        res.has_noncanonical_junction,
                        res.pos + res.read_length,
                    ))

                # Store strand-correct 3′ and 5′ ends of exonic-sense reads.
                #
                # For + strand reads: 3′ end = rightmost genomic position
                #                     5′ end = leftmost genomic position (pos)
                # For − strand reads: 3′ end = leftmost genomic position (pos)
                #                     5′ end = rightmost genomic position
                #
                # IMPORTANT: use read.reference_end (the true rightmost genomic
                # coordinate, inclusive of all intron spans) rather than
                # res.pos + res.read_length.  res.read_length is
                # query_alignment_length which counts only aligned query bases
                # (M/I/=/X); for a spliced read "100M5000N100M" it is 200, but
                # reference_end - reference_start is 5200.  Using the wrong
                # value for the rightmost end puts the computed position ~intron-
                # length away from the true 3′/5′ end, causing nearly all
                # spliced reads on the affected strand to fail the 50 bp polyA /
                # 100 bp TSS proximity check.
                if cat == ReadCategory.EXONIC_SENSE:
                    ref_end = read.reference_end or (res.pos + res.read_length)
                    three_prime = (
                        res.pos if res.is_reverse
                        else ref_end
                    )
                    five_prime = (
                        ref_end if res.is_reverse
                        else res.pos
                    )
                    _reservoir_add(
                        result.exonic_sense_three_prime,
                        (res.contig, three_prime),
                        MAX_LENGTH_SAMPLE,
                    )
                    _reservoir_add(
                        result.exonic_sense_five_prime,
                        (res.contig, five_prime),
                        MAX_LENGTH_SAMPLE,
                    )

        except ValueError as exc:
            logger.warning("Could not fetch contig '%s': %s — skipping.", contig, exc)

    if reference:
        reference.close()

    logger.debug("Contig %s: %d reads processed, %d skipped",
                 contig, result.n_reads_processed, result.n_reads_skipped)
    return result


# ---------------------------------------------------------------------------
# Main pipeline entry point
# ---------------------------------------------------------------------------

def run_pipeline(
    bam_path: str | Path,
    index: AnnotationIndex,
    meta: BamMetadata,
    *,
    whitelist: Optional[set] = None,
    cell_barcodes: Optional[set] = None,
    chimeric_distance: int = DEFAULT_CHIMERIC_DISTANCE,
    paired_end_chimeric: bool = False,
    threads: int = DEFAULT_THREADS,
    store_umis: bool = True,
    reference_path: Optional[str] = None,
    contigs: Optional[list] = None,
) -> SampleResult:
    """
    Run the full classification pipeline on *bam_path*.

    Parameters
    ----------
    bam_path:
        Input BAM (must have .bai index).
    index:
        Pre-built annotation index.
    meta:
        BAM metadata from inspect_bam().
    whitelist:
        Set of valid corrected barcode strings.
    chimeric_distance:
        SA-tag chimeric distance threshold.
    threads:
        Number of parallel worker processes.
    store_umis:
        If False, skip UMI set tracking to save memory.
    reference_path:
        Path to reference FASTA for polyA context and junction checks.
    contigs:
        Explicit list of contigs to process.  Default: all contigs in BAM
        header, excluding very small scaffolds (< 1 Mb).
    cell_barcodes:
        Set of called-cell barcode strings (after stripping the -1 suffix).
        If provided, reads whose CB tag (also -1-stripped) is not in this set
        are skipped and counted in n_reads_skipped_not_called_cell.
    """
    bam_path = str(bam_path)

    if contigs is None:
        contigs = _select_contigs(meta)

    logger.info("Processing %d contigs with %d threads …", len(contigs), threads)

    # Build worker argument list
    worker_args = [
        {
            "bam_path":            bam_path,
            "contig":              contig,
            "barcode_tag":         meta.barcode_tag,
            "umi_tag":             meta.umi_tag,
            "whitelist":           whitelist,
            "cell_barcodes":       cell_barcodes,
            "chimeric_distance":   chimeric_distance,
            "paired_end_chimeric": paired_end_chimeric,
            "store_umis":          store_umis,
            "index":               index,
            "reference_path":      reference_path,
        }
        for contig in contigs
    ]

    # Also add mitochondrial contig if present in BAM header
    mito_contigs = [c for c in meta.reference_names if c in MITO_CONTIG_NAMES]
    for mito in mito_contigs:
        if mito not in contigs:
            worker_args.append({
                **worker_args[0],
                "contig": mito,
            })

    # Run workers
    contig_results = []
    if threads == 1:
        # Single-threaded path (easier to debug)
        for args in worker_args:
            contig_results.append(_contig_worker(args))
    else:
        with Pool(processes=threads) as pool:
            contig_results = pool.map(_contig_worker, worker_args)

    # Merge
    sample = SampleResult(bam_path=Path(bam_path), meta=meta)
    for cr in contig_results:
        _merge_contig(sample, cr)

    logger.info(
        "Pipeline complete: %d reads processed, %d skipped, %d skipped (not called cell), "
        "%d unique barcodes",
        sample.n_reads_processed, sample.n_reads_skipped,
        sample.n_reads_skipped_not_called_cell,
        len(sample.read_counts),
    )
    return sample


# ---------------------------------------------------------------------------
# Public merge helper — combine two SampleResults into one (in-place)
# ---------------------------------------------------------------------------

def merge_sample_results(base: SampleResult, other: SampleResult) -> None:
    """
    Merge *other* into *base* in-place.

    Used by run-plate to aggregate per-well SampleResults into a single
    plate-level result before computing metrics.  Read counts, base counts,
    UMI sets, artifact flags, and length samples are all combined.

    *base.bam_path* and *base.meta* are unchanged; callers should set
    these to a meaningful plate-level value before calling compute_metrics.
    """
    base.n_reads_total      += other.n_reads_total
    base.n_reads_processed  += other.n_reads_processed
    base.n_reads_skipped    += other.n_reads_skipped
    base.n_reads_skipped_not_called_cell += other.n_reads_skipped_not_called_cell

    for cb, cat_counts in other.read_counts.items():
        for cat, n in cat_counts.items():
            base.read_counts[cb][cat] += n

    for cb, cat_counts in other.base_counts.items():
        for cat, n in cat_counts.items():
            base.base_counts[cb][cat] += n

    for cb, cat_umis in other.umi_sets.items():
        for cat, umi_set in cat_umis.items():
            base.umi_sets[cb][cat].update(umi_set)

    for cb, flags in other.artifact_flags.items():
        for flag_name, n in flags.items():
            base.artifact_flags[cb][flag_name] += n

    for cat, lengths in other.length_samples.items():
        existing = base.length_samples[cat]
        for L in lengths:
            _reservoir_add(existing, L, MAX_LENGTH_SAMPLE)

    for cat, bin_counts in other.length_bin_counts.items():
        for bin_idx, n in bin_counts.items():
            base.length_bin_counts[cat][bin_idx] += n

    for rec in other.intergenic_reads:
        _reservoir_add(base.intergenic_reads, rec, 500_000)

    for ins in other.insert_size_signal:
        _reservoir_add(base.insert_size_signal, ins, 500_000)
    for ins in other.insert_size_noise:
        _reservoir_add(base.insert_size_noise, ins, 500_000)

    for pos in other.exonic_sense_three_prime:
        _reservoir_add(base.exonic_sense_three_prime, pos, MAX_LENGTH_SAMPLE)
    for pos in other.exonic_sense_five_prime:
        _reservoir_add(base.exonic_sense_five_prime, pos, MAX_LENGTH_SAMPLE)

    # Propagate polyA / TSS site dicts if base doesn't have them yet
    if not getattr(base, "_polya_site_dict", None) and getattr(other, "_polya_site_dict", None):
        base._polya_site_dict = other._polya_site_dict
    if not getattr(base, "_tss_site_dict", None) and getattr(other, "_tss_site_dict", None):
        base._tss_site_dict = other._tss_site_dict


# ---------------------------------------------------------------------------
# Internal merge helper
# ---------------------------------------------------------------------------

def _merge_contig(sample: SampleResult, cr: ContigResult) -> None:
    sample.n_reads_processed += cr.n_reads_processed
    sample.n_reads_skipped   += cr.n_reads_skipped
    sample.n_reads_skipped_not_called_cell += cr.n_reads_skipped_not_called_cell

    for cb, cat_counts in cr.read_counts.items():
        for cat, n in cat_counts.items():
            sample.read_counts[cb][cat] += n

    for cb, cat_bases in cr.base_counts.items():
        for cat, n in cat_bases.items():
            sample.base_counts[cb][cat] += n

    for cb, cat_umis in cr.umi_sets.items():
        for cat, umi_set in cat_umis.items():
            sample.umi_sets[cb][cat].update(umi_set)

    for cb, flags in cr.artifact_flags.items():
        for flag_name, n in flags.items():
            sample.artifact_flags[cb][flag_name] += n

    for cat, lengths in cr.length_samples.items():
        existing = sample.length_samples[cat]
        for L in lengths:
            _reservoir_add(existing, L, MAX_LENGTH_SAMPLE)

    for cat, bin_counts in cr.length_bin_counts.items():
        for bin_idx, n in bin_counts.items():
            sample.length_bin_counts[cat][bin_idx] += n

    # Merge intergenic read records (reservoir-sampled to cap memory)
    for rec in cr.intergenic_reads:
        _reservoir_add(sample.intergenic_reads, rec, 500_000)

    # Merge insert size samples
    for ins in cr.insert_size_signal:
        _reservoir_add(sample.insert_size_signal, ins, 500_000)
    for ins in cr.insert_size_noise:
        _reservoir_add(sample.insert_size_noise, ins, 500_000)

    # Merge exonic-sense 3′ positions
    for pos in cr.exonic_sense_three_prime:
        _reservoir_add(sample.exonic_sense_three_prime, pos, MAX_LENGTH_SAMPLE)
    for pos in cr.exonic_sense_five_prime:
        _reservoir_add(sample.exonic_sense_five_prime, pos, MAX_LENGTH_SAMPLE)

    return


# ---------------------------------------------------------------------------
# Contig selection
# ---------------------------------------------------------------------------

def _select_contigs(meta: BamMetadata, min_length: int = 1_000_000) -> list:
    """
    Return the list of contigs to process in parallel.

    Excludes:
      - Mitochondrial contigs (handled separately)
      - Contigs shorter than min_length (unplaced scaffolds, alt loci)
        to avoid spawning thousands of tiny workers for hg38 alt contigs.

    The remaining contigs are sorted by length descending so the largest
    jobs start first (better load balancing).
    """
    selected = [
        name for name, length in meta.reference_lengths.items()
        if name not in MITO_CONTIG_NAMES and length >= min_length
    ]
    selected.sort(key=lambda c: meta.reference_lengths[c], reverse=True)

    if not selected:
        logger.warning(
            "No contigs ≥ %d bp found — falling back to all non-mito contigs.", min_length
        )
        selected = [c for c in meta.reference_names if c not in MITO_CONTIG_NAMES]

    return selected


# ---------------------------------------------------------------------------
# Reservoir sampling utility
# ---------------------------------------------------------------------------

def _reservoir_add(reservoir: list, value, max_size: int) -> None:
    """
    Add *value* to *reservoir* using reservoir sampling (Algorithm R).
    Keeps the list size ≤ max_size with uniform random sampling.
    """
    if len(reservoir) < max_size:
        reservoir.append(value)
    else:
        idx = random.randint(0, len(reservoir))
        if idx < max_size:
            reservoir[idx] = value
