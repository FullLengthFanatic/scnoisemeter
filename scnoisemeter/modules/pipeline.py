"""
pipeline.py
===========
Chromosome-level parallel execution of the read classifier.

Architecture
------------
One worker process per chromosome (or chromosome chunk for very large contigs).
Each worker:
  1. Opens the BAM by region using pysam.fetch(contig)
  2. Instantiates a ReadClassifier over the annotation index, which is passed
     to the pool once via the initializer rather than per task
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

UMI sets are stored per cell per category to compute UMI sequence-diversity ratios.
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
import zlib
from collections import defaultdict
from dataclasses import dataclass, field
from multiprocessing import Pool
from pathlib import Path
from typing import Optional

import pysam

from scnoisemeter.constants import (
    DEFAULT_CHIMERIC_DISTANCE,
    DEFAULT_THREADS,
    LENGTH_BIN_BREAKS,
    LOW_MAPQ_THRESHOLD,
    TSO_MIN_MATCH_LENGTH,
    ReadCategory,
)
from scnoisemeter.modules.annotation import AnnotationIndex
from scnoisemeter.modules.classifier import ReadClassifier
from scnoisemeter.utils.bam_inspector import BamMetadata

# ---------------------------------------------------------------------------
# Picklable defaultdict factory functions
# (lambdas cannot be pickled by multiprocessing on Python 3.10)
# ---------------------------------------------------------------------------

def _dd_int():
    return defaultdict(int)

def _dd_set():
    return defaultdict(set)

def _dd_flags():
    return {
        "tso": 0,
        "polya": 0,
        "noncanon": 0,
        "tso_concat": 0,
        "discordant_pair": 0,
        "low_mapq": 0,
    }


def _dd_alignment():
    return {
        "mapq_sum": 0,
        "nm_sum": 0,
        "nm_observed": 0,
        "softclip_bases": 0,
        "query_bases": 0,
    }


class _Reservoir(list):
    """
    Fixed-size list that tracks the number of items offered for true Algorithm R
    sampling. Behaves as a plain list to consumers (iteration, len, indexing).
    """
    def __init__(self, max_size: int = 0):
        super().__init__()
        self.max_size = max_size
        self.n_seen = 0


def _dd_list():
    # MAX_LENGTH_SAMPLE is defined below; resolved at call time.
    return _Reservoir(MAX_LENGTH_SAMPLE)


def _make_intergenic_reservoir():
    return _Reservoir(500_000)


def _make_insert_size_reservoir():
    return _Reservoir(500_000)


def _make_position_reservoir():
    return _Reservoir(MAX_LENGTH_SAMPLE)


def _make_endpoint_reservoir():
    return _Reservoir(MAX_LENGTH_SAMPLE)

logger = logging.getLogger(__name__)

MAX_LENGTH_SAMPLE = 5_000   # per category, per contig worker

# The annotation index is large (~75 MB pickled on GENCODE v47 with no repeats
# BED, several times that with one) and read-only once built.  It used to travel
# inside every per-contig task dict, so Pool.map serialised and re-parsed a
# fresh copy for each contig -- ~1.9 GB of transfer and ~110 s of pure pickling
# on a 25-contig run, and worse under spawn where nothing is shared by fork.
# It is now handed to the workers once through the pool initializer and read
# from this module global; task dicts may still carry an explicit "index" for
# direct callers and tests.
_WORKER_INDEX: Optional[AnnotationIndex] = None


def _init_worker(index: AnnotationIndex) -> None:
    """Pool initializer: bind the shared annotation index in each worker."""
    global _WORKER_INDEX
    _WORKER_INDEX = index


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


def _read_key(read: pysam.AlignedSegment) -> str:
    """Stable key that distinguishes the two primary mates of a template."""
    suffix = "/1" if read.is_read1 else "/2" if read.is_read2 else "/0"
    return f"{read.query_name or ''}{suffix}"


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

    # artifact_flags[cb] = {"tso": int, "polya": int, "noncanon": int, "tso_concat": int}
    artifact_flags: dict = field(default_factory=lambda: defaultdict(_dd_flags))
    alignment_sums: dict = field(default_factory=lambda: defaultdict(_dd_alignment))

    # length_samples[category] = list of read lengths (reservoir sample)
    length_samples: dict = field(default_factory=lambda: defaultdict(_dd_list))

    # length_bin_counts[category][bin_idx] = exact read count per length bin.
    # bin_idx 0–5 per _get_length_bin(); no sampling — every read is counted.
    length_bin_counts: dict = field(default_factory=lambda: defaultdict(_dd_int))

    # insert_size_signal / insert_size_noise: reservoir-sampled insert sizes
    # (abs(template_length)) for properly-paired Illumina reads, split by
    # signal (EXONIC_SENSE) vs noise (everything else).  Populated only when
    # paired_end_chimeric=True; empty list otherwise.
    insert_size_signal: list = field(default_factory=_make_insert_size_reservoir)
    insert_size_noise:  list = field(default_factory=_make_insert_size_reservoir)

    # intergenic_reads: lightweight side-table for intergenic profiler
    # Each entry additionally carries aligned bases, UMI, length and read key
    # so every downstream aggregate can be reclassified consistently.
    intergenic_reads: list = field(default_factory=_make_intergenic_reservoir)

    # exonic_sense_three_prime: reservoir-sampled 3′ end positions of
    # exonic-sense reads, used for the polyA-site-anchored fraction
    exonic_sense_three_prime: list = field(default_factory=_make_position_reservoir)

    # exonic_sense_five_prime: reservoir-sampled 5′ start positions of
    # exonic-sense reads, used for the TSS/CAGE-anchored fraction
    exonic_sense_five_prime: list = field(default_factory=_make_position_reservoir)
    # Unified per-read endpoints: (contig, strand, five_prime, three_prime, read_key)
    exonic_sense_endpoints: list = field(default_factory=_make_endpoint_reservoir)

    # Optional read_key -> (category, cell barcode), populated for compare.
    read_assignments: dict = field(default_factory=dict)

    # Total reads processed (for progress / logging)
    n_reads_processed: int = 0
    n_reads_skipped:   int = 0
    n_reads_skipped_not_called_cell: int = 0
    n_primary_mapped: int = 0
    n_secondary: int = 0
    n_supplementary: int = 0
    n_qcfail: int = 0
    n_duplicate: int = 0


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
    alignment_sums: dict = field(default_factory=lambda: defaultdict(_dd_alignment))
    length_samples: dict = field(default_factory=lambda: defaultdict(_dd_list))

    # length_bin_counts[category][bin_idx] = exact read count per length bin
    length_bin_counts: dict = field(default_factory=lambda: defaultdict(_dd_int))

    insert_size_signal: list = field(default_factory=_make_insert_size_reservoir)
    insert_size_noise:  list = field(default_factory=_make_insert_size_reservoir)

    intergenic_reads: list = field(default_factory=_make_intergenic_reservoir)
    exonic_sense_three_prime: list = field(default_factory=_make_position_reservoir)
    exonic_sense_five_prime:  list = field(default_factory=_make_position_reservoir)
    exonic_sense_endpoints:    list = field(default_factory=_make_endpoint_reservoir)
    read_assignments: dict = field(default_factory=dict)
    umi_tracking_enabled: bool = True

    n_reads_total:     int = 0
    n_reads_processed: int = 0
    n_reads_skipped:   int = 0
    n_reads_skipped_not_called_cell: int = 0
    n_reads_mapped_index: int = 0
    n_reads_unmapped: int = 0
    n_primary_mapped: int = 0
    n_secondary: int = 0
    n_supplementary: int = 0
    n_qcfail: int = 0
    n_duplicate: int = 0


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
    index               = args.get("index")
    if index is None:
        index = _WORKER_INDEX
    if index is None:
        raise RuntimeError(
            "no annotation index available in worker: pass it in the task dict "
            "or initialise the pool with _init_worker()"
        )
    reference_path      = args.get("reference_path")
    cell_barcodes       = args.get("cell_barcodes")  # set or None
    worker_seed         = args.get("seed")
    tso_sequences       = args.get("tso_sequences")
    tso_min_match       = args.get("tso_min_match", TSO_MIN_MATCH_LENGTH)
    tso_check_polyg     = args.get("tso_check_polyg", True)
    store_read_assignments = args.get("store_read_assignments", False)

    if worker_seed is not None:
        random.seed(worker_seed)

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
        tso_sequences=tso_sequences,
        tso_min_match=tso_min_match,
        tso_check_polyg=tso_check_polyg,
        reference=reference,
    )

    result = ContigResult(contig=contig)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # Handle mitochondrial contig name variants
        fetch_contig = contig
        try:
            for read in bam.fetch(fetch_contig):
                result.n_reads_processed += 1

                if read.flag & 0x200:
                    result.n_qcfail += 1
                if read.flag & 0x400:
                    result.n_duplicate += 1
                if read.is_secondary:
                    result.n_secondary += 1
                elif read.is_supplementary:
                    result.n_supplementary += 1
                elif not read.is_unmapped:
                    result.n_primary_mapped += 1

                res = classifier.classify(read)
                if res is None:
                    result.n_reads_skipped += 1
                    continue

                cb  = res.cell_barcode
                cat = res.category
                read_key = _read_key(read)

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
                if store_read_assignments:
                    result.read_assignments[read_key] = (cat, cb)

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
                if res.is_tso_concatemer:
                    result.artifact_flags[cb]["tso_concat"] += 1
                if res.is_discordant_pair:
                    result.artifact_flags[cb]["discordant_pair"] += 1
                if res.mapq < LOW_MAPQ_THRESHOLD:
                    result.artifact_flags[cb]["low_mapq"] += 1

                alignment = result.alignment_sums[cb]
                alignment["mapq_sum"] += res.mapq
                if res.edit_distance is not None:
                    alignment["nm_sum"] += res.edit_distance
                    alignment["nm_observed"] += 1
                alignment["softclip_bases"] += res.softclip_bases
                alignment["query_bases"] += res.query_bases

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

                # Store intergenic read position for the profiler.
                # Use read.reference_end (true rightmost genomic coordinate,
                # inclusive of intron spans) rather than res.pos + res.read_length,
                # and pick the strand-correct 3′ end — same reasoning as the
                # EXONIC_SENSE block below.
                if cat in (ReadCategory.INTERGENIC_SPARSE,
                           ReadCategory.INTERGENIC_REPEAT):
                    ref_end = read.reference_end or (res.pos + res.read_length)
                    three_prime = res.pos if res.is_reverse else ref_end
                    _reservoir_add(result.intergenic_reads, (
                        res.contig,
                        res.pos,
                        ref_end,
                        "-" if res.is_reverse else "+",
                        cb,
                        res.has_junction,
                        three_prime,
                        res.aligned_reference_bases,
                        res.umi,
                        res.read_length,
                        read_key,
                    ), 500_000)

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
                    _reservoir_add(
                        result.exonic_sense_endpoints,
                        (
                            res.contig,
                            "-" if res.is_reverse else "+",
                            five_prime,
                            three_prime,
                            read_key,
                        ),
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
    seed: Optional[int] = None,
    tso_sequences: Optional[list] = None,
    tso_min_match: int = TSO_MIN_MATCH_LENGTH,
    tso_check_polyg: bool = True,
    store_read_assignments: bool = False,
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
        Explicit list of contigs to process. Default: every reference with at
        least one mapped record according to the BAM index.
    cell_barcodes:
        Set of called-cell barcode strings (after stripping the -1 suffix).
        If provided, reads whose CB tag (also -1-stripped) is not in this set
        are skipped and counted in n_reads_skipped_not_called_cell.
    """
    bam_path = str(bam_path)

    mapped_index = 0
    unmapped_index = 0
    mapped_by_contig: dict[str, int] = {}
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            mapped_by_contig = {
                s.contig: int(s.mapped) for s in bam.get_index_statistics()
            }
            mapped_index = int(bam.mapped)
            unmapped_index = int(bam.unmapped)
    except Exception:
        pass
    if contigs is None:
        # Process every reference that actually has mapped records.  The old
        # 1-Mb length filter silently dropped spike-ins, viral contigs, decoys
        # and short scaffolds from both totals and artifact composition.
        try:
            contigs = [c for c, n in mapped_by_contig.items() if n > 0]
            # Longest-job-first, where the job is reads to classify rather
            # than reference span: a short contig can carry more alignments
            # than a long one, and the pool hands out work in this order.
            contigs.sort(key=lambda c: mapped_by_contig.get(c, 0), reverse=True)
        except Exception:
            contigs = list(meta.reference_names)

    logger.info("Processing %d contigs with %d threads …", len(contigs), threads)

    # If a seed is provided, seed the parent process's RNG (used during merges)
    # and derive deterministic per-worker seeds keyed on contig name so that
    # parallel workers stay independent but reproducible.
    if seed is not None:
        random.seed(seed)

    def _worker_seed(contig: str) -> Optional[int]:
        if seed is None:
            return None
        return (seed + zlib.adler32(contig.encode("utf-8"))) & 0xFFFFFFFF

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
            "reference_path":      reference_path,
            "seed":                _worker_seed(contig),
            "tso_sequences":       tso_sequences,
            "tso_min_match":       tso_min_match,
            "tso_check_polyg":     tso_check_polyg,
            "store_read_assignments": store_read_assignments,
        }
        for contig in contigs
    ]

    if not worker_args:
        logger.warning("No mapped contigs were found in the BAM index.")

    # Run workers
    contig_results = []
    if threads == 1:
        # Single-threaded path (easier to debug)
        _init_worker(index)
        for args in worker_args:
            contig_results.append(_contig_worker(args))
    else:
        with Pool(
            processes=threads, initializer=_init_worker, initargs=(index,)
        ) as pool:
            contig_results = pool.map(_contig_worker, worker_args)

    # Merge
    sample = SampleResult(bam_path=Path(bam_path), meta=meta)
    sample.n_reads_mapped_index = mapped_index
    sample.n_reads_unmapped = unmapped_index
    sample.n_reads_total = mapped_index + unmapped_index
    sample.umi_tracking_enabled = store_umis
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
    base.n_reads_mapped_index += other.n_reads_mapped_index
    base.n_reads_unmapped += other.n_reads_unmapped
    base.n_primary_mapped += other.n_primary_mapped
    base.n_secondary += other.n_secondary
    base.n_supplementary += other.n_supplementary
    base.n_qcfail += other.n_qcfail
    base.n_duplicate += other.n_duplicate
    base.umi_tracking_enabled = (
        getattr(base, "umi_tracking_enabled", True)
        and getattr(other, "umi_tracking_enabled", True)
    )

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

    for cb, values in other.alignment_sums.items():
        for key, amount in values.items():
            base.alignment_sums[cb][key] += amount

    for cat, lengths in other.length_samples.items():
        _merge_reservoir(base.length_samples[cat], lengths, MAX_LENGTH_SAMPLE)

    for cat, bin_counts in other.length_bin_counts.items():
        for bin_idx, n in bin_counts.items():
            base.length_bin_counts[cat][bin_idx] += n

    _merge_reservoir(base.intergenic_reads, other.intergenic_reads, 500_000)

    _merge_reservoir(base.insert_size_signal, other.insert_size_signal, 500_000)
    _merge_reservoir(base.insert_size_noise, other.insert_size_noise, 500_000)

    _merge_reservoir(base.exonic_sense_three_prime, other.exonic_sense_three_prime, MAX_LENGTH_SAMPLE)
    _merge_reservoir(base.exonic_sense_five_prime, other.exonic_sense_five_prime, MAX_LENGTH_SAMPLE)
    _merge_reservoir(base.exonic_sense_endpoints, other.exonic_sense_endpoints, MAX_LENGTH_SAMPLE)
    base.read_assignments.update(other.read_assignments)

    # Propagate polyA / TSS site dicts if base doesn't have them yet
    if not getattr(base, "_polya_site_dict", None) and getattr(other, "_polya_site_dict", None):
        base._polya_site_dict = other._polya_site_dict
    if not getattr(base, "_tss_site_dict", None) and getattr(other, "_tss_site_dict", None):
        base._tss_site_dict = other._tss_site_dict


def relabel_barcode(result: SampleResult, old: str, new: str) -> None:
    """Relabel a synthetic barcode across every barcode-indexed structure.

    Plate wells without CB tags are classified under ``NO_BARCODE`` inside
    their individual BAM.  Before plate aggregation that sentinel must become
    the well identity; otherwise every well collapses into one pseudo-cell.
    """
    if old == new:
        return
    for attr in (
        "read_counts", "base_counts", "umi_sets", "artifact_flags",
        "alignment_sums",
    ):
        mapping = getattr(result, attr)
        if old not in mapping:
            continue
        value = mapping.pop(old)
        if new in mapping:
            if attr == "umi_sets":
                for cat, members in value.items():
                    mapping[new][cat].update(members)
            else:
                for key, amount in value.items():
                    mapping[new][key] += amount
        else:
            mapping[new] = value

    for idx, rec in enumerate(result.intergenic_reads):
        if len(rec) >= 5 and rec[4] == old:
            mutable = list(rec)
            mutable[4] = new
            result.intergenic_reads[idx] = tuple(mutable)

    for read_key, (category, barcode) in list(result.read_assignments.items()):
        if barcode == old:
            result.read_assignments[read_key] = (category, new)


# ---------------------------------------------------------------------------
# Internal merge helper
# ---------------------------------------------------------------------------

def _merge_contig(sample: SampleResult, cr: ContigResult) -> None:
    sample.n_reads_processed += cr.n_reads_processed
    sample.n_reads_skipped   += cr.n_reads_skipped
    sample.n_reads_skipped_not_called_cell += cr.n_reads_skipped_not_called_cell
    sample.n_primary_mapped += cr.n_primary_mapped
    sample.n_secondary += cr.n_secondary
    sample.n_supplementary += cr.n_supplementary
    sample.n_qcfail += cr.n_qcfail
    sample.n_duplicate += cr.n_duplicate

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

    for cb, values in cr.alignment_sums.items():
        for key, amount in values.items():
            sample.alignment_sums[cb][key] += amount

    for cat, lengths in cr.length_samples.items():
        _merge_reservoir(sample.length_samples[cat], lengths, MAX_LENGTH_SAMPLE)

    for cat, bin_counts in cr.length_bin_counts.items():
        for bin_idx, n in bin_counts.items():
            sample.length_bin_counts[cat][bin_idx] += n

    # Merge intergenic read records (reservoir-sampled to cap memory)
    _merge_reservoir(sample.intergenic_reads, cr.intergenic_reads, 500_000)

    # Merge insert size samples
    _merge_reservoir(sample.insert_size_signal, cr.insert_size_signal, 500_000)
    _merge_reservoir(sample.insert_size_noise, cr.insert_size_noise, 500_000)

    # Merge exonic-sense 3′ positions
    _merge_reservoir(sample.exonic_sense_three_prime, cr.exonic_sense_three_prime, MAX_LENGTH_SAMPLE)
    _merge_reservoir(sample.exonic_sense_five_prime, cr.exonic_sense_five_prime, MAX_LENGTH_SAMPLE)
    _merge_reservoir(sample.exonic_sense_endpoints, cr.exonic_sense_endpoints, MAX_LENGTH_SAMPLE)
    sample.read_assignments.update(cr.read_assignments)

    return


# ---------------------------------------------------------------------------
# Reservoir sampling utility
# ---------------------------------------------------------------------------

def _reservoir_add(reservoir: list, value, max_size: int) -> None:
    """
    Add *value* to *reservoir* using reservoir sampling (Algorithm R).
    Keeps the reservoir size ≤ max_size with uniform random sampling.

    For correct Algorithm R, the running count of items offered is required.
    A `_Reservoir` carries that count; plain lists fall back to a slightly
    biased length-based heuristic (kept for backwards compatibility).
    """
    if isinstance(reservoir, _Reservoir):
        reservoir.n_seen += 1
        n = reservoir.n_seen
    else:
        n = len(reservoir) + 1

    if len(reservoir) < max_size:
        reservoir.append(value)
    else:
        idx = random.randrange(n)
        if idx < max_size:
            reservoir[idx] = value


def _merge_reservoir(base: list, other: list, max_size: int) -> None:
    """Merge two independent uniform reservoirs without contig-size bias.

    Simply offering the sampled items from each contig to Algorithm R treats a
    10-read contig and a 10-million-read contig as equal once both reservoirs
    are capped.  We instead draw how many retained items should come from each
    source using the corresponding population sizes, then sample uniformly
    within each source reservoir.
    """
    if not other:
        return
    n_a = getattr(base, "n_seen", len(base)) or len(base)
    n_b = getattr(other, "n_seen", len(other)) or len(other)
    total = n_a + n_b
    keep = min(max_size, total)

    remaining_a, remaining_b = n_a, n_b
    keep_a = 0
    # Hypergeometric draw implemented with the seeded stdlib RNG so --seed
    # controls workers and merges consistently.
    for _ in range(keep):
        remaining_total = remaining_a + remaining_b
        if remaining_total <= 0:
            break
        if random.random() < (remaining_a / remaining_total):
            keep_a += 1
            remaining_a -= 1
        else:
            remaining_b -= 1
    keep_b = keep - keep_a

    chosen_a = random.sample(list(base), keep_a) if keep_a else []
    chosen_b = random.sample(list(other), keep_b) if keep_b else []
    base[:] = chosen_a + chosen_b
    random.shuffle(base)
    if isinstance(base, _Reservoir):
        base.n_seen = total
