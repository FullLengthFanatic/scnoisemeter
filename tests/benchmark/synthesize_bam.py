#!/usr/bin/env python3
"""
synthesize_bam.py — controlled ground-truth BAMs for the scNoiseMeter benchmark.

Emits two BAMs and their ground-truth tables:

  synthetic.bam           — Exp 1 (classifier accuracy), Exp 2 (noise dose-
                            response), Exp 4 (TSO / polyA artifact-flag cells).
                            Use --platform ont to score.

  synthetic_illumina.bam  — Exp 5: paired-end CHIMERIC validation. Use
                            --platform illumina_10x to exercise the
                            paired-end chimeric branch (FLAG + RNEXT/PNEXT
                            + TLEN); long-read SA-tag chimeric is already
                            covered by the CHIMERIC cell in synthetic.bam.

Exp 1 (per-read classifier accuracy, confusion matrix):
  - "Simple" categories: one CB per true category, N reads per CB.
  - "Clustered" intergenic (HOTSPOT, NOVEL, REPEAT): multiple CBs per
    category because the intergenic profiler requires >= 3 distinct
    barcodes per locus to declare significance.

Exp 2 (per-cell noise-fraction fidelity, dose-response):
  Five CBs with known mixtures of signal (EXONIC_SENSE) and noise drawn
  from NOISE_COMPOSITION.

Exp 4 (per-cell artifact-flag counters, TSO / polyA):
  Two CBs whose reads are EXONIC_SENSE but carry a TSO 5' soft-clip or
  sit at a position with a downstream A-run. Validated against per-cell
  n_tso / n_polya in cell_metrics.tsv.

Categories covered (16 of 19):

  Simple (1 CB each):
    EXONIC_SENSE, EXONIC_ANTISENSE,
    INTRONIC_PURE, INTRONIC_BOUNDARY, INTRONIC_JXNSPAN,
    INTERGENIC_SPARSE,
    MITOCHONDRIAL, MULTIMAPPER, CHIMERIC,
    AMBIGUOUS, AMBIGUOUS_COD_COD, AMBIGUOUS_COD_NCOD,
    UNASSIGNED

  Clustered (3 CBs each, multi-locus):
    INTERGENIC_HOTSPOT, INTERGENIC_NOVEL, INTERGENIC_REPEAT

Artifact-flag cells (per-cell counter validation, layered on EXONIC_SENSE):
  TSO_INVASION (n_tso), POLYA_PRIMING (n_polya)

Side outputs:
  whitelist.txt — every legitimate CB. The runner must pass this with
  --barcode-whitelist; the UNASSIGNED CB is deliberately excluded so
  its reads trip the not-on-whitelist branch in the classifier.

Not yet covered (requires upstream feature work):
  NUMT-derived MITOCHONDRIAL sub-classification (per-read NUMT
  disambiguation is not yet implemented in scnoisemeter).

Key correctness invariants:
  * Strand-safe + exons: + strand exons with no base overlap to any - strand
    exon. Ensures EXONIC_ANTISENSE reads do not become EXONIC_SENSE vs an
    overlapping - strand gene.
  * True introns: inside a + strand gene body and outside ANY annotated exon.
  * Clean intergenic: >= 10 kb from any gene and not inside any exon.
  * JXNSPAN: read's mapped bases (both M segments) lie in intronic space, with
    one CIGAR N op — triggers the classifier's INTRONIC_PURE -> INTRONIC_JXNSPAN
    promotion regardless of canonicality.
  * AMBIGUOUS_COD_*: read is placed in a gene-body overlap region for a
    coding+coding or coding+noncoding gene pair on the same strand.
  * HOTSPOT: mono-exonic clustered reads at positions where the reference
    contains a 6+ A-run within 20 bp downstream of the read 3' end.
  * NOVEL: spliced reads with a non-canonical donor dinucleotide (NOT GT/GC/AT)
    at the junction, clustered at one locus across >= 3 CBs. The classifier
    sets has_noncanonical_junction=True, which the intergenic profiler
    interprets as splice evidence when ranking NOVEL candidates.
  * REPEAT: clustered intergenic reads whose positions are covered by the
    synthesized repeats.bed fixture emitted next to the BAM.
"""

from __future__ import annotations

import argparse
import gzip
import os
import random
import re
import sys
from bisect import bisect_right
from pathlib import Path

import pysam


WORK_CHROM = "chr22"
MITO_CHROM = "chrM"
READ_LEN = 100
N_PER_CATEGORY_EXP1 = 200
MIXTURE_CELL_N_READS = 1000
NOISE_FRACTIONS = [0.00, 0.10, 0.25, 0.40, 0.55]

NOISE_COMPOSITION = {
    "INTRONIC_PURE":     0.40,
    "INTERGENIC_SPARSE": 0.30,
    "EXONIC_ANTISENSE":  0.20,
    "INTRONIC_BOUNDARY": 0.10,
}

EXP1_SIMPLE_CATEGORIES = [
    "EXONIC_SENSE",
    "EXONIC_ANTISENSE",
    "INTRONIC_PURE",
    "INTRONIC_BOUNDARY",
    "INTRONIC_JXNSPAN",
    "INTERGENIC_SPARSE",
    "MITOCHONDRIAL",
    "MULTIMAPPER",
    "CHIMERIC",
    "AMBIGUOUS",
    "AMBIGUOUS_COD_COD",
    "AMBIGUOUS_COD_NCOD",
    "UNASSIGNED",
]

# Artifact-flag cells (per-cell n_tso / n_polya counter validation).
# These do NOT change the read-level ReadCategory; the underlying read
# is EXONIC_SENSE, the TSO/polyA flag is layered on top.
EXP4_FLAG_CATEGORIES = [
    "TSO_INVASION",
    "POLYA_PRIMING",
]

# Whitelist-failed CB used by UNASSIGNED. Excluded from whitelist.txt.
UNASSIGNED_CB = "CELL_EXP1_UNASSIGNED"

# 5' soft-clip seed for TSO_INVASION reads (must be a >=12bp substring
# of a known TSO; we use the canonical 10x TSO prefix).
TSO_SOFTCLIP_SEED = "AAGCAGTGGTAT"   # TSO_10X[:12]

EXP1_CLUSTERED_CATEGORIES = [
    "INTERGENIC_HOTSPOT",
    "INTERGENIC_NOVEL",
    "INTERGENIC_REPEAT",
]

CLUSTERED_N_CBS = 3
CLUSTERED_N_LOCI = 10
CLUSTERED_READS_PER_LOCUS_PER_CB = 7   # 3 * 7 = 21 reads/locus; ADAPTIVE_MIN_READS=5
CLUSTERED_JITTER = 100                  # bp, << INTERGENIC_LOCUS_WINDOW=500

JXNSPAN_SKIP = 300                      # bp for the CIGAR N gap
NOVEL_SKIP = 500                        # bp for the CIGAR N gap in NOVEL synthesis
NOVEL_FOOTPRINT = NOVEL_SKIP + READ_LEN  # span of a NOVEL read on the reference

CANONICAL_DONORS = {"GT", "GC", "AT"}

_NONCODING_BIOTYPES = {
    "lncRNA", "processed_transcript", "retained_intron",
    "nonsense_mediated_decay", "non_stop_decay", "TEC",
    "transcribed_unitary_pseudogene", "transcribed_processed_pseudogene",
    "transcribed_unprocessed_pseudogene", "processed_pseudogene",
    "unprocessed_pseudogene", "pseudogene", "unitary_pseudogene",
    "polymorphic_pseudogene", "misc_RNA", "miRNA", "snRNA", "snoRNA",
    "rRNA", "ribozyme", "sRNA", "scaRNA", "vault_RNA",
    "Mt_tRNA", "Mt_rRNA",
}

_CODING_BIOTYPES = {
    "protein_coding", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
    "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtf", required=True, help="GENCODE GTF (gzipped ok)")
    p.add_argument("--fasta", required=True, help="Reference FASTA with .fai (gzipped+bgzf ok)")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


# ---------------------------------------------------------------------------
# GTF parsing (now also captures biotype)
# ---------------------------------------------------------------------------

def parse_gtf(gtf_path, chrom):
    """Return (genes, exons) on `chrom`. Each gene is
    (start, end, strand, gene_id, gene_biotype); each exon is
    (start, end, strand, gene_id). Coordinates 0-based half-open.
    Missing biotype becomes the empty string."""
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    genes: list[tuple[int, int, str, str, str]] = []
    exons: list[tuple[int, int, str, str]] = []
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[0] != chrom:
                continue
            feat = cols[2]
            if feat not in ("gene", "exon"):
                continue
            start = int(cols[3]) - 1
            end = int(cols[4])
            strand = cols[6]
            attrs = cols[8]
            gene_id = None
            biotype = ""
            for kv in attrs.split(";"):
                kv = kv.strip()
                if not kv:
                    continue
                if kv.startswith("gene_id"):
                    parts = kv.split('"')
                    if len(parts) >= 2:
                        gene_id = parts[1]
                elif kv.startswith("gene_biotype") or kv.startswith("gene_type"):
                    parts = kv.split('"')
                    if len(parts) >= 2:
                        biotype = parts[1]
            if feat == "gene":
                genes.append((start, end, strand, gene_id, biotype))
            else:
                exons.append((start, end, strand, gene_id))
    return genes, exons


# ---------------------------------------------------------------------------
# Interval utilities
# ---------------------------------------------------------------------------

def _merge(intervals):
    intervals = sorted(intervals)
    out: list[list[int]] = []
    for s, e in intervals:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def _starts(merged):
    return [s for s, _ in merged]


def _overlaps(merged, starts, s, e):
    if not merged:
        return False
    i = bisect_right(starts, e - 1) - 1
    if i < 0:
        return False
    return merged[i][1] > s


def _fully_inside(merged, starts, s, e):
    if not merged:
        return False
    i = bisect_right(starts, s) - 1
    if i < 0:
        return False
    ms, me = merged[i]
    return ms <= s and e <= me


# ---------------------------------------------------------------------------
# Region pickers
# ---------------------------------------------------------------------------

def build_all_pickers(genes, exons, chrom_size, fa, buffer=10_000):
    """Return a dict of region / position lists for each synthetic category."""
    plus_exons_raw = [(s, e) for s, e, st, _ in exons if st == "+"]
    minus_exons_raw = [(s, e) for s, e, st, _ in exons if st == "-"]
    all_exons_raw = [(s, e) for s, e, _, _ in exons]

    plus_exon_union = _merge(plus_exons_raw)
    minus_exon_union = _merge(minus_exons_raw)
    all_exon_union = _merge(all_exons_raw)

    minus_starts = _starts(minus_exon_union)
    all_exon_starts = _starts(all_exon_union)

    # --- strand-safe + exons ---
    safe_plus_exons = []
    for s, e in plus_exon_union:
        if (e - s) < READ_LEN + 40:
            continue
        if _overlaps(minus_exon_union, minus_starts, s, e):
            continue
        safe_plus_exons.append((s, e))

    # --- + strand gene bodies and + strand introns ---
    plus_gene_union = _merge([(s, e) for s, e, st, *_ in genes if st == "+"])
    plus_gene_starts = _starts(plus_gene_union)

    plus_introns: list[tuple[int, int]] = []
    for gs, ge in plus_gene_union:
        i = bisect_right(all_exon_starts, ge - 1) - 1
        exons_in_gene: list[tuple[int, int]] = []
        while i >= 0:
            es, ee = all_exon_union[i]
            if ee <= gs:
                break
            if es < ge and ee > gs:
                exons_in_gene.append((max(es, gs), min(ee, ge)))
            i -= 1
        exons_in_gene.sort()
        cursor = gs
        for es, ee in exons_in_gene:
            if es > cursor:
                plus_introns.append((cursor, es))
            cursor = max(cursor, ee)
        if cursor < ge:
            plus_introns.append((cursor, ge))
    plus_introns = [(s, e) for s, e in plus_introns if (e - s) >= READ_LEN + 100]

    # --- introns usable for INTRONIC_JXNSPAN (50M{JXNSPAN_SKIP}N50M fits inside) ---
    jxnspan_footprint = 50 + JXNSPAN_SKIP + 50 + 40   # 40 bp slack at ends
    jxnspan_introns = [(s, e) for s, e in plus_introns if (e - s) >= jxnspan_footprint]

    # --- INTRONIC_BOUNDARY candidates (unchanged) ---
    boundary_positions: list[int] = []
    for s, e in safe_plus_exons:
        if e - s < 40:
            continue
        up_start = e - 30
        if up_start < s:
            continue
        if _overlaps(all_exon_union, all_exon_starts, e, e + 70):
            continue
        if not _fully_inside(plus_gene_union, plus_gene_starts, e, e + 70):
            continue
        boundary_positions.append(e)

    # --- AMBIGUOUS_COD_* overlap regions (strand-aware, biotype-aware) ---
    cod_cod_overlaps, cod_ncod_overlaps = _build_coding_overlap_regions(genes, exons)

    # --- Intergenic: >= buffer from any gene AND not inside any exon ---
    gene_union = _merge([(s, e) for s, e, *_ in genes])
    buffered = [(max(0, s - buffer), min(chrom_size, e + buffer)) for s, e in gene_union]
    buffered_union = _merge(buffered)

    # Use a stride that leaves plenty of headroom for all consumers
    intergenic_positions: list[int] = []
    stride = 2000
    prev_end = 0
    for s, e in buffered_union:
        for p in range(prev_end, s - NOVEL_FOOTPRINT - 50, stride):
            # The widest footprint (NOVEL) is NOVEL_FOOTPRINT bp; any position
            # that is clean across that window is safe for every other category.
            if not _overlaps(all_exon_union, all_exon_starts, p, p + NOVEL_FOOTPRINT):
                intergenic_positions.append(p)
        prev_end = max(prev_end, e)
    for p in range(prev_end, chrom_size - NOVEL_FOOTPRINT - 50, stride):
        if not _overlaps(all_exon_union, all_exon_starts, p, p + NOVEL_FOOTPRINT):
            intergenic_positions.append(p)

    random.shuffle(intergenic_positions)

    # --- Sub-select HOTSPOT, NOVEL, REPEAT locus seeds from the intergenic pool ---
    # Each position is used ONCE across all consumers to avoid ground-truth
    # collisions between SPARSE and the clustered categories.
    hotspot_seeds: list[int] = []
    novel_seeds: list[int] = []
    repeat_seeds: list[int] = []
    leftover: list[int] = []

    # Precompile the A-run regex once
    a_run_re = re.compile("A{6,}")

    for p in intergenic_positions:
        if len(hotspot_seeds) < CLUSTERED_N_LOCI:
            three_prime = p + READ_LEN
            try:
                ctx = fa.fetch(WORK_CHROM, three_prime, three_prime + 20).upper()
            except (ValueError, KeyError):
                ctx = ""
            if a_run_re.search(ctx):
                hotspot_seeds.append(p)
                continue
        if len(novel_seeds) < CLUSTERED_N_LOCI:
            try:
                donor = fa.fetch(WORK_CHROM, p + 50, p + 52).upper()
            except (ValueError, KeyError):
                donor = ""
            if donor and donor not in CANONICAL_DONORS:
                novel_seeds.append(p)
                continue
        if len(repeat_seeds) < CLUSTERED_N_LOCI:
            repeat_seeds.append(p)
            continue
        leftover.append(p)

    if len(hotspot_seeds) < CLUSTERED_N_LOCI:
        print(f"WARNING: only {len(hotspot_seeds)} HOTSPOT loci found "
              f"(wanted {CLUSTERED_N_LOCI})", file=sys.stderr)
    if len(novel_seeds) < CLUSTERED_N_LOCI:
        print(f"WARNING: only {len(novel_seeds)} NOVEL loci found "
              f"(wanted {CLUSTERED_N_LOCI})", file=sys.stderr)
    if len(repeat_seeds) < CLUSTERED_N_LOCI:
        print(f"WARNING: only {len(repeat_seeds)} REPEAT loci found "
              f"(wanted {CLUSTERED_N_LOCI})", file=sys.stderr)

    # --- POLYA_PRIMING seeds: positions inside safe + exons whose 20 bp
    # downstream of the read 3' end carries an A-run >= 6, matching the
    # classifier's POLYA_RUN_MIN_LENGTH / POLYA_CONTEXT_WINDOW thresholds.
    polya_positions: list[int] = []
    for s, e in safe_plus_exons:
        for try_pos in range(s + 10, e - READ_LEN - 10, 50):
            try:
                ctx = fa.fetch(WORK_CHROM, try_pos + READ_LEN,
                               try_pos + READ_LEN + 20).upper()
            except (ValueError, KeyError):
                continue
            if a_run_re.search(ctx):
                polya_positions.append(try_pos)
                if len(polya_positions) >= N_PER_CATEGORY_EXP1:
                    break
        if len(polya_positions) >= N_PER_CATEGORY_EXP1:
            break
    if len(polya_positions) < N_PER_CATEGORY_EXP1:
        print(f"WARNING: only {len(polya_positions)} POLYA_PRIMING positions found "
              f"(wanted {N_PER_CATEGORY_EXP1}); reads will recycle positions",
              file=sys.stderr)

    return {
        "safe_plus_exons": safe_plus_exons,
        "plus_introns": plus_introns,
        "jxnspan_introns": jxnspan_introns,
        "boundary_positions": boundary_positions,
        "cod_cod_overlaps": cod_cod_overlaps,
        "cod_ncod_overlaps": cod_ncod_overlaps,
        "intergenic_positions": leftover,
        "hotspot_seeds": hotspot_seeds,
        "novel_seeds": novel_seeds,
        "repeat_seeds": repeat_seeds,
        "polya_positions": polya_positions,
        "all_exon_union": all_exon_union,
    }


def _classify_overlap_type(a: str, b: str) -> str:
    a_cod = a in _CODING_BIOTYPES
    b_cod = b in _CODING_BIOTYPES
    if a_cod and b_cod:
        return "cod_cod"
    if a_cod or b_cod:
        return "cod_ncod"
    return "ncod_ncod"


def _build_coding_overlap_regions(genes, exons):
    """
    Walk same-strand exon pairs across DIFFERENT genes and emit exon-level
    overlap intervals classified by biotype. scnoisemeter computes the shared
    ambiguous regions on exons (not gene bodies) — see the note in
    annotation.py _unique_and_shared: intronic overlaps between gene bodies are
    not ambiguous for noise quantification.
    """
    biotype_by_gene: dict[str, str] = {gid: bt for _, _, _, gid, bt in genes}
    cod_cod: list[tuple[int, int, str]] = []
    cod_ncod: list[tuple[int, int, str]] = []

    for strand in ("+", "-"):
        same_strand = [
            (s, e, gid) for s, e, st, gid in exons
            if st == strand and gid is not None
        ]
        same_strand.sort(key=lambda x: x[0])
        for i, (si, ei, gi) in enumerate(same_strand):
            bi = biotype_by_gene.get(gi, "")
            for sj, ej, gj in same_strand[i + 1 :]:
                if sj >= ei:
                    break
                if gi == gj:
                    continue
                ov_start = max(si, sj)
                ov_end = min(ei, ej)
                if ov_end - ov_start < READ_LEN + 40:
                    continue
                bj = biotype_by_gene.get(gj, "")
                kind = _classify_overlap_type(bi, bj)
                if kind == "cod_cod":
                    cod_cod.append((ov_start, ov_end, strand))
                else:
                    # cod_ncod and ncod_ncod both land in AMBIGUOUS_COD_NCOD
                    # per the classifier's shared-set definition.
                    cod_ncod.append((ov_start, ov_end, strand))

    return cod_cod, cod_ncod


# ---------------------------------------------------------------------------
# Read construction
# ---------------------------------------------------------------------------

def make_read(
    name: str,
    ref_id: int,
    pos: int,
    strand: str,
    cb: str,
    umi: str,
    cigar: str = "100M",
    nh: int = 1,
    sa: str | None = None,
    mapq: int = 60,
    seq: str | None = None,
) -> pysam.AlignedSegment:
    a = pysam.AlignedSegment()
    a.query_name = name
    a.reference_id = ref_id
    a.reference_start = pos
    a.mapping_quality = mapq
    a.flag = 16 if strand == "-" else 0
    a.cigarstring = cigar
    # query length = sum of M/I/S/=/X lengths in CIGAR; we always emit Ms only
    m_len = _cigar_query_length(cigar)
    a.query_sequence = (seq if seq is not None else "N" * m_len)
    a.query_qualities = pysam.qualitystring_to_array("I" * m_len)
    tags: list[tuple[str, object]] = [("CB", cb), ("UB", umi), ("NH", nh)]
    if sa is not None:
        tags.append(("SA", sa))
    a.tags = tags
    return a


def _cigar_query_length(cigar: str) -> int:
    total = 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in ("M", "I", "S", "=", "X"):
            total += int(n)
    return total


def make_paired_read(
    name: str,
    ref_id_self: int,
    pos_self: int,
    strand_self: str,
    ref_id_mate: int,
    pos_mate: int,
    mate_unmapped: bool,
    is_read1: bool,
    tlen: int,
    cb: str,
    umi: str,
    cigar: str = "100M",
) -> pysam.AlignedSegment:
    """Build one mapped half of a paired-end alignment.

    We emit a single record per pair (always the mapped half). The
    paired-end chimeric check evaluates one record at a time by
    inspecting FLAG bits + RNEXT/PNEXT + TLEN, so a single mapped record
    is enough to trigger the branch.
    """
    a = pysam.AlignedSegment()
    a.query_name = name
    a.reference_id = ref_id_self
    a.reference_start = pos_self
    a.mapping_quality = 60

    flag = 1   # PAIRED
    if strand_self == "-":
        flag |= 16   # REVERSE
    if mate_unmapped:
        flag |= 8    # MATE_UNMAPPED
    flag |= 64 if is_read1 else 128
    a.flag = flag

    a.cigarstring = cigar
    qlen = _cigar_query_length(cigar)
    a.query_sequence = "N" * qlen
    a.query_qualities = pysam.qualitystring_to_array("I" * qlen)

    a.next_reference_id = ref_id_mate
    a.next_reference_start = pos_mate
    a.template_length = tlen
    a.tags = [("CB", cb), ("UB", umi), ("NH", 1)]
    return a


# ---------------------------------------------------------------------------
# Illumina paired-end CHIMERIC BAM (Exp 5)
# ---------------------------------------------------------------------------

ILLUMINA_N_PER_CATEGORY = 100
ILLUMINA_TLEN_LARGE = 2_000_000   # > ILLUMINA_CHIMERIC_INSERT_SIZE (1 Mb)
ILLUMINA_TLEN_NORMAL = 350         # representative library insert
ILLUMINA_CATEGORIES = [
    "ILLUMINA_CONCORDANT",          # negative control: should NOT be chimeric
    "ILLUMINA_CHIMERIC_INTERCHROM", # mate on a different chromosome
    "ILLUMINA_DISCORDANT",          # mate flagged unmapped
    "ILLUMINA_LARGETLEN",           # same chrom but |TLEN| > 1 Mb
]
# Chromosome used for the "other" mate in inter-chromosomal pairs.
ILLUMINA_OTHER_CHROM = "chr1"


def synthesize_illumina_bam(outdir, header, ref_id, picks, refs):
    """Emit synthetic_illumina.bam (Exp 5) and ground truth.

    Validates the paired-end chimeric branch reached when --platform
    illumina_10x is passed (no SA tag; classifier consults FLAG +
    RNEXT/PNEXT + TLEN). Also writes the rows for the same CBs into
    the main ground truth, and the same CBs to whitelist.txt -- caller
    handles the latter via cbs_seen aggregation.
    """
    if ILLUMINA_OTHER_CHROM not in ref_id:
        print(f"WARNING: {ILLUMINA_OTHER_CHROM} not in reference; "
              f"skipping Illumina paired-end BAM", file=sys.stderr)
        return [], []

    reads: list[pysam.AlignedSegment] = []
    ground_truth: list[tuple[str, str, str, str]] = []
    umi_counter = 0

    def umi() -> str:
        nonlocal umi_counter
        u = f"UMI{umi_counter:08d}"
        umi_counter += 1
        return u

    for cat in ILLUMINA_CATEGORIES:
        cb = f"CELL_EXP5_{cat}"
        for i in range(ILLUMINA_N_PER_CATEGORY):
            s, e = random.choice(picks["safe_plus_exons"])
            pos = random.randint(s + 10, e - READ_LEN - 10)
            name = f"{cat}_{cb}_{i:06d}"

            if cat == "ILLUMINA_CONCORDANT":
                # Both mates on chr22, normal insert, mate mapped → not chimeric.
                # The mate would sit ~ILLUMINA_TLEN_NORMAL bp downstream.
                r = make_paired_read(
                    name, ref_id[WORK_CHROM], pos, "+",
                    ref_id[WORK_CHROM], pos + ILLUMINA_TLEN_NORMAL - READ_LEN,
                    mate_unmapped=False, is_read1=True,
                    tlen=ILLUMINA_TLEN_NORMAL, cb=cb, umi=umi(),
                )
            elif cat == "ILLUMINA_CHIMERIC_INTERCHROM":
                r = make_paired_read(
                    name, ref_id[WORK_CHROM], pos, "+",
                    ref_id[ILLUMINA_OTHER_CHROM], 1_000_000 + i * 1000,
                    mate_unmapped=False, is_read1=True,
                    tlen=0, cb=cb, umi=umi(),
                )
            elif cat == "ILLUMINA_DISCORDANT":
                r = make_paired_read(
                    name, ref_id[WORK_CHROM], pos, "+",
                    ref_id[WORK_CHROM], pos,   # SAM convention: unmapped mate
                    mate_unmapped=True, is_read1=True,
                    tlen=0, cb=cb, umi=umi(),
                )
            elif cat == "ILLUMINA_LARGETLEN":
                # Far mate on same chrom (no inter-chrom signal); TLEN > 1 Mb.
                mate_pos = pos + ILLUMINA_TLEN_LARGE
                r = make_paired_read(
                    name, ref_id[WORK_CHROM], pos, "+",
                    ref_id[WORK_CHROM], mate_pos,
                    mate_unmapped=False, is_read1=True,
                    tlen=ILLUMINA_TLEN_LARGE, cb=cb, umi=umi(),
                )
            else:
                raise ValueError(f"unknown Illumina category: {cat}")

            reads.append(r)
            true_cat = "EXONIC_SENSE" if cat == "ILLUMINA_CONCORDANT" else "CHIMERIC"
            ground_truth.append((name, cb, true_cat, f"exp5;variant={cat.lower()}"))

    reads.sort(key=lambda r: (r.reference_id, r.reference_start))
    bam_path = outdir / "synthetic_illumina.bam"
    tmp_path = outdir / "synthetic_illumina.unsorted.bam"
    with pysam.AlignmentFile(str(tmp_path), "wb", header=header) as bam_out:
        for r in reads:
            bam_out.write(r)
    pysam.sort("-o", str(bam_path), str(tmp_path))
    pysam.index(str(bam_path))
    os.remove(tmp_path)

    print(f"      -> {bam_path} ({len(reads)} reads, "
          f"{len(ILLUMINA_CATEGORIES)} cells)", file=sys.stderr)
    return reads, ground_truth


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    random.seed(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[1/8] Parsing GTF for {WORK_CHROM}...", file=sys.stderr)
    genes, exons = parse_gtf(args.gtf, WORK_CHROM)
    print(f"      {len(genes)} genes, {len(exons)} exon entries", file=sys.stderr)

    print(f"[2/8] Opening FASTA...", file=sys.stderr)
    fa = pysam.FastaFile(args.fasta)
    refs = list(fa.references)
    lens = [fa.get_reference_length(r) for r in refs]
    chrom_size = lens[refs.index(WORK_CHROM)]
    mito_size = lens[refs.index(MITO_CHROM)]

    print(f"      Building picker indices for {WORK_CHROM}...", file=sys.stderr)
    picks = build_all_pickers(genes, exons, chrom_size, fa)
    print(f"      safe_plus_exons={len(picks['safe_plus_exons'])}, "
          f"plus_introns={len(picks['plus_introns'])}, "
          f"jxnspan_introns={len(picks['jxnspan_introns'])}, "
          f"boundary_positions={len(picks['boundary_positions'])}, "
          f"cod_cod={len(picks['cod_cod_overlaps'])}, "
          f"cod_ncod={len(picks['cod_ncod_overlaps'])}, "
          f"hotspot_seeds={len(picks['hotspot_seeds'])}, "
          f"novel_seeds={len(picks['novel_seeds'])}, "
          f"repeat_seeds={len(picks['repeat_seeds'])}, "
          f"polya_positions={len(picks['polya_positions'])}, "
          f"intergenic_leftover={len(picks['intergenic_positions'])}", file=sys.stderr)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": r, "LN": l} for r, l in zip(refs, lens)],
        "PG": [{"ID": "synthesize_bam", "PN": "synthesize_bam", "VN": "0.2"}],
    }
    ref_id = {r: i for i, r in enumerate(refs)}

    reads: list[pysam.AlignedSegment] = []
    ground_truth: list[tuple[str, str, str, str]] = []

    umi_counter = 0
    def next_umi() -> str:
        nonlocal umi_counter
        u = f"UMI{umi_counter:08d}"
        umi_counter += 1
        return u

    inter_cursor = 0
    def next_inter_pos() -> int:
        nonlocal inter_cursor
        if inter_cursor >= len(picks["intergenic_positions"]):
            raise RuntimeError("exhausted intergenic positions; reduce noise or use more chromosomes")
        p = picks["intergenic_positions"][inter_cursor]
        inter_cursor += 1
        return p

    def emit_simple(cat: str, cb: str, read_idx: int, notes: str = ""):
        """Emit one read for a 'simple' category (not clustered intergenic)."""
        name = f"{cat}_{cb}_{read_idx:06d}"
        u = next_umi()
        if cat == "EXONIC_SENSE":
            s, e = random.choice(picks["safe_plus_exons"])
            pos = random.randint(s + 10, e - READ_LEN - 10)
            r = make_read(name, ref_id[WORK_CHROM], pos, "+", cb, u)
        elif cat == "EXONIC_ANTISENSE":
            s, e = random.choice(picks["safe_plus_exons"])
            pos = random.randint(s + 10, e - READ_LEN - 10)
            r = make_read(name, ref_id[WORK_CHROM], pos, "-", cb, u)
        elif cat == "INTRONIC_PURE":
            s, e = random.choice(picks["plus_introns"])
            pos = random.randint(s + 20, e - READ_LEN - 20)
            r = make_read(name, ref_id[WORK_CHROM], pos, "+", cb, u)
        elif cat == "INTRONIC_BOUNDARY":
            exon_end = random.choice(picks["boundary_positions"])
            pos = exon_end - 30
            r = make_read(name, ref_id[WORK_CHROM], pos, "+", cb, u)
        elif cat == "INTRONIC_JXNSPAN":
            s, e = random.choice(picks["jxnspan_introns"])
            pos = random.randint(s + 20, e - (50 + JXNSPAN_SKIP + 50) - 20)
            r = make_read(name, ref_id[WORK_CHROM], pos, "+", cb, u,
                          cigar=f"50M{JXNSPAN_SKIP}N50M")
        elif cat == "INTERGENIC_SPARSE":
            pos = next_inter_pos()
            r = make_read(name, ref_id[WORK_CHROM], pos, "+", cb, u)
        elif cat == "MITOCHONDRIAL":
            pos = random.randint(100, mito_size - READ_LEN - 100)
            r = make_read(name, ref_id[MITO_CHROM], pos, "+", cb, u)
        elif cat == "MULTIMAPPER":
            s, e = random.choice(picks["safe_plus_exons"])
            pos = random.randint(s + 10, e - READ_LEN - 10)
            r = make_read(name, ref_id[WORK_CHROM], pos, "+", cb, u, nh=2)
        elif cat == "CHIMERIC":
            s, e = random.choice(picks["safe_plus_exons"])
            pos = random.randint(s + 10, e - READ_LEN - 10)
            other = random.randint(1_000_000, 50_000_000)
            sa = f"chr1,{other},+,100M,60,0;"
            r = make_read(name, ref_id[WORK_CHROM], pos, "+", cb, u, sa=sa)
        elif cat == "AMBIGUOUS_COD_COD":
            pool = picks["cod_cod_overlaps"]
            if not pool:
                raise RuntimeError("no coding/coding gene-body overlaps found on this chromosome")
            s, e, st = random.choice(pool)
            pos = random.randint(s + 10, e - READ_LEN - 10)
            r = make_read(name, ref_id[WORK_CHROM], pos, st, cb, u)
        elif cat == "AMBIGUOUS_COD_NCOD":
            pool = picks["cod_ncod_overlaps"]
            if not pool:
                raise RuntimeError("no coding/non-coding gene-body overlaps found on this chromosome")
            s, e, st = random.choice(pool)
            pos = random.randint(s + 10, e - READ_LEN - 10)
            r = make_read(name, ref_id[WORK_CHROM], pos, st, cb, u)
        elif cat == "AMBIGUOUS":
            # Plain AMBIGUOUS: classifier returns AMBIGUOUS when read.get_blocks()
            # is empty (no aligned bases). A pure soft-clip CIGAR satisfies this
            # while keeping the read mapped (reference_start set).
            s, e = random.choice(picks["safe_plus_exons"])
            pos = random.randint(s + 10, e - READ_LEN - 10)
            r = make_read(name, ref_id[WORK_CHROM], pos, "+", cb, u, cigar="100S")
        elif cat == "UNASSIGNED":
            # Read with a CB that is intentionally absent from whitelist.txt;
            # the classifier short-circuits to UNASSIGNED before any
            # interval-based classification runs. Genomic position is irrelevant.
            s, e = random.choice(picks["safe_plus_exons"])
            pos = random.randint(s + 10, e - READ_LEN - 10)
            r = make_read(name, ref_id[WORK_CHROM], pos, "+", cb, u)
        else:
            raise ValueError(f"unknown simple category: {cat}")
        reads.append(r)
        ground_truth.append((name, cb, cat, notes))

    def emit_flag(cat: str, cb: str, read_idx: int, notes: str = ""):
        """Emit one EXONIC_SENSE read carrying a TSO or polyA artifact flag."""
        name = f"{cat}_{cb}_{read_idx:06d}"
        u = next_umi()
        if cat == "TSO_INVASION":
            s, e = random.choice(picks["safe_plus_exons"])
            pos = random.randint(s + 10, e - READ_LEN - 10)
            # 12 bp 5' soft-clip = TSO_10X[:12], followed by READ_LEN aligned bases.
            seq = TSO_SOFTCLIP_SEED + ("N" * READ_LEN)
            r = make_read(name, ref_id[WORK_CHROM], pos, "+", cb, u,
                          cigar=f"12S{READ_LEN}M", seq=seq)
        elif cat == "POLYA_PRIMING":
            pool = picks["polya_positions"]
            if not pool:
                raise RuntimeError("no POLYA_PRIMING positions found on this chromosome")
            pos = pool[read_idx % len(pool)]
            r = make_read(name, ref_id[WORK_CHROM], pos, "+", cb, u)
        else:
            raise ValueError(f"unknown flag category: {cat}")
        reads.append(r)
        ground_truth.append((name, cb, cat, notes))

    def emit_clustered(cat: str, cb: str, seed_pos: int, read_idx: int, notes: str = ""):
        """Emit one read at a clustered intergenic locus."""
        name = f"{cat}_{cb}_{read_idx:06d}"
        u = next_umi()
        if cat == "INTERGENIC_HOTSPOT":
            # Place at exact seed so the modal 3' end lands on the verified
            # A-run downstream. Any jitter would split the modal vote.
            r = make_read(name, ref_id[WORK_CHROM], seed_pos, "+", cb, u,
                          cigar=f"{READ_LEN}M")
        elif cat == "INTERGENIC_NOVEL":
            # Place at exact seed so every read's donor dinucleotide
            # (reference[seed+50, seed+52]) is the verified non-canonical one.
            r = make_read(name, ref_id[WORK_CHROM], seed_pos, "+", cb, u,
                          cigar=f"50M{NOVEL_SKIP}N50M")
        elif cat == "INTERGENIC_REPEAT":
            # REPEAT only needs the locus bounds to overlap the BED interval;
            # jitter is fine and varies the cluster footprint slightly.
            jitter = random.randint(-CLUSTERED_JITTER, CLUSTERED_JITTER)
            pos = max(0, seed_pos + jitter)
            r = make_read(name, ref_id[WORK_CHROM], pos, "+", cb, u,
                          cigar=f"{READ_LEN}M")
        else:
            raise ValueError(f"unknown clustered category: {cat}")
        reads.append(r)
        ground_truth.append((name, cb, cat, notes))

    # ----- Exp 1 simple categories -----
    print(f"[3/8] Generating Exp 1 simple categories "
          f"({N_PER_CATEGORY_EXP1} reads × {len(EXP1_SIMPLE_CATEGORIES)} CBs)...",
          file=sys.stderr)
    for cat in EXP1_SIMPLE_CATEGORIES:
        cb = f"CELL_EXP1_{cat}"
        for i in range(N_PER_CATEGORY_EXP1):
            emit_simple(cat, cb, i, notes="exp1")

    # ----- Exp 1 clustered intergenic categories -----
    print(f"[4/8] Generating Exp 1 clustered categories "
          f"({CLUSTERED_N_CBS} CBs × {CLUSTERED_N_LOCI} loci × "
          f"{CLUSTERED_READS_PER_LOCUS_PER_CB} reads × "
          f"{len(EXP1_CLUSTERED_CATEGORIES)} categories)...", file=sys.stderr)
    for cat in EXP1_CLUSTERED_CATEGORIES:
        seeds = picks[{
            "INTERGENIC_HOTSPOT": "hotspot_seeds",
            "INTERGENIC_NOVEL":   "novel_seeds",
            "INTERGENIC_REPEAT":  "repeat_seeds",
        }[cat]]
        if not seeds:
            print(f"      SKIP {cat}: no seed positions", file=sys.stderr)
            continue
        for cb_idx in range(CLUSTERED_N_CBS):
            cb = f"CELL_EXP1_{cat}_{cb_idx}"
            read_idx = 0
            for seed in seeds:
                for _ in range(CLUSTERED_READS_PER_LOCUS_PER_CB):
                    emit_clustered(cat, cb, seed, read_idx, notes="exp1")
                    read_idx += 1

    # ----- Exp 4 artifact-flag cells (TSO_INVASION, POLYA_PRIMING) -----
    print(f"[5/8] Generating Exp 4 artifact-flag cells "
          f"({N_PER_CATEGORY_EXP1} reads × {len(EXP4_FLAG_CATEGORIES)} CBs)...",
          file=sys.stderr)
    for cat in EXP4_FLAG_CATEGORIES:
        cb = f"CELL_EXP4_{cat}"
        for i in range(N_PER_CATEGORY_EXP1):
            emit_flag(cat, cb, i, notes=f"exp4;flag={cat.lower()}")

    # ----- Exp 2 mixture cells -----
    print(f"[6/8] Generating Exp 2 mixture cells ({len(NOISE_FRACTIONS)} CBs, "
          f"{MIXTURE_CELL_N_READS} reads each)...", file=sys.stderr)
    for noise_frac in NOISE_FRACTIONS:
        cb = f"CELL_EXP2_NOISE{int(noise_frac * 100):02d}"
        n_noise = int(round(MIXTURE_CELL_N_READS * noise_frac))
        n_signal = MIXTURE_CELL_N_READS - n_noise
        for i in range(n_signal):
            emit_simple("EXONIC_SENSE", cb, i,
                        notes=f"exp2;true_noise={noise_frac:.2f}")
        remaining = n_noise
        composition_items = list(NOISE_COMPOSITION.items())
        for j, (noise_cat, frac) in enumerate(composition_items):
            n = int(round(n_noise * frac)) if j < len(composition_items) - 1 else remaining
            remaining -= n
            for i in range(n):
                emit_simple(noise_cat, cb, n_signal + j * 10_000 + i,
                            notes=f"exp2;true_noise={noise_frac:.2f};pool={noise_cat}")

    # ----- Write BAM + supporting files -----
    print(f"[7/8] Writing BAM ({len(reads)} reads)...", file=sys.stderr)
    reads.sort(key=lambda r: (r.reference_id, r.reference_start))

    bam_path = outdir / "synthetic.bam"
    tmp_path = outdir / "synthetic.unsorted.bam"
    with pysam.AlignmentFile(str(tmp_path), "wb", header=header) as bam_out:
        for r in reads:
            bam_out.write(r)
    pysam.sort("-o", str(bam_path), str(tmp_path))
    pysam.index(str(bam_path))
    os.remove(tmp_path)

    gt_path = outdir / "ground_truth.tsv"
    with open(gt_path, "w") as fh:
        fh.write("read_name\tcell_barcode\ttrue_category\tnotes\n")
        for row in ground_truth:
            fh.write("\t".join(row) + "\n")

    # Emit repeats.bed covering each REPEAT seed ± buffer. The intergenic
    # profiler's _overlaps_repeats helper assumes the per-contig interval list
    # is sorted by start and break-outs on ivl_start >= end, so we sort here.
    repeats_path = outdir / "repeats.bed"
    repeat_intervals = sorted(
        (max(0, seed - CLUSTERED_JITTER - 50),
         seed + CLUSTERED_JITTER + READ_LEN + 50)
        for seed in picks["repeat_seeds"]
    )
    with open(repeats_path, "w") as fh:
        for start, end in repeat_intervals:
            fh.write(f"{WORK_CHROM}\t{start}\t{end}\tSynth_Repeat\t0\t+\n")

    # ----- Exp 5: Illumina paired-end CHIMERIC BAM -----
    print(f"[8/8] Generating Exp 5 Illumina paired-end BAM "
          f"({len(ILLUMINA_CATEGORIES)} CBs)...", file=sys.stderr)
    _, illumina_gt = synthesize_illumina_bam(outdir, header, ref_id, picks, refs)

    if illumina_gt:
        gt_illumina_path = outdir / "ground_truth_illumina.tsv"
        with open(gt_illumina_path, "w") as fh:
            fh.write("read_name\tcell_barcode\ttrue_category\tnotes\n")
            for row in illumina_gt:
                fh.write("\t".join(row) + "\n")
        print(f"      -> {gt_illumina_path}", file=sys.stderr)

    # Emit barcode whitelist: every CB appearing in ground_truth (ONT + Illumina)
    # EXCEPT the UNASSIGNED CB. The runner must pass this via --barcode-whitelist
    # so UNASSIGNED reads trip the not-on-whitelist branch in the classifier.
    whitelist_path = outdir / "whitelist.txt"
    all_cbs = {row[1] for row in ground_truth} | {row[1] for row in illumina_gt}
    with open(whitelist_path, "w") as fh:
        for cb in sorted(all_cbs):
            if cb == UNASSIGNED_CB:
                continue
            fh.write(cb + "\n")

    fa.close()

    print(f"      -> {bam_path}", file=sys.stderr)
    print(f"      -> {bam_path}.bai", file=sys.stderr)
    print(f"      -> {gt_path}", file=sys.stderr)
    print(f"      -> {repeats_path}", file=sys.stderr)
    print(f"      -> {whitelist_path}", file=sys.stderr)
    print(f"done.", file=sys.stderr)


if __name__ == "__main__":
    main()
