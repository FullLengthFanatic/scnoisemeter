#!/usr/bin/env python3
"""
synthesize_bam.py — controlled ground-truth BAM for scNoiseMeter v1 benchmark.

For each ReadCategory we validate in v1, emit N reads with synthetic CB/UB tags.
Two experiments share one BAM:

  Exp 1 — per-read classifier accuracy (confusion matrix).
          One CB per true category, N reads per CB.
          After running scnoisemeter, the per-cell counts for each CB reveal
          the classifier's category-assignment distribution for that
          category's reads.

  Exp 2 — per-cell noise-fraction fidelity (dose-response).
          Five CBs with known mixtures of signal (EXONIC_SENSE) and noise
          (drawn from NOISE_CATEGORIES_CONSERVATIVE: INTRONIC_PURE,
          INTRONIC_BOUNDARY, INTERGENIC_SPARSE, EXONIC_ANTISENSE).

Categories covered in v1:
  EXONIC_SENSE, EXONIC_ANTISENSE, INTRONIC_PURE, INTRONIC_BOUNDARY,
  INTERGENIC_SPARSE, MITOCHONDRIAL, MULTIMAPPER, CHIMERIC.

Categories deferred to v2:
  INTRONIC_JXNSPAN, INTERGENIC_HOTSPOT, INTERGENIC_NOVEL, INTERGENIC_REPEAT,
  AMBIGUOUS_*.

Key correctness invariants used when picking positions:
  * Strand-safe + exons: + strand exon intervals with no base overlap to any
    - strand exon. Ensures EXONIC_ANTISENSE reads do not become EXONIC_SENSE
    relative to an overlapping - strand gene.
  * True introns: regions inside a + strand gene and outside the *union of
    all annotated exons* on any strand. Ensures INTRONIC_PURE and the intronic
    portion of INTRONIC_BOUNDARY reads are not exonic in any transcript.
  * Clean intergenic: ≥10 kb from any gene and not inside any exon.
"""

from __future__ import annotations

import argparse
import gzip
import os
import random
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

EXP1_CATEGORIES = [
    "EXONIC_SENSE",
    "EXONIC_ANTISENSE",
    "INTRONIC_PURE",
    "INTRONIC_BOUNDARY",
    "INTERGENIC_SPARSE",
    "MITOCHONDRIAL",
    "MULTIMAPPER",
    "CHIMERIC",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtf", required=True, help="GENCODE GTF (gzipped ok)")
    p.add_argument("--fasta", required=True, help="Reference FASTA with .fai (gzipped+bgzf ok)")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------

def parse_gtf(gtf_path, chrom):
    """Return (genes, exons) on `chrom`. Each is list of (start, end, strand, gene_id).
    Coordinates 0-based half-open."""
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    genes, exons = [], []
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
            gene_id = None
            for kv in cols[8].split(";"):
                kv = kv.strip()
                if kv.startswith("gene_id"):
                    parts = kv.split('"')
                    if len(parts) >= 2:
                        gene_id = parts[1]
                    break
            if feat == "gene":
                genes.append((start, end, strand, gene_id))
            else:
                exons.append((start, end, strand, gene_id))
    return genes, exons


# ---------------------------------------------------------------------------
# Interval utilities (sorted merged intervals + bisect)
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
    """True if [s, e) overlaps any interval in sorted non-overlapping merged."""
    if not merged:
        return False
    i = bisect_right(starts, e - 1) - 1
    if i < 0:
        return False
    return merged[i][1] > s


def _fully_inside(merged, starts, s, e):
    """True if [s, e) is fully contained within one interval of merged."""
    if not merged:
        return False
    i = bisect_right(starts, s) - 1
    if i < 0:
        return False
    ms, me = merged[i]
    return ms <= s and e <= me


# ---------------------------------------------------------------------------
# Region pickers (strand-safe, exon-aware)
# ---------------------------------------------------------------------------

def build_all_pickers(genes, exons, chrom_size, buffer=10_000):
    """Return a dict of region lists for each synthetic category."""
    plus_exons_raw = [(s, e) for s, e, st, _ in exons if st == "+"]
    minus_exons_raw = [(s, e) for s, e, st, _ in exons if st == "-"]
    all_exons_raw = [(s, e) for s, e, _, _ in exons]

    plus_exon_union = _merge(plus_exons_raw)
    minus_exon_union = _merge(minus_exons_raw)
    all_exon_union = _merge(all_exons_raw)

    minus_starts = _starts(minus_exon_union)
    all_exon_starts = _starts(all_exon_union)
    plus_exon_starts = _starts(plus_exon_union)

    # --- strand-safe + exons (no overlap with any - strand exon) ---
    safe_plus_exons = []
    for s, e in plus_exon_union:
        if (e - s) < READ_LEN + 40:
            continue
        if _overlaps(minus_exon_union, minus_starts, s, e):
            continue
        safe_plus_exons.append((s, e))

    # --- + strand introns: inside a + strand gene footprint, NOT inside any exon ---
    plus_gene_union = _merge([(s, e) for s, e, st, _ in genes if st == "+"])
    plus_gene_starts = _starts(plus_gene_union)

    # Compute the complement of all_exon_union within plus_gene_union (per + gene interval).
    plus_introns: list[tuple[int, int]] = []
    for gs, ge in plus_gene_union:
        # Find exon intervals fully or partially inside [gs, ge).
        # Walk exon intervals overlapping [gs, ge).
        i = bisect_right(all_exon_starts, ge - 1) - 1
        # Collect all exon intervals whose end > gs.
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
    # Only keep long-enough intron regions
    plus_introns = [(s, e) for s, e in plus_introns if (e - s) >= READ_LEN + 100]

    # --- INTRONIC_BOUNDARY candidates ---
    # Classifier picks the dominant category by argmax(base_counts). With a
    # 50/50 split the dict-insertion-order tie resolves to EXONIC_SENSE. Use
    # a 30/70 exon/intron split so INTRONIC_BOUNDARY wins the argmax.
    # exon_end positions on strand-safe + exons where:
    #   [exon_end-30, exon_end) is fully within this same + strand exon,
    #   [exon_end, exon_end+70) has no overlap with any annotated exon, and
    #   [exon_end, exon_end+70) is fully inside a + strand gene (so the
    #     exon-free region is intronic, not intergenic — without this the
    #     70 bp downstream would be classified INTERGENIC_SPARSE).
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

    # --- Intergenic: ≥buffer from any gene AND not inside any exon ---
    gene_union = _merge([(s, e) for s, e, _, _ in genes])
    buffered = [(max(0, s - buffer), min(chrom_size, e + buffer)) for s, e in gene_union]
    buffered_union = _merge(buffered)
    buffered_starts = _starts(buffered_union)

    intergenic_positions: list[int] = []
    stride = 2000
    prev_end = 0
    for s, e in buffered_union:
        for p in range(prev_end, s - READ_LEN - 50, stride):
            if not _overlaps(all_exon_union, all_exon_starts, p, p + READ_LEN):
                intergenic_positions.append(p)
        prev_end = max(prev_end, e)
    for p in range(prev_end, chrom_size - READ_LEN - 50, stride):
        if not _overlaps(all_exon_union, all_exon_starts, p, p + READ_LEN):
            intergenic_positions.append(p)

    random.shuffle(intergenic_positions)

    return {
        "safe_plus_exons": safe_plus_exons,
        "plus_introns": plus_introns,
        "boundary_positions": boundary_positions,
        "intergenic_positions": intergenic_positions,
        "all_exon_union": all_exon_union,
    }


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
    read_len: int = READ_LEN,
) -> pysam.AlignedSegment:
    a = pysam.AlignedSegment()
    a.query_name = name
    a.reference_id = ref_id
    a.reference_start = pos
    a.mapping_quality = mapq
    a.flag = 16 if strand == "-" else 0
    a.cigarstring = cigar
    a.query_sequence = "N" * read_len
    a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
    tags: list[tuple[str, object]] = [("CB", cb), ("UB", umi), ("NH", nh)]
    if sa is not None:
        tags.append(("SA", sa))
    a.tags = tags
    return a


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    random.seed(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[1/5] Parsing GTF for {WORK_CHROM}...", file=sys.stderr)
    genes, exons = parse_gtf(args.gtf, WORK_CHROM)
    print(f"      {len(genes)} genes, {len(exons)} exon entries", file=sys.stderr)

    print(f"[2/5] Reading FASTA index...", file=sys.stderr)
    fa = pysam.FastaFile(args.fasta)
    refs = list(fa.references)
    lens = [fa.get_reference_length(r) for r in refs]
    fa.close()
    chrom_size = lens[refs.index(WORK_CHROM)]
    mito_size = lens[refs.index(MITO_CHROM)]

    print(f"      Building picker indices for {WORK_CHROM}...", file=sys.stderr)
    picks = build_all_pickers(genes, exons, chrom_size)
    print(f"      safe_plus_exons={len(picks['safe_plus_exons'])}, "
          f"plus_introns={len(picks['plus_introns'])}, "
          f"boundary_positions={len(picks['boundary_positions'])}, "
          f"intergenic_positions={len(picks['intergenic_positions'])}", file=sys.stderr)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": r, "LN": l} for r, l in zip(refs, lens)],
        "PG": [{"ID": "synthesize_bam", "PN": "synthesize_bam", "VN": "0.1"}],
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

    def emit_read_for_category(cat: str, cb: str, read_idx: int, notes: str = ""):
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
        else:
            raise ValueError(f"unknown category: {cat}")
        reads.append(r)
        ground_truth.append((name, cb, cat, notes))

    print(f"[3/5] Generating Exp 1 reads ({N_PER_CATEGORY_EXP1} per category, {len(EXP1_CATEGORIES)} CBs)...",
          file=sys.stderr)
    for cat in EXP1_CATEGORIES:
        cb = f"CELL_EXP1_{cat}"
        for i in range(N_PER_CATEGORY_EXP1):
            emit_read_for_category(cat, cb, i, notes="exp1")

    print(f"[4/5] Generating Exp 2 mixture cells ({len(NOISE_FRACTIONS)} CBs, "
          f"{MIXTURE_CELL_N_READS} reads each)...", file=sys.stderr)
    for noise_frac in NOISE_FRACTIONS:
        cb = f"CELL_EXP2_NOISE{int(noise_frac * 100):02d}"
        n_noise = int(round(MIXTURE_CELL_N_READS * noise_frac))
        n_signal = MIXTURE_CELL_N_READS - n_noise
        for i in range(n_signal):
            emit_read_for_category("EXONIC_SENSE", cb, i,
                                   notes=f"exp2;true_noise={noise_frac:.2f}")
        remaining = n_noise
        composition_items = list(NOISE_COMPOSITION.items())
        for j, (noise_cat, frac) in enumerate(composition_items):
            n = int(round(n_noise * frac)) if j < len(composition_items) - 1 else remaining
            remaining -= n
            for i in range(n):
                emit_read_for_category(noise_cat, cb, n_signal + j * 10_000 + i,
                                       notes=f"exp2;true_noise={noise_frac:.2f};pool={noise_cat}")

    print(f"[5/5] Writing BAM ({len(reads)} reads)...", file=sys.stderr)
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

    print(f"      -> {bam_path}", file=sys.stderr)
    print(f"      -> {bam_path}.bai", file=sys.stderr)
    print(f"      -> {gt_path}", file=sys.stderr)
    print(f"done.", file=sys.stderr)


if __name__ == "__main__":
    main()
