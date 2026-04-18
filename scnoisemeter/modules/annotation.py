"""
annotation.py
=============
Builds the strand-aware genomic interval index used by the read classifier.

From a GENCODE (or compatible) GTF this module produces five pyranges
objects that together tile the genome exhaustively:

  exons          – union of all annotated exonic bases, per strand
  introns        – gene-body bases not covered by any exon, per strand
  intergenic     – complement of all gene bodies (strand-independent)
  gene_unique    – exonic/intronic bases attributable to exactly one gene
  gene_shared    – bases where two or more genes overlap (→ ambiguous)

Additionally a RepeatMasker interval set is loaded if provided, used to
sub-classify intergenic reads.

Design decisions
----------------
1.  Exon union is computed per gene (not per transcript), so a base covered
    by any exon of any transcript of a gene is "exonic".  This is the most
    permissive definition and minimises false "intronic" calls for genes with
    complex alternative splicing.

2.  Introns are gene-body minus exon-union, per gene.  This means an intron
    that is exonic in one transcript but intronic in another is treated as
    exonic (consistent with decision 1).

3.  Overlapping genes: we compute per-gene unique regions by subtracting
    all other gene bodies.  A read that falls entirely within a unique region
    is unambiguously attributable to one gene.  A read overlapping a shared
    region is flagged as ambiguous.

4.  Strand is handled explicitly throughout.  All interval sets carry a
    Strand column ('+' or '-').  Reads are matched against the strand-
    matching set first, then the anti-strand set.

5.  Pseudogenes and lncRNAs are included by default (they are in GENCODE
    and a read overlapping them is not "intergenic").  Users can filter
    biotypes via --exclude-biotypes.

6.  The mitochondrial contig is identified by name (chrM / MT) and handled
    separately by the classifier; it is excluded from these interval sets to
    avoid double-counting.

Performance
-----------
pyranges operations on GENCODE v45 (~2.9M exon records) complete in under
60 seconds on a modern laptop.  The result is cached to a compressed pickle
alongside the GTF so that subsequent runs skip the build step.
"""

from __future__ import annotations

import hashlib
import logging
import pickle
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd
import pyranges as pr

from scnoisemeter.constants import MITO_CONTIG_NAMES

logger = logging.getLogger(__name__)

# Bump this when the interval build logic changes so stale caches are invalidated
_CACHE_VERSION = "v7"


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class AnnotationIndex:
    """
    All interval sets required by the read classifier.

    All PyRanges objects carry at minimum: Chromosome, Start, End, Strand,
    gene_id, gene_name, gene_biotype columns.

    Attributes
    ----------
    exons_plus / exons_minus
        Exonic intervals on + and − strand respectively.
        These are the *union* per gene (exon bodies merged within each gene).
    introns_plus / introns_minus
        Intronic intervals (gene body minus exon union), per strand.
    gene_bodies_plus / gene_bodies_minus
        Full gene-body intervals (TSS → TES), per strand.
    gene_unique_plus / gene_unique_minus
        Gene-body bases attributable to exactly one gene, per strand.
    gene_shared
        Gene-body bases where ≥2 genes overlap (strand-aware).
    intergenic
        Complement of all gene bodies on either strand (strand-agnostic).
    repeats
        RepeatMasker intervals (strand-agnostic).  None if not provided.
    splice_sites
        Known donor/acceptor positions used for junction-spanning intronic
        sub-classification.  Dict: contig → set of (pos, strand) tuples.
    gtf_path : Path
    repeats_path : Optional[Path]
    """

    exons_plus:          pr.PyRanges
    exons_minus:         pr.PyRanges
    introns_plus:        pr.PyRanges
    introns_minus:       pr.PyRanges
    gene_bodies_plus:    pr.PyRanges
    gene_bodies_minus:   pr.PyRanges
    gene_unique_plus:    pr.PyRanges
    gene_unique_minus:   pr.PyRanges
    gene_shared_cod_cod:  pr.PyRanges   # coding vs coding overlaps (truly ambiguous)
    gene_shared_cod_ncod: pr.PyRanges   # coding vs noncoding overlaps (sub-classifiable)
    intergenic:          pr.PyRanges
    repeats:             Optional[pr.PyRanges]
    splice_sites:        dict  # contig → set[tuple[int, str]]
    gtf_path:            Path
    repeats_path:        Optional[Path]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_annotation_index(
    gtf_path: str | Path,
    *,
    repeats_path: Optional[str | Path] = None,
    exclude_biotypes: Optional[list[str]] = None,
    cache: bool = True,
) -> AnnotationIndex:
    """
    Build (or load from cache) the full annotation index.

    Parameters
    ----------
    gtf_path:
        Path to a GENCODE-format GTF (gzip-compressed or plain).
    repeats_path:
        Optional path to a RepeatMasker BED file for hg38
        (e.g. downloaded from UCSC Table Browser).
    exclude_biotypes:
        List of gene_biotype strings to exclude from the index
        (e.g. ["artifact", "TEC"]).  Default excludes nothing.
    cache:
        If True, write/read a compressed pickle next to the GTF.
    """
    gtf_path = Path(gtf_path)
    repeats_path = Path(repeats_path) if repeats_path else None
    exclude_biotypes = exclude_biotypes or []

    cache_path = _cache_path(gtf_path, exclude_biotypes)

    if cache and cache_path.exists():
        logger.info("Loading annotation index from cache: %s", cache_path)
        return _load_cache(cache_path)

    logger.info("Building annotation index from GTF: %s", gtf_path.name)
    index = _build(gtf_path, repeats_path, exclude_biotypes)

    if cache:
        _save_cache(index, cache_path)

    return index


# ---------------------------------------------------------------------------
# Build logic
# ---------------------------------------------------------------------------

def _build(
    gtf_path: Path,
    repeats_path: Optional[Path],
    exclude_biotypes: list[str],
) -> AnnotationIndex:

    # ------------------------------------------------------------------
    # 1. Load GTF
    # ------------------------------------------------------------------
    logger.info("  Reading GTF …")
    raw = pr.read_gtf(str(gtf_path))

    # Drop mitochondrial contigs — handled separately by the classifier
    mito_mask = raw.df["Chromosome"].isin(MITO_CONTIG_NAMES)
    if mito_mask.any():
        logger.debug("  Dropping %d mito records from annotation index", mito_mask.sum())
        raw = pr.PyRanges(raw.df[~mito_mask])

    # Normalise biotype column name: pyranges 0.1.x uses 'gene_type' (GENCODE original),
    # older pyranges versions used 'gene_biotype'.  We standardise to 'gene_biotype'.
    _df = raw.df.copy()
    if "gene_type" in _df.columns and "gene_biotype" not in _df.columns:
        _df = _df.rename(columns={"gene_type": "gene_biotype"})
        raw = pr.PyRanges(_df)
        logger.debug("  Renamed 'gene_type' → 'gene_biotype' for compatibility.")

    # Biotype filtering
    if exclude_biotypes and "gene_biotype" in raw.df.columns:
        before = len(raw.df)
        raw = pr.PyRanges(raw.df[~raw.df["gene_biotype"].isin(exclude_biotypes)])
        logger.debug("  Dropped %d records for excluded biotypes", before - len(raw.df))

    # ------------------------------------------------------------------
    # 2. Separate exon records
    # ------------------------------------------------------------------
    exon_df = raw.df[raw.df["Feature"] == "exon"].copy()
    gene_df = raw.df[raw.df["Feature"] == "gene"].copy()

    if exon_df.empty:
        raise ValueError("GTF contains no 'exon' feature rows — check GTF format.")
    if gene_df.empty:
        raise ValueError("GTF contains no 'gene' feature rows — check GTF format.")

    # Ensure required columns
    for col in ("gene_id", "Strand"):
        if col not in exon_df.columns:
            raise ValueError(f"GTF exon records lack column '{col}'.")

    # ------------------------------------------------------------------
    # 3. Per-gene exon union (merge overlapping exons within each gene)
    # ------------------------------------------------------------------
    logger.info("  Computing per-gene exon union …")
    exon_plus_df  = exon_df[exon_df["Strand"] == "+"]
    exon_minus_df = exon_df[exon_df["Strand"] == "-"]

    exons_plus  = _merge_exons_per_gene(exon_plus_df)
    exons_minus = _merge_exons_per_gene(exon_minus_df)

    # ------------------------------------------------------------------
    # 4. Gene bodies
    # ------------------------------------------------------------------
    # Build column list defensively — gene_name and gene_biotype are optional
    _gene_base_cols = ["Chromosome", "Start", "End", "Strand", "gene_id"]
    for _opt in ("gene_name", "gene_biotype"):
        if _opt in gene_df.columns:
            _gene_base_cols.append(_opt)
        else:
            logger.debug("GTF does not contain column '%s' — skipping.", _opt)
    gene_plus_df  = gene_df[gene_df["Strand"] == "+"][_gene_base_cols].drop_duplicates()
    gene_minus_df = gene_df[gene_df["Strand"] == "-"][_gene_base_cols].drop_duplicates()

    gene_bodies_plus  = pr.PyRanges(gene_plus_df)
    gene_bodies_minus = pr.PyRanges(gene_minus_df)

    # ------------------------------------------------------------------
    # 5. Introns = gene body minus exon union, per gene
    # ------------------------------------------------------------------
    logger.info("  Computing introns …")
    introns_plus  = _compute_introns(gene_bodies_plus,  exons_plus)
    introns_minus = _compute_introns(gene_bodies_minus, exons_minus)

    # ------------------------------------------------------------------
    # 6. Unique vs shared gene regions
    # ------------------------------------------------------------------
    logger.info("  Computing gene unique / shared regions …")
    # IMPORTANT: shared/unique is computed on EXONIC intervals, not gene bodies.
    # A read is genuinely ambiguous only when it overlaps the exons of 2+ genes.
    # Intronic overlaps between gene bodies are not ambiguous for noise quantification
    # — a read in such a region is simply intronic, regardless of which gene's
    # intron it falls in.  Computing on gene bodies inflates ambiguous 10-fold
    # because large genes (e.g. AGRN, 100 kb) have many other genes nested in
    # their introns in GENCODE, creating apparent overlap where none exists at
    # the exon level.
    gene_unique_plus,  shared_cod_cod_plus,  shared_cod_ncod_plus  = _unique_and_shared(exons_plus)
    gene_unique_minus, shared_cod_cod_minus, shared_cod_ncod_minus = _unique_and_shared(exons_minus)

    # Merge strand-specific shared sets — kept separate by type
    def _merge_strands(pr_a, pr_b):
        df = pd.concat(
            [pr_a.df if not pr_a.df.empty else pd.DataFrame(),
             pr_b.df if not pr_b.df.empty else pd.DataFrame()],
            ignore_index=True
        )
        return pr.PyRanges(df).merge() if not df.empty else pr.PyRanges()

    gene_shared_cod_cod  = _merge_strands(shared_cod_cod_plus,  shared_cod_cod_minus)
    gene_shared_cod_ncod = _merge_strands(shared_cod_ncod_plus, shared_cod_ncod_minus)

    # ------------------------------------------------------------------
    # 7. Intergenic = complement of all gene bodies (either strand)
    # ------------------------------------------------------------------
    logger.info("  Computing intergenic complement …")
    all_gene_bodies_df = pd.concat(
        [gene_plus_df, gene_minus_df], ignore_index=True
    )[["Chromosome", "Start", "End"]].drop_duplicates()

    # pyranges 0.1.x removed complement() — compute intergenic gaps manually
    intergenic = _manual_complement(all_gene_bodies_df)

    # ------------------------------------------------------------------
    # 8. Splice sites
    # ------------------------------------------------------------------
    logger.info("  Extracting splice sites …")
    splice_sites = _extract_splice_sites(exon_df)

    # ------------------------------------------------------------------
    # 9. RepeatMasker (optional)
    # ------------------------------------------------------------------
    repeats = None
    if repeats_path is not None:
        logger.info("  Loading RepeatMasker BED: %s", repeats_path.name)
        repeats = _load_repeats(repeats_path)

    logger.info("  Annotation index build complete.")

    return AnnotationIndex(
        exons_plus=exons_plus,
        exons_minus=exons_minus,
        introns_plus=introns_plus,
        introns_minus=introns_minus,
        gene_bodies_plus=gene_bodies_plus,
        gene_bodies_minus=gene_bodies_minus,
        gene_unique_plus=gene_unique_plus,
        gene_unique_minus=gene_unique_minus,
        gene_shared_cod_cod=gene_shared_cod_cod,
        gene_shared_cod_ncod=gene_shared_cod_ncod,
        intergenic=intergenic,
        repeats=repeats,
        splice_sites=splice_sites,
        gtf_path=gtf_path,
        repeats_path=repeats_path,
    )


# ---------------------------------------------------------------------------
# Helper: per-gene exon union
# ---------------------------------------------------------------------------

def _merge_exons_per_gene(exon_df: pd.DataFrame) -> pr.PyRanges:
    """
    For each gene, merge all overlapping exon intervals into a non-redundant
    set.  Returns a PyRanges with columns:
      Chromosome, Start, End, Strand, gene_id, gene_name, gene_biotype

    Implementation uses a single vectorized PyRanges merge grouped by gene
    metadata columns, which is orders of magnitude faster than a per-gene
    Python loop on large GTFs (e.g. Ensembl with ~60 k genes).
    """
    if exon_df.empty:
        return pr.PyRanges()

    keep_cols = ["Chromosome", "Start", "End", "Strand", "gene_id"]
    optional = ["gene_name", "gene_biotype"]
    for col in optional:
        if col in exon_df.columns:
            keep_cols.append(col)

    sub = exon_df[keep_cols].copy()

    # Vectorized merge: group by gene identity columns (gene_id + any optional
    # metadata that is constant per gene).  Chromosome and Strand are handled
    # by PyRanges internally (merging never crosses chromosome/strand boundaries).
    by_cols = [c for c in ["gene_id"] + optional if c in sub.columns]
    gr = pr.PyRanges(sub)
    return gr.merge(by=by_cols)


# ---------------------------------------------------------------------------
# Helper: introns
# ---------------------------------------------------------------------------

def _compute_introns(
    gene_bodies: pr.PyRanges,
    exons: pr.PyRanges,
) -> pr.PyRanges:
    """
    Introns = gene body intervals subtract exon-union intervals, per strand.
    The result retains gene_id / gene_name / gene_biotype from gene_bodies.
    """
    if gene_bodies.df.empty or exons.df.empty:
        return pr.PyRanges()

    # pyranges subtract removes exonic bases from gene bodies
    introns = gene_bodies.subtract(exons)
    return introns


# ---------------------------------------------------------------------------
# Helper: unique vs shared regions
# ---------------------------------------------------------------------------

# Non-coding biotypes: when these overlap a higher-priority gene,
# they do NOT create a truly ambiguous region — the read is sub-classified
# as AMBIGUOUS_COD_NCOD rather than being unresolvable.
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


def _classify_overlap_type(biotype_a: str, biotype_b: str) -> str:
    """
    Classify the type of overlap between two genes by their biotypes.

    Returns one of:
      "cod_cod"   — both coding → genuinely ambiguous
      "cod_ncod"  — one coding, one non-coding → sub-classifiable
      "ncod_ncod" — both non-coding → treat as ambiguous (conservative)
    """
    a_coding = biotype_a in _CODING_BIOTYPES
    b_coding = biotype_b in _CODING_BIOTYPES
    if a_coding and b_coding:
        return "cod_cod"
    if a_coding or b_coding:
        return "cod_ncod"
    return "ncod_ncod"


def _unique_and_shared(gene_bodies: pr.PyRanges) -> tuple[pr.PyRanges, pr.PyRanges, pr.PyRanges]:
    """
    Split gene body intervals into three sets:

      unique          — bases attributable to exactly one gene on that strand
      shared_cod_cod  — bases where ≥2 protein-coding genes overlap (truly ambiguous)
      shared_cod_ncod — bases where a coding and non-coding gene overlap
                        (sub-classifiable; read strand resolves in many cases)

    Key design decisions
    --------------------
    1. Strand-aware: shared regions are computed per strand.  A region shared
       between a + strand gene and a - strand gene is NOT ambiguous for a
       strand-specific read — the read's strand resolves it.

    2. Biotype-transparent: we do NOT silently assign coding-vs-noncoding
       overlaps to the coding gene.  Instead we flag them as AMBIGUOUS_COD_NCOD
       so the user can see the magnitude and decide how to handle them.

    3. The classifier uses both shared sets — AMBIGUOUS_COD_COD gets the
       highest ambiguity weight, AMBIGUOUS_COD_NCOD is reported separately.

    Returns (unique_pr, shared_cod_cod_pr, shared_cod_ncod_pr).
    """
    if gene_bodies.df.empty:
        return pr.PyRanges(), pr.PyRanges(), pr.PyRanges()

    if "gene_id" not in gene_bodies.df.columns:
        return gene_bodies, pr.PyRanges(), pr.PyRanges()

    df = gene_bodies.df.copy()
    has_biotype = "gene_biotype" in df.columns

    # Self-join to find all pairwise gene overlaps
    overlaps = gene_bodies.join(gene_bodies, suffix="_b")
    ov_df = overlaps.df

    # Keep only pairs of DIFFERENT genes
    diff = ov_df[ov_df["gene_id"] != ov_df["gene_id_b"]].copy()

    if diff.empty:
        return gene_bodies, pr.PyRanges(), pr.PyRanges()

    # Compute actual intersection interval for each pair
    diff["shared_start"] = diff[["Start", "Start_b"]].max(axis=1)
    diff["shared_end"]   = diff[["End",   "End_b"]].min(axis=1)

    # Classify each overlap by biotype
    if has_biotype and "gene_biotype_b" in diff.columns:
        diff["_overlap_type"] = diff.apply(
            lambda r: _classify_overlap_type(
                str(r.get("gene_biotype", "")),
                str(r.get("gene_biotype_b", ""))
            ), axis=1
        )
    else:
        # No biotype info — treat all overlaps as cod_cod (conservative)
        diff["_overlap_type"] = "cod_cod"

    # Build shared_cols defensively — Strand may not survive the join
    # in all pyranges 0.1.x configurations (e.g. when all entries share
    # one strand value, pyranges may treat it as implicit).
    has_strand = "Strand" in diff.columns
    shared_base = ["Chromosome", "shared_start", "shared_end"]
    if has_strand:
        shared_base.append("Strand")

    # Coding vs coding: genuinely ambiguous
    cod_cod = diff[diff["_overlap_type"] == "cod_cod"][shared_base].rename(
        columns={"shared_start": "Start", "shared_end": "End"}
    ).drop_duplicates()

    # Coding vs non-coding (and ncod vs ncod): sub-classifiable ambiguous
    ncod = diff[diff["_overlap_type"].isin(["cod_ncod", "ncod_ncod"])][shared_base].rename(
        columns={"shared_start": "Start", "shared_end": "End"}
    ).drop_duplicates()

    shared_cod_cod_pr  = pr.PyRanges(cod_cod).merge()  if not cod_cod.empty  else pr.PyRanges()
    shared_cod_ncod_pr = pr.PyRanges(ncod).merge()     if not ncod.empty     else pr.PyRanges()

    # Unique = all gene bodies minus ALL shared regions
    all_shared_df = pd.concat(
        [cod_cod, ncod], ignore_index=True
    ).drop_duplicates() if not (cod_cod.empty and ncod.empty) else pd.DataFrame()

    if all_shared_df.empty:
        unique_pr = gene_bodies
    else:
        all_shared_pr = pr.PyRanges(all_shared_df).merge()
        unique_pr = gene_bodies.subtract(all_shared_pr)

    return (
        unique_pr          if not unique_pr.df.empty          else pr.PyRanges(),
        shared_cod_cod_pr  if not shared_cod_cod_pr.df.empty  else pr.PyRanges(),
        shared_cod_ncod_pr if not shared_cod_ncod_pr.df.empty else pr.PyRanges(),
    )


# ---------------------------------------------------------------------------
# Helper: splice site extraction
# ---------------------------------------------------------------------------

def _extract_splice_sites(exon_df: pd.DataFrame) -> dict:
    """
    From the exon table, derive the set of known splice donor and acceptor
    positions for each contig.

    For a + strand gene:
      donor    = exon End   (first base of the downstream intron)
      acceptor = exon Start (last base of the upstream intron + 1)

    Returns: dict[contig_name, set[tuple[position_int, strand_str]]]
    """
    sites: dict[str, set] = {}
    for _, row in exon_df.iterrows():
        chrom  = row["Chromosome"]
        start  = int(row["Start"])
        end    = int(row["End"])
        strand = row["Strand"]
        if chrom not in sites:
            sites[chrom] = set()
        sites[chrom].add((start, strand))
        sites[chrom].add((end,   strand))
    return sites


# ---------------------------------------------------------------------------
# Helper: load RepeatMasker BED
# ---------------------------------------------------------------------------

def _load_repeats(repeats_path: Path) -> pr.PyRanges:
    """
    Load a RepeatMasker BED file.  Expects at least 6 columns:
      chrom, start, end, name, score, strand
    Returns a PyRanges with Chromosome, Start, End, Strand, repeat_name.
    """
    try:
        df = pd.read_csv(
            repeats_path,
            sep="\t",
            header=None,
            comment="#",
            usecols=[0, 1, 2, 3, 5],
            names=["Chromosome", "Start", "End", "repeat_name", "Strand"],
            dtype={"Start": int, "End": int},
            low_memory=False,
        )
        df["Strand"] = df["Strand"].fillna(".")
        return pr.PyRanges(df)
    except Exception as exc:
        logger.warning(
            "Could not load RepeatMasker BED (%s): %s — repeats layer disabled.",
            repeats_path, exc
        )
        return pr.PyRanges()


# ---------------------------------------------------------------------------
# Helper: manual complement (intergenic gaps)
# pyranges 0.1.x removed complement() so we compute gaps explicitly
# ---------------------------------------------------------------------------

def _manual_complement(gene_bodies_df: pd.DataFrame) -> pr.PyRanges:
    """
    Return the intergenic gaps between merged gene bodies as a PyRanges.

    For each chromosome, identifies the gaps between merged gene body intervals,
    including any region from position 0 to the first gene. The region after
    the last gene on each chromosome is not included because chromosome lengths
    are not available here; the read classifier handles post-last-gene reads
    via its residual fallback (any bases not in exons or introns are assigned
    INTERGENIC_SPARSE) so query correctness is unaffected.

    This object is used only as the denominator in the Poisson background model
    (compute_intergenic_bases). Excluding the post-last-gene tail introduces a
    small systematic underestimate of intergenic bases, which makes the Poisson
    test slightly conservative. The previous implementation added a sentinel
    End=2_000_000_000 per chromosome, inflating the denominator by ~24× and
    causing the test to flag almost every concentrated intergenic locus.
    """
    if gene_bodies_df.empty:
        return pr.PyRanges()

    # Merge overlapping gene bodies per chromosome first
    merged = (
        pr.PyRanges(gene_bodies_df.assign(Strand="+"))
        .merge()
        .df[["Chromosome", "Start", "End"]]
        .sort_values(["Chromosome", "Start"])
        .reset_index(drop=True)
    )

    gaps = []
    for chrom, grp in merged.groupby("Chromosome"):
        grp = grp.sort_values("Start")
        prev_end = 0
        for _, row in grp.iterrows():
            if int(row["Start"]) > prev_end:
                gaps.append({
                    "Chromosome": chrom,
                    "Start": prev_end,
                    "End":   int(row["Start"]),
                    "Strand": ".",
                })
            prev_end = max(prev_end, int(row["End"]))
        # Post-last-gene region omitted: chromosome lengths are unknown here.
        # Reads mapping there are still classified as INTERGENIC_SPARSE by the
        # classifier's residual fallback; only the Poisson denominator is affected.

    if not gaps:
        return pr.PyRanges()
    return pr.PyRanges(pd.DataFrame(gaps))


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------

def _cache_path(gtf_path: Path, exclude_biotypes: list[str]) -> Path:
    """Derive a versioned cache path next to the GTF."""
    key = f"{_CACHE_VERSION}:{':'.join(sorted(exclude_biotypes))}"
    digest = hashlib.md5(key.encode()).hexdigest()[:8]
    stem = gtf_path.stem.replace(".gtf", "").replace(".gz", "")
    return gtf_path.parent / f".scnoisemeter_{stem}_{digest}.cache.pkl.gz"


def _save_cache(index: AnnotationIndex, path: Path) -> None:
    import gzip
    logger.info("Saving annotation cache to %s", path)
    with gzip.open(path, "wb") as fh:
        pickle.dump(index, fh, protocol=5)


def _load_cache(path: Path) -> AnnotationIndex:
    import gzip
    with gzip.open(path, "rb") as fh:
        return pickle.load(fh)
