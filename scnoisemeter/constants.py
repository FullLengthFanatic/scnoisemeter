"""
constants.py
============
Single source of truth for all category names, SAM flag bits, tag names,
and tunable defaults used throughout scNoiseMeter.

Design note
-----------
Every category string and every default threshold lives here.  No other module
should hard-code these values.  This makes the taxonomy easy to audit, extend,
and document.
"""

from enum import Enum


# ---------------------------------------------------------------------------
# Read classification taxonomy
# ---------------------------------------------------------------------------

class ReadCategory(str, Enum):
    """
    Mutually-exclusive classification of every analysed, mapped primary
    alignment.  Unmapped, secondary and supplementary record totals are
    reported separately because they do not carry a comparable genomic
    location classification.

    Hierarchy (earlier entries take priority over later ones during assignment):

      UNMAPPED          – read did not align; excluded from all fractions
      SECONDARY         – SAM flag 0x100; skipped (redundant record)
      SUPPLEMENTARY     – SAM flag 0x800; routed to chimeric detector only
      MULTIMAPPER       – primary alignment with explicit multi-hit evidence
      MITOCHONDRIAL     – maps to the mitochondrial contig (chrM / MT)
      CHIMERIC          – SA tag with inter-chromosomal, strand-discordant,
                          or incompatible query/genomic order evidence
      EXONIC_SENSE      – overlaps ≥1 annotated exon base, correct strand
      EXONIC_ANTISENSE  – overlaps ≥1 annotated exon base, wrong strand
      INTRONIC_JXNSPAN  – intronic but CIGAR contains N near a splice site
                          (candidate intron-retention or NNC transcript)
      INTRONIC_PURE     – maps entirely within intron body, no junction signal
      INTRONIC_BOUNDARY – spans an exon–intron boundary without a splice op
                          (candidate incomplete RT)
      INTERGENIC_REPEAT – intergenic, overlaps RepeatMasker annotation
      INTERGENIC_HOTSPOT– intergenic, passes adaptive threshold but shows
                          internal-priming / A-rich-3′ signature
      INTERGENIC_NOVEL  – intergenic, passes adaptive threshold, strand-
                          consistent, near annotated polyA site (candidate
                          unannotated gene)
      INTERGENIC_ENRICHED – enriched intergenic locus without enough evidence
                            to call a hotspot or transcript candidate
      INTERGENIC_SPARSE – intergenic, below adaptive threshold (likely noise)
      AMBIGUOUS         – overlaps shared region of two or more genes
      UNASSIGNED        – CB tag absent or not on whitelist
    """

    UNMAPPED           = "unmapped"
    SECONDARY          = "secondary"
    SUPPLEMENTARY      = "supplementary"
    MULTIMAPPER        = "multimapper"
    CHIMERIC           = "chimeric"
    MITOCHONDRIAL      = "mitochondrial"
    EXONIC_SENSE       = "exonic_sense"
    EXONIC_ANTISENSE   = "exonic_antisense"
    INTRONIC_JXNSPAN   = "intronic_jxnspan"
    INTRONIC_PURE      = "intronic_pure"
    INTRONIC_BOUNDARY  = "intronic_boundary"
    INTERGENIC_REPEAT  = "intergenic_repeat"
    INTERGENIC_HOTSPOT = "intergenic_hotspot"
    INTERGENIC_NOVEL   = "intergenic_novel"
    INTERGENIC_ENRICHED = "intergenic_enriched"
    INTERGENIC_SPARSE  = "intergenic_sparse"
    AMBIGUOUS          = "ambiguous"           # catch-all (should be rare after fixes)
    AMBIGUOUS_COD_NCOD = "ambiguous_cod_ncod"  # protein-coding vs lncRNA/pseudogene
    AMBIGUOUS_COD_COD  = "ambiguous_cod_cod"   # two protein-coding genes
    UNASSIGNED         = "unassigned"


# Ordered list used when iterating over output columns
CATEGORY_ORDER = [
    ReadCategory.EXONIC_SENSE,
    ReadCategory.EXONIC_ANTISENSE,
    ReadCategory.INTRONIC_JXNSPAN,
    ReadCategory.INTRONIC_PURE,
    ReadCategory.INTRONIC_BOUNDARY,
    ReadCategory.INTERGENIC_SPARSE,
    ReadCategory.INTERGENIC_REPEAT,
    ReadCategory.INTERGENIC_HOTSPOT,
    ReadCategory.INTERGENIC_NOVEL,
    ReadCategory.INTERGENIC_ENRICHED,
    ReadCategory.CHIMERIC,
    ReadCategory.MITOCHONDRIAL,
    ReadCategory.MULTIMAPPER,
    ReadCategory.AMBIGUOUS,
    ReadCategory.AMBIGUOUS_COD_NCOD,
    ReadCategory.AMBIGUOUS_COD_COD,
    ReadCategory.UNASSIGNED,
    # NOTE: UNMAPPED, SECONDARY, SUPPLEMENTARY are intentionally excluded —
    # those reads are skipped before any aggregation loop that uses this list.
]

# Broad non-canonical alignment composition.  These categories are useful for
# diagnostic comparison, but they are not a causal estimate of technical
# noise: every member can contain genuine biology in an appropriate context.
BROAD_NONCANONICAL_CATEGORIES = {
    ReadCategory.EXONIC_ANTISENSE,
    ReadCategory.INTRONIC_PURE,
    ReadCategory.INTRONIC_BOUNDARY,
    ReadCategory.INTERGENIC_SPARSE,
    ReadCategory.INTERGENIC_REPEAT,
    ReadCategory.INTERGENIC_HOTSPOT,
    ReadCategory.INTERGENIC_ENRICHED,
    ReadCategory.CHIMERIC,
}

# Categories with positive alignment-level evidence for a possible artifact.
# Even these remain "candidates": chimeric alignments may be real fusions and
# A-rich hotspots can coincide with genuine transcript ends.
ARTIFACT_CANDIDATE_CATEGORIES = {
    ReadCategory.INTERGENIC_HOTSPOT,
    ReadCategory.CHIMERIC,
}

# Backwards-compatible aliases.  Output fields using the old ``noise_*`` names
# are retained through v0.7, but documentation labels them as deprecated.
NOISE_CATEGORIES = BROAD_NONCANONICAL_CATEGORIES
NOISE_CATEGORIES_CONSERVATIVE = BROAD_NONCANONICAL_CATEGORIES
NOISE_CATEGORIES_STRICT = ARTIFACT_CANDIDATE_CATEGORIES
NOISE_CATEGORIES_UNSTRANDED = (
    BROAD_NONCANONICAL_CATEGORIES - {ReadCategory.EXONIC_ANTISENSE}
)
NOISE_CATEGORIES_STRICT_UNSTRANDED = ARTIFACT_CANDIDATE_CATEGORIES

# Categories whose interpretation is ambiguous and reported separately
AMBIGUOUS_CATEGORIES = {
    ReadCategory.INTRONIC_JXNSPAN,
    ReadCategory.INTERGENIC_NOVEL,
    ReadCategory.INTERGENIC_ENRICHED,
    ReadCategory.AMBIGUOUS,
    ReadCategory.AMBIGUOUS_COD_NCOD,
    ReadCategory.AMBIGUOUS_COD_COD,
}


# ---------------------------------------------------------------------------
# SAM flag bits (bitwise masks)
# ---------------------------------------------------------------------------

class SamFlag:
    PAIRED              = 0x001
    PROPER_PAIR         = 0x002
    UNMAPPED            = 0x004
    MATE_UNMAPPED       = 0x008
    REVERSE_STRAND      = 0x010
    MATE_REVERSE        = 0x020
    READ1               = 0x040
    READ2               = 0x080
    SECONDARY           = 0x100
    QC_FAIL             = 0x200
    DUPLICATE           = 0x400
    SUPPLEMENTARY       = 0x800


# ---------------------------------------------------------------------------
# BAM tag names
# ---------------------------------------------------------------------------

class BamTag:
    # Cell barcode (corrected)
    CELL_BARCODE        = "CB"
    # Cell barcode (raw, uncorrected)
    CELL_BARCODE_RAW    = "CR"
    # Cell barcode quality
    CELL_BARCODE_QUAL   = "CY"
    # UMI (corrected)
    UMI                 = "UB"
    # UMI (raw)
    UMI_RAW             = "UR"
    # Number of alignments for this read
    NH                  = "NH"
    # Supplementary alignment loci
    SA                  = "SA"
    # Read category assigned by Cell Ranger (E/N/I)
    CELLRANGER_RE       = "RE"
    # Alignment score
    AS                  = "AS"
    # Mapping quality
    MAPQ                = "MQ"


# ---------------------------------------------------------------------------
# Platform identifiers
# ---------------------------------------------------------------------------

class Platform(str, Enum):
    ONT          = "ont"
    PACBIO       = "pacbio"
    ILLUMINA     = "illumina"      # generic Illumina (auto-selects detection path)
    ILLUMINA_10X = "illumina_10x"
    ILLUMINA_BD  = "illumina_bd"
    SMARTSEQ     = "smartseq"      # Smart-seq2/FLASH-seq/Smart-seq3 (plate-based, no CB)
    UNKNOWN      = "unknown"


# Known aligner @PG ID strings → platform mapping
# Used for auto-detection from BAM header
ALIGNER_PLATFORM_MAP = {
    # An aligner does not identify the sequencing platform.  Minimap2 is used
    # for ONT and PacBio; bare STAR is used for many bulk and single-cell
    # protocols.  Pipeline-specific @PG hints refine these UNKNOWN/generic
    # values in bam_inspector.py.
    "minimap2":   Platform.UNKNOWN,
    "pbmm2":      Platform.PACBIO,
    "STAR":       Platform.ILLUMINA,
    "STARsolo":   Platform.ILLUMINA_10X,
    "cellranger": Platform.ILLUMINA_10X,
    "bwa":        Platform.UNKNOWN,
}


# ---------------------------------------------------------------------------
# Pipeline stage identifiers
# ---------------------------------------------------------------------------

class PipelineStage(str, Enum):
    RAW          = "raw"
    PRE_FILTER   = "pre_filter"
    POST_FILTER  = "post_filter"
    CUSTOM       = "custom"


# ---------------------------------------------------------------------------
# Chemistry / barcode tag conventions
# ---------------------------------------------------------------------------

class Chemistry(str, Enum):
    TENX_V3          = "10x_v3"
    TENX_V4          = "10x_v4"
    BD_RHAPSODY_WTA  = "bd_rhapsody_wta"
    CUSTOM           = "custom"

# Barcode length (bp) per chemistry
CHEMISTRY_BARCODE_LENGTH = {
    Chemistry.TENX_V3:         16,
    Chemistry.TENX_V4:         16,
    Chemistry.BD_RHAPSODY_WTA: 27,   # 9+9+9 concatenated cell label
    Chemistry.CUSTOM:          None, # inferred from data
}


# ---------------------------------------------------------------------------
# TSO sequences (used for soft-clip TSO invasion detection)
# ---------------------------------------------------------------------------

# 10x Genomics v3/v4 TSO (same sequence for both)
TSO_10X = "AAGCAGTGGTATCAACGCAGAGTACATGGG"

# PacBio Kinnex / IsoSeq TSO (5′ SMART primer, also used in SMART-seq)
TSO_PACBIO = "AAGCAGTGGTATCAACGCAGAGT"

# Minimum matching length for TSO detection in soft-clipped bases
TSO_MIN_MATCH_LENGTH = 12

# Poly-G tail minimum length to flag as TSO-proximal
TSO_POLYG_MIN_LENGTH = 6


# ---------------------------------------------------------------------------
# Chimeric read detection thresholds
# ---------------------------------------------------------------------------

# Deprecated compatibility value. Distance alone is never sufficient for a
# chimera call because genomic distance includes introns; version 0.7 does not
# use this value in the decision rule.
DEFAULT_CHIMERIC_DISTANCE = 10_000  # bp

# For Illumina paired-end: maximum intra-chromosomal insert size (bp) before
# a read pair is called chimeric.  Illumina fragments are typically <1 kb;
# 1 Mb is a generous upper bound that catches structural variants / chimeras.
ILLUMINA_CHIMERIC_INSERT_SIZE = 1_000_000  # bp


# ---------------------------------------------------------------------------
# Intergenic adaptive threshold parameters
# ---------------------------------------------------------------------------

# Minimum number of distinct barcodes at an intergenic locus to consider it
# a "candidate" rather than sparse noise.  Scales with total barcode count.
ADAPTIVE_MIN_BARCODES_FRACTION = 0.0001   # 0.01 % of total detected barcodes
ADAPTIVE_MIN_BARCODES_ABSOLUTE = 3        # hard floor regardless of fraction

# Minimum reads at locus for Poisson significance testing
ADAPTIVE_MIN_READS = 5

# Multiple-testing-adjusted enrichment threshold for fixed intergenic windows
ADAPTIVE_PVALUE_THRESHOLD = 0.01

# Window size (bp) for locus read-depth aggregation
INTERGENIC_LOCUS_WINDOW = 500

# Require this fraction of the sampled aligned bases at an intergenic window
# to overlap RepeatMasker before assigning INTERGENIC_REPEAT.
INTERGENIC_REPEAT_MIN_FRACTION = 0.50


# ---------------------------------------------------------------------------
# Internal priming / hotspot detection
# ---------------------------------------------------------------------------

# Length of A-run immediately downstream of read 3′ end in reference
# required to flag as internal-priming candidate
POLYA_RUN_MIN_LENGTH = 6
POLYA_CONTEXT_WINDOW = 20   # bp of reference to inspect downstream of read end

# Maximum distance (bp) from annotated polyA site for a read to be considered
# "polyA-site proximal" (used to distinguish novel genes from hotspots)
POLYA_SITE_PROXIMITY = 50   # bp


# ---------------------------------------------------------------------------
# Splice site canonicality
# ---------------------------------------------------------------------------

# Canonical donor–acceptor dinucleotide pairs (GT-AG, GC-AG, AT-AC)
CANONICAL_SPLICE_SITES = {
    ("GT", "AG"),
    ("GC", "AG"),
    ("AT", "AC"),
}


# ---------------------------------------------------------------------------
# TSS proximity threshold for 5'-anchored endpoint metric
TSS_SITE_PROXIMITY = 100   # bp — CAGE peak half-width

# Cluster metadata column names (for --obs-metadata TSV)
OBS_BARCODE_COL = "cell_barcode"
OBS_CLUSTER_COL = "cluster"

# Mitochondrial contig names (handles chrM and MT naming conventions)
# ---------------------------------------------------------------------------

MITO_CONTIG_NAMES = {"chrM", "MT", "chrMT", "mitochondrion"}


# ---------------------------------------------------------------------------
# Auto-detection: fraction of reads that must carry CB tag to activate
# barcode-aware mode
# ---------------------------------------------------------------------------

BARCODE_AUTODETECT_SAMPLE_SIZE = 10_000
BARCODE_AUTODETECT_MIN_FRACTION = 0.50


# ---------------------------------------------------------------------------
# Reservoir sampling
# ---------------------------------------------------------------------------

# Seeded by default so two runs over the same BAM agree.  Several reported
# values are computed from reservoir samples, so an unseeded run moved
# broad_noncanonical_read_frac in the third decimal between identical
# invocations, which is larger than many of the differences being compared.
# Pass a different --seed to draw a different sample.
DEFAULT_SEED = 42


# ---------------------------------------------------------------------------
# `compare`: nested-pair detection
# ---------------------------------------------------------------------------

# Retention and transition tables are only meaningful when BAM B is a filtered
# subset of BAM A.  Read names are sampled from a shared genomic window of both
# BAMs; a nested pair shares nearly all of them, two independent samples share
# none.  Anything below the threshold is treated as independent, which also
# avoids materialising a per-read assignment dict for both BAMs at once.
# Below-threshold intergenic windows plotted in the report scatter.  Real
# samples produce 20k-130k windows, the vast majority never promoted; drawing
# all of them made the HTML report tens of MB without adding information.
INTERGENIC_SCATTER_MAX_SPARSE = 3_000


COMPARE_PROBE_SAMPLE_SIZE = 50_000
COMPARE_NESTED_MIN_OVERLAP = 0.50


# ---------------------------------------------------------------------------
# Parallelization
# ---------------------------------------------------------------------------

DEFAULT_THREADS = 4

# MAPQ scales differ between aligners. This threshold is therefore used only
# for a descriptive low-confidence flag, never as a classification filter.
LOW_MAPQ_THRESHOLD = 10


# ---------------------------------------------------------------------------
# Read-length stratification bins
# ---------------------------------------------------------------------------

# Break-points (bp) used with bisect_right to assign reads to bins:
#   index 0 → <150, 1 → 150–500, 2 → 500–1000,
#   index 3 → 1000–2000, 4 → 2000–5000, 5 → ≥5000
LENGTH_BIN_BREAKS = [150, 500, 1000, 2000, 5000]

# Labels when the short-read bin is active (median read length < 300 bp)
LENGTH_BIN_LABELS_LONG  = ["<150", "150–500", "500–1000", "1000–2000", "2000–5000", ">5000"]

# Standard labels (first two bins merged into "<500")
LENGTH_BIN_LABELS_SHORT = ["<500", "500–1000", "1000–2000", "2000–5000", ">5000"]

# Median read-length threshold below which the <150 bp bin is prepended
LENGTH_SHORT_READ_THRESHOLD = 300  # bp
