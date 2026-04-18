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

from enum import Enum, auto


# ---------------------------------------------------------------------------
# Read classification taxonomy
# ---------------------------------------------------------------------------

class ReadCategory(str, Enum):
    """
    Exhaustive, mutually-exclusive classification of every primary alignment.

    Hierarchy (earlier entries take priority over later ones during assignment):

      UNMAPPED          – read did not align; excluded from all fractions
      SECONDARY         – SAM flag 0x100; skipped (redundant record)
      SUPPLEMENTARY     – SAM flag 0x800; routed to chimeric detector only
      MULTIMAPPER       – primary alignment with NH > 1
      CHIMERIC          – SA tag present AND (inter-chromosomal OR
                          strand-discordant OR same-strand distance > threshold)
      MITOCHONDRIAL     – maps to the mitochondrial contig (chrM / MT)
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

# Categories counted as "noise" in summary noise-fraction metrics.
# Intronic sub-categories are reported separately; users choose interpretation.
NOISE_CATEGORIES = {
    ReadCategory.EXONIC_ANTISENSE,
    ReadCategory.INTRONIC_PURE,
    ReadCategory.INTRONIC_BOUNDARY,
    ReadCategory.INTERGENIC_SPARSE,
    ReadCategory.INTERGENIC_REPEAT,
    ReadCategory.INTERGENIC_HOTSPOT,
    ReadCategory.CHIMERIC,
}

# Conservative noise: includes INTRONIC_PURE and INTRONIC_BOUNDARY, which may
# represent genuine pre-mRNA capture rather than RT/PCR artifacts. This is the
# reported "Noise" fraction and represents an upper bound on true noise.
NOISE_CATEGORIES_CONSERVATIVE = NOISE_CATEGORIES  # alias for clarity

# Strict noise: only unambiguous RT/PCR/sequencing artifacts. Excludes
# INTRONIC_PURE and INTRONIC_BOUNDARY because these cannot be distinguished
# from genuine pre-mRNA capture at the read level. This is a lower bound on
# true noise — everything in this set is definitively an artifact.
NOISE_CATEGORIES_STRICT = {
    ReadCategory.EXONIC_ANTISENSE,
    ReadCategory.INTERGENIC_SPARSE,
    ReadCategory.INTERGENIC_REPEAT,
    ReadCategory.INTERGENIC_HOTSPOT,
    ReadCategory.CHIMERIC,
}

# Unstranded-protocol noise: excludes EXONIC_ANTISENSE, which is genuine
# cDNA signal in non-stranded libraries (Smart-seq2, FLASH-seq, Smart-seq3).
# Use this when platform == SMARTSEQ to avoid inflating the noise fraction.
NOISE_CATEGORIES_UNSTRANDED = NOISE_CATEGORIES - {ReadCategory.EXONIC_ANTISENSE}

# Categories whose interpretation is ambiguous and reported separately
AMBIGUOUS_CATEGORIES = {
    ReadCategory.INTRONIC_JXNSPAN,
    ReadCategory.INTERGENIC_NOVEL,
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
    "minimap2":   Platform.ONT,       # default assumption for minimap2
    "pbmm2":      Platform.PACBIO,
    "STAR":       Platform.ILLUMINA_10X,
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

# PacBio Kinnex adapter sequences (used by skera for segmentation)
# These appear in chimeric reads where segmentation was imperfect
KINNEX_ADAPTER = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

# Minimum matching length for TSO detection in soft-clipped bases
TSO_MIN_MATCH_LENGTH = 12

# Poly-G tail minimum length to flag as TSO-proximal
TSO_POLYG_MIN_LENGTH = 6


# ---------------------------------------------------------------------------
# Chimeric read detection thresholds
# ---------------------------------------------------------------------------

# Default maximum intra-chromosomal same-strand distance (bp) below which a
# split alignment is considered a legitimate splice, not a chimera.
# Rationale: polyA+ cDNA molecules are typically <3 kb; RT processivity
# limits capture to ~10–15 kb even for long transcripts.  10 kb is generous.
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

# Bonferroni-corrected p-value threshold for intergenic locus significance
ADAPTIVE_PVALUE_THRESHOLD = 0.01

# Window size (bp) for locus read-depth aggregation
INTERGENIC_LOCUS_WINDOW = 500


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
# TSS proximity threshold for 5'-anchored full-length metric
TSS_SITE_PROXIMITY = 100   # bp — CAGE peak half-width

# NUMT: nuclear mitochondrial DNA segments
# Reads mapping to chrM that also have a good match to a NUMT locus
# are ambiguous; scNoiseMeter flags them separately when a NUMT BED is provided.
NUMT_OVERLAP_FRACTION = 0.80   # fraction of read that must overlap NUMT

# Full-length read detection
# When a polyA site database is provided (--polya-sites), a read is considered
# full-length if its 3′ end falls within POLYA_SITE_PROXIMITY bp of a known
# polyA site.  Without the database, the fallback is a minimum length threshold.
FULL_LENGTH_MIN_LENGTH_ONT    = 500   # bp — fallback when no polyA site DB
FULL_LENGTH_MIN_LENGTH_PACBIO = 1000  # bp — PacBio HiFi reads are longer

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
# Parallelization
# ---------------------------------------------------------------------------

DEFAULT_THREADS = 4


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
