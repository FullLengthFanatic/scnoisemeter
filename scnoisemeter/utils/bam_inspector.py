"""
bam_inspector.py
================
Inspects a BAM header and a sample of reads to auto-detect:

  - Sequencing platform  (ONT / PacBio / Illumina)
  - Aligner and version  (minimap2 / STAR / STARsolo / pbmm2 / cellranger)
  - Barcode tag names    (CB/UB or custom)
  - Whether CB/UB tags are present and corrected
  - Pipeline stage hints (wf-single-cell, Cell Ranger, isoseq3, …)

All detected values are returned as a BamMetadata dataclass that is passed
into every downstream module.  No analysis logic lives here — this module
only reads and reports.

Design note
-----------
We intentionally avoid crashing on missing or malformed headers.  Real-world
BAM files from research pipelines are messy.  Every field has a sensible
fallback.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pysam

from scnoisemeter.constants import (
    ALIGNER_PLATFORM_MAP,
    BARCODE_AUTODETECT_MIN_FRACTION,
    BARCODE_AUTODETECT_SAMPLE_SIZE,
    BamTag,
    Chemistry,
    CHEMISTRY_BARCODE_LENGTH,
    Platform,
    PipelineStage,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class BamMetadata:
    """
    All auto-detected properties of a BAM file.

    Fields
    ------
    path : Path
        Absolute path to the BAM file.
    platform : Platform
        Detected or user-supplied sequencing platform.
    platform_confidence : str
        "user" | "header" | "inferred" | "unknown"
    aligner : str
        @PG ID string of the primary aligner (e.g. "minimap2").
    aligner_version : str
        @PG VN string of the primary aligner, or "" if absent.
    pipeline_hints : list[str]
        @PG ID strings of all programs recorded in the header, in order.
        Useful for detecting upstream tools (wf-single-cell, isoseq3, …).
    barcode_tag : str
        BAM tag name for the corrected cell barcode (default "CB").
    umi_tag : str
        BAM tag name for the corrected UMI (default "UB").
    barcode_tag_present : bool
        Whether the barcode_tag was found in the sampled reads.
    barcode_fraction : float
        Fraction of sampled reads that carry the barcode_tag.
    barcode_aware : bool
        True if barcode_fraction >= BARCODE_AUTODETECT_MIN_FRACTION.
    pipeline_stage : PipelineStage
        Declared or inferred processing stage.
    reference_names : list[str]
        Contig names from the BAM header (used for chromosome parallelisation).
    reference_lengths : dict[str, int]
        Mapping of contig name → length.
    warnings : list[str]
        Non-fatal issues detected during inspection (surfaced in reports).
    sort_order : str
        Value of the @HD SO tag (e.g. "coordinate", "queryname", "unsorted").
        Empty string if the tag is absent.
    """

    path: Path
    platform: Platform = Platform.UNKNOWN
    platform_confidence: str = "unknown"
    aligner: str = ""
    aligner_version: str = ""
    pipeline_hints: list[str] = field(default_factory=list)
    barcode_tag: str = BamTag.CELL_BARCODE
    umi_tag: str = BamTag.UMI
    barcode_tag_present: bool = False
    barcode_fraction: float = 0.0
    barcode_aware: bool = False
    pipeline_stage: PipelineStage = PipelineStage.CUSTOM
    reference_names: list[str] = field(default_factory=list)
    reference_lengths: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    sort_order: str = ""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def inspect_bam(
    bam_path: str | Path,
    *,
    barcode_tag: str = BamTag.CELL_BARCODE,
    umi_tag: str = BamTag.UMI,
    platform: Optional[Platform] = None,
    pipeline_stage: Optional[PipelineStage] = None,
    sample_size: int = BARCODE_AUTODETECT_SAMPLE_SIZE,
) -> BamMetadata:
    """
    Open *bam_path* (must have a .bai index), inspect the header and a
    sample of reads, and return a fully populated :class:`BamMetadata`.

    Parameters
    ----------
    bam_path:
        Path to the BAM file.
    barcode_tag:
        BAM tag to treat as the corrected cell barcode.  Overrides auto-
        detection; useful for BD Rhapsody or custom pipelines.
    umi_tag:
        BAM tag to treat as the corrected UMI.
    platform:
        If provided by the user, skip platform auto-detection.
    pipeline_stage:
        If provided by the user, skip stage auto-detection.
    sample_size:
        Number of reads to sample for barcode-tag presence check.
    """
    bam_path = Path(bam_path)
    meta = BamMetadata(path=bam_path, barcode_tag=barcode_tag, umi_tag=umi_tag)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        _parse_header(bam, meta)
        _detect_barcode_tags(bam, meta, sample_size)

    # User overrides win over auto-detection
    if platform is not None:
        meta.platform = platform
        meta.platform_confidence = "user"

    if pipeline_stage is not None:
        meta.pipeline_stage = pipeline_stage

    _log_summary(meta)
    return meta


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _parse_header(bam: pysam.AlignmentFile, meta: BamMetadata) -> None:
    """Extract platform, aligner, and pipeline hints from the BAM @PG lines."""
    header = bam.header.to_dict()

    # Sort order from @HD
    hd = header.get("HD", {})
    meta.sort_order = hd.get("SO", "")

    # Reference contigs
    sq = header.get("SQ", [])
    meta.reference_names = [s["SN"] for s in sq]
    meta.reference_lengths = {s["SN"]: s["LN"] for s in sq}

    if not meta.reference_names:
        meta.warnings.append("BAM header contains no @SQ lines — cannot build contig list.")

    # Program records
    pg_records = header.get("PG", [])
    # Collect ID and PN fields as hints, plus scan CL fields for known pipeline
    # tool names embedded in file paths (e.g. cellranger writes no own @PG record
    # but leaves its name in the CL paths of downstream samtools invocations).
    _known_pipeline_names = {
        "cellranger", "spaceranger", "starsolo", "star",
        "wf-single-cell", "guppy", "dorado", "bonito",
        "isoseq3", "pbmm2", "skera", "lima",
    }
    _hints: list[str] = []
    for p in pg_records:
        _hints.append(p.get("ID", ""))
        if p.get("PN"):
            _hints.append(p["PN"])
        cl = p.get("CL", "")
        if cl:
            cl_lower = cl.lower()
            for name in _known_pipeline_names:
                if name in cl_lower and name not in [h.lower() for h in _hints]:
                    _hints.append(name)
    meta.pipeline_hints = _hints

    # Identify the primary aligner (first @PG record whose ID matches a known aligner)
    for pg in pg_records:
        pg_id = pg.get("ID", "")
        pg_vn = pg.get("VN", "")
        for known_aligner, detected_platform in ALIGNER_PLATFORM_MAP.items():
            if known_aligner.lower() in pg_id.lower():
                meta.aligner = pg_id
                meta.aligner_version = pg_vn
                if meta.platform == Platform.UNKNOWN:
                    meta.platform = detected_platform
                    meta.platform_confidence = "header"
                break
        if meta.aligner:
            break

    # Refine platform from well-known pipeline tool names in @PG
    # (may set meta.aligner via CL-path detection, e.g. cellranger)
    _refine_platform_from_pipeline_hints(meta)

    # Only warn if aligner is still unresolved after all detection paths
    if not meta.aligner:
        meta.warnings.append(
            "Could not identify a known aligner in @PG records.  "
            "Platform auto-detection will be unreliable."
        )

    # Pipeline stage hints
    _infer_pipeline_stage(meta)


def _refine_platform_from_pipeline_hints(meta: BamMetadata) -> None:
    """
    Use @PG program names to refine platform beyond what the aligner alone tells us.

    Examples:
      - "wf-single-cell" in hints → ONT
      - "isoseq3" or "pbmm2" in hints → PacBio
      - "cellranger" in hints → Illumina 10x
    """
    hints_lower = [h.lower() for h in meta.pipeline_hints]

    ont_signals   = {"wf-single-cell", "guppy", "dorado", "bonito", "ont"}
    pacbio_signals = {"isoseq3", "pbmm2", "skera", "lima", "ccs", "kinnex"}
    illumina_signals = {"cellranger", "starsolo", "star", "spaceranger"}

    for signal in ont_signals:
        if any(signal in h for h in hints_lower):
            meta.platform = Platform.ONT
            meta.platform_confidence = "header"
            return

    for signal in pacbio_signals:
        if any(signal in h for h in hints_lower):
            meta.platform = Platform.PACBIO
            meta.platform_confidence = "header"
            return

    for signal in illumina_signals:
        if any(signal in h for h in hints_lower):
            # Distinguish 10x vs BD by checking for BD-specific keywords
            bd_signals = {"rhapsody", "bd", "seven-bridges"}
            if any(s in h for h in hints_lower for s in bd_signals):
                meta.platform = Platform.ILLUMINA_BD
            else:
                meta.platform = Platform.ILLUMINA_10X
            meta.platform_confidence = "header"
            # If aligner was not identified from @PG ID/PN fields (e.g. Cell Ranger
            # leaves no own @PG record but its name appears in downstream CL paths),
            # use the matched signal as the aligner name so the report says
            # "cellranger" rather than "unknown".
            if not meta.aligner:
                meta.aligner = signal
            return


def _infer_pipeline_stage(meta: BamMetadata) -> None:
    """
    Attempt to infer the processing stage from @PG hints.

    Heuristics:
      - "wf-single-cell" present → post_filter (wf-single-cell output BAM)
      - "cellranger" present → post_filter
      - "isoseq3" + "refine" present → pre_filter (FLNC, pre-cell-assignment)
      - Only aligner present → raw or pre_filter (ambiguous; report as pre_filter)
    """
    hints_lower = [h.lower() for h in meta.pipeline_hints]

    post_filter_signals = {"wf-single-cell", "cellranger", "starsolo"}
    for signal in post_filter_signals:
        if any(signal in h for h in hints_lower):
            meta.pipeline_stage = PipelineStage.POST_FILTER
            return

    # isoseq3 refine → FLNC, pre cell assignment
    if any("isoseq3" in h for h in hints_lower):
        if any("refine" in h for h in hints_lower):
            meta.pipeline_stage = PipelineStage.PRE_FILTER
        else:
            meta.pipeline_stage = PipelineStage.RAW
        return

    # Only an aligner — raw or very early pre-filter
    if meta.aligner:
        meta.pipeline_stage = PipelineStage.PRE_FILTER
        return

    meta.pipeline_stage = PipelineStage.CUSTOM


def _detect_barcode_tags(
    bam: pysam.AlignmentFile,
    meta: BamMetadata,
    sample_size: int,
) -> None:
    """
    Sample the first *sample_size* primary, mapped reads and check what
    fraction carry the requested barcode tag.

    Also detects whether only the raw barcode tag (CR) is present and the
    corrected tag (CB) is absent — a common situation in pre-correction BAMs.
    """
    n_sampled = 0
    n_with_tag = 0
    n_with_raw_only = 0

    for read in bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if n_sampled >= sample_size:
            break

        n_sampled += 1
        has_corrected = read.has_tag(meta.barcode_tag)
        has_raw = read.has_tag(BamTag.CELL_BARCODE_RAW)

        if has_corrected:
            n_with_tag += 1
        elif has_raw and not has_corrected:
            n_with_raw_only += 1

    if n_sampled == 0:
        meta.warnings.append("No primary mapped reads found during barcode tag sampling.")
        return

    meta.barcode_fraction = n_with_tag / n_sampled
    meta.barcode_tag_present = meta.barcode_fraction >= BARCODE_AUTODETECT_MIN_FRACTION
    meta.barcode_aware = meta.barcode_tag_present

    raw_fraction = n_with_raw_only / n_sampled
    if raw_fraction > 0.1 and not meta.barcode_tag_present:
        meta.warnings.append(
            f"{raw_fraction:.1%} of sampled reads carry the raw barcode tag "
            f"({BamTag.CELL_BARCODE_RAW}) but not the corrected tag "
            f"({meta.barcode_tag}).  This BAM may be pre-barcode-correction.  "
            f"Per-cell noise metrics will be unavailable."
        )

    if not meta.barcode_tag_present:
        meta.warnings.append(
            f"Corrected barcode tag '{meta.barcode_tag}' found in only "
            f"{meta.barcode_fraction:.1%} of sampled reads "
            f"(threshold: {BARCODE_AUTODETECT_MIN_FRACTION:.0%}).  "
            f"Running in barcode-agnostic mode: all reads will be classified "
            f"genomically and aggregated as a single sample. "
            f"Per-cell metrics will not be available.  "
            f"This is expected for PacBio FLTNC BAMs without barcode demultiplexing, "
            f"or for any BAM produced before the barcode assignment step."
        )

    logger.debug(
        "Barcode tag detection: %d/%d reads carry '%s' (%.1f%%)",
        n_with_tag, n_sampled, meta.barcode_tag, meta.barcode_fraction * 100,
    )


def _log_summary(meta: BamMetadata) -> None:
    logger.info("BAM inspection complete: %s", meta.path.name)
    logger.info("  Platform       : %s (%s)", meta.platform.value, meta.platform_confidence)
    logger.info("  Aligner        : %s %s", meta.aligner or "unknown", meta.aligner_version)
    logger.info("  Pipeline stage : %s", meta.pipeline_stage.value)
    logger.info("  Pipeline hints : %s", ", ".join(meta.pipeline_hints) or "none")
    logger.info("  Barcode-aware  : %s (%.1f%% reads tagged)",
                meta.barcode_aware, meta.barcode_fraction * 100)
    for w in meta.warnings:
        logger.warning("  WARNING: %s", w)
