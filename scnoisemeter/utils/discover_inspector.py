"""
discover_inspector.py
=====================
Extended BAM inspection for the ``discover`` subcommand.

Extends the core inspect_bam analysis with:
  - Chemistry detection from @CO / @PG / @RG fields
  - Approximate read-count from BAM index (fast, no scan)
  - Paired-end detection from read flags
  - .bai presence check
  - Formatted summary table for the CLI

None of the code here touches the network or the filesystem beyond the BAM
being inspected.  All real analysis stays in the core pipeline.
"""

from __future__ import annotations

import logging
import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pysam

from scnoisemeter.constants import (
    BARCODE_AUTODETECT_MIN_FRACTION,
    Chemistry,
    Platform,
    PipelineStage,
)
from scnoisemeter.utils.bam_inspector import BamMetadata, inspect_bam

logger = logging.getLogger(__name__)

# Max reads to sample for paired-end / chemistry detection (fast)
_SAMPLE_READS = 1_000


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class DiscoverBamInfo:
    """
    All auto-detected properties of one BAM file, for use by ``discover``.
    Extends BamMetadata with chemistry, paired-end, and approximate count.
    """
    meta: BamMetadata
    chemistry: str = "unknown"           # Chemistry value or "unknown"
    chemistry_confidence: str = "unknown"  # "header" | "rg" | "unknown"
    paired_end: bool = False
    n_reads_approx: int = 0
    has_index: bool = False
    run_issues: list[str] = field(default_factory=list)  # blocking issues

    @property
    def bam_path(self) -> Path:
        return self.meta.path

    @property
    def can_run(self) -> bool:
        """True if no blocking issues prevent running this BAM."""
        return len(self.run_issues) == 0

    @property
    def platform(self) -> Platform:
        return self.meta.platform

    @property
    def pipeline_stage(self) -> PipelineStage:
        return self.meta.pipeline_stage

    @property
    def sort_order(self) -> str:
        return self.meta.sort_order


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def inspect_bam_for_discover(bam_path: Path) -> DiscoverBamInfo:
    """
    Inspect *bam_path* and return a :class:`DiscoverBamInfo`.

    Calls :func:`inspect_bam` for core detection, then adds chemistry,
    paired-end, read count, and index checks.
    """
    bam_path = Path(bam_path).resolve()

    # Check for index before calling inspect_bam (which needs the index)
    has_index = (
        bam_path.with_suffix(".bam.bai").exists()
        or Path(str(bam_path) + ".bai").exists()
    )

    info = DiscoverBamInfo(
        meta=BamMetadata(path=bam_path),
        has_index=has_index,
    )

    if not has_index:
        info.run_issues.append("no .bai index found — run 'samtools index' first")
        return info

    # Use inspect_bam for core detection
    try:
        info.meta = inspect_bam(bam_path)
    except Exception as exc:
        info.run_issues.append(f"BAM inspection failed: {exc}")
        return info

    # Override pipeline_stage from barcode fraction (task spec: CB >50% → post_filter)
    # bam_inspector's _infer_pipeline_stage uses only @PG hints, which may not capture
    # all post-filter scenarios (e.g. EPI2ME ONT BAMs without wf-single-cell @PG entry).
    if info.meta.barcode_aware and info.meta.barcode_fraction >= BARCODE_AUTODETECT_MIN_FRACTION:
        info.meta.pipeline_stage = PipelineStage.POST_FILTER
    elif not info.meta.barcode_aware:
        if info.meta.pipeline_stage not in (PipelineStage.PRE_FILTER, PipelineStage.RAW):
            info.meta.pipeline_stage = PipelineStage.PRE_FILTER

    # Additional discovery
    _detect_chemistry(bam_path, info)
    _detect_paired_end(bam_path, info)
    _get_approx_read_count(bam_path, info)

    # Blocking issues
    if info.meta.sort_order and info.meta.sort_order.lower() not in ("coordinate", ""):
        info.run_issues.append(
            f"sort order is '{info.meta.sort_order}' — must be coordinate-sorted"
        )
    if info.meta.platform == Platform.UNKNOWN:
        info.run_issues.append("platform could not be inferred from header")

    return info


# ---------------------------------------------------------------------------
# Chemistry detection
# ---------------------------------------------------------------------------

def _detect_chemistry(bam_path: Path, info: DiscoverBamInfo) -> None:
    """
    Attempt to detect library chemistry from @CO, @PG, and @RG fields.

    Detection priority:
      1. @CO or @PG CL field containing "10x" / "chromium" → version sub-classify
      2. @RG DS field
      3. @CO or @PG containing "BD Rhapsody" / "FLASH"
    """
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            header = bam.header.to_dict()
    except Exception:
        return

    # Gather all text blobs to search
    blobs: list[str] = []
    for co in header.get("CO", []):
        blobs.append(co if isinstance(co, str) else str(co))
    for pg in header.get("PG", []):
        blobs.append(pg.get("CL", ""))
        blobs.append(pg.get("PN", ""))
        blobs.append(pg.get("ID", ""))

    combined = " ".join(blobs).lower()

    # BD Rhapsody / FLASH-seq
    # Use specific tokens to avoid false positives: "bd" matches hex UUIDs,
    # "flash" can appear in unrelated filenames.
    _bd_tokens = ("rhapsody", "flashbd", "flash-seq", "flash_seq", "bd rhapsody")
    if any(tok in combined for tok in _bd_tokens):
        info.chemistry = Chemistry.BD_RHAPSODY_WTA.value
        info.chemistry_confidence = "header"
        return

    # 10x / Chromium
    if "10x" in combined or "chromium" in combined or "gem-x" in combined:
        if "v4" in combined or "gem-x" in combined or "4.0" in combined:
            info.chemistry = Chemistry.TENX_V4.value
        elif "v3" in combined or "3.1" in combined or "3.0" in combined:
            info.chemistry = Chemistry.TENX_V3.value
        else:
            # Default 10x to v3 if version not resolvable
            info.chemistry = Chemistry.TENX_V3.value
        info.chemistry_confidence = "header"
        return

    # @RG DS field
    for rg in header.get("RG", []):
        ds = rg.get("DS", "")
        ds_lower = ds.lower()
        if "10x" in ds_lower or "chromium" in ds_lower:
            info.chemistry = Chemistry.TENX_V3.value
            info.chemistry_confidence = "rg"
            return
        if "rhapsody" in ds_lower or "bd" in ds_lower:
            info.chemistry = Chemistry.BD_RHAPSODY_WTA.value
            info.chemistry_confidence = "rg"
            return

    # If platform is illumina/pacbio/ont but chemistry unknown, leave as unknown
    info.chemistry = "unknown"
    info.chemistry_confidence = "unknown"


# ---------------------------------------------------------------------------
# Paired-end detection
# ---------------------------------------------------------------------------

def _detect_paired_end(bam_path: Path, info: DiscoverBamInfo) -> None:
    """Set info.paired_end = True if >50% of sampled reads have FLAG 0x1."""
    try:
        n_paired = 0
        n_total = 0
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                if n_total >= _SAMPLE_READS:
                    break
                n_total += 1
                if read.is_paired:
                    n_paired += 1
        if n_total > 0:
            info.paired_end = (n_paired / n_total) > 0.5
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Approximate read count from index
# ---------------------------------------------------------------------------

def _get_approx_read_count(bam_path: Path, info: DiscoverBamInfo) -> None:
    """Estimate total mapped reads from the BAM index using samtools idxstats."""
    try:
        result = subprocess.run(
            ["samtools", "idxstats", str(bam_path)],
            capture_output=True, text=True, timeout=30,
        )
        total = 0
        for line in result.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) >= 3:
                try:
                    total += int(parts[2])  # column 3: mapped reads
                except ValueError:
                    pass
        info.n_reads_approx = total
    except Exception:
        info.n_reads_approx = 0


# ---------------------------------------------------------------------------
# Platform normalisation
# ---------------------------------------------------------------------------

def _normalise_platform(platform: Platform) -> str:
    """
    Return the canonical platform string for discover display and --platform flag.

    illumina_10x and illumina_bd are both presented as "illumina" so that the
    discover table and the constructed run command are consistent with the
    --platform flag accepted by ``scnoisemeter run``.
    """
    if platform in (Platform.ILLUMINA_10X, Platform.ILLUMINA_BD):
        return Platform.ILLUMINA.value
    return platform.value


# ---------------------------------------------------------------------------
# Summary table formatting
# ---------------------------------------------------------------------------

def format_discovery_table(infos: list[DiscoverBamInfo]) -> str:
    """
    Return a formatted ASCII table summarising discovered BAMs.

    Flags problematic fields with ⚠.
    """
    if not infos:
        return "  (no BAM files found)\n"

    rows = []
    for i, info in enumerate(infos, 1):
        name = info.bam_path.name
        plat = _normalise_platform(info.platform) if info.platform != Platform.UNKNOWN else "unknown ⚠"
        chem = info.chemistry if info.chemistry != "unknown" else "unknown"
        stage = info.pipeline_stage.value if info.pipeline_stage != PipelineStage.CUSTOM else "unknown"
        reads = f"~{info.n_reads_approx / 1000:.1f}k" if info.n_reads_approx else "?"
        paired = "yes" if info.paired_end else "no"
        sort = info.sort_order or "unknown"
        if sort not in ("coordinate", ""):
            sort += " ⚠"
        if not info.has_index:
            sort = "— ⚠ no index"
        rows.append((str(i), name, plat, chem, stage, reads, paired, sort))

    # Column widths
    headers = ("#", "File", "Platform", "Chemistry", "Stage", "Reads", "Paired", "Sort")
    widths = [max(len(headers[j]), max(len(r[j]) for r in rows)) for j in range(len(headers))]

    def fmt_row(r):
        return "  " + "  ".join(f"{r[j]:<{widths[j]}}" for j in range(len(headers)))

    sep = "  " + "  ".join("-" * w for w in widths)
    lines = [fmt_row(headers), sep]
    for r in rows:
        lines.append(fmt_row(r))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Command-line prompt helpers (used by discover_cmd, factored here for testing)
# ---------------------------------------------------------------------------

def _collect_selected_indices(infos: list[DiscoverBamInfo]) -> list[int]:
    """
    Interactive prompt: ask which BAMs to run.
    Returns list of 0-based indices.
    """
    import sys
    while True:
        sys.stdout.write(
            "\nEnter the numbers of the BAMs you want to run (comma-separated),\n"
            "or 'all', or 'quit' to exit: "
        )
        sys.stdout.flush()
        raw = sys.stdin.readline().strip().lower()
        if raw == "quit":
            return []
        if raw == "all":
            return list(range(len(infos)))
        try:
            indices = []
            for token in raw.replace(" ", "").split(","):
                n = int(token)
                if 1 <= n <= len(infos):
                    indices.append(n - 1)
                else:
                    print(f"  Number {n} is out of range (1–{len(infos)}).")
                    break
            else:
                return indices
        except ValueError:
            print("  Invalid input. Enter numbers, 'all', or 'quit'.")


def _prompt_platform(info: DiscoverBamInfo) -> Optional[str]:
    """
    Ask the user to supply a platform for a BAM with unknown platform.
    Returns platform value string, or None to skip.
    """
    import sys
    name = info.bam_path.name
    print(f"\nBAM {name}: platform could not be inferred.")
    choices = [p.value for p in (Platform.ONT, Platform.PACBIO, Platform.ILLUMINA)]
    for i, c in enumerate(choices, 1):
        print(f"  {i}) {c}")
    print(f"  {len(choices)+1}) skip")
    while True:
        sys.stdout.write("Enter choice: ")
        sys.stdout.flush()
        raw = sys.stdin.readline().strip()
        try:
            n = int(raw)
            if n == len(choices) + 1:
                return None
            if 1 <= n <= len(choices):
                return choices[n - 1]
        except ValueError:
            pass
        print(f"  Enter 1–{len(choices)+1}.")


def _prompt_chemistry(info: DiscoverBamInfo) -> str:
    """
    Ask the user to supply chemistry for a BAM with unknown chemistry.
    Returns chemistry value string ("unknown" if user wants barcode-agnostic).
    """
    import sys
    name = info.bam_path.name
    print(f"\nBAM {name}: chemistry could not be inferred.")
    choices = [c.value for c in (Chemistry.TENX_V3, Chemistry.TENX_V4, Chemistry.BD_RHAPSODY_WTA)]
    for i, c in enumerate(choices, 1):
        print(f"  {i}) {c}")
    print(f"  {len(choices)+1}) unknown (run barcode-agnostic)")
    while True:
        sys.stdout.write("Enter choice: ")
        sys.stdout.flush()
        raw = sys.stdin.readline().strip()
        try:
            n = int(raw)
            if n == len(choices) + 1:
                return "unknown"
            if 1 <= n <= len(choices):
                return choices[n - 1]
        except ValueError:
            pass
        print(f"  Enter 1–{len(choices)+1}.")
