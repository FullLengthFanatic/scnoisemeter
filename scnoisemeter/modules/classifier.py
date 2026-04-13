"""
classifier.py
=============
The core per-read classification engine.

For every primary, non-secondary alignment the classifier:

  1. Applies the priority hierarchy (chimeric → mito → exonic → intronic
     → intergenic) to assign a single :class:`ReadCategory`.
  2. Simultaneously tallies the number of aligned *bases* falling into each
     category (base-level fractions — the primary benchmarking metric).
  3. Checks for TSO invasion signature in soft-clipped 5′ bases.
  4. Checks for internal oligo-dT priming signature at the read 3′ end.
  5. Checks splice-site canonicality for junction-spanning intronic reads.

Classification hierarchy
------------------------
  UNMAPPED          skip
  SECONDARY         skip
  SUPPLEMENTARY     chimeric detector only
  MULTIMAPPER       (NH > 1) → multimapper bucket, still base-classified
  CHIMERIC          SA tag + distance/strand rules
  MITOCHONDRIAL     maps to chrM / MT
  EXONIC_SENSE      ≥1 bp overlap with exon, correct strand
  EXONIC_ANTISENSE  ≥1 bp overlap with exon, wrong strand
  INTRONIC_JXNSPAN  intronic + CIGAR N at known/novel splice site
  INTRONIC_BOUNDARY spans exon-intron boundary, no splice op
  INTRONIC_PURE     entirely within intron body
  INTERGENIC_*      no gene body overlap (sub-classified later by
                    the IntegergenicProfiler module)
  AMBIGUOUS         overlaps shared region of ≥2 genes
  UNASSIGNED        CB tag absent or not on whitelist

For each read we also compute:
  - Read-level category   (single enum, used for per-cell QC)
  - Base-level breakdown  (dict category → n_bases, used for benchmarking)
  - Flags: is_tso_invasion, is_polya_priming, has_noncanonical_junction

This module is stateless: it exposes a single classify_read() function and
a ReadResult dataclass.  Aggregation is the responsibility of the caller
(the contig worker in pipeline.py).
"""

from __future__ import annotations

import re
from bisect import bisect_left
from dataclasses import dataclass, field
from typing import Optional

import pysam

from scnoisemeter.constants import (
    BamTag,
    CANONICAL_SPLICE_SITES,
    DEFAULT_CHIMERIC_DISTANCE,
    ILLUMINA_CHIMERIC_INSERT_SIZE,
    MITO_CONTIG_NAMES,
    POLYA_CONTEXT_WINDOW,
    POLYA_RUN_MIN_LENGTH,
    ReadCategory,
    SamFlag,
    TSO_10X,
    TSO_MIN_MATCH_LENGTH,
    TSO_PACBIO,
    TSO_POLYG_MIN_LENGTH,
)

# Import new sub-categories (defined in constants but also referenced directly)
_AMBIGUOUS_COD_COD  = ReadCategory.AMBIGUOUS_COD_COD
_AMBIGUOUS_COD_NCOD = ReadCategory.AMBIGUOUS_COD_NCOD
from scnoisemeter.modules.annotation import AnnotationIndex


# ---------------------------------------------------------------------------
# Result dataclass (one per classified read)
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class ReadResult:
    """
    Classification result for a single primary alignment.

    Attributes
    ----------
    query_name : str
    cell_barcode : str
        Corrected CB tag value, or "" if absent / unassigned.
    umi : str
        Corrected UB tag value, or "" if absent.
    category : ReadCategory
        Single read-level category assignment (plurality-vote rule).
    base_counts : dict[ReadCategory, int]
        Number of aligned bases assigned to each category.
        Sums to the number of aligned (non-clipped) bases.
    is_multimapper : bool
        True if NH > 1.
    is_tso_invasion : bool
        True if TSO sequence detected in 5′ soft-clip.
    is_polya_priming : bool
        True if A-rich stretch detected immediately downstream of 3′ end.
    has_noncanonical_junction : bool
        True if any CIGAR N operation falls at a non-canonical splice site
        not present in the annotation.
    contig : str
    pos : int
        0-based leftmost mapping position.
    is_reverse : bool
    read_length : int
        Length of the aligned portion (excluding soft/hard clips).
    """

    query_name:             str
    cell_barcode:           str
    umi:                    str
    category:               ReadCategory
    base_counts:            dict  # ReadCategory → int
    is_multimapper:         bool
    is_tso_invasion:        bool
    is_polya_priming:       bool
    has_noncanonical_junction: bool
    contig:                 str
    pos:                    int
    is_reverse:             bool
    read_length:            int


# ---------------------------------------------------------------------------
# Classifier
# ---------------------------------------------------------------------------

class ReadClassifier:
    """
    Stateless per-read classifier.

    Create one instance per worker process (the annotation index is shared
    read-only after fork).

    Parameters
    ----------
    index : AnnotationIndex
        Pre-built annotation intervals.
    barcode_tag : str
        BAM tag for corrected cell barcode.
    umi_tag : str
        BAM tag for corrected UMI.
    whitelist : set[str] | None
        Set of valid corrected barcodes.  None = accept all CB values.
    chimeric_distance : int
        Maximum intra-chromosomal same-strand split-alignment distance (bp)
        below which an SA-tag split is treated as a legitimate splice.
    paired_end_chimeric : bool
        When True, supplement SA-tag chimeric detection with paired-end logic:
        a pair is chimeric if the mates map to different chromosomes, same
        chromosome with insert size > ILLUMINA_CHIMERIC_INSERT_SIZE, or one
        mate is unmapped (discordant).  Activate for --platform illumina*.
    tso_sequences : list[str]
        TSO sequences to check for in soft-clipped 5′ bases.
    reference : pysam.FastaFile | None
        Reference FASTA for polyA context check and splice-site dinucleotide
        lookup.  If None, those checks are skipped.
    """

    def __init__(
        self,
        index: AnnotationIndex,
        *,
        barcode_tag: str = BamTag.CELL_BARCODE,
        umi_tag: str = BamTag.UMI,
        whitelist: Optional[set] = None,
        chimeric_distance: int = DEFAULT_CHIMERIC_DISTANCE,
        paired_end_chimeric: bool = False,
        tso_sequences: Optional[list] = None,
        reference: Optional[pysam.FastaFile] = None,
    ):
        self.index = index
        self.barcode_tag = barcode_tag
        self.umi_tag = umi_tag
        self.whitelist = whitelist
        self.chimeric_distance = chimeric_distance
        self.paired_end_chimeric = paired_end_chimeric
        self.tso_sequences = tso_sequences or [TSO_10X, TSO_PACBIO]
        self.reference = reference

        # Pre-build per-contig interval lookup tables for fast query
        # These are built lazily on first access via _get_contig_intervals()
        self._contig_cache: dict[str, dict] = {}

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def classify(self, read: pysam.AlignedSegment) -> Optional[ReadResult]:
        """
        Classify a single aligned read.

        Returns None for reads that should be silently skipped
        (unmapped, secondary, supplementary).
        Supplementary alignments are still inspected for chimeric evidence
        and their SA tag is noted, but they do not produce a ReadResult.
        """
        # --- Skip secondary entirely ---
        if read.flag & SamFlag.SECONDARY:
            return None

        # --- Skip unmapped ---
        if read.flag & SamFlag.UNMAPPED:
            return None

        # --- Supplementary: only used for chimeric signal on primary reads ---
        if read.flag & SamFlag.SUPPLEMENTARY:
            return None

        # --- Extract barcode / UMI ---
        cb, umi = self._get_tags(read)

        # --- Barcode validation ---
        # Three cases:
        # 1. CB present and on whitelist (or no whitelist) → proceed normally
        # 2. CB present but not on whitelist → UNASSIGNED (invalid barcode)
        # 3. CB absent entirely:
        #    - If whitelist provided → UNASSIGNED (expected to have barcodes)
        #    - If no whitelist → use sentinel "NO_BARCODE" and classify
        #      genomically. This enables aggregate noise metrics on BAMs
        #      that have no cell barcode tags (e.g. PacBio FLTNC without
        #      barcode demultiplexing).
        if not cb:
            if self.whitelist is not None:
                # Whitelist provided but CB absent — truly unassigned
                return ReadResult(
                    query_name=read.query_name or "",
                    cell_barcode="",
                    umi=umi,
                    category=ReadCategory.UNASSIGNED,
                    base_counts={ReadCategory.UNASSIGNED: read.query_alignment_length or 0},
                    is_multimapper=False,
                    is_tso_invasion=False,
                    is_polya_priming=False,
                    has_noncanonical_junction=False,
                    contig=read.reference_name or "",
                    pos=read.reference_start or 0,
                    is_reverse=read.is_reverse,
                    read_length=read.query_alignment_length or 0,
                )
            else:
                # No whitelist and no CB — barcode-agnostic mode
                # Use sentinel so all reads aggregate under one bucket
                cb = "NO_BARCODE"
        elif self.whitelist is not None and cb not in self.whitelist:
            # CB present but not valid
            return ReadResult(
                query_name=read.query_name or "",
                cell_barcode=cb,
                umi=umi,
                category=ReadCategory.UNASSIGNED,
                base_counts={ReadCategory.UNASSIGNED: read.query_alignment_length or 0},
                is_multimapper=False,
                is_tso_invasion=False,
                is_polya_priming=False,
                has_noncanonical_junction=False,
                contig=read.reference_name or "",
                pos=read.reference_start or 0,
                is_reverse=read.is_reverse,
                read_length=read.query_alignment_length or 0,
            )

        # --- Check flags ---
        is_multimapper = self._is_multimapper(read)

        # --- Mitochondrial ---
        contig = read.reference_name or ""
        if contig in MITO_CONTIG_NAMES:
            category = ReadCategory.MITOCHONDRIAL
            base_counts = {ReadCategory.MITOCHONDRIAL: read.query_alignment_length or 0}
            return self._make_result(read, cb, umi, category, base_counts, is_multimapper)

        # --- Chimeric ---
        is_chimeric, _ = self._check_chimeric(read)
        if is_chimeric:
            category = ReadCategory.CHIMERIC
            base_counts = {ReadCategory.CHIMERIC: read.query_alignment_length or 0}
            return self._make_result(read, cb, umi, category, base_counts, is_multimapper)

        # --- Artifact flags (TSO, polyA) ---
        is_tso   = self._check_tso_invasion(read)
        is_polya = self._check_polya_priming(read)

        # --- Genomic classification ---
        category, base_counts, has_noncano = self._classify_by_intervals(read)

        result = self._make_result(read, cb, umi, category, base_counts, is_multimapper)
        result.is_tso_invasion           = is_tso
        result.is_polya_priming          = is_polya
        result.has_noncanonical_junction = has_noncano
        return result

    # ------------------------------------------------------------------
    # Tag extraction
    # ------------------------------------------------------------------

    def _get_tags(self, read: pysam.AlignedSegment) -> tuple[str, str]:
        cb  = read.get_tag(self.barcode_tag)  if read.has_tag(self.barcode_tag)  else ""
        umi = read.get_tag(self.umi_tag)       if read.has_tag(self.umi_tag)       else ""
        return str(cb), str(umi)

    # ------------------------------------------------------------------
    # Multi-mapper check
    # ------------------------------------------------------------------

    @staticmethod
    def _is_multimapper(read: pysam.AlignedSegment) -> bool:
        """Primary criterion: NH > 1 (aligner-agnostic)."""
        if read.has_tag(BamTag.NH):
            return int(read.get_tag(BamTag.NH)) > 1
        return False

    # ------------------------------------------------------------------
    # Chimeric read detection
    # ------------------------------------------------------------------

    def _check_chimeric(self, read: pysam.AlignedSegment) -> tuple[bool, Optional[str]]:
        """
        Detect chimeric alignments.

        Returns (is_chimeric, reason_string).

        Two detection paths (evaluated in order):

        1. SA-tag path (long-read and any BAM with split alignments):
           SA tag present AND any of:
             a) SA contig ≠ primary contig  (inter-chromosomal)
             b) SA strand ≠ primary strand  (strand-discordant)
             c) Same contig + same strand, but distance > chimeric_distance

        2. Paired-end path (active only when self.paired_end_chimeric is True):
           Evaluates FLAG bits + RNEXT/PNEXT when no SA tag is present.
           A pair is chimeric if:
             a) One mate is unmapped while the other is mapped (discordant)
             b) Mates map to different chromosomes (inter-chromosomal)
             c) Same chromosome but |TLEN| > ILLUMINA_CHIMERIC_INSERT_SIZE
        """
        if read.has_tag(BamTag.SA):
            sa_str = read.get_tag(BamTag.SA)
            # SA format: rname,pos,strand,CIGAR,mapQ,NM;...
            for sa_entry in sa_str.rstrip(";").split(";"):
                parts = sa_entry.split(",")
                if len(parts) < 3:
                    continue
                sa_contig = parts[0]
                try:
                    sa_pos = int(parts[1])
                except ValueError:
                    continue
                sa_strand = parts[2]  # '+' or '-'

                primary_strand = "-" if read.is_reverse else "+"
                primary_contig = read.reference_name or ""
                primary_pos    = read.reference_start or 0

                # Inter-chromosomal
                if sa_contig != primary_contig:
                    return True, f"inter-chromosomal SA: {sa_contig}"

                # Strand-discordant
                if sa_strand != primary_strand:
                    return True, f"strand-discordant SA on {sa_contig}"

                # Same contig, same strand, but too far
                distance = abs(sa_pos - primary_pos)
                if distance > self.chimeric_distance:
                    return True, (
                        f"same-strand SA on {sa_contig} at distance {distance} bp "
                        f"> threshold {self.chimeric_distance} bp"
                    )
            return False, None

        # No SA tag — fall back to paired-end logic when enabled
        if self.paired_end_chimeric:
            return self._check_chimeric_paired(read)

        return False, None

    def _check_chimeric_paired(
        self, read: pysam.AlignedSegment
    ) -> tuple[bool, Optional[str]]:
        """
        Paired-end chimeric detection using SAM FLAG bits and RNEXT/PNEXT.

        Only meaningful when:
          - The read is marked as paired (FLAG 0x001 set)
          - RNEXT is available (not '*')

        Chimeric criteria:
          a) Mate is unmapped (FLAG 0x008) — discordant pair
          b) Mates map to different chromosomes
          c) Same chromosome, |TLEN| > ILLUMINA_CHIMERIC_INSERT_SIZE (1 Mb)
        """
        # Only evaluate paired reads
        if not (read.flag & SamFlag.PAIRED):
            return False, None

        # (a) Discordant: mate unmapped while this read is mapped
        if read.flag & SamFlag.MATE_UNMAPPED:
            return True, "discordant pair: mate is unmapped"

        # Need RNEXT to be available
        next_ref = read.next_reference_name
        if next_ref is None:
            return False, None

        primary_contig = read.reference_name or ""

        # (b) Inter-chromosomal pair
        if next_ref != primary_contig:
            return True, (
                f"inter-chromosomal pair: read on {primary_contig}, "
                f"mate on {next_ref}"
            )

        # (c) Large insert size on same chromosome
        tlen = abs(read.template_length)
        if tlen > ILLUMINA_CHIMERIC_INSERT_SIZE:
            return True, (
                f"large insert size: {tlen:,} bp "
                f"> threshold {ILLUMINA_CHIMERIC_INSERT_SIZE:,} bp"
            )

        return False, None

    # ------------------------------------------------------------------
    # TSO invasion detection
    # ------------------------------------------------------------------

    def _check_tso_invasion(self, read: pysam.AlignedSegment) -> bool:
        """
        Check for TSO sequence in the 5′ soft-clipped bases.

        A TSO invasion artifact produces a read where the 5′ end contains
        TSO sequence (the TSO primed the nascent cDNA or gDNA directly).

        We check:
          1. Whether the 5′ soft-clip contains a substring matching any known
             TSO sequence (minimum TSO_MIN_MATCH_LENGTH bp).
          2. Whether the 5′ soft-clip is a poly-G run (the TSO poly-G tail
             priming a C-rich template).
        """
        cigar = read.cigartuples
        if not cigar:
            return False

        # 5′ soft clip: cigar op 4 (S) at the start, accounting for strand
        # For a reverse-strand read, pysam reverses the sequence but NOT the
        # cigar — the first cigar tuple still corresponds to the 5′ end of
        # the read as aligned (which is the 3′ end of the molecule).
        # We check both ends for safety.
        clip_seqs = []

        if cigar[0][0] == 4:   # 5′ clip
            clip_len = cigar[0][1]
            seq = read.query_sequence
            if seq and clip_len > 0:
                clip_seqs.append(seq[:clip_len])

        if len(cigar) > 1 and cigar[-1][0] == 4:  # 3′ clip
            clip_len = cigar[-1][1]
            seq = read.query_sequence
            if seq and clip_len > 0:
                clip_seqs.append(seq[-clip_len:])

        for clip in clip_seqs:
            clip_upper = clip.upper()

            # Poly-G check (TSO poly-G tail)
            if "G" * TSO_POLYG_MIN_LENGTH in clip_upper:
                return True

            # TSO sequence check
            for tso in self.tso_sequences:
                tso_check = tso[:TSO_MIN_MATCH_LENGTH].upper()
                if tso_check in clip_upper:
                    return True

        return False

    # ------------------------------------------------------------------
    # Internal poly-A priming detection
    # ------------------------------------------------------------------

    def _check_polya_priming(self, read: pysam.AlignedSegment) -> bool:
        """
        Check for an A-rich stretch immediately downstream of the read's
        3′ mapping end in the reference genome.

        Requires self.reference (pysam.FastaFile) to be set.
        If unavailable, returns False (check is skipped gracefully).

        An A-run of >= POLYA_RUN_MIN_LENGTH within POLYA_CONTEXT_WINDOW bp
        downstream of the aligned 3′ end is flagged.
        """
        if self.reference is None:
            return False

        contig = read.reference_name
        if not contig:
            return False

        end_pos = read.reference_end  # 0-based exclusive
        if end_pos is None:
            return False

        try:
            context = self.reference.fetch(
                contig, end_pos, end_pos + POLYA_CONTEXT_WINDOW
            ).upper()
        except (ValueError, KeyError):
            return False

        a_run_re = re.compile(f"A{{{POLYA_RUN_MIN_LENGTH},}}")
        return bool(a_run_re.search(context))

    # ------------------------------------------------------------------
    # Interval-based genomic classification
    # ------------------------------------------------------------------

    def _classify_by_intervals(
        self, read: pysam.AlignedSegment
    ) -> tuple[ReadCategory, dict, bool]:
        """
        Classify the read by overlapping it against the annotation index.

        Returns
        -------
        category : ReadCategory
            Plurality-vote read-level category (most bases win).
        base_counts : dict[ReadCategory, int]
            Per-category base tally.
        has_noncanonical_junction : bool
        """
        contig   = read.reference_name or ""
        strand   = "-" if read.is_reverse else "+"
        blocks   = read.get_blocks()  # list of (start, end) 0-based half-open

        if not blocks:
            return ReadCategory.AMBIGUOUS, {ReadCategory.AMBIGUOUS: 0}, False

        # Fetch interval lookup tables for this contig
        intervals = self._get_contig_intervals(contig, strand)

        base_counts: dict[ReadCategory, int] = {}
        has_noncanonical = False

        for block_start, block_end in blocks:
            block_len = block_end - block_start

            # --- Base-level tally across all categories ---
            # For each block, compute how many bases fall into each category.
            # Shared/ambiguous bases are counted precisely; the remainder is
            # classified by the standard hierarchy.  This prevents long reads
            # that partially span a shared region from being entirely flagged
            # as ambiguous.

            # Coding-vs-coding shared: genuinely ambiguous
            n_shared_cc = _bases_in(intervals["shared_cod_cod"], block_start, block_end)
            if n_shared_cc > 0:
                _add_bases(base_counts, ReadCategory.AMBIGUOUS_COD_COD, n_shared_cc)

            # Coding-vs-noncoding shared: sub-classifiable ambiguous
            n_shared_cn = _bases_in(intervals["shared_cod_ncod"], block_start, block_end)
            if n_shared_cn > 0:
                _add_bases(base_counts, ReadCategory.AMBIGUOUS_COD_NCOD, n_shared_cn)

            n_shared = n_shared_cc + n_shared_cn

            # Remaining bases after removing shared portion
            remaining = block_len - n_shared
            if remaining <= 0:
                continue

            # For the non-shared portion, apply the standard hierarchy
            n_exon_sense = _bases_in(intervals["exon_sense"], block_start, block_end)
            # Subtract shared bases already counted (avoid double-counting)
            n_exon_sense = max(0, n_exon_sense - n_shared)
            if n_exon_sense > 0:
                _add_bases(base_counts, ReadCategory.EXONIC_SENSE, n_exon_sense)
                remaining -= n_exon_sense

            if remaining <= 0:
                continue

            n_exon_anti = _bases_in(intervals["exon_anti"], block_start, block_end)
            n_exon_anti = max(0, n_exon_anti - n_shared)
            if n_exon_anti > 0:
                _add_bases(base_counts, ReadCategory.EXONIC_ANTISENSE, n_exon_anti)
                remaining -= n_exon_anti

            if remaining <= 0:
                continue

            n_intron_sense = _bases_in(intervals["intron_sense"], block_start, block_end)
            n_intron_anti  = _bases_in(intervals["intron_anti"],  block_start, block_end)
            n_intronic = max(0, (n_intron_sense + n_intron_anti) - n_shared)
            n_intronic = min(n_intronic, remaining)
            if n_intronic > 0:
                _add_bases(base_counts, ReadCategory.INTRONIC_PURE, n_intronic)
                remaining -= n_intronic

            if remaining > 0:
                _add_bases(base_counts, ReadCategory.INTERGENIC_SPARSE, remaining)

        # --- Junction-spanning intronic sub-classification ---
        cigar_jxn_positions = _get_jxn_positions(read)
        if cigar_jxn_positions:
            known_sites = self.index.splice_sites.get(contig, set())
            for jxn_pos in cigar_jxn_positions:
                if (jxn_pos, strand) not in known_sites:
                    # Novel junction: check canonicality if reference available
                    if self.reference is not None:
                        canon = _check_junction_canonicality(
                            self.reference, contig, jxn_pos, strand
                        )
                        if not canon:
                            has_noncanonical = True

            # If intronic bases exist and read has junctions → reclassify
            # intronic bases as INTRONIC_JXNSPAN
            if ReadCategory.INTRONIC_PURE in base_counts:
                n = base_counts.pop(ReadCategory.INTRONIC_PURE)
                _add_bases(base_counts, ReadCategory.INTRONIC_JXNSPAN, n)

        # --- Exon-intron boundary check ---
        # A read that spans an exon-intron boundary without a splice op
        # (no CIGAR N at the junction) is flagged as INTRONIC_BOUNDARY
        if (
            ReadCategory.EXONIC_SENSE in base_counts
            and ReadCategory.INTRONIC_PURE in base_counts
            and not cigar_jxn_positions
        ):
            n_intronic = base_counts.pop(ReadCategory.INTRONIC_PURE)
            _add_bases(base_counts, ReadCategory.INTRONIC_BOUNDARY, n_intronic)

        if not base_counts:
            base_counts[ReadCategory.AMBIGUOUS] = read.query_alignment_length or 0

        # --- Plurality vote for read-level category ---
        category = max(base_counts, key=lambda c: base_counts[c])

        return category, base_counts, has_noncanonical

    # ------------------------------------------------------------------
    # Contig interval cache
    # ------------------------------------------------------------------

    def _get_contig_intervals(self, contig: str, strand: str) -> dict:
        """
        Build and cache a dict of sorted interval lists for *contig*.
        This is the hot path — it runs once per contig per worker process.

        Returns a dict with keys:
          exon_sense, exon_anti, intron_sense, intron_anti, shared
        each a list of (start, end) tuples sorted by start.
        """
        cache_key = f"{contig}:{strand}"
        if cache_key in self._contig_cache:
            return self._contig_cache[cache_key]

        anti = "-" if strand == "+" else "+"

        result = {
            "exon_sense":        _extract_intervals(self.index.exons_plus  if strand == "+" else self.index.exons_minus,  contig),
            "exon_anti":         _extract_intervals(self.index.exons_minus  if strand == "+" else self.index.exons_plus,   contig),
            "intron_sense":      _extract_intervals(self.index.introns_plus if strand == "+" else self.index.introns_minus, contig),
            "intron_anti":       _extract_intervals(self.index.introns_minus if strand == "+" else self.index.introns_plus, contig),
            # strand-aware shared: only same-strand overlaps are truly ambiguous
            "shared_cod_cod":    _extract_intervals(self.index.gene_shared_cod_cod,  contig),
            "shared_cod_ncod":   _extract_intervals(self.index.gene_shared_cod_ncod, contig),
        }
        self._contig_cache[cache_key] = result
        return result

    # ------------------------------------------------------------------
    # ReadResult factory
    # ------------------------------------------------------------------

    @staticmethod
    def _make_result(
        read: pysam.AlignedSegment,
        cb: str,
        umi: str,
        category: ReadCategory,
        base_counts: dict,
        is_multimapper: bool,
    ) -> ReadResult:
        return ReadResult(
            query_name=read.query_name or "",
            cell_barcode=cb,
            umi=umi,
            category=category,
            base_counts=base_counts,
            is_multimapper=is_multimapper,
            is_tso_invasion=False,
            is_polya_priming=False,
            has_noncanonical_junction=False,
            contig=read.reference_name or "",
            pos=read.reference_start or 0,
            is_reverse=read.is_reverse,
            read_length=read.query_alignment_length or 0,
        )


# ---------------------------------------------------------------------------
# Interval utility functions (pure functions, no state)
# ---------------------------------------------------------------------------

def _extract_intervals(pr_obj, contig: str) -> tuple:
    """
    Extract sorted (start, end) tuples for a contig from a PyRanges.

    Returns a tuple (intervals, prefix_max_end) where:
      - intervals is a sorted list of (start, end) tuples
      - prefix_max_end[i] = max(intervals[0..i].end), used for O(1) early exit
        in bisect-based overlap queries.
    """
    if pr_obj is None or pr_obj.df.empty:
        return ([], [])
    df = pr_obj.df
    sub = df[df["Chromosome"] == contig]
    if sub.empty:
        return ([], [])
    intervals = sorted(zip(sub["Start"].tolist(), sub["End"].tolist()))
    prefix_max_end: list[int] = []
    curr_max = 0
    for _, ivl_end in intervals:
        if ivl_end > curr_max:
            curr_max = ivl_end
        prefix_max_end.append(curr_max)
    return (intervals, prefix_max_end)


def _overlaps_any(interval_data: tuple, start: int, end: int) -> bool:
    """
    O(log n) check: does [start, end) overlap any interval in the sorted list?

    interval_data is (intervals, prefix_max_end) as returned by _extract_intervals.
    """
    intervals, prefix_max_end = interval_data
    if not intervals:
        return False
    # Find the first interval with ivl_start >= end (none of those can overlap)
    idx = bisect_left(intervals, (end,))
    if idx == 0:
        return False
    # prefix_max_end[idx-1] = max end among intervals[:idx]; if that's <= start, no overlap
    return prefix_max_end[idx - 1] > start


def _bases_in(interval_data: tuple, start: int, end: int) -> int:
    """
    Count bases in [start, end) that overlap any interval in the sorted list.

    Uses bisect to find the upper boundary in O(log n) and scans backwards with
    prefix_max_end early termination, giving O(log n + k_overlap) per call where
    k_overlap is the number of intervals that actually reach into [start, end).
    """
    intervals, prefix_max_end = interval_data
    if not intervals:
        return 0
    # First interval with ivl_start >= end cannot contribute
    idx = bisect_left(intervals, (end,))
    if idx == 0 or prefix_max_end[idx - 1] <= start:
        return 0
    total = 0
    # Scan backwards; prefix_max_end is non-decreasing left-to-right, so scanning
    # right-to-left it is non-increasing.  Once prefix_max_end[i] <= start, no
    # interval in [0..i] can reach into the query window — early exit.
    for i in range(idx - 1, -1, -1):
        if prefix_max_end[i] <= start:
            break
        ivl_start, ivl_end = intervals[i]
        ov_start = max(start, ivl_start)
        ov_end   = min(end,   ivl_end)
        if ov_end > ov_start:
            total += ov_end - ov_start
    return total


def _add_bases(counts: dict, category: ReadCategory, n: int) -> None:
    counts[category] = counts.get(category, 0) + n


def _get_jxn_positions(read: pysam.AlignedSegment) -> list:
    """
    Return the list of intron start positions (0-based) inferred from
    CIGAR N operations (splice junctions).
    """
    positions = []
    if not read.cigartuples:
        return positions
    pos = read.reference_start or 0
    for op, length in read.cigartuples:
        if op == 3:   # N = skip (intron)
            positions.append(pos)
        if op in (0, 2, 3, 7, 8):   # ops that consume reference
            pos += length
    return positions


def _check_junction_canonicality(
    reference: pysam.FastaFile,
    contig: str,
    jxn_start: int,
    strand: str,
) -> bool:
    """
    Fetch the donor dinucleotide at jxn_start and check for GT-AG / GC-AG /
    AT-AC canonicality.

    Returns True if canonical, False otherwise.
    """
    try:
        donor    = reference.fetch(contig, jxn_start,     jxn_start + 2).upper()
        # We can't easily get the acceptor without knowing intron length,
        # so we check donor only here; acceptor check can be added later.
        if strand == "-":
            # Complement the dinucleotide for reverse-strand genes
            complement = str.maketrans("ACGT", "TGCA")
            donor = donor.translate(complement)[::-1]
        for don, acc in CANONICAL_SPLICE_SITES:
            if donor == don:
                return True
        return False
    except (ValueError, KeyError):
        return True   # Can't check → don't penalise
