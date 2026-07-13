"""
classifier.py
=============
The core per-read classification engine.

For every primary, non-secondary alignment the classifier:

  1. Applies the priority hierarchy (multimapper → mito → chimeric →
     exonic → intronic → intergenic) to assign a single :class:`ReadCategory`.
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
  MITOCHONDRIAL     maps to chrM / MT
  CHIMERIC          SA tag + distance/strand rules
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
from dataclasses import dataclass
from typing import Optional

import pysam

from scnoisemeter.modules.annotation import AnnotationIndex

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


_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def _revcomp(seq: str) -> str:
    """Reverse complement of a DNA string (uppercased; non-ACGTN map to N)."""
    return seq.upper().translate(_COMPLEMENT)[::-1]


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
    is_tso_concatemer : bool
        True if the read sequence contains more than one occurrence of a TSO
        sequence or its reverse complement (template-switch concatemer).
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
    is_tso_concatemer:      bool
    has_junction:           bool
    is_discordant_pair:     bool
    contig:                 str
    pos:                    int
    is_reverse:             bool
    read_length:            int
    aligned_reference_bases: int
    mapq:                   int
    edit_distance:          Optional[int] = None
    softclip_bases:         int = 0
    query_bases:            int = 0


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
        Deprecated compatibility argument. Genomic distance alone is not
        evidence of a chimera because transcript introns can be long.
    paired_end_chimeric : bool
        When True, supplement SA-tag chimeric detection with paired-end logic:
        a pair is chimeric if the mates map to different chromosomes, same
        chromosome with insert size > ILLUMINA_CHIMERIC_INSERT_SIZE. Orphan
        and improper pairs are retained as independent discordance flags.
    tso_sequences : list[str]
        TSO sequences to check for in soft-clipped 5′ bases.
    tso_min_match : int
        Minimum prefix length (bp) of each TSO sequence required to match a
        soft-clip.  Defaults to TSO_MIN_MATCH_LENGTH.
    tso_check_polyg : bool
        If True (default), a poly-G run (≥ TSO_POLYG_MIN_LENGTH G's) in a
        soft-clip also flags TSO invasion.  Set False to require a TSO
        sequence match only.
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
        tso_min_match: int = TSO_MIN_MATCH_LENGTH,
        tso_check_polyg: bool = True,
        reference: Optional[pysam.FastaFile] = None,
    ):
        self.index = index
        self.barcode_tag = barcode_tag
        self.umi_tag = umi_tag
        self.whitelist = whitelist
        self.chimeric_distance = chimeric_distance
        self.paired_end_chimeric = paired_end_chimeric
        # ``None`` preserves the programmatic API's historical defaults;
        # an explicit empty list disables motif matching.  The CLI supplies
        # platform-aware defaults so unrelated chemistries are not screened
        # against every built-in oligo.
        self.tso_sequences = (
            [TSO_10X, TSO_PACBIO]
            if tso_sequences is None
            else list(tso_sequences)
        )
        self.tso_min_match = tso_min_match
        self.tso_check_polyg = tso_check_polyg
        self.reference = reference

        # Match prefixes for TSO-invasion detection: forward AND reverse
        # complement of each TSO, so reads whose TSO end maps antisense (the
        # clip carries revcomp(TSO)) are detected too.  Deduplicated.
        _prefixes = []
        for t in self.tso_sequences:
            up = t.upper()
            _prefixes.append(up[:tso_min_match])
            _prefixes.append(_revcomp(up)[:tso_min_match])
        self._tso_match_prefixes = list(dict.fromkeys(_prefixes))

        # Concatemer detection: non-overlapping count of full TSO sequences and
        # their reverse complements in the read.  Longest-first alternation so a
        # shorter TSO that is a substring of a longer one (e.g. PacBio ⊂ 10x) is
        # not double-counted at the same locus.
        _concat = []
        for t in self.tso_sequences:
            up = t.upper()
            _concat.append(up)
            _concat.append(_revcomp(up))
        _concat = sorted(set(_concat), key=len, reverse=True)
        self._tso_concat_pattern = (
            re.compile("|".join(re.escape(s) for s in _concat)) if _concat else None
        )

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

        aligned_reference_bases = _aligned_reference_bases(read)

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
                return self._decorate_independent_flags(
                    read,
                    self._make_result(
                        read, "", umi, ReadCategory.UNASSIGNED,
                        {ReadCategory.UNASSIGNED: aligned_reference_bases}, False,
                    ),
                )
            else:
                # No whitelist and no CB — barcode-agnostic mode
                # Use sentinel so all reads aggregate under one bucket
                cb = "NO_BARCODE"
        elif self.whitelist is not None and cb not in self.whitelist:
            # CB present but not valid
            return self._decorate_independent_flags(
                read,
                self._make_result(
                    read, cb, umi, ReadCategory.UNASSIGNED,
                    {ReadCategory.UNASSIGNED: aligned_reference_bases}, False,
                ),
            )

        # --- Multimapper (NH > 1) — highest priority per README hierarchy ---
        if self._is_multimapper(read):
            return self._decorate_independent_flags(read, self._make_result(
                read, cb, umi,
                ReadCategory.MULTIMAPPER,
                {ReadCategory.MULTIMAPPER: aligned_reference_bases},
                True,
            ))

        # --- Mitochondrial ---
        contig = read.reference_name or ""
        if contig in MITO_CONTIG_NAMES:
            category = ReadCategory.MITOCHONDRIAL
            base_counts = {ReadCategory.MITOCHONDRIAL: aligned_reference_bases}
            return self._decorate_independent_flags(read, self._make_result(
                read, cb, umi, category, base_counts, False,
            ))

        # --- Chimeric ---
        is_chimeric, _ = self._check_chimeric(read)
        if is_chimeric:
            category = ReadCategory.CHIMERIC
            base_counts = {ReadCategory.CHIMERIC: aligned_reference_bases}
            return self._decorate_independent_flags(read, self._make_result(
                read, cb, umi, category, base_counts, False,
            ))

        # --- Genomic classification ---
        category, base_counts, has_noncano = self._classify_by_intervals(read)

        result = self._make_result(read, cb, umi, category, base_counts, False)
        return self._decorate_independent_flags(
            read, result, has_noncanonical=has_noncano
        )

    def _decorate_independent_flags(
        self,
        read: pysam.AlignedSegment,
        result: ReadResult,
        *,
        has_noncanonical: Optional[bool] = None,
    ) -> ReadResult:
        """Attach evidence flags independently of the primary category."""
        result.is_tso_concatemer = self._check_tso_concatemer(read)
        result.is_tso_invasion = self._check_tso_invasion(read)
        result.is_polya_priming = self._check_polya_priming(read)
        result.is_discordant_pair = self._is_discordant_pair(read)
        result.has_noncanonical_junction = (
            self._has_noncanonical_junction(read)
            if has_noncanonical is None else has_noncanonical
        )
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
        """Detect explicit evidence that a primary alignment is non-unique.

        Only ``NH > 1`` is accepted.  Tags beginning with ``X`` are reserved
        for local use by the SAM specification, so ``XA`` cannot be assigned a
        generic meaning across producers.  MAPQ is likewise a confidence
        score, not a count of reported hits.  Low-MAPQ records are reported
        through an independent QC flag instead of replacing their genomic
        category.
        """
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
             b) SA strand ≠ primary strand  (foldback / strand-discordant)
             c) Query order is incompatible with genomic order.

           A same-contig, same-strand distance is *not* sufficient evidence:
           genomic distance includes introns and can be very large for a
           perfectly legitimate RNA molecule.

        2. Paired-end path (active only when self.paired_end_chimeric is True):
           Evaluates FLAG bits + RNEXT/PNEXT when no SA tag is present.
           A pair is chimeric if:
             a) Mates map to different chromosomes (inter-chromosomal)
             b) Same chromosome but |TLEN| > ILLUMINA_CHIMERIC_INSERT_SIZE

           Orphan and improper pairs are reported through the independent
           discordant-pair flag; they are not by themselves chimeras.
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

                # Same-contig / same-strand: use query-vs-genomic order, not
                # genomic distance. SA POS is 1-based; convert to 0-based.
                sa_start = max(sa_pos - 1, 0)
                sa_cigar = parts[3] if len(parts) > 3 else ""
                primary_cigar = read.cigarstring
                if not isinstance(primary_cigar, str):
                    primary_cigar = ""
                primary_q = _query_interval_from_cigar(primary_cigar)
                sa_q = _query_interval_from_cigar(sa_cigar)
                if primary_q and sa_q and primary_q[0] != sa_q[0]:
                    primary_first = primary_q[0] < sa_q[0]
                    genome_forward = primary_pos <= sa_start
                    expected_forward = primary_strand == "+"
                    order_ok = (
                        genome_forward == expected_forward
                        if primary_first
                        else genome_forward != expected_forward
                    )
                    if not order_ok:
                        return True, "same-contig SA with incompatible query/genomic order"
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

        # An orphan mate is discordant but is not evidence of a chimera.
        if read.flag & SamFlag.MATE_UNMAPPED:
            return False, None

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

    @staticmethod
    def _is_discordant_pair(read: pysam.AlignedSegment) -> bool:
        """Return whether a mapped paired record lacks ordinary pair support."""
        if not (read.flag & SamFlag.PAIRED):
            return False
        if read.flag & SamFlag.MATE_UNMAPPED:
            return True
        return not bool(read.flag & SamFlag.PROPER_PAIR)

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

        # The TSO is at the left soft clip for a forward alignment and the
        # right soft clip for a reverse alignment.  Checking both ends inflated
        # false positives from poly(A), adapters and random terminal sequence.
        clip_seqs = []
        seq = read.query_sequence
        if not read.is_reverse and cigar[0][0] == 4:
            clip_len = cigar[0][1]
            if seq and clip_len > 0:
                clip_seqs.append(seq[:clip_len])
        elif read.is_reverse and cigar[-1][0] == 4:
            clip_len = cigar[-1][1]
            if seq and clip_len > 0:
                clip_seqs.append(seq[-clip_len:])

        for clip in clip_seqs:
            clip_upper = clip.upper()

            # Poly-G check (TSO poly-G tail)
            if self.tso_check_polyg and "G" * TSO_POLYG_MIN_LENGTH in clip_upper:
                return True

            # TSO sequence check (forward + reverse-complement prefixes)
            for tso_check in self._tso_match_prefixes:
                max_mismatches = 1 if len(tso_check) >= 12 else 0
                if _contains_approx(clip_upper, tso_check, max_mismatches):
                    return True

        return False

    def _check_tso_concatemer(self, read: pysam.AlignedSegment) -> bool:
        """
        Check whether the read contains more than one occurrence of a TSO
        sequence or its reverse complement (a template-switch concatemer).

        Matches full TSO sequence(s) with a modest substitution tolerance,
        counting non-overlapping hits across the whole read.  Exact matching
        alone has poor sensitivity on error-prone long reads.
        """
        if self._tso_concat_pattern is None:
            return False
        seq = read.query_sequence
        if not seq:
            return False
        motifs = sorted(
            {s.upper() for s in self.tso_sequences}
            | {_revcomp(s) for s in self.tso_sequences},
            key=len,
            reverse=True,
        )
        hits: list[tuple[int, int]] = []
        upper = seq.upper()
        for motif in motifs:
            max_mm = max(1, int(round(len(motif) * 0.10)))
            for start in _approx_match_starts(upper, motif, max_mm):
                end = start + len(motif)
                if any(not (end <= a or start >= b) for a, b in hits):
                    continue
                hits.append((start, end))
                if len(hits) > 1:
                    return True
        return False

    # ------------------------------------------------------------------
    # Internal poly-A priming detection
    # ------------------------------------------------------------------

    def _check_polya_priming(self, read: pysam.AlignedSegment) -> bool:
        """
        Check for an A-rich stretch immediately downstream of the read's
        transcript 3′ end in the reference genome (strand-aware).

        For a forward read the transcript 3′ end is the rightmost coordinate
        (reference_end), and internal priming shows as an A-run just downstream
        on the + strand.  For a reverse read the transcript 3′ end is the
        leftmost coordinate (reference_start), and the same A-run in the
        transcript reads as a T-run on the + strand just upstream (lower
        coordinate).  Checking only the + strand downstream (the previous
        behaviour) under-detected minus-strand internal priming.

        Requires self.reference (pysam.FastaFile).  Returns False if unavailable.
        An A-run (forward) / T-run (reverse) of >= POLYA_RUN_MIN_LENGTH within
        POLYA_CONTEXT_WINDOW bp of the transcript 3′ end is flagged.
        """
        if self.reference is None:
            return False

        contig = read.reference_name
        if not contig:
            return False

        try:
            if read.is_reverse:
                start = read.reference_start
                if start is None:
                    return False
                lo = max(0, start - POLYA_CONTEXT_WINDOW)
                if lo >= start:
                    return False
                context = self.reference.fetch(contig, lo, start).upper()
                run_re = re.compile(f"T{{{POLYA_RUN_MIN_LENGTH},}}")
            else:
                end_pos = read.reference_end  # 0-based exclusive
                if end_pos is None:
                    return False
                context = self.reference.fetch(
                    contig, end_pos, end_pos + POLYA_CONTEXT_WINDOW
                ).upper()
                run_re = re.compile(f"A{{{POLYA_RUN_MIN_LENGTH},}}")
        except (ValueError, KeyError):
            return False

        return bool(run_re.search(context))

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
        has_noncanonical = self._has_noncanonical_junction(read)

        for block_start, block_end in blocks:
            block_counts = _partition_block(intervals, block_start, block_end)
            for cat, n_bases in block_counts.items():
                _add_bases(base_counts, cat, n_bases)

        # --- Junction-spanning intronic sub-classification ---
        cigar_junctions = _get_junctions(read)
        if cigar_junctions:
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
            and not cigar_junctions
        ):
            n_intronic = base_counts.pop(ReadCategory.INTRONIC_PURE)
            _add_bases(base_counts, ReadCategory.INTRONIC_BOUNDARY, n_intronic)

        if not base_counts:
            base_counts[ReadCategory.AMBIGUOUS] = _aligned_reference_bases(read)

        expected_bases = _aligned_reference_bases(read)
        observed_bases = sum(base_counts.values())
        if observed_bases != expected_bases:
            # This should be unreachable because _partition_block is exhaustive.
            # Keep the invariant exact even for malformed CIGAR/block output.
            delta = expected_bases - observed_bases
            if delta > 0:
                _add_bases(base_counts, ReadCategory.AMBIGUOUS, delta)
            elif delta < 0:
                raise AssertionError(
                    f"base classification exceeded aligned bases for {read.query_name}: "
                    f"{observed_bases} > {expected_bases}"
                )

        # --- Plurality vote for read-level category ---
        category = _plurality_category(base_counts)

        return category, base_counts, has_noncanonical

    def _has_noncanonical_junction(self, read: pysam.AlignedSegment) -> bool:
        """Evaluate exact annotation and both splice motifs for every junction."""
        junctions = _get_junctions(read)
        if not junctions or self.reference is None:
            return False
        contig = read.reference_name or ""
        strand = "-" if read.is_reverse else "+"
        known_junctions = getattr(self.index, "splice_junctions", {}).get(
            contig, set()
        )
        known_sites = self.index.splice_sites.get(contig, set())
        for jxn_start, jxn_end in junctions:
            exact_known = (jxn_start, jxn_end, strand) in known_junctions
            # Compatibility fallback for old cached/mocked indexes.
            if not known_junctions:
                exact_known = (
                    (jxn_start, strand) in known_sites
                    and (jxn_end, strand) in known_sites
                )
            if not exact_known and not _check_junction_canonicality(
                self.reference, contig, jxn_start, jxn_end, strand
            ):
                return True
        return False

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

        result = {
            "exon_sense":        _extract_intervals(self.index.exons_plus  if strand == "+" else self.index.exons_minus,  contig),
            "exon_anti":         _extract_intervals(self.index.exons_minus  if strand == "+" else self.index.exons_plus,   contig),
            "intron_sense":      _extract_intervals(self.index.introns_plus if strand == "+" else self.index.introns_minus, contig),
            "intron_anti":       _extract_intervals(self.index.introns_minus if strand == "+" else self.index.introns_plus, contig),
            # strand-aware shared: only same-strand overlaps are truly ambiguous
            "shared_cod_cod":    _extract_intervals(
                self.index.gene_shared_cod_cod, contig, strand=strand
            ),
            "shared_cod_ncod":   _extract_intervals(
                self.index.gene_shared_cod_ncod, contig, strand=strand
            ),
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
        is_tso_concatemer: bool = False,
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
            is_tso_concatemer=is_tso_concatemer,
            has_junction=bool(_get_junctions(read)),
            is_discordant_pair=False,
            contig=read.reference_name or "",
            pos=read.reference_start or 0,
            is_reverse=read.is_reverse,
            read_length=read.query_alignment_length or 0,
            aligned_reference_bases=_aligned_reference_bases(read),
            mapq=int(read.mapping_quality or 0),
            edit_distance=_edit_distance(read),
            softclip_bases=_softclip_bases(read),
            query_bases=int(read.query_length or read.infer_read_length() or 0),
        )


# ---------------------------------------------------------------------------
# Interval utility functions (pure functions, no state)
# ---------------------------------------------------------------------------

def _extract_intervals(pr_obj, contig: str, strand: Optional[str] = None) -> tuple:
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
    if strand is not None and "Strand" in sub.columns:
        sub = sub[sub["Strand"] == strand]
    if sub.empty:
        return ([], [])
    intervals = _merge_interval_list(
        sorted(zip(sub["Start"].tolist(), sub["End"].tolist()))
    )
    prefix_max_end: list[int] = []
    curr_max = 0
    for _, ivl_end in intervals:
        if ivl_end > curr_max:
            curr_max = ivl_end
        prefix_max_end.append(curr_max)
    return (intervals, prefix_max_end)


def _merge_interval_list(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return the union of a sorted interval list.

    GTF records are merged per gene upstream, so different genes can still
    overlap.  Counting those records independently double-counted bases and
    could make category totals exceed the aligned length.
    """
    if not intervals:
        return []
    merged: list[tuple[int, int]] = []
    cur_start, cur_end = map(int, intervals[0])
    for start, end in intervals[1:]:
        start, end = int(start), int(end)
        if end <= start:
            continue
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end
    merged.append((cur_start, cur_end))
    return merged


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


_BASE_CATEGORY_PRIORITY = [
    ("shared_cod_cod", ReadCategory.AMBIGUOUS_COD_COD),
    ("shared_cod_ncod", ReadCategory.AMBIGUOUS_COD_NCOD),
    ("exon_sense", ReadCategory.EXONIC_SENSE),
    ("exon_anti", ReadCategory.EXONIC_ANTISENSE),
    ("intron_sense", ReadCategory.INTRONIC_PURE),
    ("intron_anti", ReadCategory.INTRONIC_PURE),
]

_READ_TIEBREAK_PRIORITY = [
    ReadCategory.AMBIGUOUS_COD_COD,
    ReadCategory.AMBIGUOUS_COD_NCOD,
    ReadCategory.EXONIC_SENSE,
    ReadCategory.EXONIC_ANTISENSE,
    ReadCategory.INTRONIC_JXNSPAN,
    ReadCategory.INTRONIC_BOUNDARY,
    ReadCategory.INTRONIC_PURE,
    ReadCategory.INTERGENIC_SPARSE,
    ReadCategory.AMBIGUOUS,
]


def _overlapping_intervals(interval_data: tuple, start: int, end: int) -> list:
    """Return clipped union intervals overlapping ``[start, end)``."""
    intervals, prefix_max_end = interval_data
    if not intervals:
        return []
    idx = bisect_left(intervals, (end,))
    out = []
    for i in range(idx - 1, -1, -1):
        if prefix_max_end[i] <= start:
            break
        ivl_start, ivl_end = intervals[i]
        lo, hi = max(start, ivl_start), min(end, ivl_end)
        if hi > lo:
            out.append((lo, hi))
    return out


def _partition_block(intervals: dict, start: int, end: int) -> dict:
    """Partition one aligned reference block into disjoint categories.

    The sweep operates on atomic segments induced by all annotation interval
    boundaries.  At each segment exactly one category wins according to the
    explicit hierarchy above, guaranteeing that output bases sum to the block
    length even at complex multi-gene overlaps.
    """
    if end <= start:
        return {}
    events: dict[int, list[tuple[int, str]]] = {}
    for key, _category in _BASE_CATEGORY_PRIORITY:
        for lo, hi in _overlapping_intervals(intervals[key], start, end):
            events.setdefault(lo, []).append((1, key))
            events.setdefault(hi, []).append((-1, key))

    points = sorted({start, end, *events.keys()})
    active: dict[str, int] = {}
    counts: dict[ReadCategory, int] = {}
    for i, pos in enumerate(points[:-1]):
        # Half-open intervals ending at pos are removed before starts at pos.
        for delta, key in sorted(events.get(pos, []), key=lambda x: x[0]):
            active[key] = active.get(key, 0) + delta
            if active[key] <= 0:
                active.pop(key, None)
        next_pos = points[i + 1]
        if next_pos <= pos:
            continue
        category = ReadCategory.INTERGENIC_SPARSE
        for key, candidate in _BASE_CATEGORY_PRIORITY:
            if active.get(key, 0) > 0:
                category = candidate
                break
        _add_bases(counts, category, next_pos - pos)
    return counts


def _plurality_category(base_counts: dict[ReadCategory, int]) -> ReadCategory:
    """Return the largest base category with deterministic tie-breaking."""
    rank = {cat: i for i, cat in enumerate(_READ_TIEBREAK_PRIORITY)}
    return max(
        base_counts,
        key=lambda cat: (base_counts[cat], -rank.get(cat, len(rank))),
    )


def _aligned_reference_bases(read: pysam.AlignedSegment) -> int:
    """Number of reference bases covered by aligned blocks, consistently."""
    try:
        return sum(max(0, int(end) - int(start)) for start, end in read.get_blocks())
    except Exception:
        return int(read.query_alignment_length or 0)


def _edit_distance(read: pysam.AlignedSegment) -> Optional[int]:
    """Return NM when present; absence is distinct from zero edits."""
    try:
        return int(read.get_tag("NM")) if read.has_tag("NM") else None
    except (TypeError, ValueError):
        return None


def _softclip_bases(read: pysam.AlignedSegment) -> int:
    try:
        return sum(length for op, length in (read.cigartuples or []) if op == 4)
    except Exception:
        return 0


def _query_interval_from_cigar(cigar: str) -> Optional[tuple[int, int]]:
    """Approximate aligned query interval from a SAM CIGAR string."""
    if not cigar:
        return None
    tokens = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not tokens:
        return None
    leading = int(tokens[0][0]) if tokens[0][1] in {"S", "H"} else 0
    aligned = sum(int(n) for n, op in tokens if op in {"M", "I", "=", "X"})
    return (leading, leading + aligned)


def _approx_match_starts(sequence: str, motif: str, max_mismatches: int) -> list[int]:
    """Return substitution-tolerant motif starts (no indel model)."""
    if not motif or len(sequence) < len(motif):
        return []
    starts = []
    width = len(motif)
    for start in range(len(sequence) - width + 1):
        mismatches = 0
        for a, b in zip(sequence[start:start + width], motif):
            if a != b:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        if mismatches <= max_mismatches:
            starts.append(start)
    return starts


def _contains_approx(sequence: str, motif: str, max_mismatches: int) -> bool:
    return bool(_approx_match_starts(sequence, motif, max_mismatches))


def _get_junctions(read: pysam.AlignedSegment) -> list[tuple[int, int]]:
    """Return exact ``(intron_start, intron_end)`` CIGAR-N intervals."""
    junctions = []
    if not read.cigartuples:
        return junctions
    pos = read.reference_start or 0
    for op, length in read.cigartuples:
        if op == 3:
            junctions.append((pos, pos + length))
        if op in (0, 2, 3, 7, 8):
            pos += length
    return junctions


def _get_jxn_positions(read: pysam.AlignedSegment) -> list:
    """
    Return the list of intron start positions (0-based) inferred from
    CIGAR N operations (splice junctions).
    """
    return [start for start, _end in _get_junctions(read)]


def _check_junction_canonicality(
    reference: pysam.FastaFile,
    contig: str,
    jxn_start: int,
    jxn_end: int,
    strand: str,
) -> bool:
    """
    Fetch both donor and acceptor dinucleotides and check GT-AG / GC-AG /
    AT-AC canonicality in transcript orientation.

    Returns True if canonical, False otherwise.
    """
    try:
        if jxn_end - jxn_start < 4:
            return False
        if strand == "+":
            donor = reference.fetch(contig, jxn_start, jxn_start + 2).upper()
            acceptor = reference.fetch(contig, jxn_end - 2, jxn_end).upper()
        else:
            donor = _revcomp(reference.fetch(contig, jxn_end - 2, jxn_end))
            acceptor = _revcomp(reference.fetch(contig, jxn_start, jxn_start + 2))
        return (donor, acceptor) in CANONICAL_SPLICE_SITES
    except (ValueError, KeyError):
        return True   # Can't check → don't penalise
