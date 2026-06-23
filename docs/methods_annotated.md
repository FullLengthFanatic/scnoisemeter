# scNoiseMeter: annotated methods and code map

**A detailed, code-linked companion to [whitepaper.md](whitepaper.md).** Version 0.6.0.

This document explains what scNoiseMeter does, why each rule is written the way it is, and exactly where in the source each rule and threshold lives, so you know what to change to alter its behavior. Every code reference is `file:line` against the source tree under `scnoisemeter/`. Code blocks are quoted verbatim from that source.

It also compares the tool to existing software by name, states its strengths and limitations honestly, and ends with a list of discrepancies found between the code and its own documentation during this review.

---

## Contents

1. [Architecture and data flow](#1-architecture-and-data-flow)
2. [The single source of truth: constants.py](#2-the-single-source-of-truth-constantspy)
3. [The classification decision tree](#3-the-classification-decision-tree)
4. [Sequence-level artifact detection](#4-sequence-level-artifact-detection)
5. [The intergenic second pass](#5-the-intergenic-second-pass)
6. [The annotation model](#6-the-annotation-model)
7. [Metrics and noise definitions](#7-metrics-and-noise-definitions)
8. [Pipeline mechanics](#8-pipeline-mechanics)
9. [Platform and chemistry auto-detection](#9-platform-and-chemistry-auto-detection)
10. [CLI surface and subcommands](#10-cli-surface-and-subcommands)
11. [Output files](#11-output-files)
12. [Comparison to existing tools](#12-comparison-to-existing-tools)
13. [Strengths and limitations](#13-strengths-and-limitations)
14. [Discrepancies found during review](#14-discrepancies-found-during-review)
15. [Tuning cheat-sheet](#15-tuning-cheat-sheet)

---

## 1. Architecture and data flow

The run path is linear and easy to follow:

```
cli.py (parse flags, resolve TSO/seed/platform)
  → bam_inspector.py        detect platform, stage, barcode presence from the BAM header
  → annotation.py           build (or load cached) AnnotationIndex from the GTF
  → pipeline.py             one worker process per contig; each worker runs:
       → classifier.py      ReadClassifier.classify() on every primary alignment
  → intergenic_profiler.py  second pass: cluster + Poisson-score intergenic reads
  → metrics.py              aggregate to sample / cell / cluster, compute noise fractions
  → report.py / report_figures.py   HTML report + MultiQC JSON
```

Module sizes and roles:

| File | Role |
|---|---|
| `scnoisemeter/constants.py` | Every category name, flag bit, tag, and tunable threshold |
| `scnoisemeter/modules/classifier.py` | Per-read classification engine (`ReadClassifier`) |
| `scnoisemeter/modules/intergenic_profiler.py` | Locus clustering + Poisson scoring + promotion |
| `scnoisemeter/modules/annotation.py` | GTF parsing, interval model, caching |
| `scnoisemeter/modules/metrics.py` | Noise fractions, per-cell / per-cluster aggregation |
| `scnoisemeter/modules/pipeline.py` | Per-contig parallelism, reservoir sampling, merging |
| `scnoisemeter/modules/report.py`, `report_figures.py` | HTML report + MultiQC JSON |
| `scnoisemeter/utils/bam_inspector.py` | Platform / stage / barcode auto-detection |
| `scnoisemeter/utils/annotation_fetcher.py` | Auto-download of GTF, polyA, TSS, whitelists |
| `scnoisemeter/utils/sample_sheet.py` | Plate sample-sheet parsing |
| `scnoisemeter/utils/discover_inspector.py` | Directory scan for the `discover` subcommand |

The classifier is stateless: one `ReadClassifier` per worker, configuration only (thresholds, whitelist, reference), no per-read state. Aggregation is the caller's job (`classifier.py:39-41`).

---

## 2. The single source of truth: constants.py

Every tunable value lives in `constants.py`. The module docstring states the intent (`constants.py:7-11`): "Every category string and every default threshold lives here. No other module should hard-code these values." That is mostly true; the two exceptions found during review are noted in [section 14](#14-discrepancies-found-during-review).

| Constant | `constants.py` line | Value | What it gates | CLI flag |
|---|---|---|---|---|
| `DEFAULT_CHIMERIC_DISTANCE` | 271 | 10 000 | Max same-strand intra-chromosomal `SA` split distance before a read is chimeric | `--chimeric-distance` |
| `ILLUMINA_CHIMERIC_INSERT_SIZE` | 276 | 1 000 000 | Max paired-end insert before a pair is chimeric | none |
| `TSO_10X` | 247 | `AAGCAGTGGTATCAACGCAGAGTACATGGG` | Default 10x TSO | replaced by `--tso` |
| `TSO_PACBIO` | 250 | `AAGCAGTGGTATCAACGCAGAGT` | Default PacBio TSO | replaced by `--tso` |
| `TSO_MIN_MATCH_LENGTH` | 257 | 12 | TSO prefix length required to match a soft-clip | `--tso-min-match` |
| `TSO_POLYG_MIN_LENGTH` | 260 | 6 | Poly-G run length that flags TSO invasion | disable with `--no-polyg-tso` |
| `ADAPTIVE_MIN_BARCODES_FRACTION` | 285 | 0.0001 | Min distinct barcodes per intergenic locus, as a fraction of total | none |
| `ADAPTIVE_MIN_BARCODES_ABSOLUTE` | 286 | 3 | Hard floor on barcodes per locus | none |
| `ADAPTIVE_MIN_READS` | 289 | 5 | Min reads per locus for the Poisson test | none |
| `ADAPTIVE_PVALUE_THRESHOLD` | 292 | 0.01 | Bonferroni-corrected significance threshold | none |
| `INTERGENIC_LOCUS_WINDOW` | 295 | 500 | Single-linkage gap for clustering intergenic reads | none |
| `POLYA_RUN_MIN_LENGTH` | 304 | 6 | A-run (or T-run) length for internal-priming | none |
| `POLYA_CONTEXT_WINDOW` | 305 | 20 | bp window past the 3' end searched for the A-run | none |
| `POLYA_SITE_PROXIMITY` | 309 | 50 | Distance to an annotated polyA site that counts as "proximal" | none |
| `TSS_SITE_PROXIMITY` | 326 | 100 | Distance to a CAGE/TSS peak for the 5'-anchored metric | none |
| `MITO_CONTIG_NAMES` | 347 | `{chrM, MT, chrMT, mitochondrion}` | Mitochondrial contig names | none |
| `BARCODE_AUTODETECT_SAMPLE_SIZE` | 355 | 10 000 | Reads sampled to test for `CB` presence | none |
| `BARCODE_AUTODETECT_MIN_FRACTION` | 356 | 0.50 | `CB`-bearing fraction to enter barcode-aware mode | none |
| `DEFAULT_THREADS` | 363 | 4 | Worker processes | `--threads` |
| `LENGTH_BIN_BREAKS` | 373 | `[150,500,1000,2000,5000]` | Read-length stratification bins | none |
| `LENGTH_SHORT_READ_THRESHOLD` | 382 | 300 | Median length below which the `<150` bin is added | none |
| `CANONICAL_SPLICE_SITES` | 317-321 | GT-AG, GC-AG, AT-AC | Canonical donor dinucleotides | none |
| `NUMT_OVERLAP_FRACTION` | 331 | 0.80 | Read fraction overlapping a NUMT to flag it | `--numt-bed` enables |

Two constants worth knowing about: `MIN_READS_PER_CELL = 10` lives in `metrics.py:72` (not constants.py), and `FULL_LENGTH_THRESHOLD` (the length fallback for the full-length metric) lives in `metrics.py:79-83`. See [section 14](#14-discrepancies-found-during-review) for the consequence.

To change any "none"-flag threshold you edit `constants.py` and rerun. There is no CLI override for the intergenic and polyA thresholds by design; they are statistical and biological constants.

---

## 3. The classification decision tree

`ReadClassifier.classify()` (`classifier.py:240-358`) is the canonical priority hierarchy. The first matching rule assigns the single read-level category. Here is the dispatch in order.

### 3.1 Skip filters

```python
# classifier.py:250-259
if read.flag & SamFlag.SECONDARY:
    return None
if read.flag & SamFlag.UNMAPPED:
    return None
if read.flag & SamFlag.SUPPLEMENTARY:
    return None
```

Secondary, unmapped, and supplementary records produce no result. The supplementary record's `SA` tag is still read off the primary alignment in the chimera check, so split alignments are not lost. *To change:* this is structural; only the chimera detector consumes supplementary evidence.

### 3.2 Barcode gate

```python
# classifier.py:274-314 (abridged)
if not cb:
    if self.whitelist is not None:
        return ReadResult(..., category=ReadCategory.UNASSIGNED, ...)   # CB expected but absent
    else:
        cb = "NO_BARCODE"                                               # barcode-agnostic mode
elif self.whitelist is not None and cb not in self.whitelist:
    return ReadResult(..., category=ReadCategory.UNASSIGNED, ...)        # CB present, off-list
```

Three cases (`classifier.py:264-314`): a whitelist plus a missing or off-list `CB` gives `unassigned`; no whitelist plus no `CB` pools everything under a `NO_BARCODE` sentinel so sample-level metrics still work on barcode-free BAMs (PacBio FLNC, for example). *To change the tag names:* `--barcode-tag` / `--umi-tag` (defaults `CB` / `UB`, `constants.py:165,171`). *To supply a whitelist:* `--barcode-whitelist` or `--whitelist-db`.

Note the distinction from `--cell-barcodes`, which is a *called-cell* filter applied earlier in the worker (see [section 8](#8-pipeline-mechanics)): off-list reads there are dropped entirely and never reach this function, whereas `unassigned` reads are counted.

### 3.3 Multimapper, mitochondrial, chimeric

```python
# classifier.py:318-344 (abridged)
is_concat = self._check_tso_concatemer(read)        # independent flag, every read

if self._is_multimapper(read):                       # NH > 1
    return self._make_result(read, cb, umi, ReadCategory.MULTIMAPPER, ...)

contig = read.reference_name or ""
if contig in MITO_CONTIG_NAMES:                       # chrM / MT / chrMT / mitochondrion
    return self._make_result(read, cb, umi, ReadCategory.MITOCHONDRIAL, ...)

is_chimeric, _ = self._check_chimeric(read)
if is_chimeric:
    return self._make_result(read, cb, umi, ReadCategory.CHIMERIC, ...)
```

Multimapper is `NH > 1` (`classifier.py:373-378`), aligner-agnostic. Mitochondrial is a contig-name lookup (`classifier.py:330-336`). Note the actual order is multimapper, then mitochondrial, then chimeric; this differs from the README hierarchy table and is flagged in [section 14](#14-discrepancies-found-during-review).

Chimeric detection (`_check_chimeric`, `classifier.py:384-444`) has two paths. The `SA`-tag path:

```python
# classifier.py:424-437
if sa_contig != primary_contig:
    return True, f"inter-chromosomal SA: {sa_contig}"
if sa_strand != primary_strand:
    return True, f"strand-discordant SA on {sa_contig}"
distance = abs(sa_pos - primary_pos)
if distance > self.chimeric_distance:
    return True, (...)
```

Rationale, from the constant's own comment (`constants.py:267-271`): "polyA+ cDNA molecules are typically <3 kb; RT processivity limits capture to ~10-15 kb even for long transcripts. 10 kb is generous." *To change:* `--chimeric-distance`. Raise it for datasets with transcripts longer than 10 kb to avoid flagging legitimate long splices.

The paired-end fallback (`_check_chimeric_paired`, `classifier.py:446-491`) fires only when `paired_end_chimeric` is set (Illumina-family and Smart-seq) and there is no `SA` tag. It flags an unmapped mate, an inter-chromosomal pair, or `|TLEN| > 1 000 000` (`ILLUMINA_CHIMERIC_INSERT_SIZE`, `constants.py:276`). On standard 10x and BD BAMs only R2 is aligned, so this fallback rarely fires there; it matters for bulk-style paired Illumina and for Smart-seq / FLASH-seq where both mates are cDNA.

### 3.4 Genomic classification

For reads that survive to here, `_classify_by_intervals()` (`classifier.py:624-746`) tallies aligned bases per category against the annotation index, then takes a plurality vote for the read-level label. The base tally avoids double-counting by removing shared (ambiguous) bases first, then applying the hierarchy to the remainder:

```python
# classifier.py:662-707 (abridged)
n_shared_cc = _bases_in(intervals["shared_cod_cod"], block_start, block_end)
if n_shared_cc > 0:
    _add_bases(base_counts, ReadCategory.AMBIGUOUS_COD_COD, n_shared_cc)
n_shared_cn = _bases_in(intervals["shared_cod_ncod"], block_start, block_end)
if n_shared_cn > 0:
    _add_bases(base_counts, ReadCategory.AMBIGUOUS_COD_NCOD, n_shared_cn)
n_shared = n_shared_cc + n_shared_cn
...
n_exon_sense = max(0, _bases_in(intervals["exon_sense"], ...) - n_shared)   # then exon_anti, then intronic
...
if remaining > 0:
    _add_bases(base_counts, ReadCategory.INTERGENIC_SPARSE, remaining)
```

Interval overlap uses a bisect-based scan with a prefix-max-end early exit (`_bases_in`, `classifier.py:857-884`), so the hot path is `O(log n + k)` per block. Interval sets are cached per contig and strand (`_get_contig_intervals`, `classifier.py:752-777`), with the key insight that shared regions are only ambiguous on the same strand (`classifier.py:772`).

The plurality vote (`classifier.py:744`):

```python
category = max(base_counts, key=lambda c: base_counts[c])
```

### 3.5 Intronic subtypes

After the base tally, two reclassifications run. A read with a `CIGAR N` has its intronic bases promoted from pure to junction-spanning:

```python
# classifier.py:709-727 (abridged)
cigar_jxn_positions = _get_jxn_positions(read)
if cigar_jxn_positions:
    known_sites = self.index.splice_sites.get(contig, set())
    for jxn_pos in cigar_jxn_positions:
        if (jxn_pos, strand) not in known_sites:
            if self.reference is not None:
                canon = _check_junction_canonicality(self.reference, contig, jxn_pos, strand)
                if not canon:
                    has_noncanonical = True
    if ReadCategory.INTRONIC_PURE in base_counts:
        n = base_counts.pop(ReadCategory.INTRONIC_PURE)
        _add_bases(base_counts, ReadCategory.INTRONIC_JXNSPAN, n)
```

A read that has both exonic and intronic bases but no splice gap is the incomplete-RT signature, `intronic_boundary`:

```python
# classifier.py:732-738
if (
    ReadCategory.EXONIC_SENSE in base_counts
    and ReadCategory.INTRONIC_PURE in base_counts
    and not cigar_jxn_positions
):
    n_intronic = base_counts.pop(ReadCategory.INTRONIC_PURE)
    _add_bases(base_counts, ReadCategory.INTRONIC_BOUNDARY, n_intronic)
```

Junction positions come from `CIGAR` op 3 (`_get_jxn_positions`, `classifier.py:891-905`). Canonicality (`_check_junction_canonicality`, `classifier.py:908-933`) checks only the **donor** dinucleotide against `CANONICAL_SPLICE_SITES`; the code comment notes the acceptor is skipped because it would need the intron length (`classifier.py:922-923`), and a failed reference fetch returns `True` so a read is never penalized for a lookup error (`classifier.py:932-933`). *To enable this check at all:* supply `--reference`.

### 3.6 The 16 categories

Defined in the `ReadCategory` enum (`constants.py:21-70`), ordered for output in `CATEGORY_ORDER` (`constants.py:74-93`). `exonic_sense` is signal; the noise membership of the rest is decided in `metrics.py` (see [section 7](#7-metrics-and-noise-definitions)).

---

## 4. Sequence-level artifact detection

These flags are independent of the read category. They are computed in the classifier and carried on `ReadResult` (`classifier.py:124-137`), then counted per cell and per sample.

### 4.1 TSO invasion

The match prefixes are built once in the constructor, forward and reverse-complement, truncated to `tso_min_match`:

```python
# classifier.py:208-216
_prefixes = []
for t in self.tso_sequences:
    up = t.upper()
    _prefixes.append(up[:tso_min_match])
    _prefixes.append(_revcomp(up)[:tso_min_match])
self._tso_match_prefixes = list(dict.fromkeys(_prefixes))
```

The reverse-complement (added in v0.6) catches reads whose TSO end mapped antisense (`classifier.py:208-210`). The check itself scans both soft-clips, with a poly-G shortcut for the 10x TSO's G-rich tail:

```python
# classifier.py:533-545
for clip in clip_seqs:
    clip_upper = clip.upper()
    if self.tso_check_polyg and "G" * TSO_POLYG_MIN_LENGTH in clip_upper:
        return True
    for tso_check in self._tso_match_prefixes:
        if tso_check in clip_upper:
            return True
return False
```

*To change:* `--tso SEQ` (repeatable, **replaces** the built-in defaults), `--tso-min-match N` (default 12; note a TSO shorter than N matches at its full length, not "never"), `--no-polyg-tso` to drop the poly-G heuristic. The poly-G check is independent of the sequence list. Validation of user TSOs (ACGTN only, uppercased, short-sequence warning) is in `_resolve_tso` (`cli.py:332-364`).

### 4.2 TSO concatemer

A separate regex of full TSO sequences plus reverse-complements, longest-first so a shorter TSO that is a substring of a longer one is not double-counted (`classifier.py:218-230`). The check returns true on the second non-overlapping hit:

```python
# classifier.py:557-567
if self._tso_concat_pattern is None:
    return False
seq = read.query_sequence
if not seq:
    return False
n = 0
for _ in self._tso_concat_pattern.finditer(seq.upper()):
    n += 1
    if n > 1:
        return True
return False
```

The docstring ties the metric to its source (`classifier.py:552-555`): "Mirrors the metric in Chou et al. (bioRxiv 2025.10.06.680646): reads with > 1 TSO-or-revcomp occurrence divided by total reads." That paper ("Single-cell RNA-seq using UltraMarathonRT expands the known transcriptome") reports single-cell TSO-concatemer rates of about 2 to 6% depending on protocol, which is the benchmark this metric is meant to reproduce. *To make it meaningful:* run on untrimmed BAMs (if an upstream step removed the TSO, there is nothing to count) and pass the correct `--tso` for non-10x chemistries.

### 4.3 Internal polyA priming

Strand-aware (v0.6). A forward read looks downstream of `reference_end` for an A-run; a reverse read looks upstream of `reference_start` for a T-run:

```python
# classifier.py:597-614
if read.is_reverse:
    start = read.reference_start
    ...
    lo = max(0, start - POLYA_CONTEXT_WINDOW)
    if lo >= start:
        return False
    context = self.reference.fetch(contig, lo, start).upper()
    run_re = re.compile(f"T{{{POLYA_RUN_MIN_LENGTH},}}")
else:
    end_pos = read.reference_end  # 0-based exclusive
    ...
    context = self.reference.fetch(contig, end_pos, end_pos + POLYA_CONTEXT_WINDOW).upper()
    run_re = re.compile(f"A{{{POLYA_RUN_MIN_LENGTH},}}")
```

The docstring explains the strand logic (`classifier.py:574-588`): the previous +-strand-only check under-detected minus-strand priming. *To change:* `POLYA_RUN_MIN_LENGTH` (6) and `POLYA_CONTEXT_WINDOW` (20) in `constants.py:304-305`. *To enable:* `--reference`. Without a reference FASTA the flag is always false.

### 4.4 Reverse complement helper

```python
# classifier.py:75-80
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
def _revcomp(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT)[::-1]
```

Non-ACGTN bases map to N, a safe fallback.

---

## 5. The intergenic second pass

`intergenic_profiler.py` re-examines every read the classifier left as `intergenic_sparse`. The module docstring (`intergenic_profiler.py:1-53`) lays out the model.

### 5.1 Clustering

Single-linkage by contig, strand, and proximity, in `O(n log n)`:

```python
# intergenic_profiler.py:434-447
for idx in order:
    r = records[idx]
    if (
        r.contig != current_contig
        or r.strand != current_strand
        or r.start > current_end + INTERGENIC_LOCUS_WINDOW
    ):
        current_locus += 1
        current_end = r.end
        ...
    else:
        current_end = max(current_end, r.end)
    locus_ids[idx] = current_locus
```

*To change the merge gap:* `INTERGENIC_LOCUS_WINDOW` (500, `constants.py:295`).

### 5.2 Poisson background

The rate is genome-wide intergenic read density (`intergenic_profiler.py:210-215`). Per locus:

```python
# intergenic_profiler.py:297-305
locus_width  = max(end - start, 1)
expected     = background_rate * locus_width
raw_pvalue   = 1.0 - poisson.cdf(n_reads - 1, expected) if expected > 0 else 1.0
adj_pvalue   = min(raw_pvalue * n_tests, 1.0)   # Bonferroni
significant  = (
    adj_pvalue < alpha
    and n_reads >= ADAPTIVE_MIN_READS
    and n_barcodes >= min_barcodes
)
```

It is a one-tailed upper-tail Poisson test, Bonferroni-corrected by the number of loci (`n_tests`, `intergenic_profiler.py:230`). The barcode floor scales with sample size (`intergenic_profiler.py:220-223`):

```python
min_barcodes = max(
    ADAPTIVE_MIN_BARCODES_ABSOLUTE,                       # 3
    int(total_barcodes * ADAPTIVE_MIN_BARCODES_FRACTION), # 0.01% of total
)
```

The total intergenic base count (the denominator) comes from `compute_intergenic_bases` (`intergenic_profiler.py:551-559`). *To change:* `ADAPTIVE_MIN_READS` (5), `ADAPTIVE_PVALUE_THRESHOLD` (0.01), `ADAPTIVE_MIN_BARCODES_*` (`constants.py:285-292`).

### 5.3 Promotion rules

```python
# intergenic_profiler.py:329-344
if overlaps_repeat:
    category = ReadCategory.INTERGENIC_REPEAT
elif not significant:
    category = ReadCategory.INTERGENIC_SPARSE
elif _is_hotspot(is_monoexonic, polya_run_downstream, near_polya):
    category = ReadCategory.INTERGENIC_HOTSPOT
elif _is_novel_gene(has_junction, near_polya, strand_consistent, n_barcodes, min_barcodes):
    category = ReadCategory.INTERGENIC_NOVEL
else:
    category = ReadCategory.INTERGENIC_HOTSPOT      # significant but ambiguous → flag for review
```

The two rule functions are small and testable:

```python
# intergenic_profiler.py:367-399
def _is_hotspot(is_monoexonic, polya_run_downstream, near_polya_site):
    return is_monoexonic and polya_run_downstream

def _is_novel_gene(has_splice_evidence, near_polya_site, strand_consistent, n_barcodes, min_barcodes):
    return (
        strand_consistent
        and n_barcodes >= min_barcodes
        and (has_splice_evidence or near_polya_site)
    )
```

Strand consistency is `>= 0.80` on the dominant strand (`STRAND_CONSISTENCY_MIN`, `intergenic_profiler.py:80`, used at `:287`). The hotspot's A-run context is checked strand-aware against the modal 3' end (`_check_polya_context`, `intergenic_profiler.py:467-489`), with the same downstream-A / upstream-T logic as the per-read flag. PolyA-site proximity uses bisect (`_near_polya_site`, `intergenic_profiler.py:496-524`), and this same helper is reused for the 5' TSS metric with a 100 bp window.

Reads at a promoted locus are written back into `record_categories` (`intergenic_profiler.py:245-247`) and applied to the per-cell counts before metrics, so promoting a locus to `intergenic_novel` (ambiguous, not noise) lowers the reported noise fraction.

> **Note, flagged in [section 14](#14-discrepancies-found-during-review):** the README and this module's own docstring (`intergenic_profiler.py:30-36`) say a hotspot must also be `> 50 bp from any annotated polyA site`, but `_is_hotspot` does not enforce `near_polya`. The function docstring (`intergenic_profiler.py:376-377`) admits the distance "is not strictly required." So the behavior is the function, not the prose.

---

## 6. The annotation model

`annotation.py` turns a GTF into an `AnnotationIndex` (`annotation.py:78-124`) of PyRanges interval sets: exons and introns per strand, gene bodies, unique vs shared regions, the intergenic complement, splice sites, and optional repeats.

Design decisions, from the build docstring (`annotation.py:1-51`):

- **Exon union is per gene, not per transcript** (`_merge_exons_per_gene`, `annotation.py:335-361`), so a base that is intronic in one isoform but exonic in another is treated as exonic. This minimizes false intronic calls.
- **Introns are gene body minus exon union** (`_compute_introns`, `annotation.py:368-381`).
- **Shared vs unique is computed on exons, not gene bodies** (`_unique_and_shared`, `annotation.py:426-521`). The comment explains why (`annotation.py:262-269`): computing on gene bodies "inflates ambiguous 10-fold because large genes (e.g. AGRN, 100 kb) have many other genes nested in their introns."

Overlap type is decided by biotype:

```python
# annotation.py:408-423
def _classify_overlap_type(biotype_a, biotype_b):
    a_coding = biotype_a in _CODING_BIOTYPES
    b_coding = biotype_b in _CODING_BIOTYPES
    if a_coding and b_coding:
        return "cod_cod"
    if a_coding or b_coding:
        return "cod_ncod"
    return "ncod_ncod"
```

The coding and non-coding biotype sets are explicit (`annotation.py:391-405`); `_CODING_BIOTYPES` is `protein_coding` plus the IG/TR gene segments, everything else (lncRNA, the pseudogene family, the small RNAs) is non-coding. *To exclude biotypes from the index entirely:* `--exclude-biotypes` (repeatable).

The intergenic complement (`_manual_complement`, `annotation.py:588-637`) deliberately omits the region after the last gene on each chromosome, because chromosome lengths are not known from a GTF. The comment (`annotation.py:594-604`) notes the previous sentinel `End=2_000_000_000` inflated the Poisson denominator about 24-fold and made the test flag almost every concentrated locus. Omitting the tail makes the test slightly conservative instead.

### Caching

The index is pickled next to the GTF with a key that pins the cache version, the excluded biotypes, and the repeats path:

```python
# annotation.py:655-663
repeats_key = str(repeats_path.resolve()) if repeats_path else ""
key = (
    f"{_CACHE_VERSION}"
    f":{':'.join(sorted(exclude_biotypes))}"
    f":{repeats_key}"
)
digest = hashlib.md5(key.encode()).hexdigest()[:8]
stem = gtf_path.stem.replace(".gtf", "").replace(".gz", "")
return gtf_path.parent / f".scnoisemeter_{stem}_{digest}.cache.pkl.gz"
```

Including the repeats path means adding or removing `--repeats` invalidates the cache, so the repeats layer cannot silently disappear (`annotation.py:651-653`). *To force a rebuild:* `--no-cache`.

Auto-download (`annotation_fetcher.py`) covers the GENCODE GTF (latest or `--gtf-version N`), PolyASite 3.0 (probing GENCODE versions down to v42, `_POLYASITE_MIN_GENCODE`), PolyA_DB v4, FANTOM5 CAGE peaks, and 10x whitelists. Every URL is overridable with an `SCNM_*_URL` environment variable (`_url`, `annotation_fetcher.py:131-145`). The polyA/TSS site dictionaries are cached separately, keyed on path, mtime, size, a hash of the first 64 KB, and chromosome style. A GTF-versus-PolyASite version gap greater than 5 GENCODE releases triggers a warning (`_check_version_consistency`, `cli.py:378-399`), with the fix being either `--gtf-version 42` or `--polya-db polyadb4`.

---

## 7. Metrics and noise definitions

The noise sets are defined in `constants.py:97-128` and selected by stranding in `metrics.py`:

```python
# metrics.py:291-305
_noise_cats = NOISE_CATEGORIES_UNSTRANDED if unstranded else NOISE_CATEGORIES
sm.noise_read_frac = sum(sm.read_fracs.get(cat.value, 0.0) for cat in _noise_cats)
sm.noise_base_frac = sum(sm.base_fracs.get(cat.value, 0.0) for cat in _noise_cats)
_noise_cats_strict = NOISE_CATEGORIES_STRICT_UNSTRANDED if unstranded else NOISE_CATEGORIES_STRICT
sm.noise_read_frac_strict = sum(sm.read_fracs.get(cat.value, 0.0) for cat in _noise_cats_strict)
sm.noise_base_frac_strict = sum(sm.base_fracs.get(cat.value, 0.0) for cat in _noise_cats_strict)
```

The four sets:

| Set | Members | Use |
|---|---|---|
| `NOISE_CATEGORIES` (conservative) | exonic_antisense, intronic_pure, intronic_boundary, intergenic_sparse, intergenic_repeat, intergenic_hotspot, chimeric | Upper bound |
| `NOISE_CATEGORIES_STRICT` | the above minus intronic_pure and intronic_boundary | Lower bound (unambiguous artifacts) |
| `NOISE_CATEGORIES_UNSTRANDED` | conservative minus exonic_antisense | Smart-seq / FLASH-seq |
| `NOISE_CATEGORIES_STRICT_UNSTRANDED` | strict minus exonic_antisense | Smart-seq / FLASH-seq |

The rationale is in the comments (`constants.py:107-128`): intronic-pure and intronic-boundary "cannot be distinguished from genuine pre-mRNA capture at the read level," and antisense in an unstranded library "is genuine cDNA signal." The `unstranded` flag is set when `--platform smartseq`. Note `mitochondrial` and `multimapper` are in no noise set; they are reported as their own fractions.

Other metric formulas:

- **Strand concordance** (`metrics.py:308-310`): `exonic_sense / (exonic_sense + exonic_antisense)`.
- **Per-category fractions**: read and base, per cell and sample-wide.
- **Per-cell table**: only cells with `>= MIN_READS_PER_CELL` (10, `metrics.py:72`).
- **UMI complexity** per cell per category: `unique_UMIs / reads`, disabled by `--no-umi-dedup`.
- **Full-length fraction** (`metrics.py:346-368`): if a polyA-site dict is present, the fraction of exonic-sense 3' ends within `POLYA_SITE_PROXIMITY` (50 bp) of an annotated site; otherwise a length fallback using `FULL_LENGTH_THRESHOLD` (`metrics.py:79-83`, ONT 500, PacBio 1000).
- **TSS-anchored fraction** (5' ends within `TSS_SITE_PROXIMITY`, 100 bp): needs `--tss-sites`.
- **Length stratification**: bins from `LENGTH_BIN_BREAKS`, with the `<150` bin added only when median length `< LENGTH_SHORT_READ_THRESHOLD` (300).
- **Per-cluster** (`compute_cluster_metrics`, requires `--obs-metadata`): median and IQR of per-cell noise and of flag rates per cluster.

---

## 8. Pipeline mechanics

`pipeline.py` runs one worker process per contig (contigs above ~1 Mb, sorted longest-first for load balancing), via a process pool sized by `--threads`. Each worker builds its own `ReadClassifier`, fetches its contig, classifies, and tallies.

**Reservoir sampling** keeps memory bounded for length, insert-size, intergenic, and position samples. The reservoir is a `list` subclass that tracks `n_seen` for a correct Algorithm R (`pipeline.py:64-72`):

```python
# pipeline.py:685-705
def _reservoir_add(reservoir, value, max_size):
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
```

Length reservoirs cap at 5 000 per category (`MAX_LENGTH_SAMPLE`, `pipeline.py:108`); intergenic and insert-size reservoirs at 500 000 (`pipeline.py:80-85`). **Determinism**: with `--seed`, each contig gets its own reproducible seed (`pipeline.py:470-473`):

```python
def _worker_seed(contig):
    if seed is None:
        return None
    return (seed + zlib.adler32(contig.encode("utf-8"))) & 0xFFFFFFFF
```

Same seed and same thread count give byte-identical TSVs. The HTML report is not byte-identical because Plotly emits random div IDs.

**Intergenic collection**: every intergenic read is stored as a 7-tuple with the true `reference_end` (which includes intron spans) and a strand-correct 3' end, for the second pass. **Barcode normalization**: the `--cell-barcodes` called-cell filter strips a trailing `-1` from both the tag and the list before comparing, and drops uncalled reads entirely (they never reach the classifier). **UMI handling**: UMIs are collected into per-cell-per-category sets for complexity; no read-count deduplication is performed at classification time.

---

## 9. Platform and chemistry auto-detection

`bam_inspector.py` reads the `@PG` records. The primary aligner is matched against `ALIGNER_PLATFORM_MAP` (`constants.py:202-209`): minimap2 to ONT, pbmm2 to PacBio, STAR / STARsolo / cellranger to 10x, bwa to unknown (`bam_inspector.py:207-220`). It then refines from pipeline hints in `@PG` ID/PN/CL fields (`bam_inspector.py:246-288`): ONT tokens (`wf-single-cell`, `guppy`, `dorado`, `bonito`), PacBio tokens (`isoseq3`, `pbmm2`, `skera`, `lima`, `ccs`, `kinnex`), and a 10x-versus-BD split. One consequential branch:

```python
# bam_inspector.py (platform refinement)
if has_sc:                       # cellranger / starsolo / spaceranger present
    ... ILLUMINA_10X or ILLUMINA_BD
if has_star:                     # bare STAR, no 10x-specific tool
    meta.platform = Platform.SMARTSEQ
```

So a Smart-seq / FLASH-seq plate aligned with bare STAR can be inferred, but a STAR BAM that also carries Cell Ranger or STARsolo signatures is treated as 10x. In practice, **set `--platform smartseq` explicitly for plate data**; it is the safe path, because it also selects the unstranded noise sets and suppresses the missing-`CB` warning.

Pipeline stage is inferred from the same hints (`_infer_pipeline_stage`, `bam_inspector.py:291-322`), and barcode presence is sampled over the first 10 000 reads with a 50% threshold (`_detect_barcode_tags`, `bam_inspector.py:325-389`). Chemistry detection (10x v3/v4 16 bp, BD Rhapsody 27 bp) lives in `discover_inspector.py:147-209`. Smart-seq plate wells follow the `<PlateID>_<WellID>` convention with a well regex `[A-Pa-p]\d{1,2}` (`sample_sheet.py:79`), covering 96-well (A-H) and 384-well (A-P).

---

## 10. CLI surface and subcommands

Four subcommands, defined in `cli.py`. Shared options come from the `@_shared_options` decorator (`cli.py:92-193`).

| Subcommand | Purpose | Key required flags |
|---|---|---|
| `run` | One BAM, all outputs | `--bam`, `--output-dir` |
| `run-plate` | A Smart-seq / FLASH-seq plate, aggregated | `--plate-dir`, `--output-dir` |
| `compare` | Two BAMs, per-category chi-squared | `--bam-a`, `--bam-b`, `--output-dir` |
| `discover` | Scan a directory, infer, batch-run | `--bam-dir`, `--reference`, `--output-dir` |

Selected shared flags and defaults (full list at `cli.py:92-207`): `--platform auto`, `--pipeline-stage auto`, `--chimeric-distance 10000`, `--tso` (replaces defaults), `--tso-min-match 12`, `--no-polyg-tso`, `--threads 4`, `--no-umi-dedup`, `--reference`, `--repeats`, `--gtf` / `--gtf-version`, `--polya-db polyasite3`, `--tss-db fantom5`, `--seed`, `--no-cache`, `--offline`, `--obs-metadata`.

The `compare` subcommand uses a chi-squared contingency test per category with a Bonferroni correction (`_write_compare_outputs`, `cli.py:2566-2637`). Its own docstring flags the caveat: it is not a paired test, and BAM B is typically a subset of BAM A (post-filter is a subset of pre-filter), so the same reads appear on both sides and the independence assumption does not strictly hold. `--tso-a` / `--tso-b` let the two sides use different TSOs.

---

## 11. Output files

Verified against the example outputs under `results_examples/`.

**`read_metrics.tsv`** (two columns, metric and value): `n_reads_total`, `n_reads_classified`, `n_reads_unassigned`, `n_cells`, `noise_read_frac`, `noise_base_frac`, `noise_read_frac_strict`, `noise_base_frac_strict`, `strand_concordance`, `chimeric_read_frac`, `multimapper_read_frac`, `per_cell_noise_median`, `per_cell_noise_iqr`, the artifact-flag counts, `read_frac_<category>` and `base_frac_<category>` for all 16 categories, and `full_length_read_frac`.

**`cell_metrics.tsv`** (one row per cell with >= 10 reads): `cell_barcode`, `n_reads`, `n_bases`, the per-category read/base/UMI-complexity columns, `noise_read_frac`, `noise_base_frac`, and per-cell flag counts `n_tso`, `n_polya`, `n_noncanon`, `n_tso_concat`.

**`intergenic_loci.tsv`**: `contig`, `start`, `end`, `strand`, `n_reads`, `n_barcodes`, `has_splice_evidence`, `is_monoexonic`, `polya_run_downstream`, `near_polya_site`, `poisson_pvalue_adj`, `category`.

**`_length_stratified.tsv`**: `length_bin`, `category`, `count`, `fraction_of_bin`, `fraction_of_total`.

**`multiqc.json`**: `{id: "scnoisemeter", data: {sample: {...}}}` with the noise fractions, strand concordance, chimeric/multimapper rates, per-cell summaries, artifact-flag fractions (`tso_invasion_frac`, `polya_priming_frac`, `tso_concatemer_frac`), and the per-category read fractions.

**`report.html`**: metadata table, noise donut, per-category read and base bars, length-stratified and length-distribution charts (long-read), per-cell noise violin, artifact-flag bars, insert-size distribution (paired-end), intergenic loci bar and scatter, and per-cluster plots when `--obs-metadata` is given.

---

## 12. Comparison to existing tools

scNoiseMeter is a per-alignment, mechanistic noise classifier. It is not an ambient-RNA remover, an empty-droplet caller, a doublet detector, or an isoform-QC tool. The comparison below is by function.

### Short-read scRNA-seq

- **Cell Ranger / STARsolo built-in metrics.** Report four read classes (exonic / intronic / intergenic / antisense, by a 50%-overlap rule) as sample-level fractions, plus "reads mapped confidently to the transcriptome." scNoiseMeter splits these into 16 mechanistic categories, reports them per cell and per cluster, scores intergenic loci against a null model, and adds the TSO and internal-priming flags the aligners do not produce. (https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-gex-algorithm)
- **RNA-SeQC, Picard CollectRnaSeqMetrics.** Bulk tools that count *bases* by region (exonic/intronic/intergenic/UTR) and report strand specificity. scNoiseMeter is single-cell-aware and read-resolved, and counts both reads and bases. (https://academic.oup.com/bioinformatics/article/28/11/1530/267467 , https://gatk.broadinstitute.org/hc/en-us/articles/360037057492-CollectRnaSeqMetrics-Picard)
- **SoupX, DecontX, CellBender.** Estimate and remove ambient/background RNA, returning a single contamination fraction per cell and a corrected count matrix. scNoiseMeter does not remove anything; it classifies the *mechanism* (internal priming vs TSO invasion vs intergenic) from the BAM, not the matrix. Use them together: scNoiseMeter to diagnose, these to correct. (https://academic.oup.com/gigascience/article/9/12/giaa151/6049831 , https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6 , https://www.nature.com/articles/s41592-023-01943-7)
- **DropletUtils emptyDrops / barcodeRanks.** Decide which barcodes are real cells against a Dirichlet-multinomial ambient null. That is upstream of scNoiseMeter, which assumes droplets are already called and asks what noise is inside the real cells. (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y)
- **Scrublet, DoubletFinder.** Detect doublet cells from the count matrix. Orthogonal to scNoiseMeter, which works on the BAM at read level; a high inter-chromosomal chimera rate can hint at doublets but is not a doublet call. (https://www.sciencedirect.com/science/article/pii/S2405471218304745 , https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30073-0)
- **Internal-priming literature.** No standard automated tool comprehensively flags internal oligo-dT priming in 3' scRNA-seq; APA tools (Sierra, scAPA, SCAPTURE) infer polyA sites rather than flag artifacts. scNoiseMeter's internal-priming flag and `intergenic_hotspot` rule address this directly at read and locus level. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9142200/)

### Long-read scRNA-seq

- **SQANTI3, pigeon.** Classify *assembled isoforms* (FSM/ISM/NIC/NNC) and flag intra-priming and RT-switching post-assembly. SQANTI3's intra-priming test is a fixed downstream-A-richness threshold. scNoiseMeter works pre-assembly on individual alignments and has its own internal-priming and TSO logic; it feeds cleaner reads into these tools rather than replacing them. (https://www.nature.com/articles/s41592-024-02229-2)
- **IsoSeq, skera, lima.** PacBio segmentation and demultiplexing for Kinnex. They produce reads; scNoiseMeter classifies them. Note: scNoiseMeter has no dedicated skera-boundary detector (the `KINNEX_ADAPTER` constant is defined but unused, see [section 14](#14-discrepancies-found-during-review)); imperfectly segmented reads surface as chimeras or TSO concatemers, not as a Kinnex-specific category.
- **FLAMES, FLAIR, IsoQuant, TALON.** End-to-end isoform discovery and quantification with post-assembly QC. Different operational level; scNoiseMeter is read-level and upstream. (https://www.nature.com/articles/s41467-024-48117-3)
- **LongQC, NanoPlot, pycoQC.** Run-level QC (yield, length, quality). Orthogonal: they answer "is this run clean?", scNoiseMeter answers "is this read a technical artifact, and of what kind?" (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7144081/)

---

## 13. Strengths and limitations

**Strengths**

- One taxonomy across long-read, short-read, and plate-based assays, with platform differences switched automatically.
- Per-read and per-base resolution, aggregated to per-cell and per-cluster.
- A Poisson model for intergenic loci instead of a flat "intergenic equals noise" rule.
- Configurable TSO and internal-priming detection, including the TSO-concatemer metric from Chou et al.
- Two principled noise bounds (conservative and strict) rather than a single fragile number.
- Deterministic with `--seed`; annotation auto-downloaded and cached.

**Limitations**

- It classifies; it does not correct or remove reads. Pair it with SoupX/DecontX/CellBender for correction.
- `intronic_pure` and `intronic_boundary` cannot be separated from genuine pre-mRNA at the read level, which is why they sit in conservative but not strict noise.
- Chimeric calls mix technical chimeras and real gene fusions; interpretation needs context.
- The intergenic Poisson model assumes a homogeneous background, an approximation near highly transcribed regions.
- The junction-canonicality check inspects only the donor dinucleotide, not the acceptor (`classifier.py:922-923`).
- Human GRCh38/hg38 only; chromosome naming must match between BAM and GTF or the run aborts.
- Not benchmarked against a sorted or spike-in ground truth, a gap shared across this tool class.
- Non-MMLV chemistries need explicit configuration. A poly-T-tailed TSO (UltraMarathonRT) will not be caught by the default 10x/PacBio sequences or the poly-G heuristic; pass `--tso` with the correct sequence, and note that the internal-priming A-run heuristic is calibrated for genomic A-runs, not poly-T tails.

---

## 14. Discrepancies found during review

Surfaced here rather than smoothed over.

1. **Category priority: chimeric vs mitochondrial.** `classify()` checks mitochondrial (`classifier.py:330-336`) before chimeric (`classifier.py:338-344`), so a mitochondrial read carrying a chimeric `SA` tag is labeled `mitochondrial`. The README "Read categories" table and this file's own top docstring (`classifier.py:8`) list chimeric before mitochondrial. The effect is rare but the documented order does not match the code. Fix the docs, or reorder the checks if chimeric should win.

2. **Intergenic hotspot polyA-distance rule not enforced.** The README (and the module docstring `intergenic_profiler.py:30-36`) state a hotspot must be `> 50 bp from any annotated polyA site`. `_is_hotspot` (`intergenic_profiler.py:367-379`) only checks `is_monoexonic and polya_run_downstream`; `near_polya` is computed but not used. The function docstring admits it is "not strictly required." Either wire `near_polya` into the rule or correct both docstrings.

3. **README is missing three v0.6.0 features.** Reverse-complement TSO matching, strand-aware internal priming, and the TSO-concatemer metric are all in `docs/documentation.md` but not in `README.md`. Users who read only the README will not know about them.

4. **`KINNEX_ADAPTER` is dead code.** Defined at `constants.py:254` with a comment about skera segmentation, but never referenced in any module. There is no Kinnex-specific artifact detector. Remove the constant or implement the detector; until then, do not claim Kinnex-segmentation detection.

5. **Duplicate full-length thresholds.** `constants.py:337-338` defines `FULL_LENGTH_MIN_LENGTH_ONT`/`_PACBIO` (500 / 1000), but `metrics.py` uses its own local `FULL_LENGTH_THRESHOLD` dict (`metrics.py:79-83`) with the same numbers. Editing the constants.py values will not change behavior; edit `metrics.py:79-83`. Consolidate to one source.

6. **Example outputs are stale.** The bundled `results_examples/` files predate `n_tso_concatemer` and do not contain it, and their `multiqc.json` lacks the artifact-flag fractions that the current `to_multiqc_json` emits. Regenerate the examples with v0.6.0.

None of these change the core classification result for a typical run; items 1, 2, 4 are documentation-versus-code mismatches and items 5, 6 are maintenance issues.

---

## 15. Tuning cheat-sheet

| You want to | Change |
|---|---|
| Allow longer legitimate splits before calling chimeric | `--chimeric-distance` (default 10000) |
| Use a non-10x TSO | `--tso SEQ` (repeatable; replaces defaults) |
| Loosen/tighten TSO matching | `--tso-min-match` (default 12) |
| Stop counting poly-G as TSO | `--no-polyg-tso` |
| Enable internal-priming and non-canonical-junction flags | `--reference FASTA` |
| Treat data as unstranded (Smart-seq / FLASH-seq) | `--platform smartseq` |
| Add per-cluster noise profiles | `--obs-metadata cells.tsv` |
| Change intergenic significance (reads / barcodes / p / window) | `constants.py:285-295` |
| Change internal-priming A-run length or window | `constants.py:304-305` |
| Change polyA-site or TSS proximity | `constants.py:309` , `constants.py:326` |
| Change the full-length length fallback | `metrics.py:79-83` (not constants.py) |
| Change the per-cell minimum read count | `metrics.py:72` |
| Exclude gene biotypes from the index | `--exclude-biotypes` (repeatable) |
| Reproducible sampling | `--seed INT` (same thread count) |
| Match the polyA atlas GENCODE version | `--gtf-version 42` or `--polya-db polyadb4` |
