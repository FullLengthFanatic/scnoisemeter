# Changelog

Metric definitions are versioned deliberately: a change in how a flag is
detected can move a reported value by orders of magnitude, so `cohort` refuses
to pool result sets silently across versions. Any release that alters a
reported number says so here.

## 0.9.1 - unreleased

### Fixed

- **`run` stopped on almost every barcode-aware BAM.** The 0.9.0 invariant
  check read NaN as a value outside [0, 1]. A category that no cell has reads
  in has no UMI diversity, which is NaN: `multimapper` on any `minimap2` or
  `pbmm2` BAM, `unassigned` without a whitelist. The run ended with
  "structural invariant(s) of this result set do not hold" and wrote only
  `qc.json`. NaN is now treated as an unmeasured value, like a missing one.
  The same applied to `strand_concordance` in a sample with no exonic reads.
  A new test runs `scnoisemeter run` through the real command line, which
  nothing did before.
- **`--cell-barcodes` matching nothing finished with all-zero outputs.** The
  guard compared a counter that includes filtered reads, so it never fired.
  The run now stops with the intended message.
- **`three_prime_at_pas_site_frac` could go missing without a message.**
  `discover` and parallel `run-plate` workers loaded the polyA atlas without
  its hexamer subset and cached the result, and every later `run` was served
  that cache. The subset is now always collected and old caches are ignored.
- **`run-plate` accepted settings it did not apply.** With
  `--parallel-wells`, `--priming-criterion` fell back to `a_run` and the
  workers loaded the atlas unfiltered; `--anchor-veto-internal-priming` and the
  PAS-restricted anchoring were ignored on both paths. `run_info.json` recorded
  all of them as applied.
- `n_polya_tail`, `n_rt_switch_repeat` and `three_prime_at_pas_site_frac` were
  computed but never written to `read_metrics.tsv`.
- `run_info.json` did not contain the `priming` block that 0.9.0 documented,
  and `cohort` did not compare it. Both now happen, and `run_info.json` also
  records `intergenic_reclassification_sampled`.
- `discover` declared `--tss-sites` twice, which printed a warning on every
  call.
- The report footer and the download User-Agent named a repository URL that is
  not this project's.
- Tests wrote site-atlas caches into the real `~/.cache/scnoisemeter`;
  `tests/test_integration.py` failed to import once a dependency shipped a
  top-level `tests` package; the benchmark runner used whichever `python3` was
  first on `PATH`, and one of its ground-truth labels predated 0.7.

### Reported values that change

- **Chimeric calls on the minus strand.** For a split alignment on the same
  contig and strand, the query-order test assumed read orientation while both
  CIGARs are written in reference orientation. On the minus strand a colinear
  split (a long intron the aligner broke in two) was called `chimeric` and an
  order-inverted one was not. Verified against `minimap2` output. On chr21 of
  the Kinnex example, 591 minus-strand reads lose the call and 1,076 gain it
  (of 7,744 with an `SA` tag); the plus strand is unchanged, and the chimeric
  share of that contig moves from 3.49% to 3.58%. `chimeric_read_frac` and
  `artifact_candidate_*` move accordingly.
- **TSO poly-G shortcut on reverse alignments.** A reverse alignment stores
  the read reverse-complemented, so its poly-G is a poly-C in the clip. The
  shortcut looked for G on both strands. `n_tso_invasion` changes for
  reverse-strand reads wherever the shortcut is on (the built-in 10x motif).
- **`strand_asymmetry_pvalue` is two-sided.** The dominant strand is chosen
  after counting, so the one-sided tail reported in 0.9.0 was half the correct
  value: 3 of 4 reads is p = 0.625, not 0.31. The column still gates nothing.
- **Artifact-flag rates in `multiqc.json`** keep eight decimals. Four decimals
  wrote a measured 4.5e-05 (1,593 reads) as 0.0.

### Corrected

Numbers quoted in the 0.9.0 documentation, help text and the entry below,
re-derived from the PolyASite 3.0 GENCODE 42 atlas and GRCh38:

- Share of the genome within 50 bp of an atlas site: 37.4% unfiltered, 6.3% at
  `--polya-min-support 25`, 2.4% at 50 (was 38.6%, 7.6%, 3.0%). These are the
  values the tool itself measures and logs.
- Median A-content of the 20 bp downstream of a site: 30% at annotated polyA
  sites and 30% at random genomic positions (was 25% and 25%).
- `a_run` and `perc_a` agree on 97% of random genomic windows and 95% of
  atlas sites (was "~97%").
- PolyASite 3.0 holds 18.4 M records; the `--polya-db` help said ~569 k.
- The README described the report screenshots as a six-contig subset. They
  show the whole sample (37.4 M reads), run with 0.8.0.
- White paper: the worked read is 940 bp long and spans 1,240 bp; removing its
  `N` leaves it `exonic_sense` because a read takes the category holding most
  of its aligned bases, which the category table now says instead of
  "≥1 aligned base"; `exonic_sense` in the example run is 64.76%.

## 0.9.0 - 2026-09-10

Two settings
below change reported values; both default to 0.8.1 behaviour, so installing
this release changes no number until a flag is passed. One correction is not
behind a flag and does move numbers -- see "Reported values that change".

### Added -- new intergenic columns, promotion unchanged

- `n_molecules` per intergenic window: distinct `(cell barcode, UMI)` pairs.
  `n_reads` counts PCR copies, so a single over-amplified molecule can look
  like broad support. Counted exactly, with no mismatch collapsing -- whether
  a 1-mismatch neighbour is a sequencing error or a distinct molecule depends
  on UMI length and the library error rate, which are properties of the
  chemistry rather than of this code, so the exact count is an upper bound and
  error correction could only lower it. `-1` when no UMIs were tracked, which
  is not the same claim as zero molecules.
- `strand_asymmetry_pvalue` per intergenic window: the binomial tail
  `P(X >= dominant | n, 0.5)`. `strand_fraction` alone is depth-blind -- 3 of 4
  reads and 75 of 100 are both 0.75, and the existing consistency threshold
  accepts them equally, though the first is what a coin produces one time in
  three. `NaN` for unstranded libraries, where reads divide near 50/50 by
  design and the test would measure the protocol rather than the locus.

Neither column gates promotion in this release. Both change the same
quantities as each other and as the open question about the significance
procedure, so they ship as evidence first: the choice between the global
Poisson model and the alternatives is better made from these columns measured
on real samples than from argument.

### Not shipped, and why

- A wrong-strand **decoy null** for intergenic enrichment was implemented and
  rejected on measurement. Within any window the minority-strand count cannot
  exceed the dominant one, so a decoy distribution built from minority counts
  sits below the target distribution even under pure noise; on simulated
  windows that were entirely background the estimator returned an "FDR" of
  0.51. A biased number that looks calibrated is worse than no number.
  `strand_asymmetry_pvalue` answers the part of that question that has an
  exact null. A regression test records the property that makes the decoy
  unsound, so it is not reintroduced.

### Added — new evidence, off by default where it changes a number

Transferred from `irinakhven/flightcollapse` after a side-by-side screening of
both tools' 3'-end handling. Every setting that moves a previously reported
value defaults to 0.8.1 behaviour, so installing this changes nothing until a
flag is passed. All of them are recorded under a `priming` key in
`run_info.json` so `cohort` can refuse to pool runs that disagree.

- `--polya-min-support` filters the polyA atlas on its own per-site support
  column. Measured against PolyASite 3.0: 18,432,135 records put **38.6% of
  the primary genome** within the 50 bp proximity window, which makes
  "near an annotated polyA site" close to a constant. It vetoes
  `intergenic_hotspot` across a third of the genome and satisfies one of the
  two `intergenic_novel` evidence routes for free. 25 gives 7.6%, 50 gives
  3.0%. Every run now measures and reports the density of the atlas it loaded.
- `polya_tail_frac` — an untemplated A-run in the read's own 3' **soft** clip.
  Per-molecule evidence, where the existing internal-priming flag is a
  property of the locus: median downstream A-content is 25% both at annotated
  polyA sites and at random genomic positions. A read carrying both flags is a
  real cleavage site in A-rich context.
- `three_prime_at_pas_site_frac` — anchoring restricted to atlas records
  carrying a polyA signal hexamer, read from column 11 of a PolyASite BED.
- `rt_switch_repeat_frac` — direct repeat shared by an intron's two ends, the
  reverse-transcriptase template-switching signature. Independent of splice
  canonicality: an RT artifact can sit on a GT-AG boundary.
- `local_read_share` on every intergenic window — its share of intergenic read
  mass within 10 kb, separating an isolated locus from one window of a
  diffusely covered neighbourhood.
- `--priming-criterion {a_run,perc_a,either}`. `perc_a` is the SQANTI3/pigeon
  `perc_A_downstream_TTS` convention (>=60% A over 20 bp). The two criteria
  agree on ~97% of windows, so this is for comparability with isoform-QC
  numbers and is explicitly not a correction.
- `--anchor-veto-internal-priming` requires the absence of an internal-priming
  signature for a read to count towards 3'-end anchoring, which lowers
  `three_prime_anchored_frac`.

### Reported values that change

- **The Bonferroni denominator for the intergenic enrichment screen counted
  windows that were never tested.** `n_test_windows` was every 500 bp window of
  every non-mitochondrial contig (~6.18 M on GRCh38), but the screen only ever
  tests the intergenic complement — roughly 2.2–2.8 M windows depending on the
  annotation release. That is a 2.2–2.9x overcorrection, and it suppressed real
  promotions. The denominator is now derived from the same intergenic base
  total already computed for the background rate.

  This raises `intergenic_hotspot`, `intergenic_novel`, `intergenic_repeat` and
  `intergenic_enriched` counts, lowers `intergenic_sparse`, and therefore moves
  `broad_noncanonical_*` and `artifact_candidate_*`. It is a correction rather
  than a policy change — the previous value corrected for tests that were not
  performed — so it is not behind a flag. Result sets from earlier versions are
  not comparable in these categories, which `cohort` already flags by version.

### Added — runtime invariant checks

- Every run now checks the structural properties of its own result set and
  writes `<sample>.qc.json`: each composition sums to one, every fraction lies
  in [0, 1], `artifact_candidate_*` cannot exceed `broad_noncanonical_*`,
  `both_ends_anchored_frac` cannot exceed either single-end fraction, no
  artifact-flag count exceeds `n_reads_classified`, UMI diversity ratios stay
  in [0, 1], and `n_cells` matches the per-cell table. Observed values are
  recorded, not just verdicts. A violation stops the run; `--no-strict` writes
  the outputs and records it. Reservoir-estimated reclassification downgrades
  the conservation checks to warnings, because scaled totals are not expected
  to conserve exactly.

### Changed — runs can now stop instead of reporting a configuration error

- A barcode gate that leaves more than half of classified reads unassigned now
  **stops the run** and prints the observed counts, cell count and header-scan
  tag fraction. Previously the run completed and the failure arrived as a
  result: a sample with a handful of cells and an unremarkable composition.
  `--no-strict` downgrades this to a warning. Barcode-free BAMs are exempt,
  where 100% unassigned is the design.
- A site atlas whose coordinates fall beyond the end of their contigs is
  discarded above 20% out-of-range and warned about above 1%. Contig-name
  agreement is not assembly agreement: an hg19 atlas against an hg38 BAM
  shares every contig name, so the existing naming guard passes while every
  coordinate is shifted and every proximity test silently fails.

### Corrected

- The v0.8.1 commit `c6f50a4` and its source comment stated that PolyASite 3.0
  ships Ensembl-named contigs and that therefore "every lookup missed" for a
  UCSC-named BAM, making the endpoint-anchoring metrics and
  `intergenic_hotspot` wrong for the tool's default configuration. **That was
  incorrect.** `annotation_fetcher` requests PolyASite 3.0 in its GENCODE
  flavour (`atlas.clusters.3.0.GRCh38.GENCODE_NN.bed.gz`), whose contigs are
  `chr`-prefixed, as is the FANTOM5 CAGE atlas. The default auto-download path
  therefore matched a UCSC-named BAM directly and reported correct anchoring
  fractions; the committed report screenshots, from a UCSC-named BAM with both
  atlases auto-downloaded, show 99.62% 3'-end and 84.02% 5'-end anchoring and
  corroborate this.

  The claim was traced to a misreading of a docstring example in
  `annotation_fetcher.extract_polyasite_version_from_filename`, where
  `atlas.clusters.2.0.GRCh38.96.bed.gz -> None (Ensembl, not GENCODE)`
  documents that the *older 2.0 file* yields no GENCODE version number. It does
  not describe what the tool downloads.

  What was genuinely broken, and what the 0.8.1 change actually fixes, is
  narrower:

  - a **user-supplied Ensembl-named BED** (`--polya-sites` / `--tss-sites`) against
    a UCSC-named BAM. The old helper only ever stripped a `chr` prefix, so
    nothing rewrote in that direction and every lookup missed while the dict
    stayed non-empty — the silent-zero failure, but confined to users supplying
    their own atlas;
  - the **mitochondrion under an Ensembl-named BAM**, where stripping turned the
    atlas's `chrM` into `M` while the BAM says `MT`, silently dropping
    mitochondrial polyA/TSS proximity. Fixed in the follow-up commit.

  No reported value changes for the default auto-download path as a result of
  the 0.8.1 atlas work. `_check_site_contig_overlap` remains worthwhile as a
  guard: it turns a non-matching atlas into an absent metric plus a loud error
  rather than a measured zero.

## 0.8.1

### Reported values that change

- `unmapped_read_frac`, `low_mapq_read_frac`, `per_cell_noise_median` and
  `per_cell_noise_iqr` are now **absent** rather than `0.0` when their
  denominator does not exist — no BAM-index counters, no classified
  alignments, or a barcode-free BAM with no cells to summarise. They are
  written as blank in `<sample>.read_metrics.tsv` and omitted from
  `<sample>.multiqc.json`; `cohort` already reads a blank as a gap. A
  genuine measured zero is still reported as `0.0`.

  This is the same principle already applied to endpoint anchoring and `NM`
  summaries: reporting `0.0` for a metric that was never measurable reads as
  "measured, and there were none". Cohort figures built from 0.8.0 output will
  show these as zeros for barcode-free samples where 0.8.1 shows a gap.

No other reported value changes in this release.

### Fixed

- **Annotation index build was ~23x slower than necessary.**
  `_extract_splice_junctions` grouped by `["Chromosome", "Strand",
  "transcript_id"]` with `observed=False`. `pr.read_gtf` returns Chromosome and
  Strand as Categoricals, so that iterated the full category cartesian product:
  2,634,150 groups instead of 35,122 on chr1 alone. A full GENCODE v47 build
  went from 3551 s to 155 s, with identical output (526,926 junctions,
  1,358,754 splice sites, and unchanged interval counts and base extents).
  `_extract_splice_sites` was also converted from `iterrows` (~102k rows/s) to
  grouped set construction.

- **Repeat-overlap scoring re-normalised the whole contig interval list per
  read.** `_repeat_overlap_fraction` called `sorted()` inside the per-record
  loop on a list the caller had already merged and sorted, then scanned from
  index 0, so a read late on a chromosome walked every preceding interval. At
  hg38 chr1 RepeatMasker scale (470k intervals) a single 5-read locus took
  200 ms, which is hours across the 20k-130k windows a real sample produces.
  The union is now merged once per run by `merge_repeat_intervals()` and each
  locus binary-searches it. Verified identical on randomised disjoint,
  overlapping and unsorted inputs.

- **`requires-python` claimed 3.9, which cannot import the package.**
  `classifier.ReadResult` uses `@dataclass(slots=True)` (3.10+). Minimum is now
  3.10 and CI tests 3.10-3.13 instead of 3.9 and 3.12.

- **CI ran an unpinned `ruff`,** so an upstream default-rule change turned the
  lint job red on untouched commits. `ruff` is pinned in the `dev` extra and
  `[tool.ruff.lint]` now selects an explicit rule set (E4/E7/E9/F plus import
  ordering). Cosmetic modernisation rules are deliberately not enabled.

### Changed

- The annotation index is passed to the worker pool once through a `Pool`
  initializer instead of travelling inside every per-contig task dict. It is
  ~75 MB pickled on GENCODE v47 without a repeats BED, so a 25-contig run was
  serialising and re-parsing roughly 1.9 GB. The `pipeline` module docstring
  claimed the index was shared read-only after fork; it was not.
- Removed a `[tool.setuptools.package-data]` entry for a `scnoisemeter.data`
  package that does not exist in the tree.

### Documentation

- README and `methods_annotated.md` now state that `minimap2` and `pbmm2` do
  not emit `NH`, so the `multimapper` category is empty by construction on most
  long-read BAMs and those alignments are counted under their genomic category
  instead. `multimapper_read_frac` is therefore not comparable between a
  long-read and a STAR-derived sample. This is a **documentation-only** change;
  detection is unchanged in 0.8.1.
- README now states that the intergenic Poisson test itself runs on the
  reservoir sample, not only the totals derived from it, so promotion to
  `intergenic_hotspot` / `intergenic_novel` is depth-dependent; and that the
  Bonferroni denominator counts every 500 bp window of every non-mitochondrial
  contig, which is larger than the number of windows actually tested.
- README documents that duplicate- and QC-fail-flagged records are still
  classified, so composition from a `MarkDuplicates` BAM is duplicate-weighted.
- The `annotation` module docstring no longer claims a sub-60-second build.
