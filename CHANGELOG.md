# Changelog

Metric definitions are versioned deliberately: a change in how a flag is
detected can move a reported value by orders of magnitude, so `cohort` refuses
to pool result sets silently across versions. Any release that alters a
reported number says so here.

## Unreleased

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

  - a **user-supplied Ensembl-named BED** (`--polya-bed` / `--tss-bed`) against
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
