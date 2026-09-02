# scNoiseMeter documentation

Version 0.8.0 operational reference.

This file documents installation, inputs, commands, options, outputs, and interpretation. The scientific rationale is in [whitepaper.md](whitepaper.md); the code map is in [methods_annotated.md](methods_annotated.md).

## 1. Requirements

- Python 3.10+
- coordinate-sorted BAM
- BAM index (`.bam.bai` or `.bai`)
- matching GTF annotation
- optional indexed reference FASTA for genomic-context and splice-motif checks

Install:

```bash
pip install git+https://github.com/FullLengthFanatic/scnoisemeter.git
scnoisemeter --version
```

Prepare BAM:

```bash
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam
```

Automatic annotation downloads target human GRCh38. For other organisms, provide all resources explicitly and validate naming/coordinates.

## 2. Minimal run

```bash
scnoisemeter run \
  --bam sample.sorted.bam \
  --output-dir results/
```

Recommended reproducible run:

```bash
scnoisemeter run \
  --bam sample.sorted.bam \
  --sample-name donor_01 \
  --gtf gencode.v45.annotation.gtf.gz \
  --reference GRCh38.primary_assembly.genome.fa \
  --repeats hg38.repeatmasker.bed.gz \
  --polya-sites polyasite.bed.gz \
  --tss-sites fantom5.cage.bed.gz \
  --cell-barcodes barcodes.tsv.gz \
  --platform ont \
  --library-strand stranded \
  --seed 42 \
  --threads 8 \
  --output-dir results/
```

## 3. Input semantics

### BAM records

All references with at least one indexed mapped record are processed, including mitochondrial and short references. Unmapped/secondary/supplementary records are not assigned a genomic category but remain visible in denominators and flag totals. QC-fail and duplicate flags are counted but do not by themselves exclude a mapped primary alignment. BAM index counters are alignment-record counts, not unique templates. Historical `n_reads_*` and `read_frac_*` names are retained, but paired mates count separately; precise `n_records_*` and `n_alignments_classified` aliases are also emitted.

The default tags are:

- corrected cell barcode: `CB`;
- corrected UMI: `UB`;
- number of reported alignments: `NH`.

Override the first two with `--barcode-tag` and `--umi-tag`.

### Barcode whitelist versus called cells

`--barcode-whitelist` defines valid barcode sequences. A primary alignment with a missing/off-list barcode becomes `unassigned` and remains in sample metrics.

`--cell-barcodes` defines called cells. Reads not in this set are dropped and do not affect metrics. Cell Ranger’s `filtered_feature_bc_matrix/barcodes.tsv.gz` is accepted; a trailing `-1` is normalized for matching.

Without a whitelist or CB tag, reads pool under `NO_BARCODE`; sample-level metrics still work, but per-cell metrics do not. Plate mode converts that sentinel to a well identifier before aggregation.

### GTF

`--gtf` accepts plain or gzipped GTF. If omitted, a GENCODE annotation is cached automatically. `--gtf-version N` selects a release. Chromosome naming must match the BAM (`chr1` versus `1`). Annotation cache identity includes metadata and a content-head hash.

### Reference FASTA

`--reference` enables:

- strand-aware internal polyA-context flags;
- donor-plus-acceptor splice-motif checks;
- intergenic hotspot context.

The FASTA must match the BAM assembly and have an `.fai` index.

### Repeat, polyA, TSS, and NUMT resources

- `--repeats`: BED3+ repeat intervals used for `intergenic_repeat`.
- `--polya-sites`: repeatable BED3+; strand in column 6 is preserved.
- `--tss-sites`: repeatable BED3+; strand in column 6 is preserved.
- `--numt-bed`: records the loaded NUMT interval count as provenance only. It does not produce a NUMT read fraction.

BED coordinates are zero-based, half-open; midpoint is used for end-site resources. Unknown strand (`.` or missing column 6) can match either read orientation.

## 4. Common options

| Option | Default | Meaning |
|---|---|---|
| `--platform` | `auto` | `ont`, `pacbio`, `illumina`, `illumina_10x`, `illumina_bd`, `smartseq`, or `unknown` |
| `--library-strand` | `auto` | `stranded` or `unstranded`; explicit is preferred |
| `--pipeline-stage` | `auto` | `raw`, `pre_filter`, `post_filter`, or `custom` metadata |
| `--threads` | 4 | chromosome workers within a BAM |
| `--seed` | 42 | reservoir sampling is seeded by default, so runs agree given identical threads/input; pass another value to draw a different sample |
| `--no-cache` | off | rebuild annotation index |
| `--offline` | off | forbid downloads and require cached/provided resources |
| `--exclude-biotypes` | none | repeatable gene biotype exclusion |
| `--obs-metadata` | none | per-cell cluster metadata for cluster summaries |
| `--no-umi-tracking` | off | disable UMI-string sets to reduce memory; `--no-umi-dedup` is a legacy alias |

`--chimeric-distance` is accepted only for backward compatibility. Same-strand genomic distance is not used to call a chimera.

### Platform inference

Auto-detection uses header pipeline hints. It intentionally does not assert:

- minimap2 = ONT;
- bare STAR = Smart-seq.

Minimap2 alone becomes `unknown`; bare STAR becomes generic `illumina`. Explicit platform and strandedness avoid silent assay assumptions.

### TSO options

`--tso SEQUENCE` is repeatable and replaces the protocol-aware motif default. Only A/C/G/T/N is accepted. `--tso-min-match` sets the prefix length. `--no-polyg-tso` disables the 10x-specific poly-G shortcut.

CLI defaults:

| Platform context | Motif default |
|---|---|
| ONT or 10x | 10x TSO |
| PacBio or Smart-seq | SMART/PacBio TSO |
| Generic Illumina, BD, unknown | none; provide `--tso` if applicable |

The classifier API retains historical behavior when `tso_sequences=None`; pass `[]` to disable sequence matching programmatically.

## 5. `run`

```bash
scnoisemeter run --bam INPUT --output-dir DIRECTORY [COMMON OPTIONS]
```

Additional options:

- `--sample-name`: output prefix;
- `--tagged-bam PATH`: optional indexed BAM copy with final category in local SAM tag `sn`;
- `--cell-barcodes`: called-cell filter;
- `--barcode-whitelist` or `--whitelist-db`;
- `--barcode-tag`, `--umi-tag`;
- `--chemistry`: barcode-length expectation metadata.

The command validates coordinate sort order, BAM index counts, chromosome naming, annotation/reference compatibility, and low read counts before writing the report.

`--tagged-bam` preserves every input record and tags only classified mapped-primary alignments. It retains exact read assignments in memory and performs a second input pass. Final intergenic tags are looked up from the fixed-window locus table, including when aggregate reclassification used a reservoir.

## 6. `compare`

```bash
scnoisemeter compare \
  --bam-a raw.bam \
  --bam-b filtered.bam \
  --label-a raw \
  --label-b filtered \
  --output-dir comparison/
```

`--tso-a` and `--tso-b` override shared `--tso` independently.

Comparison outputs are dependent-data aware:

- stable read keys are query name plus `/1`, `/2`, or `/0`;
- retention is exact set intersection by category in A;
- transition counts use shared keys;
- sample fractions are descriptive;
- paired-cell median differences use cells shared by both BAMs and a fixed-seed 1,000-resample percentile bootstrap.

No independent read-count chi-square p-values are produced.

Retention, transitions and the matching summary are only written when the two BAMs really are
the same reads. Before classifying, `compare` samples read names from a shared genomic window of
both BAMs; below 50% overlap it treats them as independent, skips those three files, and warns.
That is also what keeps memory bounded: tracking per-read identity costs roughly 175 bytes per
read on each side, so a 76M-read and a 37M-read sample together approach 20 GB held at once.

Two different situations produce a low overlap, and the warning covers both. The samples may
genuinely be unrelated, in which case `cohort` is the right tool. Or a step between them may have
rewritten the read names: `isoseq dedup` replaces PacBio CCS names with `molecule/N`, so a real
pre/post-dedup pair shares no names at all. Composition metrics stay valid either way.

Override the detection with `--matched-reads` or `--no-matched-reads`.

## 6a. `cohort`

```bash
scnoisemeter cohort \
  --results results/BD46/ \
  --results results/10x_FL/ \
  --results results/PIPseq/ \
  --sample-sheet cohort.tsv \
  --output-dir cohort/
```

Compares any number of independent samples from directories `run` or `run-plate` already wrote.
No BAMs are read and nothing is re-classified, so a cohort of a dozen samples completes in
seconds. At least two `--results` directories are required.

`--sample-sheet` is an optional TSV or CSV with columns `sample`, `label`, `group` and `order`.
`sample` matches the `<sample>.read_metrics.tsv` stem; `group` drives colour in the artifact and
per-cell figures; `order` overrides the default ordering, which is cleanest sample first by
artifact-candidate composition, or by broad composition when some samples predate that metric.

Three properties of the reader matter when interpreting the output:

- A metric a sample never reported is absent, not zero. `n_tso_concatemer` does not exist before
  v0.7, and drawing it as zero would assert a measurement that was never made.
- Samples from different scnoisemeter versions or different GENCODE releases are flagged. Metric
  definitions changed between releases: the same BD46 BAM reports a TSO invasion rate of
  2.33e-03 under v0.6.1 and 3.28e-06 under v0.7.2, a 700-fold difference that is entirely
  methodological.
- Stranded and unstranded protocols are flagged and the affected heatmap cells are outlined,
  because exonic antisense and strand concordance do not mean the same thing across that boundary.

Samples classified in barcode-agnostic mode have no per-cell values; they appear in every figure
except the per-cell one, which names its exclusions.

## 7. `discover`

```bash
scnoisemeter discover \
  --bam-dir /data/bams \
  --reference GRCh38.fa \
  --output-dir batch_results/ \
  --run-all
```

Without `--run-all`, discovery displays inferred metadata and prompts for a subset. The annotation and site resources are reused across BAMs. Inference warnings should be resolved with explicit settings when protocol details are known.

## 8. `run-plate`

```bash
scnoisemeter run-plate \
  --plate-dir /data/plate \
  --sample-sheet plate.csv \
  --platform smartseq \
  --library-strand unstranded \
  --parallel-wells 8 \
  --threads 16 \
  --output-dir plate_results/
```

Expected well layout is `<PlateID>_<WellID>` with rows A–H for 96-well and A–P for 384-well plates. The sample sheet can add index and arbitrary metadata. Missing BAM indexes are reported and skipped.

`--parallel-wells` controls concurrent wells; total `--threads` is divided across workers. If a worker is killed, the command reports progress and suggests reducing parallelism.

Barcode-free well results are relabeled to `<PlateID>_<WellID>` before merging. Counts, flags, alignment summaries, UMI sets, intergenic records, endpoints, and comparison assignments follow that relabel.

## 9. Output denominators

| Quantity | Denominator |
|---|---|
| `read_frac_<category>` | classified mapped-primary alignments |
| `base_frac_<category>` | classified aligned reference bases |
| `unmapped_read_frac` | mapped + unmapped BAM-index records; absent when the index carries no counters |
| artifact flag fractions in MultiQC | classified mapped-primary alignments |
| per-cell fractions | classified reads/bases for that retained cell |
| endpoint fractions | reservoir-sampled exonic-sense reads with the required atlas |
| intergenic Poisson rate | sampled intergenic records / exact non-mitochondrial BAM-contig complement of GTF gene bodies |

`n_reads_processed` can include records later skipped by classification or called-cell filtering. Use output denominator fields, not log progress counts, for reporting.

## 10. `read_metrics.tsv`

The two-column sample table contains:

### Record counts

- `n_records_total`, `n_records_mapped`, `n_records_unmapped`, `n_alignments_classified`;
- historical aliases `n_reads_total`, `n_reads_mapped`, `n_reads_unmapped`;
- `n_primary_mapped`, `n_secondary`, `n_supplementary`;
- `n_qcfail`, `n_duplicate`;
- `n_reads_classified`, `n_reads_unassigned`, `n_cells`.

### Aggregate composition

- `broad_noncanonical_read_frac`, `broad_noncanonical_base_frac`;
- `artifact_candidate_read_frac`, `artifact_candidate_base_frac`.

The pre-v0.7 aliases `noise_read_frac`, `noise_base_frac`, `noise_read_frac_strict` and
`noise_base_frac_strict` were removed in 0.8.0. `cohort` still reads them from older result
directories, so archived runs remain usable, but nothing writes them any more.

### Alignment/location

- `strand_concordance` = sense / (sense + antisense), missing when neither exists;
- `chimeric_read_frac`, `multimapper_read_frac`, `unmapped_read_frac`;
- `low_mapq_read_frac`, `mean_mapq`, optional `mean_edit_distance`, optional `softclip_base_frac`;
- every `read_frac_<category>` and `base_frac_<category>`.

### Flags

- `n_tso_invasion`;
- `n_tso_concatemer`;
- `n_polya_priming`;
- `n_noncanon_junction`;
- `n_discordant_pair`;
- `n_low_mapq`.

### Ends and provenance

- optional `three_prime_anchored_frac`;
- optional `five_prime_anchored_frac`;
- optional `both_ends_anchored_frac`;
- deprecated `full_length_read_frac` equal to both-end anchoring;
- optional `n_numt_intervals_loaded`.

## 11. `cell_metrics.tsv`

Only non-sentinel cells with at least 10 classified reads are included. Columns include:

- `n_reads`, `n_bases`;
- category read/base fractions;
- broad and artifact-candidate fractions plus legacy aliases;
- flag counts and low-MAPQ rate;
- mean MAPQ, optional mean edit distance, softclip-base fraction;
- `umi_sequence_diversity_<category>`;
- deprecated identical `umi_complexity_<category>` aliases.

UMI sequence diversity is not gene-aware molecule complexity. Missing values can mean no reads/UMIs, tracking disabled, or exact membership unavailable after sampled intergenic reclassification.

## 12. `intergenic_loci.tsv`

Columns:

| Column | Meaning |
|---|---|
| `contig`, `start`, `end` | fixed 500 bp window |
| `strand` | dominant strand |
| `strand_fraction` | fraction on dominant strand |
| `n_reads` | sampled records in window |
| `n_barcodes` | real barcodes, excluding sentinel/missing |
| `has_splice_evidence` | any record has CIGAR `N` |
| `is_monoexonic` | no record has CIGAR `N` |
| `polya_run_downstream` | strand-aware genomic A/T-run evidence |
| `near_polya_site` | dominant-strand site proximity |
| `repeat_overlap_fraction` | aligned-span fraction overlapping repeats |
| `poisson_pvalue` | raw global-rate tail probability |
| `poisson_pvalue_adj` | Bonferroni result across all possible windows |
| `category` | final intergenic interpretation |

This table contains sampled support if more than 500,000 intergenic records were seen.

## 13. Comparison files

### `comparison.metrics.tsv`

Scalar metric values for A and B plus B−A.

### `comparison.stats.tsv`

For each category: sample fraction A/B, fraction delta, paired-cell count, paired-cell median delta, and 95% bootstrap interval.

### `comparison.retention.tsv`

For each A category: initial keys, keys present in B, retained in the same category, retained in another category, and retention fractions.

### `comparison.transitions.tsv`

Counts of category A → category B among shared keys.

### `comparison.matching.tsv`

Overall A/B key counts, shared keys, A-only, and B-only.

`comparison.retention.tsv`, `comparison.transitions.tsv` and `comparison.matching.tsv` exist only
for a nested pair. When the read names do not match they are not written at all, rather than
written full of NaN.

## 13a. Provenance and cohort files

### `<sample>.run_info.json`

Written by `run` and `run-plate` alongside the metrics: scnoisemeter version, sample name, BAM
path, platform, pipeline stage, aligner, strandedness, read and cell totals, and the GTF, polyA,
TSS and TSO sources with how each was obtained. `read_metrics.tsv` holds only numbers, so this is
what lets a later reader, or `cohort`, establish whether two result sets are comparable.

### `cohort.summary.tsv`

One row per sample: identity, group, platform, tool version, strandedness, read and cell counts,
aggregate compositions, strand concordance, chimeric rate, and each artifact-flag rate. Values a
sample never reported are empty rather than zero.

### `cohort.composition.tsv`

Samples by categories, read and base fractions, one row per sample.

### `cohort.report.html`

Composition bars excluding exonic sense, the deviation-from-median category heatmap, artifact
flag rates on a log axis, and per-cell spread for the samples that have barcodes.

## 14. Caching and reproducibility

Default cache directory is `~/.cache/scnoisemeter`. Annotation, polyA, and TSS caches are separately fingerprinted. `--offline` prevents network access. `--no-cache` affects annotation-index reuse.

`--seed` controls read-length, insert-size, endpoint, and intergenic reservoir selection, and defaults to 42. Reservoirs are merged using seen-record weights. Identical seed, inputs, software, and thread count make tabular sampling reproducible; generated HTML may contain non-deterministic Plotly element identifiers.

Seeding is on by default because several reported values derive from reservoir samples: before 0.8.0 an unseeded run moved `broad_noncanonical_read_frac` in the third decimal between identical invocations, which is larger than many of the differences a comparison is trying to resolve. A residual difference around the seventeenth significant digit remains, from the order in which parallel per-contig results are summed, and is not meaningful.

For archival work, preserve explicit input files or checksums, command line, version/commit, seed, and output denominators.

## 15. Interpretation checklist

Before calling a preparation “noisy,” ask:

1. Was `--library-strand` correct?
2. Did BAM and GTF use the same assembly/naming?
3. Were called cells supplied consistently?
4. Was the actual TSO configured?
5. Were repeat and reference resources present?
6. Is the aggregate driven by intronic biology, unresolved intergenic enrichment, or positive artifact evidence?
7. For comparison, were enough read keys and cells shared?
8. Was intergenic reclassification exact or reservoir-estimated?
9. Are endpoint atlases tissue/protocol appropriate?
10. Is an orthogonal tool needed (ambient RNA, doublet, fusion, or isoform QC)?

Use component categories and evidence flags in reports. Do not describe `broad_noncanonical_*` as a measured contamination fraction.
