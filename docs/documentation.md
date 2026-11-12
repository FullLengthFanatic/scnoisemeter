# scNoiseMeter documentation

Version 0.9.1 operational reference.

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

`--barcode-whitelist` defines valid barcode sequences. A primary alignment with a missing/off-list barcode becomes `unassigned` and remains in sample metrics. Reads with no barcode at all are left out of the per-cell table; an off-list barcode with at least 10 reads appears there as a row that is 100% `unassigned`.

`--cell-barcodes` defines called cells. Reads not in this set are dropped and do not affect metrics. Cell Ranger’s `filtered_feature_bc_matrix/barcodes.tsv.gz` is accepted; a trailing `-1` is normalized for matching.

Without a whitelist or CB tag, reads pool under `NO_BARCODE`; sample-level metrics still work, but per-cell metrics do not. Plate mode converts that sentinel to a well identifier before aggregation.

### GTF

`--gtf` accepts plain or gzipped GTF. If omitted, a GENCODE annotation is cached automatically. `--gtf-version N` selects a release. Chromosome naming must match the BAM (`chr1` versus `1`). Annotation cache identity includes metadata and a content-head hash.

### Reference FASTA

`--reference` enables:

- strand-aware internal polyA-context flags;
- donor-plus-acceptor splice-motif checks;
- intergenic hotspot context.

The FASTA must match the BAM assembly and have an `.fai` index. Without `--reference` these checks do not run, and their outputs are left empty instead of 0. `n_polya_priming`, `n_noncanon_junction` and `n_rt_switch_repeat` are blank in `read_metrics.tsv` and `per_well_metrics.tsv`. `polya_priming_frac`, `rt_switch_repeat_frac` and the `intergenic_hotspot` fractions are left out of `read_metrics.tsv` and `multiqc.json`. The per-cell columns `n_polya`, `n_noncanon` and `*_frac_intergenic_hotspot` are empty. Reads that a hotspot call would have claimed stay in the other intergenic categories, so the composition still sums to one. A FASTA that cannot be opened is treated the same way.

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
| `--polya-min-support` | 25 | drop polyA atlas records supported in fewer than this percentage of samples (column 7 of a PolyASite BED). 0 loads the atlas whole; see the note below |
| `--priming-criterion` | `a_run` | internal-priming test: `a_run` (consecutive genomic A's), `perc_a` (SQANTI3/pigeon `perc_A_downstream_TTS`, ≥60% A over 20 bp), or `either` |
| `--anchor-veto-internal-priming` | off | require the absence of an internal-priming signature for a read to count towards 3′ anchoring; lowers `three_prime_anchored_frac` |
| `--no-strict` | off | downgrade the pre-write sanity checks from errors to warnings (see §9a) |

**On `--polya-min-support`.** Two intergenic rules depend on proximity to an annotated polyA site, in opposite directions: `intergenic_hotspot` requires that no site is nearby, `intergenic_novel` accepts a nearby site as evidence. PolyASite 3.0 holds 18,432,135 records, including sites supported in 5.6% of samples, and loaded whole puts 37.4% of the genome within the 50 bp proximity window, at which point the test is nearly a constant. Filtering at 25 leaves 6.3% of the genome; 50 leaves 2.4%. Every run measures the density of the atlas it loaded and warns above 20%. The default is 25 since 0.9.1. Versions up to 0.9.0 loaded the atlas whole; pass 0 to reproduce their results.

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

The command validates coordinate sort order, BAM index counts, chromosome naming, annotation/reference compatibility, and low read counts before writing the report. Since 0.9 it also stops on two conditions that previously produced a complete but misleading result set — see §9a.

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

Without `--run-all`, discovery displays inferred metadata and prompts for a subset. The annotation and site resources are reused across BAMs. Inference warnings should be resolved with explicit settings when protocol details are known. A BAM whose result set fails the pre-write checks (§9a) is listed as failed and its outputs are not written; `--no-strict` writes them anyway.

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

## 9a. Checks that stop a run

Three checks run before any output is written. Each can be downgraded to a warning with `--no-strict`, which writes the outputs and records the violation instead. `run` applies all three, and so does `compare`, to each of its two BAMs. `discover` and `run-plate` apply the barcode gate and the result invariants, not the site-atlas coordinate check. In `discover` a BAM that fails is reported as failed and the remaining BAMs still run.

**Barcode gate.** If more than half of classified reads carry no usable cell barcode on a barcode-aware BAM, the run stops and prints the observed counts, the resulting cell count and the header-scan tag fraction. Previously the run finished and reported a sample with a handful of cells and an unremarkable composition — a configuration error arriving disguised as a finding. Barcode-free BAMs are exempt, where 100% unassigned is the design.

**Site-atlas coordinate range.** Contig names matching does not mean the assembly matches. An hg19 atlas against an hg38 BAM shares every contig name, so the existing naming check passes while every coordinate is displaced and every proximity test silently returns nothing. Positions falling beyond the end of their contig are the detectable part: above 20% out of range the atlas is discarded and its metrics reported as absent; above 1% it warns.

**Result invariants.** Every run checks the finished result set against the structural properties the tool claims for it, and writes `<sample>.qc.json` with each verdict and its observed values. Checked: each composition sums to one; every fraction lies in [0, 1]; `artifact_candidate_*` does not exceed `broad_noncanonical_*`, whose categories contain it; `both_ends_anchored_frac` does not exceed either single-end fraction, being their conjunction on the same reads; no flag count exceeds `n_reads_classified`; UNASSIGNED does not exceed the classified total; UMI diversity ratios stay in [0, 1]; `n_cells` matches the per-cell table; `strand_concordance` matches its definition recomputed from the reported composition.

A violated invariant means a number about to be written is wrong, so the run stops. Two checks are warnings rather than errors because a documented configuration breaks them: reservoir-estimated reclassification scales promoted totals from a sub-sample, so exact conservation need not survive rounding; and `--anchor-veto-internal-priming` lowers the unrestricted 3′ fraction without lowering the PAS-restricted subset. Checks whose inputs are absent are skipped, not failed.

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

Each flag is written as a count of classified alignments. The matching rates (count / `n_reads_classified`) are in `<sample>.multiqc.json`: `tso_invasion_frac`, `tso_concatemer_frac`, `polya_priming_frac`, `polya_tail_frac`, `rt_switch_repeat_frac` and `discordant_pair_frac`.

- `n_tso_invasion`;
- `n_tso_concatemer`;
- `n_polya_priming` — genomic A-run beyond the read's 3′ end. A property of the locus, not of the molecule;
- `n_polya_tail` — untemplated A-run in the read's own 3′ **soft** clip. A property of the molecule. Read with the previous flag: both set means a real 3′ end in A-rich sequence, which the genomic flag alone calls internal priming. Hard-clipped BAMs report 0 because the bases are absent from the record;
- `n_noncanon_junction`;
- `n_rt_switch_repeat` — an intron whose two ends share a ≥6 bp direct repeat, the reverse-transcriptase template-switch signature. Independent of junction canonicality, because such an artifact can occupy a GT-AG boundary. Requires a reference FASTA;
- `n_discordant_pair`;
- `n_low_mapq`.

### Ends and provenance

- optional `three_prime_anchored_frac`;
- optional `three_prime_at_pas_site_frac` — the same 3′ end restricted to atlas records naming a polyA signal hexamer (column 11 of a PolyASite BED); absent when the atlas has no such column;
- optional `five_prime_anchored_frac`;
- optional `both_ends_anchored_frac`;
- deprecated `full_length_read_frac` equal to both-end anchoring;
- optional `n_numt_intervals_loaded`.

## 11. `cell_metrics.tsv`

Only non-sentinel cells with at least 10 classified reads are included. Columns include:

- `n_reads`, `n_bases`;
- category read/base fractions;
- `broad_noncanonical_*` and `artifact_candidate_*` fractions;
- flag counts (`n_tso`, `n_polya`, `n_noncanon`, `n_tso_concat`, `n_discordant_pair`, `n_low_mapq`) and `low_mapq_read_frac`;
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
| `local_read_share` | fraction of intergenic read mass within 10 kb held by this window; near 1 is an isolated locus, low is diffuse coverage that the global test still promotes |
| `n_molecules` | distinct (barcode, UMI) pairs; `n_reads` counts PCR copies, so one over-amplified molecule can look like broad support. Exact count, no mismatch collapsing, so an upper bound. `-1` when no UMIs were tracked, which is not zero molecules |
| `strand_asymmetry_pvalue` | two-sided binomial test of the strand split against 50/50. Makes strand evidence depth-aware: 3 of 4 and 75 of 100 are both `strand_fraction` 0.75 and very different evidence. `NaN` for unstranded libraries |
| `poisson_pvalue` | raw global-rate tail probability |
| `poisson_pvalue_adj` | Bonferroni result across the intergenic windows actually tested. Before 0.9 the correction counted every window of every contig including genic sequence, a 2.2–2.9× overcorrection on GRCh38; intergenic counts are not comparable across that boundary |
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

### `<sample>.qc.json`

Every structural invariant checked over this result set, with its verdict, its statement and the observed values. Written by `run` and `discover` per sample, by `compare` per label (`<label>.qc.json`) and by `run-plate` per plate, including when all checks pass: a file listing only failures cannot be distinguished from one whose checks never ran.

### `<sample>.run_info.json`

Written by `run` and `run-plate` alongside the metrics: scnoisemeter version, sample name, BAM
path, platform, pipeline stage, aligner, strandedness, read and cell totals, and the GTF, polyA,
TSS and TSO sources with how each was obtained. Since 0.9 it also carries a `priming` block with
`polya_min_support`, the internal-priming `criterion`, and whether the anchoring veto was on --
three settings that change what a reported number means, so `cohort` can flag runs that
disagree on them, and with `intergenic_reclassification_sampled`, which is true when the
intergenic categories were estimated from the reservoir. `read_metrics.tsv` holds only numbers, so this is what lets a later reader, or
`cohort`, establish whether two result sets are comparable.

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
