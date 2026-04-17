# scNoiseMeter Documentation

Version 0.2.0

---

## 1. Overview and Purpose

scNoiseMeter quantifies technical noise in single-cell RNA-seq BAM files. It classifies every primary alignment into one of 19 mutually exclusive read categories and reports per-sample and per-cell noise fractions, strand concordance, chimeric read rates, and artifact flag counts.

The tool is platform-agnostic: it processes ONT, PacBio/Kinnex, and Illumina (10x Genomics, BD Rhapsody) BAM files using the same classification logic, with platform-specific adjustments where the underlying biology differs.

Three subcommands are provided:

- `run` — classify reads in a single BAM and produce QC metrics
- `compare` — run on two BAMs (e.g. pre- and post-filter) and produce a side-by-side comparison
- `discover` — scan a directory for BAM files, infer their parameters, and run `scnoisemeter run` on selected files

The tool requires a coordinate-sorted, indexed BAM file aligned to human GRCh38/hg38, and a GENCODE GTF annotation. Both the GTF and a PolyASite 3.0 atlas are downloaded automatically on first use if not supplied explicitly.

---

## 2. Installation

```
pip install scnoisemeter
```

Dependencies (installed automatically):

- pysam >= 0.22
- pyranges >= 0.0.129
- pandas >= 2.0
- numpy >= 1.24
- click >= 8.1
- plotly >= 5.18
- scipy >= 1.11
- tqdm >= 4.66

Development extras:

```
pip install "scnoisemeter[dev]"
```

The BAM must be coordinate-sorted and indexed before use:

```
samtools sort -o sorted.bam input.bam
samtools index sorted.bam
```

---

## 3. Read Categories

Every primary alignment receives exactly one category. The classification hierarchy is applied in the order listed; a read is assigned the first matching category.

| Category | String value | Definition |
|---|---|---|
| `UNMAPPED` | `unmapped` | Read did not align. Excluded from all fractions. |
| `SECONDARY` | `secondary` | SAM flag 0x100. Record is skipped entirely. |
| `SUPPLEMENTARY` | `supplementary` | SAM flag 0x800. Routed to chimeric detector only; not counted in fractions. |
| `MULTIMAPPER` | `multimapper` | Primary alignment with NH tag > 1. |
| `CHIMERIC` | `chimeric` | SA tag present AND the split is inter-chromosomal, strand-discordant, or the same-strand intra-chromosomal distance exceeds the chimeric distance threshold (default 10,000 bp). For Illumina paired-end BAMs, also triggered when the absolute template length exceeds 1,000,000 bp. |
| `MITOCHONDRIAL` | `mitochondrial` | Maps to the mitochondrial contig (chrM, MT, chrMT, or mitochondrion). |
| `EXONIC_SENSE` | `exonic_sense` | Overlaps at least one annotated exon base on the correct strand. |
| `EXONIC_ANTISENSE` | `exonic_antisense` | Overlaps at least one annotated exon base on the wrong strand. |
| `INTRONIC_JXNSPAN` | `intronic_jxnspan` | Maps within an intron but the CIGAR contains an N operation near a splice site (candidate intron-retention or non-consensus transcript). |
| `INTRONIC_PURE` | `intronic_pure` | Maps entirely within an intron body with no junction signal. |
| `INTRONIC_BOUNDARY` | `intronic_boundary` | Spans an exon–intron boundary without a splice operation in the CIGAR (candidate incomplete reverse transcription). |
| `INTERGENIC_REPEAT` | `intergenic_repeat` | Intergenic read overlapping a RepeatMasker interval (requires `--repeats`). |
| `INTERGENIC_HOTSPOT` | `intergenic_hotspot` | Intergenic read at a locus that passes the adaptive barcode threshold but shows an internal-priming or A-rich 3′-end signature. |
| `INTERGENIC_NOVEL` | `intergenic_novel` | Intergenic read at a locus that passes the adaptive threshold, is strand-consistent, and falls near an annotated polyA site. Candidate unannotated gene or extended 3′ UTR. |
| `INTERGENIC_SPARSE` | `intergenic_sparse` | Intergenic read at a locus below the adaptive barcode threshold. Likely noise. |
| `AMBIGUOUS` | `ambiguous` | Overlaps a region shared by two or more genes where the gene types are not clearly distinguished by the sub-categories below. |
| `AMBIGUOUS_COD_NCOD` | `ambiguous_cod_ncod` | Overlaps a shared region between a protein-coding gene and a non-coding gene (lncRNA, pseudogene, etc.). |
| `AMBIGUOUS_COD_COD` | `ambiguous_cod_cod` | Overlaps a shared region between two protein-coding genes. |
| `UNASSIGNED` | `unassigned` | CB tag absent or not on the barcode whitelist. These reads are counted in the denominator but not attributed to a cell. |

Three categories — `UNMAPPED`, `SECONDARY`, `SUPPLEMENTARY` — are excluded from the `CATEGORY_ORDER` used for output columns and fraction computation. All other 16 categories appear in per-cell and per-sample output columns.

**Noise definitions**

Two noise levels are reported:

*Conservative noise* (`noise_read_frac`, `noise_base_frac`) — includes reads that may represent genuine pre-mRNA capture. This is an upper bound on true noise:

```
EXONIC_ANTISENSE + INTRONIC_PURE + INTRONIC_BOUNDARY +
INTERGENIC_SPARSE + INTERGENIC_REPEAT + INTERGENIC_HOTSPOT + CHIMERIC
```

*Strict noise* (`noise_read_frac_strict`, `noise_base_frac_strict`) — only unambiguous RT/PCR/sequencing artifacts. Excludes `INTRONIC_PURE` and `INTRONIC_BOUNDARY`. This is a lower bound:

```
EXONIC_ANTISENSE + INTERGENIC_SPARSE + INTERGENIC_REPEAT +
INTERGENIC_HOTSPOT + CHIMERIC
```

The categories `INTRONIC_JXNSPAN`, `INTERGENIC_NOVEL`, `AMBIGUOUS`, `AMBIGUOUS_COD_NCOD`, and `AMBIGUOUS_COD_COD` are in neither noise set; their interpretation is ambiguous.

**Adaptive intergenic threshold**

Intergenic loci are evaluated using a Poisson significance test against the expected read rate across all intergenic bases. Parameters:

- Minimum distinct barcodes: max(3, 0.01% of total detected barcodes)
- Minimum reads per locus: 5
- Bonferroni-corrected p-value threshold: 0.01
- Aggregation window: 500 bp

---

## 4. Artifact Flags

Three artifact flags are computed per read and counted at the sample and per-cell level. They are not part of the classification hierarchy — a read carries the flag in addition to its category.

### TSO invasion (`n_tso_invasion`)

Detects reads where soft-clipped bases at the 5′ end match a template-switching oligonucleotide (TSO) sequence. Detection requires at least 12 bp of match.

TSO sequences checked:

- 10x Genomics v3/v4: `AAGCAGTGGTATCAACGCAGAGTACATGGG`
- PacBio Kinnex / IsoSeq: `AAGCAGTGGTATCAACGCAGAGT`

A poly-G tail of ≥ 6 bp at the 5′ soft-clip is also flagged as TSO-proximal.

### Internal polyA priming (`n_polya_priming`)

Detects reads whose 3′ end is immediately followed by an A-run in the reference genome, indicating the read likely originated from internal priming on an A-rich region rather than the true polyA tail. Detection parameters:

- Context window inspected: 20 bp downstream of the read 3′ end
- Minimum A-run length: 6 consecutive A bases within the window

### Non-canonical junction (`n_noncanon_junction`)

Detects reads with a CIGAR N operation (intron) whose donor–acceptor dinucleotide pair is not one of the three canonical splice site motifs:

- GT–AG
- GC–AG
- AT–AC

Requires a reference FASTA (`--reference`) for the dinucleotide lookup. Without a reference, this flag is not computed.

---

## 5. Reported Metrics

### Sample-wide scalar metrics (in `<sample>.read_metrics.tsv`)

| Metric | Type | Definition |
|---|---|---|
| `n_reads_total` | int | Total reads in the BAM (mapped + unmapped, from index counters). |
| `n_reads_classified` | int | Reads that received a classification category. Excludes UNMAPPED, SECONDARY, SUPPLEMENTARY. |
| `n_reads_unassigned` | int | Reads classified as UNASSIGNED (CB absent or not on whitelist). |
| `n_cells` | int | Number of distinct cell barcodes with ≥ 10 reads. Set to 1 in barcode-agnostic mode. |
| `noise_read_frac` | float | Conservative noise fraction (reads). Denominator: n_reads_classified. |
| `noise_base_frac` | float | Conservative noise fraction (aligned bases). |
| `noise_read_frac_strict` | float | Strict noise fraction (reads). |
| `noise_base_frac_strict` | float | Strict noise fraction (aligned bases). |
| `strand_concordance` | float | exonic_sense / (exonic_sense + exonic_antisense). Values < 0.95 suggest strand-switching or a non-stranded library. |
| `chimeric_read_frac` | float | Fraction of classified reads in the CHIMERIC category. |
| `multimapper_read_frac` | float | Fraction of classified reads in the MULTIMAPPER category. |
| `per_cell_noise_median` | float | Median per-cell conservative noise fraction across cells with ≥ 10 reads. |
| `per_cell_noise_iqr` | float | Interquartile range of per-cell conservative noise fraction. |
| `n_tso_invasion` | int | Number of reads with a TSO invasion flag. |
| `n_polya_priming` | int | Number of reads with an internal polyA priming flag. |
| `n_noncanon_junction` | int | Number of reads with a non-canonical splice junction flag. |
| `read_frac_<category>` | float | Fraction of classified reads in each category (one row per category in CATEGORY_ORDER). |
| `base_frac_<category>` | float | Fraction of classified aligned bases in each category. |
| `full_length_read_frac` | float | Fraction of EXONIC_SENSE reads considered full-length. Present only when a polyA site database or `--polya-sites` is provided. Without a database, computed from a length threshold (ONT: 500 bp, PacBio: 1000 bp). |
| `tss_anchored_frac` | float | Fraction of reads with 5′ end within 100 bp of an annotated TSS. Present only when `--tss-sites` is supplied. |
| `numt_read_frac` | float | Fraction of MITOCHONDRIAL reads overlapping a NUMT interval. Present only when `--numt-bed` is supplied. |

### Per-cell metrics (columns of `<sample>.cell_metrics.tsv`)

| Column | Definition |
|---|---|
| `cell_barcode` | Cell barcode string (DataFrame index). |
| `n_reads` | Total classified reads for this cell. |
| `n_bases` | Total aligned bases for this cell. |
| `read_frac_<category>` | Per-cell read fraction for each category in CATEGORY_ORDER. |
| `base_frac_<category>` | Per-cell base fraction for each category. |
| `umi_complexity_<category>` | unique UMIs / total reads for this cell and category. Only present when UMI tracking is active (default; disabled by `--no-umi-dedup`). |
| `noise_read_frac` | Per-cell conservative noise read fraction. |
| `noise_base_frac` | Per-cell conservative noise base fraction. |
| `n_tso` | TSO invasion flag count for this cell. |
| `n_polya` | Internal polyA priming flag count for this cell. |
| `n_noncanon` | Non-canonical junction flag count for this cell. |

Cells with fewer than 10 reads are excluded from the per-cell table.

---

## 6. The `run` Subcommand

Classifies reads in a single BAM and produces all output files.

### Synopsis

```
scnoisemeter run [OPTIONS]
```

### Flags

**Required:**

| Flag | Type | Description |
|---|---|---|
| `--bam PATH` | path | Input BAM file. Must be coordinate-sorted and have a `.bai` index in the same directory. |
| `--output-dir PATH` | path | Directory for output files. Created if it does not exist. |

**Optional — input:**

| Flag | Default | Description |
|---|---|---|
| `--sample-name TEXT` | BAM filename stem | Label used in output filenames and HTML report. |
| `--gtf PATH` | auto-downloaded | GENCODE GTF annotation (plain or `.gz`). Takes precedence over `--gtf-version`. If omitted, the latest GENCODE human GTF is downloaded to `~/.cache/scnoisemeter/`. |
| `--gtf-version INT` | none | GENCODE release to auto-download (e.g. `42`). Ignored when `--gtf` is set. Use `42` to match the PolyASite 3.0 atlas exactly. |
| `--barcode-whitelist PATH` | none | File of valid corrected barcodes, one per line (plain text or `.gz`). Reads whose CB tag is not in this list are classified as UNASSIGNED. Distinct from `--cell-barcodes`. |
| `--cell-barcodes PATH` | none | Called-cell barcode file (plain text or `.gz`, one per line). Reads whose CB tag is not in this list are skipped entirely and contribute to no metric. Trailing `-1` suffixes are stripped from both the file entries and the CB tags in the BAM. Compatible with Cell Ranger `filtered_feature_bc_matrix/barcodes.tsv.gz`. |
| `--barcode-tag TEXT` | `CB` | BAM tag for the corrected cell barcode. |
| `--umi-tag TEXT` | `UB` | BAM tag for the corrected UMI. |
| `--chemistry [10x_v3\|10x_v4\|bd_rhapsody_wta\|custom]` | `10x_v3` | Library chemistry. Sets the expected barcode length. |
| `--platform [ont\|pacbio\|illumina\|illumina_10x\|illumina_bd\|unknown]` | `auto` | Sequencing platform. `auto` detects from BAM header `@PG` records. |
| `--pipeline-stage [raw\|pre_filter\|post_filter\|custom]` | `auto` | Processing stage. `auto` detects from header. |
| `--repeats PATH` | none | RepeatMasker BED file (hg38). Required to classify reads as INTERGENIC_REPEAT. |
| `--reference PATH` | none | Reference FASTA (`.fa` or `.fa.gz`, with `.fai` index). Required for polyA context checks and non-canonical junction detection. |
| `--polya-sites PATH` | auto-downloaded | PolyA site BED file(s). Repeatable. Takes precedence over `--polya-db`. If omitted, a database is auto-downloaded according to `--polya-db`. |
| `--polya-db [polyasite3\|polyadb4\|both]` | `polyasite3` | PolyA site database to auto-download when `--polya-sites` is not set. `polyasite3`: PolyASite 3.0 atlas (GENCODE v42, ~569k sites). `polyadb4`: PolyA_DB v4 (not tied to a GENCODE version, works with any GTF release). `both`: load both databases simultaneously. |
| `--tss-sites PATH` | auto-downloaded | CAGE peak / TSS BED file(s). Repeatable. Takes precedence over `--tss-db`. If omitted, a database is auto-downloaded according to `--tss-db`. |
| `--tss-db [fantom5\|none]` | `fantom5` | TSS database to auto-download when `--tss-sites` is not set. `fantom5`: FANTOM5 robust CAGE peaks (hg38, ~184k peaks). `none`: skip TSS anchoring entirely. |
| `--numt-bed PATH` | none | NUMT BED file (nuclear mitochondrial DNA segments, hg38 coordinates). |
| `--obs-metadata PATH` | none | Per-cell metadata TSV with `cell_barcode` and `cluster` columns. Enables per-cluster noise profiles. |
| `--exclude-biotypes TEXT` | none | Gene biotypes to exclude from annotation. Repeatable. |

**Optional — behaviour:**

| Flag | Default | Description |
|---|---|---|
| `--chimeric-distance INT` | `10000` | Maximum same-strand intra-chromosomal SA distance (bp) below which a split alignment is not called chimeric. |
| `--threads INT` | `4` | Parallel worker processes (one per chromosome). |
| `--no-umi-dedup` | off | Skip UMI set tracking. Reduces memory for very large datasets; disables `umi_complexity_*` columns. |
| `--no-cache` | off | Do not read or write the annotation index cache. |
| `--offline` | off | Use only cached annotation files; never make network calls. Raises an error if the cache is empty. Ignored when `--gtf` and `--polya-sites` are supplied explicitly. |
| `--verbose` / `-v` | off | Enable debug logging to stderr. |

### Validation performed before processing

- BAM index (`.bai`) must exist.
- BAM must be coordinate-sorted (`@HD SO:coordinate`).
- BAM must contain at least one aligned read.
- Chromosome naming style (UCSC `chr1` vs Ensembl `1`) must match between BAM and GTF. Mismatch is fatal.
- Chromosome lengths are compared against GRCh38 expected values. Mismatch is a warning.
- If `--cell-barcodes` is supplied but the BAM has no CB tags (barcode-agnostic mode), the filter is ignored with a warning.
- If `--cell-barcodes` filters out all reads, an error is raised.

### Examples

```bash
# Minimal run (GTF and polyA atlas auto-downloaded)
scnoisemeter run \
  --bam sample.bam \
  --output-dir results/

# Full run with explicit annotation
scnoisemeter run \
  --bam sample.bam \
  --gtf gencode.v45.annotation.gtf.gz \
  --barcode-whitelist 3M-february-2018.txt \
  --cell-barcodes filtered_feature_bc_matrix/barcodes.tsv.gz \
  --platform ont \
  --pipeline-stage post_filter \
  --repeats rmsk.hg38.bed.gz \
  --reference GRCh38.fa \
  --polya-sites atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz \
  --tss-sites hg38.cage_peak_phase1and2combined_ann.bed.gz \
  --threads 16 \
  --output-dir results/ \
  --sample-name my_sample

# Restrict analysis to called cells only
scnoisemeter run \
  --bam cellranger_output/possorted_genome_bam.bam \
  --cell-barcodes cellranger_output/filtered_feature_bc_matrix/barcodes.tsv.gz \
  --output-dir results/
```

---

## 7. The `compare` Subcommand

Runs the full classification pipeline on two BAMs, computes chi-squared proportion tests for each read category, and produces a comparison report.

### Synopsis

```
scnoisemeter compare [OPTIONS]
```

### Flags

**Required:**

| Flag | Type | Description |
|---|---|---|
| `--bam-a PATH` | path | BAM A (e.g. pre-filter / raw). |
| `--bam-b PATH` | path | BAM B (e.g. post-filter). |
| `--output-dir PATH` | path | Directory for output files. |

**Optional — compare-specific:**

| Flag | Default | Description |
|---|---|---|
| `--label-a TEXT` | `sample_A` | Label for BAM A in reports and output filenames. |
| `--label-b TEXT` | `sample_B` | Label for BAM B in reports and output filenames. |

**Shared flags (same as `run`, excluding `--bam`, `--sample-name`, `--cell-barcodes`):**

`--gtf`, `--gtf-version`, `--barcode-whitelist`, `--barcode-tag`, `--umi-tag`, `--chemistry`, `--platform`, `--pipeline-stage`, `--chimeric-distance`, `--repeats`, `--reference`, `--threads`, `--no-umi-dedup`, `--no-cache`, `--exclude-biotypes`, `--obs-metadata`, `--polya-sites`, `--polya-db`, `--tss-sites`, `--tss-db`, `--numt-bed`, `--offline`, `--verbose`

A single annotation index is built once and shared between both BAM runs.

### Statistical test

For each read category, a chi-squared test is applied to the contingency table of read counts (category reads vs all other classified reads) for BAM A and BAM B. P-values are Bonferroni-corrected for the number of categories tested.

### Examples

```bash
# Pre-filter vs post-filter comparison
scnoisemeter compare \
  --bam-a raw.bam \
  --bam-b filtered.bam \
  --gtf gencode.v45.annotation.gtf.gz \
  --label-a pre_filter \
  --label-b post_filter \
  --threads 8 \
  --output-dir compare_results/

# With explicit labels and whitelist
scnoisemeter compare \
  --bam-a sample1.bam \
  --bam-b sample2.bam \
  --label-a sample1 \
  --label-b sample2 \
  --barcode-whitelist 3M-february-2018.txt \
  --platform illumina_10x \
  --output-dir compare_results/
```

---

## 8. The `discover` Subcommand

Scans a directory for `.bam` files, inspects each to infer platform and pipeline stage, presents a summary table, and runs `scnoisemeter run` on selected files.

### Synopsis

```
scnoisemeter discover [OPTIONS]
```

### Flags

**Required:**

| Flag | Type | Description |
|---|---|---|
| `--bam-dir PATH` | directory | Directory to scan for `.bam` files. |
| `--reference PATH` | path | Reference FASTA (`.fa` or `.fa.gz` with `.fai` index). |
| `--output-dir PATH` | path | Root output directory. Each BAM gets its own subdirectory named after the BAM stem. |

**Optional:**

| Flag | Default | Description |
|---|---|---|
| `--gtf PATH` | auto-downloaded | GENCODE GTF. Auto-downloaded if omitted. Takes precedence over `--gtf-version`. |
| `--gtf-version INT` | none | GENCODE release to auto-download (e.g. `42`). Ignored when `--gtf` is set. |
| `--polya-sites PATH` | auto-downloaded | PolyA site BED file(s). Repeatable. Takes precedence over `--polya-db`. |
| `--polya-db [polyasite3\|polyadb4\|both]` | `polyasite3` | PolyA site database to auto-download when `--polya-sites` is not set. |
| `--tss-sites PATH` | auto-downloaded | TSS / CAGE peak BED file(s). Repeatable. Takes precedence over `--tss-db`. |
| `--tss-db [fantom5\|none]` | `fantom5` | TSS database to auto-download when `--tss-sites` is not set. |
| `--threads INT` | `4` | Parallel worker processes per BAM. |
| `--run-all` | off | Non-interactive mode. All BAMs with fully inferable parameters are run automatically. BAMs with blocking issues (no index, wrong sort order, unknown platform) are skipped with a warning. |
| `--offline` | off | Use only cached annotation files. |
| `--verbose` / `-v` | off | Enable debug logging. |

### Behaviour

1. All `.bam` files in `--bam-dir` are inspected to infer platform, pipeline stage, and check for blocking issues (missing index, wrong sort order).
2. A summary table is printed showing each BAM, its inferred platform, pipeline stage, and any issues.
3. In interactive mode (default), the user selects which BAMs to process. BAMs with unknown platform trigger an interactive prompt.
4. In `--run-all` mode, BAMs with blocking issues are skipped silently.
5. The annotation index is built once and shared across all selected BAMs.
6. Each BAM is run inline (no subprocess). Output goes to `<output-dir>/<bam-stem>/`.
7. A summary table of successes and failures is printed at the end.

### Examples

```bash
# Interactive discovery
scnoisemeter discover \
  --bam-dir /data/bams/ \
  --reference GRCh38.fa \
  --output-dir discover_results/

# Non-interactive batch run
scnoisemeter discover \
  --bam-dir /data/bams/ \
  --reference GRCh38.fa \
  --gtf gencode.v45.annotation.gtf.gz \
  --output-dir discover_results/ \
  --threads 8 \
  --run-all
```

---

## 9. Output Files

### `run` outputs

All files are written to `--output-dir`. `<sample>` is the value of `--sample-name`.

---

#### `<sample>.read_metrics.tsv`

Two-column tab-separated file: `metric` and `value`. One row per metric.

Rows written (in order):

| Row key | Description |
|---|---|
| `n_reads_total` | Total reads in BAM (from index counters). |
| `n_reads_classified` | Reads with a category (excludes UNMAPPED, SECONDARY, SUPPLEMENTARY). |
| `n_reads_unassigned` | Reads classified as UNASSIGNED. |
| `n_cells` | Distinct barcodes with ≥ 10 reads. |
| `noise_read_frac` | Conservative noise fraction (reads). |
| `noise_base_frac` | Conservative noise fraction (bases). |
| `strand_concordance` | exonic_sense / (exonic_sense + exonic_antisense). |
| `chimeric_read_frac` | CHIMERIC read fraction. |
| `multimapper_read_frac` | MULTIMAPPER read fraction. |
| `per_cell_noise_median` | Median per-cell noise fraction. |
| `per_cell_noise_iqr` | IQR of per-cell noise fraction. |
| `n_tso_invasion` | TSO invasion flag count. |
| `n_polya_priming` | Internal polyA priming flag count. |
| `n_noncanon_junction` | Non-canonical junction flag count. |
| `read_frac_<category>` | Read fraction per category (16 rows, CATEGORY_ORDER). |
| `base_frac_<category>` | Base fraction per category (16 rows). |
| `full_length_read_frac` | Full-length read fraction (present only when polyA sites are available). |

Values are plain numbers (integers or floats with 6 decimal places for fractions).

---

#### `<sample>.cell_metrics.tsv`

Tab-separated file. Index column is `cell_barcode`. One row per cell with ≥ 10 reads.

Columns:

| Column | Description |
|---|---|
| `cell_barcode` | Corrected cell barcode string (index). |
| `n_reads` | Total classified reads for this cell. |
| `n_bases` | Total aligned bases for this cell. |
| `read_frac_exonic_sense` | Fraction of this cell's reads in each category. One column per category in CATEGORY_ORDER. |
| `read_frac_<category>` | (16 columns total, one per category in CATEGORY_ORDER) |
| `base_frac_<category>` | Base fraction per category (16 columns). |
| `umi_complexity_<category>` | unique UMIs / total reads for this cell and category (16 columns; absent when `--no-umi-dedup` is set). |
| `noise_read_frac` | Per-cell conservative noise read fraction. |
| `noise_base_frac` | Per-cell conservative noise base fraction. |
| `n_tso` | TSO invasion flag count for this cell. |
| `n_polya` | Internal polyA priming flag count for this cell. |
| `n_noncanon` | Non-canonical junction flag count for this cell. |

This file is not written when the sample is in barcode-agnostic mode and the DataFrame is empty.

---

#### `<sample>.multiqc.json`

MultiQC custom-content JSON. Contains a subset of scalar metrics formatted for ingestion by MultiQC's `custom_content` module.

---

#### `<sample>.length_distributions/`

Directory containing one TSV per read category that has at least one read with a recorded length. Filenames: `<category_value>.lengths.tsv`.

Each file has a single column:

| Column | Description |
|---|---|
| `read_length` | Length of one read (one row per sampled read). |

Read lengths are collected by reservoir sampling; the sample size is bounded in memory. Not all reads are represented.

---

#### `<sample>_length_stratified.tsv`

Tab-separated. Rows represent the cross of length bin × read category.

| Column | Description |
|---|---|
| `length_bin` | Length bin label (e.g. `<150`, `150–500`, `500–1000`, `1000–2000`, `2000–5000`, `>5000`). |
| `category` | Read category string value. |
| `count` | Number of reads in this cell of the cross-tabulation. |
| `fraction_of_bin` | Fraction of reads in this length bin that belong to this category. |

For Illumina short-read data, a comment line is prepended noting that all reads fall in the `<150 bp` bin.

Bin breaks: 150, 500, 1000, 2000, 5000 bp. When median read length is ≥ 300 bp (long-read data), all six bins are reported. When median read length is < 300 bp (short-read data), the `<150` and `150–500` bins are merged into `<500`.

---

#### `<sample>.intergenic_loci.tsv`

Tab-separated. One row per intergenic locus characterised by the intergenic profiler. Only written if intergenic reads are present.

| Column | Description |
|---|---|
| `contig` | Chromosome / contig name. |
| `start` | Locus start coordinate (0-based). |
| `end` | Locus end coordinate. |
| `strand` | `+` or `-`. |
| `n_reads` | Number of intergenic reads at this locus. |
| `n_barcodes` | Number of distinct barcodes contributing reads. |
| `has_splice_evidence` | Boolean. True if any read at this locus has a junction (N in CIGAR). |
| `is_monoexonic` | Boolean. True if no read has a junction. |
| `polya_run_downstream` | Boolean. True if an A-run ≥ 6 bp was found downstream of the locus. |
| `near_polya_site` | Boolean. True if the locus is within 50 bp of an annotated polyA site. |
| `poisson_pvalue_adj` | Bonferroni-corrected Poisson p-value for read enrichment vs background intergenic rate. |
| `category` | Category assigned to this locus: `intergenic_hotspot`, `intergenic_novel`, or `intergenic_sparse`. |

---

#### `<sample>.cluster_metrics.tsv`

Only written when `--obs-metadata` is supplied. Contains per-cluster aggregated noise metrics. The cluster column comes from the `cluster` column in the obs metadata TSV.

---

#### `<sample>.report.html`

Self-contained interactive HTML report using Plotly (pinned to version 2.35.2 from CDN, or embedded when `--offline` is set). Contains:

- Sample metadata table (platform, pipeline stage, aligner, annotation versions)
- Read category composition bar chart
- Noise fraction summary
- Per-cell noise violin (suppressed for barcode-agnostic samples)
- Read-length distributions by category (long-read platforms only)
- Noise by read length stratification (long-read platforms only)
- Insert size distribution (Illumina only, when paired reads are present)
- Per-cluster noise comparison (when `--obs-metadata` is supplied)
- Intergenic loci table (when intergenic reads are present)
- Warnings panel

---

### `compare` outputs

All written to `--output-dir`.

---

#### `comparison.metrics.tsv`

Tab-separated. Side-by-side scalar metrics for both samples.

| Column | Description |
|---|---|
| `metric` | Metric name (same keys as `<sample>.read_metrics.tsv`). |
| `<label_a>` | Value for BAM A. |
| `<label_b>` | Value for BAM B. |
| `delta` | `<label_b>` − `<label_a>`. |

Only numeric metrics (int and float fields) from `SampleMetrics` are included.

---

#### `comparison.stats.tsv`

Tab-separated. One row per read category.

| Column | Description |
|---|---|
| `category` | Read category string value. |
| `frac_<label_a>` | Read fraction for BAM A. |
| `frac_<label_b>` | Read fraction for BAM B. |
| `delta` | `frac_<label_b>` − `frac_<label_a>`. |
| `chi2` | Chi-squared statistic from the contingency table test. |
| `p_value` | Uncorrected p-value. |
| `p_adjusted` | Bonferroni-corrected p-value (clipped at 1.0). |

---

#### `comparison.report.html`

Interactive HTML comparison report. Contains:

- Side-by-side read category composition bar charts
- Noise fraction comparison bars
- Per-cell noise violin (suppressed when either sample is barcode-agnostic, i.e. `n_cells == 1`)
- Length distribution overlays
- Statistical test results table
- Warnings panel

---

## 10. Annotation Caching and Auto-Download

### Cache location

All automatically downloaded annotation files are stored in `~/.cache/scnoisemeter/`. Subsequent runs reuse cached files without any network call. Pass `--offline` to enforce cache-only mode; the tool raises an error if a required file is absent from the cache.

### GTF

On first run with no `--gtf` or `--gtf-version`, the latest GENCODE human GTF is downloaded from the GENCODE FTP. The parsed annotation index (pyranges intervals, intron complement, intergenic regions) is cached alongside the GTF as a compressed pickle; rebuilding it from a large GTF takes roughly 60 seconds, so repeated runs benefit significantly from this cache. Pass `--no-cache` to force a rebuild.

To pin a specific GENCODE release:

```bash
scnoisemeter run --bam sample.bam --gtf-version 42 --output-dir results/
```

To supply a local file (disables all auto-download for the GTF):

```bash
scnoisemeter run --bam sample.bam --gtf gencode.v45.annotation.gtf.gz --output-dir results/
```

### PolyA site databases

Two databases are supported, selected with `--polya-db`:

| Database | Flag value | Source | Genome | Notes |
|---|---|---|---|---|
| PolyASite 3.0 | `polyasite3` (default) | polyasite.unibas.ch | hg38 / GENCODE v42 | ~569k sites; tied to GENCODE v42 |
| PolyA_DB v4 | `polyadb4` | exon.njms.rutgers.edu | hg38 | Not tied to a GENCODE version; works with any GTF release |

Use `both` to load both databases simultaneously. When `--polya-sites` is provided explicitly, `--polya-db` is ignored.

The PolyASite 3.0 atlas is distributed as a BED6 file. PolyA_DB v4 is distributed as a ZIP archive; scNoiseMeter downloads and converts it to BED3 format on first use and caches the result.

**Version mismatch.** The current PolyASite 3.0 atlas is built on GENCODE v42. Auto-downloading the latest GTF (currently v49) produces a seven-version gap; the tool warns when the difference exceeds five major releases. Two ways to resolve this:

- Pass `--gtf-version 42` to auto-download GENCODE v42, matching the PolyASite 3.0 atlas exactly.
- Pass `--polya-db polyadb4` to switch to PolyA_DB v4, which is not tied to a GENCODE version and works with any GTF release.

### TSS / CAGE peak databases

Two options are supported, selected with `--tss-db`:

| Database | Flag value | Source | Notes |
|---|---|---|---|
| FANTOM5 | `fantom5` (default) | fantom.gsc.riken.jp | hg38 robust CAGE peaks, ~184k peaks, BED6 format |
| None | `none` | — | Disables TSS anchoring; `tss_anchored_frac` will not be reported |

When `--tss-sites` is provided explicitly, `--tss-db` is ignored.

---

## 11. Platform-Specific Notes

### ONT

- Platform auto-detected from `minimap2` `@PG` record in the BAM header.
- Chimeric detection uses the SA tag with the default intra-chromosomal distance threshold (10,000 bp). Split alignments within this distance on the same strand are treated as legitimate splices.
- Full-length read fraction fallback threshold: 500 bp (when no polyA site database is provided).
- Read-length distribution and noise-by-length charts are included in the HTML report.
- No insert size chart.

### PacBio / Kinnex

- Platform auto-detected from `pbmm2` `@PG` record.
- Chimeric detection uses the SA tag, same logic as ONT.
- The PacBio TSO (`AAGCAGTGGTATCAACGCAGAGT`) is used for TSO invasion detection.
- Full-length read fraction fallback threshold: 1000 bp.
- Read-length distribution and noise-by-length charts are included in the HTML report.
- No insert size chart.

### Illumina (10x Genomics, BD Rhapsody)

- Platform auto-detected from `STAR`, `STARsolo`, or `cellranger` `@PG` records.
- Chimeric detection uses paired-end mode: a read pair is chimeric if it is inter-chromosomal, strand-discordant, or has `abs(template_length)` ≥ 1,000,000 bp.
- Read-length distribution and noise-by-length charts are suppressed (all reads are the same short length).
- Insert size distribution chart is shown when properly paired reads are present (collected by reservoir sampling from read1 of each proper pair, with `0 < abs(template_length) < 2000`).
- The `_length_stratified.tsv` file contains a note that all reads fall in the `<150 bp` bin.
- `illumina_10x` and `illumina_bd` are treated identically to `illumina` in all classification logic; the distinction affects only BAM header auto-detection.

### Barcode-agnostic mode

Activated when fewer than 50% of sampled reads (10,000 reads sampled from the BAM) carry the corrected barcode tag (`CB` by default). In this mode:

- All reads are aggregated under the sentinel barcode `NO_BARCODE`.
- `n_cells` is set to 1.
- Per-cell metrics are not meaningful; the cell_metrics TSV will contain one row.
- A warning is emitted to stderr.
- If `--cell-barcodes` is supplied in barcode-agnostic mode, the filter is ignored with a warning.
- In `compare` reports, the per-cell violin is suppressed when either sample has `n_cells == 1`, and a warning is appended to the report.

---

## 12. Known Caveats and Limitations

**Genome / annotation:**

- Only human GRCh38/hg38 is supported. The chromosome length validation uses hardcoded GRCh38 expected lengths. Other species will produce length mismatch warnings.
- The GENCODE GTF and PolyASite 3.0 atlas must use the same chromosome naming convention (UCSC or Ensembl). Mismatches between BAM and GTF chromosome names cause a fatal error.
- If the GTF and polyA atlas differ by more than 5 GENCODE major versions, a warning is issued. Genes with 3′ UTRs annotated between the two versions may have reduced polyA anchoring scores. To resolve: pass `--gtf-version 42` to match the PolyASite 3.0 atlas exactly, or pass `--polya-db polyadb4` to switch to PolyA_DB v4, which is version-agnostic.

**Read classification:**

- Only primary alignments are classified. Secondary (flag 0x100) and supplementary (flag 0x800) alignments are skipped, except that supplementary records are examined by the chimeric detector.
- `MULTIMAPPER` is defined as NH tag > 1 on the primary alignment record. Reads without an NH tag are not flagged as multimappers.
- `INTRONIC_PURE` and `INTRONIC_BOUNDARY` cannot be distinguished from genuine pre-mRNA capture at the read level. They are included in conservative noise but excluded from strict noise. Their presence does not necessarily indicate an artifact.
- `INTERGENIC_NOVEL` means a locus passes the adaptive barcode threshold and is near an annotated polyA site. It does not confirm the existence of an unannotated gene.
- `INTERGENIC_REPEAT` classification requires a RepeatMasker BED file (`--repeats`). Without it, repeat-overlapping intergenic reads fall into INTERGENIC_HOTSPOT or INTERGENIC_SPARSE.
- Non-canonical junction detection requires a reference FASTA (`--reference`). Without it, `n_noncanon_junction` is 0 regardless of the actual data.

**Barcode handling:**

- The tool reads corrected barcodes from the `CB` tag (or the tag specified by `--barcode-tag`). Raw uncorrected barcodes (`CR` tag) are not used for classification.
- Trailing `-1` suffixes are stripped from both the `--cell-barcodes` file and the CB tags in the BAM, so Cell Ranger output is normalised automatically. Other suffix conventions are not handled.
- In barcode-agnostic mode, per-cell noise values are not meaningful.

**Chimeric detection:**

- The default chimeric distance threshold of 10,000 bp may flag very long transcripts (> 10 kb) that have legitimate split alignments as chimeric. For datasets with very long transcripts (e.g. PacBio full-length mRNA), consider increasing `--chimeric-distance`.
- For Illumina paired-end data, the 1,000,000 bp insert size threshold for chimeric calling is fixed and cannot be adjusted from the command line.

**Performance:**

- The annotation index is built per run. For large GTF files, this takes several minutes on first use. The cache (`--no-cache` disables it) stores the parsed index to avoid rebuilding on subsequent runs with the same GTF.
- Parallelism is at the chromosome level (`--threads`). Chromosomes are processed independently. Small contigs (alt, patch) are each dispatched as separate workers, which may cause imbalance for heavily fragmented reference assemblies.
- UMI tracking (`umi_complexity_*` columns) stores a set of UMI strings per cell per category. For very large datasets, this can consume significant memory. Use `--no-umi-dedup` to disable it.
- Read-length and insert-size sampling use reservoir sampling (Algorithm R). Not all reads are represented in the length distribution TSVs or insert size charts.

**Statistical:**

- The chi-squared test in `compare` is applied to the contingency table of read counts. It is not a paired test; it does not account for the fact that BAM B may be a strict subset of BAM A (e.g. post-filter ⊆ pre-filter). Interpret p-values accordingly.
- The adaptive intergenic threshold uses Bonferroni correction across all intergenic loci. In samples with many sparse intergenic reads, this correction is conservative and may suppress detection of low-coverage novel loci.
- Per-cell noise statistics (median, IQR) are computed only over cells with ≥ 10 reads. The threshold is fixed.
