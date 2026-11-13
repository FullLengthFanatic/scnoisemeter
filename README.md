# scNoiseMeter

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19554841.svg)](https://doi.org/10.5281/zenodo.19554841)

Current version: **0.10.0**

Quality control for single-cell RNA-seq alignments. It answers two questions about a library: where do the reads land, and which of them look like technical artifacts?

You give it an aligned BAM file. scNoiseMeter puts every read into one of 17 categories (exonic, intronic, intergenic, mitochondrial, chimeric and so on), adds independent flags for known artifact signatures such as template-switch-oligo (TSO) sequence or internal polyA priming, and writes an interactive HTML report plus plain-text tables, per sample and per cell. It only describes the library. No read is removed and no count matrix is corrected.

It works on long reads (ONT, PacBio/Kinnex) and short reads (10x, BD Rhapsody, Smart-seq/FLASH-seq plates), with or without cell barcodes. Typical uses are checking a new library preparation, comparing protocols, and seeing what a filtering step removed.

Two summary numbers lead every report:

- `broad_noncanonical_*`: antisense, unspliced intronic, most intergenic, and chimeric reads taken together. Much of this can be real biology, for example pre-mRNA in nuclei.
- `artifact_candidate_*`: the smaller part with positive evidence of an artifact (`intergenic_hotspot` and `chimeric`).

Neither is a contamination estimate. Report the individual categories next to them.

Where to read more:

- [Full documentation](docs/documentation.md): every command, option and output column.
- [White paper](docs/whitepaper.md): what the categories mean and what can be concluded from them.
- [Annotated methods](docs/methods_annotated.md): which function implements which rule.
- [Changelog](CHANGELOG.md): every release that changed a reported number.

## Installation

```bash
pip install git+https://github.com/FullLengthFanatic/scnoisemeter.git
```

For development:

```bash
git clone https://github.com/FullLengthFanatic/scnoisemeter.git
cd scnoisemeter
pip install -e ".[dev]"
pytest -q
```

Python 3.10 or newer is required. Input BAMs must be coordinate-sorted and indexed:

```bash
samtools sort -o sorted.bam input.bam
samtools index sorted.bam
```

## Quick start

```bash
scnoisemeter run \
  --bam sample.bam \
  --platform ont \
  --library-strand stranded \
  --output-dir results/
```

`--platform` is one of `ont`, `pacbio`, `illumina`, `illumina_10x`, `illumina_bd` or `smartseq`. The first run downloads a GENCODE annotation, the PolyASite 3.0 polyA-site atlas and the FANTOM5 CAGE peaks (about 330 MB in total) into `~/.cache/scnoisemeter/`; later runs reuse them. These automatic resources are human GRCh38 only. For any other genome, pass your own files (see below).

Add `--reference genome.fa` (with its `.fai` index) whenever you can. Without it, the checks that need the genome sequence do not run and their values are left blank instead of 0: internal polyA priming, non-canonical junctions, RT template-switch repeats, and the `intergenic_hotspot` category.

When the run finishes, open `results/sample.report.html`. The tables behind it are listed under [Outputs](#outputs).

Add `--tagged-bam results/sample.classified.bam` to make a coordinate-sorted, indexed copy with the final category in local SAM tag `sn`. The option retains read assignments in memory and performs a second BAM pass, so it is off by default.

For reproducible work, pass every resource explicitly:

```bash
scnoisemeter run \
  --bam sample.bam \
  --gtf gencode.v45.annotation.gtf.gz \
  --reference GRCh38.fa \
  --repeats hg38.repeatmasker.bed.gz \
  --cell-barcodes filtered_feature_bc_matrix/barcodes.tsv.gz \
  --platform ont \
  --library-strand stranded \
  --seed 42 \
  --threads 8 \
  --output-dir results/
```

Do not rely on an aligner name alone to identify the sequencer or library. `minimap2` is not proof of ONT, and bare `STAR` is not proof of Smart-seq. Ambiguous headers produce a warning; use `--platform` and `--library-strand` explicitly.

## Report preview

![scNoiseMeter classification overview and category definitions](docs/report_preview2.png)

The interactive HTML report also covers [sample metadata](docs/report_preview1.png), [read and base fractions](docs/report_preview3.png), [length-stratified and per-cell composition](docs/report_preview4.png), [read-length distributions and artifact flags](docs/report_preview5.png), and [the top intergenic loci](docs/report_preview6.png).

These screenshots come from a 10x Kinnex sample (37.4 M reads, 2,047 cells) run with version 0.10.0.

## Classification

Unmapped, secondary, and supplementary records are counted separately and excluded from the category denominator. Every remaining primary alignment gets exactly one category, in two steps:

1. Four gates are checked in order: `unassigned`, `multimapper`, `mitochondrial`, `chimeric`. The first that applies is the category.
2. Otherwise the aligned bases are compared with the annotation, and the read takes the category that holds most of its aligned bases. Ties go to the category listed first in the table.

A 1 kb unspliced read with 800 bases in an exon and 200 in the next intron is therefore `exonic_sense` as a read, while the base-level fractions record 800 exonic and 200 `intronic_boundary` bases. Read fractions and base fractions are both reported for this reason.

Counts named `n_reads_*` are alignment records, so the two mates of a pair count separately. `n_records_*` and `n_alignments_classified` are the same numbers under more precise names.

| Category | Interpretation |
|---|---|
| `unassigned` | Barcode missing/off-whitelist when a whitelist is enforced |
| `multimapper` | Explicit `NH>1` evidence. MAPQ and producer-specific `X*` tags do not replace the genomic category; low MAPQ is reported separately. **`minimap2` and `pbmm2` do not emit `NH`**, so this category is empty by construction on most long-read BAMs and those alignments are counted in their genomic category instead — see [Scope and limitations](#scope-and-limitations) |
| `mitochondrial` | Primary alignment is on `chrM`, `MT`, `chrMT`, or `mitochondrion` |
| `chimeric` | Inter-contig or strand-discordant `SA`, incompatible query/genomic split order, or extreme paired insert size |
| `ambiguous_cod_cod` | Same-strand region shared by protein-coding genes |
| `ambiguous_cod_ncod` | Same-strand region shared by coding and non-coding genes |
| `exonic_sense` | Exonic overlap on the read strand |
| `exonic_antisense` | Exonic overlap on the opposite strand, with no gene on the read's strand there |
| `intronic_jxnspan` | Intronic bases on an alignment containing a CIGAR `N` |
| `intronic_boundary` | Exonic and intronic bases without a splice operation |
| `intronic_pure` | Intronic bases without boundary or junction evidence |
| `intergenic_repeat` | At least 50% of sampled aligned span overlaps supplied repeat intervals |
| `intergenic_hotspot` | Enriched fixed window with a strand-aware genomic polyA-run signature, away from an annotated polyA site |
| `intergenic_novel` | Enriched, strand-consistent, multi-barcode window with splice or annotated polyA-end support |
| `intergenic_enriched` | Enriched window lacking evidence for a more specific interpretation |
| `intergenic_sparse` | Intergenic window below the enrichment/support thresholds |
| `ambiguous` | Defensive fallback for an otherwise unresolved alignment |

For the base-level fractions, aligned blocks are split at every annotation boundary and each piece receives exactly one category, so base counts sum exactly to the aligned reference bases even where genes overlap. Same-strand gene sharing is ambiguous; overlapping genes on opposite strands remain strand-resolvable. Where a piece lies in an intron of a gene on the read's strand and in an exon of a gene on the other strand, the read's own strand decides and the piece is intronic.

Genomic distance by itself is not used to call an RNA split chimeric: a long gap may simply be an intron. Orphan/improper pairs are reported as `discordant_pair` flags and are not automatically called chimeric.

## Intergenic enrichment

Intergenic reads are assigned by their strand-correct 3′ coordinate to predefined 500 bp windows. Using BAM contig lengths, the denominator is the exact complement of the union of gene bodies through each non-mitochondrial contig end. Each window is tested against a Poisson model with one genome-wide rate, and the Bonferroni correction counts every window of that intergenic space, including the empty ones. This is a heuristic enrichment screen, not a peak caller: it does not model local background, GC content, or mappability.

Strand consistency is measured after window assignment. A window can therefore fail strand evidence rather than being split into artificial strand-specific loci. Significant unresolved windows become `intergenic_enriched`, not `intergenic_hotspot`.

Each window also reports `local_read_share`: the fraction of intergenic read mass within 10 kb that the window itself holds. A globally significant count says a window exceeds the genome-wide rate; it does not say whether the window is the local event or one of many in a broadly covered neighbourhood. A share near 1 is an isolated locus; a low share is diffuse coverage that the global test still promotes.

**PolyA-site density matters here.** Two rules depend on proximity to an annotated polyA site in opposite directions: `intergenic_hotspot` requires the *absence* of a nearby site, and `intergenic_novel` accepts a nearby site as one of its two positive-evidence routes. PolyASite 3.0 ships 18,432,135 records, including sites supported in 5.6% of samples at a normalised score of 6e-05. Loaded unfiltered, 37.4% of the genome falls within the 50 bp proximity window, so the test is closer to a constant than a discriminator: hotspots are suppressed across a third of the genome and the novel-gene route is close to free. `--polya-min-support` filters on the atlas's own per-site support column: 25 leaves 6.3% of the genome within the window and 50 leaves 2.4%. Every run reports the measured density of the atlas it loaded. The default is 25. Pass 0 to load the atlas whole, which is what versions up to 0.9.0 did.

For datasets larger than the 500,000-record intergenic reservoir, promoted read/base totals are estimated from a uniform reservoir. UMI-diversity fields for affected categories are emitted as missing because exact category-specific UMI sets cannot be reconstructed from a sample.

The Poisson test itself also runs on that reservoir, not only the totals derived from it: both the background rate and the per-window counts come from the sampled records. A sample well above the cap is therefore tested at reduced depth and a sample below it at full depth, so promotion to `intergenic_hotspot` or `intergenic_novel` is depth-dependent. `<sample>.run_info.json` records whether reclassification was reservoir-estimated; treat these two categories as not directly rankable across samples of very different intergenic depth. Intergenic category counts from versions before 0.9 used a larger correction and are not comparable (see the [changelog](CHANGELOG.md)).

## Independent evidence flags

Flags do not change the mutually exclusive read category:

- TSO invasion: the start of the template-switch oligo found in the soft-clipped 5′ end of the read; one substitution is tolerated.
- TSO concatemer: multiple approximate full-motif occurrences.
- Internal polyA priming: downstream A-run on `+` or upstream T-run on `-`; requires a reference FASTA. This describes the *locus* in the reference, not the molecule: measured against PolyASite, median downstream A-content is 30% both at annotated polyA sites and at random genomic positions, so the genomic test alone barely separates them. `--priming-criterion perc_a` switches to the SQANTI3/pigeon `perc_A_downstream_TTS` convention (≥60% A over 20 bp); the two agree on roughly 95 to 97% of windows, so the choice is about comparability with isoform-QC numbers rather than correctness.
- Untemplated 3′ polyA tail: an A-run in the read's own 3′ **soft** clip, outside the aligned block, so it is not templated by the reference — per-molecule evidence that this read carried a tail. A read flagged for both internal priming and a tail is a real cleavage site in an A-rich context. Reported as `polya_tail_frac`. Hard-clipped BAMs report no tail, because the bases are absent from the record rather than absent from the molecule.
- Non-canonical junction: exact donor and acceptor checked in transcript orientation; exact annotated junctions are accepted.
- RT template-switch direct repeat: the two ends of one intron share a direct repeat of at least 6 bp, which reverse transcriptase can duplicate by switching template between the copies. Independent of canonicality — an RT artifact can sit on a GT-AG boundary and is invisible to the motif check. Genuine junctions in repetitive sequence also carry direct repeats, so this is evidence, never a category. Reported as `rt_switch_repeat_frac`; requires a reference FASTA.
- Discordant pair: orphan or improper paired alignment.
- Low MAPQ: descriptive `MAPQ < 10` flag; no read is filtered by this threshold.

CLI TSO defaults are protocol-aware. ONT/10x uses the 10x motif, PacBio/Smart-seq uses the SMART/PacBio motif, and generic Illumina, BD, or unknown inputs use no motif unless `--tso` is supplied. The poly-G shortcut is enabled only with the built-in 10x motif.

## End anchoring and alignment quality

“Full-length” is not inferred from read length. For the same reservoir-sampled exonic-sense read, scNoiseMeter reports:

- `three_prime_anchored_frac`: 3′ end near a strand-matched polyA site. Proximity alone counts by default; `--anchor-veto-internal-priming` additionally requires the absence of an internal-priming signature, which lowers the fraction.
- `three_prime_at_pas_site_frac`: the same 3′ end near a site the atlas annotates with a polyA signal hexamer (column 11 of a PolyASite BED) — a stricter anchor than proximity to any record;
- `five_prime_anchored_frac`: 5′ end near a strand-matched TSS/CAGE site;
- `both_ends_anchored_frac`: both conditions on that same read;
- deprecated `full_length_read_frac`: alias of `both_ends_anchored_frac`, only when both atlases are available.

These fractions are suppressed for explicitly unstranded libraries because read orientation does not identify transcript ends. Alignment-quality outputs include mean MAPQ, mean `NM` where present, soft-clipped/query-base fraction, unmapped fraction from BAM index counts, and primary/secondary/supplementary/QC-fail/duplicate totals.

`umi_sequence_diversity_<category>` is unique UMI strings divided by reads. It is not molecule complexity because UMIs are not grouped by gene/locus or error-corrected. The older `umi_complexity_*` columns remain aliases.

## Commands

### One BAM

```bash
scnoisemeter run --bam sample.bam --output-dir results/
```

### Pre/post-filter comparison

```bash
scnoisemeter compare \
  --bam-a raw.bam \
  --bam-b filtered.bam \
  --label-a pre \
  --label-b post \
  --output-dir comparison/
```

The comparison does not use an invalid independent-samples chi-square test. It produces exact read-key retention and transition tables plus descriptive composition deltas and a paired-cell bootstrap interval for the median per-cell change.

Retention and transitions require the two BAMs to be the *same reads* before and after a step. `compare` samples read names from both BAMs and, when they do not match, skips those tables and says so rather than emitting NaN. This also avoids holding a per-read map of both BAMs in memory, which for two long-read samples runs to tens of gigabytes. Note that some deduplicators rewrite read names (`isoseq dedup` replaces PacBio CCS names with `molecule/N`), so a genuine pre/post pair can still be unmatchable; the composition metrics remain valid. Force either behaviour with `--matched-reads` / `--no-matched-reads`.

### Many samples

```bash
scnoisemeter cohort \
  --results results/BD46/ \
  --results results/10x_FL/ \
  --results results/PIPseq/ \
  --sample-sheet cohort.tsv \
  --output-dir cohort/
```

Compares any number of independent samples by reading the metrics `run` already wrote, so it takes seconds rather than the hours a re-classification would need. Use it whenever the samples are separate experiments; use `compare` only for a nested pre/post pair of the same reads.

The optional sample sheet has columns `sample`, `label`, `group`, `order`, where `sample` matches the `<sample>.read_metrics.tsv` stem. Without it, samples are labelled by that stem and ordered cleanest first.

The report carries four figures: composition per sample with exonic sense excluded so the differences are not compressed into the tail of the bar; a heatmap coloured by deviation from the cohort median, which shows what is *unusual* about a method rather than restating the composition; artifact-flag rates on a log axis, since those span several orders of magnitude; and per-cell spread for whichever samples have barcodes.

Samples that came from different scnoisemeter versions, different GENCODE releases, or a mix of stranded and unstranded protocols are flagged in the report. Metrics that a sample never reported are shown as absent, never as zero.

### Directory discovery

```bash
scnoisemeter discover \
  --bam-dir /data/bams \
  --reference GRCh38.fa \
  --run-all \
  --output-dir batch_results/
```

### Plate data

```bash
scnoisemeter run-plate \
  --plate-dir /data/plate_881 \
  --sample-sheet plate_881.csv \
  --platform smartseq \
  --library-strand unstranded \
  --parallel-wells 8 \
  --output-dir results/
```

Barcode-free well BAMs are relabeled to `<plate>_<well>` before aggregation, preventing all wells from collapsing into a single `NO_BARCODE` pseudo-cell.

## Outputs

| File | Contents |
|---|---|
| `<sample>.read_metrics.tsv` | Denominators, category fractions, aggregate compositions, flags, endpoint and alignment-quality metrics |
| `<sample>.cell_metrics.tsv` | Per-cell category, aggregate, flag, UMI-diversity, and alignment-quality metrics (cells with at least 10 reads) |
| `<sample>.intergenic_loci.tsv` | Coordinates, strand/repeat evidence, raw and adjusted Poisson p-values, final category |
| `<sample>_length_stratified.tsv` | Exact category counts by length bin |
| `<sample>.length_distributions/` | Reservoir-sampled read lengths by category |
| `<sample>.cluster_metrics.tsv` | Optional summaries from `--obs-metadata` |
| `<sample>.multiqc.json` | MultiQC-compatible scalar content, including each artifact flag as a rate (`*_frac`) |
| `<sample>.report.html` | Interactive report |
| user path from `--tagged-bam` | Optional full BAM copy; classified primaries carry `sn:Z:<category>` |
| `comparison.retention.tsv` | Category-specific exact read retention |
| `comparison.transitions.tsv` | Category A→B transitions for matched read keys |
| `comparison.matching.tsv` | Overall matching counts |
| `comparison.stats.tsv` | Composition deltas and paired-cell bootstrap intervals |
| `comparison.metrics.tsv` | Scalar metrics of both BAMs side by side |
| `comparison.report.html` | Interactive comparison report |
| `<sample>.run_info.json` | Tool version, platform, strandedness, annotation sources and the settings that change what a number means |
| `<sample>.qc.json` | Every structural invariant checked over this result set, its verdict and the observed values |
| `cohort.summary.tsv` | One row per sample: provenance and headline metrics |
| `cohort.composition.tsv` | Samples x categories, read and base fractions |
| `cohort.report.html` | Cross-sample report |

`comparison.retention.tsv`, `comparison.transitions.tsv` and `comparison.matching.tsv` are written only for a nested pair.

## Result invariants

The test suite proves the code behaves on constructed inputs. Separately, `run`, `compare`, `discover` and `run-plate` each check whether *this* result set satisfies the properties the tool claims about its own output, and writes them to `<sample>.qc.json` with the observed values rather than a bare verdict.

The claims being checked are structural rather than statistical, which is what makes them checkable: each composition sums to one, every fraction lies in [0, 1], the artifact-candidate composition cannot exceed the broad non-canonical one because its categories are a subset, `both_ends_anchored_frac` cannot exceed either single-end fraction because it is their conjunction on the same reads, no artifact-flag count exceeds the denominator it is reported over, unique UMI strings cannot outnumber the reads carrying them, and `n_cells` matches the per-cell table.

A violation means a number about to be written is wrong, so the run stops; `--no-strict` writes the outputs anyway and records the violation. Two checks are warnings rather than errors because a documented configuration can break them: reservoir-estimated intergenic reclassification scales promoted totals from a uniform sub-sample, so exact conservation is not expected to survive rounding, and `--anchor-veto-internal-priming` vetoes the unrestricted 3'-anchoring fraction but not the PAS-restricted subset.

## Scope and limitations

- Classification is diagnostic; it does not remove reads or correct count matrices.
- Intronic, antisense, intergenic, mitochondrial, and chimeric observations can all reflect biology. Aggregate fields are compositions, not contamination estimates.
- `intergenic_novel` is a candidate label, not gene discovery. Validate with independent end, splice, expression, and replication evidence.
- A global Poisson background is approximate, especially near highly transcribed or poorly mappable regions.
- SA evidence can represent technical chimeras or genuine fusions.
- Multimapping is only detected where the aligner writes `NH`. STAR, STARsolo and Cell Ranger do; `minimap2` and `pbmm2` do not. Comparing `multimapper_read_frac`, or any category whose denominator it shifts, between a long-read and a short-read sample is therefore not meaningful.
- Duplicate- and QC-fail-flagged records are counted and reported, but are still classified. Composition from a BAM that has been through `MarkDuplicates` is duplicate-weighted; a Cell Ranger BAM is not.
- NUMT BED input records annotation provenance only. A NUMT read fraction is not claimed without competing-alignment evidence.
- Coordinate and annotation compatibility remain the user’s responsibility; the current automatic reference resources target human GRCh38.

## Citation and license

Please cite the archived concept DOI:

> Picelli, S. scNoiseMeter: barcode-aware alignment artifact and read-distribution QC for single-cell RNA-seq. Zenodo. https://doi.org/10.5281/zenodo.19554841

MIT licensed; see [LICENSE](LICENSE).
