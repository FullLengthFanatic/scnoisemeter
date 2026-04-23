# scNoiseMeter — validation benchmark (v1.2)

Simulation-based validation of the scnoisemeter classifier, per-cell noise
metric, artifact-flag counters, and Illumina paired-end chimeric detection.
Runs in a couple of minutes. Outputs confusion matrices and a noise-fraction
dose-response table suitable as a methods-paper supplementary figure.

## Coverage

**Exp 1 — per-read classifier accuracy (16 of 19 ReadCategory values).**

Simple categories (one CB, 200 reads each):

- EXONIC_SENSE, EXONIC_ANTISENSE
- INTRONIC_PURE, INTRONIC_BOUNDARY, INTRONIC_JXNSPAN
- INTERGENIC_SPARSE
- MITOCHONDRIAL, MULTIMAPPER, CHIMERIC
- AMBIGUOUS, AMBIGUOUS_COD_COD, AMBIGUOUS_COD_NCOD
- UNASSIGNED (CB present but not in the supplied whitelist)

Clustered intergenic categories (3 CBs, 10 loci, 21 reads per locus) — the
intergenic profiler requires >= 3 distinct barcodes per locus to call
significance, so each category ships three collaborating CBs:

- INTERGENIC_HOTSPOT (mono-exonic clusters with an A-run downstream)
- INTERGENIC_NOVEL (spliced clusters with a non-canonical donor dinucleotide)
- INTERGENIC_REPEAT (clusters whose positions are covered by `repeats.bed`)

**Exp 2 — per-cell noise-fraction fidelity.** Five synthetic cells with
known signal (EXONIC_SENSE) + noise mixtures. True noise fraction sweeps
0 → 55 %. Compared against the reported `noise_read_frac`.

**Exp 4 — artifact-flag counters.** Two dedicated cells emit reads designed
to trip the TSO-invasion and polyA-priming flags respectively (200 reads
each). Validated against per-cell `n_tso` / `n_polya` counters in
`cell_metrics.tsv`.

**Exp 5 — Illumina paired-end chimeric detection.** A small separate BAM
(`synthetic_illumina.bam`, 400 reads, 4 CBs) exercises the paired-end
chimeric branch (`_check_chimeric_paired`) via three mechanisms: mate on
a different contig, FR/RR discordance, and oversized TLEN
(> ILLUMINA_CHIMERIC_INSERT_SIZE). Scored with `--platform illumina_10x`.

## What it does NOT validate

- **Alignment-quality artifacts.** BAMs are synthesized directly with pysam;
  no aligner is invoked. Sequencing-error profiles, soft-clips, and mapQ
  distributions are stylized rather than realistic.
- **NUMT-based mitochondrial sub-classification** — requires a NUMT BED
  (feature not yet implemented).
- **384-well `run-plate`** — separate workflow scope.
- **pbsim3 + minimap2 realism** — deferred to a future pass.

## Requirements

- A GENCODE GTF and matching genome FASTA (GRCh38). Either download manually
  or reuse the files already cached by scnoisemeter under
  `~/.cache/scnoisemeter/`.
- `pysam`, `pandas` — already in scnoisemeter's install dependencies.

## Run

```bash
# 1. Synthesize ONT BAM + Illumina BAM + ground-truth tables + repeats.bed + whitelist.txt
python3 tests/benchmark/synthesize_bam.py \
    --gtf   ~/.cache/scnoisemeter/gencode.v49.annotation.gtf.gz \
    --fasta /path/to/GRCh38.primary_assembly.genome.fa.gz \
    --outdir /tmp/scnm_bench

# 2. Run scnoisemeter (ONT + Illumina passes) + score against ground truth
python3 tests/benchmark/run_and_evaluate.py \
    --bam              /tmp/scnm_bench/synthetic.bam \
    --truth            /tmp/scnm_bench/ground_truth.tsv \
    --gtf              ~/.cache/scnoisemeter/gencode.v49.annotation.gtf.gz \
    --fasta            /path/to/GRCh38.primary_assembly.genome.fa.gz \
    --repeats          /tmp/scnm_bench/repeats.bed \
    --barcode-whitelist /tmp/scnm_bench/whitelist.txt \
    --illumina-bam     /tmp/scnm_bench/synthetic_illumina.bam \
    --illumina-truth   /tmp/scnm_bench/ground_truth_illumina.tsv \
    --outdir           /tmp/scnm_bench/out
```

Pass `--skip-run` to `run_and_evaluate.py` to reuse an existing scnoisemeter
output directory (faster iteration on the evaluator). Omit `--repeats` to
skip REPEAT classification and measure only HOTSPOT / NOVEL. Omit
`--illumina-bam` / `--illumina-truth` to skip Exp 5. Omit
`--barcode-whitelist` to disable UNASSIGNED validation (Exp 1 will drop that
category).

## Outputs

| File                                          | What it contains                                         |
|-----------------------------------------------|----------------------------------------------------------|
| `synthetic.cell_metrics.tsv`                  | ONT per-cell metrics (standard scnoisemeter output)      |
| `synthetic_illumina.cell_metrics.tsv`         | Illumina per-cell metrics                                |
| `benchmark_confusion_matrix.tsv`              | Exp 1 — rows = true category, cols = assigned            |
| `benchmark_per_category_metrics.tsv`          | Exp 1 — per-category total / correct / accuracy          |
| `benchmark_noise_dose_response.tsv`           | Exp 2 — per-cell true vs reported noise fraction         |
| `benchmark_artifact_flags.tsv`                | Exp 4 — per-cell expected vs observed flag counts        |
| `benchmark_illumina_confusion_matrix.tsv`     | Exp 5 — Illumina per-read confusion matrix               |
| `benchmark_illumina_per_category_metrics.tsv` | Exp 5 — Illumina per-category total / correct / accuracy |
| `benchmark_summary.json`                      | aggregate metrics across all experiments                 |

## Interpretation

- **Per-category accuracy** should be >= 0.95 for simple categories and >= 0.90
  for clustered intergenic categories (HOTSPOT/NOVEL/REPEAT depend on Poisson
  significance thresholds that can reclassify borderline loci as SPARSE).
  Lower values indicate a classifier regression or an ambiguous synthetic
  construction.
- Ties in base counts resolve by Python dict insertion order — the synthesizer
  uses a 30/70 exon/intron split for INTRONIC_BOUNDARY reads so argmax is
  unambiguous.
- **`delta_vs_true`** in the dose-response should be within ± 0.02 at every
  mixture level. Larger deltas at higher noise fractions usually trace back
  to a single category misclassification compounding across the mixture.
- **Artifact-flag counts** (Exp 4) should match the expected count exactly
  for the target flag and stay near zero on the other flag.
- **Illumina chimeric** (Exp 5) should be >= 0.95 on each of the three
  chimeric variants (inter-chromosomal, discordant orientation, oversized
  TLEN) and on CONCORDANT → EXONIC_SENSE.

## Not yet covered

- **NUMT-derived MITOCHONDRIAL** sub-classification (needs a NUMT BED and
  the feature itself is not yet implemented).
- **pbsim3 + minimap2** for realistic ONT/PacBio error profiles and alignment
  artifacts.
- **384-well plate mode** end-to-end.
