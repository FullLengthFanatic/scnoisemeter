# scNoiseMeter — validation benchmark (v1.1)

Simulation-based validation of the scnoisemeter classifier and per-cell noise
metric. Runs in a couple of minutes. Outputs a confusion matrix and a
noise-fraction dose-response table suitable as a methods-paper supplementary
figure.

## Coverage

**Exp 1 — per-read classifier accuracy (14 of 19 ReadCategory values).**

Simple categories (one CB, 200 reads each):

- EXONIC_SENSE, EXONIC_ANTISENSE
- INTRONIC_PURE, INTRONIC_BOUNDARY, INTRONIC_JXNSPAN
- INTERGENIC_SPARSE
- MITOCHONDRIAL, MULTIMAPPER, CHIMERIC
- AMBIGUOUS_COD_COD, AMBIGUOUS_COD_NCOD

Clustered intergenic categories (3 CBs, 10 loci, 21 reads per locus) — the
intergenic profiler requires >= 3 distinct barcodes per locus to call
significance, so each category ships three collaborating CBs:

- INTERGENIC_HOTSPOT (mono-exonic clusters with an A-run downstream)
- INTERGENIC_NOVEL (spliced clusters with a non-canonical donor dinucleotide)
- INTERGENIC_REPEAT (clusters whose positions are covered by `repeats.bed`)

**Exp 2 — per-cell noise-fraction fidelity.** Five synthetic cells with
known signal (EXONIC_SENSE) + noise mixtures. True noise fraction sweeps
0 → 55 %. Compared against the reported `noise_read_frac`.

## What it does NOT validate

- **Alignment-quality artifacts.** The BAM is synthesized directly with pysam;
  no aligner is invoked. Sequencing-error profiles, soft-clips, and mapQ
  distributions are stylized rather than realistic.
- **UNASSIGNED**, polyA-priming / TSO-invasion flags, and NUMT-based
  mitochondrial sub-classification — each needs additional inputs and is
  deferred.
- **Short-read coverage (Illumina, ElemBio) and 384-well `run-plate`** —
  scope for a future pass; the classifier is platform-agnostic, so ONT-mode
  correctness implies short-read correctness for read-level classification.

## Requirements

- A GENCODE GTF and matching genome FASTA (GRCh38). Either download manually
  or reuse the files already cached by scnoisemeter under
  `~/.cache/scnoisemeter/`.
- `pysam`, `pandas` — already in scnoisemeter's install dependencies.

## Run

```bash
# 1. Synthesize BAM + ground-truth table + repeats.bed
python3 tests/benchmark/synthesize_bam.py \
    --gtf   ~/.cache/scnoisemeter/gencode.v49.annotation.gtf.gz \
    --fasta /path/to/GRCh38.primary_assembly.genome.fa.gz \
    --outdir /tmp/scnm_bench

# 2. Run scnoisemeter + score against ground truth
python3 tests/benchmark/run_and_evaluate.py \
    --bam     /tmp/scnm_bench/synthetic.bam \
    --truth   /tmp/scnm_bench/ground_truth.tsv \
    --gtf     ~/.cache/scnoisemeter/gencode.v49.annotation.gtf.gz \
    --fasta   /path/to/GRCh38.primary_assembly.genome.fa.gz \
    --repeats /tmp/scnm_bench/repeats.bed \
    --outdir  /tmp/scnm_bench/out
```

Pass `--skip-run` to `run_and_evaluate.py` to reuse an existing scnoisemeter
output directory (faster iteration on the evaluator). Omit `--repeats` to
skip REPEAT classification and measure only HOTSPOT / NOVEL.

## Outputs

| File                                    | What it contains                                   |
|-----------------------------------------|----------------------------------------------------|
| `synthetic.cell_metrics.tsv`            | scnoisemeter per-cell metrics (standard output)    |
| `benchmark_confusion_matrix.tsv`        | rows = true category, cols = assigned, values = reads |
| `benchmark_per_category_metrics.tsv`    | per-category total / correct / accuracy            |
| `benchmark_noise_dose_response.tsv`     | per-cell true vs reported noise fraction           |
| `benchmark_summary.json`                | aggregate metrics (mean/min accuracy, delta)       |

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

## Not yet covered

- **UNASSIGNED** (requires unmapped or no-CB reads and a whitelist/QC path).
- **TSO-invasion / polyA-priming flags** on the EXONIC_SENSE axis.
- **NUMT-derived MITOCHONDRIAL** sub-classification (needs a NUMT BED).
- **pbsim3 + minimap2** for realistic ONT/PacBio error profiles and alignment
  artifacts.
- **Short-read (Illumina, ElemBio) and 384-well plate mode** end-to-end.
