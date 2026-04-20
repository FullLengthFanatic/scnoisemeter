# scNoiseMeter — v1 validation benchmark

Simulation-based validation of the scnoisemeter classifier and per-cell noise
metric. Runs in under a minute. Outputs a confusion matrix and a noise-fraction
dose-response table suitable as a methods-paper supplementary figure.

## What this benchmark validates

**Exp 1 — per-read classifier accuracy.** For each of 8 read categories we
place `N=200` reads of a single true category into a dedicated synthetic cell
barcode. After running scnoisemeter, the per-cell category fractions directly
reveal how the classifier assigned those reads. Reading the results back as a
confusion matrix gives per-category accuracy.

**Exp 2 — per-cell noise-fraction fidelity.** Five synthetic cells are built
with known mixtures of signal (EXONIC_SENSE) and noise drawn from
NOISE_CATEGORIES_CONSERVATIVE. True noise fraction sweeps 0 → 55 %. We compare
the reported `noise_read_frac` to the true fraction.

## What it does NOT validate (v1 scope boundaries)

- **Alignment-quality artifacts.** The BAM is synthesized directly with
  `pysam`; no aligner is invoked. Sequencing-error profiles, soft-clips, and
  mapQ distributions are stylized, not realistic.
- **Intergenic profiler categories** (HOTSPOT, NOVEL, REPEAT). These require
  locus-level read density patterns and a RepeatMasker annotation, planned for
  v2.
- **INTRONIC_JXNSPAN.** Requires CIGAR N constructs landing on known splice
  sites; v2.
- **AMBIGUOUS_COD_COD / COD_NCOD.** Require overlapping-gene regions; v2.

## Requirements

- A GENCODE GTF and matching genome FASTA (GRCh38). Either download manually
  or reuse the files already cached by scnoisemeter under
  `~/.cache/scnoisemeter/`.
- `pysam`, `pandas` — already in scnoisemeter's install dependencies.

## Run

```bash
# 1. Synthesize BAM + ground-truth table
python3 tests/benchmark/synthesize_bam.py \
    --gtf   ~/.cache/scnoisemeter/gencode.v49.annotation.gtf.gz \
    --fasta /path/to/GRCh38.primary_assembly.genome.fa.gz \
    --outdir /tmp/scnm_bench

# 2. Run scnoisemeter + score against ground truth
python3 tests/benchmark/run_and_evaluate.py \
    --bam    /tmp/scnm_bench/synthetic.bam \
    --truth  /tmp/scnm_bench/ground_truth.tsv \
    --gtf    ~/.cache/scnoisemeter/gencode.v49.annotation.gtf.gz \
    --fasta  /path/to/GRCh38.primary_assembly.genome.fa.gz \
    --outdir /tmp/scnm_bench/out
```

Pass `--skip-run` to `run_and_evaluate.py` to reuse an existing scnoisemeter
output directory (faster iteration on the evaluator).

## Outputs

After evaluation, `--outdir` contains:

| File                                    | What it contains                                   |
|-----------------------------------------|----------------------------------------------------|
| `synthetic.cell_metrics.tsv`            | scnoisemeter per-cell metrics (standard output)    |
| `benchmark_confusion_matrix.tsv`        | rows = true category, cols = assigned, values = reads |
| `benchmark_per_category_metrics.tsv`    | per-category total / correct / accuracy            |
| `benchmark_noise_dose_response.tsv`     | per-cell true vs reported noise fraction           |
| `benchmark_summary.json`                | aggregate metrics (mean/min accuracy, delta)       |

## Interpretation

- **Per-category accuracy** should be ≥ 0.95 for every category covered in v1.
  Lower values indicate a classifier regression or an ambiguous synthetic
  construction. Ties in base counts resolve by Python dict insertion order —
  the synthesizer uses a 30/70 exon/intron split for INTRONIC_BOUNDARY reads
  so argmax is unambiguous.
- **`delta_vs_true`** in the dose-response should be within ± 0.02 at every
  mixture level. Larger deltas at higher noise fractions usually trace back
  to a single category misclassification compounding across the mixture.

## v2 plan (out of scope here)

- pbsim3 + minimap2 for realistic ONT/PacBio error profiles and alignment
  artifacts.
- INTERGENIC_HOTSPOT / NOVEL validation: spike N loci with varying read
  counts, measure sensitivity/specificity of the Poisson profiler.
- INTERGENIC_REPEAT validation: requires a RepeatMasker BED input.
- Short-read coverage (Illumina, ElemBio) via synthesized 10x-style BAMs.
- 384-well plate variant exercising `run-plate` and `--parallel-wells`.
