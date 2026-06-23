# scNoiseMeter: a white paper

**Platform-agnostic, per-read classification of technical noise in single-cell RNA-seq.**

Version 0.6.0. Written for readers who work with RNA-seq data and want to understand what the tool measures and why it is built the way it is. For the line-by-line implementation, see [methods_annotated.md](methods_annotated.md).

---

## 1. The problem

Single-cell RNA-seq libraries carry a substantial fraction of reads that do not represent a real, full-length, correctly-stranded mRNA molecule from the cell they are assigned to. The usual suspects: intronic reads from pre-mRNA or incomplete reverse transcription, antisense reads from strand-switching, reads landing in intergenic space, chimeras from template switching or ligation, internal priming on genomic A-runs, template-switching-oligo (TSO) artifacts, and ambient or multi-gene-ambiguous signal.

Most QC tools collapse all of this into one number, a single "% mapped to transcriptome" or "% intronic". That number tells you something is wrong but not what, and it gives you no handle on whether a given cell, a given read-length bin, or a given locus is the source. It also tends to be tied to one platform's aligner and one assay's conventions.

scNoiseMeter takes a different stance. It reads a coordinate-sorted BAM plus a GTF and assigns **every primary alignment to exactly one of 16 mutually-exclusive categories**, then layers independent sequence-level artifact flags on top. The same logic runs on long-read (ONT, PacBio/Kinnex) and short-read (Illumina, Element/AVITI; 10x or BD Rhapsody kits) and on plate-based Smart-seq / FLASH-seq, with the few genuinely platform-dependent rules switched automatically. The output is a mechanistic breakdown of noise, per sample, per cell, and per cluster, not a single scalar.

---

## 2. The core idea

Three design choices define the tool.

**Per-read, mutually-exclusive classification.** Each primary alignment gets one category by a fixed priority hierarchy (first matching rule wins). In parallel, the tool tallies *aligned bases* per category, so a long read that straddles an exon and an intron contributes to both at base resolution even though its single read-level label is decided by plurality vote. Read-level labels drive per-cell QC; base-level fractions are the finer benchmarking metric.

**Two noise definitions, not one.** Some categories are unambiguous artifacts (intergenic sparse, chimeric, antisense in a stranded library). Others (intronic reads with no junction) cannot be told apart from genuine pre-mRNA capture at the read level. So the tool reports a **conservative** noise fraction (upper bound, includes the intronic-pure subtypes) and a **strict** noise fraction (lower bound, only unambiguous artifacts). For unstranded protocols it drops antisense from both, because sense and antisense arrive in roughly equal proportion by design.

**Statistics where heuristics are not enough.** Intergenic reads are not all equal. A handful scattered across the genome is background; a tight cluster supported by many cells is either a real unannotated transcript or an internal-priming hotspot. scNoiseMeter resolves this with a second pass that clusters intergenic reads into loci and scores each against a Poisson background, rather than calling everything intergenic "noise".

---

## 3. The classification logic

Reads are tested in this order. The order matters: it is a priority hierarchy, and the first rule that fires assigns the category.

1. **Skip filters.** Unmapped, secondary (`0x100`), and supplementary (`0x800`) records are not classified. The supplementary record's `SA` tag is still read off the primary alignment for chimera detection.
2. **Barcode gate.** With a barcode whitelist, a read whose `CB` tag is missing or off-list becomes `unassigned`. With no whitelist and no `CB`, the tool runs in barcode-agnostic mode and pools everything under one sentinel barcode, so it still produces sample-level noise metrics on BAMs that carry no cell barcodes.
3. **Multimapper.** `NH > 1` on the primary alignment. Aligner-agnostic, so it works the same on minimap2, STAR, pbmm2.
4. **Mitochondrial.** Maps to `chrM` / `MT` / `chrMT` / `mitochondrion`.
5. **Chimeric.** From the `SA` tag: inter-chromosomal, strand-discordant, or same-strand at a distance above the threshold (default 10 kb). For paired-end short-read and Smart-seq, a fallback uses mate mapping and insert size when no `SA` tag is present.
6. **Genomic location** (the remaining reads), decided by overlapping the read's aligned blocks against the annotation:
   - `exonic_sense` / `exonic_antisense`: overlaps an annotated exon on the correct / wrong strand.
   - `intronic_jxnspan`: intronic with a `CIGAR N` (a real splice gap), candidate intron retention or a non-canonical junction.
   - `intronic_boundary`: spans an exon-intron edge with no splice gap, the signature of incomplete RT.
   - `intronic_pure`: entirely inside an intron body.
   - `intergenic_sparse`: no gene-body overlap (the default before the second pass promotes it).
   - `ambiguous_cod_cod` / `ambiguous_cod_ncod` / `ambiguous`: the read sits in a region shared by two protein-coding genes, by a coding and a non-coding gene, or by genes the tool cannot otherwise resolve.

The signal category is `exonic_sense`. Everything else is either noise, ambiguous, or (for the intronic and intergenic-novel cases) deliberately left out of the noise sets because it may be biology.

### Why the intronic split

A 3' assay should produce exonic reads ending at the polyA site. Intronic reads happen for real reasons (nuclear pre-mRNA, intron retention) and for technical reasons (RT that stalled, internal priming inside an intron). The tool cannot fully separate these per read, so instead of guessing it records the *shape* of the evidence: a junction-spanning intronic read (`intronic_jxnspan`) is treated as plausibly real and kept out of noise; a read that crosses an exon-intron boundary with no splice (`intronic_boundary`) or sits purely in an intron (`intronic_pure`) goes into conservative noise but not strict noise. The reader gets both bounds and decides.

### Why the intergenic second pass

After the first pass labels intergenic reads `intergenic_sparse`, a second pass clusters them into loci (single-linkage within 500 bp on the same strand) and tests each locus against a Poisson model whose rate is the genome-wide intergenic read density. Loci that clear a Bonferroni-corrected p-value plus minimum read and barcode support are promoted:

- **`intergenic_hotspot`**: monoexonic and sitting on a genomic A-run just past its 3' end. This is the internal-priming signature, oligo-dT or the RT priming off an A-rich stretch in the genome instead of a real tail.
- **`intergenic_novel`**: strand-consistent, multi-cell, and anchored by either a splice junction or proximity to an annotated polyA site. A candidate unannotated transcript, reported as ambiguous rather than noise.
- **`intergenic_repeat`**: overlaps a RepeatMasker interval (only when a repeats BED is supplied).

A locus that is significant but fits neither the hotspot nor the novel rule defaults to hotspot, on the principle that over-flagging an artifact is safer than over-claiming a new gene. Because promoted reads leave the sparse bucket before noise is computed, calling a locus `intergenic_novel` actually lowers the reported noise fraction.

### Sequence-level artifact flags

Separate from the category, every classified read can carry independent flags. These are not categories because a read can be, say, exonic-sense and still show a TSO artifact; the flags annotate, they do not relabel.

- **TSO invasion.** The 10x and PacBio TSOs (or a user-supplied sequence) matched against the read's soft-clips, in both forward and reverse-complement orientation (so a TSO end that mapped antisense is still caught), plus a poly-G heuristic for the 10x TSO's G-rich tail.
- **TSO concatemer.** More than one TSO occurrence in a single read, the signature of template switching running from one TSO onto another. The metric follows Chou et al. (bioRxiv 2025.10.06.680646).
- **Internal polyA priming.** A genomic A-run within 20 bp downstream of a forward read's 3' end (or a T-run upstream of a reverse read's 3' end). Strand-aware, and it needs a reference FASTA.
- **Non-canonical junction.** A `CIGAR N` whose donor dinucleotide is not GT-AG / GC-AG / AT-AC and is not in the annotation. Also needs a reference FASTA.

---

## 4. Platform-agnostic by design

The same classifier runs everywhere. What changes between platforms is auto-detected from the BAM header (`@PG` aligner plus pipeline hints) and then drives a small set of switches:

- **Strandedness.** Stranded kits (10x, BD, ONT, PacBio) count antisense as noise. Smart-seq / FLASH-seq are unstranded, so antisense is excluded from the noise sets. This one is *not* auto-detected: Smart-seq must be set with `--platform smartseq`, because a STAR-aligned plate BAM otherwise looks like 10x.
- **Chimera logic.** Long reads use the `SA`-tag distance rule; paired-end short-read and Smart-seq add the mate/insert-size fallback.
- **Read-length views.** Length-stratified charts are shown for long reads and suppressed for short-read kits where every read is sub-150 bp; insert-size distributions are shown for paired-end data instead.
- **Full-length estimate.** With a polyA-site atlas, "full length" means the 3' end sits at an annotated polyA site. Without one, it falls back to a platform-specific length threshold.

Annotation (GENCODE GTF, PolyASite 3.0 or PolyA_DB v4, FANTOM5 CAGE) is auto-downloaded and cached on first use, so the minimal invocation is just a BAM and an output directory.

---

## 5. What you get out

| Output | Content |
|---|---|
| `read_metrics.tsv` | Sample-wide fractions, conservative and strict noise, strand concordance, chimeric/multimapper rates, artifact-flag counts |
| `cell_metrics.tsv` | Per-cell breakdown for every cell with >= 10 reads, including per-category UMI complexity |
| `intergenic_loci.tsv` | Each characterized locus with its Poisson-adjusted p-value and the features that drove its call |
| `_length_stratified.tsv` | Category composition by read-length bin |
| `multiqc.json` | MultiQC-ready scalars |
| `report.html` | Interactive Plotly report |

Four subcommands cover the workflows: `run` (one BAM), `run-plate` (a 96- or 384-well Smart-seq / FLASH-seq plate), `compare` (two BAMs, e.g. pre- vs post-filter, with per-category chi-squared tests), and `discover` (scan a directory, infer parameters, batch-run). A `--obs-metadata` cluster file adds per-cluster noise profiles.

---

## 6. Where it sits relative to existing tools

scNoiseMeter is not an ambient-RNA remover, an empty-droplet caller, a doublet detector, or an isoform-QC tool. It is a per-alignment, mechanistic noise classifier, and it is most useful as a diagnostic layer alongside those tools rather than a replacement for any of them.

- **Aligner built-in metrics (Cell Ranger, STARsolo)** report four broad read classes (exonic / intronic / intergenic / antisense) as sample-level fractions. scNoiseMeter splits these into 16 mechanistic categories, reports them per cell, and adds locus scoring and sequence-artifact flags that the aligners do not.
- **Bulk RNA-seq QC (RNA-SeQC, Picard CollectRnaSeqMetrics)** count bases by region for bulk libraries. scNoiseMeter is single-cell-aware (per-cell and per-cluster) and read-resolved.
- **Ambient-RNA removal (SoupX, DecontX, CellBender)** estimate and subtract a single contamination fraction per cell from the count matrix. scNoiseMeter does not remove anything; it classifies reads so you can see the mechanism (internal priming vs TSO invasion vs intergenic), and it works from the BAM, not the matrix.
- **Empty droplets (DropletUtils emptyDrops) and doublets (Scrublet, DoubletFinder)** answer "is this barcode a real cell?" scNoiseMeter answers "within the real cells, what kind of noise is in their reads?" The two are complementary upstream and downstream of each other.
- **Long-read isoform QC (SQANTI3, pigeon, FLAMES, FLAIR, IsoQuant, TALON)** validate *assembled isoforms* and flag intra-priming and RT-switching post-assembly. scNoiseMeter operates pre-assembly on individual alignments, so it can feed cleaner reads into those pipelines, and its TSO-invasion and TSO-concatemer detection has no direct equivalent there.
- **Long-read run QC (LongQC, NanoPlot, pycoQC)** summarize yield, length, and quality. scNoiseMeter is orthogonal: it classifies artifact type, not run quality.

### Strengths

- One tool, one taxonomy, across long- and short-read and plate-based assays.
- Per-read and per-base resolution, aggregated to per-cell and per-cluster.
- A statistical model for intergenic loci instead of a flat "intergenic = noise" rule.
- Dedicated, configurable TSO and internal-priming detection, including the TSO-concatemer metric from the uMRT literature.
- Two principled noise bounds rather than one fragile number.

### Honest limitations

- It classifies; it does not correct or remove reads.
- `intronic_pure` and `intronic_boundary` cannot be separated from genuine pre-mRNA at the read level, which is exactly why they appear in conservative but not strict noise.
- Chimeric calls mix technical chimeras with real gene fusions; interpretation needs biological context.
- The intergenic Poisson model assumes a homogeneous background, which is an approximation near highly transcribed regions.
- Human GRCh38/hg38 only; other genomes warn but are not supported.
- It has not been benchmarked against a sorted or spike-in ground truth, an open problem shared by every tool in this space.

---

## 7. Tuning

Every threshold lives in one file, `scnoisemeter/constants.py`, and the common ones are exposed as CLI flags (`--chimeric-distance`, `--tso`, `--tso-min-match`, `--no-polyg-tso`, `--platform`, and the annotation flags). The detailed companion document maps each rule to its exact code location and says what to change to alter its behavior. See [methods_annotated.md](methods_annotated.md).
