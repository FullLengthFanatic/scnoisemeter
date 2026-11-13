# scNoiseMeter: white paper

**A barcode-aware alignment artifact and read-distribution QC framework for single-cell RNA-seq**

Version 0.10.0. This document explains the scientific position of scNoiseMeter: what the
categories mean, why they are drawn where they are, and what a reader should and should not
conclude from the numbers. For implementation details see
[methods_annotated.md](methods_annotated.md); for command and output reference see
[documentation.md](documentation.md).

## Abstract

RNA-seq alignments outside an annotated, sense exonic model are frequently described as “noise.”
That shorthand is useful operationally but scientifically unsafe. Intronic RNA supports
transcriptional dynamics; antisense transcription may be regulated; intergenic signal may
represent incomplete annotation; split alignments may represent real fusions; and mitochondrial
signal varies with biology as well as sample quality. A BAM alone rarely identifies a causal
mechanism.

scNoiseMeter therefore treats technical-noise analysis as an evidence-labeling and composition
problem. Every mapped primary alignment receives one mutually exclusive location/alignment
category, its aligned reference bases are partitioned without double-counting, and independent
context flags record evidence such as internal priming, an untemplated polyA tail, TSO sequence,
non-canonical splice boundaries, reverse-transcriptase template-switch repeats, low mapping
confidence, or pair discordance. Where a flag describes the genome rather than the molecule, a
companion flag describing the molecule is reported beside it, because the two disagree in
precisely the cases that matter. Results are aggregated by sample, cell, cluster, read length,
or plate well.

The tool reports a broad non-canonical composition and a narrower artifact-candidate composition.
Neither is presented as an estimate of contamination. Every run additionally checks its own
finished output against the structural properties claimed here and refuses to write a result set
that violates them. The framework is designed to be auditable
and comparable across long-read, droplet short-read, and plate protocols while remaining explicit
about assay-specific assumptions.

## 1. Scientific motivation

No single alignment category is synonymous with a technical artifact.

- Intronic signal can contain nascent RNA. The distinction between spliced and unspliced RNA is
  the basis of [RNA velocity](https://www.nature.com/articles/s41586-018-0414-6), and bulk work
  has also shown that intronic and exonic reads carry distinct regulatory information.
- Internal oligo(dT) priming is a real, systematic bias across bulk and single-cell data,
  especially near genomic A-rich sequence
  ([Svoboda et al., 2022](https://academic.oup.com/nargab/article/4/2/lqac035/6592171)).
- Long-read transcript models require post-alignment curation of junctions, ends, and artifact
  evidence; [SQANTI3](https://www.nature.com/articles/s41592-024-02229-2) operates at the
  transcript-model level.
- Ambient RNA tools such as
  [SoupX](https://academic.oup.com/gigascience/article/9/12/giaa151/6049831) estimate
  contamination from count-matrix and empty-droplet evidence. An alignment-category tool cannot
  reproduce that inference from mapped location alone.
- Large benchmarking efforts such as
  [LRGASP](https://www.nature.com/articles/s41592-024-02298-3) demonstrate that long-read
  transcript reconstruction and quantification depend on both protocol and analysis method.

This motivates a narrower claim: scNoiseMeter describes alignment composition and flags candidate
artifacts. It does not infer that every non-exonic read is technical, nor does it remove
molecules.

## 2. Design principles

### 2.1 Separate observation from interpretation

The primary categories describe what is observed: mitochondrial mapping, multimapping, exonic
orientation, intronic structure, same-strand gene ambiguity, intergenic enrichment, or
split-alignment evidence. Independent flags describe additional evidence. Causal terms are
reserved for cases with a specific context signature and are still qualified as candidates.

### 2.2 Preserve exact denominators

The BAM index supplies mapped and unmapped totals. Primary, secondary, supplementary, QC-fail,
and duplicate counts are reported separately. Genomic category fractions use classified
mapped-primary alignments. Unmapped fraction uses the full mapped-plus-unmapped BAM-index
denominator.

At base level, every aligned reference block is split at annotation boundaries. Exactly one
category wins each atomic interval. This prevents the former failure mode where overlapping
exons, introns, or genes counted the same base multiple times.

### 2.3 Keep strand as evidence

Opposite-strand genes are not automatically ambiguous. Ambiguity is formed from same-strand
shared gene regions, while exonic sense and antisense are classified relative to read
orientation. Users can state `stranded` or `unstranded`; automatic inference is intentionally
conservative because aligner names do not identify library strandedness.

### 2.4 Avoid pseudo-replication

Pre/post-filter BAMs are dependent, and individual reads are nested within molecules and cells.
scNoiseMeter does not apply an independent-samples chi-square test to read counts. It matches
read keys exactly for retention/transitions and uses cells as paired units for a bootstrap
interval on median per-cell differences.

### 2.5 Make missing evidence explicit

Endpoint anchoring is absent when the required atlas is absent. `NM` summaries are absent when
the tag is unavailable. UMI sequence diversity is absent for categories whose exact UMI
membership cannot be reconstructed after reservoir-based intergenic reclassification. Per-cell
summaries are absent for a barcode-free BAM, and the unmapped fraction is absent when the index
carries no counters. Missing is preferable to a fabricated proxy.

The same rule applies to a resource that is present but unusable. An atlas whose contig names
share nothing with the BAM would make every proximity test fail while leaving the metric
non-empty, so it is discarded with an error and the derived fractions are reported as absent
rather than as measured zeros.

## 3. The classification model

Unmapped, secondary, and supplementary records are counted but never assigned a genomic
category: they carry no comparable location. Every mapped primary alignment then passes through
one hierarchy, and the first rule that fires wins.

### 3.1 One read, end to end

The hierarchy is easiest to read by following a single alignment through it. Consider a 940 bp
PacBio read aligned across 1,240 bp of the genome, `chr1:1,000,100-1,001,340`, with CIGAR
`620M300N320M` and tags `NH:i:1`,
`CB:Z:AAACCTGAGAAACCAT`, `UB:Z:TTTGGCCAAT`, no `SA`.

1. **Barcode gate.** `CB` is present and on the whitelist, so the read is attributed to that
   cell. A missing `CB`, or one absent from the whitelist, ends here as `unassigned` — still
   counted in the denominator, because it is a real read that simply cannot be attributed.
2. **Multimapper evidence.** `NH:i:1`, so this is a single-hit alignment. `NH:i:3` would end
   here as `multimapper`, held out of the location categories entirely.
3. **Mitochondrial contig.** `chr1` is not `chrM`/`MT`, so the read continues. The mitochondrial
   check sits ahead of annotation, which matters because the mitochondrial genome is
   almost entirely exonic: tested later, nearly every mitochondrial read would be absorbed into
   `exonic_sense` and the fraction people watch for lysis would disappear.
4. **Chimera evidence.** No `SA` tag and no discordant mate, so the read continues.
5. **Atomic annotation overlap.** 620 aligned bases fall in one exon of a `+`-strand gene, 320
   in the next exon of the same gene, and the 300 `N` bases are skipped rather than aligned. Read
   orientation matches gene orientation, so all 940 aligned bases are sense-exonic and the
   category is `exonic_sense`. When an alignment mixes several kinds of bases, the kind holding
   most of them names the read.
6. **Base partition.** 940 aligned reference bases are attributed to `exonic_sense`. The 300 `N`
   bases are not aligned reference bases and enter no category — this is why base fractions and
   read fractions differ, and why they differ more for long reads than short ones.
7. **Flags, computed independently of the category.** The `N` interval is checked as a
   `(donor, acceptor)` pair in transcript orientation, found to be GT–AG and present as an exact
   junction of an annotated transcript, so no non-canonical flag is raised. The 5′ soft clip is
   screened for the TSO motif, and the reference immediately downstream of the strand-correct 3′
   end is screened for a genomic A-run.

Small changes to the same read change the outcome in instructive ways. Remove the `N` and the
300 intronic bases are now aligned. They are counted as `intronic_boundary` bases, but the read
stays `exonic_sense`, because 940 of its 1,240 aligned bases are exonic. Shift that unspliced
read until most of it lies in the intron and it becomes `intronic_boundary`; place it wholly
inside the intron and it is `intronic_pure`. Flip the read orientation and it becomes `exonic_antisense`. Add an `SA` tag pointing to
another contig and it becomes `chimeric` before annotation is ever consulted — a chimeric
alignment has no single meaningful location.

### 3.2 The categories and what they usually mean

Seventeen categories are reported. The rule is what the code tests; the interpretation is what a
reader should have in mind, and it is never certain from a BAM alone. The rules for the exonic,
intronic and intergenic categories describe kinds of aligned bases: a read takes the category
that holds most of its aligned bases (§3.1), and the base-level fractions count each kind
exactly.

**Exonic**

| Category | Rule | What it usually means |
|---|---|---|
| `exonic_sense` | Aligned bases in an annotated exon, read on the gene's strand | The expected signal. In a healthy stranded library this dominates. |
| `exonic_antisense` | Same, opposite strand | In a **stranded** library: candidate strand-switching during reverse transcription or template switching, or genuine regulated antisense transcription. In an **unstranded** library it is expected cDNA signal and is excluded from the broad composition. |

**Intronic** — these three separate the spliced/unspliced question that RNA velocity depends on.

| Category | Rule | What it usually means |
|---|---|---|
| `intronic_jxnspan` | Intronic bases on an alignment containing CIGAR `N` | A spliced molecule that still covers intronic sequence: intron retention, or a novel junction. Junction canonicality is flagged separately rather than folded in. |
| `intronic_pure` | Entirely within an intron body | Nascent/pre-mRNA. Rises sharply and expectedly in single-nuclei preparations. Also where internal priming lands when intronic sequence is A-rich. |
| `intronic_boundary` | Spans an exon–intron boundary with no `N` | Incomplete reverse transcription, or a genuinely unspliced molecule caught mid-intron. |

**Intergenic** — assigned in a second pass, because a single read cannot tell you whether a locus
is enriched. See §3.5.

| Category | Rule | What it usually means |
|---|---|---|
| `intergenic_sparse` | Intergenic, below the enrichment screen | Diffuse background: mismapping, residual genomic DNA, or annotation that does not cover a real transcript. |
| `intergenic_repeat` | ≥50% of sampled aligned span overlaps supplied repeat intervals | Multi-copy sequence, where placement is intrinsically uncertain. Note that Alu poly-A tails make repeats also the densest source of genomic A-runs, so this category and internal priming overlap physically. |
| `intergenic_hotspot` | Enriched, monoexonic, A-rich strand-aware 3′ context, and *not* near a strand-matched annotated polyA site | The internal oligo(dT) priming candidate. This is one of only two categories in the artifact-candidate set. |
| `intergenic_novel` | Enriched, dominant-strand, multiple real barcodes, plus splice or annotated-polyA support | Candidate unannotated transcript. A label for follow-up, not a discovery claim. |
| `intergenic_enriched` | Enriched, none of the above | Deliberately named for the statistic rather than a cause. Calling this a hotspot would convert enrichment into a mechanism. |

**Alignment-level and unattributable**

| Category | Rule | What it usually means |
|---|---|---|
| `chimeric` | `SA` evidence that is inter-contig, strand-discordant, or incompatible in query-vs-genomic order; extreme paired insert size as a fallback | Technical chimera or a real fusion. These are not separable from an `SA` tag, which is why the category is not called "artifact". |
| `mitochondrial` | Maps to `chrM`/`MT` | Both biology and sample quality. High values are conventionally read as stress or lysis, but tissue identity matters more than a single threshold. |
| `multimapper` | Explicit `NH > 1` | Ambiguous placement, held out of the location categories. **`minimap2` and `pbmm2` do not emit `NH`**, so this category is empty by construction on most long-read BAMs and those alignments are counted under their genomic category instead. This shifts the denominator relative to a STAR-derived BAM and is the single largest obstacle to comparing platforms. |
| `ambiguous_cod_cod`, `ambiguous_cod_ncod` | Aligned bases fall where two or more **same-strand** genes overlap | Genuinely unattributable to one gene. Split by biotype pairing because a coding/coding overlap and a coding/lncRNA overlap warrant different scepticism. Measured at exon level, which keeps it small — on GENCODE v47 the two shared sets are 0.9% and 1.5% of exonic sequence. |
| `unassigned` | No usable cell barcode | A real read that cannot be attributed to a cell. Counted in the denominator. Reads with no barcode at all are left out of the per-cell table. |

### 3.3 Evidence flags are independent of category

A read has exactly one category but any number of flags, and the flags are what carry mechanistic
evidence:

- **TSO invasion** — the template-switch oligo motif in the soft-clipped end where it is
  molecularly plausible, with one substitution tolerated.
- **TSO concatemer** — the motif appearing more than once, indicating oligo-oligo product.
- **Internal polyA priming** — a genomic A-run in the strand-correct reference context
  immediately beyond the read's 3′ end. Read this carefully: it is a statement about the
  *place in the genome*, not about the molecule. It says an oligo-dT primer *could* have
  annealed here. It cannot say whether this particular molecule was primed there or arrived
  by genuine cleavage and polyadenylation.
- **Untemplated 3′ polyA tail** — an A-run in the read's own soft-clipped end, outside the
  aligned block. Those bases are in the read but not in the reference, so the molecule brought
  them: this is evidence about the molecule rather than the locus. The two flags are most
  informative together, and the table in §3.3.1 explains why.
- **Non-canonical junction** — a CIGAR `N` interval whose donor/acceptor pair is neither an
  annotated junction nor GT–AG/GC–AG/AT–AC.
- **RT template-switch repeat** — the two ends of one intron share a short identical sequence.
  Reverse transcriptase can jump between two copies of the same short sequence, stitching
  together fragments that were never adjacent; the join then looks like a splice junction.
  This is deliberately independent of the canonicality test above, because such an artifact can
  land on a perfectly canonical GT–AG boundary and be invisible to a motif check. Real junctions
  in repetitive sequence also carry direct repeats, so this is evidence and never a category.
- **Low MAPQ** — reported descriptively, never used to filter or reclassify, because MAPQ scales
  differ between aligners.
- **Discordant pair** — a missing or improper mate, which is a flag and not an automatic chimera.

Keeping these orthogonal to the taxonomy is what allows a read to be, for example,
`exonic_sense` *and* internally primed — a combination that a single-label scheme would have to
discard.

#### 3.3.1 Why two polyA flags rather than one

Internal oligo-dT priming is the artifact most often blamed for spurious 3′ signal, and the usual
test for it is genomic: look at the reference sequence just past the read's 3′ end and ask whether
it is A-rich. scNoiseMeter does that, but the test is weaker than its popularity suggests.
Measured against the PolyASite atlas (each site taken at its representative cleavage position),
the median A-content in the 20 bp downstream window is
**30% at annotated polyA sites and 30% at random genomic positions**. The genomic context
barely distinguishes the two, because A-rich stretches are simply common.

What does distinguish them is whether the A's the read carries are in the genome. At a real
cleavage site the polyA tail is added enzymatically after transcription, so those bases appear in
the read but not in the reference, and the aligner pushes them into a soft clip. At an internal
priming event the A's the read carries are the genome's own, and they align.

Reading the two flags together therefore separates four situations that a single flag collapses
into two:

| genomic A-run downstream | untemplated tail in the soft clip | what it means |
|---|---|---|
| yes | no | internal priming — the classical case, now with corroboration |
| yes | yes | a real 3′ end that happens to sit in A-rich sequence; the genomic test alone calls this an artifact |
| no | yes | a real 3′ end |
| no | no | no evidence either way |

Only soft clips are read. If an aligner hard-clips, the bases are removed from the record
entirely, and scNoiseMeter reports no tail — which means "not observable here", not "no tail on
this molecule". The distinction is recorded rather than papered over.

A second choice is worth stating because it is easy to get wrong. The tail is decided by the
length of the A-run alone, not by how A-rich the clip is overall. A genuine 12 bp tail sitting
inside a 60 bp clip that also holds adapter or TSO sequence is only about 20% A, so requiring
both a run and a high fraction would quietly demand a tail of roughly 45 bp and miss most real
ones.

### 3.4 Chimeric and splice evidence

The `SA` tag is interpreted structurally. Inter-contig and strand-discordant splits are
candidates. Same-contig, same-strand splits are candidates only when aligned query order is
incompatible with genomic order. Genomic distance alone is not evidence, because RNA alignments
routinely span long introns — a distance cutoff would reclassify ordinary splicing as chimerism.
This is consistent with the general role of splice-aware aligners such as
[STAR](https://academic.oup.com/bioinformatics/article/29/1/15/272537), which detect both splice
junctions and chimeric transcripts, and with recent work on long-read foldback artifacts
([Breakinator](https://link.springer.com/article/10.1186/s12864-025-12492-y)).

For splice evidence, CIGAR `N` operations are converted to exact
`(donor boundary, acceptor boundary)` pairs. Exact transcript junctions from the GTF are
accepted; otherwise both dinucleotides are checked in transcript orientation. A donor-only test
would overcall canonicality, and a position-only splice-site set could accept a novel acceptor
merely because some other transcript has a donor there.

### 3.5 Intergenic evidence and why it needs a second pass

Whether an intergenic read matters depends on how many others land beside it, which no per-read
rule can see. The second pass therefore assigns strand-correct 3′ coordinates to predefined
500 bp windows across the full non-mitochondrial intergenic complement — defined from BAM
reference lengths and the union of GTF gene bodies, so it includes sequence after the last
annotated gene — and applies a global intergenic read-density Poisson model.

Two properties of this screen should be understood before its output is used:

- It is a **global-rate** enrichment screen, not a peak caller. It has no local-background, GC,
  or mappability correction, so it is least reliable exactly where the local rate departs most
  from the genome mean — near highly transcribed regions and in poorly mappable sequence.
- The correction spans every window in the intergenic complement rather than only windows
  containing a read, and the test is performed on the intergenic reservoir sample rather than
  the full read set. Both make it conservative, and the second makes promotion depth-dependent:
  a library far above the reservoir cap is tested at lower effective depth than one below it.
  (Before 0.9 the correction counted every window of every contig, including all the genic
  sequence where no test is ever performed — roughly a 2.2–2.9× overcorrection on GRCh38 that
  suppressed real enrichment. Intergenic category counts are not comparable across that
  boundary.)

Each window also reports **`local_read_share`**: how much of the intergenic read mass within
10 kb the window itself holds. This answers a question significance cannot. A window can clear
the genome-wide rate because it is a sharp isolated pile-up, or because it sits inside a broad
smear of coverage that happens to exceed the average everywhere. A share near 1 is the first; a
low share is the second. They warrant different follow-up, and the p-value is identical.

Two further columns describe each window without changing how it is
classified. **`n_molecules`** counts distinct barcode/UMI pairs rather than
reads, because thirty reads of one over-amplified molecule is one observation
and not thirty; it is reported as an exact count, since deciding whether a
one-base-different UMI is an error or a second molecule depends on the
chemistry rather than on anything visible here. **`strand_asymmetry_pvalue`**
asks whether a window is more one-sided than a coin would produce, which the
existing strand fraction cannot: three reads out of four and seventy-five out
of a hundred are the same fraction, and only the second is evidence. Neither
column promotes anything yet. They exist so that the open question — whether
promotion should count molecules, and which multiple-testing procedure to use —
can be settled from measurements on real samples rather than by argument.

**A caution about the polyA atlas.** Two of the rules
below turn on proximity to an annotated polyA site, and they turn on it in opposite directions:
`intergenic_hotspot` requires that a site is *not* nearby, while `intergenic_novel` accepts a
nearby site as positive evidence. That makes the density of the atlas a silent parameter.
PolyASite 3.0 contains 18.4 million records, many supported in only a few percent of samples,
and loaded whole it puts **37.4% of the genome within the 50 bp proximity window**. At
that density the test is closer to a constant than a discriminator: hotspots are vetoed across
a third of the genome, and the novel-gene evidence route is nearly free. Filtering the atlas on
its own per-site support column restores it: 25% support leaves 6.3% of the genome within the
window and 50% leaves 2.4%. Every run measures and reports the density of the atlas it loaded.
Since 0.9.1 the default is 25% support. Earlier versions loaded the atlas whole, and in the
synthetic benchmark that setting reported 6 of 10 internal-priming loci as `intergenic_novel`.

Enriched windows are then interpreted by the independent evidence in §3.2. The ordering matters:
repeat overlap is tested before enrichment, so a window that is both repeat-derived and A-rich is
reported as `intergenic_repeat`. Because Alu poly-A tails are the densest genomic A-runs, some
genuine internal-priming loci are therefore labelled by their sequence context rather than their
mechanism. The per-locus table carries `repeat_overlap_fraction` alongside the category so this
is recoverable.

## 4. Aggregate metrics

The categories above are the primitives; the aggregates below are conveniences built from them,
and each is a composition rather than an estimate.

### Broad non-canonical composition

The broad set contains exonic antisense (unless unstranded), intronic pure and boundary,
sparse/repeat/enriched/hotspot intergenic reads, and chimeric alignments. Note what it excludes:
`intronic_jxnspan`, because a spliced molecule covering intronic sequence is ordinary
intron retention; `intergenic_novel`, because its evidence points at annotation incompleteness
rather than technical failure; `mitochondrial`, because it is biology as often as quality; and
`multimapper` and `ambiguous`, because uncertain placement is not the same as a wrong molecule.

It is useful for comparing preparations or filters on the same protocol. It is not a
contamination rate, and its components should be reported alongside it.

### Artifact-candidate composition

The narrow set contains only `intergenic_hotspot` and `chimeric` — the two categories with
positive alignment or context evidence for a specific mechanism. Even this is not ground truth:
A-rich contexts coincide with real 3′ ends, and `SA` evidence represents real fusions as well as
technical ones.

### UMI sequence diversity

Unique UMI strings divided by reads, within a cell and category. This is **not** molecule
complexity: UMIs are not grouped by gene or genomic molecule and sequencing errors are not
corrected. A deprecated `umi_complexity_*` alias is emitted alongside the accurate name.

### Endpoint anchoring

Read length alone is not evidence of transcript completeness — a long alignment may start late,
end early, retain introns, or cross a fusion. scNoiseMeter instead measures strand-matched
proximity to a TSS/CAGE atlas and to a polyA atlas **on the same read**, and reports
`both_ends_anchored_frac` only when both resources exist. The deprecated
`full_length_read_frac` is exactly that same-read both-end fraction, never a length fallback.
This is consistent with work emphasising the difficulty of assigning trustworthy long-read
transcript ends, and with cap- and 3′-enriched approaches such as
[CapTrap-seq](https://www.nature.com/articles/s41467-024-49523-3).

Proximity to an atlas record is a weak anchor when the atlas is dense (§3.5), so two stricter
readings are available alongside it. `three_prime_at_pas_site_frac` counts only proximity to
records the atlas itself annotates with a polyadenylation signal hexamer — a real cleavage site
usually has one 10–40 nt upstream, though a strict hexamer also occurs at 9.4% of random genomic
positions, so it is corroboration rather than proof. And `--anchor-veto-internal-priming`
additionally requires that a read carry no internal-priming signature before it counts as
anchored, on the reasoning that proximity alone should not outvote contrary evidence on the same
read. Both are reported, and the veto is off by default because it lowers the headline fraction.

## 5. Reading the output

The numbers below are from the example PacBio Kinnex run whose report is shown in the README
screenshots: 37.4 M reads, 2,047 cells, aligned with `pbmm2`, annotated with GENCODE v45, with
PolyASite 3.0 and FANTOM5 CAGE auto-downloaded, run with version 0.8.0. They are one sample on one protocol, not
reference values — but they show what the fields look like together, which no field description
does on its own.

| Field | Value | How to read it |
|---|---|---|
| `exonic_sense` (reads / bases) | 64.76% / 63.06% | The bulk of the library, as expected for a stranded long-read cDNA protocol. |
| `broad_noncanonical_read_frac` | 23.22% | Not a contamination estimate. Dominated here by `intronic_pure` at 14.00%. |
| `artifact_candidate_read_frac` | 3.70% | Almost entirely the chimeric rate: `intergenic_hotspot` adds only 0.002% in this sample. |
| `strand_concordance` | 98.86% | Near 1 confirms a stranded library. Near 0.5 would mean unstranded, and would invalidate antisense comparisons. |
| `intronic_pure` | 14.00% | Substantial, and expected for whole-cell long-read data; would be higher still for single nuclei. |
| `three_prime_anchored_frac` | 99.62% | Almost every read's 3′ end sits at an annotated polyA site — the signature of an oligo(dT)-primed protocol working as designed. |
| `five_prime_anchored_frac` | 84.02% | Lower than the 3′ figure, which is the normal asymmetry: 5′ completeness depends on reverse transcription reaching the cap. |
| `both_ends_anchored_frac` | 83.84% | Candidate complete molecules. Not a transcript-compatibility claim. |
| `polya_priming_frac` (`multiqc.json`) | 2.31% | The internal-priming flag. Sits alongside an `intergenic_hotspot` of 0.002%, i.e. flagged reads are not forming enriched intergenic loci here. |
| `tso_invasion_frac` (`multiqc.json`) | 0.0045% | Low. This metric is highly sensitive to detection settings — see §7. |
| `mitochondrial` | 3.77% | Interpret against the tissue, not a universal threshold. |

Patterns that should attract attention, rather than any single threshold:

- **`strand_concordance` near 0.5 with a `stranded` declaration.** The declaration is wrong, or
  the library is not what it is thought to be. Every antisense and endpoint number downstream is
  then measuring something else.
- **`intergenic_hotspot` rising with `polya_priming_frac`.** The two agreeing is the coherent
  internal-priming signature; the flag rising alone is likelier to be threshold sensitivity.
- **`three_prime_anchored_frac` far below the 5′ figure** in an oligo(dT) protocol. That
  inverts the expected asymmetry and points at the atlas or the reference, not the library.
- **`multimapper` exactly 0 on a long-read BAM.** Expected, not reassuring — see §3.2.
- **`chimeric` climbing between a pre- and post-filter pair.** Filtering should not create split
  alignments; if it appears to, check whether read names were rewritten and the pairing is real.
- **A category that is absent rather than zero.** This means the evidence was never available.
  Absent and zero are different claims and the output distinguishes them deliberately.
- **`polya_priming_frac` high while `polya_tail_frac` is comparably high in the same reads.**
  Read together (§3.3.1), that is a population of real 3′ ends sitting in A-rich sequence, not
  an internal-priming problem. The genomic flag alone would have said the opposite.
- **`rt_switch_repeat_frac` elevated while non-canonical junctions are not.** Template-switch
  artifacts can occupy canonical GT–AG boundaries, so the junction-motif check will not see them.
  This is the pairing that makes the direct-repeat flag worth having.
- **A warning that the loaded polyA atlas covers a large share of the genome.** Reported by every
  run. Above roughly a fifth, treat `intergenic_hotspot` and `intergenic_novel` as compromised in
  opposite directions rather than as measurements (§3.5).

## 6. Protocol handling

Sequencer, aligner, library chemistry, and strandedness are different concepts. Minimap2 is used
beyond ONT; STAR is used beyond Smart-seq; a 10x library can be sequenced on several short-read
instruments. Auto-detection uses explicit pipeline hints where available and otherwise returns
generic or unknown states with warnings.

TSO screening follows protocol evidence: a 10x motif for ONT/10x contexts, a SMART/PacBio motif
for PacBio/Smart-seq, and no motif for generic or unknown inputs unless `--tso` is supplied.
Approximate matching tolerates one substitution at the molecularly appropriate soft-clipped end.
The poly-G shortcut is restricted to the built-in 10x motif; on two-colour Illumina chemistry a
terminal G-homopolymer is a basecalling artifact rather than a TSO, so this option should be used
with care there.

Unstranded libraries suppress orientation-dependent endpoint metrics and remove exonic antisense
from the broad composition. Explicit `--library-strand` is recommended over inference.

## 7. Comparing across samples

Two comparisons are possible and they are not interchangeable.

A **nested** comparison asks what a processing step did to a fixed set of reads, and is only
defined when both files contain the same reads. Read-level retention and category transitions
belong here. This is narrower than "before and after a step": some deduplicators rewrite read
names, so a genuine pre/post pair can be unmatchable at the read level even though both files
describe the same library.

A **cohort** comparison asks how independent samples differ, and can only be distributional.
Four limits apply, and each is a way a cross-sample figure can invent a result.

**Metric definitions change between tool versions.** The same BD46 alignment file reports a TSO
invasion rate of 2.33e-03 under v0.6.1 and 3.28e-06 under v0.7.2 — a 700-fold difference caused
entirely by a change in how the flag is detected. Pooling result sets from different versions
produces differences that look biological and are not. `cohort` flags version mixtures loudly
for this reason.

**Annotation completeness changes what "intergenic" means.** Samples annotated against different
GENCODE releases are not comparable in any category that depends on gene coverage.

**Strandedness changes what a category measures.** In a non-stranded protocol, exonic-antisense
reads are expected cDNA signal and strand concordance near 50% is the design, so neither quantity
can be ranked against a stranded sample.

**Aligner capability changes the denominator.** A BAM without `NH` has an empty multimapper
category, so its location fractions are computed over a slightly different read set than a
STAR-derived sample's.

A metric that a sample never recorded is a gap in the record, not a measurement of zero, and is
reported as absent.

## 8. Novelty and relationship to existing tools

The individual ingredients are established: genomic partitioning, splice motif checks,
internal-priming context, endpoint atlases, per-cell summaries, and split-alignment evidence.
Bulk genomic-region composition in particular has been available for years through RSeQC's
`read_distribution.py` and Picard's `CollectRnaSeqMetrics`. The useful novelty here is narrower
and should be stated as such: barcode-aware, base-exact, protocol-aware integration of those
ingredients at the individual-alignment level, under one taxonomy, across droplet, long-read and
plate workflows.

scNoiseMeter complements rather than replaces:

- aligner summaries (mapping rate and broad genomic regions);
- ambient-RNA models such as SoupX, DecontX and CellBender;
- cell-calling and doublet detection;
- isoform-model QC such as SQANTI3 and pigeon, which share several of the same underlying signals
  (RT-switching evidence, intra-priming, junction canonicality) at transcript-model level.
  Where a signal is shared, the definition here is deliberately aligned with theirs so the
  numbers can be read side by side: `--priming-criterion perc_a` is their
  `perc_A_downstream_TTS` convention, and the direct-repeat test is the same RT-switch signature
  they apply to transcript models, applied instead to individual alignments;
- read-structure artifact classifiers, which work from the read rather than the alignment and can
  therefore see products — such as TSO-TSO concatemers — that align poorly or not at all;
- pre-alignment chimera removal, which sees adapter-bridged reads that never surface as `SA`
  evidence;
- dedicated internal-priming filters such as polyAfilter, which act on the count matrix rather
  than describing the alignment;
- fusion callers, which apply gene-level evidence and filtering well beyond an `SA` tag.

The fixed-window intergenic screen, same-read two-end anchoring, exact pre/post read retention,
and the separation of descriptive composition from artifact-candidate evidence are the most
distinctive elements. They should still be externally benchmarked.

## 9. Validation status and limitations

The test suite covers category invariants, overlapping annotation, strand-aware ambiguity, exact
splice boundaries and motifs, reservoir merging, site-strand loading, plate relabeling,
comparison matching, and an end-to-end synthetic BAM/GTF path. CI runs supported Python versions.

The synthetic benchmark constructs alignments directly rather than through an aligner, so its
ground truth is generated by the same assumptions being tested. Per-category accuracy against it
measures self-consistency, and should be read that way.

A test suite, however, only proves the code behaves on inputs someone thought to construct. Every
run therefore also checks the finished result set against the properties this document claims for
it, and writes the verdicts and the observed numbers to `<sample>.qc.json`. Those properties are
structural rather than statistical, which is what makes them checkable at all: each composition
sums to one; every fraction lies between zero and one; the artifact-candidate composition cannot
exceed the broad non-canonical one, because its categories are a subset of them; both-end
anchoring cannot exceed either single end, because it is their conjunction on the same reads; no
flag count exceeds the denominator it is reported over; unique UMI strings cannot outnumber the
reads carrying them. A violation means a number about to be written is wrong, so the run stops
rather than emitting it. The point is the failure mode being guarded against: a broken invariant
would not announce itself as a crash — it would look like a slightly wrong figure in a table
somebody later cites.

Two checks are warnings rather than errors, because a documented configuration legitimately
breaks them. Reservoir-estimated intergenic reclassification scales promoted totals from a
uniform sub-sample, so exact conservation is not expected to survive rounding; and the anchoring
veto lowers the unrestricted 3′ fraction without lowering the hexamer-restricted subset.

Two configuration errors that previously produced a complete, plausible-looking result set now
stop the run instead. A barcode whitelist from the wrong chemistry leaves nearly every read
unassigned, and the run used to finish and report a sample with a handful of cells and an
unremarkable composition — a failure that arrives disguised as a finding. And a site atlas built
on a different genome assembly shares every contig name with the BAM while every coordinate is
displaced, so the existing name check passes and every proximity test silently returns nothing;
positions falling beyond the end of their contig are the detectable part of that, and are now
tested for. Both can be downgraded to warnings with `--no-strict`.

Important remaining limitations:

- no experimentally sorted positive/negative ground truth across protocols;
- no local-background or mappability model for intergenic enrichment;
- multimapping undetectable where the aligner omits `NH`;
- no gene-aware molecule deduplication, and intergenic promotion counts reads rather than
  UMI-collapsed molecules, so one heavily amplified molecule can carry a window over threshold;
- the untemplated-tail flag is unavailable on hard-clipped BAMs, where the bases are absent from
  the record;
- the RT direct-repeat flag cannot separate a template-switch artifact from a genuine junction in
  repetitive sequence;
- no read-level NUMT inference from competing alignments;
- no distinction between real fusions and technical chimeras;
- duplicate- and QC-fail-flagged records are counted but still classified;
- no correction or removal of reads;
- automatic reference downloads currently focus on GRCh38.

Accordingly, scNoiseMeter should be used for QC, method comparison, and hypothesis generation.
Candidate artifact reads or loci should be validated with orthogonal evidence before exclusion or
biological interpretation.

## 10. Reproducible reporting

A defensible methods section should state:

- scNoiseMeter version/commit;
- BAM generation and filtering stage;
- GTF, genome, repeat, polyA, and TSS resource versions;
- platform and library strandedness supplied to the CLI;
- barcode/UMI tags, whitelist, and called-cell filter;
- custom TSO sequence and whether poly-G checking was enabled;
- the polyA-atlas support threshold (`--polya-min-support`), because it changes what
  "near an annotated polyA site" means and therefore both intergenic rules that depend on it;
- the internal-priming criterion (`--priming-criterion`), and whether the anchoring veto was on;
- seed, thread count, and whether intergenic reclassification was reservoir-estimated;
- which denominator and aggregate metric was reported.

Do not relabel `broad_noncanonical_read_frac` as percent contamination. Report its component
categories and the narrower artifact-candidate composition alongside it.

Every run writes `<sample>.run_info.json` next to the metrics, carrying the version, platform,
strandedness, annotation sources, and the priming settings above. The metrics table itself holds only numbers, so this file is
the record that makes an archived result set interpretable later, and it is what the cross-sample
comparison uses to decide whether two samples may be placed side by side at all.
