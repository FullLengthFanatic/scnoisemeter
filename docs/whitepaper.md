# scNoiseMeter: white paper

**A barcode-aware alignment artifact and read-distribution QC framework for single-cell RNA-seq**

Version 0.8.0. This document explains the scientific position of scNoiseMeter. For implementation details see [methods_annotated.md](methods_annotated.md); for command and output reference see [documentation.md](documentation.md).

## Abstract

RNA-seq alignments outside an annotated, sense exonic model are frequently described as “noise.” That shorthand is useful operationally but scientifically unsafe. Intronic RNA supports transcriptional dynamics; antisense transcription may be regulated; intergenic signal may represent incomplete annotation; split alignments may represent real fusions; and mitochondrial signal varies with biology as well as sample quality. A BAM alone rarely identifies a causal mechanism.

scNoiseMeter therefore treats technical-noise analysis as an evidence-labeling and composition problem. Every mapped primary alignment receives one mutually exclusive location/alignment category, its aligned reference bases are partitioned without double-counting, and independent context flags record evidence such as internal priming, TSO sequence, non-canonical splice boundaries, low mapping confidence, or pair discordance. Results are aggregated by sample, cell, cluster, read length, or plate well.

The tool reports a broad non-canonical composition and a narrower artifact-candidate composition. Neither is presented as an estimate of contamination. The framework is designed to be auditable and comparable across long-read, droplet short-read, and plate protocols while remaining explicit about assay-specific assumptions.

## 1. Scientific motivation

No single alignment category is synonymous with a technical artifact.

- Intronic signal can contain nascent RNA. The distinction between spliced and unspliced RNA is the basis of [RNA velocity](https://www.nature.com/articles/s41586-018-0414-6), and bulk work has also shown that intronic and exonic reads carry distinct regulatory information.
- Internal oligo(dT) priming is a real, systematic bias across bulk and single-cell data, especially near genomic A-rich sequence ([Svoboda et al., 2022](https://academic.oup.com/nargab/article/4/2/lqac035/6592171)).
- Long-read transcript models require post-alignment curation of junctions, ends, and artifact evidence; [SQANTI3](https://www.nature.com/articles/s41592-024-02229-2) operates at the transcript-model level.
- Ambient RNA tools such as [SoupX](https://academic.oup.com/gigascience/article/9/12/giaa151/6049831) estimate contamination from count-matrix and empty-droplet evidence. An alignment-category tool cannot reproduce that inference from mapped location alone.
- Large benchmarking efforts such as [LRGASP](https://www.nature.com/articles/s41592-024-02298-3) demonstrate that long-read transcript reconstruction and quantification depend on both protocol and analysis method.

This motivates a narrower claim: scNoiseMeter describes alignment composition and flags candidate artifacts. It does not infer that every non-exonic read is technical, nor does it remove molecules.

## 2. Design principles

### 2.1 Separate observation from interpretation

The primary categories describe what is observed: mitochondrial mapping, multimapping, exonic orientation, intronic structure, same-strand gene ambiguity, intergenic enrichment, or split-alignment evidence. Independent flags describe additional evidence. Causal terms are reserved for cases with a specific context signature and are still qualified as candidates.

### 2.2 Preserve exact denominators

The BAM index supplies mapped and unmapped totals. Primary, secondary, supplementary, QC-fail, and duplicate counts are reported separately. Genomic category fractions use classified mapped-primary alignments. Unmapped fraction uses the full mapped-plus-unmapped BAM-index denominator.

At base level, every aligned reference block is split at annotation boundaries. Exactly one category wins each atomic interval. This prevents the former failure mode where overlapping exons, introns, or genes counted the same base multiple times.

### 2.3 Keep strand as evidence

Opposite-strand genes are not automatically ambiguous. Ambiguity is formed from same-strand shared gene regions, while exonic sense and antisense are classified relative to read orientation. Users can state `stranded` or `unstranded`; automatic inference is intentionally conservative because aligner names do not identify library strandedness.

### 2.4 Avoid pseudo-replication

Pre/post-filter BAMs are dependent, and individual reads are nested within molecules and cells. scNoiseMeter does not apply an independent-samples chi-square test to read counts. It matches read keys exactly for retention/transitions and uses cells as paired units for a bootstrap interval on median per-cell differences.

### 2.5 Make missing evidence explicit

Endpoint anchoring is absent when the required atlas is absent. `NM` summaries are absent when the tag is unavailable. UMI sequence diversity is absent for categories whose exact UMI membership cannot be reconstructed after reservoir-based intergenic reclassification. Missing is preferable to a fabricated proxy.

## 3. Classification model

Unmapped, secondary, and supplementary records are counted but not assigned genomic categories. A mapped primary read then passes through:

1. barcode/whitelist gate;
2. multimapper evidence;
3. mitochondrial contig check;
4. split/pair chimera evidence;
5. atomic annotation overlap;
6. intergenic second-pass refinement.

The resulting 17 categories are listed in the README.

### Chimeric evidence

The `SA` tag is interpreted structurally. Inter-contig and strand-discordant splits are candidates. Same-contig, same-strand splits are candidates only when aligned query order is incompatible with genomic order. Genomic distance alone is not evidence because RNA alignments routinely span long introns. This is consistent with the general role of splice-aware aligners such as [STAR](https://academic.oup.com/bioinformatics/article/29/1/15/272537), which detect both splice junctions and chimeric transcripts.

For paired cDNA alignments, inter-contig mates or an extreme template length can produce a chimeric candidate. A missing or improper mate produces a discordance flag, not an automatic chimera. Real fusions and technical chimeras remain inseparable without additional biological evidence.

Recent work on long-read foldback and chimera artifacts further supports treating orientation/order as evidence rather than using a fixed genomic gap ([Breakinator](https://link.springer.com/article/10.1186/s12864-025-12492-y)).

### Splice evidence

CIGAR `N` operations are converted to exact `(donor boundary, acceptor boundary)` pairs. Exact transcript junctions extracted from the GTF are accepted. Otherwise both donor and acceptor dinucleotides are checked in transcript orientation against GT–AG, GC–AG, and AT–AC. A donor-only test would overcall canonicality; a position-only set could incorrectly combine donor and acceptor from different transcripts.

### Intergenic evidence

The second pass assigns strand-correct 3′ coordinates to predefined 500 bp windows. BAM reference lengths and the union of GTF gene bodies define the complete non-mitochondrial intergenic complement, including sequence after the final annotated gene. A global intergenic read-density model supplies a Poisson tail probability. The correction includes every possible window, not only windows observed to contain a read. The model is deliberately described as a global-rate enrichment screen; unlike a mature peak caller, it has no local-background, GC, or mappability correction.

Enriched windows are interpreted by independent evidence:

- repeat-derived: at least half of sampled aligned span overlaps repeat intervals;
- internal-priming hotspot candidate: monoexonic, strand-aware genomic polyA-run context, and not near a strand-matched annotated polyA site;
- candidate novel transcript: dominant-strand support, multiple real barcodes, and either splice or annotated polyA-end support;
- enriched unresolved: significant but lacking the evidence above.

Calling the last group “enriched” rather than “hotspot” avoids converting statistical enrichment into a causal artifact label.

## 4. Aggregate metrics

### Broad non-canonical composition

The broad set contains exonic antisense (unless unstranded), intronic pure/boundary, sparse/repeat/enriched/hotspot intergenic reads, and chimeric alignments. It is useful for comparing preparations or filters but can include real biology.

### Artifact-candidate composition

The narrow set contains internal-priming hotspot candidates and chimeric candidates. Even this is not ground truth: A-rich contexts can coincide with real ends and chimeric evidence can represent real fusions.

The old `noise_*` and `noise_*_strict` fields remain compatibility aliases during the 0.7 transition. New work should use the descriptive names.

### UMI sequence diversity

Unique UMI strings divided by reads measures sequence diversity within a cell/category. It is not molecule complexity because scNoiseMeter does not group UMIs by gene or genomic molecule and does not correct sequencing errors. Both the accurate name and a deprecated `umi_complexity_*` alias are emitted.

### Endpoint anchoring

Read length alone is not evidence of transcript completeness. A long genomic alignment may start late, end early, include retained introns, or cross a fusion. scNoiseMeter instead measures strand-matched proximity to a TSS/CAGE atlas and a polyA atlas on the same read. A `both_ends_anchored_frac` is reported only when both resources exist. The deprecated `full_length_read_frac` is exactly this same-read both-end fraction, never a length fallback.

This conservative position is consistent with work emphasizing the difficulty of assigning trustworthy long-read transcript starts and ends and with cap/3′-end-enriched approaches such as [CapTrap-seq](https://www.nature.com/articles/s41467-024-49523-3).

## 5. Protocol handling

Sequencer, aligner, library chemistry, and strandedness are different concepts. Minimap2 is used beyond ONT; STAR is used beyond Smart-seq; a 10x library can be sequenced on different short-read instruments. Auto-detection uses explicit pipeline hints when available and otherwise returns generic or unknown states with warnings.

TSO sequence screening also follows protocol evidence. The CLI uses a 10x motif for ONT/10x contexts, a SMART/PacBio motif for PacBio/Smart-seq, and no motif for generic/unknown inputs unless the user supplies `--tso`. Approximate matching tolerates one substitution at the molecularly appropriate soft-clipped end. The poly-G shortcut is restricted to the built-in 10x motif.

Unstranded libraries suppress orientation-dependent endpoint metrics and remove exonic antisense from the broad composition. Explicit `--library-strand` is recommended.

## 6. Novelty and relationship to existing tools

The individual ingredients are established: genomic partitioning, splice motif checks, internal-priming context, endpoint atlases, per-cell summaries, and split-alignment evidence. The useful novelty is their integration at the individual-alignment level with a common taxonomy and exact read/base accounting across droplet, long-read, and plate workflows.

scNoiseMeter complements rather than replaces:

- aligner summaries (mapping and broad genomic regions);
- ambient-RNA models such as SoupX, DecontX, and CellBender;
- cell-calling and doublet detection;
- isoform-model QC such as SQANTI3 and pigeon;
- run-level long-read QC;
- fusion callers, which apply gene-level evidence and filtering beyond an `SA` tag.

The fixed-window intergenic screen, same-read two-end anchoring, exact pre/post read retention, and separation of descriptive composition from artifact-candidate evidence are the most distinctive elements of version 0.7. They should still be externally benchmarked.

## 7. Validation status and limitations

The test suite covers category invariants, overlapping annotation, strand-aware ambiguity, exact splice boundaries and motifs, reservoir merging, site-strand loading, plate relabeling, comparison matching, and an end-to-end synthetic BAM/GTF path. CI runs supported Python versions.

Important remaining limitations are:

- no experimentally sorted positive/negative ground truth across protocols;
- no local/mappability model for intergenic enrichment;
- no gene-aware molecule deduplication;
- no read-level NUMT inference from competing alignments;
- no distinction between real fusions and technical chimeras;
- no correction/removal of reads;
- automatic reference downloads currently focus on GRCh38.

Accordingly, scNoiseMeter should be used for QC, method comparison, and hypothesis generation. Candidate artifact reads or loci should be validated with orthogonal evidence before exclusion or biological interpretation.

## 8. Reproducible reporting

A defensible methods section should state:

- scNoiseMeter version/commit;
- BAM generation and filtering stage;
- GTF, genome, repeat, polyA, and TSS resource versions;
- platform and library strandedness supplied to the CLI;
- barcode/UMI tags, whitelist, and called-cell filter;
- custom TSO sequence and whether poly-G checking was enabled;
- seed, thread count, and whether intergenic reclassification was reservoir-estimated;
- which denominator and aggregate metric was reported.

Do not relabel `broad_noncanonical_read_frac` as percent contamination. Report its component categories and the narrower artifact-candidate composition alongside it.

Every run writes `<sample>.run_info.json` next to the metrics, carrying the version, platform,
strandedness and annotation sources. The metrics table itself holds only numbers, so this file is
the record that makes an archived result set interpretable later, and it is what the cross-sample
comparison uses to decide whether two samples may be placed side by side at all.

### Comparing across samples

Two comparisons are possible and they are not interchangeable.

A **nested** comparison asks what a processing step did to a fixed set of reads, and is only
defined when both files contain the same reads. Read-level retention and category transitions
belong here. Note that this is narrower than "before and after a step": some deduplicators
rewrite read names, so a genuine pre/post pair can be unmatchable at the read level even though
both files describe the same library.

A **cohort** comparison asks how independent samples differ, and can only be distributional.
Three limits apply, and each of them is a way a cross-sample figure can invent a result.

Metric definitions change between tool versions. The same BD46 alignment file reports a TSO
invasion rate of 2.33e-03 under v0.6.1 and 3.28e-06 under v0.7.2, a 700-fold difference caused
entirely by a change in how the flag is detected. Pooling result sets from different versions
therefore produces differences that look biological and are not.

Annotation completeness changes what "intergenic" means. Samples annotated against different
GENCODE releases are not comparable in the categories that depend on gene coverage.

Strandedness changes what a category measures. In a non-stranded protocol, exonic-antisense reads
are expected cDNA signal and strand concordance near 50% is the design, so neither quantity can be
ranked against a stranded sample.

A metric that a sample never recorded is a gap in the record, not a measurement of zero, and is
reported as absent.
