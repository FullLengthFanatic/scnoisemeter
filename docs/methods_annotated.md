# scNoiseMeter annotated methods

**Code-linked implementation companion for version 0.7.1**

This document maps the scientific rules to modules and functions. It avoids fixed line numbers because they become stale during refactoring. See [whitepaper.md](whitepaper.md) for rationale and literature context, and [documentation.md](documentation.md) for the CLI/output reference.

## 1. Taxonomy and aggregate sets

Source: `scnoisemeter/constants.py`.

`ReadCategory` defines the output taxonomy. `CATEGORY_ORDER` is the canonical output/report order and intentionally excludes the three record states that are not genomic classifications: unmapped, secondary, and supplementary.

`BROAD_NONCANONICAL_CATEGORIES` is a descriptive composition, not a causal noise estimator. `ARTIFACT_CANDIDATE_CATEGORIES` is the smaller hotspot/chimeric set. The older `NOISE_*` constants alias these sets for compatibility.

Important constants:

| Constant | Default | Role |
|---|---:|---|
| `ILLUMINA_CHIMERIC_INSERT_SIZE` | 1,000,000 bp | Extreme paired-template fallback |
| `LOW_MAPQ_THRESHOLD` | 10 | Descriptive flag only; never filters |
| `TSO_MIN_MATCH_LENGTH` | 12 bp | Motif prefix used at soft-clipped end |
| `TSO_POLYG_MIN_LENGTH` | 6 | 10x-specific poly-G shortcut |
| `INTERGENIC_LOCUS_WINDOW` | 500 bp | Fixed 3′-end window width |
| `ADAPTIVE_MIN_READS` | 5 | Minimum reads for enrichment promotion |
| `ADAPTIVE_PVALUE_THRESHOLD` | 0.01 | Bonferroni-adjusted threshold |
| `ADAPTIVE_MIN_BARCODES_ABSOLUTE` | 3 | Candidate-novel minimum barcode support |
| `INTERGENIC_REPEAT_MIN_FRACTION` | 0.50 | Repeat-overlap fraction |
| `POLYA_RUN_MIN_LENGTH` | 6 | Internal-priming context run |
| `POLYA_CONTEXT_WINDOW` | 20 bp | Genomic sequence context window |
| `POLYA_SITE_PROXIMITY` | 50 bp | 3′ atlas proximity |
| `TSS_SITE_PROXIMITY` | 100 bp | 5′ atlas proximity |

`DEFAULT_CHIMERIC_DISTANCE` and `--chimeric-distance` remain accepted for API/CLI compatibility but do not determine same-strand split classification.

## 2. Annotation construction

Source: `scnoisemeter/modules/annotation.py`.

### Parsing and interval sets

`AnnotationIndex._build()` parses genes, transcripts, and exons into PyRanges structures. Exons and inferred introns are separated by strand. Mitochondrial contigs are not inserted into nuclear annotation sets because mitochondrial mapping is classified earlier.

`_unique_and_shared()` computes same-strand shared gene intervals and splits them into coding–coding and coding–non-coding sets. Opposite-strand overlap remains resolvable as sense versus antisense and is not inserted into a shared bucket.

The shared sets preserve strand. A global unstranded merge would erase the property needed by `ReadClassifier._get_contig_intervals()`.

### Exact junctions

`_extract_splice_junctions()` groups exons by transcript and records exact `(intron_start, intron_end, strand)` tuples. The older independent splice-site set is retained as a compatibility fallback for old/mocked indexes.

### Cache identity

The annotation cache key includes file path, size, nanosecond mtime, and a hash of the first 64 KiB for both GTF and repeats input. This reduces stale-cache risk after an in-place edit that preserves a filename.

## 3. Record filters and denominators

Sources: `classifier.py`, `pipeline.py`, and `metrics.py`.

`ReadClassifier.classify()` returns `None` for unmapped, secondary, or supplementary records. `_contig_worker()` still counts record flags before classification, producing:

- `n_primary_mapped`;
- `n_secondary`;
- `n_supplementary`;
- `n_qcfail`;
- `n_duplicate`.

`run_pipeline()` reads `.bai` statistics for mapped/unmapped totals and processes every reference with mapped records. No 1 Mb contig-size exclusion is used, so mitochondrial, viral, spike-in, decoy, and small-scaffold alignments are not silently dropped.

`SampleMetrics.n_records_total` is mapped plus unmapped from the index. Category fractions use classified mapped-primary alignments. `unmapped_read_frac` uses unmapped/total. Historical `n_reads_*` and `read_frac_*` labels remain compatibility aliases; paired mates are separate alignment records.

## 4. Barcode handling

`ReadClassifier._get_tags()` reads corrected tags, default `CB` and `UB`.

- With a whitelist, missing or off-list CB becomes `unassigned`.
- Without a whitelist, missing CB is pooled under `NO_BARCODE`, allowing sample-level BAM QC.
- `--cell-barcodes` is different: the pipeline drops uncalled cells before aggregation.

`pipeline.relabel_barcode()` updates read/base/UMI/flag/alignment structures plus intergenic records and read assignments. `run-plate` uses it to convert each barcode-free well’s sentinel into the well identity before plate merging.

## 5. Priority classification

The early priority in `ReadClassifier.classify()` is:

1. barcode gate;
2. multimapper;
3. mitochondrial;
4. chimeric;
5. genomic interval plurality.

Independent artifact/context checks run for every classified primary alignment, including early categories.

### Multimapper

`_is_multimapper()` uses `NH > 1` when `NH` exists. If `NH` is absent, a non-empty `XA` tag is accepted as explicit alternative-hit evidence. MAPQ alone does not establish how many alignments were reported: it is a mapping-confidence score, and aligner conventions differ. Low-MAPQ records retain their genomic category and are reported through the independent low-MAPQ metric.

### Chimeric alignment

`_check_chimeric()` parses each `SA` entry:

- different contig → candidate chimera;
- different strand → candidate foldback/chimera;
- same contig/strand → `_query_interval_from_cigar()` checks whether query and genomic segment order are compatible.

Distance alone is ignored. `_check_chimeric_paired()` adds inter-contig mates and extreme `TLEN` when paired detection is enabled. `_is_discordant_pair()` separately records missing/improper mates.

### Genomic base partition

`_classify_by_intervals()` obtains aligned blocks from `get_blocks()`. `_partition_block()` creates an event sweep across every interval boundary. At each atomic segment one category wins in this order:

1. shared coding–coding;
2. shared coding–non-coding;
3. sense exon;
4. antisense exon;
5. sense intron;
6. antisense intron;
7. intergenic fallback.

The invariant is:

```text
sum(base_counts.values()) == sum(end - start for start, end in read.get_blocks())
```

If this invariant is ever exceeded, classification raises rather than hiding double-counting. The read label is a deterministic plurality vote over the exact base partition.

### Intronic refinements

`_get_junctions()` advances reference position through CIGAR operations and returns both boundaries for every `N`.

- Intronic bases plus any junction become `intronic_jxnspan`.
- Sense-exonic plus intronic bases without `N` convert the intronic component to `intronic_boundary`.
- Remaining intronic bases are `intronic_pure`.

`_check_junction_canonicality()` fetches donor and acceptor and reverses both correctly for minus-strand transcript orientation. It accepts GT–AG, GC–AG, and AT–AC. Exact annotated transcript junctions bypass motif penalties.

## 6. Sequence and alignment flags

### TSO invasion

`_check_tso_invasion()` examines only the molecularly appropriate terminal soft clip: left for a forward alignment and right for a reverse alignment. `_approx_prefix_match()` tolerates one substitution in the configured prefix. This avoids testing both read ends against a motif that should be directional.

The CLI resolves protocol-aware motifs with `_default_tso_sequences()` and `_effective_tso_sequences()`. Supplying `--tso` replaces the default. The classifier API treats `None` as its historical two-motif default and an explicit empty list as disabled.

### TSO concatemer

`_check_tso_concatemer()` screens the query for more than one approximate non-overlapping occurrence of a configured full motif or its reverse complement. The tolerance is about 10% substitutions. This is a sequence heuristic; it does not model indels.

### Internal polyA priming

`_check_polya_priming()` uses reference sequence after the 3′ end of `+` alignments and before the 3′ end of `-` alignments, checking A and T runs respectively. It returns no flag when reference context cannot be fetched.

### Alignment quality

`ReadResult` carries MAPQ, optional `NM`, soft-clipped bases, query bases, and aligned reference bases. `_contig_worker()` accumulates these per cell. `metrics.compute_metrics()` reports mean MAPQ, mean edit distance among reads carrying `NM`, the softclip/query-base fraction, and low-MAPQ rate.

## 7. Intergenic profiler

Source: `intergenic_profiler.py`.

### Input and sampling

The pipeline stores a uniform reservoir of up to 500,000 intergenic records. Each record includes contig, genomic span, strand, barcode, junction flag, strand-correct 3′ coordinate, aligned bases, UMI, read length, and stable read key.

`_merge_reservoir()` combines contig or well reservoirs using the number of records seen, not equal weighting by contig. This prevents small contigs from being overrepresented.

### Windows and null model

`_cluster_reads()` assigns each record to `(contig, three_prime // 500)`. Strand is intentionally excluded from the key so that `strand_fraction` is measured rather than predetermined.

`compute_intergenic_bases(by_contig=True, contig_lengths=...)` computes the exact complement of the union of gene bodies through each non-mitochondrial BAM contig end. The annotation-gap object is only a compatibility fallback when lengths are unavailable. CLI multiple-testing correction uses every fixed non-mitochondrial BAM window, including gene-only and empty windows, as a conservative pre-defined upper bound. `_score_locus()` uses `scipy.stats.poisson.sf(n_reads - 1, expected)` and Bonferroni adjustment.

The minimum-read and adjusted-p thresholds determine enrichment. Barcode support is required for `intergenic_novel`, not for technical hotspots, allowing barcode-free sample-level BAMs to reveal internal-priming candidates.

### Evidence hierarchy

1. At least 50% repeat-span overlap → `intergenic_repeat`.
2. Not enriched/supported → `intergenic_sparse`.
3. Monoexonic + strand-aware polyA run + not near matched atlas site → `intergenic_hotspot`.
4. Dominant strand fraction ≥0.80 + minimum real barcodes + splice or polyA-site evidence → `intergenic_novel`.
5. Otherwise → `intergenic_enriched`.

The repeat label is an annotation class and does not require Poisson enrichment.

### Reclassification consistency

`cli._apply_intergenic_reclassification()` moves read and aligned-base counts, length-bin counts, read assignments, and exact UMI sets when the reservoir is exhaustive. For sampled datasets it estimates category counts from the uniform reservoir and marks affected UMI sequence-diversity values missing.

## 8. End anchoring

For exonic-sense reads, `_contig_worker()` stores the strand-correct 5′ and 3′ coordinates together with one read key in `exonic_sense_endpoints`. It uses `reference_end`, not `reference_start + query_alignment_length`, so spliced alignments are not displaced by intron length.

`compute_metrics()` calls the strand-aware `_near_polya_site()` for both atlases:

- `three_prime_anchored_frac` when polyA sites exist;
- `five_prime_anchored_frac` when TSS sites exist;
- `both_ends_anchored_frac` only when both exist and both match the same read.

Unstranded libraries suppress these orientation-dependent values. `full_length_read_frac` is a deprecated alias of the both-end result. The old platform-specific length fallback is not used.

BED loaders preserve `+`, `-`, or unknown `.` in dictionary keys. An unknown-strand site is available to either orientation; a known opposite-strand site is not.

## 9. UMI sequence diversity

The pipeline stores sets of raw corrected UMI strings by cell and category. Metrics calculate:

```text
unique UMI strings / category reads
```

The result is named `umi_sequence_diversity_*`. It is not gene-aware molecule complexity and does not error-correct UMIs. `umi_complexity_*` is an output compatibility alias only. `--no-umi-tracking` disables the sets; `--no-umi-dedup` is retained as a legacy alias.

## 10. Comparisons

`pipeline._read_key()` combines query name with mate suffix `/1`, `/2`, or `/0`. When `store_read_assignments=True`, `compare` retains key → `(category, barcode)`.

`cli._write_compare_outputs()` writes:

- exact key-matching summary;
- retention by category from A;
- category transition counts among shared keys;
- sample composition deltas;
- median paired-cell category-fraction delta and a deterministic 1,000-resample percentile bootstrap interval.

It deliberately emits no read-level independent-samples p-value.

`cli._write_category_tagged_bam()` is the optional audit trail for `run`. It copies all records in coordinate order and writes local tag `sn:Z:<category>` on classified mapped-primary alignments. Intergenic tags use final fixed-window calls rather than sampled aggregate assignments; the output BAM is indexed after writing.

## 11. Platform inference and strandedness

Source: `utils/bam_inspector.py`.

Explicit platform hints have priority. Minimap2 alone maps to `unknown`; bare STAR maps to generic `illumina`; Cell Ranger/STARsolo and BD-specific pipeline hints refine the library family. Warnings ask for explicit input when the header is insufficient.

`cli._is_unstranded_library()` uses an explicit CLI choice first. In auto mode only `smartseq` is assumed unstranded. Generic/unknown data are treated as stranded for reporting but receive a warning.

## 12. Testing map

- `tests/test_core.py`: taxonomy, TSO, chimera, exact junctions, canonical motifs, interval invariants, metrics helpers.
- `tests/test_intergenic.py`: fixed windows, evidence rules, reclassification, strand-aware annotation and site loaders.
- `tests/test_integration.py`: synthetic GTF/BAM through annotation, inspection, pipeline, and metrics.
- `tests/test_adversarial.py`: empty/tiny input and validation failures.
- `tests/test_cell_barcodes.py`: called-cell filtering.
- `tests/test_discover.py` and `tests/test_sample_sheet.py`: platform discovery and plate metadata.
- `.github/workflows/tests.yml`: Python 3.9 and 3.12 CI with Ruff and pytest.

## 13. Safe tuning guidance

| Goal | Preferred action |
|---|---|
| Correct protocol orientation | Set `--library-strand` explicitly |
| Correct sequencer/library context | Set `--platform`; do not infer from aligner alone |
| Use another TSO | Supply the exact `--tso`; repeat for multiple motifs |
| Require exact TSO evidence | Use `--no-polyg-tso` and increase `--tso-min-match` |
| Enable context/junction checks | Supply `--reference` |
| Interpret repeats | Supply a matching RepeatMasker BED |
| Reproduce reservoirs | Supply `--seed` and keep thread count fixed |
| Alter intergenic thresholds | Change named constants, document them, and revalidate |

Do not tune `--chimeric-distance`; it is deprecated and has no scientific decision role in version 0.7.
