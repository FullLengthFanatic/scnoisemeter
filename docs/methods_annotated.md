# scNoiseMeter annotated methods

**Code-linked implementation companion for version 0.9.0**

This document maps the scientific rules to modules and functions. It avoids fixed line numbers because they become stale during refactoring. See [whitepaper.md](whitepaper.md) for rationale and literature context, and [documentation.md](documentation.md) for the CLI/output reference.

## 1. Taxonomy and aggregate sets

Source: `scnoisemeter/constants.py`.

`ReadCategory` defines the output taxonomy. `CATEGORY_ORDER` is the canonical output/report order and intentionally excludes the three record states that are not genomic classifications: unmapped, secondary, and supplementary.

`BROAD_NONCANONICAL_CATEGORIES` is a descriptive composition, not a causal noise estimator. `ARTIFACT_CANDIDATE_CATEGORIES` is the smaller hotspot/chimeric set. The older `NOISE_*` constants remain as internal aliases of these sets; the corresponding `noise_*` output fields were removed in 0.8.0, so nothing in the metrics tables or MultiQC JSON carries the old vocabulary. `modules/cohort.py` maps the old names forward when reading result directories written before that.

`ARTIFACT_CANDIDATE_CATEGORIES` contains no orientation-dependent category, so unlike the broad set it needs no stranded/unstranded variant. The benchmark evaluator imports both sets from `constants.py` rather than restating them.

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

`_is_multimapper()` uses only `NH > 1`. MAPQ alone does not establish how many alignments were reported: it is a mapping-confidence score, and aligner conventions differ. `XA` is not interpreted generically because the SAM specification reserves tags beginning with `X` for local use, so its meaning is producer-specific. Low-MAPQ records retain their genomic category and are reported through the independent low-MAPQ metric. When an aligner does not emit `NH`, scNoiseMeter does not infer primary-read multimapping from unrelated fields. `minimap2` and `pbmm2` fall in that group: they report secondary alignments through the `tp` tag and score secondaries through `s1`/`s2`, and write no `NH`. On those BAMs the multimapper category is empty by construction, and the affected alignments are counted under their genomic category rather than held out, which shifts the classification denominator relative to a STAR-derived BAM.

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

`priming_criterion` selects the test: `a_run` requires `POLYA_RUN_MIN_LENGTH` consecutive matching bases anywhere in `POLYA_CONTEXT_WINDOW`; `perc_a` requires `POLYA_PERC_A_MIN` percent over the first `POLYA_PERC_A_WINDOW` bases, the SQANTI3/pigeon convention on a 0–100 scale; `either` accepts one or both. The fraction is counted on `T` for a reverse alignment rather than reverse-complementing the slice. `ReadClassifier.__init__` rejects any other value, so the option cannot be accepted and silently ignored. Threaded CLI → `run_pipeline` → task dict → `_contig_worker` → `ReadClassifier`.

### Untemplated 3′ polyA tail

`_polya_tail_stats()` returns `(soft_clip_len, terminal run, A-fraction)` for the transcript-3′ clip: the trailing clip counted on `A` for a forward alignment, the leading clip counted on `T` for a reverse one, reversed so index 0 is always the base adjacent to the cleavage site. Only CIGAR op 4 qualifies; op 5 returns `(0, 0, 0.0)` because hard-clipped bases are absent from `query_sequence` and indexing the aligned sequence would read alignment bases as a tail.

`_has_polya_tail()` requires `run >= POLYA_TAIL_MIN_LENGTH`. `POLYA_TAIL_MIN_FRAC` defaults to 0 and is therefore not part of the decision: a genuine 12 bp tail inside a 60 bp clip carrying adapter or TSO sequence scores about 0.2, so requiring both would implicitly demand a ~45 bp tail.

### RT template-switch direct repeat

`_direct_repeat_len(contig, donor, acceptor)` takes the longest common suffix of the `RT_SWITCH_REPEAT_SCAN` bases preceding each intron boundary, where `donor` is the first intronic base and `acceptor` one past the last — the CIGAR-N interval convention used by `_get_junctions()`. `N` terminates the run. `_has_rt_switch_repeat()` returns True when any junction reaches `RT_SWITCH_REPEAT_MIN`, and False without a reference FASTA. Independent of `_has_noncanonical_junction()` by design.

### Alignment quality

`ReadResult` carries MAPQ, optional `NM`, soft-clipped bases, query bases, and aligned reference bases. `_contig_worker()` accumulates these per cell. `metrics.compute_metrics()` reports mean MAPQ, mean edit distance among reads carrying `NM`, the softclip/query-base fraction, and low-MAPQ rate.

## 7. Intergenic profiler

Source: `intergenic_profiler.py`.

### Input and sampling

The pipeline stores a uniform reservoir of up to 500,000 intergenic records. Each record includes contig, genomic span, strand, barcode, junction flag, strand-correct 3′ coordinate, aligned bases, UMI, read length, and stable read key.

`_merge_reservoir()` combines contig or well reservoirs using the number of records seen, not equal weighting by contig. This prevents small contigs from being overrepresented.

### Windows and null model

`_cluster_reads()` assigns each record to `(contig, three_prime // 500)`. Strand is intentionally excluded from the key so that `strand_fraction` is measured rather than predetermined.

`compute_intergenic_bases(by_contig=True, contig_lengths=...)` computes the exact complement of the union of gene bodies through each non-mitochondrial BAM contig end. The annotation-gap object is only a compatibility fallback when lengths are unavailable. CLI multiple-testing correction counts the windows actually tested: the intergenic base total divided by `INTERGENIC_LOCUS_WINDOW`, derived from the same figure that feeds the background rate. Empty intergenic windows are still counted, which keeps the correction conservative. Before 0.9 the denominator was every window of every non-mitochondrial contig, including all genic sequence where no test is performed — 6,176,550 windows on GRCh38 against 2.16M–2.78M tested, a 2.2–2.9× overcorrection whose effect ran one way, suppressing real enrichment. `_score_locus()` uses `scipy.stats.poisson.sf(n_reads - 1, expected)` and Bonferroni adjustment.

`_assign_local_read_share()` sets `IntergenicLocus.local_read_share` after every window is scored, since it needs the whole window set. Per contig it sorts by start and builds a prefix sum, so each window's ±`LOCAL_SHARE_FLANK` neighbourhood total is a subtraction rather than a scan — one sort per contig instead of a neighbourhood walk per window.

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

## 7a. Site-atlas loading and pre-write guards

`_load_polya_sites(paths, chrom_style, min_support, pas_out)` parses column 7 as the percentage of samples supporting a record and drops those below `min_support`; a record without a parseable column 7 is kept, so a plain BED3 is unaffected. `min_support` is part of the pickle cache key, so a cached unfiltered dict is never served to a filtered run. When `pas_out` is supplied it is filled with the records whose column 11 names a hexamer, harmonised and pruned to the same keys as the main dict.

`_polya_atlas_density()` merges ±`POLYA_SITE_PROXIMITY` windows per contig and returns the covered fraction of the non-mitochondrial reference; `_warn_if_atlas_too_dense()` records a warning above `POLYA_ATLAS_DENSITY_WARN`. The density is a property of the file, so it is measured once per run rather than assumed.

`_check_site_coordinate_range()` runs after `_check_site_contig_overlap()` at all four load sites. Positions at or beyond their contig length are counted; above `SITE_COORD_OUT_OF_RANGE_DROP` the atlas is returned as `None` so its metrics take the absent path, above `SITE_COORD_OUT_OF_RANGE_WARN` it warns. Contigs absent from the BAM are skipped rather than counted against the atlas.

`_check_barcode_gate()` compares `n_reads_unassigned / n_reads_classified` against `UNASSIGNED_GATE_WARN`, only when `meta.barcode_aware is True`. Counters that are not numbers are skipped via `_as_number()` rather than raised on: the guard runs before every report and must not be what fails a valid run.

`modules/selfcheck.py` runs the result invariants and writes `<sample>.qc.json`. `run_invariants()` returns passing and failing checks alike; `failures()` filters by severity; the whole block is wrapped at the call site so a defect in the checking code cannot fail an otherwise valid analysis. Severity is `ERROR` except where a documented configuration legitimately breaks the property — reservoir-estimated reclassification for the conservation checks, and the anchoring veto for the PAS-subset check.

## 8. End anchoring

For exonic-sense reads, `_contig_worker()` stores the strand-correct 5′ and 3′ coordinates together with one read key in `exonic_sense_endpoints`. It uses `reference_end`, not `reference_start + query_alignment_length`, so spliced alignments are not displaced by intron length.

`compute_metrics()` calls the strand-aware `_near_polya_site()` for both atlases:

- `three_prime_anchored_frac` when polyA sites exist;
- `three_prime_at_pas_site_frac` against `result._polya_pas_dict`, the subset of atlas records whose column 11 names a hexamer;
- `five_prime_anchored_frac` when TSS sites exist;
- `both_ends_anchored_frac` only when both exist and both match the same read.

`result._anchor_veto_internal_priming` makes the 3′ hit additionally require that the read was not flagged for internal priming. The endpoint reservoir record gained a sixth field in 0.9 for that per-read flag; `compute_metrics()` indexes defensively so five-field records written by 0.8 still load.

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

The first three of those are conditional. `cli._probe_nested_overlap()` opens both BAMs before the
pipeline runs, takes up to `COMPARE_PROBE_SAMPLE_SIZE` primary read names from the first shared
contig of B, and measures how many occur in A over the same genomic window. Below
`COMPARE_NESTED_MIN_OVERLAP` the pair is treated as independent: `store_read_assignments` stays
False, the matching outputs are skipped, and the run warns. Iteration is over A testing membership
in B's name set, so probe memory is bounded by the sample size rather than by A's density.

This is a memory control as much as a correctness one. `read_assignments` costs roughly 175 bytes
per read, measured with real PacBio names, so retaining it for both sides of a 76M-read and a
37M-read pair needs about 20 GB simultaneously. `--matched-reads` and `--no-matched-reads`
override the probe.

`modules/cohort.py` covers the other comparison shape: N independent samples read from finished
result directories rather than BAMs. It shares no code with the paired path because it can share
no assumptions, and it refuses to silently pool samples from different tool versions, different
GENCODE releases, or mixed strandedness.

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
- `tests/test_cohort.py`: reading legacy and current result directories, absent-versus-zero
  handling, the version / annotation / strandedness comparability warnings, ranking fallback when
  a metric predates some samples, sample-sheet labelling and ordering, and that the deviation
  heatmap leaves a near-zero row neutral instead of painting it as the strongest signal.
- `tests/test_priming_evidence.py`: the untemplated-tail measurement including hard clips,
  reverse-strand orientation and the adapter-bearing clip; atlas support filtering, hexamer
  extraction and BED3 passthrough; atlas density and its warning; and that the option defaults
  preserve pre-0.9 behaviour.
- `tests/test_robustness_transfers.py`: both guards in both directions -- firing on the failure
  they exist for, silent on the legitimate lookalike (healthy sample, barcode-free BAM, in-range
  atlas, contigs absent from the BAM) -- plus direct-repeat measurement and local read share.
- `tests/test_selfcheck.py`: one corruption per invariant, each asserted to fire by name, and a
  healthy result set that trips none of them.
- `.github/workflows/tests.yml`: Python 3.10-3.13 CI with Ruff (pinned) and pytest.

`tests/benchmark/` builds a synthetic BAM with known per-read truth and scores the classifier
against it. It covers 16 of the 17 classified categories; `INTERGENIC_ENRICHED` has no ground
truth because constructing a window that clears the enrichment screen while failing every
positive-evidence rule is fragile to threshold changes. The evaluator does score the column, so
reads landing there are counted rather than dropped from the confusion matrix, and it imports the
aggregate sets from `constants.py` so they cannot drift from the classifier again.

## 13. Safe tuning guidance

| Goal | Preferred action |
|---|---|
| Correct protocol orientation | Set `--library-strand` explicitly |
| Correct sequencer/library context | Set `--platform`; do not infer from aligner alone |
| Use another TSO | Supply the exact `--tso`; repeat for multiple motifs |
| Require exact TSO evidence | Use `--no-polyg-tso` and increase `--tso-min-match` |
| Enable context/junction checks | Supply `--reference` |
| Interpret repeats | Supply a matching RepeatMasker BED |
| Reproduce reservoirs | Seeded by default (42); keep thread count fixed |
| Make polyA-site proximity discriminating again | Set `--polya-min-support 25` (or 50). Unfiltered, PolyASite 3.0 puts 38.6% of the genome within the proximity window; the run reports the density it measured |
| Compare internal priming with isoform-QC numbers | Use `--priming-criterion perc_a`, the SQANTI3/pigeon convention. It disagrees with the default on ~3% of windows, so switch for comparability, not correction |
| Stop proximity outvoting contrary evidence at the 3′ end | Use `--anchor-veto-internal-priming`; it lowers `three_prime_anchored_frac` |
| Detect RT template-switch artifacts | Supply `--reference`; the direct-repeat flag needs it and is silent without it |
| Continue past a failed pre-write check | `--no-strict`. Prefer fixing the cause: each of the three checks fires on a condition that otherwise produces a complete but wrong result set |
| Alter intergenic thresholds | Change named constants, document them, and revalidate |

Every option in the first block above changes what a reported number means, which is why each is recorded in `run_info.json` and defaults to pre-0.9 behaviour. Changing one mid-study splits the cohort.

Do not tune `--chimeric-distance`; it is deprecated and has no scientific decision role in version 0.7.
