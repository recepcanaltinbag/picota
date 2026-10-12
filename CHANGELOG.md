# PICOTA — Changelog & Development Notes

All notable changes are documented here. Follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [1.1.0] — 2026-10-13

Results written by 1.0.0-rc1 are not comparable with this release: defaults,
scoring and several output fields changed. See **Output changes** under
*Changed* before comparing old and new runs.

### Added
- **`picota.py` command-line interface** — five subcommands (`sra_download`,
  `assembly`, `analysis`, `scoring`, `all`), each a thin front end over the same
  stage functions `picota_testv3.py` runs, so the CLI and the config-driven
  pipeline execute the same code. `--config`, `--threads` and `--output` are
  accepted after the subcommand; `all` takes one of `--sra`, `--fastq` or `--gfa`
  (mutually exclusive) or, with none, runs the config's sample list; `--long`
  names a long-read accession for validation; `--version` prints the version.
  The tool check before any work starts is derived from the tool paths in the
  config. Previously `picota/picota.py` was an empty file in rc1, and the first
  version (4fdf66f) could not run any subcommand (see *Fixed*).
- **`picota/__init__.py`**, so `setup.py`'s `find_packages()` finds the package.
  Version is 1.1.0 in `setup.py`, `picota.py` and `picota/__init__.py`.
- **`output_formatter.py`** — new module that generates `picota_enriched.csv` alongside
  `picota_final_tab`. Each detected composite transposon receives a unique CT tag
  (CT001, CT002, …). Rows are expanded one-per-antibiotic-class. Includes
  antibiotic-class inference from 23 gene-name patterns (Beta-lactam,
  Aminoglycoside, Tetracycline, MLS, …).
- **`picota_enriched.csv`** output, columns `CT_Tag`, `Category` (Novel/Known),
  `CycleID`, `SRA_ID`, `SRA_Organism`, `CT_Length_bp`, `Score`, `Score2`, `NumIS`,
  `IS_Group`, `IS_Family`, `IS_Length_bp`, `IS_Names`, `NumAMR`,
  `Antibiotic_Class`, `Resistance_Gene`, `ARO_ID`, `ARO_Name`, `ARO_Drug_Class`,
  `NumXeno`, `Xenobiotic_Functions`, `NumCompTN`, `Known_CompTN`.
- **IS-family lookup from the ISFinder taxonomy** — `infer_is_family()` resolves
  names in three tiers: `DBs/ISes/isfinder_taxonomy.tsv` (5570 elements, built
  from the bundled `IS.fna` by `scripts/build_isfinder_taxonomy.py`), a curated
  exact-name table for elements absent from `IS.fna`, then regex rules for
  variant suffixes. Database prefixes (`ISFinder_`, `ISFDB_`, `TNcentral_`) and
  accession suffixes (`IS26-MH257753`) are stripped first; `Tn` + digits is
  reported as `Composite`.
- **`SRA_Organism`** — organism looked up from NCBI eUtils (esearch for the
  UID, then esummary).
- **`test_complete_pipeline.py` as a real pipeline runner** — replaced
  `MockAnalysis` with real pipeline calls; `--gfa_mode` runs without SRA download
  or assembly (`test_sra_ids.csv` is auto-created in that mode); `--config`
  selects a YAML, from which `outdir`, `sra_id_file` and threads are read unless
  given on the command line; every step reads its parameters from the config
  rather than hardcoded values; completed steps are skipped when their outputs
  exist (resume), and download is skipped when a GFA already exists; raw FASTQs
  are deleted after a successful assembly when `delete_fastq_files: true`;
  long-read validation downloads the long-read FASTQ when missing, then runs
  `bam_file_analyze` and `analyze_blocks` after mapping, as `picota_testv3.py`
  does.
- **`picota_summary.csv`** (`test_complete_pipeline.py`) — one table over all
  samples with a reduced column set, written by `write_summary_csv()`.
- **Cross-sample novel CT clustering** (`src/ct_cluster.py`,
  `test_complete_pipeline.py`) — novel CT sequences from all samples are compared
  all-vs-all with blastn (identity ≥90%, coverage ≥80%) and grouped into
  connected components, written to `picota_ct_clusters.tsv` with
  `CTCluster001`, … tags (largest first) and `Unique` for singletons. Runs only
  when more than one sample is processed; skipped when blastn/makeblastdb are
  absent.
- **`score3`** (`total_score_type: 3`) — a 0-100 score that gates on IS hit
  quality (`0.5 + 0.5 × quality`; no IS scores zero) and on cycle length, then
  weights a structural term against cargo, where `score2` put ninety of its
  hundred points on whether the cargo hits a database. Cargo uses the best hit's
  quality rather than the hit count, so cargo absent from CARD and the
  xenobiotic set no longer caps the score at 10. Written as its own
  `picota_final_tab` column. The default from this release (see *Changed*); `tests/test_score3.py` pins its behaviour.
- **Bounds for unattended runs** — `max_path_expansions` (options) caps the DFS
  expansions one path search may spend, and reports how many searches were
  abandoned; deterministic, so the same graph gives the same cycles on any
  machine. `spades_memory_gb` and `spades_timeout_s` (paths) pass `-m` and a wall
  clock to SPAdes, which otherwise assumes 250 GB and swaps rather than fails on
  a smaller machine. All default to `None`, the previous unbounded behaviour.
  `cycle_timeout_s` and `spades_memory_gb_retry` are accepted by the config
  loader but not read by the shipped runners.
- **SPAdes failure reason** — the last fifteen lines of SPAdes' output now reach
  the log and the raised error; previously only the exit status was logged.
- **`blast_batch`** option (default `true`) — see *Performance*.
- **Output provenance sidecars** (`src/stale_output.py`) — a BLAST result is
  reused only when `<output>.picota-sig.json` records the same inputs (hashes of
  the query and reference FASTA, program, output format). The sidecar is written
  only after the search succeeds.
- **`scripts/build_is_db.py`** — rebuilds the IS reference set from the merged
  TnCentral + ISfinder FASTA and ISfinder's `IS.fna`: repairs glued records,
  adds missing ISfinder elements, and stops on any record it cannot parse
  rather than guessing. Deterministic. `tests/test_build_is_db.py` covers each
  glued-record shape found in the distributed file. See *Output changes*.
- **Curated xenobiotic additions** — `build_xenobiotic_db.py --stage curated`
  appends `picota/src/xenobiotics_curated.fasta` (tracked in git) to the KEGG
  set, once each, with `PATH:curated` in the header. It holds one enzyme,
  QIH10984.1, the aromatic ring-hydroxylating dioxygenase alpha subunit of
  *Pseudomonas* sp. BIOMIG1BAC, which KEGG assigns no KO. 16,539 → 16,540
  sequences. `tests/test_xenobiotic_curated.py`.
- **`src/build_xenobiotic_db.py`** — rebuilds the xenobiotic database from KEGG:
  19 biodegradation pathways, keeping only KOs that BRITE br08901 places in
  "Xenobiotics biodegradation and metabolism" and in no other metabolic
  category, sequences from UniProt by EC, CARD overlap removed by homology.
  16,539 sequences against the previous 82,164, with zero mobile-element and zero
  non-enzymatic entries. `DBs/` is gitignored, so the script is shipped and the
  database must be built (see *Output changes*).
- **`docs/ROADMAP.md`** — development roadmap: table of six confirmed defects
  (D1-D6) in cycle detection and deduplication, the "no phase changes existing
  output silently" principle, and a phased plan with per-phase exit criteria.
- **`docs/VALIDATION.md`** — honest record of what has and has not been
  validated, and the design for a closed-reference-genome benchmark. Notes that
  PICOTA has never been measured against a closed genome and that no recall or
  precision figure exists for any organism.
- **`tests/synthetic_gfa.py`** — deterministic synthetic assembly graphs
  reproducing the composite-transposon topology (one repeated IS node with N
  distinct cargo branches), with LN/KC/dp depth tags.
- **`tests/test_known_defects.py`** — characterization tests for D1-D6.
  Assertions of the correct behaviour are marked `xfail(strict=True)`, so a fix
  turns them into XPASS failures. Measured today: when N distinct composite
  transposons share one IS node, the reported output saturates at two candidates
  regardless of N (2/2, 2/3, 2/4, 2/5).
- **`src/cycle_dedup.py`** (roadmap phase 2) — deduplication that cannot delete
  a distinct composite transposon. Paths are duplicates only when they are the
  same cycle: `canonical_path_key()` is invariant under rotation and
  reverse-complement traversal, so bubbles through a shared IS node all survive.
  Sequences are compared with `multiset_jaccard()` over canonical circular
  k-mers -- counting k-mers keeps copy number visible, the symmetric denominator
  removes the traversal-order dependence, and circular extraction makes two
  rotations of one cycle compare equal.
- **`dedup_mode` option** (`legacy` | `strict`, default `legacy`) in
  `config.yaml`, plumbed through `cycle_analysis()`. Measured on synthetic
  graphs where N composite transposons share one IS node: legacy reports 2
  candidates for every N (2/2, 2/3, 2/4, 2/5, 2/8), strict reports N/N. On the
  bundled `testNitro.gfa` strict returns legacy's five cycles plus one
  34,953 bp candidate legacy was deleting.
- **`scripts/simulate_ct_genome.py`** (roadmap phase 0.5) — builds a synthetic
  chromosome whose composite transposons are known exactly, using real IS
  elements from ISfinder and real cargo genes from CARD. `--shared-is` controls
  how many CTs are built on the SAME IS element with different cargo, which is
  the case PICOTA had never been tested on. `--backbone-fasta` implants into a
  real chromosome (e.g. *E. coli* K-12 MG1655) instead of random filler, at the
  cost of the host's own IS elements appearing as apparent false positives.
  Emits `ground_truth.tsv`, `ground_truth_cts.fasta` and `ground_truth.json`.
  Later options: `--novel-cts` (cargo built from random ORFs with no database
  homology; ground truth gains `Cargo_Type`), `--cargo-repeats`,
  `--cargo-is-mode different|same` (an IS inserted into the cargo), and
  `--is-copies-per-element` (free-standing copies of each element's own IS, so
  copy number varies without sharing).
- **`tests/test_e2e_simulated_ct.py`** — end-to-end test of the claim the tool
  makes: implant composite transposons, sequence, assemble with SPAdes, and
  check PICOTA recovers them, including that two CTs built on the same IS
  element come back as separate candidates. ~13 s; reads come from ART, or
  wgsim where ART is absent; skipped when neither simulator, spades.py or blastn
  is available.
- **D7, found by that test**: `min_size_of_cycle: 2000` silently drops compact
  composite transposons. The graph cycle is IS + cargo, not IS + cargo + IS, so
  an IS26-sized element with a single resistance gene produces a cycle under
  2 kb. A simulated CT with a 1877 bp cycle is missed at the default and
  recovered at 1000.
- **`scripts/run_ct_benchmark.py`** — sweeps the whole benchmark over several
  host chromosomes, CT counts, `min_size_of_cycle` values and seeds, writing one
  row per case to `summary.tsv`. Resumable: a case already in the summary is
  skipped. `--assemblers` takes SPAdes, MEGAHIT or both (reads simulated once
  per case and reused), with older summaries migrated to the new `Assembler`
  column; `--score` also runs the scoring stage and adds a `Stage`
  (detection / scored) row per case. Rows carry `METRIC_VERSION`, and
  `summarize_benchmark.py` refuses a file stamped with another one.
- **Reads are simulated with ART** (Huang et al. 2012), profile `HSXt`
  (HiSeqX TruSeq 150 bp PE). ART carries empirical per-cycle quality profiles;
  wgsim applies one uniform error rate and invents quality strings — measured
  mean Q per cycle is a flat 30 across all 150 cycles, against ART's
  36 → 41 → 37. wgsim was first kept as `--simulator wgsim` and has since been
  removed from the benchmark driver; ART is the only path.
- **`--placement random`** (now the default in `simulate_ct_genome.py`) —
  implants are scattered subject to `--spacing` rather than evenly spaced.
  Even spacing put every CT the same distance apart, a systematic artefact.
- **`scripts/score_picota_benchmark.py`** — scores reported cycles against that
  ground truth: CT recall, precision, and copy-distinctness recall (how many CTs
  sharing an IS came back as *separate* cycles). A candidate counts as a true
  positive only when it covers at least `min_coverage` of an element, the bar
  recall uses; fragments no longer count.
- **Benchmark harness and report** — `summarize_benchmark.py` (TP/FN/FP,
  sensitivity and PPV per stage and by cargo type; specificity deliberately not
  reported), `run_scenarios.py` (graded scenarios changing one structural
  property at a time; `--min-cycle-size` defaults to 1000 with a test that it
  matches `config.yaml`), `run_copy_sweep.py` / `score_copy_sweep.py` (IS copy
  number 2, 3, 4, 5, 6, 8 on any architecture, graph and contig routes for both
  assemblers; incomplete cases are skipped), `compare_contig_recovery.py`,
  `compare_tools.py` (MobileElementFinder on the true chromosome and on SPAdes
  and MEGAHIT contigs, scored by a single-alignment rule; PICOTA's collapsed
  cycles keep the union rule), `funnel_analysis.py` (which stage lost each missed
  element), `catalogue_elements.py`, `run_manifest.py` (per-case
  `run_parameters.json` with run settings and tool versions),
  `package_supplementary.py`, `score_one_case.py` (cases scored in parallel) and
  `build_report.sh`, which produces every table from one invocation.
  Scenarios use four elements per genome and one genome per seed.
- **`tests/test_simulate_ct_genome.py`** — 74 tests, the strictest being that
  `genome[CT_Start-1:CT_End]` really is the CT sequence; every recall number
  depends on it.
- **`estimated_ani()`** (roadmap phase 3) — Mash transform converting k-mer
  Jaccard into estimated nucleotide identity, so the deduplication threshold is
  a property of the sequences rather than of `k`. The same 0.5%-divergent pair
  estimates 99.49% identity at k=21 and 99.47% at k=80, where the raw Jaccard
  swings from 0.81 to 0.48. Used with a length-ratio criterion (`min_length_ratio`,
  default 0.95): both must hold, because identity alone merges a composite
  transposon with the fragment inside it — IS-cargo-IS against IS-cargo
  estimates 99.65% identity at k=80. Thresholds configurable as `dedup_min_ani`
  (default 99.0) and `dedup_min_jaccard` (default 0.85), with a dedup-specific
  `DEDUP_KMER_SIZE = 21`.
- **Deduplication is deliberately conservative**, and docs/ROADMAP.md §3 records
  why: no k-mer statistic separates "one cycle with assembly noise" from "two
  CTs sharing an IS with different cargo". Measured across k, the gap is at best
  +0.06 and inverts from k=31 — pairs that must stay separate score higher than
  pairs that should merge. Strict mode therefore keeps both when in doubt;
  reporting a near-duplicate is recoverable, deleting a CT is not.
- **`GraphWork.truncated_searches`** — `path_limit` truncation is counted and
  `cycle_analysis()` warns that the enumeration was not exhaustive. Previously
  the truncation assigned a dead local and reported nothing, so a truncated
  search was indistinguishable from a complete one. On the bundled
  `test_data/testNitro.gfa` the default `path_limit: 25` truncates 96 searches;
  raising it to 50 truncates none and returns the same 35 paths, so nothing was
  lost there — but establishing that required a second run, which is the point
  of the warning.
- **`tests/test_cycle_dedup.py`** — 61 tests pinning rotation and
  reverse-complement invariance, multiset counting, Jaccard symmetry, and that
  genuine duplicates are still removed, plus the identity estimate's stability
  across k and the two-criterion rule.
- **Node depth parsing** (roadmap phase 1) — `parse_segment_depth()` in
  `cycle_finderv2.py` reads coverage from all four encodings PICOTA can meet:
  the SPAdes `_cov_` field in node names, MEGAHIT `multi=`, the GFA1 `dp:f:`/
  `DP:f:` tag, and `KC:i:` divided by `LN:i:`. `parse_gfa()` now records it as
  `Depth` on both strands of each node; unknown depth is `None`, never `0`.
- **`Cycle.depth_ratio`** — most-covered node depth over least-covered, which
  approximates the copy number of the repeated IS node collapsed into the bubble.
- **`Cycle.node_spans`, `depth_at()` and `junction_depth_ratio()`** — depth over
  the annotated IS divided by depth over the cargo. Not a copy-number estimate:
  a genome with 16 copies of the shared element gives ratios of 1.8 to 3.8. It
  carries a direction only (medians of roughly 2.3 for multi-copy against 1.6 for
  two-copy IS across the scenarios) and feeds no scoring decision.
- **`<cycles>.depths.tsv` sidecar** — per-cycle depth report, including node
  spans, written next to the cycle FASTA by `cycle_analysis()`. **Report-only**:
  no detection, scoring or filtering decision reads it, and cycle output is
  unchanged.
- **`scripts/select_benchmark_strains.py`** (roadmap phase 0.5, stage 1) —
  shortlists benchmark strains reproducibly: closed RefSeq genomes that also
  have public Illumina paired-end WGS runs from the same BioSample, with the
  read depth each run would give over that genome. Network access is isolated
  in `EntrezClient` so the parsing and ranking logic is testable offline.
  Deliberately not restricted to NCBI reference genomes: that filter reduces
  E. coli from 5350 closed assemblies to 2, and the survivor is K-12 MG1655,
  a lab strain with essentially no resistance composite transposons.
  `--reference-only` restores the restriction. `--resolve-strains` fills the
  strain name from BioSample (the only record that carries it); the `Instrument`
  column reports the sequencer model, which implies read length and therefore
  how well the assembly resolves repeats.
- **`tests/test_select_benchmark_strains.py`** — 45 offline tests using record
  shapes taken from live NCBI esummary responses.
- **`tests/test_depth.py`** — 34 tests covering every depth encoding, the
  `depth_ratio` edge cases, the sidecar report, and a guard that depth tags do
  not change which cycles are reported.
- **Step-progress logging** in `picota_testv3.py` — five clearly labelled steps with
  elapsed time, step banners, and a per-sample summary block printed to stdout/log.
- **`test_e2e_gfa.py`** — real end-to-end integration tests using bundled
  `testNitro.gfa`. Covers cycle detection, BLAST scoring, enriched CSV validation,
  and unit tests for `output_formatter.py` (31 tests, no mock data).
- **`tests/test_cli.py`** — 21 tests. Stage functions are patched with
  `autospec`, so a CLI call that does not bind to the real signature fails, plus
  an end-to-end `analysis` run on `testNitro.gfa`.
- **Tests for modules that had none** — `tests/test_split_cycle_coords.py` and
  `tests/test_analyze_blocks.py`, each verified against the pre-fix module.

### Changed
#### Output changes
An existing run re-executed with this release can write different results:
- **`total_score_type` defaults to 3** (`score3`), was 0, in `config.yaml` and in
  `scoring_main()`'s own default. `threshold_final_score` now applies to `score3`;
  all four scores are still written. `score0` is unbounded and sums hit scores,
  so the same threshold meant something different in every run.
- **`picota_final_tab` gains a `score3` column** after `score2`.
- **`min_size_of_cycle` now defaults to 1000**, was 2000. The graph cycle of a
  composite transposon is IS + cargo rather than IS + cargo + IS, so a compact
  element yields a cycle far shorter than itself. Twelve implanted elements of
  2,470–3,095 bp produce cycles of 1,761–2,280 bp: the old threshold recovered
  8 of 12, the new one all 12, at a cost on a wild-type genome of 17 candidate
  cycles rising to 22. IS26 with a single resistance gene sits in this range.
- **`assembler_type` now defaults to `spades`** in all shipped configs, was
  `megahit`. On simulated genomes with known composite transposons, PICOTA
  recovers 40/40 elements from the SPAdes assembly graph and 19/40 from the
  MEGAHIT graph of the same reads. MEGAHIT simplifies repeats harder, so the
  bubbles are not in its graph to be found; its perfect precision reflects
  reporting less, not reporting better.
- **One insertion sequence per IS locus, not per cycle** (52ec265). `NumIS`
  counts IS loci; it could previously only be 0 or 1, and was 1 for every
  reported candidate. `ISproducts`/`IScoords` list every non-overlapping locus,
  best-scoring first; overlapping hits still compete for one locus. `score0`,
  which sums IS scores, changes; `score3` uses the best hit and does not. In a
  12,877-candidate survey 7,036 carry two or more IS at distinct coordinates, and
  the candidates "carrying IS26" at score3 > 50 grow from 751 to 1,180.
  Protein databases keep the single best hit.
- **Tied IS references are all named** (18e597b). IS26, IS26.var1, IS15DII and
  IS15DIV score identically at many loci; keeping one made the name depend on
  BLAST row order (22 of 1,158 loci named differently between two batchings).
  Every reference within `TIE_EPS` of the best score at a locus is reported,
  joined with `|` and sorted, e.g.
  `ISFinder_IS15DII|ISFinder_IS15DIV|ISFinder_IS26|ISFinder_IS26.var1`.
- **`IS_Group`/`IS_Family` for tied names** in `picota_enriched.csv`: each name
  in a `|`-joined tie is resolved separately and disagreeing answers are joined
  with `|`, sorted. The tie above was reported as `IS15`/`IS15` and is now `IS6`.
  Affects every locus where references tie.
- **Known-transposon match requires ≥95% identity over ≥80% of both sequences**
  (772180c). The old criterion, reference coverage × identity > 80, accepted a
  candidate that merely contained a known element (a 30 kb cycle covering the
  1,243 bp Tn5466_partial record was reported as Tn5466). Candidates with a
  `Known_CompTN` name fall from 522 to 48 in a 13,390-candidate survey and from
  1,696 to 243 in a 68-genome benchmark; there, name recall moves 68 → 52 of 89,
  false names 112 → 46, PPV 0.378 → 0.531. Changes `Category` (Known/Novel).
- **Transposon/cargo split uses one IS locus** (3d0ca5c). `split_cycles_from_picota`
  bounds the transposon block by the best-scoring locus (`IScoords[0]`) instead
  of the first start to the last end of all loci, which after 52ec265 swallowed
  the cargo. Output is identical wherever a cycle has a single IS.
- **Long-read gap acceptance** (b98f021) uses the configured floor per pair;
  some gaps previously accepted are now rejected, so `analyze_blocks` pattern
  counts are not directly comparable. They can be regenerated from existing
  `_full`/`_partial` files without re-mapping.
- **Xenobiotic database and its names** (d5fc75a, d8866b5). The database built by
  `build_xenobiotic_db.py` replaces the keyword-matched set, in which 464
  mobile-element proteins were 0.6% of the file but 29.3% of its xenobiotic
  annotations and a further 28.3% were antibiotic resistance genes already in
  CARD. Its headers spell description and organism with underscores, and
  `Xenoproducts` reads `description|KO:…|EC:…` with the KO as the gene; the older
  classified names pass through unchanged. Xenobiotic hits, and the scores built
  from them, differ from rc1. The database was built in September but the
  shipped configs kept pointing at `Xenobiotics_classified.fasta`; all three
  (`config.yaml`, `config_IS26.yaml`, `config_random.yaml`) now name
  `Xenobiotics_kegg.fasta`, as do `test_complete_pipeline.py`, the benchmark
  drivers and the tests. On the bundled `testNitro.gfa` this takes the report
  from two candidates to one. `Cycle_3`, whose only xenobiotic hit was a
  recombinase-family protein, falls below the threshold. `Cycle_2` lost its
  aromatic ring-hydroxylating dioxygenase hit, which has no KO, and dropped
  from score3 85.5 to 55.7. That enzyme is now in the set as a curated
  addition (below), which brings `Cycle_2` back to 84.4; the remaining
  difference is its MFS transporter, a transporter rather than a degradation
  enzyme, which stays out.
- **IS reference set** — `path_to_ises` now names `picota_IS_db.fasta`, built
  by the new `scripts/build_is_db.py`. In the merged
  `_tncentral_nointegrall_isfinder-TNs.fasta` 21 records had lost the line
  breaks around their header, so BLAST searched each as the tail of the record
  before it, under that record's name; they are split back out. 495 ISfinder
  elements present in the bundled `IS.fna` and absent from the merged set, by
  name and by sequence on either strand, are added as `ISFinder_<name>`
  (5,409 → 5,904 records). Of the elements added, 298 had no BLAST hit at all
  in the old set, so a cycle carrying one could not pass the IS gate. MITEs
  are not added (`--include-mites` adds them): they encode no transposase and
  cannot flank a composite transposon, but counted as IS hits they lift a
  cycle through the IS gate. On the simulated benchmark, with them in, a
  non-transposon cycle went from score3 17 to 29 on `MITEEc1` alone. The merged
  set's own 2 MITEs stay. The set covers all eight insertion sequences the
  benchmark generator draws that the old set lacked: seven by name, and
  IS200G through `ISFinder_IS1541A`, which has the same sequence. On
  `testNitro.gfa` the report is unchanged. `config_IS26.yaml` and
  `config_random.yaml` named `tncentral_integrall_isfinder-TNs.fasta`, which
  does not exist; they now name `picota_IS_db_integrall.fasta`, the same repair
  applied to the INTEGRALL-inclusive `tncentral_integrall_isfinder.fa`
  (6,094 → 6,589 records).
- **IS families in `picota_enriched.csv`** come from the ISFinder taxonomy, and
  numbered ISCR and ISEcp members (ISCR4 to ISCR30, ISEcp2, ISEcp3), previously
  `Unknown`, are assigned.
- **Existing GFAs are recomputed**, not reused: `assembly_main()` removes and
  rebuilds a GFA it finds, since it may describe older reads or parameters.

#### Other changes
- **`config.yaml` no longer carries one machine's paths.** `outdir` and
  `log_file` pointed at `/media/lin-bio/back2/picota_project_RandomSRA`, and
  `sra_id_file` at that project's sample sheet. They are now
  `picota_results`, `picota_results/picota.log` and `picota/samples.csv`.
  `config_random.yaml` keeps the RandomSRA project's settings.
- **BLAST threads default to 8**, was 4; `PICOTA_BLAST_THREADS` overrides it, and
  `picota.py --threads` sets it. Hits are identical at every thread count.
- **First run after upgrading repeats BLAST searches**: outputs without a
  `.picota-sig.json` sidecar are recomputed.
- `assembly_main()` returns the path of the GFA it selected; Bandage scores are
  cached in `gfa_scores.json` and reused only when the cached set matches the
  input set.
- `run_sra_down()` takes `expected_files` (default 2); pass 1 for single-end runs.
- `scoring_main()` calls `write_enriched_csv()` at the end, writing
  `picota_enriched.csv` next to `picota_final_tab` (effective from 126423a, see
  *Fixed*). `test_complete_pipeline.py` always rewrites it from an existing
  `picota_final_tab`.
- README updated with enriched output column reference and E2E test instructions.

### Fixed
- **Every `picota.py` subcommand failed at run time** (4fdf66f): the CLI imported
  `run_assembly_pipeline`, `run_scoring_pipeline` and `analyze_gfa_file`, which do
  not exist, and read a nonexistent `args.sra_file`.
- **Assembly deleted the reads it was given.** With `skip_filtering: true`
  (the shipped setting) SPAdes and MEGAHIT read the caller's FASTQs directly,
  and cleanup removed them unless `keep_temp_files` was set: right for reads the
  pipeline downloaded, destructive for reads a user pointed it at.
  `assembly_main(delete_inputs=)` now decides that, default true, so the SRA
  pipeline behaves as before. `picota assembly` and `picota all --fastq` pass
  false. Reads fastp wrote are intermediates and are still removed.
- **Single-end assembly passed SPAdes `-1` without `-2`**, which SPAdes rejects.
  A single file now goes in as `-s`. Single-end downloads were fixed in d58663e,
  but their assembly still failed.
- **`picota` subcommands, found in pre-release review:**
  - **Missing database:** a run with a missing reference database finished with
    exit 0 and an empty table, because scoring only warns. The CLI now stops
    first and names the missing path.
  - **Relative database paths:** these were resolved against the working
    directory only, so every one was "missing" when `picota` ran outside the
    repository root. They now also resolve against the repository root.
  - **Log file:** logging went to the config's `log_file`, a path on the
    config author's machine, and crashed when its directory did not exist. A CLI
    run now logs to `<output>/picota.log`.
  - **Samples sharing one `--output`:** two samples shared one `scoring/`
    directory, and every stage skips work whose output exists. A second sample
    therefore got the first one's `picota_final_tab`, split against its own
    cycles. Each sample now has `scoring/<sample>/` and `annot/<sample>/`, and
    `all --fastq` assembles into `assembly/<sample>/`.
  - **Sample name:** cycle files were named after the GFA. SPAdes names its best
    graph by k (`99.gfa`), so `SRAID` came out as `99`. Cycle files are now
    `<sample>_<graph>.fasta`, as `picota_testv3.py` names them. `--name`
    overrides the sample name; `_` and `.`, which would truncate it in
    `SRAID`, become `-`.
  - **`--threads`** set `PICOTA_BLAST_THREADS` after the scoring module had
    already read it at import, so BLAST kept 8 threads. It now sets the thread
    count directly.
  - **Nested output:** `all --fastq` failed on the nested output directory,
    because `assembly_main` creates only one level.
  - **`all` without input:** with no input and no sample sheet, `all` raised a
    traceback; it now prints an error naming `sra_id_file`.
- **Tied IS names with an unresolvable member**: `IS26|<unknown>` reported
  `IS6|Unknown` and a bare `|` reported an empty group. Unresolvable names are
  now ignored when any other name resolves, and the result is `Unknown` when
  none does.
- **`pysam` and `matplotlib` were missing from `requirements.txt`, `setup.py`
  and the README's conda command.** `picota_testv3.py` imports both at start-up,
  for long-read validation, so a fresh install failed on the first command even
  without long reads.
- **`picota_enriched.csv` was never written** by a normal run: `scoring_main`
  imported `output_formatter` as a top-level module, which fails when scoring is
  imported as `src.scoringv4ProtBlast`, and a broad `except` logged it as a
  warning. A missing formatter now fails loudly; a write error is still
  tolerated.
- **Config keys with no effect**: `max_path_expansions` now reaches
  `cycle_analysis()` and `spades_memory_gb`/`spades_timeout_s` reach
  `assembly_main()`, in both `picota_testv3.py` and `test_complete_pipeline.py`.
  Defaults remain `None`, so default output is unchanged.
- **`test_complete_pipeline.py` ignored parts of the config**: it read the
  prodigal/blast tool paths from `paths` while `config.yaml` keeps them under
  `options` (custom tool paths were replaced by bare names), hardcoded the four
  database paths instead of `path_to_*`, and never passed `dedup_mode`,
  `dedup_min_ani` or `blast_batch`. With the shipped config the result is
  unchanged.
- **`score3` wiring and design** — the dispatch list stopped at `score2`, so
  `score3` raised `IndexError`, and the benchmark scripts' `--score-type`
  rejected 3 (af53f46). Shared-IS
  candidates scored 49.3 to 49.9 against a threshold of 50; component tolerance
  now scales with multi-copy evidence and the IS gate no longer halves a weak
  hit (7a7943a). Large wild-type cycles carrying one IS and one CARD hit about
  30 kb apart scored 68.6 and 65.8; length now gates the score and they score
  31.1 and 26.6 (e5c1776). The depth term never ran (every candidate took the
  unknown-depth default, capping scores at 85) and no monotonic function of the
  depth ratio separates real from spurious cycles, so it was removed and the
  remaining weights renormalised (4d0b41c).
- **A truncated or stale BLAST output was accepted as a result**, including after
  the database was rebuilt (4b60831).
- **Single-end SRA runs failed their own completeness check** (d58663e): two
  FASTQ files were required, so every ONT and PacBio download was rejected as
  incomplete after it finished.
- **`analyze_blocks` gap floor** (b98f021): `UnboundLocalError` when no earlier
  pair in a read was near the mode (10 of 98 samples in one run), and otherwise a
  value leaked from an earlier gap (up to 5000 where the floor is 1000). The
  module logger is no longer `None` before `analyze_blocks()` runs.
- **Multi-k SPAdes assembly** could not work: input FASTQs were deleted after the
  first k. Deleting inputs is now the last k's job.
- **`fastq-dump` failures** — a four-step fallback (parallel-fastq-dump,
  fasterq-dump on the local `.sra`, fasterq-dump by accession, prefetch +
  fasterq-dump), and fasterq-dump temp files go to the output path rather than
  `/tmp`.
- **`test_complete_pipeline.py` runner errors**: `assembly_main` called with
  keyword arguments it does not take (`TypeError`); long-read `run_sra_down` call
  missing required arguments; `split_min_score` hardcoded to 0; long-read FASTQ
  re-downloaded when present; a finished sample skipped wholesale, including
  long-read validation; `ModuleNotFoundError` when run outside `picota/`;
  `--sra_list` default broken from the script's own directory; bare tool names
  (`megahit`) and `./picota/tools/...` paths resolved to nonexistent files
  (exit 127).
- **IS-family inference** — prefix rules matched IS10 as IS1 and IS21/IS26 as
  IS3; replaced by the taxonomy lookup. `SRA_Organism` lookup passed the
  accession to esummary, which needs a numeric UID.
- **Benchmark ground truth** — IS copy number omitted the copy `--cargo-is-mode
  same` inserts into the cargo; the interrupting IS in `different` mode could
  equal a flanking one; `run_scenarios.py` hardcoded `min_size_of_cycle` 2000 and
  did not copy the `.depths.tsv` sidecar beside `cycles.fasta`; half-finished
  cases entered the copy-sweep aggregate. Benchmark precision fell from 96/97
  (99.0%) to 92/97 (94.8%) once fragments stopped counting (3ac646b).
- `output_formatter._TAB_COLS` lists `score3`, matching the header
  `picota_final_tab` writes.

### Performance
- **One BLAST search per database, not one per cycle** (`blast_batch: true`).
  Twenty cycles against the four bundled databases: 23.5 s over 80 invocations
  against 13.1 s over 4. `picota_final_tab` is byte-identical, checked on 4 and
  40 cycles; the transposon parser's best hit is now resolved per cycle rather
  than per result file.
- **Strict deduplication**: k-mer reverse complements computed once per window
  with `str.translate`, `multiset_jaccard` walks only the smaller multiset, and
  the length criterion is checked first. A graph that hit the 45-minute timeout
  (ERR3053723) finishes in 1.8 minutes; results identical.
- **Cycle duplicate check** uses a set of sequences instead of comparing each new
  cycle with every accepted one, and reverse paths are deduplicated before their
  sequences are assembled. `testNitro.gfa` gives a byte-identical cycles FASTA.
- **Resumed runs skip BLAST** when the sidecar matches: 40 cycles took 16.8 s on
  a first run and 5.8 s on a second.
- **Bandage scores cached** in `gfa_scores.json` (see *Changed*).

---

## [1.0.0-rc1] — 2026-04-02

### Core pipeline (`picota_testv3.py`)
End-to-end workflow orchestrated per sample:
1. SRA download (skip if FASTQ already present)
2. Genome assembly (MEGAHIT or SPAdes → GFA)
3. Cycle detection on assembly graph
4. BLAST scoring (CARD + IS-finder/TnCentral + KEGG)
5. Annotation split (transposon/cargo boundaries)
6. Long-read validation (minimap2 + BAM analysis)

### Bug fixes applied during development
| Location | Issue | Fix |
|----------|-------|-----|
| `cycle_finderv2.py` | `NameError: name 'cycle'` in f-strings (lines 417, 428) | Changed `{cycle}` → `{path}` |
| `cycle_finderv2.py` | `overlap_length` potentially unbound | Added `overlap_length = 0` default |
| `cycle_finderv2.py` | Duplicate cycles from reverse-strand paths | Added strand-agnostic `frozenset(n[:-1])` comparison in `cycle_match_based_on_contig_id` |
| `scoringv4ProtBlast.py` | Bare `except:` | Changed to `except (IndexError, AttributeError):` |
| `scoringv4ProtBlast.py` | `split('|')` index out of range | Added `len(parts) >= N` guards |
| `scoringv4ProtBlast.py` | `IndexError` at line 611 on files without `_` in name | Safe `_basename_parts[1] if len > 1` guard |
| `assembly.py` | Shell injection risk in Bandage call | Added `shlex.quote()` |
| `long_read_checkIS26.py` | `KeyError` on missing FASTQ read | Added `if qname not in fastq_dict: continue` |

### Performance optimisations
| Module | Optimisation | Impact |
|--------|-------------|--------|
| `cycle_kmer_hash.py` | Replaced O(n²) pairwise k-mer comparison with inverted index | Linear scaling for large cycle sets |
| `cycle_finderv2.py` | Replaced `+=` string concatenation with list+join | Avoids O(n²) memory copies for long sequences |
| `cycle_finderv2.py` | `_searched_pairs` set in `dfs_iterative` to skip redundant path expansions | Removes duplicate cycle edges |
| `scoringv4ProtBlast.py` | Pre-computed `fwd_key`/`rev_key` in `parse_gfa` | Halves `reverse_sign()` calls |

### Test suite
- `test_cycle_finder.py` — 56 unit tests for graph cycle detection
- `test_scoring.py` — 37 unit tests for scoring functions
- `test_kmer_hash.py` — 22 unit tests for k-mer filtering
- `test_integration_gfa.py` — 19 integration tests using `testNitro.gfa`
- `test_smoke.py` — 45 smoke tests across 8 pipeline steps

---

## Algorithm Overview

```
Raw reads (FASTQ)
      │
      ▼
[1] fastp             Quality filtering (optional)
      │
      ▼
[2] MEGAHIT / SPAdes  Genome assembly  ──►  assembly.gfa  (GFA format)
      │
      ▼
[3] cycle_finderv2    Parse GFA → directed graph (contig+ / contig- nodes)
                      Iterative DFS → back-edges = cycle candidates
                      find_paths() → enumerate simple paths per cycle edge
                      Sequence reconstruction (overlap-aware)
                      Deduplication:
                        (a) strand-agnostic frozenset(node[:-1])
                        (b) weighted path overlap
                        (c) k-mer Jaccard (inverted index)
                      ──►  cycles.fasta
      │
      ▼
[4] Prodigal (meta)   ORF prediction on each cycle sequence
                      ──►  .faa (proteins), .fasta (nucleotides), .gbk
      │
      ▼
[5] BLAST             blastp vs CARD (antibiotics)
                      blastp vs KEGG (xenobiotics)
                      blastn vs IS-finder / TnCentral (insertion sequences)
                      blastn vs CompTn DB (known composite transposons)
      │
      ▼
[6] scoring_main      Z-score normalised composite score:
                        z = |len - μ_ct| / σ_ct   (μ=5850, σ=2586 bp)
                        score = (Σ_ant + Σ_is + Σ_xeno) ^ (1 - z/max_z)
                      Four scoring modes (score0/1/2/3; default score3)
                      Threshold filter → ranked hit list
                      GenBank annotation
                      ──►  picota_final_tab, picota_enriched.csv
      │
      ▼
[7] split_cycle_coords  Annotate transposon / cargo boundaries
                        ──►  annot/<accession>/*.fasta
      │
      ▼
[8] minimap2 + samtools  Map ONT long reads to candidate CT FASTA (optional)
    bam_analyse          Coverage analysis, circular-read detection
    analyze_blocksv3     Block-pattern scoring (TCT, CTCT, …)
                         ──►  mapping/<accession>/
```

### Key data structures

```python
class CodingRegion:
    start, end, strand     # genomic coordinates
    fullname               # raw BLAST sseqid
    r_type                 # 'Antibiotics' | 'InsertionSequences' | 'Xenobiotics' | 'CompTNs'
    score                  # (match_len / slen) * pident
    product, gene          # parsed from fullname

class GeneticInfo:
    seq_id, seq_acc, seq_description
    feature_list: List[CodingRegion]
    nuc_seq: str
    score0, score1, score2, score3: float
```

---

## Known limitations / future work

- IS-family annotation relies on name prefixes; full taxonomy lookup against ISFinder
  API would be more accurate.
- Antibiotic-class inference uses regex on gene names; integrating CARD ARO ontology
  would cover edge cases.
- Organism field in enriched CSV is currently the SRA accession; NCBI metadata lookup
  would populate species names.
- `find_all_path=True` can be exponential for highly connected graph regions;
  `max_path_expansions` bounds it, but is off by default.
- `samtools` output requires `libncurses.so.5` on some Linux systems — install via
  `apt install libncurses5` or use the conda-forge build.
