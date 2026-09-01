# PICOTA Development Roadmap

Status: 2026-09-01. Companion document: [VALIDATION.md](VALIDATION.md).

---

## 1. Where we are

### Implemented

- **v1.0 pipeline** — SRA download → assembly (SPAdes / MEGAHIT) → GFA cycle
  detection → Prodigal ORF calling → BLAST scoring against CARD, ISfinder +
  TnCentral and a curated xenobiotic set → transposon/cargo boundary annotation →
  optional Oxford Nanopore long-read evidence.
- **v1.1 (unreleased)** — `output_formatter.py` and the CT-tagged
  `picota_enriched.csv`, cross-sample novel-CT clustering (`ct_cluster.py`),
  step-progress logging, E2E tests over the bundled `test_data/testNitro.gfa`.
- **Test suite** — 125 unit and smoke tests passing.

### Not implemented

- **No ground-truth benchmark.** PICOTA has never been measured against a closed
  reference genome. There is no recall or precision figure for any organism.
  See [VALIDATION.md](VALIDATION.md).
- **Deduplication deletes real composite transposons.** Quantified below.

---

## 2. Confirmed defects

Every row is reproduced by `tests/test_known_defects.py`. Tests that assert the
correct behaviour are marked `xfail(strict=True)`, so a fix turns them into
XPASS failures and forces the marker to be removed.

| ID | Location | Defect | Measured impact |
|----|----------|--------|-----------------|
| ~~**D1**~~ | [`cycle_finderv2.py`](../picota/src/cycle_finderv2.py) | ~~`parse_gfa` keeps only `Name` and `Sequence`, dropping every depth encoding~~ **Resolved in phase 1** | Node depth now parsed from SPAdes `cov_`, MEGAHIT `multi=`, `dp:f:`/`DP:f:` and `KC:i:`/`LN:i:` |
| **D2** | [`cycle_finderv2.py:348-387`](../picota/src/cycle_finderv2.py) | `cycle_match_based_on_contig_id` strips strand for its exact-duplicate test but keeps strand for its similarity test, and treats >70% shared node length as duplicate. **Fixed by `dedup_mode: strict`** | Legacy output saturates at 2 candidates regardless of ground truth (2/2, 2/3, 2/4, 2/5, 2/8). Strict recovers N/N |
| **D3** | [`cycle_kmer_hash.py:12`](../picota/src/cycle_kmer_hash.py) | `get_kmer_hashes` returns a set, not a multiset, so repeat copy number is invisible. **Fixed by `dedup_mode: strict`** | Legacy collapses `IS-cargo` and `IS-cargo-IS` into one candidate, discarding the *complete* CT. Strict keeps both |
| **D4** | [`cycle_kmer_hash.py:118-119`](../picota/src/cycle_kmer_hash.py) | Similarity denominator is `len(new_cycle)`, a containment measure used to make an identity decision. **Fixed by `dedup_mode: strict`** | Legacy result depends on DFS traversal order: the same two cycles give 1 or 2 outputs. Strict is symmetric |
| **D5** | `k_mer_sim: 80` in `config.yaml` | Exact k-mer matching is a cliff whose position depends entirely on k. **Fixed by `dedup_mode: strict`** via the Mash transform (`estimated_ani`) plus a length criterion | Legacy: shared k-mers fall 86% → 64% → 34% across 0.1% → 0.5% → 1% divergence. Strict: the same pair estimates 99.49% and 99.47% identity at k=21 and k=80, so the threshold is a sequence identity rather than an artefact of k |
| ~~**D6**~~ | [`cycle_finderv2.py`](../picota/src/cycle_finderv2.py) | ~~`path_limit` truncation assigns a dead local and reports nothing~~ **Resolved** | `GraphWork.truncated_searches` counts them and `cycle_analysis` warns. First run of the fix found **96 truncated searches on the bundled `testNitro.gfa`** at the default `path_limit: 25` — results on that graph were incomplete and nothing said so |

D2 is the defect that matters most for the whole-genome benchmark: it is
precisely the "one IS, several copies, different cargo" case that a real
*E. coli* or *K. pneumoniae* genome presents.

---

## 3. Architectural principle

**No phase changes existing output silently.**

Phase 2 is the first phase with a behavioural flag: `dedup_mode` in
`config.yaml`, defaulting to `legacy`. On the bundled `testNitro.gfa`, `strict`
returns the same five cycles as `legacy` plus one 34,953 bp candidate that
legacy was deleting.

Every behavioural change ships behind a config flag, defaulting to the current
behaviour until the benchmark shows the new path is better. Each run writes an
`algo_version` plus its resolved parameters alongside the results, so the
question "which of my previous analyses must be redone?" is answered by a diff,
not by guesswork.

---

## 4. Phases

| Phase | Work | Exit criterion | Invalidates prior analyses? |
|-------|------|----------------|------------------------------|
| **0** ✅ | Characterization tests over synthetic assembly graphs (`tests/synthetic_gfa.py`, `tests/test_known_defects.py`). No production code touched. | Every defect in §2 reproduced by a test | No |
| **0.5** | Closed-genome benchmark harness. Stage 1 ✅ `scripts/select_benchmark_strains.py` shortlists closed genomes with matching Illumina runs. Stage 2 (IS annotation → ground-truth CT catalogue) and stage 3 (run + metrics) outstanding. See [VALIDATION.md](VALIDATION.md). | A baseline recall/precision number exists | No |
| **1** ✅ | Parse node depth in `parse_gfa`; expose `depth_ratio` (max node depth / min node depth) via `Cycle.depth_ratio` and a `<cycles>.depths.tsv` sidecar. **Report-only** — nothing in detection, scoring or filtering reads it. | D1 test passes; cycle output byte-identical | No — sidecar file only |
| **2** ✅ | `src/cycle_dedup.py`: paths are duplicates only when they are the same cycle (rotation- and reverse-complement-invariant key); sequences compared by multiset Jaccard over canonical circular k-mers. Behind `dedup_mode: legacy \| strict`, default `legacy`. | D2, D3, D4 pass in strict mode | Only when `dedup_mode: strict` is set |
| **3** ✅ | `estimated_ani()` converts k-mer Jaccard to nucleotide identity (Mash transform), paired with a length-ratio criterion; threshold `dedup_min_ani` (default 99.0%). Both criteria are needed — identity alone merges a CT with the fragment inside it. | D5 tests pass at k=21, 31 and 80 | Only when `dedup_mode: strict` |
| **4** | Deterministic superbubble enumeration and IS-centric search: mark IS-like nodes first (high depth, high degree, IS BLAST hit), then collect paths through them. D6 is now *reported* rather than resolved — the 96 truncated searches on `testNitro.gfa` are exactly what this phase removes the need for. | `path_limit` removed | Yes — larger refactor |

Phases 2-4 are only meaningful once phase 0.5 exists; without a benchmark there
is no way to tell an improvement from a regression.

---

## 5. Immediate next steps

1. **Phase 0.5 stage 2** — run `scripts/select_benchmark_strains.py`, then
   annotate IS elements on the shortlisted genomes (ISEScan) and build the
   ground-truth CT catalogue described in [VALIDATION.md](VALIDATION.md) §3.2.
   Stage 1 is done; stage 2 needs a download budget and an IS annotator in the
   environment.
2. **Raise `path_limit`, or start phase 4.** The D6 counter reports 96 truncated
   searches on the bundled `testNitro.gfa` at the default `path_limit: 25`, so
   candidate cycles are being missed on a graph already used for testing. Phase 4
   removes the limit entirely by enumerating superbubbles deterministically.
3. **Flip `dedup_mode` to `strict` by default** once the phase 0.5 benchmark
   confirms on real data what the synthetic graphs already show.

Known limitation of phase 1: when depth comes from `KC:i:`/`LN:i:` the value is a
k-mer count per base, which underestimates true depth by `(length - k + 1)/length`.
The assembly k is not recorded in the GFA, so this is left uncorrected; it biases
nodes shorter than a few times k and is harmless for the rest. On the bundled
`testNitro.gfa` the resulting `DepthRatio` spans 1.32 to 7.01 across five cycles.

Phase 1 shipped a `depths.tsv` sidecar but nothing consumes it yet. Before
`depth_ratio` is allowed to influence scoring, its distribution has to be
characterised on real runs: single-copy bubbles should sit near 1 and genuine
composite transposons at or above 2, and that separation needs to be observed,
not assumed.
