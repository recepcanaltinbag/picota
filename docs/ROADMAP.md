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
- **Read depth is discarded.** Assemblers report per-node coverage and PICOTA
  throws it away, so the one signal that distinguishes a multi-copy IS from a
  single-copy region is unused.
- **Deduplication deletes real composite transposons.** Quantified below.

---

## 2. Confirmed defects

Every row is reproduced by `tests/test_known_defects.py`. Tests that assert the
correct behaviour are marked `xfail(strict=True)`, so a fix turns them into
XPASS failures and forces the marker to be removed.

| ID | Location | Defect | Measured impact |
|----|----------|--------|-----------------|
| **D1** | [`cycle_finderv2.py:238-244`](../picota/src/cycle_finderv2.py) | `parse_gfa` keeps only `Name` and `Sequence`; `dp:f:`, `KC:i:`, `LN:i:` and the SPAdes `cov_` field in node names are all dropped | Copy number of an IS node is unknowable downstream |
| **D2** | [`cycle_finderv2.py:348-387`](../picota/src/cycle_finderv2.py) | `cycle_match_based_on_contig_id` strips strand for its exact-duplicate test but keeps strand for its similarity test, and treats >70% shared node length as duplicate | **Output saturates at 2 candidates** regardless of ground truth: recall 2/2, 2/3, 2/4, 2/5 for N distinct CTs sharing one IS |
| **D3** | [`cycle_kmer_hash.py:12`](../picota/src/cycle_kmer_hash.py) | `get_kmer_hashes` returns a set, not a multiset, so repeat copy number is invisible | `IS-cargo` and `IS-cargo-IS` collapse to one candidate — and it is the *complete* CT that gets discarded |
| **D4** | [`cycle_kmer_hash.py:118-119`](../picota/src/cycle_kmer_hash.py) | Similarity denominator is `len(new_cycle)`, a containment measure used to make an identity decision | Result depends on DFS traversal order: same two cycles give 1 or 2 outputs depending on input order |
| **D5** | `k_mer_sim: 80` in `config.yaml` | Exact 80-mer matching is a cliff, not a gradient | Shared k-mers fall 86% → 64% → 34% across 0.1% → 0.5% → 1% divergence, crossing the 80% threshold with no usable middle ground. IS copies within one genome sit exactly in this range |
| **D6** | [`cycle_finderv2.py:157-158`](../picota/src/cycle_finderv2.py) | `path_limit` truncation assigns a dead local (`LIMIT = True`) and reports nothing | Number of paths lost to truncation is unknown and unlogged |

D2 is the defect that matters most for the whole-genome benchmark: it is
precisely the "one IS, several copies, different cargo" case that a real
*E. coli* or *K. pneumoniae* genome presents.

---

## 3. Architectural principle

**No phase changes existing output silently.**

Every behavioural change ships behind a config flag, defaulting to the current
behaviour until the benchmark shows the new path is better. Each run writes an
`algo_version` plus its resolved parameters alongside the results, so the
question "which of my previous analyses must be redone?" is answered by a diff,
not by guesswork.

---

## 4. Phases

| Phase | Work | Exit criterion | Invalidates prior analyses? |
|-------|------|----------------|------------------------------|
| **0** | Characterization tests over synthetic assembly graphs (`tests/synthetic_gfa.py`, `tests/test_known_defects.py`). No production code touched. | Every defect in §2 reproduced by a test | No |
| **0.5** | Closed-genome benchmark harness: ground-truth CT catalogue, PICOTA run, metrics. See [VALIDATION.md](VALIDATION.md). | A baseline recall/precision number exists | No |
| **1** | Parse node depth in `parse_gfa`; expose `depth_ratio` (IS node / cargo node) as a **report-only** column. | D1 test passes; benchmark numbers unchanged | No — new column only |
| **2** | Replace sequence-similarity deduplication with node-multiset identity; keep containment as a recorded parent/child relation instead of a deletion. | D2, D3, D4 tests pass; benchmark recall improves | Yes, once the flag is enabled |
| **3** | Replace raw k-mer identity with IS-annotation-driven deduplication (two candidates sharing an IS family but carrying different cargo are two CTs, never duplicates). | D5 test passes; benchmark recall improves | Yes |
| **4** | Deterministic superbubble enumeration and IS-centric search: mark IS-like nodes first (high depth, high degree, IS BLAST hit), then collect paths through them. | D6 resolved; `path_limit` removed | Yes — larger refactor |

Phases 2-4 are only meaningful once phase 0.5 exists; without a benchmark there
is no way to tell an improvement from a regression.

---

## 5. Immediate next steps

1. Phase 0.5 — select benchmark strains (criteria in [VALIDATION.md](VALIDATION.md)),
   build the ground-truth CT catalogue, record baseline metrics.
2. Phase 1 — depth parsing, report-only.
3. Phase 2 — deduplication rewrite, measured against the phase 0.5 baseline.
