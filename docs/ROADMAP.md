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
- **Test suite** — unit, smoke, and an end-to-end test
  (`tests/test_e2e_simulated_ct.py`) that implants composite transposons,
  sequences them with wgsim, assembles with SPAdes and checks PICOTA recovers
  them. About 13 seconds; skipped when the tools are absent.

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
| **D5** | `k_mer_sim: 80` in `config.yaml` | Exact k-mer matching is a cliff whose position depends entirely on k. **Mitigated, not fixed.** `estimated_ani()` gives a k-stable identity estimate, but no k-mer statistic can carry the merge decision — see §3.1 | Legacy: shared k-mers fall 86% → 64% → 34% across 0.1% → 0.5% → 1% divergence. Strict adds a Jaccard floor and a length criterion and errs toward keeping candidates |
| ~~**D7**~~ | `min_size_of_cycle` in `config.yaml` — **resolved, default now 1000** | The graph cycle for a composite transposon is IS + cargo, not IS + cargo + IS, so the default threshold silently drops compact CTs. IS26 at 820 bp plus one resistance gene yields a cycle under 2 kb | Found by `tests/test_e2e_simulated_ct.py`: an implanted CT whose cycle is 1877 bp is missed at the default and recovered at `min_size_of_cycle: 1000`. Compact IS26-bounded single-gene units are among the most common CTs in clinical isolates, so the default is not neutral |
| ~~**D6**~~ | [`cycle_finderv2.py`](../picota/src/cycle_finderv2.py) | ~~`path_limit` truncation assigns a dead local and reports nothing~~ **Resolved** | `GraphWork.truncated_searches` counts them and `cycle_analysis` warns that the enumeration was not exhaustive. On the bundled `testNitro.gfa` the default `path_limit: 25` truncates 96 searches; raising it to 50 truncates none and returns the **same 35 paths** (longest 14 nodes), so nothing was actually lost there. The counter flags *unverified completeness*, not proven loss |

D2 is the defect that matters most for the whole-genome benchmark: it is
precisely the "one IS, several copies, different cargo" case that a real
*E. coli* or *K. pneumoniae* genome presents.

---

## 3. A limitation worth stating plainly

**k-mer statistics cannot tell scattered differences from clustered ones.**

Deduplication has to separate two situations: one cycle assembled twice with a
little noise (should merge) and two composite transposons sharing an IS element
but carrying different cargo (must not merge). Measured on simulated assemblies,
multiset Jaccard for those two cases is:

| k | same cycle, 0.5% noise | different cargo (real pairs) | gap |
|---|---|---|---|
| 11 | 0.896 | 0.432 – 0.835 | +0.061 |
| 21 | 0.814 | 0.406 – 0.787 | +0.027 |
| 31 | 0.740 | 0.385 – 0.748 | **−0.008** |
| 80 | 0.484 | 0.318 – 0.634 | **−0.150** |

There is no safe threshold at any k, and from k=31 the ordering inverts: pairs
that must stay separate score *higher* than pairs that should merge. This is not
a tuning problem. A k-mer set records how many k-mers differ, never whether the
differences are spread through the sequence or sit in one block, and that is
exactly the distinction required.

The Mash transform does not rescue it. It assumes k-mer loss comes from
scattered substitutions, so it reads a 20% block substitution as about 1% of
uniform divergence: cycles differing only in cargo estimate 99.0–99.7% identity
while sharing 52–79% of their k-mers. Merging on identity alone deletes exactly
the candidates the pipeline exists to find.

So strict deduplication is deliberately conservative — high Jaccard floor, and
it keeps both candidates when in doubt. Reporting a near-duplicate is visible
and can be collapsed downstream (`ct_cluster.py` already does cross-sample
clustering); deleting a composite transposon is invisible and unrecoverable.
Doing better needs alignment, not k-mers, and the pipeline already runs BLAST
later — folding that back into deduplication is future work.

---

## 4. Architectural principle

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

## 5. Phases

| Phase | Work | Exit criterion | Invalidates prior analyses? |
|-------|------|----------------|------------------------------|
| **0** ✅ | Characterization tests over synthetic assembly graphs (`tests/synthetic_gfa.py`, `tests/test_known_defects.py`). No production code touched. | Every defect in §2 reproduced by a test | No |
| **0.5** | Closed-genome benchmark harness. Stage 1 ✅ `scripts/select_benchmark_strains.py` shortlists closed genomes with matching Illumina runs. Stage 2 (IS annotation → ground-truth CT catalogue) and stage 3 (run + metrics) outstanding. See [VALIDATION.md](VALIDATION.md). | A baseline recall/precision number exists | No |
| **1** ✅ | Parse node depth in `parse_gfa`; expose `depth_ratio` (max node depth / min node depth) via `Cycle.depth_ratio` and a `<cycles>.depths.tsv` sidecar. **Report-only** — nothing in detection, scoring or filtering reads it. | D1 test passes; cycle output byte-identical | No — sidecar file only |
| **2** ✅ | `src/cycle_dedup.py`: paths are duplicates only when they are the same cycle (rotation- and reverse-complement-invariant key); sequences compared by multiset Jaccard over canonical circular k-mers. Behind `dedup_mode: legacy \| strict`, default `legacy`. | D2, D3, D4 pass in strict mode | Only when `dedup_mode: strict` is set |
| **3** ◐ | `estimated_ani()` plus a Jaccard floor (`dedup_min_jaccard`, default 0.85), a length-ratio criterion and a smaller dedup-specific k (`DEDUP_KMER_SIZE = 21`). All three criteria are load-bearing. D5 is mitigated rather than fixed — see §3. | Benchmark precision improves; recall unchanged | Only when `dedup_mode: strict` |
| **4** | Deterministic superbubble enumeration and IS-centric search: mark IS-like nodes first (high depth, high degree, IS BLAST hit), then collect paths through them. D6 is now *reported* rather than resolved — this phase removes the need to check truncation by hand at all. | `path_limit` removed | Yes — larger refactor |

Phases 2-4 are only meaningful once phase 0.5 exists; without a benchmark there
is no way to tell an improvement from a regression.

---

## 6. Assembler choice

`assembler_type` now defaults to `spades` in all shipped configs; it was
`megahit`. Measured on simulated genomes with known composite transposons, same
reads and same ground truth:

| | recall | PPV |
|---|---|---|
| SPAdes graph | **40/40 (100%)** | 46/48 (96%) |
| MEGAHIT graph | **19/40 (48%)** | 19/19 (100%) |

MEGAHIT's perfect precision is not an advantage — it reports fewer candidates
because its graph holds fewer bubbles, having simplified the repeats harder. The
trade is one-directional for this task: a structure missing from the graph
cannot be recovered by any scoring downstream, whereas a spurious candidate can
be filtered. The same effect appears in the contigs: elements intact in a single
contig were 3/8 under SPAdes and 0/8 under MEGAHIT.

---

## 7. Scoring modes

`total_score_type` selects among four. They are not interchangeable and the
choice changes what the tool can find at all.

| mode | form | novel cargo | range | ranks? |
|---|---|---|---|---|
| 0 | sum of homology hits, unbounded | yes | 65–300 | weakly |
| 1 | presence x 100 per class | yes | — | no |
| 2 | `90 x [cargo] x [IS] + 10^z` | **never** | 96–100 | no |
| 3 | `gate x (0.40 structure + 0.30 multi-copy + 0.30 cargo)` | yes | 50–100 | yes |

score2 cannot report cargo absent from CARD and the xenobiotic set at any
threshold: with no such hit its first term is zero, leaving a maximum of 10
against a threshold of 50. A composite transposon carrying a metabolic operon is
still one, so defining the structure by its cargo answers a different question
from the one the tool claims to. The same imbalance costs it ranking — ninety of
its hundred points are one yes/no question, so a cycle of eighteen components
scores within four points of a clean two-component one.

score3 inverts that balance and is the mode this benchmark recommends. Two
design points were established by measurement rather than assumed:

**Component count and multi-copy evidence are not independent.** A first draft
penalised many components uniformly, which put every shared-IS candidate at
49.3–49.9 against a threshold of 50 — reporting none of them, on the case the
tool exists for. A two-copy element should assemble into a clean bubble, so many
components there is real evidence against it; an element whose IS sits in dozens
of copies necessarily threads many nodes, and charging it the same penalty bills
the candidate for the structure being detected. Component tolerance now scales
with the multi-copy evidence.

**The depth term is load-bearing, not decorative.** Removing it was measured
against keeping it: a genuine shared-IS element scores 72.4 against 54.0 for an
unexplained complex cycle with the term, and 62.0 against 78.0 without — the
ordering inverts, because component count is then penalised uniformly and the
better IS hit wins. This supersedes the earlier finding in §3 that per-cycle
depth added nothing: that was true of `depth_ratio` (max node over min node),
which tracks cycle length, and false of `junction_depth_ratio`, which compares
the annotated IS against the cargo.

---

## 8. Negative control

Wild-type *E. coli* K-12 MG1655, no implanted elements, same reads, assembly and
pipeline: 17 cycles in the graph, **0 reported under score0 and 2 under score2**.

The graph of an untouched genome does contain cycles — MG1655 carries dozens of
native IS elements — and scoring rejects essentially all of them. The two
surviving under score2 are most likely its native *ampC*-family genes, which
CARD recognises, flanked by IS. This is not a false-positive rate against zero
truth, since those may be real IS-bounded structures; what it establishes is
that the pipeline does not manufacture composite transposons where none were
placed.

---

## 9. Minimum cycle size

`min_size_of_cycle` now defaults to 1000; it was 2000. The graph cycle of a
composite transposon is IS + cargo, not IS + cargo + IS, because both flanking
copies collapse into a single node — so a compact element produces a cycle far
shorter than itself and the old threshold discarded it before anything scored it.

| threshold | compact elements recovered | cycles on a wild-type genome |
|---|---|---|
| 2000 | 8/12 | 17 |
| 1500 | 12/12 | 21 |
| **1000** | **12/12** | **22** |
| 500 | 12/12 | 23 |

Twelve implanted elements of 2,470–3,095 bp produce cycles of 1,761–2,280 bp.
Recovery is complete from 1500 down, and the cost measured on wild-type MG1655
with nothing implanted rises only from 17 candidate cycles to 22. 1000 was
chosen over 1500 for headroom: IS26 at 820 bp with a single resistance gene —
the most common composite transposon in clinical isolates — yields a cycle near
1,400 bp and would fall under 1500.

This is a threshold moved to a measured place, not removed. An element more
compact still would fall under 1000 too. Replacing the length cut with a
structural test — does the cycle contain a complete IS and at least one ORF —
would remove the class of failure rather than relocate it, and belongs with the
scoring work rather than here.

---

## 10. Immediate next steps

1. **Phase 0.5 stage 2** — run `scripts/select_benchmark_strains.py`, then
   annotate IS elements on the shortlisted genomes (ISEScan) and build the
   ground-truth CT catalogue described in [VALIDATION.md](VALIDATION.md) §3.2.
   Stage 1 is done; stage 2 needs a download budget and an IS annotator in the
   environment.
2. **Raise the default `path_limit`, or start phase 4.** The D6 counter reports
   96 truncated searches on the bundled `testNitro.gfa` at `path_limit: 25`;
   raising it to 50 truncates none and returns the same 35 paths, so on that
   graph nothing was lost — but that could only be established by re-running.
   On a real genome the same warning would need the same check, every time.
   Phase 4 removes the limit entirely by enumerating superbubbles
   deterministically, which removes the question rather than the symptom.
3. **Decide what to do about D7.** Lowering `min_size_of_cycle` recovers compact
   composite transposons but admits more short spurious cycles; the benchmark
   can measure that trade directly rather than leaving the default to habit.
4. **Flip `dedup_mode` to `strict` by default** once the benchmark sweep
   confirms on more genomes what the current cases already show.

Known limitation of phase 1: when depth comes from `KC:i:`/`LN:i:` the value is a
k-mer count per base, which underestimates true depth by `(length - k + 1)/length`.
The assembly k is not recorded in the GFA, so this is left uncorrected; it biases
nodes shorter than a few times k and is harmless for the rest. On the bundled
`testNitro.gfa` the resulting `DepthRatio` spans 1.32 to 7.01 across five cycles.

**Superseded:** `junction_depth_ratio` — the annotated IS against the cargo,
rather than max node over min node — is now a weighted term in score3 and
carries the discrimination described in §7. The finding below stands only for
the cruder `depth_ratio`.

Phase 1 shipped a `depths.tsv` sidecar; score3 consumes it. Before
`depth_ratio` is allowed to influence scoring, its distribution has to be
characterised on real runs: single-copy bubbles should sit near 1 and genuine
composite transposons at or above 2, and that separation needs to be observed,
not assumed.
