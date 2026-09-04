<p align="center">
  <img src="logo/picota_logo.png" alt="PICOTA" width="180">
</p>

<h1 align="center">PICOTA</h1>

<p align="center">
  <b>Pipeline for Identification of Composite Transposons from Assembly graphs</b><br>
  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/python-3.8%2B-blue" alt="Python 3.8+"></a>
  <a href="./LICENSE"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="MIT"></a>
  <a href="https://doi.org/10.5281/zenodo.21769753"><img src="https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21769753-blue" alt="DOI"></a>
</p>

---

A composite transposon is cargo — antibiotic resistance genes, a degradation
pathway — flanked by two copies of the same insertion sequence. Short-read
assemblers collapse those two copies into one contig, so the element is rarely
present as a contiguous sequence to be found. What survives is its shape: the
collapsed IS and its cargo close a **cycle in the assembly graph**.

PICOTA enumerates those cycles and scores them, so composite transposons can be
recovered from draft assemblies that never resolved them into contigs.

```
     [ IS ] ──── [ cargo ] ──── [ IS ]          in the genome
                    ↓
     both IS copies collapse into one node       in the assembly graph
                    ↓
              ┌── IS ──┐
              │        │                          a cycle: IS + cargo,
              └─ cargo ┘                          shorter than the element
```

Measured on 480 implanted elements across 120 simulated genomes: **95.0%
sensitivity at 90.9% precision across five of six scenarios**, and on ten
wild-type genomes with nothing implanted, none of the 180 candidate cycles clear
the threshold. The sixth scenario — where the flanking IS also occurs inside the
cargo — recovers 5%, and the loss is in detection rather than scoring. Per
scenario numbers are in [Scoring](#scoring).

---

## Pipeline

| | Stage | Tool | Produces |
|---|---|---|---|
| 1 | Read filtering | fastp | trimmed FASTQ |
| 2 | Assembly | SPAdes (or MEGAHIT) | GFA graph |
| 3 | Cycle detection | DFS over the graph | candidate FASTA |
| 4 | Gene prediction | Prodigal (`-p meta`) | ORFs |
| 5 | Homology search | BLAST against four reference sets | hits |
| 6 | Scoring | length, structure, cargo quality | ranked table |
| 7 | Boundary annotation | IS/cargo split | GenBank |
| 8 | Long-read validation *(optional)* | minimap2 | circular-read evidence |

Stages 3–7 can be run alone on an existing GFA.

---

## Installation

```bash
git clone https://github.com/recepcanaltinbag/picota.git
cd picota

conda create -n picota -c bioconda -c conda-forge \
    python=3.9 prodigal blast spades megahit fastp \
    minimap2 samtools biopython pandas pyyaml tqdm
conda activate picota
pip install -r requirements.txt
```

The pipeline runs from the source directory; there is nothing further to install.

### Requirements

| | Tool | Version | Needed for |
|---|---|---|---|
| required | BLAST+ | ≥ 2.12 | scoring |
| required | Prodigal | ≥ 2.6 | gene prediction |
| required | SPAdes | ≥ 3.13 | assembly |
| optional | MEGAHIT | ≥ 1.2 | alternative assembler |
| optional | fastp | — | read filtering |
| optional | minimap2, samtools | — | long-read validation |

Linux, Python ≥ 3.8, ≥ 8 GB RAM (16 GB for large genomes), and disk for the
reference sets. macOS and Windows are untested.

---

## Quick start

### From a GFA you already have

```bash
cd picota/
conda run -n picota python3 test_complete_pipeline.py --gfa_mode
```

Runs detection, scoring and annotation over the bundled *P. nitroreducens*
graph (`picota/test_data/testNitro.gfa`) — no download, no assembly. Results
land in `picota_results/`.

### From SRA accessions

**1. Configure** `picota/config.yaml`:

```yaml
paths:
  outdir: "/path/to/output"
  sra_id_file: "picota/sra_ids.csv"
  path_to_antibiotics: "picota/DBs/Antibiotics/protein_fasta_protein_homolog_model.fasta"
  path_to_xenobiotics: "picota/DBs/Xenobiotics/Xenobiotics_kegg.fasta"
  path_to_ises: "picota/DBs/ISes/_tncentral_nointegrall_isfinder-TNs.fasta"
  path_to_TNs: "picota/DBs/CompTns/Known_Tns.fasta"
  assembler_type: "spades"

options:
  min_size_of_cycle: 1000
  total_score_type: 3          # bounded 0-100 score; see Scoring
  threshold_final_score: 50
```

**2. List the samples** in `picota/sra_ids.csv`:

```csv
sra_short_id,sra_long_id
SRR12345678,SRR98765432
SRR11111111,-
```

`sra_long_id` is an ONT accession for optional validation; use `-` when there
is none.

**3. Run:**

```bash
cd picota/
python3 picota/picota_testv3.py -c picota/config.yaml
```

### Individual stages

`picota/picota.py` exposes the stages separately:

```bash
python picota/picota.py sra_download --sra SRR11362851 --output data/
python picota/picota.py assembly     --fastq reads.fq  --output assembly/
python picota/picota.py analysis     --gfa assembly.gfa --output results/
python picota/picota.py scoring      --cycle_folder cycles/ --output results/
python picota/picota.py all          --sra SRR11362851 --output results/
```

### Cycle detection as a library

```python
from src.cycle_finderv2 import cycle_analysis

cycle_analysis(
    "assembly.gfa", "cycles.fasta",
    find_all_path=True, path_limit=25,
    min_size_of_cycle=1000, max_size_of_cycle=40000,
    name_prefix_cycle="Cycle",
    min_component_number=1, max_component_number=25,
    k_mer_sim=80, threshold_sim=80,
)
```

---

## Reference databases

Four reference sets are searched, configured by path in `config.yaml`.

1. **CARD** — antibiotic resistance genes.
   [card.mcmaster.ca](https://card.mcmaster.ca/), protein FASTA
   (`protein_fasta_protein_homolog_model.fasta`).

2. **Insertion sequences** — ISfinder and TnCentral merged, nucleotide FASTA
   (`_tncentral_nointegrall_isfinder-TNs.fasta`). Finding one is a necessary
   condition for a candidate to score at all.

3. **Known composite transposons** — nucleotide FASTA (`Known_Tns.fasta`). A
   cycle matching one is reported as `Known` rather than `Novel`.

4. **Xenobiotic degradation genes** (KEGG-derived)
   - Built by `picota/src/build_xenobiotic_db.py` from KEGG KO identifiers
     resolved to UniProt sequences, then clustered
   - Format: FASTA (proteins), 16,539 sequences
   - Header: `accession|KO:K16045|EC:1.1.1.145|PATH:map00984|description|organism`

   **Headers must not contain spaces.** BLAST's `sseqid` is the header up to its
   first space, so a description written with spaces is cut mid-word and the
   KO/EC provenance never reaches a result table — which is what happened to all
   16,539 headers of the first build. The builder writes descriptions and
   organism names with underscores for that reason. A custom set can use any
   header; PICOTA reads the KEGG fields when they are present and otherwise
   passes the name through unchanged.

---

## Output

```
output/
├── assembly/<accession>/     GFA assembly graphs
├── cycles/                   Candidate cycle sequences (FASTA)
├── scoring/<accession>/
│   ├── picota_final_tab      Tab-delimited raw results
│   ├── picota_enriched.csv   Enriched CT-tagged results
│   ├── genbank/              GenBank annotations
│   └── Pico_Temp/            Intermediate BLAST / Prodigal files
├── annot/<accession>/        Split transposon/cargo FASTA
└── mapping/<accession>/      BAM files and long-read analysis
```

### Raw results table (`picota_final_tab`)

| Column | Description |
|--------|-------------|
| `CycleID` | Cycle identifier: `Cycle_N-lenXXXX-compY-` |
| `score0/1/2/3` | Four scoring strategies, all always written (see Scoring section) |
| `NumIS` / `ISproducts` / `IScoords` | IS element count, names, coordinates |
| `NumAnt` / `Antproducts` / `Antcoords` | AMR gene count, names, coordinates |
| `NumXeno` / `Xenoproducts` / `Xenocoords` | Xenobiotic gene count, names, coordinates. With the KEGG-derived set a name reads `description\|KO:K16045\|EC:1.1.1.145` |
| `NumCompTN` / `CompTN` | Known composite transposon matches |

### Enriched output (`picota_enriched.csv`)

A tidy, analysis-ready CSV automatically generated alongside `picota_final_tab`.
Each CT gets a unique tag (`CT001`, `CT002`, …). When a CT carries genes from
multiple antibiotic classes it is expanded to one row per class.

| Column | Description |
|--------|-------------|
| `CT_Tag` | Unique composite transposon tag (`CT001`, …) |
| `Category` | `Novel` (not in CompTn DB) or `Known` |
| `CycleID` | Original cycle identifier |
| `SRA_ID` | Source sample accession |
| `CT_Length_bp` | Total cycle length in base pairs |
| `Score` | score0 — **not** the score the gate uses. `total_score_type` ships as 3, but the enriched CSV reads score0 and does not carry score3 at all. Read score3 from `picota_final_tab` until this is reconciled |
| `NumIS` | Number of IS elements detected |
| `IS_Group` | IS superfamily group (e.g., `IS6`, `IS3`) |
| `IS_Family` | IS family name (e.g., `IS26`, `ISEcp`) |
| `IS_Length_bp` | Length of the longest IS element found |
| `IS_Names` | Semicolon-separated IS element names |
| `NumAMR` | Number of AMR genes detected |
| `Antibiotic_Class` | Antibiotic class (one row per class) |
| `Resistance_Gene` | Representative resistance gene for this class |
| `NumXeno` | Number of xenobiotic genes detected |
| `Xenobiotic_Functions` | Semicolon-separated xenobiotic gene names |
| `NumCompTN` | Number of known CompTn database matches |
| `Known_CompTN` | Matched known composite transposon name(s) |

---

## Scoring

Detection returns every cycle the graph closes, which on a real chromosome means
roughly 20 candidates per genome before anything is filtered — most of them the
host's own repeat structures rather than composite transposons. Scoring is what
separates the two.

### Choosing a scoring mode

`total_score_type` selects the formula. **Use `3`.** The others are kept because
published results reference them.

| Mode | Range | Behaviour |
|------|-------|-----------|
| `0` | unbounded | Sums every hit score, then applies a length/component exponent. Ranks well but the value is not comparable between runs, and many weak hits outweigh one strong hit. |
| `1` | unbounded | Presence/absence variant of `0`. |
| `2` | unbounded | Requires cargo in a database, so an element whose cargo is novel is unreportable at any threshold. |
| `3` | **0–100** | Bounded, quality-based, and each term corresponds to a property the definition of a composite transposon asserts. |

### How score3 is computed

```
score3 = 100 × IS_gate × length_gate × (0.50 · component_fit + 0.50 · cargo_quality)
```

Two **gates**, which multiply and can zero the score outright:

- `IS_gate` — `0.5 + 0.5 × IS_hit_quality`, or **0** when no insertion sequence
  is found. An IS is a necessary condition, not a bonus. A hit at 59%
  identity-coverage is still unambiguously an IS, which is why the gate does not
  fall all the way to zero for a weak one.
- `length_gate` — how far the cycle sits from the expected composite transposon
  length (`mean_of_CompTns`, `std_of_CompTns`). Under `dist_type: 1` anything
  shorter than the mean pays nothing, so the gate only penalises oversized
  cycles. This gate carries most of the specificity: without it, precision on
  the benchmark falls from 69% to 47%.

Two **weighted terms**, equally split:

- `component_fit` — `1 / (1 + |components − 2| / 12)`. A cycle threading many
  graph components is wandering through host repeats. Note what this does *not*
  measure: it is not a structural score. Real elements thread 2 components when
  their flanking IS has two genomic copies and 16 when it has sixteen, so the
  count tracks how widespread the IS is. It is a host-repeat filter, and it is
  weighted accordingly.
- `cargo_quality` — best hit against CARD or the xenobiotic set, as
  (alignment length / subject length) × percent identity. When no database hit
  is found it falls back to a floor of **0.30** rather than zero, so an element
  carrying cargo absent from every database is still reportable. This is the one
  thing mode `2` cannot do.

The 0.50/0.50 split is measured, not chosen: swept across 1272 benchmark
candidates it gives the best F1, and its own optimal threshold lands on 50 —
the value `threshold_final_score` ships with.

### What the benchmark says

Measured on 480 implanted elements across 120 simulated genomes, plus 10
wild-type controls (see [docs/VALIDATION.md](docs/VALIDATION.md)):

| Scenario | Sensitivity | Precision |
|---|---|---|
| `baseline` — unique IS, CARD cargo | 92.5% | 97.4% |
| `novel_cargo` — cargo in no database | 90.0% | 100% |
| `compact` — small elements | 97.5% | 100% |
| `shared_is` — several CTs on one 16-copy IS | 100% | 88.9% |
| `cargo_is_diff` — a different IS inside the cargo | 95.0% | 74.5% |
| `cargo_is_same` — flanking IS repeated in the cargo | 5.0% | 2.5% |
| **All but `cargo_is_same`** | **95.0%** | **90.9%** |

On ten wild-type genomes with nothing implanted, 180 candidate cycles are
detected and **none** clear the threshold.

Two limits worth knowing before relying on the score:

- **`cargo_is_same` is not solved.** When the flanking IS also occurs inside the
  cargo, the element usually never becomes a candidate cycle at all — the loss
  is in detection, not scoring, and no threshold recovers it.
- **The IS gate did not discriminate on this benchmark.** Every one of the 1272
  candidates carried an IS hit, so the gate rejected nothing. It is kept because
  it is definitionally required, not because it was shown to filter.

---

## Performance

### BLAST

Scoring used to search a cycle at a time: for every cycle, one BLAST invocation
per reference database. BLAST's per-invocation setup — opening the database and
walking its index — does not amortise over a query of five predicted proteins,
so the cost was paid once per cycle per database. `blast_batch` (default true)
sends every cycle in one query instead.

Measured on 20 real cycles against the four bundled databases, BLAST at four
threads and the databases already built:

| search | per cycle | batched |
|---|---|---|
| blastn, IS elements | 1.96 s | 0.20 s |
| blastn, known transposons | 1.81 s | 0.20 s |
| blastp, CARD | 3.74 s | 1.52 s |
| blastp, xenobiotics | 14.01 s | 10.24 s |
| **total** | **23.5 s** (80 invocations) | **13.1 s** (4) |

`picota_final_tab` is byte-identical either way, checked on 4 and on 40 cycles.
An E-value depends on the database and on the individual query, not on what else
shares the query file, and the score is computed from alignment length, subject
length and identity — none of which the batching touches.

Threads follow from that: one large search parallelises where 160 tiny ones did
not. On the batched query of 40 cycles against the 82k-sequence xenobiotic set,
62.0 s at one thread, 15.4 s at four, 9.7 s at eight, 9.1 s at sixteen, 8.5 s at
twenty-four. The default is **8**, where the curve flattens and there is still
room for an assembler running beside it. `PICOTA_BLAST_THREADS` overrides it for
callers that score several samples at once.

The remaining cost is the xenobiotic search, which gains only 1.4× from batching
because it is dominated by real alignment work rather than by setup. A faster
aligner, not fewer invocations, is what would pay there.

### Re-running

Work is skipped when its output already exists, which is what makes a long run
resumable. Existence alone is the wrong test: an interrupted search leaves a
truncated table, and re-running after a database is rebuilt or the scoring code
changes would report the previous answer.

So each output carries a sidecar naming the inputs it came from:

```
Blast_Out/blast_files/SRR123_Xenobiotics.batch.out
Blast_Out/blast_files/SRR123_Xenobiotics.batch.out.picota-sig.json
```

The sidecar is written only after the work succeeds, so an interrupted run
leaves an output without one and the next run redoes it. Inputs are hashed
rather than stat'ed, so a checkout or an rsync that moves mtime without changing
a byte does not force a recomputation. Re-running 40 cycles with unchanged
inputs takes 5.8 s against 16.8 s for the first pass.

Nothing has to be deleted to pick up a rebuilt database or new code — the
signature no longer matches and the work is redone.

---

## Algorithm details

### Cycle detection

1. Parse GFA → directed graph where each contig becomes two nodes (`contig+`, `contig-`)
2. Run iterative DFS to detect back-edges (cycles) and reverse-complement cycles
3. For each detected cycle edge, find simple paths from source to destination
4. Reconstruct full cycle sequences, accounting for contig overlaps
5. Filter duplicates: strand-agnostic node-set comparison, then k-mer similarity (Jaccard-based with inverted index)

### Deduplication strategy

| Level | Method | Catches |
|-------|--------|---------|
| Path (strand-agnostic) | Frozenset of node IDs without `+/-` | Reverse-traversal duplicates |
| Path (weighted) | Shared node length / max length | Nested subpath duplicates |
| Sequence (k-mer) | Inverted index + Jaccard similarity | Near-identical sequences from different paths |

---

## Configuration reference

Every default below is the one in the shipped `config.yaml`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_size_of_cycle` | 1000 | Minimum cycle length (bp). The graph cycle is IS + cargo, not IS + cargo + IS, so a compact element yields a shorter cycle than itself — see [ROADMAP §9](docs/ROADMAP.md) |
| `max_size_of_cycle` | 40000 | Maximum cycle length (bp) |
| `min_component_number` | 1 | Minimum contigs in cycle |
| `max_component_number` | 25 | Maximum contigs in cycle |
| `k_mer_sim` | 80 | K-mer size for similarity filtering |
| `threshold_sim` | 80 | K-mer similarity threshold (%) |
| `find_all_path` | true | Enumerate all paths (exponential — use with caution) |
| `path_limit` | 25 | Maximum nodes in one enumerated path. Searches that hit it are cut short and counted; the run warns that the enumeration is not exhaustive |
| `dedup_mode` | legacy | `legacy` reproduces the historical deduplication; `strict` never discards a candidate for merely sharing a repeat node — see [ROADMAP §2](docs/ROADMAP.md) |
| `dedup_min_ani` | 99.0 | Identity (%) above which two same-sized candidates are one sequence. `strict` only |
| `mean_of_CompTns` | 5850 | Mean composite transposon length (bp) |
| `std_of_CompTns` | 2586 | Standard deviation of composite transposon length |
| `max_z` | 20 | Length z-score at which the length term reaches zero |
| `dist_type` | 1 | Length-fit distribution |
| `total_score_type` | 3 | Which score the threshold applies to (0-3). All four are always written — see [Scoring](#scoring) |
| `threshold_final_score` | 50 | Minimum score to report a cycle |
| `split_min_score` | 92 | Minimum score for boundary annotation |
| `assembler_type` | spades | Assembly tool (`spades` or `megahit`). SPAdes by measurement: 40/40 implanted elements recovered against 19/40 from the MEGAHIT graph of the same reads — see [ROADMAP §6](docs/ROADMAP.md) |
| `blast_batch` | true | One BLAST search per database over every cycle at once, rather than one per cycle per database. Same hits either way — see [Performance](#performance) |

---

## Testing

```bash
conda run -n picota python3 -m pytest tests/ -q          # everything
conda run -n picota python3 -m pytest tests/test_smoke.py -q   # module health
conda run -n picota python3 -m pytest tests/test_e2e_gfa.py -q # full post-assembly run
```

`tests/test_e2e_gfa.py` runs detection through scoring to enriched CSV over the
bundled graph and checks every output file and column. `tests/test_known_defects.py`
holds characterisation tests for the open defects in
[docs/ROADMAP.md](docs/ROADMAP.md), marked `xfail(strict=True)` so a fix turns
them into failures and forces the marker to be removed.

One test asserts a wall-clock budget for the k-mer filter and will fail on a
loaded machine without anything being wrong.


---

## Documentation

| | |
|---|---|
| [docs/ROADMAP.md](docs/ROADMAP.md) | Open defects, phase plan, and what each measurement settled |
| [docs/VALIDATION.md](docs/VALIDATION.md) | Benchmark design and ground truth |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Development setup and pull request process |

Bug reports: [GitHub issues](https://github.com/recepcanaltinbag/picota/issues).
Please include the Python version, OS, tool versions and a minimal example.

---

## Citation

> *A de novo computational pipeline for discovery and analysis of composite
> transposons from assembly graphs reveals hidden antibiotic resistance gene
> mobilization networks.* Manuscript in preparation.

[10.5281/zenodo.21769753](https://doi.org/10.5281/zenodo.21769753)

---

## License

MIT. See [LICENSE](LICENSE).
