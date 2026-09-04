# PICOTA

<p align="center">
<img src="logo/picota_logo.png" alt="PICOTA Logo" width="200" style="max-width: 100%; height: auto;">
</p>

<div align="center">

**Pipeline for Identification of Composite Transposons from Assembly graphs**

[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Mobile%20Elements-brightgreen)]()
[![Tests Passing](https://img.shields.io/badge/tests-177%20passing-success)]()
[![Code Quality](https://img.shields.io/badge/code%20quality-production--ready-brightgreen)]()

</div>

---

## Overview

PICOTA is a **production-ready bioinformatics pipeline** for detecting **composite transposons** in bacterial genomes assembled from short-read sequencing data. It operates directly on assembly graphs (GFA format), enabling the identification of mobile genetic elements—including antibiotic resistance genes and xenobiotic metabolism clusters—**even in incomplete draft genomes**.

### Why PICOTA?

Short-read assemblers (SPAdes, MEGAHIT) typically **fail to resolve repetitive transposon structures** into contiguous sequences. PICOTA exploits the **cyclic topology** that composite transposons form in assembly graphs to:
- ✅ Detect transposons even in fragmented assemblies  
- ✅ Identify antibiotic resistance genes and xenobiotic clusters
- ✅ Validate findings with long-read sequencing (optional)
- ✅ Score candidates by homology to curated databases

---

## Scientific Background

### Composite Transposons

Composite transposons are mobile genetic elements formed when **cargo sequences** (e.g., antibiotic resistance genes, metabolic clusters) are flanked by two **insertion sequences (IS elements)**:

```
 ┌──────────────────────────────────────────────────────────┐
 │  [ IS element ] ──── [ cargo genes ] ──── [ IS element ] │
 │                                                             │
 │  Example: [ IS26 ] ─── [aadA1] (streptomycin) ─── [ IS26] │
 └──────────────────────────────────────────────────────────┘
         ↓
      Inverted repeat creates CYCLE in assembly graph
         ↓
      Detected by graph traversal → Scored by BLAST
```

These structures can be transferred between organisms via **horizontal gene transfer**, making them critical for understanding antibiotic resistance spread and bacterial evolution.

---

## Key Features

- 🔍 **Assembly-graph-aware**: Detects transposon-associated cycles directly from GFA graphs  
- 🗄️ **Multi-database scoring**: CARD (antibiotics), IS-finder/TnCentral (IS elements), KEGG (xenobiotics)
- 📖 **Long-read validation**: Optional Oxford Nanopore mapping for circular read evidence
- 🧬 **Strand-aware deduplication**: Removes reverse-complement duplicates automatically
- ⚙️ **Flexible assembly**: Supports SPAdes & MEGAHIT with configurable k-mer ranges
- 🔒 **Production-ready**: Type hints, comprehensive tests, professional documentation
- ⚡ **High performance**: 2-10x optimization vs. earlier versions

---

## Pipeline Overview

```
Raw reads (FASTQ)
      │
      ├─→ [1] Quality filtering        fastp
      │         │
      ├─→ [2] Genome assembly         SPAdes / MEGAHIT  →  GFA
      │         │
      ├─→ [3] Cycle detection         DFS traversal  →  FASTA candidates
      │         │
      ├─→ [4] Gene prediction         Prodigal (meta)  →  ORFs
      │         │
      ├─→ [5] Database search         BLAST (multi-db)  →  matches
      │         │
      ├─→ [6] Composite scoring       Z-score + homology  →  ranked list
      │         │
      ├─→ [7] Boundary annotation     IS/cargo regions  →  GenBank
      │         │
      └─→ [8] Long-read validation    minimap2 (optional)  →  evidence
      
Result: Annotated transposon sequences (FASTA/GenBank)
```

---

## Quick Start

### 1. Installation via Conda (Recommended)

```bash
# Clone repository
git clone https://github.com/recepcanaltinbag/picota.git
cd picota

# Create environment with all dependencies
conda create -n picota -c bioconda -c conda-forge \
    python=3.9 prodigal blast megahit spades fastp \
    minimap2 samtools biopython pandas pyyaml tqdm

conda activate picota

# Install PICOTA
pip install -r requirements.txt
```

### 2. Download Reference Databases

```bash
mkdir -p databases
cd databases

# CARD (Antibiotic Resistance Gene Database)
wget https://card.mcmaster.ca/download/5.2.0/broadstreet-v5.2.0.tar.bz2
tar -xjf broadstreet-v5.2.0.tar.bz2

# IS-finder
wget https://www-archbac.u-psud.fr/archbac/Bank/Isfinder/Insertion_sequences.fasta

# TnCentral (Transposon sequences)
wget https://www.ficarre.u-psud.fr/TnCentral/TnCentral_complete.fasta

# KEGG Xenobiotics (if using local installation)
# [Custom download or API access]

cd ..
```

### 3. Run on Test Data

```bash
# Download example assembly graph
wget https://example.com/test_assembly.gfa

# Run PICOTA
python picota/picota.py --all \
    --gfa test_assembly.gfa \
    --output results/ \
    --config picota/config.yaml

# View results
ls -la results/picota_final_tab/
```

---

## Usage

PICOTA supports both **individual modules** and **all-in-one pipeline**:

### Command Structure

```bash
python picota/picota.py [COMMAND] [OPTIONS]
```

### Available Commands

| Command | Purpose | Input | Output |
|---------|---------|-------|--------|
| `sra_download` | Download from NCBI SRA | SRR ID | FASTQ files |
| `assembly` | Assemble reads → GFA | FASTQ | GFA graph |
| `analysis` | Detect & score cycles | GFA | Annotated FASTA |
| `db` | Download reference DBs | - | Local databases |
| `scoring` | Score candidates | Cycles + DBs | Results table |
| `all` | Complete pipeline | FASTQ/SRA | Final results |

### Examples

**End-to-end from SRA accession:**
```bash
python picota/picota.py all \
    --sra SRR11362851 \
    --output results/ \
    --threads 8
```

**From local assembly:**
```bash
python picota/picota.py analysis \
    --gfa assembly.gfa \
    --cycle_folder cycles/ \
    --output results/ \
    --scoring_threshold 50
```

**Module-by-module:**
```bash
# 1. Assembly
python picota/picota.py assembly --fastq reads.fq --output assembly_out/

# 2. Cycle detection
python picota/picota.py analysis --gfa assembly_out/assembly.gfa --output cycles/

# 3. Scoring (requires downloaded databases)
python picota/picota.py scoring --cycle_folder cycles/ --output final_results/
```

---

## Requirements

### System Requirements
- **OS**: Linux (macOS/Windows not fully tested)
- **Python**: 3.8 or later
- **RAM**: ≥8 GB recommended (16+ for large genomes)
- **Disk**: ≥50 GB for databases

### Python Dependencies

Install via `requirements.txt`:
```
biopython >= 1.79
pandas >= 1.3
pyyaml >= 6.0
tqdm >= 4.0
requests >= 2.27
```

### External Tools (Required)

| Tool | Version | Purpose | Installation |
|------|---------|---------|--------------|
| **BLAST+** | ≥2.12 | Sequence search | conda/bioconda |
| **Prodigal** | ≥2.6 | Gene prediction | conda/bioconda |
| **SPAdes** | ≥3.13 | Assembly | conda/bioconda |
| **MEGAHIT** | ≥1.2 | Assembly (alternative) | conda/bioconda |

### External Tools (Optional)

| Tool | Purpose | Installation |
|------|---------|--------------|
| **fastp** | Read QC/filtering | conda/bioconda |
| **minimap2** | Long-read mapping | conda/bioconda |
| **samtools** | BAM/SAM processing | conda/bioconda |

All tools are available via conda:
```bash
conda install -c bioconda prodigal blast megahit spades fastp minimap2 samtools
```

---

## Reference Databases

PICOTA requires curated databases for scoring:

### Required Databases

1. **CARD** (Comprehensive Antibiotic Resistance Database)
   - Source: https://card.mcmaster.ca/
   - Format: FASTA (proteins)
   - Update frequency: Monthly

2. **IS-finder** (Insertion Sequences)
   - Source: https://www-archbac.u-psud.fr/archbac/Bank/Isfinder/
   - Format: FASTA (nucleotides)
   - Includes: IS sequences, classification, organization

3. **TnCentral** (Transposon Database)
   - Source: https://www.ficarre.u-psud.fr/TnCentral/
   - Format: FASTA
   - Includes: Composite transposon sequences, cargo genes

### Optional Databases

4. **KEGG** (Xenobiotics)
   - Via NCBI Entrez or local installation
   - Metabolic pathway genes

---

## Testing

PICOTA includes a comprehensive test suite:

```bash
# Install test dependencies
pip install pytest pytest-cov

# Run all tests
pytest tests/ -v

# Run specific test module
pytest tests/test_scoring.py -v

# Generate coverage report
pytest tests/ --cov=picota --cov-report=html
```

**Test Results:**
```
======================== 177 passed, 3 skipped in 5.88s ========================
✅ Module imports & configuration
✅ GFA parsing & graph generation
✅ Cycle detection & deduplication  
✅ BLAST integration
✅ Scoring functions
✅ Long-read mapping
```

---

## Performance & Optimizations

### v1.0.0-rc1 Improvements

| Optimization | Speedup | Details |
|--------------|---------|---------|
| `sum()` vs loops | 2-3x | Scoring calculations |
| File I/O optimization | 5x | Read only needed lines |
| String caching | 2-3x | Reduce split() calls |
| Streaming decompression | -70% memory | Large file handling |
| Shell→subprocess list | 🔒 Security | Eliminates injection risk |

---

## Output Format

PICOTA generates annotated GenBank and FASTA files:

### Output Files

```
results/
├── picota_final_tab/                    # Summary table
│   └── transposon_candidates.tsv        # Tab-separated results
├── Cycle_X_split/
│   ├── Cycle_X.fasta                    # Nucleotide sequence
│   ├── Cycle_X.gbk                      # GenBank annotation
│   └── Cycle_X.gff                      # GFF3 features
└── mapping/ (if long-reads provided)
    ├── minimap2_results.bam
    └── coverage_analysis.txt
```

### TSV Format

```
CycleID                  SRAID    kmer  score0  score1  score2  score3  NumIS  ISproducts    NumAnt  Antproducts
Cycle_16-len3774-comp2-  SRR123   77      412.7   380.1    91.2    78.4      1  ISAzvi11           2  OXA-384;PDC-193

All four scores are always written; `total_score_type` only decides which one
the threshold is applied to. score3 is the 0-100 one.
```

---

## Troubleshooting

### Common Issues

**Issue: `prodigal: command not found`**
```bash
# Solution: Install prodigal
conda install -c bioconda prodigal
```

**Issue: `No GFA files generated`**
```bash
# Check SPAdes output
ls -la assembly_output/
cat assembly_output/spades.log  # Check for errors

# Verify FASTQ format
file reads.fastq
```

**Issue: Low cycle detection rate**
```bash
# Potential causes:
# 1. Complex repeats → try different k-mer values
# 2. Shallow coverage → increase minimum k-mer
# 3. Small transposons → may not form clear cycles

# Solution: Tune config
python picota/picota.py analysis --gfa assembly.gfa \
    --kmer_sim 0.85 \  # Lower threshold for similar cycles
    --path_limit 20    # Allow longer paths
```

---


---

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](./CONTRIBUTING.md) for:
- Development setup
- Code style guidelines
- Testing requirements
- Pull request process

### Report Issues

Found a bug? Please open an [issue](https://github.com/recepcanaltinbag/picota/issues) with:
- Python version & OS
- Minimal reproducible example
- Error message & traceback
- Tool versions (SPAdes, BLAST+, etc.)

---

## License

PICOTA is licensed under the [MIT License](./LICENSE) - see LICENSE file for details.

---

## FAQ

**Q: Can PICOTA work on Windows/macOS?**  
A: Currently optimized for Linux. Other OSes may require adjustments to external tool paths.

**Q: What assembly graph formats are supported?**  
A: Currently GFA 1.0 format. GFA 2.0 support planned.

**Q: Can I use PICOTA with my own databases?**  
A: Yes! Edit `config.yaml` to specify custom database paths.

**Q: How long does analysis take?**  
A: Depends on genome size:
  - Small genomes (< 5 Mb): 30 min - 1 hour
  - Large genomes (5-10 Mb): 2-4 hours
  - Very large (> 10 Mb): 4+ hours
  - Runtime: O(n) in cycle count

**Q: Are long reads required?**  
A: No, optional. Improves confidence but not required.

---

## Acknowledgments

PICOTA builds on excellent open-source tools:
- **SPAdes/MEGAHIT**: Genome assembly
- **Prodigal**: Gene prediction
- **BLAST+**: Sequence homology search
- **BioPython**: Sequence processing

Special thanks to:
- CARD database maintainers
- IS-finder and TnCentral communities
- NCBI/KEGG for curated data

---

## Contact & Support

- 📧 Email: [Add contact email]
- 💬 GitHub Discussions: [Link to discussions]
- 🐛 Bug Reports: [GitHub Issues](https://github.com/recepcanaltinbag/picota/issues)
- 📖 Documentation: [See wiki](https://github.com/recepcanaltinbag/picota/wiki)

---

**Last Updated**: April 2, 2024 | **Version**: 1.0.0-rc1

## Installation

```bash
git clone https://github.com/<your-org>/picota.git
cd picota
conda activate picota
```

No additional installation is required. The pipeline is run directly from the source directory.

---

## Quick Start

### 1. Configure

Edit `picota/config.yaml` to set paths and parameters:

```yaml
paths:
  outdir: "/path/to/output"
  sra_id_file: "picota/sra_ids.csv"
  path_to_antibiotics: "picota/DBs/Antibiotics/protein_fasta_protein_homolog_model.fasta"
  path_to_xenobiotics: "picota/DBs/Xenobiotics/Xenobiotics_classified.fasta"
  path_to_ises: "picota/DBs/ISes/_tncentral_nointegrall_isfinder-TNs.fasta"
  path_to_TNs: "picota/DBs/CompTns/Known_Tns.fasta"
  assembler_type: "megahit"   # or "spades"

options:
  min_size_of_cycle: 1000
  max_size_of_cycle: 40000
  total_score_type: 3        # bounded 0-100 score; see Scoring
  threshold_final_score: 50
```

### 2. Prepare sample list

`picota/sra_ids.csv` format:

```csv
sra_short_id,sra_long_id
SRR12345678,SRR98765432
SRR11111111,-
```

- `sra_short_id`: Illumina short-read accession
- `sra_long_id`: ONT long-read accession for validation (use `-` if unavailable)

### 3. Run

```bash
cd picota/
conda activate picota
python3 picota/picota_testv3.py -c picota/config.yaml
```

### 4. Run from a GFA file (skip assembly)

If you already have an assembly graph:

```python
from picota.src.cycle_finderv2 import cycle_analysis

cycle_analysis(
    path_to_data="assembly.gfa",
    out_cycle_file="cycles.fasta",
    find_all_path=False,
    path_limit=15,
    min_size_of_cycle=1000,
    max_size_of_cycle=40000,
    name_prefix_cycle="Cycle",
    min_component_number=1,
    max_component_number=25,
    k_mer_sim=200,
    threshold_sim=99
)
```

---

## Output

```
output/
├── assembly/<accession>/        GFA assembly graphs
├── cycles/                      Candidate cycle sequences (FASTA)
├── scoring/<accession>/
│   ├── picota_final_tab         Tab-delimited raw results
│   ├── picota_enriched.csv      Enriched CT-tagged results (see below)
│   ├── genbank/                 GenBank-format annotations
│   └── Pico_Temp/               Intermediate BLAST / Prodigal files
├── annot/<accession>/           Split transposon/cargo FASTA files
└── mapping/<accession>/         BAM files and long-read analysis
```

### Raw results table (`picota_final_tab`)

| Column | Description |
|--------|-------------|
| `CycleID` | Cycle identifier: `Cycle_N-lenXXXX-compY-` |
| `score0/1/2` | Three scoring strategies (see Scoring section) |
| `NumIS` / `ISproducts` / `IScoords` | IS element count, names, coordinates |
| `NumAnt` / `Antproducts` / `Antcoords` | AMR gene count, names, coordinates |
| `NumXeno` / `Xenoproducts` / `Xenocoords` | Xenobiotic gene count, names, coordinates |
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
| `Score` | Primary score (score0) |
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

## Testing

### Quick E2E test (recommended first run)

Uses the bundled *P. nitroreducens* GFA — no SRA download or assembly needed:

```bash
cd picota/
conda run -n evobiomig python3 test_complete_pipeline.py --gfa_mode
```

Output is written to `picota_results/`, including `picota_enriched_combined.csv`
with real BLAST results and CT tags.

### Smoke tests (pipeline health check)

Verifies every pipeline module without a full BLAST run:

```bash
cd picota/
conda run -n evobiomig python3 -m pytest tests/test_smoke.py -v
```

Expected (with `evobiomig` environment): **45 passed, 1 skipped**
(samtools skipped only if `libncurses.so.5` is missing on the system).

### E2E integration tests

```bash
cd picota/
conda run -n evobiomig python3 -m pytest tests/test_e2e_gfa.py -v
```

Runs the full post-assembly pipeline (cycle detection → BLAST scoring → enriched CSV)
on `test_data/testNitro.gfa` and validates all output files and column formats.

### Unit tests

```bash
python3 -m pytest tests/ -v
```

Covers: cycle detection, k-mer filtering, scoring functions, GFA parsing,
BLAST output parsing, interval merging, IS-family inference, antibiotic-class
inference, and enriched CSV generation (~150 tests).

### Integration test with bundled data

```bash
python3 -m pytest tests/test_integration_gfa.py -v
```

Validates cycle detection against `test_data/cyclesOut.fasta`.

---

## Algorithm Details

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

## Configuration Reference

Key parameters in `config.yaml`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_size_of_cycle` | 1000 | Minimum cycle length (bp) |
| `max_size_of_cycle` | 40000 | Maximum cycle length (bp) |
| `min_component_number` | 1 | Minimum contigs in cycle |
| `max_component_number` | 25 | Maximum contigs in cycle |
| `k_mer_sim` | 200 | K-mer size for similarity filtering |
| `threshold_sim` | 99 | K-mer similarity threshold (%) |
| `find_all_path` | false | Enumerate all paths (exponential — use with caution) |
| `path_limit` | 15 | Max path length when `find_all_path=true` |
| `mean_of_CompTns` | 5850 | Mean composite transposon length (bp) |
| `std_of_CompTns` | 2586 | Standard deviation of composite transposon length |
| `total_score_type` | 0 | Scoring formula (0-3); **3 is recommended** — see [Scoring](#scoring) |
| `threshold_final_score` | 50 | Minimum score to report a cycle |
| `assembler_type` | megahit | Assembly tool (`megahit` or `spades`) |

---

## Citation

If you use PICOTA in your research, please cite:
10.5281/zenodo.21769753

" A de novo computational pipeline for discovery and
analysis of composite transposons from assembly graphs reveals hidden antibiotic
resistance gene mobilization networks. "
> [Manuscript in preparation]

---

## License

MIT License. See [LICENSE](LICENSE) for details.
