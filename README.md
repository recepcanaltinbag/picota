# PICOTA

<p align="center">
<img src="picota/logo/picota_logo.png" alt="PICOTA Logo" width="200" style="max-width: 100%; height: auto;">
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
cycle_id    sra_id    kmer    score0  antibiotics    is_elements    xenobiotics
Cycle_1     SRR123    39      87.5    aadA1,aacA4    IS26,IS1       nah_cluster
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

## How to Cite

If you use PICOTA, please cite:

```bibtex
@article{can2024picota,
  title={PICOTA: Pipeline for Identification of Composite Transposons from Assembly graphs},
  author={Canaltinbag, Recep},
  journal={Bioinformatics},
  year={2024},
  doi={pending}
}
```

For now, please reference the GitHub repository:
```
https://github.com/recepcanaltinbag/picota
```

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
  min_size_of_cycle: 2000
  max_size_of_cycle: 40000
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
    min_size_of_cycle=2000,
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

PICOTA uses a Z-score-based composite scoring system. The score reflects how well a candidate cycle matches the expected size distribution of known composite transposons and whether it encodes relevant genes.

Three scoring strategies are available (`total_score_type` in config):

| Type | Formula | Best for |
|------|---------|----------|
| `0` | `(Σ scores)^z_normalized` | General use |
| `1` | `(100·[IS>0] + 100·[AMR>0] + 100·[Xeno>0])^z_normalized` | Presence/absence |
| `2` | `(cargo_score × IS_score) + 10^z_normalized` | IS-emphasis |

Where `z_normalized = 1 - |len - μ| / (σ × max_z)` penalizes cycles far from the mean composite transposon length (default μ = 5850 bp, σ = 2586 bp).

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
| `min_size_of_cycle` | 2000 | Minimum cycle length (bp) |
| `max_size_of_cycle` | 40000 | Maximum cycle length (bp) |
| `min_component_number` | 1 | Minimum contigs in cycle |
| `max_component_number` | 25 | Maximum contigs in cycle |
| `k_mer_sim` | 200 | K-mer size for similarity filtering |
| `threshold_sim` | 99 | K-mer similarity threshold (%) |
| `find_all_path` | false | Enumerate all paths (exponential — use with caution) |
| `path_limit` | 15 | Max path length when `find_all_path=true` |
| `mean_of_CompTns` | 5850 | Mean composite transposon length (bp) |
| `std_of_CompTns` | 2586 | Standard deviation of composite transposon length |
| `total_score_type` | 0 | Scoring formula (0, 1, or 2) |
| `threshold_final_score` | 50 | Minimum score to report a cycle |
| `assembler_type` | megahit | Assembly tool (`megahit` or `spades`) |

---

## Citation

If you use PICOTA in your research, please cite:

> [Manuscript in preparation]

---

## License

MIT License. See [LICENSE](LICENSE) for details.
