# PICOTA — Changelog & Development Notes

All notable changes are documented here. Follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased] — v1.1.0 (in progress)

### Added
- **`output_formatter.py`** — new module that generates `picota_enriched.csv` alongside
  `picota_final_tab`. Each detected composite transposon receives a unique CT tag
  (CT001, CT002, …). Rows are expanded one-per-antibiotic-class. Includes IS-family
  inference from 30+ IS superfamily prefixes and antibiotic-class inference from
  23 gene-name patterns (Beta-lactam, Aminoglycoside, Tetracycline, MLS, …).
- **`picota_enriched.csv`** output: `CT_Tag`, `Category` (Novel/Known), `IS_Group`,
  `IS_Family`, `IS_Length_bp`, `Antibiotic_Class`, `Resistance_Gene`, `CT_Length_bp`,
  `Score`, `SRA_ID`, `Xenobiotic_Functions`, `Known_CompTN`.
- **Step-progress logging** in `picota_testv3.py` — five clearly labelled steps with
  elapsed time, step banners, and a per-sample summary block printed to stdout/log.
- **`test_e2e_gfa.py`** — real end-to-end integration tests using bundled
  `testNitro.gfa`. Covers cycle detection, BLAST scoring, enriched CSV validation,
  and unit tests for `output_formatter.py` (~25 tests, no mock data).
- **`test_complete_pipeline.py`** rewritten — replaced `MockAnalysis` with real
  pipeline calls (`cycle_analysis`, `scoring_main`, `output_formatter`). Supports
  `--gfa_mode` flag to run without SRA download or assembly.
- `test_sra_ids.csv` auto-created in `--gfa_mode` to avoid CLI errors.

### Fixed
- **`scoringv4ProtBlast.py:611`** — `IndexError: list index out of range` when the
  cycle FASTA basename has no underscore (e.g., `cycles.fasta`). Now uses a safe
  `_basename_parts[1] if len(...) > 1 else ''` guard.

### Changed
- `scoring_main()` now calls `write_enriched_csv()` automatically at the end.
  The enriched CSV is written to the same directory as `picota_final_tab`.
- README updated with enriched output column reference and E2E test instructions.

---

## [1.0.0] — 2026-03 (initial release)

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
- `test_e2e_gfa.py` — 25 E2E tests with real BLAST execution (new in v1.1)

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
                      Three scoring modes (score0/1/2)
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
    score0, score1, score2: float
```

---

## Known limitations / future work

- IS-family annotation relies on name prefixes; full taxonomy lookup against ISFinder
  API would be more accurate.
- Antibiotic-class inference uses regex on gene names; integrating CARD ARO ontology
  would cover edge cases.
- Organism field in enriched CSV is currently the SRA accession; NCBI metadata lookup
  would populate species names.
- `find_all_path=True` can be exponential for highly connected graph regions; a
  heuristic path-budget would help.
- `samtools` output requires `libncurses.so.5` on some Linux systems — install via
  `apt install libncurses5` or use the conda-forge build.
