# PICOTA Validation Status and Benchmark Design

Status: 2026-09-01. Companion document: [ROADMAP.md](ROADMAP.md).

---

## 1. What has actually been validated

- **Unit and smoke tests** — 125 tests over parsing, cycle construction,
  deduplication helpers, scoring and output formatting.
- **End-to-end tests** on the bundled `test_data/testNitro.gfa`, checking that the
  pipeline runs and that output files and columns are well formed.
- **Long-read evidence, per candidate** — optional minimap2 mapping of Oxford
  Nanopore reads onto a reported cycle, used as supporting evidence that the
  circular structure exists. This is per-candidate corroboration, not a
  genome-wide measurement.
- **Known/Novel labelling** by homology against a small curated transposon set
  (Tn4352, Tn6010, Tn6023, Tn6309, Tn6925, TnPMLUA4).

The largest analysis result currently in the repository is
`picota_results/picota_enriched_combined.csv`: **three composite transposons from
one run** (SRR20032745, *Enterococcus faecium*).

## 2. What has not been validated

**PICOTA has never been benchmarked against a closed reference genome.** There is
no ground-truth CT catalogue, no recall or precision figure, and no comparison
against a completed assembly anywhere in this repository. Any statement of the
form "strain X carries N composite transposons and PICOTA recovered M of them"
would be unsupported today.

### A clarification worth recording

PICOTA does **not** operate on pre-selected, CT-carrying reads. `sra_download.py`
retrieves the complete SRA run with no filtering (`skip_filtering: true`), and
`assembly.py` assembles the whole genome; cycles are then sought in the
whole-genome assembly graph. The input has always been whole-genome short reads.
What is missing is the ground-truth comparison, not the genome-wide input.

---

## 3. Benchmark design (phase 0.5)

### 3.1 Strain selection

Required properties, in order:

1. A **closed** genome (long-read or hybrid assembly, single circular contig plus
   plasmids).
2. Short-read Illumina data **from the same isolate**, available in SRA. Same
   isolate is non-negotiable: a different strain makes the ground truth drift.
3. **Several copies of one IS family**, at least two of which flank *different*
   cargo. This is the discriminating case; a genome with one CT tests nothing
   that `testNitro.gfa` does not already cover.

Candidate groups to screen (to be confirmed, not yet selected): *E. coli* ST131
and *K. pneumoniae* lineages carrying IS26-bounded resistance regions,
*Salmonella* carrying IS26/Tn6029 arrays, and *A. baumannii* AbGRI islands
built on ISAba1 and IS26.

Strain selection is scripted so the criteria are reproducible rather than
anecdotal. `scripts/select_benchmark_strains.py` screens properties 1 and 2
against NCBI (closed RefSeq assemblies, then public Illumina paired-end WGS runs
from the same BioSample) and writes a shortlist TSV with the estimated read
depth each run would give.

It deliberately does **not** restrict to NCBI reference genomes. That filter
takes E. coli from 5350 closed assemblies down to 2, and the survivor is K-12
MG1655 -- a lab strain with essentially no resistance composite transposons,
which is the opposite of what this benchmark needs. Clinical MDR isolates are
where multi-copy IS elements with differing cargo live, and none of them are
reference genomes. `--reference-only` restores the restriction if it is ever
wanted.

```
python scripts/select_benchmark_strains.py \
    --email you@example.org \
    --species "Escherichia coli" "Klebsiella pneumoniae" \
    --limit 200 --min-coverage 40 --resolve-strains \
    --out benchmark_candidates.tsv
```

A run over 240 closed assemblies of these two species yielded 83 candidates
(33 *E. coli*, 50 *K. pneumoniae*) at 41x to 1894x estimated depth, median 97x.

`--resolve-strains` is optional because the strain designation lives only in the
BioSample record -- assembly `Biosource` comes back empty for most modern
submissions and SRA reports species only -- so filling it costs two extra
throttled requests per candidate. The assembly and BioSample accessions identify
the isolate without it; the strain name is for human-readable reporting. The
`Instrument` column is free and worth reading: model implies read length, and
read length drives how well a short-read assembly resolves the very repeats this
benchmark measures.

Property 3 -- several copies of one IS family with at least two carrying
different cargo -- is not screened there, because it needs the genome downloaded
and annotated. That is stage 2, described below.

### 3.2 Ground truth construction

1. Annotate IS elements on the closed genome (ISEScan, cross-checked against
   ISfinder).
2. Enumerate same-family IS pairs in the appropriate orientation separated by
   less than a chosen maximum cargo distance.
3. For each candidate CT record: IS family, IS copy number genome-wide, cargo
   coordinates, cargo gene content, and total length.
4. Group CTs that share an IS family, and flag those whose cargo differs — these
   are the cases the benchmark exists to measure.

### 3.3 Execution

Assemble the matching short reads with SPAdes, take the GFA, and run PICOTA
exactly as a user would. Repeat at subsampled depths (10x, 30x, 100x) to
characterise depth dependence.

### 3.4 Metrics

| Metric | Definition |
|--------|------------|
| **CT recall** | Fraction of ground-truth CTs recovered as a reported cycle at ≥95% identity over ≥95% of length |
| **Precision** | Fraction of reported cycles that map to a ground-truth CT |
| **Copy-distinctness recall** | Of the N ground-truth CTs sharing one IS family, how many are reported as *separate* candidates. **The headline metric** |
| **Cargo recall** | Fraction of ground-truth cargo genes present in at least one reported candidate |
| **Boundary accuracy** | Offset between reported and true transposon/cargo boundaries |

### 3.5 Simulated benchmark (implemented)

Downloading closed genomes gives realism but never certainty: the true CT
catalogue of a real isolate is itself an annotation with its own errors.
`scripts/simulate_ct_genome.py` inverts that trade. It builds a chromosome whose
composite transposons are known by construction, using **real sequences** — IS
elements from the bundled ISfinder set, cargo from CARD — so what is synthetic
is the arrangement, not the biology. `scripts/score_picota_benchmark.py` scores
PICOTA's output against that ground truth.

```
python scripts/simulate_ct_genome.py --out-dir sim/ \
    --backbone-length 400000 --n-cts 5 --shared-is 4 \
    --is-copies-outside 5 --is-divergence 0.5 --cargo-genes 1 --seed 3
art_illumina -ss HSXt -i sim/genome.fasta -p -l 150 -f 50 -m 350 -s 50 \
    -o sim/r  # empirical HiSeq X quality profile, not a flat error rate
spades.py -1 sim/r1.fq -2 sim/r2.fq -o sim/sp -k 55,77,99 --only-assembler
# run cycle_analysis on sim/sp/assembly_graph_with_scaffolds.gfa, then:
python scripts/score_picota_benchmark.py \
    --ground-truth sim/ground_truth.tsv \
    --ground-truth-fasta sim/ground_truth_cts.fasta \
    --cycles sim/cycles_strict.fasta
```

Because a random backbone contains no repeats other than the ones the script
implants, the ground truth is exact but the genome is easier to assemble than a
real one. **Recall measured that way is an upper bound** and should be reported
as such.

`--backbone-fasta` implants the same CTs into a real chromosome instead, which
removes that objection. The implanted CTs remain exact ground truth, but the
host genome brings its own IS elements -- K-12 MG1655 carries dozens -- and
those form cycles of their own. They are counted as unexplained by the precision
metric, so **precision against a real backbone is a lower bound**: it charges
PICOTA for finding real structures that simply are not in our ground-truth
list.

#### Results, 50x simulated Illumina, SPAdes k=55,77,99

**Real backbone: *E. coli* K-12 MG1655 (GCF_000005845.2), 6 CTs implanted, four
of them sharing one IS element present in 14 genome copies.** Genome 4,681,301 bp;
assembly graph 300 segments / 410 links.

| mode | CT recall | precision | copy-distinctness |
|---|---|---|---|
| legacy | 6/6 | 6/19 | 4/4 |
| strict | 6/6 | 7/24 | 4/4 |

Every implanted composite transposon was recovered, and all four sharing
ISNaoc2 came back as separate cycles. The low precision is the MG1655 native-IS
effect described above, not a failure to discriminate.

**Random backbone**

| case | IS / cargo | CT recall | precision | copy-distinctness |
|---|---|---|---|---|
| IS 1221 bp, cargo 2154 bp | legacy | 6/6 | 6/12 | 3/3 |
| | strict | 6/6 | 7/13 | 3/3 |
| IS 1655 bp, cargo 810 bp | legacy | 5/5 | 7/14 | 4/4 |
| | strict | 5/5 | 21/28 | 4/4 |
| IS 2496 bp, cargo 810 bp | legacy | 5/5 | 8/15 | 4/4 |
| | strict | 5/5 | 11/18 | 4/4 |

Two findings, both worth stating against expectation:

**PICOTA already handles the reviewer's scenario.** Composite transposons
sharing one IS element (present in 11–13 genome copies) with different cargo
were recovered separately in every case, in legacy mode as well as strict.

**The D2 recall failure did not reproduce here.** On idealised two-node graphs
legacy saturates at two candidates whatever the ground truth; on real assembly
graphs it did not, because the assembler fragments each cycle into 5–24 nodes
and the shared IS is therefore rarely more than 70% of the total node length —
the condition the legacy threshold needs to misfire. The defect is real, and
the conditions for it are narrower than the synthetic graphs implied.

What strict mode actually buys on real data is **precision**: a higher fraction
of reported cycles correspond to genuine composite transposons, with recall
unchanged.

#### Structural scenarios, twelve elements each

Six chromosomes, each implanting twelve elements into MG1655 and varying **one**
structural property while holding everything else fixed, so a drop between rows
is attributable to the structure rather than to chance. Scored under score3 at a
threshold of 50.

| scenario | structure | in graph | reported | TP | FN | FP | score |
|---|---|---|---|---|---|---|---|
| baseline | unique IS, CARD cargo | 12 | 12 | **12** | 0 | 0 | 79–94 |
| novel_cargo | cargo in no database | 12 | 12 | **12** | 0 | 0 | 58–73 |
| cargo_is_diff | a different IS inside the cargo | 15 | 14 | **12** | 0 | 1 | 80–91 |
| compact | cycle below `min_size_of_cycle` | 8 | 7 | 7 | 5 | 0 | 68–96 |
| shared_is | twelve elements on one 48-copy IS | 27 | 6 | 4 | 8 | 2 | 63–66 |
| cargo_is_same | flanking IS also inside the cargo | 22 | 22 | 0 | 12 | 22 | — |

Three things this table is meant to show, including the ones that do not
flatter the tool:

**Scoring is not the bottleneck anywhere.** Across all six scenarios and 72
elements the scoring stage loses exactly one candidate. Every other loss happens
earlier — at assembly, at graph traversal, or at the size filter. `shared_is`
reports 4 of 12 not because scoring rejects the other eight but because only
four survive the graph as recoverable cycles; all three scoring modes return the
same four.

**novel_cargo is the case score2 cannot do at all.** Twelve elements whose cargo
matches nothing in CARD or the xenobiotic set are all recovered at 58–73. Under
score2 every one of them scores 10 against a threshold of 50.

**cargo_is_same is a genuine failure, not a scoring artefact.** When the
flanking IS also occurs inside the cargo, the cycle that PICOTA reports is a
sub-element, and none of the twelve implanted elements is recovered whole. This
is a limit of reading cycles out of the graph and is not addressed by any
scoring change.

### 3.6 Negative control

Wild-type *E. coli* K-12 MG1655 with nothing implanted, run through the same
reads, assembly and pipeline: 17 cycles detected in the graph, **0 reported
under score0, 2 under score2 and 0 under score3**.

Read as a false-positive rate this would be misleading, because MG1655 is not a
composite-transposon-free genome — it carries dozens of native IS elements and
some may genuinely flank cargo. The two surviving candidates are most likely its
*ampC*-family genes, which CARD recognises. What the control establishes is
narrower and still worth having: the pipeline does not manufacture composite
transposons where none were placed.

Both are 35–38 kb cycles whose IS and CARD hit sit some 30 kb apart — an IS and
a resistance gene in one chromosomal neighbourhood, not an element bounded by
two copies of one IS. score3 scores them 31.1 and 26.6, because length gates
there rather than contributing 20% of the total; this is the case that motivated
that change (ROADMAP §7).

### 3.7 Does the assembly already answer the question?

If the elements can be read off the contigs, neither the assembly graph nor
PICOTA is needed. Measured on the same assemblies, counting an element as
recovered only when a *single* contig carries at least 95% of it in one
contiguous alignment:

| route | intact | on a 14-copy IS | on a 2-copy IS |
|---|---|---|---|
| SPAdes contigs | 3/8 | 0/4 | 3/4 |
| MEGAHIT contigs | 0/8 | 0/4 | 0/4 |
| PICOTA, SPAdes graph | **8/8** | **4/4** | **4/4** |

Across the eight-run sweep, 26 of 80 elements are intact in a SPAdes contig.
Assembly fails exactly where the IS is multi-copy: two copies leave two possible
pairings that paired-end links usually settle, fourteen do not. The one two-copy
element SPAdes could not resolve has the longest IS of the four, 1,828 bp against
1,304, 1,449 and 1,602 — consistent with repeat length as a second limit, on
four cases.

Summing several alignments from one contig would not do: a contig holding only
IS+cargo aligns twice to an IS-cargo-IS element, once at each flanking copy, and
appears complete while missing an entire IS.

### 3.8 Predicted baseline

Copy-distinctness recall is expected to be poor before phases 2 and 3 land. On
synthetic graphs the current deduplication caps output at two candidates
irrespective of ground truth (2/3, 2/4, 2/5 — defect D2 in
[ROADMAP.md](ROADMAP.md)). The benchmark exists to measure that on real data and
to demonstrate the improvement afterwards.
