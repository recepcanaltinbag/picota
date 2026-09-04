# PICOTA Validation Status and Benchmark Design

Status: 2026-09-01. Companion document: [ROADMAP.md](ROADMAP.md).

---

## 1. What has actually been validated

- **Unit and smoke tests** — 530 tests over parsing, cycle construction,
  deduplication helpers, scoring and output formatting.
- **Simulated benchmark with exact ground truth** — 520 composite transposons
  implanted into 130 simulated genomes across six architectures and seven IS
  copy levels, plus ten wild-type controls. Sensitivity 95.0% and precision
  90.9% excluding the one architecture the method cannot solve; zero false
  positives on the wild-type controls (§3.5–3.7).
- **Head-to-head against MobileElementFinder** on the same genomes and the same
  matching criterion (§3.9).
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

**PICOTA has never been benchmarked against a real closed genome.** The figures
in §3.5–3.7 come from simulated chromosomes, where ground truth is exact by
construction. That is the right design for comparing detection methods — on a
real isolate a missed element and a mis-annotated reference are
indistinguishable — but it does not establish performance on real data. Any
statement of the form "strain X carries N composite transposons and PICOTA
recovered M of them" would still be unsupported today.

**One host genome.** Every simulated case is built on *E. coli* K-12 MG1655.
Generalisation across species and GC content is untested; the pipeline accepts
any backbone, so this is a run that has not happened rather than a design limit.

**Inverted flanks and diverged flanks.** Every implanted element is in direct
orientation with 0.5% divergence between its flanking copies. Neither axis has
been varied.

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

Counts below are **elements**, forty per scenario: four implanted into each of
ten independently simulated genomes. Ten genomes rather than one crowded genome
matters — elements sharing a chromosome share its assembly and its read set, so
they are not independent measurements, and forty elements in one genome makes it
several times more repetitive than any real one.

Scored under score3 at a threshold of 50.

| scenario | structure | in graph | reported | lost detecting | lost scoring | FP | sensitivity | PPV |
|---|---|---|---|---|---|---|---|---|
| baseline | unique IS, CARD cargo | 40/40 | **37** | 0 | 3 | 1 | 92.5% | 97.4% |
| novel_cargo | cargo in no database | 40/40 | **36** | 0 | 4 | 0 | 90.0% | 100% |
| compact | small elements, one cargo gene | 40/40 | **39** | 0 | 1 | 0 | 97.5% | 100% |
| cargo_is_diff | a different IS inside the cargo | 38/40 | **38** | 2 | 0 | 13 | 95.0% | 74.5% |
| shared_is | four elements per genome on one 16-copy IS | 40/40 | **40** | 0 | 0 | 5 | 100% | 88.9% |
| cargo_is_same | flanking IS also inside the cargo | 2/40 | 2 | 38 | 0 | 77 | 5.0% | 2.5% |
| **all six** | | 200/240 | **192** | 40 | 8 | 96 | 80.0% | 66.7% |
| **excluding cargo_is_same** | | 198/200 | **190** | 2 | 8 | 19 | **95.0%** | **90.9%** |

"In graph" counts elements covered at ≥95% by some candidate cycle before any
scoring, so the two loss columns separate what detection lost from what scoring
lost. Four things this table is meant to show, including the ones that do not
flatter the tool:

**Scoring is not the bottleneck.** Of 240 implanted elements, scoring loses
**eight**. The other 40 losses happen before scoring sees them, and 38 of those
are the single `cargo_is_same` architecture. Outside it, detection loses 2 of
200 and scoring loses 8 of 198. Tuning the score would move eight elements; the
larger constraint is upstream, in what the assembly graph preserves.

**novel_cargo is the case score2 cannot do at all.** Cargo matching nothing in
CARD or the xenobiotic set is still recovered at 90%, because score3 applies a
floor of 0.30 when no database hit is found rather than zeroing the term. Under
score2 these elements score in single digits against a threshold of 50 and are
unreportable at any threshold.

**shared_is is now fully recovered.** Four elements per genome built on one
insertion sequence present in sixteen genomic copies: 40 of 40, at 100%
sensitivity. This is the case a contig-based method loses first — SPAdes contigs
recover none of these forty (§3.7).

**cargo_is_same is a genuine failure, and it is a detection failure.** When the
flanking IS also occurs inside the cargo, only 2 of 40 elements are covered by
any candidate cycle at all. The 77 false positives are unrelated cycles the
scorer accepted; no scoring change recovers an element detection never proposed.
PICOTA should not be trusted on elements whose flanking IS recurs within their
own cargo.

### 3.6 Negative control

Ten wild-type *E. coli* K-12 MG1655 genomes with nothing implanted, one per
seed, through identical reads, assembly and pipeline.

| stage | count |
|---|---|
| candidate cycles detected | 12–22 per genome, median 19 (180 total) |
| reported under score3 at threshold 50 | **0 of 180** |

Running ten genomes rather than one mattered: the single genome originally used
yielded 22 cycles, which turns out to be the top of the range rather than
typical.

Read as a false-positive rate this would mislead, because MG1655 is not a
composite-transposon-free genome — it carries dozens of native IS elements and
some may genuinely flank cargo. What the control establishes is narrower and
still worth having: the pipeline does not manufacture composite transposons
where none were placed.

The candidates are the same native repeat structures that make the unmatched
cycle rate in the implanted runs unremarkable — a genome with no implants at all
still closes 12–22 cycles, so those unmatched candidates are overwhelmingly the
host's, not artefacts of the implants.

Most of the rejection is done by the length gate. MG1655's native repeats close
large cycles, a median of 21 kb against 3.9 kb for implanted elements, and
score3 gates on length rather than weighting it, so an oversized candidate is
rejected outright however good its homology. Removing that gate drops precision
on the implanted benchmark from 91% to 47%.

### 3.7 Does the assembly already answer the question?

If the elements can be read off the contigs, neither the assembly graph nor
PICOTA is needed. Measured on the same assemblies, across 520 implanted elements
spanning six architectures and seven copy levels, counting an element as
recovered only when a *single* contig or cycle carries at least 95% of it in one
contiguous alignment:

| route | recovered | missed | |
|---|---|---|---|
| MEGAHIT contigs | 0/520 | 520 | 0% |
| SPAdes contigs | 94/520 | 426 | 18% |
| PICOTA, MEGAHIT graph | 188/520 | 332 | 36% |
| PICOTA, SPAdes graph | **480/520** | 40 | **92%** |

The gap is not uniform, and the conditions under which it opens are the point:

| architecture | IS copies | SPAdes contigs | SPAdes graph |
|---|---|---|---|
| baseline | 2 | 52% | 100% |
| baseline | 3 | 5% | 100% |
| baseline | 8 | 0% | 100% |
| compact | 2 | **78%** | 100% |
| compact | 3 | 7.5% | 100% |
| cargo_is_diff | 2 | 0% | 95% |
| shared_is | 16 | 0% | 100% |

**Copy number decides, not element size.** `compact` implants the smallest
elements in the benchmark — a 700–900 bp IS around one cargo gene — which is the
easiest thing a contig can carry whole, and at two copies it gives the best
contig recovery anywhere at 78%. One additional genomic copy of the flanking IS
takes it to 7.5%, the same factor `baseline` loses going from 52% to 5%.

**An internal repeat is fatal on its own.** `cargo_is_diff` carries a two-copy
flanking IS, the easy case, but a second different IS inside its cargo. Contig
recovery is zero.

Two copies is therefore close to the ceiling for a contig-based method, not a
typical case: insertion sequences in real genomes are rarely present in only two
copies. A contig-based tool run on real data should be expected to perform below
these figures rather than at them.

Summing several alignments from one contig would not do: a contig holding only
IS+cargo aligns twice to an IS-cargo-IS element, once at each flanking copy, and
appears complete while missing an entire IS.

### 3.8 Deduplication, measured rather than predicted

An earlier version of this section predicted that copy-distinctness recall would
be poor, because deduplication capped output at two candidates irrespective of
ground truth (defect D2). That prediction is superseded: under `dedup_mode:
strict` the benchmark recovers four distinct elements per genome in every
scenario except `cargo_is_same`, including `shared_is`, where four elements
share one insertion sequence present in sixteen genomic copies and all four are
returned (§3.5). The cap no longer applies.

### 3.9 Against MobileElementFinder

MobileElementFinder 1.1.2 (Johansson et al. 2021) predicts composite transposons
from assembled sequence, so both tools can be scored against one ground truth
without a format translation favouring either. Run on the same genomes, with the
same 95% coverage and identity criterion, across 240 elements at six copy levels:

| route | 2 copies | 3 | 4 | 8 |
|---|---|---|---|---|
| PICOTA, SPAdes graph | **100%** | **100%** | **100%** | **100%** |
| MEF on the true chromosome | 98% | 98% | 98% | 98% |
| MEF on SPAdes contigs | 52% | 5% | 10% | 0% |
| MEF on MEGAHIT contigs | 0% | 0% | 0% | 0% |

The second row is what makes the comparison fair rather than flattering. Run on
the true chromosome — a perfect assembly — MobileElementFinder recovers 98% at
every copy level. Copy number does not degrade the tool; it degrades the assembly
the tool has to read. Without that row, "the tool missed these" could not be
separated from "the assembly lost these".

Its detection reads assembled contigs — fastq input only annotates depth in its
GFF output — so it inherits whatever the assembler collapsed. The claim
generalises past this one tool to any contig-based detector.

Two findings from this comparison are worth recording separately:

**Matching rule.** An earlier version of this table credited MobileElementFinder
with six recoveries from MEGAHIT contigs. Those were an artefact of unioning
alignment positions: each prediction matched the element in two separate places,
once at each flanking IS, and the union made a fragment look complete. Their
longest single alignments cover 64–81% of the element. MobileElementFinder emits
the element as a linear sequence, so it is scored by the single-alignment rule.
PICOTA is scored on the union, and only because a cycle is the collapsed form —
it carries one copy of the flanking IS where the element has two, so it cannot
align to the whole element in one piece. Measured over 240 elements its longest
single alignment spans a median 75% against 73% predicted for exactly "element
minus one IS copy", with 205 of 240 within two points.

**Database coverage bounds both tools, asymmetrically.** PICOTA's single
`compact` miss is an element whose flanking ISHne6 is absent from its IS
reference; MobileElementFinder's consistent 39/40 on perfect genomes is ISVa13,
which PICOTA's reference does contain and which PICOTA recovers. Neither is an
algorithmic failure and neither is fixed by better traversal or scoring.

