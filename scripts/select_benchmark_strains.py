#!/usr/bin/env python3
"""
select_benchmark_strains.py
---------------------------
Stage 1 of the PICOTA validation benchmark (docs/ROADMAP.md phase 0.5).

Finds strains that can serve as ground truth for measuring how well PICOTA
reconstructs composite transposons. A usable benchmark strain needs:

  1. A CLOSED genome -- assembly level "Complete Genome", so every composite
     transposon and every IS copy is resolved and countable.
  2. Illumina paired-end WGS reads from the SAME BioSample. Same isolate is
     non-negotiable: reads from a different strain make the ground truth drift.
  3. Enough read depth that a short-read assembly is meaningful.

Property 3 of the benchmark design -- several copies of one IS family, at least
two flanking different cargo -- is NOT screened here. That requires downloading
the genome and running IS annotation, which is stage 2. This script produces the
shortlist that stage 2 screens.

Usage:
  python scripts/select_benchmark_strains.py \
      --email you@example.org \
      --species "Escherichia coli" "Klebsiella pneumoniae" \
      --limit 200 --min-coverage 40 \
      --out benchmark_candidates.tsv

  # NCBI allows 10 requests/second with a key, 3 without:
  python scripts/select_benchmark_strains.py --email ... --api-key $NCBI_API_KEY ...

Output TSV columns:
  Species, Strain, AssemblyAccession, AssemblyLevel, GenomeLength, BioSample,
  RunAccession, Platform, LibraryLayout, TotalBases, EstimatedCoverage

The network layer is isolated in EntrezClient so the parsing and ranking logic
can be tested offline (see tests/test_select_benchmark_strains.py).
"""

import argparse
import os
import re
import sys
import time
from collections import namedtuple

AssemblyCandidate = namedtuple(
    "AssemblyCandidate",
    "accession organism strain level genome_length biosample")

SraRun = namedtuple(
    "SraRun",
    "accession platform layout total_bases total_spots")

BenchmarkCandidate = namedtuple(
    "BenchmarkCandidate",
    "assembly run estimated_coverage")

OUTPUT_COLUMNS = [
    "Species", "Strain", "AssemblyAccession", "AssemblyLevel", "GenomeLength",
    "BioSample", "RunAccession", "Platform", "LibraryLayout", "TotalBases",
    "EstimatedCoverage",
]

# Species whose resistance regions are commonly built on multi-copy IS elements,
# which is what makes them informative for this benchmark. Not exhaustive, and
# deliberately overridable: the point of this script is that strain selection is
# reproducible, not that this list is authoritative.
DEFAULT_SPECIES = [
    "Escherichia coli",
    "Klebsiella pneumoniae",
    "Salmonella enterica",
    "Acinetobacter baumannii",
]

_STRAIN_RE = re.compile(r"strain=([^;]+)")


def build_assembly_query(species):
    """Search term for closed, current RefSeq genomes of one species."""
    return (f'"{species}"[Organism] AND "complete genome"[Assembly Level] '
            'AND "latest refseq"[filter] AND "reference genome"[filter] NOT anomalous[filter]')


def build_assembly_query_broad(species):
    """As build_assembly_query, without the reference-genome restriction."""
    return (f'"{species}"[Organism] AND "complete genome"[Assembly Level] '
            'AND "latest refseq"[filter] NOT anomalous[filter]')


def build_sra_query(biosample):
    """Search term for public Illumina paired-end WGS runs of one BioSample."""
    return (f'{biosample}[BioSample] AND illumina[Platform] AND wgs[Strategy] '
            'AND paired[Layout] AND public[Access]')


def _to_int(value, default=0):
    try:
        return int(str(value).replace(",", ""))
    except (TypeError, ValueError):
        return default


def parse_assembly_summary(record):
    """
    Turn one Entrez assembly esummary record into an AssemblyCandidate.

    Returns None when the record is not a closed genome or is missing the
    BioSample link that stage 2 needs to pair reads with the genome.
    """
    level = record.get("AssemblyStatus", "")
    if level != "Complete Genome":
        return None

    biosample = record.get("BioSampleAccn", "").strip()
    if not biosample:
        return None

    strain = ""
    for entry in record.get("Biosource", {}).get("InfraspeciesList", []) or []:
        if entry.get("Sub_type") == "strain":
            strain = entry.get("Sub_value", "").strip()
            break
    if not strain:
        match = _STRAIN_RE.search(record.get("Biosource", {}).get("Isolate", "") or "")
        strain = match.group(1).strip() if match else ""

    # Genome length lives in the assembly stats blob under different keys
    # depending on the record vintage; total_length is the one we want.
    genome_length = _to_int(record.get("Meta_total_length"))
    if not genome_length:
        meta = record.get("Meta", "") or ""
        match = re.search(r'category="total_length"[^>]*>(\d+)<', meta)
        genome_length = _to_int(match.group(1)) if match else 0

    return AssemblyCandidate(
        accession=record.get("AssemblyAccession", ""),
        organism=record.get("SpeciesName", "") or record.get("Organism", ""),
        strain=strain,
        level=level,
        genome_length=genome_length,
        biosample=biosample,
    )


def parse_sra_summary(record):
    """
    Extract runs from one Entrez SRA esummary record.

    The useful fields are buried in two XML fragments: ExpXml carries the
    platform and library layout, Runs carries one entry per sequencing run.
    """
    exp_xml = record.get("ExpXml", "") or ""
    runs_xml = record.get("Runs", "") or ""

    platform_match = re.search(r"<Platform[^>]*>([^<]+)</Platform>", exp_xml)
    platform = platform_match.group(1).strip() if platform_match else ""
    if not platform:
        instrument = re.search(r'instrument_model="([^"]+)"', exp_xml)
        platform = instrument.group(1).strip() if instrument else ""

    layout = "PAIRED" if "<PAIRED" in exp_xml else ("SINGLE" if "<SINGLE" in exp_xml else "")

    runs = []
    for run_match in re.finditer(r"<Run\s+([^>]+)/?>", runs_xml):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', run_match.group(1)))
        accession = attrs.get("acc", "")
        if not accession:
            continue
        runs.append(SraRun(
            accession=accession,
            platform=platform,
            layout=layout,
            total_bases=_to_int(attrs.get("total_bases")),
            total_spots=_to_int(attrs.get("total_spots")),
        ))
    return runs


def estimate_coverage(total_bases, genome_length):
    """Read depth a run would give over a genome of this size, or None."""
    if not total_bases or not genome_length:
        return None
    return total_bases / genome_length


def pick_best_run(runs, genome_length, min_coverage):
    """
    Choose the run that best suits a short-read assembly benchmark.

    Prefers paired Illumina runs that clear `min_coverage`, then the smallest
    such run: excess depth costs download time and assembly memory without
    making the graph more informative.
    """
    eligible = []
    for run in runs:
        if run.layout and run.layout != "PAIRED":
            continue
        coverage = estimate_coverage(run.total_bases, genome_length)
        if coverage is None or coverage < min_coverage:
            continue
        eligible.append((coverage, run))
    if not eligible:
        return None, None
    coverage, run = min(eligible, key=lambda pair: pair[0])
    return run, coverage


class EntrezClient:
    """Thin wrapper over Bio.Entrez, isolated so the logic above stays testable."""

    def __init__(self, email, api_key=None, sleep=None):
        from Bio import Entrez  # imported lazily: only the live path needs it

        self._entrez = Entrez
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        # NCBI permits 10 requests/second with an API key, 3 without.
        self._sleep = sleep if sleep is not None else (0.11 if api_key else 0.34)

    def _throttle(self):
        time.sleep(self._sleep)

    def esearch(self, db, term, retmax):
        self._throttle()
        with self._entrez.esearch(db=db, term=term, retmax=retmax) as handle:
            # validate=False: NCBI ships esummary tags that are absent from the
            # published DTD (AnnotRptUrl, for one), and strict parsing rejects
            # the whole response over them.
            return self._entrez.read(handle, validate=False).get("IdList", [])

    def esummary(self, db, ids):
        if not ids:
            return []
        self._throttle()
        with self._entrez.esummary(db=db, id=",".join(ids)) as handle:
            result = self._entrez.read(handle, validate=False)
        # The assembly db nests its records; sra returns a plain list.
        if isinstance(result, dict):
            return result.get("DocumentSummarySet", {}).get("DocumentSummary", [])
        return result


def collect_candidates(client, species_list, limit, min_coverage,
                       broad=False, progress=None):
    """
    Find one benchmark candidate per closed genome that has usable reads.

    `client` needs esearch(db, term, retmax) and esummary(db, ids); anything
    satisfying that works, which is what makes this testable without network.
    """
    candidates = []
    for species in species_list:
        query = (build_assembly_query_broad if broad else build_assembly_query)(species)
        assembly_ids = client.esearch("assembly", query, limit)
        if progress:
            progress(f"{species}: {len(assembly_ids)} closed assemblies")

        for chunk_start in range(0, len(assembly_ids), 50):
            chunk = assembly_ids[chunk_start:chunk_start + 50]
            for record in client.esummary("assembly", chunk):
                assembly = parse_assembly_summary(record)
                if assembly is None:
                    continue

                sra_ids = client.esearch("sra", build_sra_query(assembly.biosample), 20)
                if not sra_ids:
                    continue

                runs = []
                for sra_record in client.esummary("sra", sra_ids):
                    runs.extend(parse_sra_summary(sra_record))

                run, coverage = pick_best_run(runs, assembly.genome_length, min_coverage)
                if run is None:
                    continue

                candidates.append(BenchmarkCandidate(assembly, run, coverage))
                if progress:
                    progress(f"  {assembly.accession} {assembly.strain or '-'} "
                             f"-> {run.accession} ({coverage:.0f}x)")
    return candidates


def write_candidates(path, candidates):
    """Write the shortlist as TSV, deepest-coverage-last for stable diffs."""
    ordered = sorted(candidates, key=lambda c: (c.assembly.organism,
                                                c.assembly.accession))
    with open(path, "w") as handle:
        handle.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for candidate in ordered:
            assembly, run = candidate.assembly, candidate.run
            handle.write("\t".join([
                assembly.organism,
                assembly.strain or "NA",
                assembly.accession,
                assembly.level,
                str(assembly.genome_length),
                assembly.biosample,
                run.accession,
                run.platform or "NA",
                run.layout or "NA",
                str(run.total_bases),
                f"{candidate.estimated_coverage:.1f}",
            ]) + "\n")
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Shortlist closed genomes with matching Illumina reads "
                    "for the PICOTA benchmark (docs/ROADMAP.md phase 0.5).")
    parser.add_argument("--email", required=True,
                        help="Contact address; NCBI requires one for E-utilities.")
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"),
                        help="NCBI API key, raising the rate limit to 10 req/s.")
    parser.add_argument("--species", nargs="+", default=DEFAULT_SPECIES,
                        help="Species to search. Default: %(default)s")
    parser.add_argument("--limit", type=int, default=100,
                        help="Max assemblies to inspect per species.")
    parser.add_argument("--min-coverage", type=float, default=40.0,
                        help="Reject runs below this estimated depth.")
    parser.add_argument("--broad", action="store_true",
                        help="Search all closed genomes, not only NCBI reference genomes.")
    parser.add_argument("--out", default="benchmark_candidates.tsv",
                        help="Output TSV path.")
    args = parser.parse_args(argv)

    client = EntrezClient(args.email, args.api_key)
    candidates = collect_candidates(
        client, args.species, args.limit, args.min_coverage,
        broad=args.broad,
        progress=lambda message: print(message, file=sys.stderr, flush=True))

    if not candidates:
        print("No candidates found. Try --broad, a lower --min-coverage, "
              "or a higher --limit.", file=sys.stderr)
        return 1

    write_candidates(args.out, candidates)
    print(f"\n{len(candidates)} candidate(s) written to {args.out}", file=sys.stderr)
    print("Next: stage 2 screens these for multi-copy IS families with "
          "differing cargo (docs/VALIDATION.md 3.2).", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
