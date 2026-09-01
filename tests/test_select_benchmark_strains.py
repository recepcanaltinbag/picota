"""
Tests for scripts/select_benchmark_strains.py (roadmap phase 0.5, stage 1).

Record shapes below mirror what NCBI actually returns -- they were taken from
live esummary responses for Escherichia coli assemblies and their SRA runs, so
the parsers are exercised against the real field layout rather than an idealised
one. No network access is needed.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from select_benchmark_strains import (  # noqa: E402
    OUTPUT_COLUMNS,
    build_assembly_query,
    build_biosample_query,
    build_assembly_query_reference_only,
    build_sra_query,
    collect_candidates,
    estimate_coverage,
    parse_assembly_summary,
    parse_biosample_summary,
    parse_sra_summary,
    pick_best_run,
    resolve_strain,
    write_candidates,
)

# Real Meta blob shape: assembly stats are an XML fragment inside a string field.
META_BLOB = ('<Stats><Stat category="total_length" sequence_tag="all">5147771</Stat>'
             '<Stat category="contig_count" sequence_tag="all">1</Stat></Stats>')


def assembly_record(**overrides):
    record = {
        "AssemblyAccession": "GCF_982530155.1",
        "SpeciesName": "Escherichia coli",
        "AssemblyStatus": "Complete Genome",
        "BioSampleAccn": "SAMEA122163316",
        "Biosource": {"InfraspeciesList": [{"Sub_type": "strain", "Sub_value": "ST131-A"}],
                      "Isolate": ""},
        "Meta": META_BLOB,
    }
    record.update(overrides)
    return record


def sra_record(acc="ERR17004303", total_bases="184339160", total_spots="634551",
               layout="PAIRED", platform="ILLUMINA", instrument="Illumina NovaSeq X"):
    layout_tag = f"<{layout}/>" if layout else ""
    return {
        "ExpXml": f'<Summary><Platform instrument_model="{instrument}">{platform}'
                  f'</Platform></Summary><Library_descriptor><LIBRARY_LAYOUT>'
                  f'{layout_tag}</LIBRARY_LAYOUT></Library_descriptor>',
        "Runs": f'<Run acc="{acc}" total_spots="{total_spots}" '
                f'total_bases="{total_bases}" is_public="true"/>',
    }


# Real BioSample SampleData shape: attributes live in an escaped XML blob.
def biosample_record(strain="JNEC-KY203", extra=""):
    attrs = ""
    if strain:
        attrs += f'&lt;Attribute attribute_name="strain"&gt;{strain}&lt;/Attribute&gt;'
    attrs += extra
    return {"SampleData": f"&lt;BioSample&gt;&lt;Attributes&gt;{attrs}"
                          f"&lt;/Attributes&gt;&lt;/BioSample&gt;"}


class TestQueries:
    def test_assembly_query_restricts_to_closed_current_genomes(self):
        query = build_assembly_query("Escherichia coli")
        assert '"Escherichia coli"[Organism]' in query
        assert '"complete genome"[Assembly Level]' in query
        assert '"latest refseq"[filter]' in query

    def test_default_query_is_not_restricted_to_reference_genomes(self):
        """
        The reference-genome filter cuts E. coli from 5350 closed assemblies to
        2, and the survivor is K-12 MG1655 -- no resistance CTs, so useless here.
        """
        assert '"reference genome"[filter]' not in build_assembly_query("Escherichia coli")

    def test_reference_only_query_adds_the_restriction(self):
        query = build_assembly_query_reference_only("Escherichia coli")
        assert '"reference genome"[filter]' in query
        assert '"complete genome"[Assembly Level]' in query

    def test_sra_query_pins_platform_strategy_and_layout(self):
        query = build_sra_query("SAMEA122163316")
        for token in ("SAMEA122163316[BioSample]", "illumina[Platform]",
                      "wgs[Strategy]", "paired[Layout]", "public[Access]"):
            assert token in query


class TestParseAssemblySummary:
    def test_parses_a_complete_genome(self):
        candidate = parse_assembly_summary(assembly_record())
        assert candidate.accession == "GCF_982530155.1"
        assert candidate.organism == "Escherichia coli"
        assert candidate.strain == "ST131-A"
        assert candidate.biosample == "SAMEA122163316"

    def test_genome_length_comes_from_the_meta_blob(self):
        """NCBI has no Meta_total_length field; the number is inside Meta."""
        assert parse_assembly_summary(assembly_record()).genome_length == 5147771

    def test_draft_assemblies_are_rejected(self):
        """A draft genome cannot serve as ground truth -- its CTs are unresolved."""
        assert parse_assembly_summary(assembly_record(AssemblyStatus="Scaffold")) is None
        assert parse_assembly_summary(assembly_record(AssemblyStatus="Contig")) is None

    def test_missing_biosample_is_rejected(self):
        """Without a BioSample there is no way to pin reads to this same isolate."""
        assert parse_assembly_summary(assembly_record(BioSampleAccn="")) is None

    def test_absent_strain_is_tolerated(self):
        """Many records carry an empty InfraspeciesList; that is not fatal."""
        record = assembly_record(Biosource={"InfraspeciesList": [], "Isolate": ""})
        assert parse_assembly_summary(record).strain == ""

    def test_strain_falls_back_to_isolate_field(self):
        record = assembly_record(
            Biosource={"InfraspeciesList": [], "Isolate": "strain=K-12;other=x"})
        assert parse_assembly_summary(record).strain == "K-12"

    def test_unparseable_meta_yields_zero_length(self):
        record = assembly_record(Meta="<Stats></Stats>")
        assert parse_assembly_summary(record).genome_length == 0


class TestParseSraSummary:
    def test_parses_platform_layout_and_run(self):
        runs = parse_sra_summary(sra_record())
        assert len(runs) == 1
        assert runs[0].accession == "ERR17004303"
        assert runs[0].platform == "ILLUMINA"
        assert runs[0].layout == "PAIRED"
        assert runs[0].total_bases == 184339160
        assert runs[0].total_spots == 634551

    def test_single_end_layout_detected(self):
        assert parse_sra_summary(sra_record(layout="SINGLE"))[0].layout == "SINGLE"

    def test_multiple_runs_in_one_experiment(self):
        record = sra_record()
        record["Runs"] += '<Run acc="ERR17004304" total_spots="10" total_bases="20"/>'
        assert [r.accession for r in parse_sra_summary(record)] == ["ERR17004303", "ERR17004304"]

    def test_run_without_accession_skipped(self):
        record = sra_record()
        record["Runs"] = '<Run total_spots="10" total_bases="20"/>'
        assert parse_sra_summary(record) == []

    def test_empty_record_is_not_an_error(self):
        assert parse_sra_summary({}) == []


class TestCoverageAndSelection:
    def test_estimate_coverage(self):
        assert estimate_coverage(184339160, 5147771) == pytest.approx(35.81, abs=0.01)

    def test_coverage_is_none_when_inputs_missing(self):
        assert estimate_coverage(0, 5147771) is None
        assert estimate_coverage(184339160, 0) is None

    def test_run_below_threshold_is_rejected(self):
        """The live ERR17004303 case: 35.8x against a 40x floor."""
        runs = parse_sra_summary(sra_record())
        assert pick_best_run(runs, 5147771, 40.0) == (None, None)

    def test_run_above_threshold_is_accepted(self):
        runs = parse_sra_summary(sra_record())
        run, coverage = pick_best_run(runs, 5147771, 30.0)
        assert run.accession == "ERR17004303"
        assert coverage == pytest.approx(35.81, abs=0.01)

    def test_single_end_runs_are_rejected(self):
        runs = parse_sra_summary(sra_record(total_bases="900000000", layout="SINGLE"))
        assert pick_best_run(runs, 5147771, 30.0) == (None, None)

    def test_smallest_sufficient_run_wins(self):
        """Excess depth costs download and assembly time without adding signal."""
        runs = (parse_sra_summary(sra_record(acc="ERR1", total_bases="900000000"))
                + parse_sra_summary(sra_record(acc="ERR2", total_bases="300000000")))
        run, _ = pick_best_run(runs, 5147771, 30.0)
        assert run.accession == "ERR2"


class FakeClient:
    """Stands in for EntrezClient; records the queries it was asked to run."""

    def __init__(self, assembly_ids, assembly_records, sra_ids, sra_records,
                 biosample_records=None):
        self.assembly_ids = assembly_ids
        self.assembly_records = assembly_records
        self.sra_ids = sra_ids
        self.sra_records = sra_records
        self.biosample_records = biosample_records or []
        self.queries = []

    def esearch(self, db, term, retmax):
        self.queries.append((db, term))
        if db == "assembly":
            return self.assembly_ids
        if db == "biosample":
            return ["9"] if self.biosample_records else []
        return self.sra_ids

    def esummary(self, db, ids):
        if db == "assembly":
            return self.assembly_records
        if db == "biosample":
            return self.biosample_records
        return self.sra_records


class TestCollectCandidates:
    def _client(self, **kwargs):
        defaults = dict(assembly_ids=["1"], assembly_records=[assembly_record()],
                        sra_ids=["2"], sra_records=[sra_record()])
        defaults.update(kwargs)
        return FakeClient(**defaults)

    def test_returns_one_candidate_per_usable_genome(self):
        found = collect_candidates(self._client(), ["Escherichia coli"], 10, 30.0)
        assert len(found) == 1
        assert found[0].assembly.accession == "GCF_982530155.1"
        assert found[0].run.accession == "ERR17004303"
        assert found[0].estimated_coverage == pytest.approx(35.81, abs=0.01)

    def test_genome_without_matching_reads_is_skipped(self):
        assert collect_candidates(self._client(sra_ids=[]), ["Escherichia coli"], 10, 30.0) == []

    def test_genome_whose_reads_are_too_shallow_is_skipped(self):
        assert collect_candidates(self._client(), ["Escherichia coli"], 10, 100.0) == []

    def test_draft_genomes_never_reach_the_sra_lookup(self):
        client = self._client(assembly_records=[assembly_record(AssemblyStatus="Contig")])
        assert collect_candidates(client, ["Escherichia coli"], 10, 30.0) == []
        assert [db for db, _ in client.queries] == ["assembly"]

    def test_each_species_is_searched(self):
        client = self._client()
        collect_candidates(client, ["Escherichia coli", "Klebsiella pneumoniae"], 10, 30.0)
        assembly_terms = [term for db, term in client.queries if db == "assembly"]
        assert len(assembly_terms) == 2
        assert any("Klebsiella pneumoniae" in term for term in assembly_terms)


class TestWriteCandidates:
    def test_output_has_the_documented_columns(self, tmp_path):
        found = collect_candidates(
            FakeClient(["1"], [assembly_record()], ["2"], [sra_record()]),
            ["Escherichia coli"], 10, 30.0)
        out = write_candidates(str(tmp_path / "c.tsv"), found)
        rows = [line.rstrip("\n").split("\t") for line in open(out)]
        assert rows[0] == OUTPUT_COLUMNS
        assert rows[1][rows[0].index("RunAccession")] == "ERR17004303"
        assert rows[1][rows[0].index("EstimatedCoverage")] == "35.8"
        assert rows[1][rows[0].index("GenomeLength")] == "5147771"

    def test_missing_strain_written_as_na(self, tmp_path):
        record = assembly_record(Biosource={"InfraspeciesList": [], "Isolate": ""})
        found = collect_candidates(FakeClient(["1"], [record], ["2"], [sra_record()]),
                                   ["Escherichia coli"], 10, 30.0)
        out = write_candidates(str(tmp_path / "c.tsv"), found)
        rows = [line.rstrip("\n").split("\t") for line in open(out)]
        assert rows[1][rows[0].index("Strain")] == "NA"

    def test_empty_input_writes_header_only(self, tmp_path):
        out = write_candidates(str(tmp_path / "c.tsv"), [])
        assert open(out).read().strip() == "\t".join(OUTPUT_COLUMNS)


class TestInstrument:
    def test_instrument_model_captured_separately_from_platform(self):
        """Read length drives repeat resolution, so the model matters, not just ILLUMINA."""
        run = parse_sra_summary(sra_record())[0]
        assert run.platform == "ILLUMINA"
        assert run.instrument == "Illumina NovaSeq X"

    def test_platform_falls_back_to_instrument_when_absent(self):
        record = sra_record()
        record["ExpXml"] = record["ExpXml"].replace(">ILLUMINA<", "><")
        assert parse_sra_summary(record)[0].platform == "Illumina NovaSeq X"

    def test_missing_instrument_is_empty_not_an_error(self):
        record = sra_record()
        record["ExpXml"] = record["ExpXml"].replace(
            ' instrument_model="Illumina NovaSeq X"', "")
        assert parse_sra_summary(record)[0].instrument == ""


class TestBiosampleStrain:
    def test_strain_extracted_from_sample_data_blob(self):
        assert parse_biosample_summary(biosample_record()) == "JNEC-KY203"

    def test_isolate_used_when_strain_absent(self):
        record = biosample_record(
            strain="",
            extra='&lt;Attribute attribute_name="isolate"&gt;ST131-A&lt;/Attribute&gt;')
        assert parse_biosample_summary(record) == "ST131-A"

    def test_strain_preferred_over_isolate(self):
        record = biosample_record(
            extra='&lt;Attribute attribute_name="isolate"&gt;other&lt;/Attribute&gt;')
        assert parse_biosample_summary(record) == "JNEC-KY203"

    def test_no_strain_attribute_returns_empty(self):
        assert parse_biosample_summary(biosample_record(strain="")) == ""

    def test_empty_record_returns_empty(self):
        assert parse_biosample_summary({}) == ""

    def test_biosample_query_targets_the_accession(self):
        assert build_biosample_query("SAMD01706138") == "SAMD01706138[Accession]"

    def test_resolve_strain_round_trip(self):
        client = FakeClient(["1"], [assembly_record()], ["2"], [sra_record()],
                            biosample_records=[biosample_record()])
        assert resolve_strain(client, "SAMD01706138") == "JNEC-KY203"

    def test_resolve_strain_without_match_is_empty(self):
        client = FakeClient(["1"], [assembly_record()], ["2"], [sra_record()])
        assert resolve_strain(client, "SAMD01706138") == ""


class TestResolveStrainsFlag:
    def _client(self):
        record = assembly_record(Biosource={"InfraspeciesList": [], "Isolate": ""})
        return FakeClient(["1"], [record], ["2"], [sra_record()],
                          biosample_records=[biosample_record()])

    def test_strain_left_empty_by_default(self):
        """The default path must not pay for two extra requests per candidate."""
        client = self._client()
        found = collect_candidates(client, ["Escherichia coli"], 10, 30.0)
        assert found[0].assembly.strain == ""
        assert "biosample" not in [db for db, _ in client.queries]

    def test_strain_filled_when_requested(self):
        client = self._client()
        found = collect_candidates(client, ["Escherichia coli"], 10, 30.0,
                                   resolve_strains=True)
        assert found[0].assembly.strain == "JNEC-KY203"

    def test_existing_strain_is_not_re_fetched(self):
        client = FakeClient(["1"], [assembly_record()], ["2"], [sra_record()],
                            biosample_records=[biosample_record()])
        found = collect_candidates(client, ["Escherichia coli"], 10, 30.0,
                                   resolve_strains=True)
        assert found[0].assembly.strain == "ST131-A"
        assert "biosample" not in [db for db, _ in client.queries]

    def test_instrument_reaches_the_output(self, tmp_path):
        found = collect_candidates(self._client(), ["Escherichia coli"], 10, 30.0)
        out = write_candidates(str(tmp_path / "c.tsv"), found)
        rows = [line.rstrip("\n").split("\t") for line in open(out)]
        assert rows[1][rows[0].index("Instrument")] == "Illumina NovaSeq X"
