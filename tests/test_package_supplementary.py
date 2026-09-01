"""
Tests for scripts/package_supplementary.py.

The package is what a reader receives, so the invariant that matters is that it
is self-describing and complete: every case indexed, every ground-truth file
carried across, and the parameters recorded so the case can be rebuilt.
"""

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from package_supplementary import (  # noqa: E402
    MANIFEST_COLUMNS,
    case_dirs,
    main,
    summarise,
)


def make_case(root, name, n_cts=2, novel=0, shared=2, backbone="NC_000913.3", seed=1):
    case = os.path.join(root, name)
    os.makedirs(case, exist_ok=True)
    cts = []
    for i in range(n_cts):
        cts.append({
            "CT_ID": "CT%03d" % (i + 1),
            "IS_Name": "ISfoo%d" % i,
            "Cargo_Type": "novel" if i < novel else "AMR",
            "CT_Length": 4000,
            "Sequence": "ACGT" * 1000,
        })
    payload = {
        "parameters": {"seed": seed, "shared_is": shared, "cargo_is_mode": "none",
                       "n_cts": n_cts},
        "backbone": backbone,
        "genome_length": 500000,
        "composite_transposons": cts,
    }
    with open(os.path.join(case, "ground_truth.json"), "w") as handle:
        json.dump(payload, handle)
    with open(os.path.join(case, "genome.fasta"), "w") as handle:
        handle.write(">chr\n" + "ACGT" * 100 + "\n")
    with open(os.path.join(case, "ground_truth.tsv"), "w") as handle:
        handle.write("CT_ID\tIS_Name\n"
                     + "".join("%s\t%s\n" % (c["CT_ID"], c["IS_Name"]) for c in cts))
    with open(os.path.join(case, "ground_truth_cts.fasta"), "w") as handle:
        for c in cts:
            handle.write(">%s\n%s\n" % (c["CT_ID"], c["Sequence"]))
    return case


class TestCaseDiscovery:
    def test_finds_cases_with_ground_truth(self, tmp_path):
        run = tmp_path / "run"
        make_case(str(run), "case_a")
        make_case(str(run), "case_b")
        assert len(list(case_dirs([str(run)]))) == 2

    def test_ignores_directories_without_ground_truth(self, tmp_path):
        run = tmp_path / "run"
        make_case(str(run), "case_a")
        os.makedirs(str(run / "not_a_case"))
        assert len(list(case_dirs([str(run)]))) == 1

    def test_reads_several_run_directories(self, tmp_path):
        make_case(str(tmp_path / "bench"), "c1")
        make_case(str(tmp_path / "scenarios"), "c2")
        found = list(case_dirs([str(tmp_path / "bench"), str(tmp_path / "scenarios")]))
        assert len(found) == 2


class TestSummarise:
    def test_counts_novel_cargo(self, tmp_path):
        case = make_case(str(tmp_path), "c", n_cts=4, novel=2)
        payload = json.load(open(os.path.join(case, "ground_truth.json")))
        assert summarise(payload)["NovelCargo"] == 2

    def test_reports_ct_count_and_backbone(self, tmp_path):
        case = make_case(str(tmp_path), "c", n_cts=3, backbone="NC_002695.2")
        payload = json.load(open(os.path.join(case, "ground_truth.json")))
        summary = summarise(payload)
        assert summary["NumCTs"] == 3
        assert summary["Backbone"] == "NC_002695.2"

    def test_missing_fields_do_not_raise(self):
        assert summarise({})["NumCTs"] == 0


class TestPackaging:
    @pytest.fixture
    def packaged(self, tmp_path):
        run = tmp_path / "bench"
        make_case(str(run), "case_a", n_cts=2, novel=1)
        make_case(str(run), "case_b", n_cts=3)
        out = tmp_path / "supp"
        main(["--runs", str(run), "--out", str(out)])
        return out

    def test_every_case_contributes_four_files(self, packaged):
        for label in ("bench_case_a", "bench_case_b"):
            for suffix in ("genome.fasta", "ground_truth.tsv",
                           "ground_truth.fasta", "parameters.json"):
                assert (packaged / ("%s_%s" % (label, suffix))).exists(), suffix

    def test_manifest_indexes_every_case(self, packaged):
        rows = [l.rstrip("\n").split("\t") for l in open(packaged / "manifest.tsv")]
        assert rows[0] == MANIFEST_COLUMNS
        assert len(rows) == 3
        assert {r[0] for r in rows[1:]} == {"bench_case_a", "bench_case_b"}

    def test_parameters_are_recoverable(self, packaged):
        params = json.load(open(packaged / "bench_case_a_parameters.json"))
        assert params["seed"] == 1
        assert params["n_cts"] == 2

    def test_readme_explains_the_coordinate_convention(self, packaged):
        """A reviewer will check this; the package has to state it."""
        readme = (packaged / "README.md").read_text()
        assert "CT_Start-1" in readme
        assert "1-based" in readme

    def test_readme_says_how_to_regenerate_reads(self, packaged):
        readme = (packaged / "README.md").read_text()
        assert "art_illumina" in readme and "spades.py" in readme

    def test_reads_and_assemblies_are_not_shipped(self, packaged):
        names = os.listdir(packaged)
        assert not any(n.endswith((".fq", ".fastq", ".gfa")) for n in names)
