"""
test_e2e_gfa.py – Real end-to-end integration test for PICOTA.

Uses the bundled testNitro.gfa assembly graph to exercise the full
post-assembly pipeline (cycle detection → BLAST scoring → enriched CSV)
without needing an SRA download or fresh genome assembly.

Run with the evobiomig conda environment:

    conda run -n evobiomig python3 -m pytest tests/test_e2e_gfa.py -v

The test writes all output to a pytest tmp_path directory and cleans
up automatically on success.
"""

import csv
import os
import shutil
import subprocess
import sys

import pytest

# ─── Path setup (mirrors conftest.py) ───────────────────────────────────────
HERE         = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)                        # picota/
INNER_ROOT   = os.path.join(PROJECT_ROOT, 'picota')         # picota/picota/
TEST_DATA    = os.path.join(INNER_ROOT, 'test_data')
GFA_FILE     = os.path.join(TEST_DATA, 'testNitro.gfa')

DB_BASE = os.path.join(INNER_ROOT, 'DBs')


def db_path(*parts):
    return os.path.join(DB_BASE, *parts)


def tool_available(name: str) -> bool:
    return shutil.which(name) is not None


# ─── Skip guard ──────────────────────────────────────────────────────────────
REQUIRED_TOOLS = ['prodigal', 'blastn', 'blastp', 'blastx', 'makeblastdb']
MISSING_TOOLS  = [t for t in REQUIRED_TOOLS if not tool_available(t)]

skip_no_tools = pytest.mark.skipif(
    bool(MISSING_TOOLS),
    reason=f"Required tools missing: {', '.join(MISSING_TOOLS)}"
)

skip_no_dbs = pytest.mark.skipif(
    not os.path.exists(db_path('Antibiotics', 'protein_fasta_protein_homolog_model.fasta')),
    reason="BLAST databases not found — set up picota/picota/DBs/ first"
)


# ─── Helper: count FASTA records ─────────────────────────────────────────────
def count_fasta(path: str) -> int:
    if not os.path.exists(path):
        return 0
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith('>'))


# ═══════════════════════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestE2E_CycleDetection:
    """Step 1 – cycle detection directly from testNitro.gfa."""

    def test_gfa_file_present(self):
        assert os.path.exists(GFA_FILE), f"Test GFA not found: {GFA_FILE}"

    def test_cycle_analysis_produces_fasta(self, tmp_path):
        from src.cycle_finderv2 import cycle_analysis

        out_fasta = str(tmp_path / 'cycles.fasta')
        cycle_analysis(
            GFA_FILE, out_fasta,
            find_all_path=False, path_limit=15,
            min_size_of_cycle=3000, max_size_of_cycle=100_000,
            name_prefix_cycle='Cycle',
            min_component_number=1, max_component_number=25,
            k_mer_sim=200, threshold_sim=99
        )
        assert os.path.exists(out_fasta), "cycle_analysis produced no output FASTA"
        n = count_fasta(out_fasta)
        assert n > 0, "No cycles detected in testNitro.gfa"

    def test_cycle_fasta_has_length_in_header(self, tmp_path):
        from src.cycle_finderv2 import cycle_analysis

        out_fasta = str(tmp_path / 'cycles.fasta')
        cycle_analysis(
            GFA_FILE, out_fasta,
            find_all_path=False, path_limit=15,
            min_size_of_cycle=3000, max_size_of_cycle=100_000,
            name_prefix_cycle='Cycle',
            min_component_number=1, max_component_number=25,
            k_mer_sim=200, threshold_sim=99
        )
        with open(out_fasta) as fh:
            headers = [l.strip() for l in fh if l.startswith('>')]
        # Every header should look like ">Cycle_N-lenXXX-compY-"
        import re
        for h in headers:
            assert re.search(r'-len\d+-', h), f"Header missing length: {h}"


class TestE2E_Scoring:
    """Step 2 – BLAST scoring on real cycle FASTA."""

    @pytest.fixture(scope='class')
    def cycle_fasta(self, tmp_path_factory):
        """Generate cycle FASTA once for the whole class."""
        from src.cycle_finderv2 import cycle_analysis

        out = str(tmp_path_factory.mktemp('cycles') / 'cycles.fasta')
        cycle_analysis(
            GFA_FILE, out,
            find_all_path=False, path_limit=15,
            min_size_of_cycle=3000, max_size_of_cycle=100_000,
            name_prefix_cycle='Cycle',
            min_component_number=1, max_component_number=25,
            k_mer_sim=200, threshold_sim=99
        )
        return out

    @skip_no_tools
    @skip_no_dbs
    def test_scoring_creates_final_tab(self, tmp_path, cycle_fasta):
        from src.scoringv4ProtBlast import scoring_main

        if not os.path.exists(cycle_fasta) or os.path.getsize(cycle_fasta) == 0:
            pytest.skip("Cycle FASTA empty")

        scoring_main(
            cycle_fasta, str(tmp_path),
            db_path('Antibiotics', 'protein_fasta_protein_homolog_model.fasta'),
            db_path('Xenobiotics', 'Xenobiotics_classified.fasta'),
            db_path('ISes', '_tncentral_nointegrall_isfinder-TNs.fasta'),
            db_path('CompTns', 'Known_Tns.fasta'),
            str(tmp_path),
            mean_of_CompTns=5850, std_of_CompTns=2586,
            total_score_type=0, threshold_final_score=0,
            max_z=20, dist_type=1,
            path_of_prodigal='prodigal',
            path_of_blastn='blastn',
            path_of_makeblastdb='makeblastdb',
            path_of_blastx='blastx',
            path_of_blastp='blastp',
            logger_name='picota_e2e'
        )
        tab = str(tmp_path / 'picota_final_tab')
        assert os.path.exists(tab), "picota_final_tab not created"
        assert os.path.getsize(tab) > 0, "picota_final_tab is empty"

    @skip_no_tools
    @skip_no_dbs
    def test_scoring_tab_has_header(self, tmp_path, cycle_fasta):
        from src.scoringv4ProtBlast import scoring_main

        if not os.path.exists(cycle_fasta) or os.path.getsize(cycle_fasta) == 0:
            pytest.skip("Cycle FASTA empty")

        scoring_main(
            cycle_fasta, str(tmp_path),
            db_path('Antibiotics', 'protein_fasta_protein_homolog_model.fasta'),
            db_path('Xenobiotics', 'Xenobiotics_classified.fasta'),
            db_path('ISes', '_tncentral_nointegrall_isfinder-TNs.fasta'),
            db_path('CompTns', 'Known_Tns.fasta'),
            str(tmp_path),
            mean_of_CompTns=5850, std_of_CompTns=2586,
            total_score_type=0, threshold_final_score=0,
            max_z=20, dist_type=1,
            path_of_prodigal='prodigal',
            path_of_blastn='blastn',
            path_of_makeblastdb='makeblastdb',
            path_of_blastx='blastx',
            path_of_blastp='blastp',
            logger_name='picota_e2e'
        )
        tab = str(tmp_path / 'picota_final_tab')
        with open(tab) as fh:
            header = fh.readline().strip().split('\t')
        expected = ['CycleID', 'SRAID', 'kmer', 'score0', 'score1', 'score2',
                    'NumIS', 'ISproducts', 'IScoords']
        for col in expected:
            assert col in header, f"Missing column '{col}' in picota_final_tab"


class TestE2E_EnrichedCSV:
    """Step 3 – enriched CSV generated from picota_final_tab."""

    @pytest.fixture(scope='class')
    def scored_outdir(self, tmp_path_factory):
        """Run full scoring and return the output directory."""
        from src.cycle_finderv2 import cycle_analysis
        from src.scoringv4ProtBlast import scoring_main

        if MISSING_TOOLS:
            pytest.skip(f"Tools missing: {', '.join(MISSING_TOOLS)}")
        if not os.path.exists(db_path('Antibiotics', 'protein_fasta_protein_homolog_model.fasta')):
            pytest.skip("BLAST databases not found")

        out_dir   = tmp_path_factory.mktemp('scoring')
        cyc_dir   = tmp_path_factory.mktemp('cycles')
        cycle_fasta = str(cyc_dir / 'cycles.fasta')

        cycle_analysis(
            GFA_FILE, cycle_fasta,
            find_all_path=False, path_limit=15,
            min_size_of_cycle=3000, max_size_of_cycle=100_000,
            name_prefix_cycle='Cycle',
            min_component_number=1, max_component_number=25,
            k_mer_sim=200, threshold_sim=99
        )

        scoring_main(
            cycle_fasta, str(out_dir),
            db_path('Antibiotics', 'protein_fasta_protein_homolog_model.fasta'),
            db_path('Xenobiotics', 'Xenobiotics_classified.fasta'),
            db_path('ISes', '_tncentral_nointegrall_isfinder-TNs.fasta'),
            db_path('CompTns', 'Known_Tns.fasta'),
            str(out_dir),
            mean_of_CompTns=5850, std_of_CompTns=2586,
            total_score_type=0, threshold_final_score=0,
            max_z=20, dist_type=1,
            path_of_prodigal='prodigal',
            path_of_blastn='blastn',
            path_of_makeblastdb='makeblastdb',
            path_of_blastx='blastx',
            path_of_blastp='blastp',
            logger_name='picota_e2e'
        )
        return str(out_dir)

    def test_enriched_csv_created(self, scored_outdir):
        csv_path = os.path.join(scored_outdir, 'picota_enriched.csv')
        assert os.path.exists(csv_path), "picota_enriched.csv was not generated"

    def test_enriched_csv_has_expected_columns(self, scored_outdir):
        csv_path = os.path.join(scored_outdir, 'picota_enriched.csv')
        if not os.path.exists(csv_path):
            pytest.skip("picota_enriched.csv not found")
        with open(csv_path, newline='') as fh:
            reader = csv.DictReader(fh)
            header = reader.fieldnames or []
        required = ['CT_Tag', 'Category', 'CycleID', 'CT_Length_bp',
                    'IS_Group', 'IS_Family', 'Antibiotic_Class', 'Resistance_Gene', 'Score']
        for col in required:
            assert col in header, f"Missing column '{col}' in picota_enriched.csv"

    def test_enriched_csv_ct_tags_sequential(self, scored_outdir):
        csv_path = os.path.join(scored_outdir, 'picota_enriched.csv')
        if not os.path.exists(csv_path):
            pytest.skip("picota_enriched.csv not found")
        with open(csv_path, newline='') as fh:
            reader = csv.DictReader(fh)
            tags = [row['CT_Tag'] for row in reader]
        assert len(tags) > 0, "No rows in enriched CSV"
        # Tags must all match CT\d{3}
        import re
        for tag in tags:
            assert re.match(r'^CT\d{3,}$', tag), f"Unexpected CT tag format: {tag}"

    def test_enriched_csv_category_valid(self, scored_outdir):
        csv_path = os.path.join(scored_outdir, 'picota_enriched.csv')
        if not os.path.exists(csv_path):
            pytest.skip("picota_enriched.csv not found")
        with open(csv_path, newline='') as fh:
            reader = csv.DictReader(fh)
            categories = {row['Category'] for row in reader}
        assert categories.issubset({'Novel', 'Known'}), \
            f"Unexpected Category values: {categories - {'Novel', 'Known'}}"

    def test_enriched_csv_ct_length_positive(self, scored_outdir):
        csv_path = os.path.join(scored_outdir, 'picota_enriched.csv')
        if not os.path.exists(csv_path):
            pytest.skip("picota_enriched.csv not found")
        with open(csv_path, newline='') as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)
        for row in rows:
            assert int(row['CT_Length_bp']) >= 0, \
                f"Negative CT_Length_bp for {row['CT_Tag']}"


class TestE2E_OutputFormatter:
    """Unit tests for output_formatter.py (no BLAST needed)."""

    def _make_tab(self, tmp_path, content: str) -> str:
        p = str(tmp_path / 'picota_final_tab')
        with open(p, 'w') as fh:
            fh.write(content)
        return p

    def test_empty_tab_returns_empty(self, tmp_path):
        from src.output_formatter import build_enriched_rows

        tab = self._make_tab(tmp_path,
            'CycleID\tSRAID\tkmer\tscore0\tscore1\tscore2\t'
            'NumIS\tISproducts\tIScoords\t'
            'NumAnt\tAntproducts\tAntcoords\t'
            'NumXeno\tXenoproducts\tXenocoords\t'
            'NumCompTN\tCompTN\tCompTNscoords\n'
        )
        rows = build_enriched_rows(tab)
        assert rows == []

    def test_single_novel_ct(self, tmp_path):
        from src.output_formatter import build_enriched_rows

        tab = self._make_tab(tmp_path,
            'CycleID\tSRAID\tkmer\tscore0\tscore1\tscore2\t'
            'NumIS\tISproducts\tIScoords\t'
            'NumAnt\tAntproducts\tAntcoords\t'
            'NumXeno\tXenoproducts\tXenocoords\t'
            'NumCompTN\tCompTN\tCompTNscoords\n'
            'Cycle_1-len8500-comp4-\tSRR123\t99\t85.0\t80.0\t70.0\t'
            '2\tIS26;IS26\t100-1500;7000-8500\t'
            '1\tblaKPC-2\t2000-2800\t'
            '0\t\t\t'
            '0\tNovel\t\n'
        )
        rows = build_enriched_rows(tab)
        assert len(rows) == 1
        assert rows[0]['CT_Tag'] == 'CT001'
        assert rows[0]['Category'] == 'Novel'
        assert rows[0]['CT_Length_bp'] == 8500
        assert rows[0]['IS_Family'] == 'IS26'
        assert rows[0]['Antibiotic_Class'] == 'Beta-lactam'

    def test_known_ct_has_category_known(self, tmp_path):
        from src.output_formatter import build_enriched_rows

        tab = self._make_tab(tmp_path,
            'CycleID\tSRAID\tkmer\tscore0\tscore1\tscore2\t'
            'NumIS\tISproducts\tIScoords\t'
            'NumAnt\tAntproducts\tAntcoords\t'
            'NumXeno\tXenoproducts\tXenocoords\t'
            'NumCompTN\tCompTN\tCompTNscoords\n'
            'Cycle_2-len6000-comp3-\tSRR456\t99\t75.0\t60.0\t55.0\t'
            '2\tIS1;IS1\t50-1150;4800-5900\t'
            '0\t\t\t'
            '0\t\t\t'
            '1\tTn3\t0-6000\n'
        )
        rows = build_enriched_rows(tab)
        assert rows[0]['Category'] == 'Known'
        assert rows[0]['Known_CompTN'] == 'Tn3'

    def test_multiple_amr_classes_expand_rows(self, tmp_path):
        from src.output_formatter import build_enriched_rows

        tab = self._make_tab(tmp_path,
            'CycleID\tSRAID\tkmer\tscore0\tscore1\tscore2\t'
            'NumIS\tISproducts\tIScoords\t'
            'NumAnt\tAntproducts\tAntcoords\t'
            'NumXeno\tXenoproducts\tXenocoords\t'
            'NumCompTN\tCompTN\tCompTNscoords\n'
            'Cycle_3-len10000-comp5-\tSRR789\t99\t90.0\t88.0\t80.0\t'
            '2\tIS26;IS26\t100-1500;8500-10000\t'
            '2\tblaKPC-2;tetA\t2000-2800;4000-4800\t'
            '0\t\t\t'
            '0\tNovel\t\n'
        )
        rows = build_enriched_rows(tab)
        # Two distinct AMR classes → two rows for the same CT
        assert len(rows) == 2
        tags = {r['CT_Tag'] for r in rows}
        assert tags == {'CT001'}
        classes = {r['Antibiotic_Class'] for r in rows}
        assert 'Beta-lactam' in classes
        assert 'Tetracycline' in classes

    def test_write_enriched_csv_creates_file(self, tmp_path):
        from src.output_formatter import write_enriched_csv

        tab = self._make_tab(tmp_path,
            'CycleID\tSRAID\tkmer\tscore0\tscore1\tscore2\t'
            'NumIS\tISproducts\tIScoords\t'
            'NumAnt\tAntproducts\tAntcoords\t'
            'NumXeno\tXenoproducts\tXenocoords\t'
            'NumCompTN\tCompTN\tCompTNscoords\n'
            'Cycle_1-len5000-comp2-\tSRR001\t99\t70.0\t65.0\t60.0\t'
            '2\tIS26;IS26\t50-1500;3500-5000\t'
            '1\tblaOXA-48\t1600-2400\t'
            '0\t\t\t'
            '0\tNovel\t\n'
        )
        out_csv = str(tmp_path / 'enriched.csv')
        n = write_enriched_csv(tab, out_csv)
        assert n == 1
        assert os.path.exists(out_csv)
        with open(out_csv, newline='') as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 1
        assert rows[0]['CT_Tag'] == 'CT001'

    def test_is_family_inference_is26(self):
        from src.output_formatter import infer_is_family
        group, family = infer_is_family('IS26')
        assert group == 'IS6'
        assert family == 'IS26'

    def test_is_family_inference_isecp(self):
        # ISFinder: ISEcp1 is in IS1380 superfamily
        from src.output_formatter import infer_is_family
        group, family = infer_is_family('ISEcp1')
        assert group == 'IS1380'

    def test_is_family_inference_is10(self):
        # IS10 is in IS4 superfamily, NOT IS1
        from src.output_formatter import infer_is_family
        group, family = infer_is_family('IS10')
        assert group == 'IS4'
        assert family == 'IS10'

    def test_is_family_inference_is10_variant(self):
        from src.output_formatter import infer_is_family
        group, family = infer_is_family('IS10_B')
        assert group == 'IS4'

    def test_is_family_inference_is21(self):
        # IS21 is its own superfamily, NOT IS3
        from src.output_formatter import infer_is_family
        group, family = infer_is_family('IS21')
        assert group == 'IS21'

    def test_is_family_inference_is26_variant(self):
        from src.output_formatter import infer_is_family
        group, family = infer_is_family('IS26_1')
        assert group == 'IS6'

    def test_is_family_inference_is2_superfamily(self):
        # IS2 element belongs to IS3 superfamily
        from src.output_formatter import infer_is_family
        group, family = infer_is_family('IS2')
        assert group == 'IS3'

    def test_is_family_inference_iscr(self):
        from src.output_formatter import infer_is_family
        group, family = infer_is_family('ISCR3')
        assert group == 'IS91'

    def test_is_family_unknown(self):
        from src.output_formatter import infer_is_family
        group, family = infer_is_family('XXXXXXXX')
        assert group == 'Unknown'

    def test_antibiotic_class_bla(self):
        from src.output_formatter import infer_antibiotic_class
        assert infer_antibiotic_class('blaKPC-2') == 'Beta-lactam'

    def test_antibiotic_class_tet(self):
        from src.output_formatter import infer_antibiotic_class
        assert infer_antibiotic_class('tetA') == 'Tetracycline'

    def test_antibiotic_class_erm(self):
        from src.output_formatter import infer_antibiotic_class
        assert infer_antibiotic_class('erm(B)') == 'Macrolide-Lincosamide-Streptogramin'

    def test_antibiotic_class_unknown(self):
        from src.output_formatter import infer_antibiotic_class
        assert infer_antibiotic_class('hypothetical_protein') == 'Other'
