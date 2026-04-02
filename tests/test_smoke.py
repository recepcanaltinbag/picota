"""
PICOTA Pipeline Smoke Tests
============================
Her pipeline adımının çalışıp çalışmadığını hızlıca kontrol eder.
Assembly ve SRA download adımları atlanır (external tool / ağ bağımlılığı).

Çalıştırmak için:
    cd picota/
    python3 -m pytest tests/test_smoke.py -v

Renkli özet:
    python3 -m pytest tests/test_smoke.py -v --tb=short
"""

import sys
import os
import shutil
import subprocess
import tempfile
import textwrap

import pytest

# ─── path setup ──────────────────────────────────────────────────────────────
PICOTA_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # picota/
SRC_ROOT    = os.path.join(PICOTA_ROOT, 'picota')                           # picota/picota/
sys.path.insert(0, SRC_ROOT)

TEST_DATA   = os.path.join(SRC_ROOT, 'test_data')
GFA_FILE    = os.path.join(TEST_DATA, 'testNitro.gfa')
DB_ROOT     = os.path.join(SRC_ROOT, 'DBs')

# ─── helpers ─────────────────────────────────────────────────────────────────

def tool_available(name):
    return shutil.which(name) is not None


def skip_if_missing_tool(*tools):
    missing = [t for t in tools if not tool_available(t)]
    if missing:
        return pytest.mark.skipif(True, reason=f"Tool(s) eksik: {', '.join(missing)}")
    return pytest.mark.skipif(False, reason="")


def db_path(relative):
    return os.path.join(DB_ROOT, relative)


# ─────────────────────────────────────────────────────────────────────────────
# ADIM 0 — Import kontrolü
# ─────────────────────────────────────────────────────────────────────────────

class TestStep0_Imports:
    """Tüm modüller import edilebiliyor mu?"""

    def test_import_cycle_finderv2(self):
        from src.cycle_finderv2 import GraphWork, Cycle, cycle_analysis
        assert GraphWork is not None

    def test_import_scoringv4(self):
        from src.scoringv4ProtBlast import scoring_main, calculate_total_score, merge_intervals
        assert scoring_main is not None

    def test_import_cycle_kmer_hash(self):
        from src.cycle_kmer_hash import filter_cycles_with_kmer, get_kmer_hashes
        assert filter_cycles_with_kmer is not None

    def test_import_config_loader(self):
        from src.config_loader import load_config, Config
        assert Config is not None

    def test_import_split_cycle(self):
        from src.split_cycle_coords_for_is import split_cycles_from_picota
        assert split_cycles_from_picota is not None

    def test_import_logger_setup(self):
        from src.logger_setup import setup_logger_from_config
        assert setup_logger_from_config is not None

    def test_import_assembly(self):
        from src.assembly import assembly_main
        assert assembly_main is not None

    def test_import_sra_download(self):
        from src.sra_download import run_sra_down
        assert run_sra_down is not None


# ─────────────────────────────────────────────────────────────────────────────
# ADIM 1 — Config yükleme
# ─────────────────────────────────────────────────────────────────────────────

class TestStep1_Config:
    """config.yaml doğru parse ediliyor mu?"""

    CONFIG_PATH = os.path.join(SRC_ROOT, '..', 'picota', 'config.yaml')

    def test_config_loads(self):
        from src.config_loader import load_config
        cfg = load_config(os.path.join(SRC_ROOT, 'config.yaml'))
        assert cfg is not None

    def test_config_has_required_fields(self):
        from src.config_loader import load_config
        cfg = load_config(os.path.join(SRC_ROOT, 'config.yaml'))
        assert cfg.paths.path_to_antibiotics != ''
        assert cfg.options.min_size_of_cycle > 0
        assert cfg.options.max_size_of_cycle > cfg.options.min_size_of_cycle
        assert cfg.options.k_mer_sim > 0

    def test_config_tolerances(self):
        from src.config_loader import load_config
        cfg = load_config(os.path.join(SRC_ROOT, 'config.yaml'))
        assert cfg.tolerances.tolerance > 0
        assert cfg.tolerances.max_reasonable_gap > 0

    def test_config_patterns_defined(self):
        from src.config_loader import load_config
        cfg = load_config(os.path.join(SRC_ROOT, 'config.yaml'))
        assert 'TCT' in cfg.patterns   # transposon-cargo-transposon temel yapı
        assert len(cfg.patterns) >= 3


# ─────────────────────────────────────────────────────────────────────────────
# ADIM 2 — GFA parse + graph oluşturma
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.skipif(not os.path.exists(GFA_FILE), reason="testNitro.gfa bulunamadı")
class TestStep2_GfaParse:
    """GFA dosyası okunuyor ve graph oluşturuluyor mu?"""

    def test_parse_produces_nodes_and_edges(self):
        from src.cycle_finderv2 import GraphWork
        gw = GraphWork()
        nodes, edges = gw.parse_gfa(GFA_FILE)
        assert len(nodes) > 0, "Node bulunamadı"
        assert len(edges) > 0, "Edge bulunamadı"

    def test_graph_generated(self):
        from src.cycle_finderv2 import GraphWork
        gw = GraphWork()
        nodes, edges = gw.parse_gfa(GFA_FILE)
        G = gw.generate_genome_graph(nodes, edges)
        assert G is not None
        assert len(G.adj) > 0

    def test_both_strands_present(self):
        from src.cycle_finderv2 import GraphWork
        gw = GraphWork()
        nodes, _ = gw.parse_gfa(GFA_FILE)
        sample_id = list(nodes.keys())[0][:-1]
        assert sample_id + '+' in nodes
        assert sample_id + '-' in nodes

    def test_sequences_not_empty(self):
        from src.cycle_finderv2 import GraphWork
        gw = GraphWork()
        nodes, _ = gw.parse_gfa(GFA_FILE)
        for key, val in list(nodes.items())[:10]:
            assert len(val['Sequence']) > 0, f"{key} sekansı boş"


# ─────────────────────────────────────────────────────────────────────────────
# ADIM 3 — Cycle detection
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture(scope='module')
def detected_cycles():
    """Modül genelinde bir kez çalışır; cycle listesi döner."""
    from src.cycle_finderv2 import GraphWork, cycle_info_optimized, cycle_match_based_on_contig_id
    from src.cycle_kmer_hash import filter_cycles_with_kmer

    if not os.path.exists(GFA_FILE):
        pytest.skip("testNitro.gfa bulunamadı")

    gw = GraphWork()
    gw.find_all_path = False
    gw.path_limit = 15
    nodes, edges = gw.parse_gfa(GFA_FILE)
    gw.dfs_iterative(gw.generate_genome_graph(nodes, edges))

    node_lengths = {k: len(v['Sequence']) for k, v in nodes.items()}
    unique_paths = []
    for p in gw.paths:
        if cycle_match_based_on_contig_id(p, node_lengths, unique_paths):
            unique_paths.append(p)

    cycle_list = []
    for p in unique_paths:
        obj = cycle_info_optimized(p, nodes, edges, cycle_list)
        if obj and obj != 'Pass' and 3000 <= obj.length <= 100000:
            obj.name = f'Cycle_{len(cycle_list)+1}'
            cycle_list.append(obj)

    return filter_cycles_with_kmer(cycle_list, 200, 99, 'Cycle')


@pytest.mark.skipif(not os.path.exists(GFA_FILE), reason="testNitro.gfa bulunamadı")
class TestStep3_CycleDetection:
    """Cycle tespiti ve filtreleme çalışıyor mu?"""

    def test_cycles_detected(self, detected_cycles):
        assert len(detected_cycles) > 0, "Hiç cycle bulunamadı"

    def test_no_strand_duplicates(self, detected_cycles):
        """Aynı fiziksel node setini içeren iki cycle olmamalı."""
        seen = []
        for c in detected_cycles:
            stripped = frozenset(n[:-1] for n in c.path)
            assert stripped not in seen, f"Duplicate path: {c.name}"
            seen.append(stripped)

    def test_sizes_within_bounds(self, detected_cycles):
        for c in detected_cycles:
            assert 3000 <= c.length <= 100000, f"{c.name}: {c.length}bp sınır dışı"

    def test_sequences_valid_dna(self, detected_cycles):
        valid = set('ATCGatcg')
        for c in detected_cycles:
            bad = set(c.sequence) - valid
            assert not bad, f"{c.name}: geçersiz nükleotid {bad}"

    def test_sequence_length_consistent(self, detected_cycles):
        for c in detected_cycles:
            assert len(c.sequence) == c.length

    def test_cycle_fasta_written(self, detected_cycles, tmp_path):
        """cycle_analysis() FASTA dosyası yazıyor mu?"""
        from src.cycle_finderv2 import cycle_analysis
        out = str(tmp_path / 'cycles_smoke.fasta')
        cycle_analysis(
            GFA_FILE, out,
            find_all_path=False, path_limit=15,
            min_size_of_cycle=3000, max_size_of_cycle=100000,
            name_prefix_cycle='Cycle',
            min_component_number=1, max_component_number=25,
            k_mer_sim=200, threshold_sim=99
        )
        assert os.path.exists(out), "FASTA dosyası oluşturulmadı"
        assert os.path.getsize(out) > 0, "FASTA dosyası boş"


# ─────────────────────────────────────────────────────────────────────────────
# ADIM 4 — Veritabanı dosyaları mevcut mu?
# ─────────────────────────────────────────────────────────────────────────────

class TestStep4_Databases:
    """Scoring için gerekli DB dosyaları erişilebilir mi?"""

    def test_antibiotics_db_exists(self):
        p = db_path('Antibiotics/protein_fasta_protein_homolog_model.fasta')
        assert os.path.exists(p), f"Antibiotics DB bulunamadı: {p}"

    def test_antibiotics_db_not_empty(self):
        p = db_path('Antibiotics/protein_fasta_protein_homolog_model.fasta')
        if not os.path.exists(p):
            pytest.skip("DB yok")
        assert os.path.getsize(p) > 1000

    def test_is_db_exists(self):
        p = db_path('ISes/_tncentral_nointegrall_isfinder-TNs.fasta')
        assert os.path.exists(p), f"IS DB bulunamadı: {p}"

    def test_xenobiotics_db_exists(self):
        p = db_path('Xenobiotics/Xenobiotics_classified.fasta')
        assert os.path.exists(p), f"Xenobiotics DB bulunamadı: {p}"

    def test_tn_db_exists(self):
        p = db_path('CompTns/Known_Tns.fasta')
        assert os.path.exists(p), f"CompTns DB bulunamadı: {p}"


# ─────────────────────────────────────────────────────────────────────────────
# ADIM 5 — External tool varlığı
# ─────────────────────────────────────────────────────────────────────────────

class TestStep5_ExternalTools:
    """Pipeline'ın gerektirdiği araçlar sistemde var mı?"""

    # Zorunlu araçlar — olmadan scoring çalışmaz
    @pytest.mark.parametrize("tool", ["blastn", "blastp", "blastx", "makeblastdb"])
    def test_blast_tools(self, tool):
        assert tool_available(tool), f"{tool} bulunamadı (BLAST+ kurulu değil?)"

    # Önerilen araçlar — skip ile geç
    @pytest.mark.parametrize("tool", ["prodigal", "minimap2", "samtools"])
    def test_recommended_tools(self, tool):
        if not tool_available(tool):
            pytest.skip(f"{tool} sistemde yok — bu adım çalışmayacak")

    @pytest.mark.parametrize("tool", ["spades.py", "megahit"])
    def test_assembly_tools(self, tool):
        if not tool_available(tool):
            pytest.skip(f"{tool} sistemde yok — assembly çalışmayacak")

    def test_blast_version_runs(self):
        """BLAST kurulu ve çalışır durumda mı?"""
        r = subprocess.run(['blastn', '-version'], capture_output=True, text=True)
        assert r.returncode == 0, "blastn -version başarısız"
        assert 'blastn' in r.stdout.lower() or 'blast' in r.stdout.lower()


# ─────────────────────────────────────────────────────────────────────────────
# ADIM 6 — Scoring (prodigal + blast)
# ─────────────────────────────────────────────────────────────────────────────

class TestStep6_Scoring:
    """Scoring modülü temel fonksiyonlar çalışıyor mu?"""

    def test_calculate_total_score_all_types(self):
        from src.scoringv4ProtBlast import calculate_total_score
        for stype in (0, 1, 2):
            score = calculate_total_score(
                total_score_type=stype, dist_type=1, max_z=20,
                mean_of_CompTns=5850, std_of_CompTns=2586,
                len_of_cycle=8000,
                lst_ant=[80.0], lst_is=[90.0], lst_xe=[],
                comp_number=2
            )
            assert isinstance(score, float), f"type={stype}: float beklendi"
            assert score >= 0, f"type={stype}: negatif skor"

    def test_merge_intervals_basic(self):
        from src.scoringv4ProtBlast import merge_intervals
        result = merge_intervals([(100, 300), (250, 500), (800, 1000)])
        assert result == [(100, 500), (800, 1000)]

    def test_parsing_blast_file_empty(self, tmp_path):
        from src.scoringv4ProtBlast import parsing_blast_file
        empty = str(tmp_path / 'empty.txt')
        open(empty, 'w').close()
        result = parsing_blast_file(empty, 'Antibiotics', 50.0, {})
        assert result == []

    @pytest.mark.skipif(
        not tool_available('prodigal') or not tool_available('makeblastdb') or not tool_available('blastn'),
        reason="prodigal veya BLAST eksik — scoring E2E atlandı"
    )
    def test_scoring_main_e2e(self, tmp_path):
        """scoring_main() gerçek cycle FASTA üzerinde çalışıyor mu?"""
        from src.scoringv4ProtBlast import scoring_main
        from src.cycle_finderv2 import cycle_analysis

        # Önce cycle FASTA oluştur
        cycle_fasta = str(tmp_path / 'cycles.fasta')
        cycle_analysis(
            GFA_FILE, cycle_fasta,
            find_all_path=False, path_limit=15,
            min_size_of_cycle=3000, max_size_of_cycle=100000,
            name_prefix_cycle='Cycle',
            min_component_number=1, max_component_number=25,
            k_mer_sim=200, threshold_sim=99
        )
        if not os.path.exists(cycle_fasta) or os.path.getsize(cycle_fasta) == 0:
            pytest.skip("Cycle FASTA boş — scoring test atlandı")

        scoring_main(
            cycle_fasta, str(tmp_path),
            db_path('Antibiotics/protein_fasta_protein_homolog_model.fasta'),
            db_path('Xenobiotics/Xenobiotics_classified.fasta'),
            db_path('ISes/_tncentral_nointegrall_isfinder-TNs.fasta'),
            db_path('CompTns/Known_Tns.fasta'),
            str(tmp_path),
            mean_of_CompTns=5850, std_of_CompTns=2586,
            total_score_type=0, threshold_final_score=0,   # 0 → tüm cycle'lar çıktsın
            max_z=20, dist_type=1,
            path_of_prodigal='prodigal',
            path_of_blastn='blastn',
            path_of_makeblastdb='makeblastdb',
            path_of_blastx='blastx',
            path_of_blastp='blastp',
            logger_name='picota_analysis'
        )
        tab = str(tmp_path / 'picota_final_tab')
        assert os.path.exists(tab), "picota_final_tab oluşturulmadı"
        assert os.path.getsize(tab) > 0, "picota_final_tab boş"


# ─────────────────────────────────────────────────────────────────────────────
# ADIM 7 — split_cycle_coords_for_is
# ─────────────────────────────────────────────────────────────────────────────

class TestStep7_SplitCycles:
    """split_cycles_from_picota() minimal tabloyla çalışıyor mu?"""

    def _write_minimal_tab(self, path, cycle_id, seq_len):
        import textwrap
        header = "CycleID\tscore2\tIScoords\tAntcoords\tXenocoords\n"
        row = (
            f"{cycle_id}\t95\t100-800\t900-1200\t\n"
        )
        with open(path, 'w') as f:
            f.write(header + row)

    def _write_minimal_fasta(self, path, cycle_id, seq_len):
        with open(path, 'w') as f:
            f.write(f">{cycle_id}\n")
            f.write("A" * seq_len + "\n")

    def test_split_creates_fasta(self, tmp_path):
        from src.split_cycle_coords_for_is import split_cycles_from_picota

        cycle_id = "Cycle_1-len2000-comp2-"
        tab  = str(tmp_path / 'picota_final_tab')
        fasta = str(tmp_path / 'cycles.fasta')
        outdir = str(tmp_path / 'annot')

        self._write_minimal_tab(tab, cycle_id, 2000)
        self._write_minimal_fasta(fasta, cycle_id, 2000)

        results = split_cycles_from_picota(tab, fasta, outdir, split_min_score=50)
        assert len(results) == 1, f"1 dosya beklendi, {len(results)} üretildi"
        assert os.path.exists(results[0])

    def test_split_low_score_skipped(self, tmp_path):
        from src.split_cycle_coords_for_is import split_cycles_from_picota

        cycle_id = "Cycle_1-len2000-comp2-"
        tab  = str(tmp_path / 'picota_final_tab')
        fasta = str(tmp_path / 'cycles.fasta')
        outdir = str(tmp_path / 'annot')

        # score2=30, split_min_score=50 → atlanmalı
        with open(tab, 'w') as f:
            f.write("CycleID\tscore2\tIScoords\tAntcoords\tXenocoords\n")
            f.write(f"{cycle_id}\t30\t100-800\t\t\n")
        self._write_minimal_fasta(fasta, cycle_id, 2000)

        results = split_cycles_from_picota(tab, fasta, outdir, split_min_score=50)
        assert results == [], "Düşük skorlu cycle atlanmalıydı"

    def test_split_missing_is_coords_skipped(self, tmp_path):
        from src.split_cycle_coords_for_is import split_cycles_from_picota

        cycle_id = "Cycle_1-len2000-comp2-"
        tab  = str(tmp_path / 'picota_final_tab')
        fasta = str(tmp_path / 'cycles.fasta')
        outdir = str(tmp_path / 'annot')

        with open(tab, 'w') as f:
            f.write("CycleID\tscore2\tIScoords\tAntcoords\tXenocoords\n")
            f.write(f"{cycle_id}\t95\t\t\t\n")   # IScoords boş
        self._write_minimal_fasta(fasta, cycle_id, 2000)

        results = split_cycles_from_picota(tab, fasta, outdir, split_min_score=50)
        assert results == [], "IS koordinatları olmayan cycle atlanmalıydı"


# ─────────────────────────────────────────────────────────────────────────────
# ADIM 8 — minimap2 + samtools (long-read mapping)
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.skipif(
    not tool_available('minimap2') or not tool_available('samtools'),
    reason="minimap2 veya samtools eksik"
)
class TestStep8_LongReadMapping:
    """minimap2 + samtools pipeline çalışıyor mu?"""

    def _make_tiny_fasta(self, path):
        """Kısa bir referans FASTA yazar."""
        with open(path, 'w') as f:
            f.write(">ref\n" + "ATCGATCGATCG" * 50 + "\n")

    def _make_tiny_fastq(self, path):
        """Tek read'lik bir FASTQ yazar."""
        with open(path, 'w') as f:
            f.write("@read1\n")
            f.write("ATCGATCGATCG" * 10 + "\n")
            f.write("+\n")
            f.write("I" * 120 + "\n")

    def test_minimap2_runs(self, tmp_path):
        ref   = str(tmp_path / 'ref.fasta')
        reads = str(tmp_path / 'reads.fastq')
        sam   = str(tmp_path / 'out.sam')
        self._make_tiny_fasta(ref)
        self._make_tiny_fastq(reads)
        r = subprocess.run(
            ['minimap2', '-ax', 'map-ont', ref, reads],
            capture_output=True, text=True
        )
        assert r.returncode == 0, f"minimap2 hata verdi:\n{r.stderr}"

    def test_samtools_version_runs(self, tmp_path):
        r = subprocess.run(['samtools', 'version'], capture_output=True, text=True)
        if r.returncode != 0:
            pytest.skip(
                f"samtools çalışmıyor (sistem kütüphanesi eksik olabilir — "
                f"'sudo apt install libncurses5' deneyin):\n{r.stderr.strip()}"
            )
        assert 'samtools' in r.stdout.lower()
