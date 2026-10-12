"""
test_cli.py – the picota.py command-line front end.

Each subcommand must reach the real stage functions with arguments their
signatures accept. The stage functions are patched with autospec=True, so a
call that does not bind to the real signature fails here instead of at run
time. Dispatch tests use a config whose databases are empty placeholder files,
so they run without the real databases. One end-to-end test runs
`picota analysis` on the bundled testNitro.gfa when BLAST and the databases are
available.
"""

import importlib.util
import os
import shutil
import subprocess
import sys
from unittest import mock

import pytest
import yaml

HERE         = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)                        # picota/
INNER_ROOT   = os.path.join(PROJECT_ROOT, 'picota')         # picota/picota/
CLI_PATH     = os.path.join(INNER_ROOT, 'picota.py')
CONFIG       = os.path.join(INNER_ROOT, 'config.yaml')
GFA_FILE     = os.path.join(INNER_ROOT, 'test_data', 'testNitro.gfa')


@pytest.fixture(scope='module')
def cli():
    spec = importlib.util.spec_from_file_location('picota_cli', CLI_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_config(tmp_path, **overrides):
    """The shipped config with placeholder databases and the given overrides.

    overrides: {'options.max_path_expansions': 1234, ...}
    """
    with open(CONFIG) as fh:
        data = yaml.safe_load(fh)
    for attr in ('path_to_antibiotics', 'path_to_xenobiotics', 'path_to_ises', 'path_to_TNs'):
        db = tmp_path / 'dbs' / f'{attr}.fasta'
        db.parent.mkdir(exist_ok=True)
        db.write_text('>x\nACGT\n')
        data['paths'][attr] = str(db)
    for key, value in overrides.items():
        section, name = key.split('.')
        data[section][name] = value
    path = tmp_path / 'config.yaml'
    path.write_text(yaml.safe_dump(data))
    return str(path)


@pytest.fixture
def config(tmp_path):
    return write_config(tmp_path)


def touch_fastqs(tmp_path, *names):
    paths = []
    for n in names:
        p = tmp_path / 'reads' / n
        p.parent.mkdir(exist_ok=True)
        p.write_text('@r\nACGT\n+\nIIII\n')
        paths.append(str(p))
    return paths


def fake_assembly(name, raw, out, *a, **k):
    open(os.path.join(out, '99.gfa'), 'w').close()


def fake_cycles(gfa, out_file, *a, **k):
    with open(out_file, 'w') as fh:
        fh.write('>c\nACGT\n')


def fake_scoring(cycle_file, out_folder, *a, **k):
    with open(os.path.join(out_folder, 'picota_final_tab'), 'w') as fh:
        fh.write(cycle_file + '\n')


class TestParser:

    @pytest.mark.parametrize('argv', [
        ['sra_download', '--sra', 'SRR1', '-o', 'x'],
        ['assembly', '--fastq', 'a_1.fastq', 'a_2.fastq', '-o', 'x', '-t', '8'],
        ['analysis', '--gfa', 'g.gfa', '-o', 'x'],
        ['scoring', '--cycle_folder', 'c/', '-o', 'x'],
        ['all', '--sra', 'SRR1', '--long', 'SRR2', '-o', 'x'],
        ['all', '-o', 'x'],
    ])
    def test_documented_invocations_parse(self, cli, argv):
        args = cli.build_parser().parse_args(argv)
        assert args.command == argv[0]
        assert args.command in cli.COMMANDS

    def test_all_inputs_are_exclusive(self, cli):
        with pytest.raises(SystemExit):
            cli.build_parser().parse_args(['all', '--sra', 'S', '--gfa', 'g', '-o', 'x'])

    def test_no_command_prints_help(self, cli, capsys):
        assert cli.main([]) == 1
        assert 'COMMAND' in capsys.readouterr().out

    @pytest.mark.parametrize('fastq, name', [
        ('SRR123_1.fastq', 'SRR123'),
        ('/data/SRR123_2.fastq.gz', 'SRR123'),
        ('sample_R1_001.fastq.gz', 'sample'),
        ('reads.fq', 'reads'),
    ])
    def test_sample_name_from_fastq(self, cli, fastq, name):
        assert cli.fastq_sample_name(fastq) == name

    def test_sample_name_survives_scoring(self, cli):
        # Scoring cuts SRAID at the first '_' or '.'.
        assert cli.sample_name('my.sample_v2') == 'my-sample-v2'


class TestDispatch:
    """Every subcommand reaches the real stage functions with valid arguments."""

    def test_sra_download(self, cli, tmp_path, config):
        def fake_download(acc, out_path, *a, **k):
            for i in (1, 2):
                open(os.path.join(out_path, f'{acc}_{i}.fastq'), 'w').close()
        with mock.patch.object(cli.pipeline, 'run_sra_down', autospec=True,
                               side_effect=fake_download) as down:
            assert cli.main(['sra_download', '--sra', 'SRR1', '-o', str(tmp_path / 'o'),
                             '-c', config]) == 0
        assert down.call_args.args[0] == 'SRR1'

    def test_assembly(self, cli, tmp_path, config):
        fastqs = touch_fastqs(tmp_path, 'S_1.fastq', 'S_2.fastq')
        with mock.patch.object(cli.pipeline, 'assembly_main', autospec=True,
                               side_effect=fake_assembly) as asm:
            assert cli.main(['assembly', '--fastq', *fastqs, '-o', str(tmp_path / 'o'),
                             '-t', '3', '-c', config]) == 0
        args, kwargs = asm.call_args
        assert args[0] == 'S'
        assert args[1] == fastqs
        assert args[3] == 3                                     # threads override
        assert kwargs['delete_inputs'] is False                 # the user's reads

    def test_assembly_missing_fastq_exits(self, cli, tmp_path, config):
        with pytest.raises(SystemExit):
            cli.main(['assembly', '--fastq', str(tmp_path / 'nope_1.fastq'),
                      '-o', str(tmp_path / 'o'), '-c', config])

    def test_scoring(self, cli, tmp_path, config):
        cycles = tmp_path / 'cycles.fasta'
        cycles.write_text('>c\nACGT\n')
        out = tmp_path / 'out'
        with mock.patch.object(cli.pipeline, 'scoring_main', autospec=True) as score:
            assert cli.main(['scoring', '--cycle_folder', str(cycles), '-o', str(out),
                             '-c', config]) == 0
        assert score.call_args.args[0] == str(cycles)
        assert score.call_args.args[1] == str(out)

    def test_analysis_runs_detection_scoring_and_split(self, cli, tmp_path, config):
        out = tmp_path / 'o'
        with mock.patch.object(cli.pipeline, 'cycle_analysis', autospec=True,
                               side_effect=fake_cycles) as cyc, \
             mock.patch.object(cli.pipeline, 'scoring_main', autospec=True,
                               side_effect=fake_scoring) as score, \
             mock.patch.object(cli.pipeline, 'split_cycles_from_picota', autospec=True,
                               return_value=[]) as split:
            assert cli.main(['analysis', '--gfa', GFA_FILE, '-o', str(out), '-c', config]) == 0
        assert cyc.call_args.args[0] == GFA_FILE
        # <sample>_<graph>.fasta, which scoring reads SRAID and kmer back out of
        assert os.path.basename(cyc.call_args.args[1]) == 'testNitro_testNitro.fasta'
        assert score.call_args.args[1] == str(out / 'scoring' / 'testNitro')
        assert split.called

    def test_analysis_without_cycles_stops_before_scoring(self, cli, tmp_path, config):
        with mock.patch.object(cli.pipeline, 'cycle_analysis', autospec=True), \
             mock.patch.object(cli.pipeline, 'scoring_main', autospec=True) as score:
            assert cli.main(['analysis', '--gfa', GFA_FILE, '-o', str(tmp_path / 'o'),
                             '-c', config]) == 0
        assert not score.called

    def test_second_sample_does_not_reuse_the_first(self, cli, tmp_path, config):
        # Each stage skips work whose output exists, so two samples sharing one
        # output directory must not share their stage directories.
        other = tmp_path / 'other.gfa'
        shutil.copy(GFA_FILE, other)
        out = tmp_path / 'o'
        with mock.patch.object(cli.pipeline, 'cycle_analysis', autospec=True,
                               side_effect=fake_cycles), \
             mock.patch.object(cli.pipeline, 'scoring_main', autospec=True,
                               side_effect=fake_scoring) as score, \
             mock.patch.object(cli.pipeline, 'split_cycles_from_picota', autospec=True,
                               return_value=[]):
            cli.main(['analysis', '--gfa', GFA_FILE, '-o', str(out), '-c', config])
            cli.main(['analysis', '--gfa', str(other), '-o', str(out), '-c', config])
        assert score.call_count == 2
        tab = (out / 'scoring' / 'other' / 'picota_final_tab').read_text()
        assert 'other_other.fasta' in tab

    def test_all_with_fastq_keeps_reads_and_sample_name(self, cli, tmp_path, config):
        fastqs = touch_fastqs(tmp_path, 'SRR9_1.fastq', 'SRR9_2.fastq')
        out = tmp_path / 'o'
        with mock.patch.object(cli.pipeline, 'assembly_main', autospec=True,
                               side_effect=fake_assembly) as asm, \
             mock.patch.object(cli.pipeline, 'cycle_analysis', autospec=True,
                               side_effect=fake_cycles) as cyc, \
             mock.patch.object(cli.pipeline, 'scoring_main', autospec=True,
                               side_effect=fake_scoring), \
             mock.patch.object(cli.pipeline, 'split_cycles_from_picota', autospec=True,
                               return_value=[]):
            assert cli.main(['all', '--fastq', *fastqs, '-o', str(out), '-c', config]) == 0
        assert asm.call_args.kwargs['delete_inputs'] is False
        # The best SPAdes graph is named after its k, not the sample; the cycle
        # file must still carry the sample name.
        assert os.path.basename(cyc.call_args.args[1]) == 'SRR9_99.fasta'

    def test_all_with_sra_runs_the_config_pipeline(self, cli, tmp_path, config):
        with mock.patch.object(cli.pipeline, 'process_accession', autospec=True) as proc:
            assert cli.main(['all', '--sra', 'SRR1', '--long', 'SRR2', '-o', str(tmp_path / 'o'),
                             '-c', config]) == 0
        short, long_, cfg = proc.call_args.args
        assert (short, long_) == ('SRR1', 'SRR2')
        assert cfg.paths.outdir == str(tmp_path / 'o')

    def test_all_without_input_or_sample_sheet_exits(self, cli, tmp_path):
        config = write_config(tmp_path, **{'paths.sra_id_file': str(tmp_path / 'none.csv')})
        with pytest.raises(SystemExit, match='sra_id_file'):
            cli.main(['all', '-o', str(tmp_path / 'o'), '-c', config])

    def test_missing_config_exits(self, cli, tmp_path):
        with pytest.raises(SystemExit):
            cli.main(['analysis', '--gfa', GFA_FILE, '-o', str(tmp_path),
                      '-c', str(tmp_path / 'nope.yaml')])

    def test_missing_database_stops_before_any_work(self, cli, tmp_path):
        # Scoring itself only warns, so the run would "succeed" with an empty table.
        config = write_config(tmp_path, **{'paths.path_to_ises': str(tmp_path / 'none.fasta')})
        with mock.patch.object(cli.pipeline, 'cycle_analysis', autospec=True) as cyc, \
             pytest.raises(SystemExit, match='path_to_ises'):
            cli.main(['analysis', '--gfa', GFA_FILE, '-o', str(tmp_path / 'o'), '-c', config])
        assert not cyc.called

    def test_relative_database_paths_resolve_from_any_directory(self, cli, tmp_path, monkeypatch):
        # The shipped paths are repo-relative ("picota/DBs/..."); run elsewhere,
        # they used to be skipped with a warning.
        rel = os.path.relpath(GFA_FILE, PROJECT_ROOT)   # any repo file stands in for a DB
        config = write_config(tmp_path, **{'paths.path_to_ises': rel})
        monkeypatch.chdir(tmp_path)
        with mock.patch.object(cli.pipeline, 'process_accession', autospec=True) as proc:
            cli.main(['all', '--sra', 'S', '-o', str(tmp_path / 'o'), '-c', config])
        assert proc.call_args.args[2].paths.path_to_ises == os.path.join(PROJECT_ROOT, rel)

    def test_log_goes_to_the_output_directory(self, cli, tmp_path, config):
        # The config's log_file is the config author's; a missing directory there
        # used to crash the run before it started.
        with mock.patch.object(cli.pipeline, 'process_accession', autospec=True) as proc:
            cli.main(['all', '--sra', 'S', '-o', str(tmp_path / 'o'), '-c', config])
        assert proc.call_args.args[2].logging.log_file == str(tmp_path / 'o' / 'picota.log')

    def test_threads_reach_blast(self, cli, tmp_path, config, monkeypatch):
        monkeypatch.setattr(cli.scoring, 'BLAST_THREADS', cli.scoring.BLAST_THREADS)
        with mock.patch.object(cli.pipeline, 'process_accession', autospec=True):
            cli.main(['all', '--sra', 'S', '-o', str(tmp_path / 'o'), '-t', '3', '-c', config])
        assert cli.scoring.BLAST_THREADS == 3

    def test_traversal_and_spades_caps_reach_the_stages(self, cli, tmp_path):
        # These config keys were loaded but never passed on, so setting them
        # had no effect.
        config = write_config(tmp_path, **{'options.max_path_expansions': 1234,
                                           'paths.spades_memory_gb': 7,
                                           'paths.spades_timeout_s': 60})
        out = tmp_path / 'o'
        with mock.patch.object(cli.pipeline, 'cycle_analysis', autospec=True) as cyc:
            cli.main(['analysis', '--gfa', GFA_FILE, '-o', str(out), '-c', config])
        assert cyc.call_args.kwargs['max_expansions'] == 1234

        fastqs = touch_fastqs(tmp_path, 'S_1.fastq')
        with mock.patch.object(cli.pipeline, 'assembly_main', autospec=True,
                               side_effect=fake_assembly) as asm:
            cli.main(['assembly', '--fastq', *fastqs, '-o', str(out / 'asm'), '-c', config])
        assert asm.call_args.kwargs['spades_memory_gb'] == 7
        assert asm.call_args.kwargs['spades_timeout_s'] == 60


class TestEntryPoint:

    @pytest.mark.parametrize('cwd', [PROJECT_ROOT, INNER_ROOT, '/'])
    def test_help_runs_from_any_directory(self, cwd):
        p = subprocess.run([sys.executable, CLI_PATH, '--help'], cwd=cwd,
                           capture_output=True, text=True)
        assert p.returncode == 0, p.stderr
        assert 'analysis' in p.stdout

    def test_version_matches_setup(self, cli):
        with open(os.path.join(PROJECT_ROOT, 'setup.py')) as fh:
            assert f'version="{cli.__version__}"' in fh.read()


_TOOLS = ['prodigal', 'blastn', 'blastp', 'makeblastdb']
_DBS = [os.path.join(INNER_ROOT, 'DBs', *p) for p in (
    ('Antibiotics', 'protein_fasta_protein_homolog_model.fasta'),
    ('Xenobiotics', 'Xenobiotics_kegg.fasta'),
    ('ISes', 'picota_IS_db.fasta'),
    ('CompTns', 'Known_Tns.fasta'))]


@pytest.mark.skipif(any(shutil.which(t) is None for t in _TOOLS)
                    or not all(os.path.exists(d) for d in _DBS),
                    reason='BLAST tools or databases missing')
def test_analysis_end_to_end_from_another_directory(cli, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)        # the README says repo root; it must not matter
    assert cli.main(['analysis', '--gfa', GFA_FILE, '-o', str(tmp_path / 'o'), '-c', CONFIG]) == 0
    out = tmp_path / 'o'
    assert (out / 'cycles' / 'testNitro_testNitro.fasta').exists()
    tab = (out / 'scoring' / 'testNitro' / 'picota_final_tab').read_text().splitlines()
    assert (out / 'scoring' / 'testNitro' / 'picota_enriched.csv').exists()
    assert len(tab) > 1, 'no candidate scored: are the databases being found?'
    assert tab[1].split('\t')[1] == 'testNitro'     # SRAID
