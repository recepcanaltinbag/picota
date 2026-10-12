"""
test_cli.py – the picota.py command-line front end.

Each subcommand must reach the real stage functions with arguments their
signatures accept. The stage functions are patched with autospec=True, so a
call that does not bind to the real signature fails here instead of at run
time. One end-to-end test runs `picota analysis` on the bundled testNitro.gfa
when BLAST and the databases are available.
"""

import importlib.util
import os
import shutil
import subprocess
import sys
from unittest import mock

import pytest

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


def run(cli, *argv):
    return cli.main(list(argv) + ['-c', CONFIG])


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


class TestDispatch:
    """Every subcommand reaches the real stage functions with valid arguments."""

    def test_sra_download(self, cli, tmp_path):
        def fake_download(acc, out_path, *a, **k):
            for i in (1, 2):
                open(os.path.join(out_path, f'{acc}_{i}.fastq'), 'w').close()
        with mock.patch.object(cli.pipeline, 'run_sra_down', autospec=True,
                               side_effect=fake_download) as down:
            assert run(cli, 'sra_download', '--sra', 'SRR1', '-o', str(tmp_path)) == 0
        assert down.call_args.args[0] == 'SRR1'

    def test_assembly(self, cli, tmp_path):
        def fake_assembly(name, raw, out, *a, **k):
            open(os.path.join(out, 'assembly_graph.gfa'), 'w').close()
        with mock.patch.object(cli.pipeline, 'assembly_main', autospec=True,
                               side_effect=fake_assembly) as asm:
            assert run(cli, 'assembly', '--fastq', 'S_1.fastq', 'S_2.fastq',
                       '-o', str(tmp_path), '-t', '3') == 0
        args = asm.call_args.args
        assert args[0] == 'S'
        assert args[1] == ['S_1.fastq', 'S_2.fastq']
        assert args[3] == 3                                     # threads override

    def test_scoring(self, cli, tmp_path):
        cycles = tmp_path / 'cycles.fasta'
        cycles.write_text('>c\nACGT\n')
        out = tmp_path / 'out'
        with mock.patch.object(cli.pipeline, 'scoring_main', autospec=True) as score:
            assert run(cli, 'scoring', '--cycle_folder', str(cycles), '-o', str(out)) == 0
        assert score.call_args.args[0] == str(cycles)
        assert score.call_args.args[1] == str(out)

    def test_analysis_runs_detection_scoring_and_split(self, cli, tmp_path):
        def fake_cycles(gfa, out_file, *a, **k):
            with open(out_file, 'w') as fh:
                fh.write('>c\nACGT\n')

        def fake_scoring(cycle_file, out_folder, *a, **k):
            open(os.path.join(out_folder, 'picota_final_tab'), 'w').close()

        with mock.patch.object(cli.pipeline, 'cycle_analysis', autospec=True,
                               side_effect=fake_cycles) as cyc, \
             mock.patch.object(cli.pipeline, 'scoring_main', autospec=True,
                               side_effect=fake_scoring) as score, \
             mock.patch.object(cli.pipeline, 'split_cycles_from_picota', autospec=True,
                               return_value=[]) as split:
            assert run(cli, 'analysis', '--gfa', GFA_FILE, '-o', str(tmp_path)) == 0
        assert cyc.call_args.args[0] == GFA_FILE
        assert score.called and split.called

    def test_analysis_without_cycles_stops_before_scoring(self, cli, tmp_path):
        with mock.patch.object(cli.pipeline, 'cycle_analysis', autospec=True), \
             mock.patch.object(cli.pipeline, 'scoring_main', autospec=True) as score:
            assert run(cli, 'analysis', '--gfa', GFA_FILE, '-o', str(tmp_path)) == 0
        assert not score.called

    def test_all_with_sra_runs_the_config_pipeline(self, cli, tmp_path):
        with mock.patch.object(cli.pipeline, 'process_accession', autospec=True) as proc:
            assert run(cli, 'all', '--sra', 'SRR1', '--long', 'SRR2', '-o', str(tmp_path)) == 0
        short, long_, cfg = proc.call_args.args
        assert (short, long_) == ('SRR1', 'SRR2')
        assert cfg.paths.outdir == str(tmp_path)

    def test_traversal_and_spades_caps_reach_the_stages(self, cli, tmp_path):
        # These config keys were loaded but never passed on, so setting them
        # had no effect.
        import yaml
        with open(CONFIG) as fh:
            data = yaml.safe_load(fh)
        data['options']['max_path_expansions'] = 1234
        data['paths']['spades_memory_gb'] = 7
        data['paths']['spades_timeout_s'] = 60
        cfg_path = tmp_path / 'config.yaml'
        cfg_path.write_text(yaml.safe_dump(data))
        out = tmp_path / 'out'

        with mock.patch.object(cli.pipeline, 'cycle_analysis', autospec=True) as cyc:
            cli.main(['analysis', '--gfa', GFA_FILE, '-o', str(out), '-c', str(cfg_path)])
        assert cyc.call_args.kwargs['max_expansions'] == 1234

        with mock.patch.object(cli.pipeline, 'assembly_main', autospec=True) as asm:
            with pytest.raises(SystemExit):        # the mock writes no GFA
                cli.main(['assembly', '--fastq', 'S_1.fastq', '-o', str(out),
                          '-c', str(cfg_path)])
        assert asm.call_args.kwargs['spades_memory_gb'] == 7
        assert asm.call_args.kwargs['spades_timeout_s'] == 60

    def test_missing_config_exits(self, cli, tmp_path):
        with pytest.raises(SystemExit):
            cli.main(['analysis', '--gfa', GFA_FILE, '-o', str(tmp_path),
                      '-c', str(tmp_path / 'nope.yaml')])


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
_DB = os.path.join(INNER_ROOT, 'DBs', 'Antibiotics', 'protein_fasta_protein_homolog_model.fasta')


@pytest.mark.skipif(any(shutil.which(t) is None for t in _TOOLS) or not os.path.exists(_DB),
                    reason='BLAST tools or databases missing')
def test_analysis_end_to_end(cli, tmp_path, monkeypatch):
    monkeypatch.chdir(PROJECT_ROOT)            # config DB paths are repo-relative
    assert run(cli, 'analysis', '--gfa', GFA_FILE, '-o', str(tmp_path)) == 0
    assert (tmp_path / 'cycles' / 'testNitro.fasta').exists()
    assert (tmp_path / 'scoring' / 'picota_final_tab').exists()
    assert (tmp_path / 'scoring' / 'picota_enriched.csv').exists()
