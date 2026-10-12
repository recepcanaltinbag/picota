"""
test_assembly_inputs.py – what assembly does with the reads it is given.

SPAdes is replaced by a stub that records its command line and writes a GFA,
so these run without an assembler.
"""

import os
from unittest import mock

import pytest

from src import assembly


def fake_spades(calls):
    def run(cmd, *a, **k):
        calls.append(cmd)
        out = cmd[cmd.index('-o') + 1]
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, 'assembly_graph.gfa'), 'w') as fh:
            fh.write('H\tVN:Z:1.0\n')
    return run


def reads(tmp_path, *names):
    paths = []
    for n in names:
        p = tmp_path / 'reads' / n
        p.parent.mkdir(exist_ok=True)
        p.write_text('@r\nACGT\n+\nIIII\n')
        paths.append(str(p))
    return paths


def assemble(tmp_path, raw, delete_inputs):
    calls = []
    with mock.patch.object(assembly.subprocess, 'run', side_effect=fake_spades(calls)), \
         mock.patch.object(assembly, 'process_gfa_files', side_effect=lambda gfas, _: gfas[0]):
        assembly.assembly_main(
            'S', raw, str(tmp_path / 'out'), 1, '99',
            True, False,                     # quiet, keep_temp_files=False
            'spades.py', 'fastp', True,      # skip_filtering: SPAdes reads raw files
            delete_inputs=delete_inputs)
    return calls


def test_caller_reads_survive_when_not_to_be_deleted(tmp_path):
    raw = reads(tmp_path, 'S_1.fastq', 'S_2.fastq')
    assemble(tmp_path, raw, delete_inputs=False)
    assert all(os.path.exists(f) for f in raw)


def test_downloaded_reads_are_still_cleaned_up_by_default(tmp_path):
    # The SRA pipeline relies on this to keep disk use bounded.
    raw = reads(tmp_path, 'S_1.fastq', 'S_2.fastq')
    assemble(tmp_path, raw, delete_inputs=True)
    assert not any(os.path.exists(f) for f in raw)


@pytest.mark.parametrize('names, flags', [
    (['S_1.fastq', 'S_2.fastq'], ['-1', '-2']),
    (['S.fastq'], ['-s']),            # SPAdes rejects -1 without its -2 mate
])
def test_spades_read_flags(tmp_path, names, flags):
    calls = assemble(tmp_path, reads(tmp_path, *names), delete_inputs=False)
    cmd = calls[0]
    assert [a for a in cmd if a in ('-1', '-2', '-s')] == flags
