"""
test_xenobiotic_curated.py – curated additions to the KEGG xenobiotic set.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), 'picota', 'src')
sys.path.insert(0, SRC)

import build_xenobiotic_db as bx  # noqa: E402
from src.scoringv4ProtBlast import xenobiotic_names  # noqa: E402

CURATED = os.path.join(SRC, 'xenobiotics_curated.fasta')


def test_curated_headers_read_as_their_description():
    for name, seq in bx.read_fasta(CURATED):
        assert ' ' not in name, 'BLAST cuts sseqid at the first space'
        assert seq and seq.isalpha()
        product, gene = xenobiotic_names(name)
        assert product == name.split('|')[4]
        assert gene == name.split('|')[0]


def test_dioxygenase_is_curated():
    names = [n for n, _ in bx.read_fasta(CURATED)]
    assert any(n.startswith('QIH10984.1|') for n in names)


def test_curated_stage_adds_once(tmp_path):
    final = tmp_path / 'Xenobiotics_kegg.fasta'
    final.write_text('>P1|KO:K1|EC:1.1.1.1|PATH:map00984|desc|org\nMKV\n')
    bx.stage_curated(str(final), CURATED)
    once = final.read_text()
    bx.stage_curated(str(final), CURATED)
    assert final.read_text() == once
    assert once.count('>') == 1 + sum(1 for _ in bx.read_fasta(CURATED))
