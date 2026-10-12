"""
test_build_is_db.py – repair and completion of the IS reference set.

The glued-record shapes below are the ones found in the distributed
_tncentral_nointegrall_isfinder-TNs.fasta, shortened.
"""

import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'scripts'))

import build_is_db as b  # noqa: E402


def lines(text):
    return text.strip('\n').split('\n')


def test_clean_file_is_unchanged():
    recs, repaired = b.parse_merged(lines(">A-X1\nacgt\nacgt\n>ISFinder_IS1\nGGCC\n"))
    assert recs == [('A-X1', 'acgtacgt'), ('ISFinder_IS1', 'GGCC')]
    assert repaired == []


def test_header_glued_mid_line_with_its_sequence():
    recs, repaired = b.parse_merged(lines(
        ">In1-X1\naaaacgaca>In1223-KX784502tgtcattt\nggcc\n>Tn1-Y2\nacgt\n"))
    assert recs == [('In1-X1', 'aaaacgaca'),
                    ('In1223-KX784502', 'tgtcatttggcc'),
                    ('Tn1-Y2', 'acgt')]
    assert repaired == ['In1223-KX784502']


def test_header_split_across_the_line_wrap():
    recs, _ = b.parse_merged(lines(">A-X1\nggaactgaca>T\nn6023-GU562437.2ggcactg\nttgc\n"))
    assert recs[1] == ('Tn6023-GU562437.2', 'ggcactgttgc')


def test_accession_split_across_the_line_wrap():
    recs, _ = b.parse_merged(lines(">A-X1\nacga>In781_p-CP0\n03505acgt\n"))
    assert recs[1] == ('In781_p-CP003505', 'acgt')


def test_lowercase_letters_in_the_name_are_kept():
    recs, _ = b.parse_merged(lines(">A-X1\nacga>Tn7CVna1-CP009977.1tacg\n"))
    assert recs[1] == ('Tn7CVna1-CP009977.1', 'tacg')


def test_consecutive_glued_records():
    recs, repaired = b.parse_merged(lines(">A-X1\nac>B-Y1gg\ntt>C-Z1aa\n"))
    assert recs == [('A-X1', 'ac'), ('B-Y1', 'ggtt'), ('C-Z1', 'aa')]
    assert repaired == ['B-Y1', 'C-Z1']


def test_unrecognisable_glued_record_stops_the_build():
    with pytest.raises(b.BuildError):
        b.parse_merged(lines(">A-X1\nac>nothinglikeaheaderacgt\n"))


def test_is_fna_names_drop_family_and_group():
    recs = b.parse_is_fna(lines(">IS26_IS6_unknown\nACGT\n>ISBj2_B_IS5_IS5\nGG\n"))
    assert recs == [('IS26', 'ACGT'), ('ISBj2_B', 'GG')]


def test_missing_elements_are_added_and_present_ones_are_not():
    merged = lines(">ISFinder_IS26\nacgtacgtaa\n>Tn6023-GU562437.2\nggggcccc\n")
    is_fna = lines(
        ">IS26_IS6_unknown\nTTTTTTTTTT\n"     # present by name
        ">IS15DII_IS6_unknown\nTTACGTACGT\n"  # present by sequence (reverse complement)
        ">IS100_IS21_unknown\nAAACCCGGG\n")   # missing
    records, report = b.build(merged, is_fna)
    assert report['added'] == ['ISFinder_IS100']
    assert records[-1] == ('ISFinder_IS100', 'AAACCCGGG')
    assert report['total'] == 3


def test_written_file_round_trips(tmp_path):
    out = tmp_path / 'db.fasta'
    b.write_fasta([('A-X1', 'a' * 130)], str(out))
    text = out.read_text().split('\n')
    assert text[0] == '>A-X1' and len(text[1]) == b.WRAP
    recs, _ = b.parse_merged(b.read_fasta_lines(str(out)))
    assert recs == [('A-X1', 'a' * 130)]
