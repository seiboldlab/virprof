"""
Unit tests for fasta module
"""
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring
# pylint: disable=protected-access
# pylint: disable=relative-beyond-top-level

from tempfile import NamedTemporaryFile

from ..regionlist import RegionList
from ..fasta import FastaFile, scaffold_contigs

subject =  b"AAAAATTTTTCCCCCGGGGG"
subject_comp = b"TTTTTAAAAAGGGGGCCCCC"

contigs = FastaFile(None, mode="")
contigs.sequences = {
    b"first10bp": subject[0:10],
    # contig matching first 10 bp, but reverse
    b"first10bp-revcomp": subject_comp[9::-1],
    b"second10bp": subject[10:20],
    b"second10bp-revcomp": subject_comp[19:9:-1],
    b"middle10bp": subject[5:15],
    b"middle10bp-revcomp": subject_comp[14:4:-1],
}


def test_FastaFile():
    tmp = NamedTemporaryFile(delete=False)
    fa = FastaFile(tmp, "w")
    for acc, seq in contigs.sequences.items():
        fa.put(acc.decode("ASCII"), seq)
    fa.close()
    print(fa.name)
    fa = FastaFile(tmp, "r")
    for acc, seq in contigs.sequences.items():
        seq2 = fa.get(acc.decode("ASCII"))
        assert seq == seq2


def test_scaffold_contigs_simple():
    rl = RegionList()
    rl.add(1, 10, ("first10bp", 1, 10, 1, 10))
    sequence = scaffold_contigs(rl, contigs)
    assert sequence["first10bp"] == subject[0:10]


def test_scaffold_contigs_simple_reversed():
    rl = RegionList()
    rl.add(10, 1, ("first10bp-revcomp", 10, 1, 1, 10))
    sequence = scaffold_contigs(rl, contigs)
    assert sequence["first10bp-revcomp"] == subject[0:10]


def test_scaffold_contigs_split():
    rl = RegionList()
    rl.add(1, 10, ("first10bp", 1, 10, 1, 10))
    rl.add(11, 20, ("second10bp", 11, 20, 1, 10))
    sequence = scaffold_contigs(rl, contigs)
    assert len(sequence) == 1
    assert next(iter(sequence.values())) == subject


def test_scaffold_contigs_overlap():
    rl = RegionList()
    rl.add(1, 10, ("first10bp", 1, 10, 1, 10))
    rl.add(11, 20, ("second10bp", 11, 20, 1, 10))
    rl.add(6, 15, ("middle10bp", 6, 15, 1, 10))
    sequence = scaffold_contigs(rl, contigs)
    assert len(sequence) == 1
    assert next(iter(sequence.values())) == subject
