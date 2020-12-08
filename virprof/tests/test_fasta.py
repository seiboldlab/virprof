"""
Unit tests for fasta module
"""
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring
# pylint: disable=protected-access
# pylint: disable=relative-beyond-top-level

from tempfile import NamedTemporaryFile

from ..regionlist import RegionList
from ..fasta import FastaFile, merge_contigs

subject = b'AAAAAGGGGGCCCCCTTTTT'
subject_comp = bytes(reversed(subject))
contigs = FastaFile(None, mode='')
contigs.sequences = {
    b'acc1': subject[0:10],
    b'acc1r': subject_comp[9::-1],
    b'acc2': subject[10:20],
    b'acc2r': subject_comp[19:9:-1],
    b'acc3': subject[5:15],
    b'acc3r': subject_comp[14:4:-1],
}


def test_FastaFile():
    tmp = NamedTemporaryFile(delete=False)
    fa = FastaFile(tmp, "w")
    for acc, seq in contigs.sequences.items():
        fa.put(acc.decode('ASCII'), seq)
    fa.close()
    print(fa.name)
    fa = FastaFile(tmp, "r")
    for acc, seq in contigs.sequences.items():
        seq2 = fa.get(acc.decode('ASCII'))
        assert seq == seq2
