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

subject = b"AGAAA" b"TGTTT" b"CACCC" b"GAGGG"
subject_comp = b"TCTTT" b"ACAAA" b"GTGGG" b"CTCCC"
insert = b"CGCA"

contigs = FastaFile(None, mode="")
contigs.sequences = {
    b"identity": subject,
    b"first10bp": subject[0:10],
    # contig matching first 10 bp, but reverse
    b"first10bp-revcomp": subject_comp[9::-1],
    b"second10bp": subject[10:20],
    b"second10bp-revcomp": subject_comp[19:9:-1],
    b"middle10bp": subject[5:15],
    b"middle10bp-revcomp": subject_comp[14:4:-1],
    b"deletion": subject[0:10] + subject[15:20],
    b"insertion": subject[0:10] + insert + subject[10:20],
    b"replacement": subject[0:5] + insert + subject[10:20],
    b"replacement-revcomp": subject_comp[19:9:-1] + insert + subject_comp[4::-1],
    b"replacement-mixed": subject[10:20] + insert + subject_comp[4::-1],
    b"left_overhang": insert + subject,
    b"right_overhang": subject + insert,
}


def test_FastaFile():
    tmp = NamedTemporaryFile(delete=False)
    fa = FastaFile(tmp, "w")
    for acc, seq in contigs.sequences.items():
        fa.put(acc.decode("ASCII"), seq)
    fa.close()
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


def test_scaffold_gap():
    """Two contigs mapped with gap on reference in between should have
    missing piece filled with Ns"""
    rl = RegionList()
    rl.add(1, 10, ("first10bp", 1, 10, 1, 10))
    rl.add(16, 20, ("second10bp", 16, 20, 6, 10))
    sequence = scaffold_contigs(rl, contigs)
    assert len(sequence) == 1
    assert next(iter(sequence.values())) == subject[0:10] + b"n" * 5 + subject[15:20]


def test_scaffold_deletion_in_contig():
    """A contig with a deletion w.r.t. the reference should stay
    intact. No N's reprecenting missing reference bases should be
    inserted.
    """
    rl = RegionList()
    rl.add(1, 10, ("deletion", 1, 10, 1, 10))
    rl.add(16, 20, ("deletion", 16, 20, 11, 15))
    sequence = scaffold_contigs(rl, contigs)
    assert sequence["deletion"] == contigs.get("deletion")


def test_scaffold_split_contig_inserted():
    """A contig is split onto two parts of the reference, and a second
    contig mapped into that split area.
    """
    rl = RegionList()
    # Mapping deletion contig with 5bp gap
    rl.add(1, 10, ("deletion", 1, 10, 1, 10))
    rl.add(16, 20, ("deletion", 16, 20, 11, 15))
    # Mapping something else into the middle
    rl.add(12, 14, ("identity", 12, 14, 12, 14))
    sequence = scaffold_contigs(rl, contigs)
    expected = subject[0:10] + b"n" + subject[11:14] + b"n" + subject[15:20]
    assert len(sequence) == 1
    assert next(iter(sequence.values())) == expected


def test_scaffold_insertion_in_contig():
    """A contig with an insertion w.r.t. the reference should keep that
    insertion.
    """
    rl = RegionList()
    rl.add(1, 10, ("deletion", 1, 10, 1, 10))
    rl.add(16, 20, ("deletion", 16, 20, 11, 15))
    sequence = scaffold_contigs(rl, contigs)
    assert sequence["deletion"] == contigs.get("deletion")


def test_scaffold_replacement_in_contig():
    """A contig in which bases were replaced w.r.t. the reference should
    stay intact.
    """
    rl = RegionList()
    # First match of 5 bp
    rl.add(1, 5, ("replacement", 1, 5, 1, 5))
    # Second match after replacement of 4bp for 5bp in reference
    rl.add(11, 20, ("replacement", 11, 20, 10, 19))
    sequence = scaffold_contigs(rl, contigs)
    assert sequence["replacement"] == contigs.get("replacement")


def test_scaffold_replacement_in_contig_reversed():
    """A contig in which bases were replaced w.r.t. the reference should
    stay intact. Reverse complemented contig
    """
    rl = RegionList()
    # First match of 5 bp
    rl.add(5, 1, ("replacement-revcomp", 5, 1, 15, 20))
    # Second match after replacement of 4bp for 5bp in reference
    rl.add(20, 11, ("replacement-revcomp", 20, 11, 1, 10))
    sequence = scaffold_contigs(rl, contigs)
    assert sequence["replacement-revcomp"] == contigs.get("replacement")


def test_scaffold_replacement_in_contig_mixed():
    """A contig in which bases were replaced w.r.t. the reference should
    stay intact. Second hit is reversed.
    """
    rl = RegionList()
    # First match of 5 bp (reversed)
    rl.add(5, 1, ("replacement-mixed", 5, 1, 15, 20))
    # Second match after replacement of 4bp for 5bp in reference
    rl.add(11, 20, ("replacement-mixed", 11, 20, 1, 10))
    sequence = scaffold_contigs(rl, contigs)
    assert sequence["replacement-mixed"] == contigs.get("replacement")
