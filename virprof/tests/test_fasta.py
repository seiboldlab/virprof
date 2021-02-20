"""
Unit tests for fasta module
"""
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring
# pylint: disable=protected-access
# pylint: disable=relative-beyond-top-level

from tempfile import NamedTemporaryFile

from ..regionlist import RegionList
from ..fasta import FastaFile, Btop, scaffold_contigs

subject = b"AGAAA" b"TGTTT" b"CACCC" b"GAGGG"
subject_comp = b"TCTTT" b"ACAAA" b"GTGGG" b"CTCCC"

mutated = b"AGACC" b"TGTTT" b"CACCC" b"GAGGG"
mutated_btop = "3CACA15"
insertion = b"AGAAAG" b"TGTTT" b"CACCC" b"GAGGG"
insertion_subj = b"AGAAA-" b"TGTTT" b"CACCC" b"GAGGG"
insertion_btop = "5G-15"
deletion = b"AGAA" b"TGTTT" b"CACCC" b"GAGGG"
deletion_alig = b"AGAA-" b"TGTTT" b"CACCC" b"GAGGG"
deletion_btop = "4-A15"

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


def test_Btop_identity():
    """Simple BTOP of exact match"""
    btop = Btop(str(len(subject)))
    assert btop.get_aligned_query(subject) == subject
    assert btop.get_aligned_subject(subject) == subject


def test_Btop_mutation():
    """Simple BTOP with 2bp mutated"""
    btop = Btop(mutated_btop)
    assert btop.get_aligned_query(mutated) == mutated
    assert btop.get_aligned_subject(mutated) == subject


def test_Btop_insertion():
    """Simple BTOP with one insertion in query"""
    btop = Btop(insertion_btop)
    assert btop.get_aligned_query(insertion) == insertion
    assert btop.get_aligned_subject(insertion) == insertion_subj


def test_Btop_deletion():
    """Simple BTOP with one deletion in query"""
    btop = Btop(deletion_btop)
    assert btop.get_aligned_query(deletion) == deletion_alig
    assert btop.get_aligned_subject(deletion) == subject


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
    """A contig with two hits to the reference where one is reversed.

    There is no way to know which orientation the middle bases should
    have. This situation really should be avoided when choosing the
    reference.

    Expecting N's in the middle.
    """
    rl = RegionList()
    # First match of 5 bp - REVERSED
    rl.add(5, 1, ("replacement-mixed", 5, 1, 15, 20))
    # Second match after replacement of 4bp for 5bp in reference - FORWARD
    rl.add(11, 20, ("replacement-mixed", 11, 20, 1, 10))
    sequence = scaffold_contigs(rl, contigs)
    expected = (
        subject[0:5]  # first part
        + b"n" * 5  # gap not replaced
        + subject[10:]  # second part
        + contigs.get("replacement-mixed")[10:]  # right overhang
    )
    assert sequence["replacement-mixed"] == expected


def test_scaffold_left_overhanging_contig():
    """A contig overlapping the reference on its left end should stay intact"""
    rl = RegionList()
    rl.add(1, 20, ("left_overhang", 1, 20, 5, 24))
    sequence = scaffold_contigs(rl, contigs)
    assert sequence["left_overhang"] == contigs.get("left_overhang")


def test_scaffold_right_overhanging_contig():
    """A contig overlapping the reference on its left end should stay intact"""
    rl = RegionList()
    rl.add(1, 20, ("right_overhang", 1, 20, 1, 20))
    sequence = scaffold_contigs(rl, contigs)
    assert sequence["right_overhang"] == contigs.get("right_overhang")
