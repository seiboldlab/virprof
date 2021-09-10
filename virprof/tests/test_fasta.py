"""
Unit tests for fasta module
"""
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring
# pylint: disable=protected-access
# pylint: disable=relative-beyond-top-level

from tempfile import NamedTemporaryFile
from itertools import product

from ..regionlist import RegionList
from ..fasta import FastaFile, Btop, scaffold_contigs, revcomp, combine_inserts

subject = b"AGAAA" b"TGTTT" b"CACCC" b"GAGGG"
subject_comp = b"TCTTT" b"ACAAA" b"GTGGG" b"CTCCC"

mutated = b"AGACC" b"TGTTT" b"CACCC" b"GAGGG"
mutated_btop = "3CACA15"
mutated1 = b"AGAAC" b"TGTTT" b"CACCC" b"GAGGG"
mutated1_comp = b"TCTTG" b"ACAAA" b"GTGGG" b"CTCCC"
mutated1_comp_btop = "15GT4"
deletion = b"AGAA" b"TGTTT" b"CACCC" b"GAGGG"
deletion_alig = b"AGAA-" b"TGTTT" b"CACCC" b"GAGGG"
deletion_btop = "4-A15"

insert = b"CGCA"
insert_comp = b"GCGT"

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
    b"replacement-revcomp": subject_comp[19:9:-1]
    + insert_comp[::-1]
    + subject_comp[4::-1],
    b"replacement-mixed": subject[10:20] + insert + subject_comp[4::-1],
    b"mutated1-revcomp": mutated1_comp[::-1],
    b"left_overhang": insert + subject,
    b"right_overhang": subject + insert,
    b"left_revcomp": subject_comp[9::-1] + subject[10:20],
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
    btop = Btop(1, 20, 1, str(len(subject)))
    assert btop.get_aligned_query(subject) == (subject, [], 1, len(subject))
    assert btop.get_aligned_subject(subject) == (subject, [], 1, len(subject))
    assert btop.get_aligned_query(subject, 1, len(subject)) == (
        subject,
        [],
        1,
        len(subject),
    )
    assert btop.get_aligned_subject(subject, 1, len(subject)) == (
        subject,
        [],
        1,
        len(subject),
    )
    assert btop.get_aligned_query(subject, 10, 15) == (subject[9:15], [], 10, 15)
    assert btop.get_aligned_subject(subject, 10, 15) == (subject[9:15], [], 10, 15)


def test_Btop_mutation():
    """Simple BTOP with 2bp mutated (4 and 5)"""
    btop = Btop(1, 20, 1, mutated_btop)
    assert btop.get_aligned_query(mutated) == (mutated, [], 1, len(mutated))
    assert btop.get_aligned_subject(mutated) == (subject, [], 1, len(mutated))


def test_Btop_mutation1_comp():
    """Simple BTOP with 2bp mutated (4 and 5)"""
    btop = Btop(1, 20, 1, mutated1_comp_btop)
    query = mutated1_comp[::-1]
    subj = subject_comp[::-1]
    assert btop.get_aligned_query(query) == (query, [], 1, len(query))
    assert btop.get_aligned_subject(query) == (subj, [], 1, len(query))
    for lf, ri in product(*[range(1, 6)] * 2):
        if ri < lf:
            continue
        assert btop.get_aligned_query(query, lf, ri) == (query[lf - 1 : ri], [], lf, ri), (
            lf,
            ri,
        )
        assert btop.get_aligned_subject(query, lf, ri) == (subj[lf - 1 : ri], [], lf, ri), (
            lf,
            ri,
        )


def test_Btop_insertion():
    """Simple BTOP with one insertion in query"""
    insertion = b"AGAAAG" b"TGTTT" b"CACCC" b"GAGGG"
    insertion_subject = b"AGAAA-" b"TGTTT" b"CACCC" b"GAGGG"
    insertion_btop = "5G-15"

    btop = Btop(1, 20, 1, insertion_btop)
    assert btop.get_aligned_query(insertion) == (insertion, [5], 1, len(insertion))
    assert btop.get_aligned_subject(insertion) == (
        insertion_subject,
        [5],
        1,
        len(insertion),
    )
    # Insertion between subject 5 and 6 is prepended on position 6
    assert btop.get_aligned_query(insertion, 6, 6) == (insertion[5:7], [0], 6, 7)
    # Expecting gap on aligned subject
    assert btop.get_aligned_subject(insertion, 6, 6) == (
        insertion_subject[5:7],
        [0],
        6,
        7,
    )
    # Full match is 11 bp long
    assert btop.get_aligned_query(insertion, 1, 20) == (
        insertion,
        [5],
        1,
        len(insertion),
    )
    assert btop.get_aligned_subject(insertion, 1, 20) == (
        insertion_subject,
        [5],
        1,
        len(insertion),
    )
    # Part before is normal
    assert btop.get_aligned_query(insertion, 1, 5) == (insertion[0:5], [], 1, 5)
    assert btop.get_aligned_subject(insertion, 1, 5) == (subject[0:5], [], 1, 5)
    # Part after has offset in aligned query template
    assert btop.get_aligned_query(insertion, 11, 15) == (insertion[11:16], [], 12, 16)
    # But not for subject
    assert btop.get_aligned_subject(insertion, 11, 15) == (subject[10:15], [], 12, 16)


def test_Btop_deletion():
    """Simple BTOP with one deletion in query

    Base 5 is deleted in query
    """
    btop = Btop(1, 20, 1, deletion_btop)
    assert btop.get_aligned_query(deletion) == (deletion_alig, [], 1, len(deletion))
    assert btop.get_aligned_subject(deletion) == (subject, [], 1, len(deletion))
    # At subject position 5 we get a gap character
    assert btop.get_aligned_query(deletion, 5, 5) == (b"-", [], 5, 5)
    assert btop.get_aligned_subject(deletion, 5, 5) == (subject[4:5], [], 5, 5)
    # Part before
    assert btop.get_aligned_query(deletion, 1, 3) == (subject[0:3], [], 1, 3)
    assert btop.get_aligned_subject(deletion, 1, 3) == (subject[0:3], [], 1, 3)
    # Part spanning
    assert btop.get_aligned_query(deletion, 5, 6) == (deletion_alig[4:6], [], 5, 5)
    assert btop.get_aligned_subject(deletion, 5, 6) == (subject[4:6], [], 5, 5)
    assert btop.get_aligned_query(deletion, 6, 7) == (deletion_alig[5:7], [], 5, 6)
    assert btop.get_aligned_subject(deletion, 6, 7) == (subject[5:7], [], 5, 6)
    assert btop.get_aligned_query(deletion, 5, 7) == (deletion_alig[4:7], [], 5, 6)
    assert btop.get_aligned_subject(deletion, 5, 7) == (subject[4:7], [], 5, 6)
    # Part after
    assert btop.get_aligned_query(deletion, 11, 15) == (subject[10:15], [], 10, 14)
    assert btop.get_aligned_subject(deletion, 11, 15) == (subject[10:15], [], 10, 14)


def test_scaffold_contigs_simple():
    rl = RegionList()
    rl.add(1, 10, ("first10bp", Btop(1, 10, 1, "10")))
    _, sequence = scaffold_contigs(rl, contigs)
    assert sequence["first10bp"] == subject[0:10]


def test_scaffold_contigs_simple_reversed():
    rl = RegionList()
    rl.add(10, 1, ("first10bp-revcomp", Btop(10, 1, 1, "10")))
    _, sequence = scaffold_contigs(rl, contigs)
    assert sequence["first10bp-revcomp"] == subject[0:10]


def test_scaffold_contigs_split():
    rl = RegionList()
    rl.add(1, 10, ("first10bp", Btop(1, 10, 1, "10")))
    rl.add(11, 20, ("second10bp", Btop(11, 20, 1, "10")))
    _, sequence = scaffold_contigs(rl, contigs)
    assert len(sequence) == 1
    assert next(iter(sequence.values())) == subject


def test_scaffold_contigs_overlap():
    rl = RegionList()
    rl.add(1, 10, ("first10bp", Btop(1, 10, 1, "10")))
    rl.add(11, 20, ("second10bp", Btop(11, 20, 1, "10")))
    rl.add(6, 15, ("middle10bp", Btop(6, 15, 1, "10")))
    _, sequence = scaffold_contigs(rl, contigs)
    assert len(sequence) == 1
    assert next(iter(sequence.values())) == subject


def test_scaffold_gap():
    """Two contigs mapped with gap on reference in between should have
    missing piece filled with Ns or be split if too long"""
    rl = RegionList()
    rl.add(1, 10, ("first10bp", Btop(1, 10, 1, "10")))
    rl.add(16, 20, ("second10bp", Btop(16, 20, 6, "10")))
    _, sequence = scaffold_contigs(rl, contigs)
    assert len(sequence) == 1
    assert next(iter(sequence.values())) == subject[0:10] + b"n" * 5 + subject[15:20]
    _, sequence = scaffold_contigs(rl, contigs, 5)
    assert len(sequence) == 1
    _, sequence = scaffold_contigs(rl, contigs, 4)
    assert len(sequence) == 2
    assert sequence["first10bp+second10bp"] == subject[0:10]
    assert sequence["first10bp+second10bp.2"] == subject[15:20]


def test_scaffold_deletion_in_contig():
    """A contig with a deletion w.r.t. the reference should stay
    intact. No N's reprecenting missing reference bases should be
    inserted.
    """
    rl = RegionList()
    rl.add(1, 10, ("deletion", Btop(1, 10, 1, "20")))
    rl.add(16, 20, ("deletion", Btop(16, 20, 11, "20")))
    _, sequence = scaffold_contigs(rl, contigs)
    assert sequence["deletion"] == contigs.get("deletion")


def test_scaffold_split_contig_inserted():
    """A contig is split onto two parts of the reference, and a second
    contig mapped into that split area.
    """
    rl = RegionList()
    # Mapping deletion contig with 5bp gap
    rl.add(1, 10, ("deletion", Btop(1, 10, 1, "10")))
    rl.add(16, 20, ("deletion", Btop(16, 20, 11, "20")))
    # Mapping something else into the middle
    rl.add(12, 14, ("identity", Btop(12, 14, 12, "20")))
    _, sequence = scaffold_contigs(rl, contigs)
    expected = subject[0:10] + b"n" + subject[11:14] + b"n" + subject[15:20]
    assert len(sequence) == 1
    assert next(iter(sequence.values())) == expected


def test_scaffold_insertion_in_contig():
    """A contig with an insertion w.r.t. the reference should keep that
    insertion.
    """
    rl = RegionList()
    rl.add(1, 10, ("deletion", Btop(1, 10, 1, "10")))
    rl.add(16, 20, ("deletion", Btop(16, 20, 11, "20")))
    _, sequence = scaffold_contigs(rl, contigs)
    assert sequence["deletion"] == contigs.get("deletion")


def test_scaffold_replacement_in_contig():
    """A contig in which bases were replaced w.r.t. the reference should
    stay intact.
    """
    rl = RegionList()
    # First match of 5 bp
    rl.add(1, 5, ("replacement", Btop(1, 5, 1, "10")))
    # Second match after replacement of 4bp for 5bp in reference
    rl.add(11, 20, ("replacement", Btop(11, 20, 11, "20")))
    _, sequence = scaffold_contigs(rl, contigs)
    assert sequence["replacement"] == contigs.get("replacement")


def test_scaffold_replacement_in_contig_reversed():
    """A contig in which bases were replaced w.r.t. the reference should
    stay intact. Reverse complemented contig
    """
    rl = RegionList()
    # First match of 5 bp
    rl.add(5, 1, ("replacement-revcomp", Btop(5, 1, 15, "20")))
    # Second match after replacement of 4bp for 5bp in reference
    rl.add(20, 11, ("replacement-revcomp", Btop(20, 11, 1, "20")))
    _, sequence = scaffold_contigs(rl, contigs)
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
    rl.add(5, 1, ("replacement-mixed", Btop(5, 1, 15, "5")))
    # Second match after replacement of 4bp for 5bp in reference - FORWARD
    rl.add(11, 20, ("replacement-mixed", Btop(11, 20, 1, "10")))
    _, sequence = scaffold_contigs(rl, contigs)
    expected = (
        subject[0:5]  # first part
        + b"n" * 5  # gap not replaced
        + subject[10:]  # second part
    )
    assert sequence["replacement-mixed"] == expected


def test_scaffold_mutation_revcomp():
    """Reverse complemented contig with 1 bp changed"""
    rl = RegionList()
    # Whole sequence
    rl.add(20, 1, ("mutated1-revcomp", Btop(20, 1, 1, mutated1_comp_btop)))
    _, sequence = scaffold_contigs(rl, contigs)
    assert sequence["mutated1-revcomp"] == mutated1
    rl = RegionList()
    # Disjoint alignment with 2 hits
    # First match of 5 bp
    rl.add(5, 1, ("mutated1-revcomp", Btop(5, 1, 15, mutated1_comp_btop)))
    # Second match after replacement of 4bp for 5bp in reference
    rl.add(20, 11, ("mutated1-revcomp", Btop(20, 11, 1, mutated1_comp_btop)))
    _, sequence = scaffold_contigs(rl, contigs)
    assert sequence["mutated1-revcomp"] == mutated1


def test_scaffold_left_overhanging_contig():
    """A contig overlapping the reference on its left end should stay intact"""
    rl = RegionList()
    rl.add(1, 20, ("left_overhang", Btop(1, 20, 5, "24")))
    _, sequence = scaffold_contigs(rl, contigs)
    assert sequence["left_overhang"] == contigs.get("left_overhang")


def test_scaffold_right_overhanging_contig():
    """A contig overlapping the reference on its left end should stay intact"""
    rl = RegionList()
    rl.add(1, 20, ("right_overhang", Btop(1, 20, 1, "20")))
    _, sequence = scaffold_contigs(rl, contigs)
    assert sequence["right_overhang"] == contigs.get("right_overhang")


def test_scaffold_mixed():
    rl = RegionList()
    rl.add(1, 10, ("left_revcomp", Btop(10, 1, 1, "10")))
    rl.add(11, 20, ("left_revcomp", Btop(11, 20, 11, "10")))
    _, sequence = scaffold_contigs(rl, contigs)
    assert sequence["left_revcomp"] == subject


def test_scaffold_bug():
    """Two contigs eacho aligned to reference need to match up right"""

    # The subject region is 14 bp
    slen = 14
    sstart = 5
    send = sstart + slen - 1
    subject = b"x" * (sstart - 1) + b"AATAATAATAATAA"
    assert len(subject) == send

    # The first contig has 1 'G' insertion between subject 1 and 2
    btop_1 = Btop(sstart, send, 1, "1G-13")
    assert btop_1.subject_length == slen
    assert btop_1.query_length == slen + 1
    contig_1 = subject[sstart - 1 : sstart] + b"G" + subject[sstart:]
    assert len(contig_1) == slen + 1

    alig_1, ins_1, start_1, end_1 = btop_1.get_aligned_query(
        contig_1, sstart, send, subject_coordinates=True
    )
    assert (start_1, end_1) == (1, slen + 1)
    assert len(alig_1) == slen + 1
    subj_1, sins_1, _, _ = btop_1.get_aligned_subject(
        contig_1, sstart, send, subject_coordinates=True
    )

    # The second contig has  1 'GG' insertion between subject 1 and 2,
    # an GG insertion between subject 7 and 8,
    # and a deleted A at subject 13
    btop_2 = Btop(sstart, send, 1, "1G-G-6G-G-5-A1")
    assert btop_2.subject_length == slen
    assert btop_2.query_length == slen + 3
    contig_2 = (
        subject[sstart - 1 : sstart]
        + b"GG"
        + subject[sstart : sstart + 6]
        + b"GG"
        + subject[sstart + 6 : sstart + 12]
        + subject[sstart + 14 :]
    )
    assert len(contig_2) == slen + 3

    alig_2, ins_2, start_2, end_2 = btop_2.get_aligned_query(
        contig_2, sstart, send, subject_coordinates=True
    )
    assert (start_2, end_2) == (1, slen + 3)
    assert len(alig_2) == slen + 4
    subj_2, sins_2, _, _ = btop_2.get_aligned_subject(
        contig_2, sstart, send, subject_coordinates=True
    )

    btop_2r = Btop(send, sstart, 1, "1-T5C-C-6C-C-1")
    contig_2r = revcomp(contig_2)
    alig_2r, ins_2r, start_2r, end_2r = btop_2r.get_aligned_query(
        contig_2r, sstart, send, subject_coordinates=True
    )
    assert alig_2r == alig_2
    assert ins_2r == ins_2
    subj_2r, sins_2r, _, _ = btop_2r.get_aligned_subject(
        contig_2r, sstart, send, subject_coordinates=True
    )
    assert subj_2r == subj_2
    assert sins_2r == sins_2

    merged_1, merged_2 = combine_inserts([ins_1, ins_2r], [alig_1, alig_2r])
    assert len(merged_1) == len(merged_2) == slen + 4

    msubj_1, msubj_2 = combine_inserts([sins_1, sins_2], [subj_1, subj_2r])
    assert (
        bytes([x for x in msubj_1 if x != 45])
        == bytes([x for x in msubj_2 if x != 45])
        == subject[sstart - 1 : send]
    )

    fasta = FastaFile(None, mode="")
    fasta.sequences = {b"c1": contig_1, b"c2": contig_2r}
    rl = RegionList()
    rl.add(sstart, send, ("c1", btop_1))
    rl.add(sstart, send, ("c2", btop_2r))
    _, sequences = scaffold_contigs(rl, fasta)
    sequence = next(iter(sequences.values()))
    assert sequence == subject[sstart - 1 : sstart] + b"G" + subject[sstart:send]
