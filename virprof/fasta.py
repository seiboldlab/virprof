"""Helper methods for fasta"""

import logging
import subprocess as sp
import re
from collections import Counter
from typing import Iterator, Sequence, Collection, BinaryIO, Optional, Dict, Tuple, List

import Bio.Seq

LOG = logging.getLogger(__name__)


def read_from_command(args: Sequence[str]) -> Iterator[bytes]:
    """Runs external command yielding output lines"""
    proc = sp.Popen(args, stdout=sp.PIPE)
    assert proc.stdout is not None
    while True:
        line = proc.stdout.readline()
        if not line:
            break
        yield line
    proc.stdout.close()
    proc.wait()


def get_accs_from_fasta(fileobj: BinaryIO) -> Iterator[bytes]:
    """Reads accession numbers from (gzipped) FASTA file"""
    for line in read_from_command(["zgrep", "^>", fileobj.name]):
        if line:
            yield line[1:].split(maxsplit=1)[0]


def filter_fasta(
    filein: BinaryIO, fileout: BinaryIO, accs: Collection[str], remove: bool
) -> Tuple[int, int]:
    """Creates filtered copy of gzipped FASTA file

    Args:
        filein: Object with name property indicating input file
        fileout: Object with name property indicating output file
        accs: Accession numbers
        remove: Whether ``accs`` lists sequences to be removed or sequences
                to be kept.
    """
    unzip = read_from_command(["gzip", "-dc", filein.name])
    outzip = sp.Popen(["gzip", "-c"], stdout=fileout, stdin=sp.PIPE)
    assert outzip.stdin is not None
    skip = True
    fasta_header = b">"[0]
    accs_b = set(acc.encode("ascii") for acc in accs)
    n_seqs_in = n_seqs_out = 0
    for line in unzip:
        if line[0] == fasta_header:
            n_seqs_in += 1
            acc = line[1:].split(maxsplit=1)[0]
            if remove:
                skip = acc in accs_b
            else:
                skip = acc not in accs_b
            if not skip:
                n_seqs_out += 1
        if not skip:
            outzip.stdin.write(line)
    outzip.stdin.close()
    outzip.wait()
    return n_seqs_in, n_seqs_out


class FastaFile:
    """Handles access to GZipp'ed FASTA format file"""

    def __init__(self, iofile: Optional[BinaryIO] = None, mode="r") -> None:
        self.iofile = iofile
        self.mode = mode
        self.sequences = self._try_load_all()
        self._written = set()

        if "w" in mode:
            self.outzip: Optional[sp.Popen[bytes]] = sp.Popen(
                ["gzip", "-c"], stdout=iofile, stdin=sp.PIPE
            )
        else:
            self.outzip = None

    @property
    def name(self) -> str:
        """File name"""
        return self.iofile.name

    def __str__(self) -> str:
        return f"{self.__class__.__name__} {self.name}"

    def close(self) -> None:
        """Close potentially open file handles"""
        if self.outzip is not None and self.outzip.stdin is not None:
            self.outzip.stdin.close()
            self.outzip.wait()

    def __len__(self) -> int:
        return len(self.sequences) if self.sequences else 0

    def _try_load_all(self) -> Optional[Dict[bytes, bytes]]:
        if "r" not in self.mode:
            return {}
        sequences = {}
        fastafile = read_from_command(["gunzip", "-dc", self.iofile.name])
        fasta_header = b">"[0]
        acc = None
        lines = []
        for line in fastafile:
            if line[0] == fasta_header:
                if acc is not None:
                    sequences[acc] = b"".join(lines)
                    lines = []
                acc = line[1:].split(maxsplit=1)[0]
            else:
                lines.append(line.strip())
        if acc is not None:
            sequences[acc] = b"".join(lines)
        return sequences

    def __iter__(self):
        for acc in self.sequences:
            yield acc.decode("utf-8")

    def get(self, acc: str, start: int = 1, stop: int = None) -> bytes:
        """Retrieve a sequence or a part of a sequence

        Params:
          acc: Unique sequence identifier (first word on header line)
          start: Start position of subsequence (1 indexed)
          stop: End position of subsequence (``None`` means until the end)
        """
        if not self.sequences:
            raise IndexError("Empty or write only FASTA")
        seq = self.sequences[acc.encode("utf-8")]
        if stop is None:
            stop = len(seq)
        return seq[start - 1 : stop]

    def put(self, acc: str, sequence: bytes, comment: str = None):
        """Write a sequence

        Params:
          acc: Unique sequence identifier
          sequence: Sequence data
          comment: Additional data to add to seuqence header
        """
        if not self.outzip or not self.outzip.stdin:
            raise IOError("FastaFile not writeable")

        # Make sure we have unique id
        if acc in self._written:
            n = 1
            while f"{acc}_{n}" in self._written:
                n = n + 1
            acc = f"{acc}_{n}"
        self._written.add(acc)

        if comment is not None:
            header = ">{} {}".format(acc, comment).encode("utf-8")
        else:
            header = ">{}".format(acc).encode("utf-8")
        self.outzip.stdin.write(b"\n".join((header, sequence, b"")))


class Btop:
    """BLAST alignment

    All parameters are 1-indexed.

    Args:
       sstart: First base of alignment on subject (reference) sequence
       send: Last base of alignment on subject (reference) sequence
       qstart: First base of alignment on query (contig) sequence
       btop: Blast trace-back operations string (CIGAR like)
    """

    def __init__(self, sstart: int, send: int, qstart: int, btop: str) -> None:
        self._sstart = sstart
        self._send = send
        self._qstart = qstart
        self._btop = btop
        self._subject_length, self._query_length, self._ops = self._parse_btop(btop)

    def __str__(self):
        return self._btop

    def __repr__(self):
        return f"{self.__class__.__name__}{self._sstart, self._send, self._qstart, self._btop}"

    @property
    def subject_length(self):
        return self._subject_length

    @property
    def query_length(self):
        return self._query_length

    def is_forward(self):
        """Check if query and subject have same orientation"""
        return self._send > self._sstart

    @staticmethod
    def _parse_btop(btop: str):
        ops = list()
        subject_length = 0
        query_length = 0
        for match in re.finditer("([0-9]+)?([A-Z-]{2}|$)", btop):
            digits, letters = match.groups()
            matched = int(digits) if digits else 0
            subject_length += matched
            query_length += matched
            if letters:
                query, subject = letters.encode("ascii")
                if subject != 45:  # '-':
                    subject_length += 1
                if query != 45:  # '-':
                    query_length += 1
            else:
                query, subject = None, None
            ops.append((matched, query, subject))
        return subject_length, query_length, ops

    def _get_aligned(
        self,
        sequence: bytes,
        start: int = None,
        end: int = None,
        get_query: bool = True,
    ) -> Tuple[bytes, List[Tuple[int, int]], int, int]:
        """Get aligned query/subject in subject coordinates"""
        if start is None:
            start = 0
        elif start < 0:
            raise IndexError(f"{repr(self)}: out of bounds (start={start} < 0")
        else:
            start -= 1

        sequence = sequence[self._qstart - 1 :]

        if end is None:
            end = self._subject_length
        elif end > self._subject_length:
            raise IndexError(
                f"{repr(self)}: out of bounds (end={end}) > len={self._subject_length})"
            )

        if start >= end:
            raise IndexError(f"{repr(self)}: empty or negative range {start} to {end}")

        aligned = []
        insertions = []
        query_ptr = 0
        subject_ptr = 0
        query_first_base = None
        query_last_base = None
        for matched, query, subject in self._ops:
            offset = query_ptr - subject_ptr
            # Process matches
            startpos = max(start + offset, query_ptr)
            endpos = min(end + offset, query_ptr + matched)
            if endpos > startpos:
                aligned += sequence[startpos:endpos]
                if query_first_base is None:
                    query_first_base = startpos
            query_ptr += matched
            subject_ptr += matched

            if endpos >= end + offset:
                query_last_base = endpos - 1
                break

            if not query:
                raise RuntimeError("Failed to process sequence alignment (b)")

            # Process mismatch/indel
            if start <= subject_ptr < end:
                if query_first_base is None:
                    query_first_base = query_ptr
                if subject == 45:
                    insertions.append(subject_ptr - start)
                if get_query:
                    aligned.append(query)
                else:
                    aligned.append(subject)
            if subject != 45:  # '-'
                subject_ptr += 1
            if query != 45:  # '-'
                assert sequence[query_ptr] == query
                query_ptr += 1
            if subject_ptr == end:
                query_last_base = query_ptr
                break
        else:
            raise RuntimeError("Failed to process sequence alignment")
        return (
            bytes(aligned),
            insertions,
            query_first_base + self._qstart,
            query_last_base + self._qstart,
        )

    def get_aligned_query(
        self,
        sequence: bytes,
        start: int = None,
        end: int = None,
        subject_coordinates: bool = False,
    ) -> bytes:
        """Return aligned query sequence (with gaps)

        All coordinates are 1-indexed closed intervals (matching BLAST).

        Args:
          sequence: Entire query sequence
          start: Position of first base to extract
          end: Position of last base to extract
          subject_coordinates: If true, coordinates are relative to
            the start of the alignment.

        Returns:
          Three-tuple of aligned sequence, and start and end position
            of the extracted sequence œin query coordinates.
        """
        if subject_coordinates:
            if self.is_forward():
                start = start - self._sstart + 1
                end = end - self._sstart + 1
            else:
                start = start - self._send + 1
                end = end - self._send + 1
        aligned, insertions, left, right = self._get_aligned(
            sequence, start, end, get_query=True
        )
        if subject_coordinates and not self.is_forward():
            aligned = revcomp(aligned)
            insertions = [end - start + 1 - pos for pos in reversed(insertions)]
        return aligned, insertions, left, right

    def get_aligned_subject(
        self,
        sequence: bytes,
        start: int = None,
        end: int = None,
        subject_coordinates: bool = False,
    ) -> bytes:
        """Return aligned subject sequence (with gaps)

        See `get_aligned_query`
        """
        if subject_coordinates:
            if self.is_forward():
                start = start - self._sstart + 1
                end = end - self._sstart + 1
            else:
                start = start - self._send + 1
                end = end - self._send + 1
        aligned, insertions, left, right = self._get_aligned(
            sequence, start, end, get_query=False
        )
        if subject_coordinates and not self.is_forward():
            aligned = revcomp(aligned)
            insertions = [end - start + 1 - pos for pos in reversed(insertions)]
        return aligned, insertions, left, right


def revcomp(sequence: bytes) -> bytes:
    """Reverse complement a bytes sequence"""
    seq = Bio.Seq.Seq(sequence.decode("ASCII"))
    return str(seq.reverse_complement()).encode("ASCII")


def consensus(sequences: List[bytes], strip_gaps: bool = True) -> bytes:
    """Compute consensus for a set of sequences

    Uncertain bases are filled with 'n'.

    Raises `ValueError` if the ``sequences`` vary in length.
    """
    if any(len(sequences[0]) != len(seq) for seq in sequences):
        LOG.error("Consensus on sequences with differing lengths:")
        for num, seq in enumerate(sequences):
            LOG.error(f"  {num}: {seq}")
        raise ValueError(
            f"section sizes differ: {', '.join(str(len(seq)) for seq in sequences)}"
        )
    ## Overlapping piece - fill with consensus
    bases = []
    for base_counts in map(Counter, zip(*sequences)):
        if len(base_counts) == 1:
            best, best_count = next(iter(base_counts.items()))
            second_count = 0
        else:
            (best, best_count), (_, second_count) = base_counts.most_common(2)

        if best_count == second_count:
            bases.append(110)  # 110 == 'n'
        elif best != 45 or not strip_gaps:  # '-'
            bases.append(best)
    return bytes(bases)


def combine_inserts(inserts: List[List[int]], sequences: List[bytes]) -> List[bytes]:
    if len(sequences) < 2 or not any(inserts):
        return sequences
    if len(inserts) != len(sequences):
        raise ValueError("Must have the same number of sequences as inserts")
    inserts = [ins.copy() for ins in inserts]
    offsets = [0] * len(sequences)
    result = [b""] * len(sequences)
    last_pos = 0
    while any(inserts):
        pos = min(min(x) for x in inserts if x)
        for i, (ins, seq) in enumerate(zip(inserts, sequences)):
            result[i] += seq[offsets[i] + last_pos : offsets[i] + pos]
            if pos in ins:
                ins.remove(pos)
                result[i] += seq[offsets[i] + pos : offsets[i] + pos + 1]
                offsets[i] += 1
            else:
                result[i] += b"-"
        last_pos = pos
    for i, (ins, seq) in enumerate(zip(inserts, sequences)):
        result[i] += seq[offsets[i] + last_pos :]

    return result


def scaffold_contigs(
    regs: "RegionList", contigs: FastaFile, max_fill_length: int = -1
) -> Dict[str, bytes]:
    """Scaffold sequence from BLAST hits

        TODO: describe exact output decisions

        Args:

          regs: A RegionList containing the blast hits. The data field
            must contain tuples of the identifyer string used in
            ``contigs`` and a `Btop` object
    .
          contigs: FastaFile object containing the sequences
    """
    sequence = []  # sequence fragments
    qaccs = set()  # seen contig names
    result = {}
    # Iterate over each disjoined piece
    mappings = []
    for section_num, (section_start, section_end, hits) in enumerate(regs):
        section_seqs = []  # aligned contig fragments
        section_inserts = []  # matching insert positions
        section_subjs = []  # aligned subject fragments
        section_subjins = []  # matching insert positions
        mappings.append({"sstart": section_start, "send": section_end, "qacc": {}})

        # Iterate over each hit overlapping piece
        for qacc, btop in hits:
            qaccs.add(qacc)
            seq = contigs.get(qacc)
            aligned, inserts, start, end = btop.get_aligned_query(
                seq, section_start, section_end, subject_coordinates=True
            )
            section_seqs.append(aligned)
            section_inserts.append(inserts)
            ref, inserts, _, _ = btop.get_aligned_subject(
                seq, section_start, section_end, subject_coordinates=True
            )
            section_subjs.append(ref)
            section_subjins.append(inserts)
            is_forward = btop.is_forward()
            # Remember piece used
            mappings[-1]["qacc"].setdefault(qacc, []).append((start, end, is_forward))

            # Handle split contig
            if (
                len(mappings) >= 2
                and not mappings[-2]["qacc"]
                and qacc in mappings[-3]["qacc"]
                and len(mappings[-3]["qacc"][qacc]) == 1
            ):
                print("  fixing split")
                l_start, l_end, l_is_forward = mappings[-3]["qacc"][qacc][0]
                if is_forward and l_is_forward:
                    sequence[-1] = seq[l_end : start - 1]
                    if sequence[-1]:
                        mappings[-2]["qacc"][qacc] = [(l_end + 1, start - 1, is_forward)]
                elif not is_forward and not l_is_forward:
                    sequence[-1] = revcomp(seq[end : l_start - 1])
                    if sequence[-1]:
                        mappings[-2]["qacc"][qacc] = [(end + 1, l_start - 1, is_forward)]

        if len(section_seqs) == 0:
            # No contig covering this piece of reference.
            section_len = section_end - section_start + 1
            sequence.append(b"n" * section_len)
            mappings[-1]['scaffold'] = "gap"
        elif len(section_seqs) == 1:
            sequence.append(consensus(section_seqs))
            mappings[-1]['scaffold'] = "unique"
        else:
            section_subjs = combine_inserts(section_subjins, section_subjs)
            subj = consensus(section_subjs, strip_gaps=False)
            section_seqs = combine_inserts(section_inserts, section_seqs)
            sequence.append(consensus(section_seqs + [subj]))
            mappings[-1]['scaffold'] = "consensus"

    # Handle edge overhang
    for left in (True, False):
        # Get the respective outside mapping
        mapping = mappings[0] if left else mappings[-1]
        if len(mapping["qacc"]) != 1:
            continue  # two overlapping contigs on edge, can't handle
        qacc = next(iter(mapping["qacc"]))
        if len(mapping["qacc"][qacc]) != 1:
            continue  # same contig aligned twice on edge, can't handle
        start, end, is_forward = mapping["qacc"][qacc][0]

        # Extract the remaining contig piece
        seq = contigs.get(qacc)
        qstart, qend = (1, start - 1) if left == is_forward else (end + 1, len(seq))
        edge = seq[qstart - 1 : qend]
        if not edge:
            continue
        if not is_forward:
            edge = revcomp(edge)

        # Check if edge overlaps with any existing mapping.
        other = [
            m3
            for m2 in mappings
            if qacc in m2["qacc"]
            for m3 in m2["qacc"][qacc]
        ]
        if any(
            ostart <= qstart <= oend or ostart <= qend <= oend
            for ostart, oend, ois_forward in other
        ):
            continue

        # Add and register overhang
        data = {"qacc": {qacc: [(qstart, qend, is_forward)]}, "scaffold": "edge"}
        if left:
            sequence.insert(0, edge)
            data["sstart"] = mapping["sstart"] - len(edge)
            data["send"] = mapping["sstart"] - 1
            mappings.insert(0, data)
        else:
            sequence.append(edge)
            data["sstart"] = mapping["send"] + 1
            data["send"] = mapping["send"] + len(edge)
            mappings.append(data)


    map_dicts = {}
    parts = []
    base_acc = "+".join(sorted(qaccs))
    acc = f"{base_acc}"
    for seq, mapping in zip(sequence, mappings):
        if 0 < max_fill_length < len(seq) and len(seq) == seq.count(b"n"):
            result[acc] = b"".join(parts)
            acc = f"{base_acc}.{len(result)+1}"
            parts = []
            continue
        parts.append(seq)
        map_dicts.setdefault(acc, []).extend(
            {
                "sstart": mapping["sstart"],
                "send": mapping["send"],
                "qacc": qacc,
                "qstart": qstart,
                "qend": qend,
                "qlen": len(contigs.get(qacc)),
                "reversed": not is_forward,
                "bp": sum(seq.count(base) for base in (b"A", b"G", b"C", b"T")),
                "scaffold": mapping["scaffold"],

            }
            for qacc in mapping["qacc"]
            for qstart, qend, is_forward in mapping["qacc"][qacc]
        )
    result[acc] = b"".join(parts)

    return map_dicts, result
