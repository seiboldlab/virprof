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
            return None
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
    """BLAST Trace-back operations (alignment) string"""

    def __init__(self, btop: str, qstart: int = 1) -> None:
        self._btop = btop
        self._qstart = qstart
        self._length, self._ops = self._parse_btop(btop)

    def __str__(self):
        return self._btop

    def __repr__(self):
        return f"{self.__class__.__name__}({self._btop, self._qstart})"

    @staticmethod
    def _parse_btop(btop: str):
        ops = list()
        length = 0
        for match in re.finditer("([0-9]+)?([A-Z-]{2}|$)", btop):
            digits, letters = match.groups()
            matched = int(digits) if digits else 0
            length += matched
            if letters:
                query, subject = letters.encode("ascii")
                if subject != 45:  # '-':
                    length += 1
            else:
                query, subject = None, None
            ops.append((matched, query, subject))
        return length, ops

    def _get_aligned(
        self,
        sequence: bytes,
        start: int = None,
        end: int = None,
        get_query: bool = True,
    ) -> bytes:
        """Get aligned query/subject in subject coordinates"""
        if start is None:
            start = 0
        elif start < 0:
            raise IndexError(f"{repr(self)}: out of bounds (start={start} < 0")
        else:
            start -= 1

        sequence = sequence[self._qstart - 1 :]

        if end is None:
            end = self._length
        elif end > self._length:
            raise IndexError(
                f"{repr(self)}: out of bounds (end={end}) > len={self._length})"
            )

        aligned = []
        query_ptr = 0
        subject_ptr = 0
        for matched, query, subject in self._ops:
            offset = query_ptr - subject_ptr
            # Process matches
            startpos = max(start + offset, query_ptr)
            endpos = min(end + offset, query_ptr + matched)
            if endpos > startpos:
                aligned += sequence[startpos:endpos]
            query_ptr += matched
            subject_ptr += matched
            if endpos == end + offset:
                break

            if not query:
                LOG.error("ARGH")
                break

            # Process mismatch/indel
            if start <= subject_ptr < end:
                if get_query:
                    aligned.append(query)
                else:
                    aligned.append(subject)
            if subject != 45:  # '-'
                subject_ptr += 1
            if query != 45:  # '-'
                query_ptr += 1
        return bytes(aligned)

    def get_aligned_query(
        self,
        sequence: bytes,
        start: int = None,
        end: int = None,
    ) -> bytes:
        """Return aligned query sequence (with gaps)

        Sequence is identical to passed ``sequence``, with "-"
        characters inserted for gaps.

        Args:
          sequence: query sequence

        """
        return self._get_aligned(sequence, start, end, get_query=True)

    def get_aligned_subject(
        self,
        sequence: bytes,
        start: int = None,
        end: int = None,
    ) -> bytes:
        """Return aligned subject sequence (with gaps)

        Sequence is mutated from passed ``sequence`` and has had gaps
        ("-") inserted.

        Args:
          sequence: query sequence

        """
        return self._get_aligned(sequence, start, end, get_query=False)


def revcomp(sequence: bytes) -> bytes:
    """Reverse complement a bytes sequence"""
    seq = Bio.Seq.Seq(sequence.decode("ASCII"))
    return str(seq.reverse_complement()).encode("ASCII")


def consensus(sequences: List[bytes]) -> bytes:
    """Compute consensus for a set of sequences

    Uncertain bases are filled with 'n'.
    """
    if any(len(sequences[0]) != len(seq) for seq in sequences):
        LOG.error(
            "section sizes differ: %s", " ".join(str(len(seq)) for seq in sequences)
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
        elif best != 45:  # '-'
            bases.append(best)
    return bytes(bases)


def scaffold_contigs(regs: "RegionList", contigs: FastaFile) -> Dict[str, bytes]:
    """Scaffold sequence from BLAST hits

    Args:
      regs: A RegionList containing the blast hits. The data field
        must be qacc, sstart, send, qstart, qend.
      contigs: FastaFile object containing the sequences referenced
        in ``regs``.
    """
    sequence = []  # sequence fragments
    qaccs = set()  # seen contig names
    last_hits = {}
    last_section_from_reference = False
    # Iterate over each disjoined piece
    for section_num, (section_start, section_end, hits) in enumerate(regs):
        section_len = section_end - section_start + 1
        section_seqs = []
        current_hits = {}

        # Iterate over each hit overlapping piece
        for qacc, sstart, send, qstart, _qend, btop in hits:
            qaccs.add(qacc)
            seq = contigs.get(qacc)
            is_forward = send > sstart

            # Compute start and end within contig:
            if is_forward:
                start = qstart + section_start - sstart
            else:
                start = qstart + section_start - send
            end = start + section_len - 1

            # Handle split contig
            current_hits[qacc] = (start, end, is_forward)
            if last_section_from_reference and qacc in last_hits:
                l_start, l_end, l_is_forward = last_hits[qacc]
                if is_forward and l_is_forward:
                    sequence[-1] = seq[l_end : start - 1]
                elif not is_forward and not l_is_forward:
                    sequence[-1] = revcomp(seq[end : l_start - 1])

            if is_forward:
                aligned = btop.get_aligned_query(
                    seq, section_start - sstart + 1, section_end - sstart + 1
                )
            else:
                aligned = btop.get_aligned_query(
                    seq, section_start - send + 1, section_end - send + 1
                )
                aligned = revcomp(aligned)

            # Add outside overhang:
            if len(hits) == 1:

                def get_edge(
                    left, seq=seq, start=start, end=end, is_forward=is_forward
                ):
                    res = seq[: start - 1] if left == is_forward else seq[end:]
                    if is_forward:
                        return res
                    return revcomp(res)

                if section_num == 0:
                    aligned = get_edge(True) + aligned
                if section_num == len(regs) - 1:
                    aligned = aligned + get_edge(False)
            section_seqs.append(aligned)

        # if hits:
        #    section_seqs.append(btop.get_aligned_subject(seq, start + 1, end))

        last_section_from_reference = False
        if len(section_seqs) == 0:
            # No contig covering this piece of reference.
            last_section_from_reference = True
            sequence.append(b"n" * section_len)
        else:
            sequence.append(consensus(section_seqs))
            last_hits = current_hits
    return {"+".join(qaccs): b"".join(sequence)}
