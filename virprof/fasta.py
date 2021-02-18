"""Helper methods for fasta"""

import logging
import subprocess as sp
from collections import Counter
from typing import Iterator, Sequence, Collection, BinaryIO, Optional, Dict, Tuple

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


def revcomp(sequence: bytes) -> bytes:
    seq = Bio.Seq.Seq(sequence.decode("ASCII"))
    return str(seq.reverse_complement()).encode("ASCII")


def scaffold_contigs(regs: "RegionList", contigs: FastaFile) -> Dict[str, bytes]:
    sequence = []  # sequence fragments
    qaccs = set()  # seen contig names
    last_hits = {}
    last_section_from_reference = False
    # Iterate over each disjoined piece
    for section_start, section_end, hits in regs:
        section_len = section_end - section_start + 1
        section_seqs = []

        # Iterate over each hit overlapping piece
        for qacc, sstart, send, qstart, qend in hits:
            seq = contigs.get(qacc)
            qaccs.add(qacc)

            # Remove/replace reference only sections with the same
            # contig mapped to either side.  This addresses larger
            # insertion/deletion/replacement events that BLAST will
            # align as separate hits.
            if qacc in last_hits and last_section_from_reference:
                (_, lend), (rstart, _) = sorted((last_hits[qacc], (qstart, qend)))
                sequence[-1] = seq[lend : rstart - 1]

            # Reverse contig if blast hit was on negative strand
            if sstart > send:
                seq = revcomp(seq)
                qstart = len(seq) - qstart + 1

            # Extract sequence matching reference section and add to list
            offset = qstart + section_start - sstart - 1
            section_seqs.append(seq[offset : offset + section_len])

        if hits:
            last_hits = {qacc: (qstart, qend) for qacc, _, _, qstart, qend in hits}

        last_section_from_reference = False
        if len(section_seqs) == 0:
            # No contig covering this piece of reference.
            last_section_from_reference = True
            sequence.append(b"n" * section_len)
        elif len(section_seqs) == 1:
            ## Singleton piece - fill with sequence
            sequence.append(section_seqs[0])
        else:
            if any(len(section_seqs[0]) != len(seq) for seq in section_seqs):
                LOG.error("section sizes differ!")
                LOG.error(
                    "section: start=%i end=%i len=%i lens=%s",
                    section_start,
                    section_end,
                    section_len,
                    ", ".join(str(len(seq)) for seq in section_seqs),
                )
            ## Overlapping piece - fill with consensus
            section_consensus = []
            for base_counts in map(Counter, zip(*section_seqs)):
                if len(base_counts) == 1:
                    best, best_count = next(iter(base_counts.items()))
                    second_count = 0
                else:
                    (best, best_count), (_, second_count) = base_counts.most_common(2)

                if best_count == second_count:
                    section_consensus.append(110)  # 110 == 'n'
                else:
                    section_consensus.append(best)
            sequence.append(bytes(section_consensus))
    return {"+".join(qaccs): b"".join(sequence)}
