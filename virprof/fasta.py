"""Helper methods for fasta"""

import subprocess as sp

from typing import Iterator, Sequence, Collection, BinaryIO


def read_from_command(args: Sequence[str]) -> Iterator[str]:
    """Runs external command yielding output lines"""
    proc = sp.Popen(args, stdout=sp.PIPE)
    while True:
        line = proc.stdout.readline()
        if not line:
            break
        yield line
    proc.stdout.close()
    proc.wait()


def get_accs_from_fasta(fileobj: BinaryIO) -> Iterator[str]:
    """Reads accession numbers from (gzipped) FASTA file"""
    for line in read_from_command(["zgrep", "^>", fileobj.name]):
        if line:
            yield line[1:].split(maxsplit=1)[0]


def filter_fasta(filein: BinaryIO, fileout: BinaryIO,
                 accs: Collection[str], remove: bool) -> None:
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
    skip = True
    fasta_header = b'>'[0]
    accs_b = set(acc.encode('ascii') for acc in accs)
    for line in unzip:
        if line[0] == fasta_header:
            acc = line[1:].split(maxsplit=1)[0]
            if remove:
                skip = acc in accs_b
            else:
                skip = acc not in accs_b
        if not skip:
            outzip.stdin.write(line)
    outzip.stdin.close()
    outzip.wait()


class FastaFile:
    """Handles access to GZipp'ed FASTA format file

    """
    def __init__(self, iofile: BinaryIO, mode='r') -> None:
        self.iofile = iofile
        self.mode = mode

        if 'r' in mode:
            self.sequences = self._load_all()
        else:
            self.sequences = None

        if 'w' in mode:
            self.outzip = sp.Popen(["gzip", "-c"], stdout=iofile, stdin=sp.PIPE)
        else:
            self.outzip = None

    def close(self) -> None:
        """Close potentially open file handles"""
        if self.outzip is not None:
            self.outzip.stdin.close()
            self.outzip.wait()

    def __len__(self) -> int:
        return len(self.sequences)

    def _load_all(self) -> dict:
        sequences = dict()
        fastafile = read_from_command(["gunzip", "-dc", self.iofile.name])
        fasta_header = b'>'[0]
        acc = None
        lines = []
        for line in fastafile:
            if line[0] == fasta_header:
                if acc is not None:
                    sequences[acc] = b''.join(lines)
                acc = line[1:].split(maxsplit=1)[0]
            else:
                lines.append(line.strip())
        sequences[acc] = b''.join(lines)
        return sequences

    def get(self, acc: str, start: int = 1, stop: int = None) -> bytes:
        """Retrieve a sequence or a part of a sequence

        Params:
          acc: Unique sequence identifier (first word on header line)
          start: Start position of subsequence (1 indexed)
          stop: End position of subsequence (``None`` means until the end)
        """
        seq = self.sequences[acc.encode('utf-8')]
        if stop is None:
            stop = len(seq)
        return seq[start-1:stop]

    def put(self, acc: str, sequence: bytes, comment: str = None):
        """Write a sequence

        Params:
          acc: Unique sequence identifier
          sequence: Sequence data
          comment: Additional data to add to seuqence header
        """
        if comment is not None:
            header = ">{} {}".format(acc, comment).encode('utf-8')
        else:
            header = ">{}".format(acc).encode('utf-8')
        self.outzip.stdin.write(b"\n".join((header, sequence, b"")))
