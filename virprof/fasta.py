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
    def __init__(self, fname: str) -> None:
        self.fname = fname
        self.sequences = dict()
        self._load_all()
        

    def _load_all(self) -> None:
        fastafile = read_from_command(["gunzip", "-dc", self.fname])
        fasta_header = b'>'[0]
        acc = None
        lines = []
        for line in fastafile:
            if line[0] == fasta_header:
                if acc is not None:
                    self.sequences[acc] = b''.join(lines)
                acc = line[1:].split(maxsplit=1)[0]
            else:
                lines.append(line.strip())
        self.sequences[acc] = b''.join(lines)
    
    def get(self, acc: str, start: int, stop: int) -> str:
        seq = self.sequences[acc.encode('utf-8')]
        return seq[start-1:stop]
