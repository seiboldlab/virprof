"""
Code from YMP

Parsers for blast output formats 6 (CSV) and 7 (CSV with comments between queries).
"""

import logging

from collections import namedtuple
from typing import List, NamedTuple, Iterator, Mapping, Callable

LOG = logging.getLogger(__name__)  # pylint: disable=invalid-name


def reader(fileobj, fmt: int = 7) -> "BlastParser":
    """
    Creates a reader for files in BLAST format

    >>> with open(blast_file) as infile:
    >>>    reader = blast.reader(infile)
    >>>    for hit in reader:
    >>>       print(hit)

    Args:
      fileobj: iterable yielding lines in blast format
      fmt: number of blast format type
    """
    if fmt == 7:
        return Fmt7Parser(fileobj)
    if fmt == 6:
        return Fmt6Parser(fileobj)
    raise NotImplementedError()


def writer(fileobj, fmt: int = 7) -> "BlastWriter":
    """
    Creates a writer for files in BLAST format

    >>> with open(blast_file) as outfile:
    >>>    writer = blast.writer(outfile)
    >>>    for hit in hits:
    >>>       writer.write_hit(hit)
    """
    if fmt == 7:
        return Fmt7Writer(fileobj)
    raise NotImplementedError()


class BlastBase:
    "Base class for BLAST readers and writers"

    @staticmethod
    def tupleofint(text):
        """Converts semicolon separated character string into tuple of int"""
        if text == "N/A":
            return tuple()
        try:
            return tuple(int(i) for i in text.split(";"))
        except ValueError:
            LOG.warning("Error parsing BLAST file at line='%s'", text)
            return tuple()

    #: Map between field short and long names
    FIELD_MAP = {
        # "": "qseqid",  # Query Seq-id
        # "": "qgi",  # Query GI
        "query acc.": "qacc",  # Query accession
        # "": "qaccver",  # Query accession.version
        "query length": "qlen",  # Query sequence length
        # "": "sseqid",  # Subject Seq-id
        # "": "sallseqid",  # All subject Seq-id(s), separated by ';'
        # "": "sgi",  # Subject GI
        # "": "sallgi",  # All subject GIs
        "subject acc.": "sacc",  # Subject accession
        # "": "saccver",  # Subject accession.version
        # "": "sallacc",  # All subject accessions
        # "": "slen",  # Subject sequence length
        "q. start": "qstart",  # Start of alignment in query
        "q. end": "qend",  # End of alignment in query
        "s. start": "sstart",  # Start of alignment in subject
        "s. end": "send",  # Start of alignment in query
        # "": "qseq",  # Aligned part of query sequence
        # "": "sseq",  # Aligned part of subject sequence
        "evalue": "evalue",  # Expect value
        "bit score": "bitscore",  # Bit score
        "score": "score",  # Raw score
        "alignment length": "length",  # Alignment length
        "% identity": "pident",  # Percentage of identical matches
        "mismatches": "mismatch",  # Number of mismatches
        # "": "positive", # Number of positive-scoring matches
        "gap opens": "gapopen",  # Number of gap openings
        # "": "gaps",  # Total number of gaps
        # "": "ppos",  # Percentage of positive-soring matches
        # "": "frames",  # Query and subject frames separated by a '/'
        "query frame": "qframe",  # Query frame
        "sbjct frame": "sframe",  # Subject frame
        "BTOP": "btop",  # Blast traceback operations (BTOP)
        # "": "staxid",  # Subject Taxonomy ID
        # "": "scciname",  # Subject Scientifi Name
        # "": "scomname",  # Subject Common Name
        # "": "sblastname",  # Subject Blast Name
        # "": "sskingdom",  # Subject Super Kingdom
        "subject tax ids": "staxids",  # sorted unique ';'-separated
                                       # Subject Taxonomy ID(s)
        # "": "sscinames",  # unique Subject Scientific Name(s)
        # "": "scomnames",  # unique Subject Common Name(s)
        # "": "sblastnames",  # unique Subject Blast Name(s)
        # "": "sskingdoms",  # unique Subject Super Kingdom(s)
        "subject title": "stitle",  # Subject Title
        # "": "sakktutkes",  # All Subject Title(s) separated by '<>'
        "subject strand": "sstrand",  # Subject Strand
        # "": "qcovs",  # Query Coverage per Subject
        # "": "qcovhsp",  # Query Coverage per HSP
        # "": "qcovus",  # Query Coverage per Unique Subject (blastn only)
    }

    #: Reversed map from short to long name
    FIELD_REV_MAP = {value: key for key, value in FIELD_MAP.items()}

    #: Map defining types of fields
    FIELD_TYPE: Mapping[str, Callable] = {
        "pident": float,
        "length": int,
        "mismatch": int,
        "gapopen": int,
        "qstart": int,
        "qend": int,
        "qlen": int,
        "sstart": int,
        "send": int,
        "evalue": float,
        "bitscore": float,
        "score": float,
        "sframe": int,
        "qframe": int,
        "stitle": str,
        "staxids": tupleofint,
    }


class BlastHit(NamedTuple):
    # pylint: disable=too-few-public-methods
    """Base type for a BLAST hit

    This class is only used for type checking.
    """
    qacc: str
    score: float
    sacc: str
    send: int
    sstart: int
    qstart: int
    qend: int
    qlen: int
    pident: float
    length: int
    stitle: str
    staxids: List[int]
    bitscore: int
    btop: str


class BlastParser(BlastBase):
    """Base class for BLAST readers"""

    def get_fields(self):
        raise NotImplementedError()

    def __iter__(self) -> Iterator[BlastHit]:
        raise NotImplementedError()


class BlastWriter(BlastBase):
    """Base class for BLAST writers"""

    def write_hit(self, hit: BlastHit):
        raise NotImplementedError()


class Fmt7Parser(BlastParser):
    """
    Parses BLAST results in format '7' (CSV with comments)
    """

    FIELDS = "# Fields: "
    QUERY = "# Query: "
    DATABASE = "# Database: "
    HITSFOUND = " hits found"

    def __init__(self, fileobj):
        self.fileobj = fileobj
        self.fields = None
        self.query = "undefined"
        self.database = "undefined"
        self.Hit = None
        if "BLAST" not in fileobj.readline():
            raise ValueError("not a BLAST7 formatted file")

    def get_fields(self) -> List[str]:
        """Returns list of available field names

        Format 7 specifies which columns it contains in comment lines, allowing
        this parser to be agnostic of the selection of columns made when running
        BLAST.

        Returns:
          List of field names (e.g. ``['sacc', 'qacc', 'evalue']``)
        """
        return self.fields

    def __iter__(self) -> Iterator[BlastHit]:
        for line in self.fileobj:
            if line.startswith(self.FIELDS):
                self.fields = [
                    self.FIELD_MAP[field] if field in self.FIELD_MAP else field
                    for field in line[len(self.FIELDS) :].strip().split(", ")
                ]
                self.Hit = namedtuple("BlastHit", self.fields)
            elif line.startswith(self.QUERY):
                self.query = line[len(self.QUERY) :].strip()
            elif line.startswith(self.DATABASE):
                self.database = line[len(self.DATABASE) :].strip()
            elif line.strip().endswith(self.HITSFOUND):
                self.hits = int(line.split()[1])
                self.hit = 0
            elif line[0] == "#":
                continue
            else:
                self.hit += 1
                yield self.Hit(
                    *[
                        self.FIELD_TYPE[key](value) if key in self.FIELD_TYPE else value
                        for key, value in zip(self.fields, line.strip().split("\t"))
                    ]
                )

    def isfirsthit(self) -> bool:
        """Returns `True` if the current hit is the first hit for the current
        query"""
        return self.hit == 1


class Fmt6Parser(BlastParser):
    """Parser for BLAST format 6 (CSV)"""

    #: Default field types
    fields = (
        "qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore"
    ).split()
    field_types = [BlastParser.FIELD_TYPE.get(n, None) for n in fields]
    Hit = namedtuple("BlastHit", fields)

    def __init__(self, fileobj):
        self.fileobj = fileobj

    def get_fields(self):
        return self.fields

    def __iter__(self):
        for line in self.fileobj:
            yield self.Hit(
                *[t(v) if t else v for v, t in zip(line.split("\t"), self.field_types)]
            )


class Fmt7Writer(BlastWriter):
    def __init__(self, fileobj):
        self.fileobj = fileobj
        self.toolname = "YMP writer"
        self.query = None
        self.database = "undefined"
        self.fields = "undefined"
        self.hits = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, ext_value, tb):
        self.write_hitset()

    def write_header(self):
        """Writes BLAST7 format header"""
        self.fileobj.write(
            f"# {self.toolname}\n"
            f"# Query: {self.query}\n"
            f"# Database: {self.database}\n"
            f"# Fields: {self.fields}\n"
        )

    def write_hitset(self):
        self.query = self.hits[0].qacc
        self.fields = self.hits[0]._fields
        self.write_header()
        self.fileobj.write(f"# {len(self.hits)} found")

    def write_hit(self, hit):
        if self.hits and hit.qacc != self.hits[0].qacc:
            self.write_hitset()
        self.hits.append(hit)
