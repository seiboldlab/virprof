"""Calculate sequence entropy"""

from collections import Counter
from math import log2
import re

def sequence_entropy(sequence: bytes, klen: int = 5):
    """Calculates sequence entropy

    The base entropy formula is H(X) = sum(p(x) * log(p(x)) for x in X).

    For the purpose of this simple estimation, we assume that the
    states of X, the words of the sequence, are independent. That's
    untrue, but hopefully a sufficient approximation and also used
    commonly. Instead of using just a single basepair for each word,
    which would only look at how even the base distribution is, we
    consider k-mers.

    To make values with different k, and thus different number of
    categories more comparable, we divide by log(|categories|), which,
    since we actually use log_2 and assume 4 different bases, ends
    up being 2k (log_2 4^k == 2k).
    """

    # bail out if too short
    if len(sequence) < klen:
        return 0

    # do counting
    counts = Counter(
        sequence[i:i+klen]
        for i in range(len(sequence) - klen + 1)
    )

    # aggregate kmers containing N if any
    n_count = sum(n for kmer, n in counts.items() if b"N" in kmer)
    if n_count > 0:
        counts = {kmer:n for kmer,n in counts.items() if not b"N" in kmer}
        counts["N" * klen] = n_count
    words = sum(counts.values())
    logp = lambda x: x * log2(x)
    entropy = -sum(logp(c/words) for c in counts.values())
    return entropy / (2*klen)


def homopolymer_ratio(sequence: bytes, kmin: int = 5):
    """Calculates fraction of sequences composed of homopolymers

    A consecutive strech of identical bases of at least kmin length is
    considered a homopolymer.
    """

    # remove homopolymers
    # (regex will be faster than manual code in python)
    no_hp = re.sub(f"([A-Z])\\1{{{kmin},}}".encode(), b"", sequence)

    return 1-len(no_hp)/len(sequence)


def normalize_sequence(sequence: bytes):
    """Normalizes sequence for calculating entropy and hp count

    - remove n* sequences
    - uppercase
    - map wobbles to N
    """

    # remove scaffold padding 'n'* (assembled is N)
    sequence = re.sub(b"n+", b"", sequence)
    # make sure we have no lowercase left
    sequence = sequence.upper()
    # map unresolved IUPAC to N
    sequence = re.sub(b"[^AGCT]", b"n", sequence)

    return sequence
