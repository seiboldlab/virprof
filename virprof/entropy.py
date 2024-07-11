"""Calculate sequence entropy"""

from collections import Counter
from math import log2
import re

def sequence_entropy(sequence: bytes, klen: int = 5):
    """Calculates sequence entropy

    The entropy formula is H(X) = sum(p(x) * log(p(x)) for x in X).

    For the purpose of this simple estimation, we assume that the
    states of X, the words of the sequence, are independent. That's
    untrue, but hopefully a sufficient approximation and also used
    commonly. Instead of using just a single basepair for each word,
    which would only look at how even the base distribution is, we
    consider k-mers.
    """

    # bail out if too short
    if len(sequence) < klen:
        return 0

    # remove 'n'*
    sequence = re.sub(b"n+", b"", sequence)
    # ignore case
    sequence = sequence.upper()

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
    return -sum(logp(c/words) for c in counts.values())


def homopolymer_ratio(sequence: bytes, kmin: int = 5):
    """Calculates fraction of sequences composed of homopolymers

    A consecutive strech of identical bases of at least kmin length is
    considered a homopolymer.
    """
    # remove 'n'*
    sequence = re.sub(b"n+", b"", sequence)

    # ignore case
    sequence = sequence.upper()

    # remove homopolymers
    # (regex will be faster than manual code in python)
    #print(f"([A-Z])\\1{{{kmin},}}")
    no_hp = re.sub(f"([A-Z])\\1{{{kmin},}}".encode(), b"", sequence)

    return 1-len(no_hp)/len(sequence)
