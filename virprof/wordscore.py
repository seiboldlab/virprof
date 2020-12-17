"""This module handles summarizing blast hit description lines"""

import string

from collections import defaultdict
from typing import Optional, Set, Iterable, List, Dict

from .blastbin import HitChain


class WordScorer:
    """Creates best-scoring-word summary from sequence titles

    This class takes a set of `HitChain`s and generates a string
    from the ``stitle`` fields containing the sequence description
    using word order, frequency and HitChain score.

    Args:
      stopwords: List of words to omit from output
      maxwordlen: Maximum length of words to include in output
      keepwords: Maximum number of words to include in output
    """

    def __init__(self,
                 stopwords: Optional[Set[str]] = None,
                 maxwordlen: int = 27,
                 orderweight: int = 2,
                 keepwords: int = 4) -> None:
        if stopwords is None:
            self.stopwords = set([
                'genome', 'complete', 'sequence', 'strain',
                'isolate', 'and', 'or', 'sapiens', 'genes',
                'predicted', 'human', 'homo', 'assembly',
                'mus', 'musculus', 'virus', 'polyprotein'
            ])
        else:
            self.stopwords = stopwords
        self.maxwordlen = maxwordlen
        self.orderweight = orderweight
        self.keepwords = keepwords

    @staticmethod
    def split_words(chain: HitChain) -> Iterable[str]:
        """Split title of ``chain`` into words"""
        return chain.stitle.split()

    @staticmethod
    def strip_punctuation(words: Iterable[str]) -> Iterable[str]:
        """Remove punctuation from the outside of each word

        Args:
          words: sequence of words
        Returns:
          Generator over words
        """
        for word in words:
            yield word.strip(string.punctuation)

    @staticmethod
    def uncapitalize(words: Iterable[str]) -> Iterable[str]:
        """Remove capitalization for words without internal caps

        Args:
          words: sequence of words
        Returns:
          Generator over words
        """
        for word in words:
            if word[1:].islower():
                yield word.lower()
            else:
                yield word

    @staticmethod
    def combine_shortwords(words: Iterable[str]) -> Iterable[str]:
        """Combine short words with their predecessor

        Args:
          words: sequence of words
        Returns:
          Generator over words
        """
        worditer = iter(words)
        try:
            stack = [next(worditer)]
        except StopIteration:
            return
        for word in words:
            if len(word) <= 2:
                stack += [' '.join((stack[-1], word))]
            else:
                yield from reversed(stack)
                stack = [word]
        yield from reversed(stack)

    def filter_stopwords(self, words: Iterable[str]) -> Iterable[str]:
        """Filter out words occurring in our stopword list

        Args:
          words: sequence of words
        Returns:
          Generator over words
        """
        for word in words:
            if word.lower() not in self.stopwords:
                yield word

    def score(self, chains: List[HitChain]) -> str:
        """Generate summary using word scoring

        Breaks the title of each passed in HitChain into words,
        removes outside punctuation, removes captitalization, merges
        short words with their preceding word(s) and removes
        (singleton) stopwords.

        Each word occurring first in any input title receives a score
        equal to the sum of the log-evalues (exponent) of the
        respective HitChain. The words and scores are added with each
        subsequent word, dividing the log-evalue by ``orderweight`` as
        we move through the title.

        The resulting words are sorted by score and emitted as
        result. Up to ``keepwords`` words may comprise the
        result. Words composed from multiple words are split for this
        accounting, and no (sub)word may occur twice.
        """
        word_scores: Dict[str, float] = defaultdict(float)
        for chain in chains:
            words = self.split_words(chain)
            words = self.strip_punctuation(words)
            words = self.uncapitalize(words)
            words = self.combine_shortwords(words)
            words = self.filter_stopwords(words)

            for idx, word in enumerate(words):
                word_scores[word] += (
                    chain.log10_evalue / self.orderweight ** idx
                )

        sorted_words = sorted(word_scores.items(), key=lambda k: k[1])

        result_words = []
        seen = set()
        for word, _score in sorted_words:
            parts = [part for part in word.split()
                     if part not in seen]
            seen.update(parts)
            for part in parts:
                result_words.append(part)
            if len(result_words) > self.keepwords:
                break
        return " ".join(result_words[:self.keepwords])

    @staticmethod
    def score_taxids(chains: List[HitChain]) -> int:
        """Selects best scoring taxid from chains

        Each input HitChain may have multiple NCBI taxonomy IDs
        assigned (NCBI taxonomy may assign multiple to a single
        reference sequence). For each found ID, the respective
        log-evalues are summed and the best score selected.
        """
        taxid_scores: Dict[int, float] = defaultdict(float)
        for chain in chains:
            for taxid in chain.staxids:
                taxid_scores[taxid] += chain.log10_evalue
                # FIXME! This is bad for chains with many hits
        staxids = [
            taxid for taxid, score
            in sorted(taxid_scores.items(), key=lambda k: k[1])
        ]
        return staxids[0]
