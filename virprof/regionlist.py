"""Implements RegionList"""

from typing import List, Any, Tuple, Generator


class RegionList:
    """A container of ranges with associated data

    >>> rl = RegionList()
    >>> rl.add(10, 20, "ten-to-twenty")
    >>> rl.get(15)
    ['ten-to-twenty']
    """

    def __init__(self) -> None:
        self._region_starts: List[int] = []
        self._region_data: List[List] = []

    def __len__(self) -> int:
        return 0 if not self._region_starts else len(self._region_starts) - 1

    def __str__(self) -> str:
        return f"{self.__class__.__name__}: starts={self._region_starts} data={self._region_data}"

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        return (
            self._region_starts == other._region_starts
            and self._region_data == other._region_data
        )

    def __iter__(self) -> Generator[Tuple[int, int, Any], None, None]:
        zipped = zip(self._region_starts, self._region_data)
        try:
            start, data = next(zipped)
        except StopIteration:
            return
        for nextstart, nextdata in zipped:
            yield start, nextstart - 1, data
            start = nextstart
            data = nextdata

    def __copy__(self):
        cpy = object.__new__(type(self))
        cpy._region_starts = list(self._region_starts)
        cpy._region_data = list(self._region_data)
        return cpy

    def _find(self, start: int) -> int:
        """Find index of region containing start

        Returns -1 if before any region
        """
        if not self._region_starts:
            return -1
        for idx, ostart in enumerate(self._region_starts):
            if start < ostart:
                return idx - 1
        return len(self._region_starts) - 1

    def get(self, position) -> List:
        """Access data for ``position``"""
        idx = self._find(position)
        try:
            return self._region_data[idx]
        except IndexError:
            return []

    def _ensure_region(self, start: int) -> int:
        """Finds or creates region for ``start``"""

        index = self._find(start)
        try:
            oldstart = self._region_starts[index]
        except IndexError:
            # The new start point was out of bounds
            if index == -1:
                # Before everything, need to shift before inserting
                index = 0
            self._region_starts.insert(index, start)
            # Initialize with empty list
            self._region_data.insert(index, [])
            return index

        if oldstart != start:
            # We are in the middle of a region, copy it and return the
            data = list(self._region_data[index])
            index += 1
            self._region_starts.insert(index, start)
            self._region_data.insert(index, data)
        return index

    def add(self, start: int, stop: int, data) -> None:
        """Adds a new region ranging from ``stop`` to ``start`` inclusively"""
        if stop < start:
            stop, start = start, stop
        stop += 1  # make 'stop' exclusive
        start_idx = self._ensure_region(start)
        stop_idx = self._ensure_region(stop)
        for idx in range(start_idx, stop_idx):
            self._region_data[idx].append(data)

    def remove(self, start: int, stop: int, data) -> None:
        """Removes a region data entry"""
        if stop < start:
            stop, start = start, stop
        start_idx = self._find(start)
        stop_idx = self._find(stop)
        for idx in range(start_idx, stop_idx + 1):
            self._region_data[idx].remove(data)
        dups = set()
        for idx in range(1, len(self._region_data)):
            if self._region_data[idx - 1] == self._region_data[idx]:
                dups.add(idx)
        if self._region_data[0] == []:
            dups.add(0)
        if dups:
            self._region_data = [
                self._region_data[i]
                for i in range(len(self._region_data))
                if i not in dups
            ]
            self._region_starts = [
                self._region_starts[i]
                for i in range(len(self._region_starts))
                if i not in dups
            ]

    def end_pos(self) -> int:
        """Retrieve last position covered by a region"""
        return self._region_starts[-1] - 1
