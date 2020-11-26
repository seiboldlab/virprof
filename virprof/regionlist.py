from typing import List

class RegionList:
    def __init__(self) -> None:
        self._region_starts: List[int] = []
        self._region_data: List[List] = []

    def __len__(self) -> int:
        return len(self._region_starts)
        
    def __str__(self) -> str:
        return f"{self.__class__.__name__}: starts={self._region_starts} data={self._region_data}"

    def _find(self, start: int) -> int:
        """Find index of region containing start
        
        Returns -1 if before any region
        """
        if not self._region_starts:
            return -1
        for idx, ostart in enumerate(self._region_starts):
            if start < ostart:
                return idx - 1
        return len(self._region_starts)

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
        stop += 1  # make 'stop' exclusive
        start_idx = self._ensure_region(start)
        stop_idx = self._ensure_region(stop)
        for idx in range(start_idx, stop_idx):
            self._region_data[idx].append(data)
       
            
