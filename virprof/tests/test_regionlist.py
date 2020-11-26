import copy

from ..regionlist import RegionList


def test_empty():
    rl = RegionList()
    assert len(rl) == 0
    assert rl.get(0) == []


def test_nonoverlapping():
    """
    --- ---
    """
    rl = RegionList()
    rl.add(10, 20, "10-20")
    rl.add(30, 40, "30-40")
    assert rl.get(0)  == []
    assert rl.get(10) == ["10-20"]
    assert rl.get(20) == ["10-20"]
    assert rl.get(21) == []
    assert rl.get(29) == []
    assert rl.get(30) == ["30-40"]
    assert rl.get(40) == ["30-40"]
    assert rl.get(41) == []


def test_overlapping():
    """
    --- ---
    -------
    """
    rl = RegionList()
    rl.add(10, 20, "10-20")
    rl.add(30, 40, "30-40")
    rl.add(10, 40, "10-40")
    assert rl.get(15) == ["10-20", "10-40"]
    assert rl.get(25) == ["10-40"]
    assert rl.get(35) == ["30-40", "10-40"]


def test_sub():
    """
    --------
      ---
    """
    rl = RegionList()
    rl.add(10, 40, "10-40")
    rl.add(20, 30, "20-30")
    assert rl.get(15) == ["10-40"]
    assert rl.get(25) == ["10-40", "20-30"]
    assert rl.get(35) == ["10-40"]


def test_super():
    """
      ---
    --------
    """
    rl = RegionList()
    rl.add(20, 30, "20-30")
    rl.add(10, 40, "10-40")
    assert rl.get(15) == ["10-40"]
    assert rl.get(25) == ["20-30", "10-40"]
    assert rl.get(35) == ["10-40"]


def test_right_match():
    """
    ------
      ----
    """
    rl = RegionList()
    rl.add(10, 40, "10-40")
    rl.add(20, 40, "20-40")
    assert rl.get(15) == ["10-40"]
    assert rl.get(25) == ["10-40", "20-40"]
    assert rl.get(40) == ["10-40", "20-40"]
    assert rl.get(41) == []


def test_right_overlap():
    """
    ----
      ----
    """
    rl = RegionList()
    rl.add(10, 30, "10-30")
    rl.add(20, 40, "20-40")
    assert rl.get(15) == ["10-30"]
    assert rl.get(25) == ["10-30", "20-40"]
    assert rl.get(35) == ["20-40"]
    assert rl.get(41) == []


def test_left_overlap():
    """
      ----
    ----
    """
    rl = RegionList()
    rl.add(20, 40, "20-40")
    rl.add(10, 30, "10-30")
    assert rl.get(15) == ["10-30"]
    assert rl.get(25) == ["20-40", "10-30"]
    assert rl.get(35) == ["20-40"]
    assert rl.get(41) == []


def test_span():
    """
      --- ---
    -----------
    """
    rl = RegionList()
    rl.add(10, 20, "10-20")
    rl.add(30, 40, "30-40")
    rl.add(0, 50, "0-50")
    assert rl.get(-1) == []
    assert rl.get(0) == ["0-50"]
    assert rl.get(5) == ["0-50"]
    assert rl.get(10) == ["10-20", "0-50"]
    assert rl.get(20) == ["10-20", "0-50"]
    assert rl.get(21) == ["0-50"]
    assert rl.get(29) == ["0-50"]
    assert rl.get(30) == ["30-40", "0-50"]
    assert rl.get(40) == ["30-40", "0-50"]
    assert rl.get(41) == ["0-50"]
    assert rl.get(50) == ["0-50"]
    assert rl.get(51) == []


def test_intersect():
    """
    --- ----
      ---
    """
    rl = RegionList()
    rl.add(10, 20, "10-20")
    rl.add(30, 40, "30-40")
    rl.add(15, 35, "15-35")
    assert rl.get(9) == []
    assert rl.get(10) == ["10-20"]
    assert rl.get(14) == ["10-20"]
    assert rl.get(15) == ["10-20", "15-35"]
    assert rl.get(20) == ["10-20", "15-35"]
    assert rl.get(21) == ["15-35"]
    assert rl.get(29) == ["15-35"]
    assert rl.get(30) == ["30-40", "15-35"]
    assert rl.get(35) == ["30-40", "15-35"]
    assert rl.get(36) == ["30-40"]
    assert rl.get(40) == ["30-40"]
    assert rl.get(41) == []


def test_remove():
    rl = RegionList()
    rl.add(10, 20, "10-20")
    rl.add(30, 40, "30-40")
    rl.add(15, 35, "15-35")
    rl.remove(30, 40, "30-40")

    rl2 = RegionList()
    rl2.add(10, 20, "10-20")
    rl2.add(15, 35, "15-35")

    assert rl._region_starts == rl2._region_starts
    assert rl._region_data == rl2._region_data
    assert rl == rl2


def test_remove_all():
    rl = RegionList()
    rl.add(10, 20, "10-20")
    rl.add(30, 40, "30-40")
    rl.add(15, 35, "15-35")
    rl.remove(30, 40, "30-40")
    rl.remove(10, 20, "10-20")
    rl.remove(15, 35, "15-35")

    assert RegionList() == rl


def test_iter_empty():
    rl = RegionList()
    assert list(rl) == list()


def test_iter():
    """
    --- ----
      ---
    """
    rl = RegionList()
    rl.add(10, 20, "10-20")
    rl.add(30, 40, "30-40")
    rl.add(15, 35, "15-35")

    expected = [
        (10, 14, ["10-20"]),
        (15, 20, ["10-20", "15-35"]),
        (21, 29, ["15-35"]),
        (30, 35, ["30-40", "15-35"]),
        (36, 40, ["30-40"])
    ]

    assert expected == list(rl)


def test_reversed():
    rl = RegionList()
    rl2 = RegionList()
    rl.add(10, 20, "10-20")
    rl2.add(20, 10, "10-20")
    assert rl == rl2

    rl.add(30, 40, "30-40")
    rl2.add(30, 40, "30-40")
    assert rl == rl2

    rl.remove(20, 10, "10-20")
    rl2.remove(10, 20, "10-20")
    assert rl == rl2


def test_copy():
    rl = RegionList()
    rl.add(10, 20, "10-20")

    rl2 = copy.copy(rl)
    assert rl2 == rl

    rl2.add(30, 40, "30-40")
    assert rl2 != rl
    assert rl2.get(35) == ["30-40"]
    assert rl.get(35) == []
