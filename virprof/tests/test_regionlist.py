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

    
