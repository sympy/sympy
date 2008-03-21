from sympy.functions.elementary.miscellaneous import min_, max_

def test_min():
    assert min_(5, 4) == 4

def test_max():
    assert max_(5, 4) == 5
