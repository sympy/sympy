from sympy.mpmath import *

def test_diff():
    assert diff(cos, 1).ae(-sin(1))
    assert diff(abs, 0).ae(0)
    assert diff(abs, 0, 1).ae(1)
    assert diff(abs, 0, -1).ae(-1)

def test_diffc():
    assert diffc(exp, 1).ae(e)
    assert diffc(exp, 1, 5).ae(e)
