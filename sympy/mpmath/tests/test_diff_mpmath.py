from sympy.mpmath import *

def test_diff():
    assert diff(log, 2.0, n=0).ae(log(2))
    assert diff(cos, 1.0).ae(-sin(1))
    assert diff(abs, 0.0) == 0
    assert diff(abs, 0.0, direction=1) == 1
    assert diff(abs, 0.0, direction=-1) == -1
    assert diff(exp, 1.0).ae(e)
    assert diff(exp, 1.0, n=5).ae(e)
    assert diff(exp, 2.0, n=5, direction=3*j).ae(e**2)
    assert diff(lambda x: x**2, 3.0, method='quad').ae(6)
    assert diff(lambda x: 3+x**5, 3.0, n=2, method='quad').ae(540)
    assert diff(lambda x: 3+x**5, 3.0, n=2, method='step').ae(540)
    assert diffun(sin)(2).ae(cos(2))
    assert diffun(sin, n=2)(2).ae(-sin(2))

def test_taylor():
    # Easy to test since the coefficients are exact in floating-point
    assert taylor(sqrt, 1, 4) == [1, 0.5, -0.125, 0.0625, -0.0390625]
