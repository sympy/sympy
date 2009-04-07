from sympy import sin, cos, exp, E, series
from sympy.abc import x

def test_sin():
    e1 = sin(x).series(x, 0)
    e2 = series(sin(x), x, 0)
    assert e1 == e2

def test_cos():
    e1 = cos(x).series(x, 0)
    e2 = series(cos(x), x, 0)
    assert e1 == e2

def test_exp():
    e1 = exp(x).series(x, 0)
    e2 = series(exp(x), x, 0)
    assert e1 == e2

def test_exp2():
    e1 = exp(cos(x)).series(x, 0)
    e2 = series(exp(cos(x)),x,0)
    assert e1 == e2
