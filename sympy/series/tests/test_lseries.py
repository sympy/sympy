from sympy import sin, cos, exp, E, S
from sympy.abc import x

def test_sin():
    e = sin(x).lseries(x)
    assert e.next() == x
    assert e.next() == -x**3/6
    assert e.next() == x**5/120

def test_cos():
    e = cos(x).lseries(x)
    assert e.next() == 1
    assert e.next() == -x**2/2
    assert e.next() == x**4/24

def test_exp():
    e = exp(x).lseries(x)
    assert e.next() == 1
    assert e.next() == x
    assert e.next() == x**2/2
    assert e.next() == x**3/6

def test_exp2():
    e = exp(cos(x)).lseries(x)
    assert e.next() == E
    assert e.next() == -E*x**2/2
    assert e.next() == E*x**4/6
    assert e.next() == -31*E*x**6/720

def test_simple():
    assert [t for t in x.lseries(x)] == [x]
    assert [t for t in S.One.lseries(x)] == [1]
