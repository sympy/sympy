from sympy import sin, cos, exp, E, S, Order
from sympy.abc import x, y

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
    assert [t for t in x.lseries()] == [x]
    assert [t for t in S.One.lseries(x)] == [1]
    assert not ((x/(x + y)).lseries(y).next()).has(Order)

def test_issue_2084():
    s = (x + 1/x).lseries()
    assert [si for si in s] == [1/x, x]
    assert (x + x**2).lseries().next() == x
    assert ((1+x)**7).lseries(x).next() == 1
    assert (sin(x + y)).series(x, n=3).lseries(y).next() == x
    # it would be nice if all terms were grouped, but in the
    # following case that would mean that all the terms would have
    # to be known since, for example, every term has a constant in it.
    s = ((1+x)**7).series(x, 1, n=None)
    assert [s.next() for i in range(2)] == [128, -448 + 448*x]
