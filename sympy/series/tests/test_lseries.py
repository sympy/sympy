from sympy import sin, cos, exp, E, S, Order
from sympy.abc import x, y


def test_sin():
    e = sin(x).lseries(x)
    assert next(e) == x
    assert next(e) == -x**3/6
    assert next(e) == x**5/120


def test_cos():
    e = cos(x).lseries(x)
    assert next(e) == 1
    assert next(e) == -x**2/2
    assert next(e) == x**4/24


def test_exp():
    e = exp(x).lseries(x)
    assert next(e) == 1
    assert next(e) == x
    assert next(e) == x**2/2
    assert next(e) == x**3/6


def test_exp2():
    e = exp(cos(x)).lseries(x)
    assert next(e) == E
    assert next(e) == -E*x**2/2
    assert next(e) == E*x**4/6
    assert next(e) == -31*E*x**6/720


def test_simple():
    assert [t for t in x.lseries()] == [x]
    assert [t for t in S.One.lseries(x)] == [1]
    assert not next((x/(x + y)).lseries(y)).has(Order)


def test_issue_2084():
    s = (x + 1/x).lseries()
    assert [si for si in s] == [1/x, x]
    assert next((x + x**2).lseries()) == x
    assert next(((1 + x)**7).lseries(x)) == 1
    assert next((sin(x + y)).series(x, n=3).lseries(y)) == x
    # it would be nice if all terms were grouped, but in the
    # following case that would mean that all the terms would have
    # to be known since, for example, every term has a constant in it.
    s = ((1 + x)**7).series(x, 1, n=None)
    assert [next(s) for i in range(2)] == [128, -448 + 448*x]
