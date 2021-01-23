
from sympy import (ask, cos, Float, I, log, oo,
    pi, Q, Rational, sin, symbols)
from sympy.printing import sstr, pretty, latex

x,y,z = symbols('x y z')

def test_printing():
    eq = Q.eq(x,y)
    assert sstr(eq) == "x = y"
    assert pretty(eq) == "x = y"
    assert latex(eq) == "x = y"

def test_doit():
    assert Q.eq(x, 0).doit() == Q.eq(x, 0)
    assert Q.eq(1, 1).doit() == Q.eq(1, 1)

def test_ask():
    assert ask(Q.eq(0, 0)) is True
    assert ask(Q.eq(1, 0)) is False
    assert ask(Q.eq(I, 2)) is False
    a = Float('.000000000000000000001', '')
    b = Float('.0000000000000000000001', '')
    assert ask(Q.eq(pi+a, pi+b)) is False
    assert ask(Q.eq(x, x)) is True

    assert ask(Q.eq(log(cos(2)**2 + sin(2)**2), 0)) is True

    assert ask(Q.eq(x, 1), Q.negative(x)) is False

def test_reversed():
    assert Q.eq(x, y).reversed == Q.eq(y, x)

def test_issue_10304():
    d = cos(1)**2 + sin(1)**2 - 1
    assert d.is_comparable is False  # if this fails, find a new d
    e = 1 + d*I
    assert Q.eq(e, 0).simplify() == Q.eq(0, 1)
    assert ask(Q.eq(e, 0)) is False

def test_issue_18412():
    d = (Rational(1, 6) + z / 4 / y)
    assert Q.eq(x, pi * y**3 * d).replace(y**3, z) == Q.eq(x, pi * z * d)

# XXX : Migrate test_issue_10401 from test_relational after
# inequality handlers are implemented

def test_issue_10633():
    assert ask(Q.eq(True, False)) is False
    assert ask(Q.eq(False, True)) is False
    assert ask(Q.eq(True, True)) is True
    assert ask(Q.eq(False, False)) is True

def test_issue_10927():
    assert ask(Q.eq(x, oo)) is None
    assert ask(Q.eq(x, -oo)) is None
