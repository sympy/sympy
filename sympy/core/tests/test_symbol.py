from sympy import Symbol, Inequality, StrictInequality
import py


def test_Symbol():
    a = Symbol("a")
    x1 = Symbol("x")
    x2 = Symbol("x")
    xdummy1 = Symbol("x", dummy=True)
    xdummy2 = Symbol("x", dummy=True)

    assert a != x1
    assert a != x2
    assert x1 == x2
    assert x1 != xdummy1
    assert xdummy1 != xdummy2

    assert Symbol("x") == Symbol("x")
    assert Symbol("x", dummy=True) != Symbol("x", dummy=True)


def test_lt_gt():
    x, y = Symbol('x'), Symbol('y')

    assert (x <= y) == Inequality(x, y)
    assert (x >= y) == Inequality(y, x)
    assert (x <= 0) == Inequality(x, 0)
    assert (x >= 0) == Inequality(0, x)

    assert (x < y) == StrictInequality(x, y)
    assert (x > y) == StrictInequality(y, x)
    assert (x < 0) == StrictInequality(x, 0)
    assert (x > 0) == StrictInequality(0, x)

    assert (x**2+4*x+1 > 0) == StrictInequality(0, x**2+4*x+1)
