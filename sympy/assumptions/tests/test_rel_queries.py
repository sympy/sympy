from sympy.assumptions.lra_satask import lra_satask
from sympy.assumptions.ask import Q, ask

from sympy.core import symbols, Symbol
from sympy.matrices.expressions.matexpr import MatrixSymbol
from sympy.core.numbers import I

from sympy.testing.pytest import raises

def test_lra_satask():
    x = symbols("x", real=True)
    im = Symbol('im', imaginary=True)

    # test preprocessing of unequalities is working correctly
    assert lra_satask(Q.eq(x, 1), ~Q.ne(x, 0)) is False
    assert lra_satask(Q.eq(x, 0), ~Q.ne(x, 0)) is True
    assert lra_satask(~Q.ne(x, 0), Q.eq(x, 0)) is True
    assert lra_satask(~Q.eq(x, 0), Q.eq(x, 0)) is False
    assert lra_satask(Q.ne(x, 0), Q.eq(x, 0)) is False

    assert lra_satask(Q.gt(x, 0), Q.gt(x, 1)) is True

    # check that True/False are handled
    assert lra_satask(Q.gt(x, 0), True) is None
    assert raises(ValueError, lambda: lra_satask(Q.gt(x, 0), False))

    # check imaginary numbers are correctly handled
    assert lra_satask(Q.gt(im * I, 0), Q.gt(im * I, 0)) is None  # (im * I).is_real returns True so this is an edge case

    # check matrix inputs
    X = MatrixSymbol("X", 2, 2)
    assert lra_satask(Q.lt(X, 2) & Q.gt(X, 3)) is None


def test_rel_queries():
    x, y, z = symbols("x y z", real=True)
    assert ask(Q.lt(x, 2) & Q.gt(x, 3)) is False
    assert ask(Q.positive(x - z), (x > y) & (y > z)) is True
    assert ask(x + y > 2, (x < 0) & (y <0)) is False
    assert ask(x > z, (x > y) & (y > z)) is True


def test_unhandled_queries():
    X = MatrixSymbol("X", 2, 2)
    assert ask(Q.lt(X, 2) & Q.gt(X, 3)) is None
