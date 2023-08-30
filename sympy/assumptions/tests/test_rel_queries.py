from sympy.assumptions.lra_satask import lra_satask
from sympy.assumptions.ask import Q, ask

from sympy.core import symbols

def test_lra_satask():
    x = symbols("x", real=True)

    # test preprocessing of unequalities is working correctly
    assert lra_satask(Q.eq(x, 1), ~Q.ne(x, 0)) is False
    assert lra_satask(Q.eq(x, 0), ~Q.ne(x, 0)) is True
    assert lra_satask(~Q.ne(x, 0), Q.eq(x, 0)) is True
    assert lra_satask(~Q.eq(x, 0), Q.eq(x, 0)) is False
    assert lra_satask(Q.ne(x, 0), Q.eq(x, 0)) is False

    assert lra_satask(Q.gt(x, 0), Q.gt(x, 1)) is True

    # check that True/False are handled
    assert lra_satask(Q.gt(x, 0), True) is None
    assert lra_satask(Q.gt(x, 0), False) is None


def test_rel_queries():
    x, y, z = symbols("x y z", real=True)
    assert ask(Q.positive(x - z), (x > y) & (y > z)) is True

    assert ask(x + y > 2, (x < 0) & (y <0)) is False
    assert ask(x > z, (x > y) & (y > z)) is True
