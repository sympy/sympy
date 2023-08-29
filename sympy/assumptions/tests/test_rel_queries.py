from sympy.assumptions.lra_satask import lra_satask
from sympy.assumptions.ask import Q

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
