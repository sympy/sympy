from sympy.assumptions.lra_satask import lra_satask
from sympy.assumptions.ask import Q

from sympy.core import symbols

def test_lra_satask():
    x = symbols("x", finite=True)
    assert lra_satask(Q.gt(x, 1), Q.gt(x, 0)) is True
