
from sympy import *

x, y = symbols('xy')

def test_erf():
    assert erf(nan) == nan

    assert erf(oo) == 1
    assert erf(-oo) == -1

    assert erf(0) == 0

    assert erf(-2) == -erf(2)
    assert erf(-x*y) == -erf(x*y)

def test_erf_series():
    assert erf(x).series(x, 7) == 2*x/sqrt(pi) - \
        2*x**3/3/sqrt(pi) + x**5/5/sqrt(pi) + O(x**7)
