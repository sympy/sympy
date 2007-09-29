from sympy import Symbol, exp, log
from sympy.series.limits2 import compare, mrv
from sympy.utilities.pytest import XFAIL

"""
This test suite is testing the limit algoritm using the bottom up approach.
See the documentation in limits2.py. The algorithm itself is highly recursive
by nature, so "compare" is logically the lowest part of the algorithm, yet in
some sense it's the most complex part, because it needs to calculate a limit to
return the result. 
"""

x = Symbol('x')
m = Symbol('m')

def test_compare1():
    """test_compare1, test_compare2, ... test compare, from easy cases to the
    most difficult"""

    assert compare(2, x, x) == "<"
    assert compare(x, exp(x), x) == "<"
    assert compare(exp(x), exp(x**2), x) == "<"
    assert compare(exp(x**2),exp(exp(x)), x) == "<"
    assert compare(1,exp(exp(x)), x) == "<"

    assert compare(x, 2, x) == ">"
    assert compare(exp(x), x, x) == ">"
    assert compare(exp(x**2), exp(x), x) == ">"
    assert compare(exp(exp(x)), exp(x**2), x) == ">"
    assert compare(exp(exp(x)), 1, x) == ">"

    assert compare(2, 3, x) == "="
    assert compare(3, -5, x) == "="
    assert compare(2, -5, x) == "="

    assert compare(x, x**2, x) == "="
    assert compare(x**2, x**3, x) == "="
    assert compare(x**3, 1/x, x) == "="
    assert compare(1/x, x**m, x) == "="
    assert compare(x**m, -x, x) == "="

    assert compare(exp(x), exp(-x), x) == "="
    assert compare(exp(-x), exp(2*x), x) == "="
    assert compare(exp(2*x), exp(x)**2, x) == "="
    assert compare(exp(x)**2, exp(x+exp(-x)), x) == "="
    assert compare(exp(x), exp(x+exp(-x)), x) == "="

    assert compare(exp(x**2), 1/exp(x**2), x) == "="

def test_mrv():
    assert mrv(x+1/x, x) == set([x])
    assert mrv(x**2, x) == set([x])
    assert mrv(log(x), x) == set([x])
    assert mrv(exp(x), x) == set([exp(x)])
    assert mrv(exp(x**2), x) == set([exp(x**2)])
    assert mrv(exp(x+1/x), x) == set([exp(x+1/x)])
    assert mrv(exp(x+exp(-exp(x))),x) == set([exp(-exp(x))])
