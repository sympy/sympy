from sympy import Symbol, exp, log, Mul, sqrt
from sympy.calculus.singularities import singularities
from sympy.abc import x, y

from sympy.utilities.pytest import raises, XFAIL


def test_singularities():

    # results should be sorted
    v = sqrt(3) - sqrt(7), sqrt(3)
    assert singularities(1/Mul(*[(x - i) for i in v])) == tuple(v)

    eq = 1/(x - 1) + 1/y
    raises(ValueError, lambda: singularities(eq))
    assert singularities(eq, y) == (0,)
    assert singularities(eq, x) == (1,)

    assert singularities(x**2, x) == ()
    assert singularities(x/(x**2 + 3*x + 2), x) == (-2, -1)
    assert singularities(1/(1/x - 1)) == (1,)
    assert singularities(1/(1/x - 1), strict=True) == (0, 1,)
    assert singularities(1/(exp(x) - exp(0)), poles=False) == (0,)
    assert singularities(exp(1/x), x, poles=False) == (0,)
    assert singularities(log((x - 2)**2), x, poles=False) == ()
    assert singularities(log(1 - x)/(x - 2), poles=False) == (2,)

@XFAIL
def test_singularity_poles():
    singularities(log(x)) == (0,)
