from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.trigonometric import cos, sin
from sympy.abc import x
from sympy.series.pade_approximant import pade_approximant
from sympy.simplify import simplify

def test_sin_52():
    assert pade_approximant(sin(x), x, 0, 5, 2) == \
        (x - x**3/7 + (11*x**5)/2520)/(1 + x**2/42)

def test_exp_22():
    assert pade_approximant(exp(x), x, 0, 2, 2) == \
        (x**2/12 + x/2 + 1)/(x**2/12 - x/2 + 1)

def test_cos_24():
    assert pade_approximant(cos(x), x, 0, 2, 4) == \
        (1 - 61*x**2/150)/(x**4/200 + 7*x**2/75 + 1)

def test_exp_cos_33():
    assert pade_approximant(exp(cos(x)), x, 0, 3, 3) == \
        (exp(1) - exp(1)*x**2/6)/(x**2/3 + 1)

def test_log_32():
    """
    The below assertion evaluates to False, despite being correct.
    assert pade_approximant(log(1 + x), x, 1, 3, 2)\
      == ((x - 1)**3/240 + (x - 1)**2*(3*log(2)/40 + 7/40)\
      + (x - 1)*(3*log(2)/5 + 1/2) + log(2))/(3*x/5 + 3*(x - 1)**2/40 + 2/5)
    This assertion is mathematically equivalent to the above """

    assert 0==\
        simplify(pade_approximant(log(1 + x), x, 1, 3, 2)\
            - ((x - 1)**3/240 + (x - 1)**2*(3*log(2)/40 + 7/40)\
                + (x - 1)*(3*log(2)/5 + 1/2) + log(2))\
                /(3*x/5 + 3*(x - 1)**2/40 + 2/5))

def test_rational_02():
    assert pade_approximant(1/(1 + x**2), x, 0, 0, 2) == 1/(1 + x**2)

def test_rational_20():
    assert pade_approximant(2 + x**2, x, 0, 2, 0) == 2 + x**2
