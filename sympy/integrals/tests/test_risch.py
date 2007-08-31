
from sympy import *

x = Symbol('x')

def test_risch_norman_polynomials():
    assert risch_norman(1, x) == x
    assert risch_norman(x, x) == x**2/2
    assert risch_norman(x**17, x) == x**18/18

def test_risch_norman_fractions():
    assert risch_norman(1/x, x) == log(x)
    assert risch_norman(1/(2+x), x) == log(x + 2)

    assert risch_norman(1/x**2, x) == -1/x
    assert risch_norman(-1/x**5, x) == 1/(4*x**4)

def test_risch_norman_log():
    assert risch_norman(log(x), x) == x*log(x) - x
    assert risch_norman(log(x**2), x) == x*log(x**2) - 2*x

def test_risch_norman_exp():
    assert risch_norman(exp(x), x) == exp(x)
    assert risch_norman(exp(-x), x) == -exp(-x)
    assert risch_norman(exp(17*x), x) == exp(17*x) / 17
    assert risch_norman(x*exp(x), x) == x*exp(x) - exp(x)
    assert risch_norman(x*exp(x**2), x) == exp(x**2) / 2

    # here will come trigs