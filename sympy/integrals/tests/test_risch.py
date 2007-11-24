
from sympy import *
from sympy.utilities.pytest import XFAIL

x, y = symbols('xy')

def test_risch_norman_polynomials():
    assert risch_norman(1, x) == x
    assert risch_norman(x, x) == x**2/2
    assert risch_norman(x**17, x) == x**18/18

def test_risch_norman_fractions():
    assert risch_norman(1/x, x) == log(x)
    assert risch_norman(1/(2 + x), x) == log(x + 2)

    assert risch_norman(5*x**5/(2*x**6 + 5), x) == 5*log(5 + 2*x**6) / 12

    assert risch_norman(1/x**2, x) == -1/x
    assert risch_norman(-1/x**5, x) == 1/(4*x**4)

def test_risch_norman_log():
    assert risch_norman(log(x), x) == x*log(x) - x
    assert risch_norman(log(3*x), x) == x*log(3*x) - x
    assert risch_norman(log(x**2), x) in [x*log(x**2) - 2*x, 2*x*log(x) - 2*x]

def test_risch_norman_exp():
    assert risch_norman(exp(x), x) == exp(x)
    assert risch_norman(exp(-x), x) == -exp(-x)
    assert risch_norman(exp(17*x), x) == exp(17*x) / 17
    assert risch_norman(x*exp(x), x) == x*exp(x) - exp(x)
    assert risch_norman(x*exp(x**2), x) == exp(x**2) / 2

    assert risch_norman(exp(-x**2), x) is None

def test_risch_norman_trigonometric():
    assert risch_norman(sin(x), x) == -cos(x)
    assert risch_norman(cos(x), x) == sin(x)

    assert risch_norman(sin(x)*sin(y), x) == -cos(x)*sin(y)
    assert risch_norman(sin(x)*sin(y), y) == -cos(y)*sin(x)

    assert risch_norman(sin(x)*cos(x), x) == sin(x)**2 / 2
    assert risch_norman(cos(x)/sin(x), x) == log(sin(x))

    assert risch_norman(x*sin(7*x), x) == sin(7*x) / 49 - x*cos(7*x) / 7
    assert risch_norman(x**2*cos(x), x) == x**2*sin(x) - 2*sin(x) + 2*x*cos(x)

def test_risch_norman_hyperbolic():
    assert risch_norman(sinh(x), x) == cosh(x)
    assert risch_norman(cosh(x), x) == sinh(x)

    assert risch_norman(x*sinh(x), x) == x*cosh(x) - sinh(x)
    assert risch_norman(x*cosh(x), x) == x*sinh(x) - cosh(x)

def test_risch_norman_mixed():
    assert risch_norman(sin(x)*exp(x), x) == exp(x)*sin(x)/2 - exp(x)*cos(x)/2

def test_risch_norman_special():
    assert risch_norman(erf(x), x) == x*erf(x) + exp(-x**2)/sqrt(pi)
    assert risch_norman(exp(-x**2)*erf(x), x) == sqrt(pi)*erf(x)**2 / 4

def test_risch_norman_issue442():
    assert risch_norman(1/(x+y), x)         == log(x+y)
#   assert risch_norman(1/(x**2+y), x)      == ?


@XFAIL
def test_risch_norman_issue442_0():
    assert risch_norman(1/(x+sqrt(2)), x)   == log(x+sqrt(2))

@XFAIL
def test_resch_norman_issue442_1():
    assert risch_norman(1/(x+sin(y)), x)    == log(x+sin(y))



### These are examples from the Poor Man's Integrator
### http://www-sop.inria.fr/cafe/Manuel.Bronstein/pmint/examples/
#
# NB: correctness assured as ratsimp(diff(g,x) - f) == 0 in maxima
# SymPy is unable to do it :(

@XFAIL
def test_pmint_rat():
    f = (x**7-24*x**4-4*x**2+8*x-8) / (x**8+6*x**6+12*x**4+8*x**2)
    g = (4 + 8*x**2 + 6*x + 3*x**3) / (x*(x**4 + 4*x**2 + 4))  +  log(x)

    assert risch_norman(f, x) == g


@XFAIL
def test_pmint_trig():
    f = (x-tan(x)) / tan(x)**2  +  tan(x)
    g = (-x - tan(x)*x**2 / 2) / tan(x)  +  log(1+tan(x)**2) / 2

    assert risch_norman(f, x) == g


@XFAIL
def test_pmint_logexp():
    f = (1+x+x*exp(x))*(x+log(x)+exp(x)-1)/(x+log(x)+exp(x))**2/x
    g = 1/(x+log(x)+exp(x)) + log(x + log(x) + exp(x))

    assert risch_norman(f, x) == g


@XFAIL
def test_pmint_erf():
    f = exp(-x**2)*erf(x)/(erf(x)**3-erf(x)**2-erf(x)+1)
    g = sqrt(pi)/4 * (-1/(erf(x)-1) - log(erf(x)+1)/2 + log(erf(x)-1)/2)

    assert risch_norman(f, x) == g


# TODO: convert the rest of PMINT tests:
# - Airy
# - Bessel
# - Whittaker
# - LambertW
# - Wright omega

