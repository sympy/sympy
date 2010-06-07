from sympy import Rational, sqrt, symbols, sin, exp, log, sinh, cosh, cos, pi, \
    I, S, erf, tan, asin, asinh, acos, acosh, Function, Derivative, diff, simplify
from sympy.integrals.heurisch import heurisch, components
from sympy.utilities.pytest import XFAIL, skip

x, y, z = symbols('x,y,z')
f = Function('f')

def test_components():
    assert components(x*y, x) == set([x])
    assert components(1/(x+y), x) == set([x])
    assert components(sin(x), x) == set([sin(x), x])
    assert components(sin(x)*sqrt(log(x)), x) == \
       set([log(x), sin(x), sqrt(log(x)), x])
    assert components(x*sin(exp(x)*y), x) == \
       set([sin(y*exp(x)), x, exp(x)])
    assert components(x**Rational(17,54)/sqrt(sin(x)), x) == \
       set([sin(x), x**Rational(1,54), sqrt(sin(x)), x])

    assert components(f(x), x) == \
        set([x, f(x)])
    assert components(Derivative(f(x),x), x) == \
        set([x, f(x), Derivative(f(x), x)])
    assert components(f(x)*diff(f(x), x),  x) == \
        set([x, f(x), Derivative(f(x), x), Derivative(f(x), x)])

def test_heurisch_polynomials():
    assert heurisch(1, x) == x
    assert heurisch(x, x) == x**2/2
    assert heurisch(x**17, x) == x**18/18

def test_heurisch_fractions():
    assert heurisch(1/x, x) == log(x)
    assert heurisch(1/(2 + x), x) == log(x + 2)
    assert heurisch(1/(x+sin(y)), x) == log(x+sin(y))

    # Up to a constant, where C = 5*pi*I/12, Mathematica gives identical
    # result in the first case. The difference is because sympy changes
    # signs of expressions without any care.
    # XXX ^ ^ ^ is this still correct?
    assert heurisch(5*x**5/(2*x**6 - 5), x) in [5*log(2*x**6 - 5) / 12, 5*log(-2*x**6 + 5) / 12]
    assert heurisch(5*x**5/(2*x**6 + 5), x) == 5*log(2*x**6 + 5) / 12

    assert heurisch(1/x**2, x) == -1/x
    assert heurisch(-1/x**5, x) == 1/(4*x**4)

def test_heurisch_log():
    assert heurisch(log(x), x) == x*log(x) - x
    assert heurisch(log(3*x), x) == -x + x*log(3) + x*log(x)
    assert heurisch(log(x**2), x) in [x*log(x**2) - 2*x, 2*x*log(x) - 2*x]

def test_heurisch_exp():
    assert heurisch(exp(x), x) == exp(x)
    assert heurisch(exp(-x), x) == -exp(-x)
    assert heurisch(exp(17*x), x) == exp(17*x) / 17
    assert heurisch(x*exp(x), x) == x*exp(x) - exp(x)
    assert heurisch(x*exp(x**2), x) == exp(x**2) / 2

    assert heurisch(exp(-x**2), x) is None

    assert heurisch(2**x, x) == 2**x/log(2)
    assert heurisch(x*2**x, x) == x*2**x/log(2) - 2**x*log(2)**(-2)

def test_heurisch_trigonometric():
    assert heurisch(sin(x), x) == -cos(x)
    assert heurisch(pi*sin(x)+1, x) == x-pi*cos(x)

    assert heurisch(cos(x), x) == sin(x)
    assert heurisch(tan(x), x) in [
        log(1 + tan(x)**2)/2,
        log(tan(x) + I) + I*x,
        log(tan(x) - I) - I*x,
    ]

    assert heurisch(sin(x)*sin(y), x) == -cos(x)*sin(y)
    assert heurisch(sin(x)*sin(y), y) == -cos(y)*sin(x)

    # gives sin(x) in answer when run via setup.py and cos(x) when run via py.test
    assert heurisch(sin(x)*cos(x), x) in [sin(x)**2 / 2, -cos(x)**2 / 2]
    assert heurisch(cos(x)/sin(x), x) == log(sin(x))

    assert heurisch(x*sin(7*x), x) == sin(7*x) / 49 - x*cos(7*x) / 7
    assert heurisch(1/pi/4 * x**2*cos(x), x) == 1/pi/4*(x**2*sin(x) - 2*sin(x) + 2*x*cos(x))

    assert heurisch(acos(x/4) * asin(x/4), x) == 2*x - ((16-x**2)**Rational(1,2))*asin(x/4) \
        + ((16 - x**2)**Rational(1,2))*acos(x/4) + x*asin(x/4)*acos(x/4)

def test_heurisch_hyperbolic():
    assert heurisch(sinh(x), x) == cosh(x)
    assert heurisch(cosh(x), x) == sinh(x)

    assert heurisch(x*sinh(x), x) == x*cosh(x) - sinh(x)
    assert heurisch(x*cosh(x), x) == x*sinh(x) - cosh(x)

    assert heurisch(x*asinh(x/2), x) == x**2*asinh(x/2)/2 + asinh(x/2) - x*(4+x**2)**Rational(1,2)/4

def test_heurisch_mixed():
    assert heurisch(sin(x)*exp(x), x) == exp(x)*sin(x)/2 - exp(x)*cos(x)/2

def test_heurisch_radicals():
    assert heurisch(x**Rational(-1,2), x) == 2*x**Rational(1,2)
    assert heurisch(x**Rational(-3,2), x) == -2*x**Rational(-1,2)
    assert heurisch(x**Rational(3,2), x) == 2*x**Rational(5,2) / 5

    assert heurisch(sin(x)*sqrt(cos(x)), x) == -2*cos(x)**Rational(3,2) / 3
    assert heurisch(sin(y*sqrt(x)), x) == 2*y**(-2)*sin(y*x**S.Half) - \
                                          2*x**S.Half*cos(y*x**S.Half)/y

def test_heurisch_special():
    assert heurisch(erf(x), x) == x*erf(x) + exp(-x**2)/sqrt(pi)
    assert heurisch(exp(-x**2)*erf(x), x) == sqrt(pi)*erf(x)**2 / 4

def test_heurisch_symbolic_coeffs():
    assert heurisch(1/(x+y), x)         == log(x+y)
    assert heurisch(1/(x+sqrt(2)), x)   == log(x+sqrt(2))
    assert simplify(diff(heurisch(log(x+y+z), y), y)) == log(x+y+z)

def test_heurisch_symbolic_coeffs_1130():
    assert heurisch(1/(x**2+y), x) in [I*y**(-S.Half)*log(x + (-y)**S.Half)/2 - \
    I*y**(-S.Half)*log(x - (-y)**S.Half)/2, I*log(x + I*y**Rational(1,2)) / \
    (2*y**Rational(1,2)) - I*log(x - I*y**Rational(1,2))/(2*y**Rational(1,2))]

def test_heurisch_hacking():
    assert heurisch(sqrt(1 + 7*x**2), x, hints=[]) == \
        x*sqrt(1+7*x**2)/2 + sqrt(7)*asinh(sqrt(7)*x)/14
    assert heurisch(sqrt(1 - 7*x**2), x, hints=[]) == \
        x*sqrt(1-7*x**2)/2 + sqrt(7)*asin(sqrt(7)*x)/14

    assert heurisch(1/sqrt(1 + 7*x**2), x, hints=[]) == \
        sqrt(7)*asinh(sqrt(7)*x)/7
    assert heurisch(1/sqrt(1 - 7*x**2), x, hints=[]) == \
        sqrt(7)*asin(sqrt(7)*x)/7

    assert heurisch(exp(-7*x**2),x,hints=[]) == \
        sqrt(7*pi)*erf(sqrt(7)*x)/14

    assert heurisch(1/sqrt(9 - 4*x**2), x, hints=[]) == \
        asin(2*x/3)/2

    assert heurisch(1/sqrt(9 + 4*x**2), x, hints=[]) == \
        asinh(2*x/3)/2

def test_heurisch_function():
    df = diff(f(x), x)

    assert heurisch(f(x), x)            == None
    assert heurisch(f(x)*df, x)         == f(x)**2/2
    assert heurisch(f(x)**2 * df, x)    == f(x)**3/3
    assert heurisch(df / f(x), x)       == log(f(x))

def test_issue510():
    assert heurisch(1/(x * (1 + log(x)**2)), x) == I*log(log(x) + I)/2 - \
                                                   I*log(log(x) - I)/2

### These are examples from the Poor Man's Integrator
### http://www-sop.inria.fr/cafe/Manuel.Bronstein/pmint/examples/
#
# NB: correctness assured as ratsimp(diff(g,x) - f) == 0 in maxima
# SymPy is unable to do it :(

# Besides, they are skipped(), because they take too much time to execute.

@XFAIL
def test_pmint_rat():
    skip('takes too much time')
    f = (x**7-24*x**4-4*x**2+8*x-8) / (x**8+6*x**6+12*x**4+8*x**2)
    g = (4 + 8*x**2 + 6*x + 3*x**3) / (x*(x**4 + 4*x**2 + 4))  +  log(x)

    assert heurisch(f, x) == g


@XFAIL
def test_pmint_trig():
    skip('takes too much time')
    f = (x-tan(x)) / tan(x)**2  +  tan(x)
    g = (-x - tan(x)*x**2 / 2) / tan(x)  +  log(1+tan(x)**2) / 2

    assert heurisch(f, x) == g


@XFAIL
def test_pmint_logexp():
    skip('takes too much time')
    f = (1+x+x*exp(x))*(x+log(x)+exp(x)-1)/(x+log(x)+exp(x))**2/x
    g = 1/(x+log(x)+exp(x)) + log(x + log(x) + exp(x))

    assert heurisch(f, x) == g


@XFAIL
def test_pmint_erf():
    skip('takes too much time')
    f = exp(-x**2)*erf(x)/(erf(x)**3-erf(x)**2-erf(x)+1)
    g = sqrt(pi)/4 * (-1/(erf(x)-1) - log(erf(x)+1)/2 + log(erf(x)-1)/2)

    assert heurisch(f, x) == g


# TODO: convert the rest of PMINT tests:
# - Airy
# - Bessel
# - Whittaker
# - LambertW
# - Wright omega

