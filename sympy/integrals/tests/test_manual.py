from sympy import (sin, cos, tan, sec, csc, cot, log, exp, atan, asin, acos,
                   Symbol, Integral, integrate, pi, Dummy, Derivative,
                   diff, I, sqrt, erf, Piecewise, Eq, symbols,
                   And, Heaviside, S, asinh, acosh)
from sympy.integrals.manualintegrate import manualintegrate, find_substitutions, \
    _parts_rule

x, y, u, n, a, b = symbols('x y u n a b')


def test_find_substitutions():
    assert find_substitutions((cot(x)**2 + 1)**2*csc(x)**2*cot(x)**2, x, u) == \
        [(cot(x), 1, -u**6 - 2*u**4 - u**2)]
    assert find_substitutions((sec(x)**2 + tan(x) * sec(x)) / (sec(x) + tan(x)),
                              x, u) == [(sec(x) + tan(x), 1, 1/u)]
    assert find_substitutions(x * exp(-x**2), x, u) == [(-x**2, -S.Half, exp(u))]


def test_manualintegrate_polynomials():
    assert manualintegrate(y, x) == x*y
    assert manualintegrate(exp(2), x) == x * exp(2)
    assert manualintegrate(x**2, x) == x**3 / 3
    assert manualintegrate(3 * x**2 + 4 * x**3, x) == x**3 + x**4

    assert manualintegrate((x + 2)**3, x) == (x + 2)**4 / 4
    assert manualintegrate((3*x + 4)**2, x) == (3*x + 4)**3 / 9

    assert manualintegrate((u + 2)**3, u) == (u + 2)**4 / 4
    assert manualintegrate((3*u + 4)**2, u) == (3*u + 4)**3 / 9


def test_manualintegrate_exponentials():
    assert manualintegrate(exp(2*x), x) == exp(2*x) / 2
    assert manualintegrate(2**x, x) == (2 ** x) / log(2)

    assert manualintegrate(1 / x, x) == log(x)
    assert manualintegrate(1 / (2*x + 3), x) == log(2*x + 3) / 2
    assert manualintegrate(log(x)**2 / x, x) == log(x)**3 / 3


def test_manualintegrate_parts():
    assert manualintegrate(exp(x) * sin(x), x) == \
        (exp(x) * sin(x)) / 2 - (exp(x) * cos(x)) / 2
    assert manualintegrate(2*x*cos(x), x) == 2*x*sin(x) + 2*cos(x)
    assert manualintegrate(x * log(x), x) == x**2*log(x)/2 - x**2/4
    assert manualintegrate(log(x), x) == x * log(x) - x
    assert manualintegrate((3*x**2 + 5) * exp(x), x) == \
        -6*x*exp(x) + (3*x**2 + 5)*exp(x) + 6*exp(x)
    assert manualintegrate(atan(x), x) == x*atan(x) - log(x**2 + 1)/2

    # Make sure _parts_rule doesn't pick u = constant but can pick dv =
    # constant if necessary, e.g. for integrate(atan(x))
    assert _parts_rule(cos(x), x) == None
    assert _parts_rule(exp(x), x) == None
    assert _parts_rule(x**2, x) == None
    result = _parts_rule(atan(x), x)
    assert result[0] == atan(x) and result[1] == 1


def test_manualintegrate_trigonometry():
    assert manualintegrate(sin(x), x) == -cos(x)
    assert manualintegrate(tan(x), x) == -log(cos(x))

    assert manualintegrate(sec(x), x) == log(sec(x) + tan(x))
    assert manualintegrate(csc(x), x) == -log(csc(x) + cot(x))

    assert manualintegrate(sin(x) * cos(x), x) in [sin(x) ** 2 / 2, -cos(x)**2 / 2]
    assert manualintegrate(-sec(x) * tan(x), x) == -sec(x)
    assert manualintegrate(csc(x) * cot(x), x) == -csc(x)
    assert manualintegrate(sec(x)**2, x) == tan(x)
    assert manualintegrate(csc(x)**2, x) == -cot(x)

    assert manualintegrate(x * sec(x**2), x) == log(tan(x**2) + sec(x**2))/2
    assert manualintegrate(cos(x)*csc(sin(x)), x) == -log(cot(sin(x)) + csc(sin(x)))


def test_manualintegrate_trigpowers():
    assert manualintegrate(sin(x)**2 * cos(x), x) == sin(x)**3 / 3
    assert manualintegrate(sin(x)**2 * cos(x) **2, x) == \
        x / 8 - sin(4*x) / 32
    assert manualintegrate(sin(x) * cos(x)**3, x) == -cos(x)**4 / 4
    assert manualintegrate(sin(x)**3 * cos(x)**2, x) == \
        cos(x)**5 / 5 - cos(x)**3 / 3

    assert manualintegrate(tan(x)**3 * sec(x), x) == sec(x)**3/3 - sec(x)
    assert manualintegrate(tan(x) * sec(x) **2, x) == sec(x)**2/2

    assert manualintegrate(cot(x)**5 * csc(x), x) == \
        -csc(x)**5/5 + 2*csc(x)**3/3 - csc(x)
    assert manualintegrate(cot(x)**2 * csc(x)**6, x) == \
        -cot(x)**7/7 - 2*cot(x)**5/5 - cot(x)**3/3


def test_manualintegrate_inversetrig():
    # atan
    assert manualintegrate(exp(x) / (1 + exp(2*x)), x) == atan(exp(x))
    assert manualintegrate(1 / (4 + 9 * x**2), x) == atan(3 * x/2) / 6
    assert manualintegrate(1 / (16 + 16 * x**2), x) == atan(x) / 16
    assert manualintegrate(1 / (4 + x**2), x) == atan(x / 2) / 2
    assert manualintegrate(1 / (1 + 4 * x**2), x) == atan(2*x) / 2
    assert manualintegrate(1/(a + b*x**2), x) == \
        Piecewise(((sqrt(a/b)*atan(x*sqrt(b/a))/a), And(a > 0, b > 0)))
    assert manualintegrate(1/(4 + b*x**2), x) == \
        Piecewise((sqrt(1/b)*atan(sqrt(b)*x/2)/2, b > 0))
    assert manualintegrate(1/(a + 4*x**2), x) == \
        Piecewise((atan(2*x*sqrt(1/a))/(2*sqrt(a)), a > 0))
    assert manualintegrate(1/(4 + 4*x**2), x) == atan(x) / 4

    # asin
    assert manualintegrate(1/sqrt(1-x**2), x) == asin(x)
    assert manualintegrate(1/sqrt(4-4*x**2), x) == asin(x)/2
    assert manualintegrate(3/sqrt(1-9*x**2), x) == asin(3*x)
    assert manualintegrate(1/sqrt(4-9*x**2), x) == asin(3*x/2)/3

    # asinh
    assert manualintegrate(1/sqrt(x**2 + 1), x) == \
        asinh(x)
    assert manualintegrate(1/sqrt(x**2 + 4), x) == \
        asinh(x/2)
    assert manualintegrate(1/sqrt(4*x**2 + 4), x) == \
        asinh(x)/2
    assert manualintegrate(1/sqrt(4*x**2 + 1), x) == \
        asinh(2*x)/2
    assert manualintegrate(1/sqrt(a*x**2 + 1), x) == \
        Piecewise((sqrt(-1/a)*asin(x*sqrt(-a)), a < 0), (sqrt(1/a)*asinh(sqrt(a)*x), a > 0))
    assert manualintegrate(1/sqrt(a + x**2), x) == \
        Piecewise((asinh(x*sqrt(1/a)), a > 0), (acosh(x*sqrt(-1/a)), a < 0))

    # acosh
    assert manualintegrate(1/sqrt(x**2 - 1), x) == \
        acosh(x)
    assert manualintegrate(1/sqrt(x**2 - 4), x) == \
        acosh(x/2)
    assert manualintegrate(1/sqrt(4*x**2 - 4), x) == \
        acosh(x)/2
    assert manualintegrate(1/sqrt(9*x**2 - 1), x) == \
        acosh(3*x)/3
    assert manualintegrate(1/sqrt(a*x**2 - 4), x) == \
        Piecewise((sqrt(1/a)*acosh(sqrt(a)*x/2), a > 0))
    assert manualintegrate(1/sqrt(-a + 4*x**2), x) == \
        Piecewise((asinh(2*x*sqrt(-1/a))/2, -a > 0), (acosh(2*x*sqrt(1/a))/2, -a < 0))

    # piecewise
    assert manualintegrate(1/sqrt(a-b*x**2), x) == \
        Piecewise((sqrt(a/b)*asin(x*sqrt(b/a))/sqrt(a), And(-b < 0, a > 0)),
                  (sqrt(-a/b)*asinh(x*sqrt(-b/a))/sqrt(a), And(-b > 0, a > 0)),
                  (sqrt(a/b)*acosh(x*sqrt(b/a))/sqrt(-a), And(-b > 0, a < 0)))
    assert manualintegrate(1/sqrt(a + b*x**2), x) == \
        Piecewise((sqrt(-a/b)*asin(x*sqrt(-b/a))/sqrt(a), And(a > 0, b < 0)),
                  (sqrt(a/b)*asinh(x*sqrt(b/a))/sqrt(a), And(a > 0, b > 0)),
                  (sqrt(-a/b)*acosh(x*sqrt(-b/a))/sqrt(-a), And(a < 0, b > 0)))


def test_manualintegrate_trig_substitution():
    assert manualintegrate(sqrt(16*x**2 - 9)/x, x) == \
        Piecewise((sqrt(16*x**2 - 9) - 3*acos(3/(4*x)),
                   And(x < 3*S.One/4, x > -3*S.One/4)))
    assert manualintegrate(1/(x**4 * sqrt(25-x**2)), x) == \
        Piecewise((-sqrt(-x**2/25 + 1)/(125*x) -
                   (-x**2/25 + 1)**(3*S.Half)/(15*x**3), And(x < 5, x > -5)))
    assert manualintegrate(x**7/(49*x**2 + 1)**(3 * S.Half), x) == \
        ((49*x**2 + 1)**(5*S.Half)/28824005 -
         (49*x**2 + 1)**(3*S.Half)/5764801 +
         3*sqrt(49*x**2 + 1)/5764801 + 1/(5764801*sqrt(49*x**2 + 1)))


def test_manualintegrate_rational():
    assert manualintegrate(1/(4 - x**2), x) == -log(x - 2)/4 + log(x + 2)/4
    assert manualintegrate(1/(-1 + x**2), x) == log(x - 1)/2 - log(x + 1)/2


def test_manualintegrate_derivative():
    assert manualintegrate(pi * Derivative(x**2 + 2*x + 3), x) == \
        pi * ((x**2 + 2*x + 3))
    assert manualintegrate(Derivative(x**2 + 2*x + 3, y), x) == \
        x * Derivative(x**2 + 2*x + 3, y)
    assert manualintegrate(Derivative(sin(x), x, x, y, x), x) == \
        Derivative(sin(x), x, x, y)


def test_manualintegrate_Heaviside():
    assert manualintegrate(Heaviside(x), x) == x*Heaviside(x)
    assert manualintegrate(x*Heaviside(2), x) == x**2/2
    assert manualintegrate(x*Heaviside(-2), x) == 0
    assert manualintegrate(x*Heaviside( x), x) == x**2*Heaviside( x)/2
    assert manualintegrate(x*Heaviside(-x), x) == x**2*Heaviside(-x)/2
    assert manualintegrate(Heaviside(2*x + 4), x) == (x+2)*Heaviside(2*x + 4)
    assert manualintegrate(x*Heaviside(x), x) == x**2*Heaviside(x)/2
    assert manualintegrate(Heaviside(x + 1)*Heaviside(1 - x)*x**2, x) == \
        ((x**3/3 + S(1)/3)*Heaviside(x + 1) - S(2)/3)*Heaviside(-x + 1)

    y = Symbol('y')
    assert manualintegrate(sin(7 + x)*Heaviside(3*x - 7), x) == \
            (- cos(x + 7) + cos(S(28)/3))*Heaviside(3*x - S(7))

    assert manualintegrate(sin(y + x)*Heaviside(3*x - y), x) == \
            (cos(4*y/3) - cos(x + y))*Heaviside(3*x - y)


def test_issue_6799():
    r, x, phi = map(Symbol, 'r x phi'.split())
    n = Symbol('n', integer=True, positive=True)

    integrand = (cos(n*(x-phi))*cos(n*x))
    limits = (x, -pi, pi)
    assert manualintegrate(integrand, x).has(Integral)
    assert r * integrate(integrand.expand(trig=True), limits) / pi == r * cos(n * phi)
    assert not integrate(integrand, limits).has(Dummy)


def test_issue_3796():
    assert manualintegrate(diff(exp(x + x**2)), x) == exp(x + x**2)
    assert integrate(x * exp(x**4), x, risch=False) == -I*sqrt(pi)*erf(I*x**2)/4


def test_manual_true():
    assert integrate(exp(x) * sin(x), x, manual=True) == \
        (exp(x) * sin(x)) / 2 - (exp(x) * cos(x)) / 2
    assert integrate(sin(x) * cos(x), x, manual=True) in \
        [sin(x) ** 2 / 2, -cos(x)**2 / 2]


def test_issue_6746():
    y = Symbol('y')
    n = Symbol('n')
    assert manualintegrate(y**x, x) == \
        Piecewise((x, Eq(log(y), 0)), (y**x/log(y), True))
    assert manualintegrate(y**(n*x), x) == \
        Piecewise(
            (x, Eq(n, 0)),
            (Piecewise(
                (n*x, Eq(log(y), 0)),
                (y**(n*x)/log(y), True))/n, True))
    assert manualintegrate(exp(n*x), x) == \
        Piecewise((x, Eq(n, 0)), (exp(n*x)/n, True))

    y = Symbol('y', positive=True)
    assert manualintegrate((y + 1)**x, x) == (y + 1)**x/log(y + 1)
    y = Symbol('y', zero=True)
    assert manualintegrate((y + 1)**x, x) == x
    y = Symbol('y')
    n = Symbol('n', nonzero=True)
    assert manualintegrate(y**(n*x), x) == \
        Piecewise((n*x, Eq(log(y), 0)), (y**(n*x)/log(y), True))/n
    y = Symbol('y', positive=True)
    assert manualintegrate((y + 1)**(n*x), x) == \
        (y + 1)**(n*x)/(n*log(y + 1))
    a = Symbol('a', negative=True)
    assert manualintegrate(1 / (a + b*x**2), x) == \
        Integral(1/(a + b*x**2), x)


def test_issue_2850():
    assert manualintegrate(asin(x)*log(x), x) == -x*asin(x) - sqrt(-x**2 + 1) \
            + (x*asin(x) + sqrt(-x**2 + 1))*log(x) - Integral(sqrt(-x**2 + 1)/x, x)
    assert manualintegrate(acos(x)*log(x), x) == -x*acos(x) + sqrt(-x**2 + 1) + \
        (x*acos(x) - sqrt(-x**2 + 1))*log(x) + Integral(sqrt(-x**2 + 1)/x, x)
    assert manualintegrate(atan(x)*log(x), x) == -x*atan(x) + (x*atan(x) - \
            log(x**2 + 1)/2)*log(x) + log(x**2 + 1)/2 + Integral(log(x**2 + 1)/x, x)/2


def test_constant_independent_of_symbol():
    assert manualintegrate(Integral(y, (x, 1, 2)), x) == x*Integral(y, (x, 1, 2))
