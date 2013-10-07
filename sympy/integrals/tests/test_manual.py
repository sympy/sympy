from sympy import (sin, cos, tan, sec, csc, cot, log, exp, atan,
                   Symbol, Mul, Integral, integrate, pi, Dummy,
                   Derivative, diff, I, sqrt, erf, Piecewise,
                   Eq, Ne, Q, assuming, symbols, And)
from sympy.integrals.manualintegrate import manualintegrate, find_substitutions, \
    integral_steps, _parts_rule

x, y, u, n, a, b = symbols('x y u n a b')

def test_find_substitutions():
    assert find_substitutions((cot(x)**2 + 1)**2*csc(x)**2*cot(x)**2, x, u) == \
        [(cot(x), 1, -u**6 - 2*u**4 - u**2)]
    assert find_substitutions((sec(x)**2 + tan(x) * sec(x)) / (sec(x) + tan(x)),
                              x, u) == [(sec(x) + tan(x), 1, 1/u)]

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
    assert manualintegrate(exp(x) / (1 + exp(2*x)), x) == atan(exp(x))
    assert manualintegrate(1 / (4 + 9 * x**2), x) == atan(3 * x/2) / 6
    assert manualintegrate(1 / (16 + 16 * x**2), x) == atan(x) / 16
    assert manualintegrate(1 / (4 + x**2), x) == atan(x / 2) / 2
    assert manualintegrate(1 / (1 + 4 * x**2), x) == atan(2*x) / 2
    assert manualintegrate(1/(a + b*x**2), x) == \
        Piecewise((atan(x*sqrt(b/a))/(a*sqrt(b/a)), And(a > 0, b > 0)))
    assert manualintegrate(1/(4 + b*x**2), x) == \
        Piecewise((atan(sqrt(b)*x/2)/(2*sqrt(b)), b > 0))
    assert manualintegrate(1/(a + 4*x**2), x) == \
        Piecewise((atan(2*x*sqrt(1/a))/(2*a*sqrt(1/a)), a > 0))
    assert manualintegrate(1/(4 + 4*x**2), x) == atan(x) / 4

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

def test_issue_3700():
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

def test_issue_3647():
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

    with assuming(~Q.zero(log(y))):
        assert manualintegrate(y**x, x) == y**x/log(y)
    with assuming(Q.zero(log(y))):
        assert manualintegrate(y**x, x) == x
    with assuming(~Q.zero(n)):
        assert manualintegrate(y**(n*x), x) == \
            Piecewise((n*x, Eq(log(y), 0)), (y**(n*x)/log(y), True))/n
    with assuming(~Q.zero(n) & ~Q.zero(log(y))):
        assert manualintegrate(y**(n*x), x) == \
            y**(n*x)/(n*log(y))
    with assuming(Q.negative(a)):
        assert manualintegrate(1 / (a + b*x**2), x) == \
            Integral(1/(a + b*x**2), x)
