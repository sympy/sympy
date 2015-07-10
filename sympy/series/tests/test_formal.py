from sympy import (symbols, factorial, sqrt, Rational, atan, I, log, fps, O,
                   Sum, oo, S, pi, cos, sin, Function, exp, Derivative, asin,
                   airyai, acos, acosh, gamma)
from sympy.series.formal import (rational_algorithm, FormalPowerSeries,
                                 rational_independent, simpleDE, exp_re,
                                 hyper_re)
from sympy.utilities.pytest import raises, XFAIL

x, y, z = symbols('x y z')
n, m, k = symbols('n m k')
f, r = Function('f'), Function('r')


def test_rational_algorithm():
    f = 1 / ((x - 1)**2 * (x - 2))
    assert rational_algorithm(f, x, k) == \
        (-2**(-k - 1) + 1 - (factorial(k + 1) / factorial(k)), 0, 0)

    f = (1 + x + x**2 + x**3) / ((x - 1) * (x - 2))
    assert rational_algorithm(f, x, k) == \
        (-15*2**(-k - 1) + 4, x + 4, 0)

    f = z / (y*m - m*x - y*x + x**2)
    assert rational_algorithm(f, x, k) == \
        (((-y**(-k - 1)*z) / (y - m)) + ((m**(-k - 1)*z) / (y - m)), 0, 0)

    f = x / (1 - x - x**2)
    assert rational_algorithm(f, x, k) is None
    assert rational_algorithm(f, x, k, full=True) == \
        (((-Rational(1, 2) + sqrt(5)/2)**(-k - 1) *
         (-sqrt(5)/10 + Rational(1, 2))) +
         ((-sqrt(5)/2 - Rational(1, 2))**(-k - 1) *
         (sqrt(5)/10 + Rational(1, 2))), 0, 0)

    f = 1 / (x**2 + 2*x + 2)
    assert rational_algorithm(f, x, k) is None
    assert rational_algorithm(f, x, k, full=True) == \
        ((I*(-1 + I)**(-k - 1)) / 2 - (I*(-1 - I)**(-k - 1)) / 2, 0, 0)

    f = log(1 + x)
    assert rational_algorithm(f, x, k) == \
        (-(-1)**(-k) / k, 0, 1)

    f = atan(x)
    assert rational_algorithm(f, x, k) is None
    assert rational_algorithm(f, x, k, full=True) == \
        (((I*I**(-k)) / 2 - (I*(-I)**(-k)) / 2) / k, 0, 1)

    f = x*atan(x) - log(1 + x**2) / 2
    assert rational_algorithm(f, x, k) is None
    assert rational_algorithm(f, x, k, full=True) == \
        (((I*I**(-k + 1)) / 2 - (I*(-I)**(-k + 1)) / 2) /
         (k*(k - 1)), 0, 2)

    f = log((1 + x) / (1 - x)) / 2 - atan(x)
    assert rational_algorithm(f, x, k) is None
    assert rational_algorithm(f, x, k, full=True) == \
        ((-(-1)**(-k) / 2 - (I*I**(-k)) / 2 + (I*(-I)**(-k)) / 2 +
          Rational(1, 2)) / k, 0, 1)

    assert rational_algorithm(cos(x), x, k) is None


def test_rational_independent():
    ri = rational_independent
    assert ri([], x) == []
    assert ri([cos(x), sin(x)], x) == [cos(x), sin(x)]
    assert ri([x**2, sin(x), x*sin(x), x**3], x) == \
        [x**3 + x**2, x*sin(x) + sin(x)]
    assert ri([S.One, x*log(x), log(x), sin(x)/x, cos(x), sin(x), x], x) == \
        [x + 1, x*log(x) + log(x), sin(x)/x + sin(x), cos(x)]


def test_simpleDE():
    assert simpleDE(exp(x), x, f) == -f(x) + Derivative(f(x), x)
    assert simpleDE(sin(x), x, f) == f(x) + Derivative(f(x), x, x)
    assert simpleDE(log(1 + x), x, f) == \
        (x + 1)*Derivative(f(x), x, 2) + Derivative(f(x), x)
    assert simpleDE(asin(x), x, f) == \
        x*Derivative(f(x), x) + (x**2 - 1)*Derivative(f(x), x, x)
    assert simpleDE(exp(x)*sin(x), x, f) == \
        2*f(x) - 2*Derivative(f(x)) + Derivative(f(x), x, x)
    assert simpleDE(((1 + x)/(1 - x))**n, x, f) == \
        2*n*f(x) + (x**2 - 1)*Derivative(f(x), x)
    assert simpleDE(airyai(x), x, f) == -x*f(x) + Derivative(f(x), x, x)


def test_exp_re():
    d = -f(x) + Derivative(f(x), x)
    assert exp_re(d, r, k) == -r(k) + r(k + 1)

    d = f(x) + Derivative(f(x), x, x)
    assert exp_re(d, r, k) == r(k) + r(k + 2)

    d = f(x) + Derivative(f(x), x) + Derivative(f(x), x, x)
    assert exp_re(d, r, k) == r(k) + r(k + 1) + r(k + 2)

    d = Derivative(f(x), x) + Derivative(f(x), x, x)
    assert exp_re(d, r, k) == r(k) + r(k + 1)

    d = Derivative(f(x), x, 3) + Derivative(f(x), x, 4) + Derivative(f(x))
    assert exp_re(d, r, k) == r(k) + r(k + 2) + r(k + 3)


def test_hyper_re():
    d = f(x) + Derivative(f(x), x, x)
    assert hyper_re(d, r, k) == r(k) + (k+1)*(k+2)*r(k + 2)

    d = -x*f(x) + Derivative(f(x), x, x)
    assert hyper_re(d, r, k) == (k + 2)*(k + 3)*r(k + 3) - r(k)

    d = 2*f(x) - 2*Derivative(f(x), x) + Derivative(f(x), x, x)
    assert hyper_re(d, r, k) == \
        (-2*k - 2)*r(k + 1) + (k + 1)*(k + 2)*r(k + 2) + 2*r(k)

    d = 2*n*f(x) + (x**2 - 1)*Derivative(f(x), x)
    assert hyper_re(d, r, k) == \
        k*r(k) + 2*n*r(k + 1) + (-k - 2)*r(k + 2)

    d = (x**10 + 4)*Derivative(f(x), x) + x*(x**10 - 1)*Derivative(f(x), x, x)
    assert hyper_re(d, r, k) == \
        (k*(k - 1) + k)*r(k) + (4*k - (k + 9)*(k + 10) + 40)*r(k + 10)

    d = ((x**2 - 1)*Derivative(f(x), x, 3) + 3*x*Derivative(f(x), x, x) +
         Derivative(f(x), x))
    assert hyper_re(d, r, k) == \
        ((k*(k - 2)*(k - 1) + 3*k*(k - 1) + k)*r(k) +
         (-k*(k + 1)*(k + 2))*r(k + 2))


def test_fps():
    assert fps(1) == 1
    assert fps(2, x) == 2
    assert fps(2, x, dir='+') == 2
    assert fps(2, x, dir='-') == 2
    assert fps(x**2 + x + 1) == x**2 + x + 1
    assert fps(1/x + 1/x**2) == 1/x + 1/x**2
    assert fps(log(1 + x), hyper=False, rational=False) == log(1 + x)

    f = fps(log(1 + x))
    assert isinstance(f, FormalPowerSeries)
    assert f.function == log(1 + x)
    assert f.subs(x, y) == f
    assert f[:5] == [0, x, -x**2/2, x**3/3, -x**4/4]
    assert f.as_leading_term(x) == x
    assert f.polynomial(6) == x - x**2/2 + x**3/3 - x**4/4 + x**5/5

    k = f.ak.variables[0]
    assert f.infinite == Sum((-(-1)**(-k)*x**k)/k, (k, 1, oo))

    ft, s = f.truncate(n=None), f[:5]
    for i, t in enumerate(ft):
        if i == 5:
            break
        assert s[i] == t

    raises(NotImplementedError, lambda: fps(y*x))
    raises(ValueError, lambda: fps(x, dir=0))


def test_fps__rational():
    assert fps(1/x) == (1/x)
    assert fps((x**2 + x + 1) / x**3, dir=-1) == (x**2 + x + 1) / x**3

    f = 1 / ((x - 1)**2 * (x - 2))
    assert fps(f, x).truncate() == \
        (-Rational(1, 2) - 5*x/4 - 17*x**2/8 - 49*x**3/16 - 129*x**4/32
         - 321*x**5/64 + O(x**6))

    f = (1 + x + x**2 + x**3) / ((x - 1) * (x - 2))
    assert fps(f, x).truncate() == \
        (Rational(1, 2) + 5*x/4 + 17*x**2/8 + 49*x**3/16 + 113*x**4/32
         + 241*x**5/64 + O(x**6))

    f = x / (1 - x - x**2)
    assert fps(f, x, full=True).truncate() == \
        x + x**2 + 2*x**3 + 3*x**4 + 5*x**5 + O(x**6)

    f = 1 / (x**2 + 2*x + 2)
    assert fps(f, x, full=True).truncate() == \
        Rational(1, 2) - x/2 + x**2/4 - x**4/8 + x**5/8 + O(x**6)

    f = log(1 + x)
    assert fps(f, x).truncate() == \
        x - x**2/2 + x**3/3 - x**4/4 + x**5/5 + O(x**6)
    assert fps(f, x, dir=1).truncate() == fps(f, x, dir=-1).truncate()
    assert fps(f, x, 2).truncate() == \
        (log(3) - Rational(2, 3) - (x - 2)**2/18 + (x - 2)**3/81
         - (x - 2)**4/324 + (x - 2)**5/1215 + x/3 + O((x - 2)**6, (x, 2)))
    assert fps(f, x, 2, dir=-1).truncate() == \
        (log(3) - Rational(2, 3) - (-x + 2)**2/18 - (-x + 2)**3/81
         - (-x + 2)**4/324 - (-x + 2)**5/1215 + x/3 + O((x - 2)**6, (x, 2)))

    f = atan(x)
    assert fps(f, x, full=True).truncate() == x - x**3/3 + x**5/5 + O(x**6)
    assert fps(f, x, full=True, dir=1).truncate() == \
        fps(f, x, full=True, dir=-1).truncate()
    assert fps(f, x, 2, full=True).truncate() == \
        (atan(2) - Rational(2, 5) - 2*(x - 2)**2/25 + 11*(x - 2)**3/375
         - 6*(x - 2)**4/625 + 41*(x - 2)**5/15625 + x/5 + O((x - 2)**6, (x, 2)))
    assert fps(f, x, 2, full=True, dir=-1).truncate() == \
        (atan(2) - Rational(2, 5) - 2*(-x + 2)**2/25 - 11*(-x + 2)**3/375
         - 6*(-x + 2)**4/625 - 41*(-x + 2)**5/15625 + x/5 + O((x - 2)**6,
                                                              (x, 2)))
    assert fps(f, x, oo, full=True).truncate() == \
        -1/(5*x**5) + 1/(3*x**3) - 1/x + pi/2 + O(1/x**6, (x, oo))
    assert fps(f, x, -oo, full=True).truncate() == \
        -1/(5*x**5) + 1/(3*x**3) - 1/x - pi/2 + O(1/x**6, (x, oo))

    f = x*atan(x) - log(1 + x**2) / 2
    assert fps(f, x, full=True).truncate() == x**2/2 - x**4/12 + O(x**6)

    f = log((1 + x) / (1 - x)) / 2 - atan(x)
    assert fps(f, x, full=True).truncate(n=10) == 2*x**3/3 + 2*x**7/7 + O(x**10)


def test_fps__hyper():
    f = sin(x)
    assert fps(f, x).truncate() == x - x**3/6 + x**5/120 + O(x**6)

    f = cos(x)
    assert fps(f, x).truncate() == 1 - x**2/2 + x**4/24 + O(x**6)

    f = exp(x)
    assert fps(f, x).truncate() == \
        1 + x + x**2/2 + x**3/6 + x**4/24 + x**5/120 + O(x**6)

    f = atan(x)
    assert fps(f, x).truncate() == x - x**3/3 + x**5/5 + O(x**6)

    f = exp(acos(x))
    assert fps(f, x).truncate() == \
        (exp(pi/2) - x*exp(pi/2) + x**2*exp(pi/2)/2 - x**3*exp(pi/2)/3 +
         5*x**4*exp(pi/2)/24 - x**5*exp(pi/2)/6 + O(x**6))

    f = exp(acosh(x))
    assert fps(f, x).truncate() == I + x - I*x**2/2 - I*x**4/8 + O(x**6)

    f = atan(1/x)
    assert fps(f, x).truncate() == pi/2 - x + x**3/3 - x**5/5 + O(x**6)

    f = x*atan(x) - log(1 + x**2) / 2
    assert fps(f, x, rational=False).truncate() == x**2/2 - x**4/12 + O(x**6)

    f = log(1 + x)
    assert fps(f, x, rational=False).truncate() == \
        x - x**2/2 + x**3/3 - x**4/4 + x**5/5 + O(x**6)

    f = x*exp(x)*sin(2*x)  # TODO: solved using rsolve, improve simpleDE
    assert fps(f, x).truncate() == 2*x**2 + 2*x**3 - x**4/3 - x**5 + O(x**6)

    f = airyai(x**2)
    assert fps(f, x).truncate() == \
        (3**Rational(5, 6)*gamma(Rational(1, 3))/(6*pi) -
         3**Rational(2, 3)*x**2/(3*gamma(Rational(1, 3))) + O(x**6))

    f = exp(x)*sin(x)
    assert fps(f, x).truncate() == x + x**2 + x**3/3 - x**5/30 + O(x**6)


def test_fps_shift():
    f = x**-5*sin(x)
    assert fps(f, x).truncate() == \
        1/x**4 - 1/(6*x**2) + S.One/120 - x**2/5040 + x**4/362880 + O(x**6)

    f = x**2*atan(x)
    assert fps(f, x, rational=False).truncate() == \
        x**3 - x**5/3 + O(x**6)


def test_fps__Add_expr():
    f = x*atan(x) - log(1 + x**2) / 2
    assert fps(f, x).truncate() == x**2/2 - x**4/12 + O(x**6)

    f = sin(x) + cos(x) - exp(x) + log(1 + x)
    assert fps(f, x).truncate() == x - 3*x**2/2 - x**4/4 + x**5/5 + O(x**6)

    f = 1/x + sin(x)
    assert fps(f, x).truncate() == 1/x + x - x**3/6 + x**5/120 + O(x**6)


@XFAIL
def test_xfail_fps__hyper():
    f = exp(x)
    assert fps(f, x, -oo).truncate() == O(1/x**6, (x, oo))

    f = sin(sqrt(x)) / x
    assert fps(f, x).truncate() == \
        (1/sqrt(x) - sqrt(x)/6 + x**Rational(3, 2)/120 - x**Rational(5, 2)/5040
         + x**Rational(7, 2)/362880 - x**Rational(9, 2)/39916800
         + x**Rational(11, 2)/6227020800 + O(x**6))

    f = x**n*sin(x**2)
    assert fps(f, x).truncate() == x**n*(x**2 + O(x**6))

    f = exp(x)*sin(x)/x
    assert fps(f, x).truncate() == 1 + x + x**2/3 - x**4/30 - x**5/90 + O(x**6)
