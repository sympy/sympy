from sympy import sin, exp, asin, Derivative, Symbol
from sympy.core import Function
from sympy.functions import airyai
from sympy.series.formal import simpleDE, DEtoRE, FormalSeries
from sympy.abc import x, k, m, n


def test_simpleDE():
    f = Function('f')

    assert simpleDE(x, sin(x), f) == f(x) + Derivative(f(x), x, x)
    assert simpleDE(x, exp(x), f) == -f(x) + Derivative(f(x), x)
    assert simpleDE(x, sin(x)*exp(x), f) == 2*f(x) - 2*Derivative(f(x), x) + Derivative(f(x), x, x)
    assert simpleDE(x, airyai(x), f) == -x*f(x) + Derivative(f(x), x, x)
    assert simpleDE(x, asin(x), f) == x*Derivative(f(x), x) + (x**2 - 1)*Derivative(f(x), x, x)
    assert simpleDE(x, asin(x**5), f) == (x**10 + 4)*Derivative(f(x), x) + x*(x**10 - 1)*Derivative(f(x), x, x)
    assert simpleDE(x, asin(x)**3, f) == x**4*Derivative(f(x), x, x, x, x) + 6*x**3*Derivative(f(x), x, x, x) + \
            7*x**2*Derivative(f(x), x, x) - 2*x**2*Derivative(f(x), x, x, x, x) + x*Derivative(f(x), x) - \
            6*x*Derivative(f(x), x, x, x) - 4*Derivative(f(x), x, x) + Derivative(f(x), x, x, x, x)
    assert simpleDE(x, x**n*exp(m*x), f) == x*Derivative(f(x), x) - (m*x + n)*f(x)
    assert simpleDE(x, ((1+x)/(1-x))**n, f) == 2*n*f(x) + (x**2 - 1)*Derivative(f(x), x)
    assert simpleDE(x, exp(m*x)*sin(n*x), f) == -2*m*Derivative(f(x), x) + (m**2 + n**2)*f(x) + Derivative(f(x), x, x)


def test_DEtoRE():
    f = Function('f')
    r = Function('r')

    DE = f(x) + Derivative(f(x), x, x)
    assert DEtoRE(x, DE, r, k) == (k + 1)*(k + 2)*r(k + 2) + r(k)
    DE = 2*f(x) - 2*Derivative(f(x), x) + Derivative(f(x), x, x)
    assert DEtoRE(x, DE, r, k) == (-2*k - 2)*r(k + 1) + (k + 1)*(k + 2)*r(k + 2) + 2*r(k)
    DE = 2*n*f(x) + (x**2 - 1)*Derivative(f(x), x)
    assert DEtoRE(x, DE, r, k) == 2*n*r(k + 1) + (-k - 2)*r(k + 2) + k*r(k)
    DE = (x**10 + 4)*Derivative(f(x), x) + x*(x**10 - 1)*Derivative(f(x), x, x)
    assert DEtoRE(x, DE, r, k) == k*(k - 1)*r(k) + k*r(k) - (k + 9)*(k + 10)*r(k + 10) + (4*k + 40)*r(k + 10)
    DE = -x*f(x) + Derivative(f(x), x, x)
    assert DEtoRE(x, DE, r, k) == (k + 2)*(k + 3)*r(k + 3) - r(k)


def test_solveRE():
    pass


def test_series():
    s = FormalSeries(x, function=sin(x))
    assert s.as_series() == x**5/120 - x**3/6 + x
    s = FormalSeries(x, function=exp(x))
    assert s.as_series() == x**6/720 + x**5/120 + x**4/24 + x**3/6 + x**2/2 + x + 1
    s = FormalSeries(x, function=1/(1-x))
    assert s.as_series() == x**6 + x**5 + x**4 + x**3 + x**2 + x + 1
    s = FormalSeries(x, function=1/(1+x))
    assert s.as_series() == x**6 - x**5 + x**4 - x**3 + x**2 - x + 1
    s = FormalSeries(x, sequence=(1, 2, 3))
    assert s.as_series() == x**6 + 3*x**5 + 2*x**4 + 1*x**3 + 3*x**2 + 2*x + 1


def test_series_add():
    s1 = FormalSeries(x, function=sin(x))
    s2 = FormalSeries(x, function=exp(x))
    s = s1 + s2
    assert s.as_series() == x**6/720 + x**5/60 + x**4/24 + x**2/2 + 2*x + 1
    s1 = FormalSeries(x, function=1/(1-x))
    s2 = FormalSeries(x, sequence=(1, 2, 3))
    s = s1 + s2
    assert s.as_series() == 2*x**6 + 4*x**5 + 3*x**4 + 2*x**3 + 4*x**2 + 3*x + 2


def test_series_mul():
    s1 = FormalSeries(x, function=sin(x))
    s2 = FormalSeries(x, function=exp(x))
    s = s1 * s2
    assert s.as_series() == -x**6/90 + -x**5/30 + x**3/3 + x**2 + x
    s1 = FormalSeries(x, function=1/(1-x))
    s2 = FormalSeries(x, function=1/(1+x))
    s = s1 * s2
    assert s.as_series() == x**6 + x**4 + x**2 + 1


def test_series_scale_shift():
    s = FormalSeries(x, sequence=(1, 2, 3))
    assert s.scale(-1).as_series() == -x**6 - 3*x**5 - 2*x**4 - 1*x**3 - 3*x**2 - 2*x - 1
    assert s.shift(2).as_series() == x**8 + 3*x**7 + 2*x**6 + 1*x**5 + 3*x**4 + 2*x**3 + x**2
    s = FormalSeries(x, function=sin(x))
    assert s.scale(0).as_series() == 0
    assert s.shift(1).scale(1).as_series() == x**6/120 -x**4/6 + x**2


def test_series_div():
    s1 = FormalSeries(x, function=sin(x))
    s2 = FormalSeries(x, function=exp(x))
    s = s1 / s2
    assert s.as_series() == x**6/90 - x**5/30 + x**3/3 - x**2 + x
    s1 = FormalSeries(x, function=1/(1-x))
    s2 = FormalSeries(x, function=1/(1+x))
    s = s1 / s2
    assert s.as_series() == 2*x**6 + 2*x**5 + 2*x**4 + 2*x**3 + 2*x**2 + 2*x + 1


def test_series_inverse():
    s = FormalSeries(x, function=1/(1-x))
    assert s.invert().as_series() == -x + 1
    s = FormalSeries(x, function=exp(x))
    assert s.invert().as_series() == x**6/720 - x**5/120 + x**4/24 - x**3/6 + x**2/2 - x + 1
    s = FormalSeries(x, sequence=(1, 2, 3))
    assert s.invert().as_series() == 9*x**6 + 9*x**5 - 9*x**4 + 3*x**3 + x**2 - 2*x + 1


def test_series_compose():
    s1 = FormalSeries(x, function=sin(x))
    s2 = FormalSeries(x, function=exp(x))
    s = s2.compose(s1)
    assert s.as_series() == -x**6/240 - x**5/15 - x**4/8 + x**2/2 + x + 1
