from sympy import sin, cos, exp, asin, atan, asinh, sqrt, O, Derivative, Symbol
from sympy.functions import factorial, RisingFactorial, airyai
from sympy.core import Function
from sympy.series.formal import simpleDE, DEtoRE, FormalSeries, solveRE
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
    assert DEtoRE(DE, r, k) == (k + 1)*(k + 2)*r(k + 2) + r(k)
    DE = 2*f(x) - 2*Derivative(f(x), x) + Derivative(f(x), x, x)
    assert DEtoRE(DE, r, k) == (-2*k - 2)*r(k + 1) + (k + 1)*(k + 2)*r(k + 2) + 2*r(k)
    DE = (x**10 + 4)*Derivative(f(x), x) + x*(x**10 - 1)*Derivative(f(x), x, x)
    assert DEtoRE(DE, r, k) == k*(k - 1)*r(k) + k*r(k) - (k + 9)*(k + 10)*r(k + 10) + (4*k + 40)*r(k + 10)
    DE = -x*f(x) + Derivative(f(x), x, x)
    assert DEtoRE(DE, r, k) == (k + 2)*(k + 3)*r(k + 3) - r(k)
    DE = 2*n*f(x) + (x**2 - 1)*Derivative(f(x), x)
    assert DEtoRE(DE, r, k) == k*r(k) + 2*n*r(k + 1) + (-k - 2)*r(k + 2)

def test_solveRE():
    r = Function('r')

    RE = (k + 1)*r(k + 1) - r(k)
    assert solveRE(RE, r, k, exp(x)) == [(1/factorial(k), k)]
    RE = (k + 1)*(k + 2)*r(k + 2) + r(k)
    assert solveRE(RE, r, k, sin(x)) == [((-1/4)**k/(RisingFactorial(3/2, k)*factorial(k)), 2*k + 1)]
    assert solveRE(RE, r, k, cos(x)) == [((-1/4)**k/(RisingFactorial(1/2, k)*factorial(k)), 2*k)]
    RE = k*(k - 1)*r(k) + k*r(k) - (k + 1)*(k + 2)*r(k + 2)
    assert solveRE(RE, r, k, asin(x)) == [(RisingFactorial(1/2, k)**2/(RisingFactorial(3/2, k)*factorial(k)), 2*k + 1)]
    RE = k*(k - 1)*r(k) + 2*k*r(k) + (k + 1)*(k + 2)*r(k + 2)
    assert solveRE(RE, r, k, atan(x)) == [((-1)**k*RisingFactorial(1/2, k)/RisingFactorial(3/2, k), 2*k + 1)]
    RE = k*(k - 1)*r(k) + k*r(k) + (k + 1)*(k + 2)*r(k + 2) - r(k)
    assert solveRE(RE, r, k, exp(asinh(x))) == [((-1)**k*RisingFactorial(-1/2, k)/factorial(k), 2*k),
            ((-1)**k*RisingFactorial(0, k)/RisingFactorial(3/2, k), 2*k + 1)]
    RE = 4*k*(k + 1)*r(k + 1) + 10*(k + 1)*r(k + 1) + r(k) + 2*r(k + 1)
    assert solveRE(RE, r, k, sin(sqrt(x))/x) == [((-1/4)**k/(RisingFactorial(3/2, k)*factorial(k)), k - 1/2)]


def test_series():
    s = FormalSeries(x, function=sin(x))
    assert s.as_series() == x - x**3/6 + x**5/120 + O(x**6)
    s = FormalSeries(x, function=exp(x))
    assert s.as_series() == 1 + x + x**2/2 + x**3/6 + x**4/24 + x**5/120 + O(x**6)
    s = FormalSeries(x, function=x + 1/(1-x))
    assert s.as_series() == 1 + 2*x + x**2 + x**3 + x**4 + x**5 + O(x**6)
    s = FormalSeries(x, function=sin(x) + cos(x))
    assert s.as_series() == 1 + x - x**2/2 - x**3/6 + x**4/24 + x**5/120 + O(x**6)
    s = FormalSeries(x, sequence=(1, 2, 3))
    assert s.as_series() == 1 + 2*x + 3*x**2 + x**3 + 2*x**4 + 3*x**5 + O(x**6)


def test_series_add():
    s1 = FormalSeries(x, function=sin(x))
    s2 = FormalSeries(x, function=exp(x))
    s = s1 + s2
    assert s.as_series() == 1 + 2*x + x**2/2 + x**4/24 + x**5/60 + O(x**6)
    s1 = FormalSeries(x, function=1/(1-x))
    s2 = FormalSeries(x, sequence=(1, 2, 3))
    s = s1 + s2
    assert s.as_series() == 2 + 3*x + 4*x**2 + 2*x**3 + 3*x**4 + 4*x**5 + O(x**6)


def test_series_mul():
    s1 = FormalSeries(x, function=sin(x))
    s2 = FormalSeries(x, function=exp(x))
    s = s1 * s2
    assert s.as_series() == x + x**2 + x**3/3 - x**5/30 + O(x**6)
    s1 = FormalSeries(x, function=1/(1-x))
    s2 = FormalSeries(x, function=1/(1+x))
    s = s1 * s2
    assert s.as_series() == 1 + x**2 + x**4 + O(x**6)
    s = FormalSeries(x, sequence=(1, 2, 3))
    assert (s*(-1)).as_series() == -1 - 2*x - 3*x**2 - x**3 - 2*x**4 - 3*x**5 + O(x**6)
    assert (s*x**2).as_series() == x**2 + 2*x**3 + 3*x**4 + x**5 + O(x**6)


def test_series_div():
    s1 = FormalSeries(x, function=sin(x))
    s2 = FormalSeries(x, function=exp(x))
    s = s1 / s2
    assert s.as_series() == x - x**2 + x**3/3 - x**5/30 + O(x**6)
    s1 = FormalSeries(x, function=1/(1-x))
    s2 = FormalSeries(x, function=1/(1+x))
    s = s1 / s2
    assert s.as_series() == 1 + 2*x + 2*x**2 + 2*x**3 + 2*x**4 + 2*x**5 + O(x**6)


def test_series_inverse():
    s = FormalSeries(x, function=1/(1-x))
    assert s.inverse().as_series() == 1 - x + O(x**6)
    s = FormalSeries(x, function=exp(x))
    assert s.inverse().as_series() == 1 - x + x**2/2 - x**3/6 + x**4/24 - x**5/120 + O(x**6)
    s = FormalSeries(x, sequence=(1, 2, 3))
    assert s.inverse().as_series() == 1 - 2*x + x**2 + 3*x**3 - 9*x**4 + 9*x**5 + O(x**6)


def test_series_compose():
    s1 = FormalSeries(x, function=sin(x))
    s2 = FormalSeries(x, function=exp(x))
    s = s2.compose(s1)
    assert s.as_series() == 1 + x + x**2/2 - x**4/8 - x**5/15 + O(x**6)
