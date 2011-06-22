from sympy import Symbol, gamma, oo, nan, zoo, factorial, sqrt, Rational, log,\
        polygamma, EulerGamma, pi, uppergamma, S, expand_func, loggamma, sin, cos, \
        O, cancel, lowergamma, exp,  erf
from sympy.utilities.randtest import (test_derivative_numerically as td,
                                      random_complex_number as randcplx,
                                      test_numerically as tn)

x = Symbol('x')
y = Symbol('y')
n = Symbol('n', integer=True)

def test_gamma():
    assert gamma(nan) == nan
    assert gamma(oo) == oo

    assert gamma(-100) == zoo
    assert gamma(0) == zoo

    assert gamma(1) == 1
    assert gamma(2) == 1
    assert gamma(3) == 2

    assert gamma(102) == factorial(101)

    assert gamma(Rational(1,2)) == sqrt(pi)

    assert gamma(Rational(3, 2)) == Rational(1, 2)*sqrt(pi)
    assert gamma(Rational(5, 2)) == Rational(3, 4)*sqrt(pi)
    assert gamma(Rational(7, 2)) == Rational(15, 8)*sqrt(pi)

    assert gamma(Rational(-1, 2)) == -2*sqrt(pi)
    assert gamma(Rational(-3, 2)) == Rational(4, 3)*sqrt(pi)
    assert gamma(Rational(-5, 2)) == -Rational(8, 15)*sqrt(pi)

    assert gamma(Rational(-15, 2)) == Rational(256, 2027025)*sqrt(pi)

    assert gamma(x).diff(x) == gamma(x)*polygamma(0, x)

    assert gamma(x - 1).expand(func=True) == gamma(x)/(x-1)
    assert gamma(x + 2).expand(func=True, mul=False) == x*(x+1)*gamma(x)

    assert expand_func(gamma(x + Rational(3, 2))) == \
        (x + Rational(1, 2))*gamma(x + Rational(1, 2))

    assert expand_func(gamma(x - Rational(1, 2))) == \
        gamma(Rational(1, 2) + x)/(x - Rational(1, 2))

def test_gamma_series():
    assert gamma(x + 1).series(x, 0, 3) == \
        1 - x*EulerGamma + x**2*EulerGamma**2/2 + pi**2*x**2/12 + O(x**3)
    assert gamma(x).series(x, -1, 3) == \
        -1/x + EulerGamma - 1 + EulerGamma*x - EulerGamma**2*x/2 - pi**2*x/12 \
        - x + EulerGamma*x**2 + EulerGamma*pi**2*x**2/12 - \
        x**2*polygamma(2, 1)/6 + EulerGamma**3*x**2/6 - EulerGamma**2*x**2/2 \
        - pi**2*x**2/12 - x**2 + O(x**3)

def test_lowergamma():
    from sympy import meijerg
    assert lowergamma(x, y).diff(y) == y**(x-1)*exp(-y)
    assert td(lowergamma(randcplx(), y), y)
    assert lowergamma(x, y).diff(x) == \
           gamma(x)*polygamma(0, x) - uppergamma(x, y)*log(y) \
           + meijerg([], [1, 1], [0, 0, x], [], y)

    assert lowergamma(S.Half, x) == sqrt(pi)*erf(sqrt(x))
    assert not lowergamma(S.Half - 3, x).has(lowergamma)
    assert not lowergamma(S.Half + 3, x).has(lowergamma)
    assert lowergamma(S.Half, x, evaluate=False).has(lowergamma)
    assert tn(lowergamma(S.Half + 3, x, evaluate=False),
              lowergamma(S.Half + 3, x), x)
    assert tn(lowergamma(S.Half - 3, x, evaluate=False),
              lowergamma(S.Half - 3, x), x)

def test_uppergamma():
    from sympy import meijerg
    assert uppergamma(4, 0) == 6
    assert uppergamma(x, y).diff(y) == -y**(x-1)*exp(-y)
    assert td(uppergamma(randcplx(), y), y)
    assert uppergamma(x, y).diff(x) == \
           uppergamma(x, y)*log(y) + meijerg([], [1, 1], [0, 0, x], [], y)
    assert td(uppergamma(x, randcplx()), x)

    assert uppergamma(S.Half, x) == sqrt(pi)*(1 - erf(sqrt(x)))
    assert not uppergamma(S.Half - 3, x).has(uppergamma)
    assert not uppergamma(S.Half + 3, x).has(uppergamma)
    assert uppergamma(S.Half, x, evaluate=False).has(uppergamma)
    assert tn(uppergamma(S.Half + 3, x, evaluate=False),
              uppergamma(S.Half + 3, x), x)
    assert tn(uppergamma(S.Half - 3, x, evaluate=False),
              uppergamma(S.Half - 3, x), x)

def test_polygamma():

    assert polygamma(n, nan) == nan

    assert polygamma(0, oo) == oo
    assert polygamma(1, oo) == 0
    assert polygamma(5, oo) == 0

    assert polygamma(0, -9) == zoo

    assert polygamma(0, -9) == zoo
    assert polygamma(0, -1) == zoo

    assert polygamma(0, 0) == zoo

    assert polygamma(0, 1) == -EulerGamma
    assert polygamma(0, 7) == Rational(49, 20) - EulerGamma

    assert polygamma(1, 1) == pi**2/6
    assert polygamma(1, 2) == pi**2/6 - 1
    assert polygamma(1, 3) == pi**2/6 - Rational(5, 4)
    assert polygamma(3, 1) == pi**4 / 15
    assert polygamma(3, 5) == 6*(Rational(-22369,20736) + pi**4/90)
    assert polygamma(5, 1) == 8 * pi**6 / 63

    assert polygamma(3, 7*x).diff(x) == 7*polygamma(4, 7*x)

def test_polygamma_expand_func():
    assert polygamma(0, x).expand(func=True) == polygamma(0, x)
    assert polygamma(0, 2*x).expand(func=True) == \
           polygamma(0, x)/2 + polygamma(0, Rational(1, 2) + x)/2 + log(2)
    assert polygamma(1, 2*x).expand(func=True) == \
           polygamma(1, x)/4 + polygamma(1, Rational(1, 2) + x)/4
    assert polygamma(2, x).expand(func=True) == \
           polygamma(2, x)
    assert polygamma(0, -1 + x).expand(func=True) == \
           polygamma(0, x) - 1/(x - 1)
    assert polygamma(0, 1 + x).expand(func=True) == \
           1/x + polygamma(0, x )
    assert polygamma(0, 2 + x).expand(func=True) == \
           1/x + 1/(1 + x) + polygamma(0, x)
    assert polygamma(0, 3 + x).expand(func=True) == \
           polygamma(0, x) + 1/x + 1/(1 + x) + 1/(2 + x)
    assert polygamma(0, 4 + x).expand(func=True) == \
           polygamma(0, x) + 1/x + 1/(1 + x) + 1/(2 + x) + 1/(3 + x)
    assert polygamma(1, 1 + x).expand(func=True) == \
           polygamma(1, x) - 1/x**2
    assert polygamma(1, 2 + x).expand(func=True, multinomial=False) == \
           polygamma(1, x) - 1/x**2 - 1/(1 + x)**2
    assert polygamma(1, 3 + x).expand(func=True, multinomial=False) == \
           polygamma(1, x) - 1/x**2 - 1/(1 + x)**2 - 1/(2 + x)**2
    assert polygamma(1, 4 + x).expand(func=True, multinomial=False) == \
           polygamma(1, x) - 1/x**2 - 1/(1 + x)**2 - \
           1/(2 + x)**2 - 1/(3 + x)**2
    assert polygamma(0, x + y).expand(func=True) == \
           polygamma(0, x + y)
    assert polygamma(1, x + y).expand(func=True) == \
           polygamma(1, x + y)
    assert polygamma(1, 3 + 4*x + y).expand(func=True, multinomial=False) == \
           polygamma(1, y + 4*x) - 1/(y + 4*x)**2 - \
           1/(1 + y + 4*x)**2 - 1/(2 + y + 4*x)**2
    assert polygamma(3, 3 + 4*x + y).expand(func=True, multinomial=False) == \
           polygamma(3, y + 4*x) - 6/(y + 4*x)**4 - \
           6/(1 + y + 4*x)**4 - 6/(2 + y + 4*x)**4
    assert polygamma(3, 4*x + y + 1).expand(func=True, multinomial=False) == \
           polygamma(3, y + 4*x) - 6/(y + 4*x)**4
    e = polygamma(3, 4*x + y + S(3)/2)
    assert e.expand(func=True) == e
    e = polygamma(3, x + y + S(3)/4)
    assert e.expand(func = True, basic = False) == e

def test_loggamma():
    s1 = loggamma(1/(x+sin(x))+cos(x)).nseries(x,n=4)
    s2 = (-log(2*x)-1)/(2*x) - log(x/pi)/2 + (4-log(2*x))*x/24 + O(x**2)
    assert (s1 - s2).expand(force=True).removeO() == 0
    s1 = loggamma(1/x).series(x)
    s2 = (1/x-S(1)/2)*log(1/x) - 1/x  + log(2*pi)/2 + \
         x/12 - x**3/360 + x**5/1260 +  O(x**7)
    assert ((s1 - s2).expand(force=True)).removeO() == 0

    def tN(N, M):
        assert loggamma(1/x)._eval_nseries(x,n=N,logx=None).getn() == M
    tN(0, 0)
    tN(1, 1)
    tN(2, 3)
    tN(3, 3)
    tN(4, 5)
    tN(5, 5)

def test_polygamma_expansion():
    # A. & S., pa. 259 and 260
    assert polygamma(0, 1/x).nseries(x, n=3) \
           == -log(x) - x/2 - x**2/12 + O(x**4)
    assert polygamma(1, 1/x).series(x, n=5) \
           == x + x**2/2 + x**3/6 + O(x**5)
    assert polygamma(3, 1/x).nseries(x, n=8) \
           == 2*x**3 + 3*x**4 + 2*x**5 - x**7 + 4*x**9/3 + O(x**11)
