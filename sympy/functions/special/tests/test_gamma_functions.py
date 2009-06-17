from sympy import Symbol, gamma, oo, nan, zoo, factorial, sqrt, Rational, log,\
        polygamma, EulerGamma, pi, uppergamma, S, expand_func

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

    assert expand_func(gamma(x + Rational(3, 2))) ==\
    (x + Rational(1, 2))*gamma(x + Rational(1, 2))

def test_lowergamma():
    pass

def test_uppergamma():
    assert uppergamma(4, 0) == 6

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
    assert polygamma(0, 2*x).expand(func=True) == log(2) + polygamma(0, 2*x)
    assert polygamma(2, x).expand(func=True, basic=False) == polygamma(2, x)
    #assert polygamma(2, 3*x).expand(func=True) == polygamma(2, 3*x)/9
    #assert polygamma(3, 4*x).expand(func=True,basic=False) == polygamma(3, 4*x)/64
    assert polygamma(0, 1 + x).expand(func=True, basic=False) == 1 + x + polygamma(0, x )
    assert polygamma(0, 2 + x).expand(func=True, basic=False) == 5 + 2*x + polygamma(0, x)
    assert polygamma(0, 3 + x).expand(func=True, basic=False) == 12 + 3*x + polygamma(0, x)
    assert polygamma(0, 4 + x).expand(func=True, basic=False) == 22 + 4*x + polygamma(0, x)
    assert polygamma(1, 1 + x).expand(func=True,basic=False) == -1 - x + polygamma(1, x)
    assert polygamma(1, 2 + x).expand(func=True,basic=False) == -5 - 2*x + polygamma(1, x)
    assert polygamma(1, 3 + x).expand(func=True,basic=False) == -12 - 3*x + polygamma(1, x)
    assert polygamma(1, 4 + x).expand(func=True,basic=False) == -22 - 4*x + polygamma(1, x)
    assert polygamma(0, x + y).expand(func= True, basic = False) == polygamma(0, x + y)
    assert polygamma(1, x + y).expand(func = True, basic = False) == polygamma(1, x + y)
    assert polygamma(1, 3 + 4*x + y).expand(func = True,basic = False) == -12 - 12*x - 3*y + polygamma(1, y + 4*x)
    assert polygamma(3, 3 + 4*x + y).expand(func = True,basic = False) == -72 - 72*x -  18*y + polygamma(3, y + 4*x)
    assert polygamma(3, 4*x + y + 1).expand(func = True, basic = False) == -6 - 24*x - 6*y + polygamma(3, y + 4*x)
    assert polygamma(3, 4*x + y + S(3)/2).expand(func = True, basic = False) == polygamma(3, S(3)/2 + y + 4*x)
    assert polygamma(3, x + y + S(3)/4).expand(func = True, basic = False) == polygamma(3, S(3)/4 + x + y)

def test_loggamma():
    pass
