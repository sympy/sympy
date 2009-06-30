from sympy import (Symbol, Sum, oo, Real, Rational, sum, pi, cos, zeta,
Catalan, exp, log, factorial, sqrt, E, sympify, binomial, EulerGamma, Function)
from sympy.concrete.summations import getab
from sympy.utilities.pytest import XFAIL

a, b, c, d, m, k = map(Symbol, 'abcdmk')
n = Symbol('n', integer=True)

def test_arithmetic_sums():
    assert sum(1, (n, a, b)) == b-a+1
    assert sum(1, (n, 1, 10)) == 10
    assert sum(2*n, (n, 0, 10**10)) == 100000000010000000000
    assert sum(4*n*m, (n, a, 1), (m, 1, d)).expand() == \
        2*d + 2*d**2 + a*d + a*d**2 - d*a**2 - a**2*d**2
    assert sum(cos(n), (n, -2, 1)) == cos(-2)+cos(-1)+cos(0)+cos(1)

def test_polynomial_sums():
    assert sum(n**2, (n, 3, 8)) == 199
    assert sum(n, (n, a, b)) == \
        ((a+b)*(b-a+1)/2).expand()
    assert sum(n**2, (n, 1, b)) == \
        ((2*b**3+3*b**2+b)/6).expand()
    assert sum(n**3, (n, 1, b)) == \
        ((b**4+2*b**3+b**2)/4).expand()
    assert sum(n**6, (n, 1, b)) == \
        ((6*b**7+21*b**6+21*b**5-7*b**3+b)/42).expand()

def test_geometric_sums():
    assert sum(pi**n, (n, 0, b)) == (1-pi**(b+1)) / (1-pi)
    assert sum(2 * 3**n, (n, 0, b)) == 3**(b+1) - 1
    assert sum(Rational(1,2)**n, (n, 1, oo)) == 1
    assert sum(2**n, (n, 0, b)) == 2**(b+1) - 1
    assert sum(2**n, (n, 1, oo)) == oo
    assert sum(2**(-n), (n, 1, oo)) == 1
    assert sum(3**(-n), (n, 4, oo)) == Rational(1,54)
    assert sum(2**(-4*n+3), (n, 1, oo)) == Rational(8,15)
    assert sum(2**(n+1), (n, 1, b)).expand() == 4*(2**b-1)

def test_composite_sums():
    f = Rational(1,2)*(7 - 6*n + Rational(1,7)*n**3)
    s = sum(f, (n, a, b))
    assert not isinstance(s, Sum)
    A = 0
    for i in range(-3, 5):
        A += f.subs(n, i)
    B = s.subs(a,-3).subs(b,4)
    assert A == B

fac = factorial

def NS(e, n=15, **options):
    return str(sympify(e).evalf(n, **options))

def test_evalf_fast_series():
    # Euler transformed series for sqrt(1+x)
    assert NS(Sum(fac(2*n+1)/fac(n)**2/2**(3*n+1), (n, 0, oo)), 100) == NS(sqrt(2), 100)

    # Some series for exp(1)
    estr = NS(E, 100)
    assert NS(Sum(1/fac(n), (n, 0, oo)), 100) == estr
    assert NS(1/Sum((1-2*n)/fac(2*n), (n, 0, oo)), 100) == estr
    assert NS(Sum((2*n+1)/fac(2*n), (n, 0, oo)), 100) == estr
    assert NS(Sum((4*n+3)/2**(2*n+1)/fac(2*n+1), (n, 0, oo))**2, 100) == estr

    pistr = NS(pi, 100)
    # Ramanujan series for pi
    assert NS(9801/sqrt(8)/Sum(fac(4*n)*(1103+26390*n)/fac(n)**4/396**(4*n), (n, 0, oo)), 100) == pistr
    assert NS(1/Sum(binomial(2*n,n)**3 * (42*n+5)/2**(12*n+4), (n, 0, oo)), 100) == pistr
    # Machin's formula for pi
    assert NS(16*Sum((-1)**n/(2*n+1)/5**(2*n+1), (n, 0, oo)) - \
        4*Sum((-1)**n/(2*n+1)/239**(2*n+1), (n, 0, oo)), 100) == pistr

    # Apery's constant
    astr = NS(zeta(3), 100)
    P = 126392*n**5 + 412708*n**4 + 531578*n**3 + 336367*n**2 + 104000*n + 12463
    assert NS(Sum((-1)**n * P / 24 * (fac(2*n+1)*fac(2*n)*fac(n))**3 / fac(3*n+2) / fac(4*n+3)**3, (n, 0, oo)), 100) == astr
    assert NS(Sum((-1)**n * (205*n**2 + 250*n + 77)/64 * fac(n)**10 / fac(2*n+1)**5, (n, 0, oo)), 100) == astr

def test_evalf_fast_series_issue998():
    # Catalan's constant
    assert NS(Sum((-1)**(n-1)*2**(8*n)*(40*n**2-24*n+3)*fac(2*n)**3*\
        fac(n)**2/n**3/(2*n-1)/fac(4*n)**2, (n, 1, oo))/64, 100) == \
        NS(Catalan, 100)
    astr = NS(zeta(3), 100)
    assert NS(5*Sum((-1)**(n-1)*fac(n)**2 / n**3 / fac(2*n), (n, 1, oo))/2, 100) == astr
    assert NS(Sum((-1)**(n-1)*(56*n**2-32*n+5) / (2*n-1)**2 * fac(n-1)**3 / fac(3*n), (n, 1, oo))/4, 100) == astr

def test_evalf_slow_series():
    assert NS(Sum((-1)**n / n, (n, 1, oo)), 15) == NS(-log(2), 15)
    assert NS(Sum((-1)**n / n, (n, 1, oo)), 50) == NS(-log(2), 50)
    assert NS(Sum(1/n**2, (n, 1, oo)), 15) == NS(pi**2/6, 15)
    assert NS(Sum(1/n**2, (n, 1, oo)), 100) == NS(pi**2/6, 100)
    assert NS(Sum(1/n**2, (n, 1, oo)), 500) == NS(pi**2/6, 500)
    assert NS(Sum((-1)**n / (2*n+1)**3, (n, 0, oo)), 15) == NS(pi**3/32, 15)
    assert NS(Sum((-1)**n / (2*n+1)**3, (n, 0, oo)), 50) == NS(pi**3/32, 50)

def test_euler_maclaurin():
    # Exact polynomial sums with E-M
    def check_exact(f, a, b, m, n):
        A = Sum(f, (k, a, b))
        s, e = A.euler_maclaurin(m, n)
        assert (e == 0) and (s.expand() == A.doit())
    check_exact(k**4, a, b, 0, 2)
    check_exact(k**4 + 2*k, a, b, 1, 2)
    check_exact(k**4 + k**2, a, b, 1, 5)
    check_exact(k**5, 2, 6, 1, 2)
    check_exact(k**5, 2, 6, 1, 3)
    # Not exact
    assert Sum(k**6, (k, a, b)).euler_maclaurin(0, 2)[1] != 0
    # Numerical test
    for m, n in [(2, 4), (2, 20), (10, 20), (18, 20)]:
        A = Sum(1/k**3, (k, 1, oo))
        s, e = A.euler_maclaurin(m, n)
        assert abs((s-zeta(3)).evalf()) < e.evalf()

def test_evalf_euler_maclaurin():
    assert NS(Sum(1/k**k, (k, 1, oo)), 15) == '1.29128599706266'
    assert NS(Sum(1/k**k, (k, 1, oo)), 50) == '1.2912859970626635404072825905956005414986193682745'
    assert NS(Sum(1/k-log(1+1/k), (k, 1, oo)), 15) == NS(EulerGamma, 15)
    assert NS(Sum(1/k-log(1+1/k), (k, 1, oo)), 50) == NS(EulerGamma, 50)
    assert NS(Sum(log(k)/k**2, (k, 1, oo)), 15) == '0.937548254315844'
    assert NS(Sum(log(k)/k**2, (k, 1, oo)), 50) == '0.93754825431584375370257409456786497789786028861483'
    assert NS(Sum(1/k, (k, 1000000, 2000000)), 15) == '0.693147930560008'
    assert NS(Sum(1/k, (k, 1000000, 2000000)), 50) == '0.69314793056000780941723211364567656807940638436025'

@XFAIL
def test_simple_products():
    assert Product(2, (n, a, b)) == 2**(b-a+1)
    assert Product(n, (n, 1, b)) == factorial(b)
    assert Product(n**3, (n, 1, b)) == factorial(b)**3
    assert Product(3**(2+n), (n, a, b)) \
           == 3**(2*(1-a+b)+b/2+(b**2)/2+a/2-(a**2)/2)
    assert Product(cos(n), (n, 3, 5)) == cos(3)*cos(4)*cos(5)
    # If Product managed to evaluate this one, it most likely got it wrong!
    assert isinstance(Product(n**n, (n, 1, b)), Product)

@XFAIL
def test_rational_products():
    assert Product(1+1/n, (n, a, b)) == (1+b)/a
    assert Product(n+1, (n, a, b)) == factorial(1+b)/factorial(a)
    assert Product((n+1)/(n-1), (n, a, b)) == b*(1+b)/(a*(a-1))
    assert Product(n/(n+1)/(n+2), (n, a, b)) \
           == a*factorial(a+1)/(b+1)/factorial(b+2)
    assert Product(n*(n+1)/(n-1)/(n-2), (n, a, b)) \
           == b**2*(b-1)*(1+b)/(a-1)**2/(a*(a-2))

@XFAIL
def test_wallis_product():
    # Wallis product, given in two different forms to ensure that Product
    # can factor simple rational expressions
    A = Product(4*n**2 / (4*n**2-1), (n, 1, b))
    B = Product((2*n)*(2*n)/(2*n-1)/(2*n+1), (n, 1, b))
    half = Rational(1,2)
    R = pi/2 * factorial(b)**2 / factorial(b-half) / factorial(b+half)
    assert A == R
    assert B == R
    # This one should eventually also be doable (Euler's product formula for sin)
    # assert Product(1+x/n**2, (n, 1, b)) == ...

def test_telescopic_sums():
    #checks also input 2 of comment 1 issue 1028
    assert Sum(1/k - 1/(k+1),(k,1,n)).doit() == 1 - 1/(1 + n)
    f = Function("f")
    assert Sum(f(k)-f(k+2),(k,m,n)).doit() == -f(1+n) - f(2+n) + f(m) + f(1+m)
    assert Sum(cos(k)-cos(k+3),(k,1,n)).doit() == -cos(1 + n) - cos(2 + n) - \
                                           cos(3 + n) + cos(1) + cos(2) + cos(3)

def test_Sum_limit_subs():
    assert Sum(a*exp(a), (a, -2, 2)) == Sum(a*exp(a), (a, -b, b)).subs(b,2)
