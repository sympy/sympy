from sympy import (binomial, Catalan, cos, Derivative, E, exp, EulerGamma,
                   factorial, Function, harmonic, Integral, log, nan, oo, pi,
                   Product, product, Rational, S, sqrt, Sum, summation, Symbol,
                   sympify, zeta, oo)
from sympy.abc import a, b, c, d, k, m, n, x, y, z
from sympy.concrete.summations import telescopic
from sympy.utilities.pytest import XFAIL, raises

n = Symbol('n', integer=True)

def test_arithmetic_sums():
    assert summation(1, (n, a, b)) == b-a+1
    assert Sum(S.NaN, (n, a, b)) is S.NaN
    assert Sum(x, (n, a, a)).doit() == x
    assert Sum(x, (x, a, a)).doit() == a
    assert Sum(x, (n, 1, a)).doit() == a*x
    lo, hi = 1, 2
    s1 = Sum(n, (n, lo, hi))
    s2 = Sum(n, (n, hi, lo))
    assert s1 != s2
    assert s1.doit() == s2.doit() == 3
    lo, hi = x, x + 1
    s1 = Sum(n, (n, lo, hi))
    s2 = Sum(n, (n, hi, lo))
    assert s1 != s2
    assert s1.doit() == s2.doit() == 2*x + 1
    assert Sum(Integral(x, (x, 1, y)) + x, (x, 1, 2)).doit() == \
        y**2 + 2
    assert summation(1, (n, 1, 10)) == 10
    assert summation(2*n, (n, 0, 10**10)) == 100000000010000000000
    assert summation(4*n*m, (n, a, 1), (m, 1, d)).expand() == \
        2*d + 2*d**2 + a*d + a*d**2 - d*a**2 - a**2*d**2
    assert summation(cos(n), (n, -2, 1)) == cos(-2) + cos(-1) + cos(0) + cos(1)
    assert summation(cos(n), (n, x, x + 2)) == cos(x) + cos(x + 1) + cos(x + 2)
    assert isinstance(summation(cos(n), (n, x, x + S.Half)), Sum)

def test_polynomial_sums():
    assert summation(n**2, (n, 3, 8)) == 199
    assert summation(n, (n, a, b)) == \
        ((a+b)*(b-a+1)/2).expand()
    assert summation(n**2, (n, 1, b)) == \
        ((2*b**3+3*b**2+b)/6).expand()
    assert summation(n**3, (n, 1, b)) == \
        ((b**4+2*b**3+b**2)/4).expand()
    assert summation(n**6, (n, 1, b)) == \
        ((6*b**7+21*b**6+21*b**5-7*b**3+b)/42).expand()

def test_geometric_sums():
    assert summation(pi**n, (n, 0, b)) == (1-pi**(b+1)) / (1-pi)
    assert summation(2 * 3**n, (n, 0, b)) == 3**(b+1) - 1
    assert summation(Rational(1,2)**n, (n, 1, oo)) == 1
    assert summation(2**n, (n, 0, b)) == 2**(b+1) - 1
    assert summation(2**n, (n, 1, oo)) == oo
    assert summation(2**(-n), (n, 1, oo)) == 1
    assert summation(3**(-n), (n, 4, oo)) == Rational(1,54)
    assert summation(2**(-4*n+3), (n, 1, oo)) == Rational(8,15)
    assert summation(2**(n+1), (n, 1, b)).expand() == 4*(2**b-1)

def test_harmonic_sums():
    assert summation(1/k, (k, 0, n)) == Sum(1/k, (k, 0, n))
    assert summation(1/k, (k, 1, n)) == harmonic(n)
    assert summation(n/k, (k, 1, n)) == n*harmonic(n)
    assert summation(1/k, (k, 5, n)) == harmonic(n) - harmonic(4)

def test_composite_sums():
    f = Rational(1,2)*(7 - 6*n + Rational(1,7)*n**3)
    s = summation(f, (n, a, b))
    assert not isinstance(s, Sum)
    A = 0
    for i in range(-3, 5):
        A += f.subs(n, i)
    B = s.subs(a,-3).subs(b,4)
    assert A == B

def test_hypergeometric_sums():
    assert summation(binomial(2*k, k)/4**k, (k, 0, n)) == (1 + 2*n)*binomial(2*n, n)/4**n

def test_other_sums():
    f = m**2 + m*exp(m)
    g = 3*exp(S(3)/2)/2 + exp(S(1)/2)/2 - exp(-S(1)/2)/2 - 3*exp(-S(3)/2)/2 + 5

    assert summation(f, (m, -S(3)/2, S(3)/2)).expand() == g
    assert summation(f, (m, -1.5, 1.5)).evalf().epsilon_eq(g.evalf(), 1e-10)

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

def test_simple_products():
    assert Product(S.NaN, (x, 1, 3)) is S.NaN
    assert product(S.NaN, (x, 1, 3)) is S.NaN
    assert Product(x, (n, a, a)).doit() == x
    assert Product(x, (x, a, a)).doit() == a
    assert Product(x, (y, 1, a)).doit() == x**a
    lo, hi = 1, 2
    s1 = Product(n, (n, lo, hi))
    s2 = Product(n, (n, hi, lo))
    assert s1 != s2
    assert s1.doit() == s2.doit() == 2
    lo, hi = x, x + 1
    s1 = Product(n, (n, lo, hi))
    s2 = Product(n, (n, hi, lo))
    assert s1 != s2
    assert s1.doit() == s2.doit() == x*(x + 1)
    assert Product(Integral(2*x, (x,1,y)) + 2*x, (x,1,2)).doit() == \
        (y**2 + 1)*(y**2 + 3)
    assert product(2, (n, a, b)) == 2**(b-a+1)
    assert product(n, (n, 1, b)) == factorial(b)
    assert product(n**3, (n, 1, b)) == factorial(b)**3
    assert product(3**(2+n), (n, a, b)) \
           == 3**(2*(1-a+b)+b/2+(b**2)/2+a/2-(a**2)/2)
    assert product(cos(n), (n, 3, 5)) == cos(3)*cos(4)*cos(5)
    assert product(cos(n), (n, x, x + 2)) == cos(x)*cos(x + 1)*cos(x + 2)
    assert isinstance(product(cos(n), (n, x, x + S.Half)), Product)
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

    # dummy variable shouldn't matter
    assert telescopic(1/m, -m/(1+m),(m, n-1, n)) == \
           telescopic(1/k, -k/(1+k),(k, n-1, n))

    assert Sum(1/x/(x - 1), (x, a, b)).doit() == -((a - b - 1)/(b*(a - 1)))

def test_sum_reconstruct():
    s = Sum(n**2, (n, -1, 1))
    assert s == Sum(*s.args)
    raises(ValueError, lambda: Sum(x, x))
    raises(ValueError, lambda: Sum(x, (x, 1)))

def test_Sum_limit_subs():
    assert Sum(a*exp(a), (a, -2, 2)) == Sum(a*exp(a), (a, -b, b)).subs(b, 2)
    assert Sum(a, (a, Sum(b, (b, 1, 2)), 4)).subs(Sum(b, (b, 1, 2)), c) == \
           Sum(a, (a, c, 4))

@XFAIL
def test_issue2166():
    assert Sum(x, (x, 1, x)).subs(x, a) == Sum(x, (x, 1, a))

def test_Sum_doit():
    assert Sum(n*Integral(a**2), (n, 0, 2)).doit() == a**3
    assert Sum(n*Integral(a**2), (n, 0, 2)).doit(deep = False) == \
        3*Integral(a**2)
    assert summation(n*Integral(a**2), (n, 0, 2)) == 3*Integral(a**2)

def test_Product_doit():
    assert Product(n*Integral(a**2), (n, 1, 3)).doit() == 2 * a**9 / 9
    assert Product(n*Integral(a**2), (n, 1, 3)).doit(deep = False) == \
        6*Integral(a**2)**3
    assert product(n*Integral(a**2), (n, 1, 3)) == 6*Integral(a**2)**3

def test_Sum_interface():
    assert isinstance(Sum(0, (n, 0, 2)), Sum)
    assert Sum(nan, (n, 0, 2)) is nan
    assert Sum(nan, (n, 0, oo)) is nan
    assert Sum(0, (n, 0, 2)).doit() == 0
    assert isinstance(Sum(0, (n, 0, oo)), Sum)
    assert Sum(0, (n, 0, oo)).doit() == 0
    raises(ValueError, lambda: Sum(1))
    raises(ValueError, lambda: summation(1))

def test_eval_diff():
    assert Sum(x, (x, 1, 2)).diff(x) == 0
    assert Sum(x*y, (x, 1, 2)).diff(x) == 0
    assert Sum(x*y, (y, 1, 2)).diff(x) == Sum(y, (y, 1, 2))
    e = Sum(x*y, (x, 1, a))
    assert e.diff(a) == Derivative(e, a)
    assert Sum(x*y, (x, 1, 3), (a, 2, 5)).diff(y) == \
           Sum(x*y, (x, 1, 3), (a, 2, 5)).doit().diff(y) == \
           24

def test_hypersum():
    from sympy import simplify, sin, hyper
    assert simplify(summation(x**n/fac(n), (n, 1, oo))) == -1 + exp(x)
    assert summation((-1)**n * x**(2*n) / fac(2*n), (n, 0, oo)) == cos(x)
    assert simplify(summation((-1)**n*x**(2*n+1)/factorial(2*n+1),
                              (n, 3, oo))) \
           == -x + sin(x) + x**3/6 - x**5/120

    assert summation(1/(n+2)**3, (n, 1, oo)) == \
           -S(9)/8 + zeta(3)
    assert summation(1/n**4, (n, 1, oo)) == pi**4/90

    s = summation(x**n*n, (n, -oo, 0))
    assert s.is_Piecewise
    assert s.args[0].args[0] == -1/(x*(1 - 1/x)**2)
    assert s.args[0].args[1] == (abs(1/x) < 1)

    m = Symbol('n', integer=True, positive=True)
    assert summation(binomial(m, k), (k, 0, m)) == 2**m

def test_issue_1071():
    assert summation(1/factorial(k), (k, 0, oo)) == E

def test_is_zero():
    for func in [Sum, Product]:
        assert func(0, (x, 1, 1)).is_zero is True
        assert func(x, (x, 1, 1)).is_zero is None

def test_is_commutative():
    from sympy.physics.secondquant import NO, F, Fd
    m = Symbol('m', commutative=False)
    for f in (Sum, Product, Integral):
        assert f(z, (z, 1, 1)).is_commutative is True
        assert f(z*y, (z, 1, 6)).is_commutative is True
        assert f(m*x, (x, 1, 2)).is_commutative is False

        assert f(NO(Fd(x)*F(y))*z, (z, 1, 2)).is_commutative is False

def test_is_number():
    assert Sum(1, (x, 1, 1)).is_number is True
    assert Sum(1, (x, 1, x)).is_number is False
    assert Sum(0, (x, y, z)).is_number is True
    assert Sum(x, (y, 1, 2)).is_number is False
    assert Sum(x, (y, 1, 1)).is_number is False
    assert Sum(x, (x, 1, 2)).is_number is True
    assert Sum(x*y, (x, 1, 2), (y, 1, 3)).is_number is True

    assert Product(2, (x, 1, 1)).is_number is True
    assert Product(2, (x, 1, y)).is_number is False
    assert Product(0, (x, y, z)).is_number is True
    assert Product(1, (x, y, z)).is_number is True
    assert Product(x, (y, 1, x)).is_number is False
    assert Product(x, (y, 1, 2)).is_number is False
    assert Product(x, (y, 1, 1)).is_number is False
    assert Product(x, (x, 1, 2)).is_number is True

def test_free_symbols():
    for func in [Sum, Product]:
        assert func(1, (x, 1, 2)).free_symbols == set()
        assert func(0, (x, 1, y)).free_symbols == set()
        assert func(2, (x, 1, y)).free_symbols == set([y])
        assert func(x, (x, 1, 2)).free_symbols == set()
        assert func(x, (x, 1, y)).free_symbols == set([y])
        assert func(x, (y, 1, y)).free_symbols == set([x, y])
        assert func(x, (y, 1, 2)).free_symbols == set([x])
        assert func(x, (y, 1, 1)).free_symbols == set([x])
        assert func(x, (y, 1, z)).free_symbols == set([x, z])
        assert func(x, (x, 1, y), (y, 1, 2)).free_symbols == set()
        assert func(x, (x, 1, y), (y, 1, z)).free_symbols == set([z])
        assert func(x, (x, 1, y), (y, 1, y)).free_symbols == set([y])
        assert func(x, (y, 1, y), (y, 1, z)).free_symbols == set([x, z])
    assert Sum(1, (x, 1, y)).free_symbols == set([y])
    assert Product(1, (x, 1, y)).free_symbols == set()

@XFAIL
def test_issue_1072() :
    k = Symbol("k")
    assert summation(factorial(2*k + 1)/factorial(2*k), (k, 0, oo)) == oo

@XFAIL
def test_issue_3174():
    # when this passes, the doctests involving Sum in
    # is_constant can be unskipped
    assert Sum(x, (x, 1, n)).n(2, subs={n: 0}) == 1

@XFAIL
def test_issue_3175():
    assert Sum(x, (x, 1, 0)).doit() == 0
