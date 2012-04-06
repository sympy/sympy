from sympy import bernoulli, Symbol, symbols, Sum, harmonic, Rational, oo, \
                    zoo, pi, I, bell, fibonacci, lucas, euler, catalan, \
                    binomial, gamma, sqrt, hyper, log, polygamma, diff, \
                    Expr, sympify

from sympy.utilities.pytest import XFAIL

x = Symbol('x')

def test_bernoulli():
    assert bernoulli(0) == 1
    assert bernoulli(1) == Rational(-1,2)
    assert bernoulli(2) == Rational(1,6)
    assert bernoulli(3) == 0
    assert bernoulli(4) == Rational(-1,30)
    assert bernoulli(5) == 0
    assert bernoulli(6) == Rational(1,42)
    assert bernoulli(7) == 0
    assert bernoulli(8) == Rational(-1,30)
    assert bernoulli(10) == Rational(5,66)
    assert bernoulli(1000001) == 0

    assert bernoulli(0, x) == 1
    assert bernoulli(1, x) == x-Rational(1,2)
    assert bernoulli(2, x) == x**2-x+Rational(1,6)
    assert bernoulli(3, x) == x**3 - (3*x**2)/2 + x/2

    # Should be fast; computed with mpmath
    b = bernoulli(1000)
    assert b.p % 10**10  == 7950421099
    assert b.q == 342999030

    b = bernoulli(10**6, evaluate=False).evalf()
    assert str(b) == '-2.23799235765713e+4767529'


def test_fibonacci():
    assert [fibonacci(n) for n in range(-3, 5)] == [2, -1, 1, 0, 1, 1, 2, 3]
    assert fibonacci(100) == 354224848179261915075
    assert [lucas(n) for n in range(-3, 5)] == [-4, 3, -1, 2, 1, 3, 4, 7]
    assert lucas(100) == 792070839848372253127

    assert fibonacci(1, x) == 1
    assert fibonacci(2, x) == x
    assert fibonacci(3, x) == x**2 + 1
    assert fibonacci(4, x) == x**3 + 2*x

def test_bell():
    assert [bell(n) for n in range(8)] == [1, 1, 2, 5, 15, 52, 203, 877]

    assert bell(0, x) == 1
    assert bell(1, x) == x
    assert bell(2, x) == x**2 + x
    assert bell(5, x) == x**5 + 10*x**4 + 25*x**3 + 15*x**2 + x

    X = symbols('x:6')
    # X = (x0, x1, .. x5)
    # at the same time: X[1] = x1, X[2] = x2 for standard readablity.
    # but we must supply zero-based indexed object X[1:] = (x1, .. x5)

    assert bell(6, 2, X[1:]) == 6*X[5]*X[1] + 15*X[4]*X[2] + 10*X[3]**2
    assert bell(6, 3, X[1:]) == 15*X[4]*X[1]**2 + 60*X[3]*X[2]*X[1] + 15*X[2]**3

    X = (1, 10, 100, 1000, 10000)
    assert bell(6, 2, X) == (6 + 15 + 10)*10000

    X = (1, 2, 3, 3, 5)
    assert bell(6, 2, X) == 6*5 + 15*3*2 + 10*3**2

    X = (1, 2, 3, 5)
    assert bell(6, 3, X) == 15*5 + 60*3*2 + 15*2**3

def test_harmonic():
    assert harmonic(1,1) == 1
    assert harmonic(2,1) == Rational(3,2)
    assert harmonic(3,1) == Rational(11,6)
    assert harmonic(4,1) == Rational(25,12)
    assert harmonic(3,1) == harmonic(3)
    assert harmonic(3,5) == 1 + Rational(1,2**5) + Rational(1,3**5)
    assert harmonic(10,0) == 10
    assert harmonic(oo,1) == zoo
    assert harmonic(oo,2) == (pi**2)/6

def test_euler():
    assert euler(0) == 1
    assert euler(1) == 0
    assert euler(2) == -1
    assert euler(3) == 0
    assert euler(4) == 5
    assert euler(6) == -61
    assert euler(8) == 1385

    assert euler(20, evaluate=False) != 370371188237525

    n = Symbol('n', integer=True)
    assert euler(n) != -1
    assert euler(n).subs(n, 2) == -1

    assert euler(20).evalf() == 370371188237525.0
    assert euler(20, evaluate=False).evalf() == 370371188237525.0

    assert euler(n).rewrite(Sum) == euler(n)
    # XXX: Not sure what the guy who wrote this test was trying to do with the _j and _k stuff
    assert euler(2*n+1).rewrite(Sum) == 0

@XFAIL
def test_euler_failing():
    # depends on dummy variables being implemented http://code.google.com/p/sympy/issues/detail?id=2566
    assert euler(2*n).rewrite(Sum) ==  I*Sum(Sum((-1)**_j*2**(-_k)*I**(-_k)*(-2*_j + _k)**(2*n + 1)*binomial(_k, _j)/_k, (_j, 0, _k)), (_k, 1, 2*n + 1))

def test_catalan():
    assert catalan(1) == 1
    assert catalan(2) == 2
    assert catalan(3) == 5
    assert catalan(4) == 14

    assert catalan(x) == catalan(x)
    assert catalan(2*x).rewrite(binomial) == binomial(4*x, 2*x)/(2*x + 1)
    assert catalan(Rational(1,2)).rewrite(gamma) == 8/(3*pi)
    assert catalan(3*x).rewrite(gamma) == 4**(3*x)*gamma(3*x + Rational(1,2))/(sqrt(pi)*gamma(3*x + 2))
    assert catalan(x).rewrite(hyper) == hyper((-x + 1, -x), (2,), 1)

    assert diff(catalan(x),x) == (polygamma(0, x + Rational(1,2)) - polygamma(0, x + 2) + 2*log(2))*catalan(x)

    c = catalan(0.5).evalf()
    assert str(c) == '0.848826363156775'
