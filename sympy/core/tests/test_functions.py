from sympy import *


def test_func():
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    p = Rational(5)
    e = a*b + sin(b**p)
    assert e == a*b + sin(b**5)
    assert e.diff(a) == b
    assert e.diff(b) == a+5*b**4*cos(b**5)
    e = tan(c)
    assert e == tan(c)
    assert e.diff(c) in [cos(c)**(-2),1 + sin(c)**2/cos(c)**2, 1 + tan(c)**2]
    e = c*log(c)-c
    assert e == -c+c*log(c)
    assert e.diff(c) == log(c)
    e = log(sin(c))
    assert e == log(sin(c))
    assert e.diff(c) == sin(c)**(-1)*cos(c)
    assert e.diff(c) != cos(c)**(-1)*sin(c)
    assert e.diff(c) != sin(c)**(-2)*cos(c)
    assert e.diff(c) != sin(c)**(-3)*cos(c)
    t = Rational(2)
    e = (t**a/log(t))
    assert e == 2**a*log(Rational(2))**(-1)
    assert e.diff(a) == 2**a

def test_log():
    assert log(2) > 0
    assert log(1).is_zero
    assert log(0.5).is_negative == True

def test_exp_log():
    x = Symbol("x", real=True)
    assert log(exp(x)) == x
    assert exp(log(x)) == x

def test_log_expansion():
    x = Symbol("x", real=True)
    y = Symbol("y", real=True)
    #assert log(x*y) != log(x)+log(y)
    #assert log(x**2) != 2*log(x)
    assert log(x*y).expand() == log(x)+log(y)
    assert log(x**2).expand() == 2*log(x)
    assert (log(x**-5)**-1).expand() == -1/log(x)/5

def test_log_hashing_bug():
    x = Symbol("y")
    assert x != log(log(x))
    assert hash(x) != hash(log(log(x)))
    assert log(x) != log(log(log(x)))

    e = 1/log(log(x)+log(log(x)))
    assert isinstance(e.base, Basic.ApplyLog)
    e = 1/log(log(x)+log(log(log(x))))
    assert isinstance(e.base, Basic.ApplyLog)

    x = Symbol("x")
    e = log(log(x))
    assert isinstance(e, Basic.ApplyLog)
    assert not isinstance(x, Basic.ApplyLog)
    assert hash(log(log(x))) != hash(x)
    assert e != x

def test_sign():
    assert sign(log(2)) == 1

def test_exp_bug():
    x = Symbol("x")
    assert exp(1*log(x)) == x

def test_exp_expand():
    x = Symbol("x")
    y = Symbol("y")

    e = exp(log(Rational(2))*(1+x)-log(Rational(2))*x)
    assert e.expand() == 2
    # XXX The RHS is automatically combined so that it equals the LHS
    #assert exp(x+y) != exp(x)*exp(y)
    assert exp(x+y).expand() == exp(x)*exp(y)

def test_pi():
    assert cos(pi)==-1
    assert cos(2*pi)==1
    assert sin(pi)==0
    assert sin(2*pi)==0

def test_bug1():
    x = Symbol("x")
    w = Symbol("w")

    e = sqrt(-log(w))
    assert e.subs(log(w),-x) == sqrt(x)

    e = sqrt(-5*log(w))
    assert e.subs(log(w),-x) == sqrt(5*x)

def test_Derivative():
    x = Symbol("x")
    e = Derivative(log(x),x)
    assert e == 1/x

def _test_invtrig(): # XXX No inverse trig yet
    x = Symbol("x")
    assert atan(0) == 0
    assert atan(x).diff(x) == 1/(1+x**2)

def test_general_function():
    nu = Function('nu', nofargs=1)
    x = Symbol("x")
    y = Symbol("y")

    e = nu(x)
    edx = e.diff(x)
    edy = e.diff(y)
    edxdx = e.diff(x).diff(x)
    edxdy = e.diff(x).diff(y)
    assert e == nu(x)
    assert edx != nu(x)
    assert edx == Derivative(nu(x), x)
    assert edy == 0
    assert edxdx == Derivative(Derivative(nu(x), x), x)
    assert edxdy == 0

def test_derivative_subs_bug():
    x = Symbol("x")
    l = Function('l', nofargs=1)
    n = Function('n', nofargs=1)

    e = Derivative(n(x), x)
    assert e.subs(n, l) != e
    assert e.subs(n, l) == Derivative(l(x), x)
    assert e.subs(n, -l) == Derivative(-l(x), x)

def test_derivative_linearity():
    x = Symbol("x")
    y = Symbol("y")
    n = Function('n', nofargs=1)

    assert Derivative(-n(x), x) == -Derivative(n(x), x)
    assert Derivative(8*n(x), x) == 8*Derivative(n(x), x)
    assert Derivative(8*n(x), x) != 7*Derivative(n(x), x)
    assert Derivative(8*n(x)*x, x) == 8*n(x) + 8*x*Derivative(n(x), x)
    assert Derivative(8*n(x)*y*x, x) == 8*y*n(x) + 8*y*x*Derivative(n(x), x)

def _test_combine():
    # XXX combine no longer exists
    x = Symbol("x")
    y = Symbol("y")
    assert exp(x)*exp(-x) != 1
    assert (exp(x)*exp(-x)).combine() == 1

    assert exp(x)**2 != exp(2*x)
    assert (exp(x)**2).combine() == exp(2*x)

    assert exp(x)*exp(-x/2)*exp(-x/2) != 1
    assert (exp(x)*exp(-x/2)*exp(-x/2)).combine() == 1

    assert (2*log(x)).combine() == log(x**2)
    assert exp(2*log(x)) != x**2
    assert exp(2*log(x)).combine() == x**2

    assert exp(x)*exp(-x)-1 !=0
    assert (exp(x)*exp(-x)-1).combine() == 0

    assert (2*exp(x)*exp(-x)).combine() == 2
    assert (x/exp(x)*exp(-x)).combine() == x*exp(-2*x)
