from sympy import Symbol, Rational, cos, sin, exp

def testdiff():
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    p = Rational(5)
    e = a*b + b**p
    assert e.diff(a) == b
    assert e.diff(b) == a + 5*b**4
    assert e.diff(b).diff(a) == Rational(1)
    e = a*(b + c)
    assert e.diff(a) == b + c
    assert e.diff(b) == a
    assert e.diff(b).diff(a) == Rational(1)
    e = c**p
    assert e.diff(c, 6) == Rational(0)
    assert e.diff(c, 5) == Rational(120)
    e = c**Rational(2)
    assert e.diff(c) == 2*c
    e = a*b*c
    assert e.diff(c) == a*b

def testdiff2():
    n3 = Rational(3)
    n2 = Rational(2)
    n6 = Rational(6)
    x = Symbol("x")
    c = Symbol("c")
    e = n3*(-n2 + x**n2)*cos(x) + x*(-n6 + x**n2)*sin(x)
    assert e == 3*(-2 + x**2)*cos(x) + x*(-6 + x**2)*sin(x)
    assert e.diff(x).expand() == x**3*cos(x)

    e = (x + 1)**3
    assert e.diff(x) == 3*(x + 1)**2
    e = x*(x + 1)**3
    assert e.diff(x) == (x + 1)**3 + 3*x*(x + 1)**2
    e = 2*exp(x*x)*x
    assert e.diff(x) == 2*exp(x**2) + 4*x**2*exp(x**2)
