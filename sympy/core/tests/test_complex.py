from sympy import *


def test_complex():
    a = Symbol("a")
    b = Symbol("b")
    e = (a+I*b)*(a-I*b)
    assert e.expand() == a**2+b**2
    assert (a+I*b).conjugate() == conjugate(a)-I*conjugate(b)
    assert str(abs(a)) == "abs(a)"

def test_abs1():
    a=Symbol("a", real=True)
    b=Symbol("b", real=True)
    assert abs(a) == abs(a)
    assert abs(-a) == abs(a)
    # XXX this isn't working
    #assert abs(a+I*b) == sqrt(a*a+b*b)

def test_abs2():
    a=Symbol("a", real=False)
    b=Symbol("b", real=False)
    assert abs(a) != a
    assert abs(-a) != a
    assert abs(a+I*b) != sqrt(a*a+b*b)

def _test_evalc():
    # XXX evalc() methods currently do not exist
    x=Symbol("x", real=True)
    y=Symbol("y", real=True)
    assert ((x+I*y)**2).evalc() == x**2+2*I*x*y - y**2

    assert exp(I*x) != cos(x)+I*sin(x)
    assert exp(I*x).evalc() == cos(x)+I*sin(x)

    assert exp(I*x+y).evalc() == exp(y)*cos(x)+I*sin(x)*exp(y)

    assert sin(I*x).evalc() == I * (exp(x)/2-exp(-x)/2)
    assert sin(x+I*y).evalc() == sin(x)*(exp(y)/2+exp(-y)/2) + \
            I * (exp(y)/2-exp(-y)/2) * cos(x)

    assert cos(I*x).evalc() == (exp(x)/2+exp(-x)/2)
    assert cos(x+I*y).evalc() == cos(x)*(exp(y)/2+exp(-y)/2) - \
            I * (exp(y)/2-exp(-y)/2) * sin(x)


def test_pythoncomplex():
    x = Symbol("x")
    assert 4j*x == 4*x*I
    assert 4j*x != 4.0*x*I
    assert 4.1j*x != 4*x*I

def _test_rootcomplex():
    # XXX evalc() methods currently do not exist
    assert ((1+I)**Rational(1,2)).evalc() == 2**Rational(1,4)*cos(Rational(1,2)*atan(1))+2**(Rational(1,4))*sin(Rational(1,2)*atan(1))*I
    assert (sqrt(-10)*I).get_re_im() == (-sqrt(10), 0)
