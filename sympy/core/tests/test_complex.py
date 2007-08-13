from sympy import *


def test_complex():
    a=Symbol("a", real=True)
    b=Symbol("b", real=True)
    x=Symbol('x')
    e = (a+I*b)*(a-I*b)
    assert e.expand() == a**2+b**2
    assert (a+I*b).conjugate() == a-I*b
    assert exp(a+I*b).conjugate() == exp(a-I*b)
    assert exp(I*x).conjugate() == exp(-I*conjugate(x))
    assert sqrt(I) == (-1)**Rational(1,4)
    assert str(abs(a)) == "abs(a)"

def test_abs1():
    a=Symbol("a", real=True)
    b=Symbol("b", real=True)
    assert abs(a) == abs(a)
    assert abs(-a) == abs(a)
    assert abs(a+I*b) == sqrt(a**2+b**2)

def test_abs2():
    a=Symbol("a", real=False)
    b=Symbol("b", real=False)
    assert abs(a) != a
    assert abs(-a) != a
    assert abs(a+I*b) != sqrt(a**2+b**2)

def test_evalc():
    x=Symbol("x", real=True)
    y=Symbol("y", real=True)
    assert ((x+I*y)**2).expand(complex=True) == x**2+2*I*x*y - y**2

    assert exp(I*x) != cos(x)+I*sin(x)
    assert exp(I*x).expand(complex=True) == cos(x)+I*sin(x)

    assert exp(I*x+y).expand(complex=True) == exp(y)*cos(x)+I*sin(x)*exp(y)

    assert sin(I*x).expand(complex=True) == I * sinh(x)
    assert sin(x+I*y).expand(complex=True) == sin(x)*cosh(y) + \
            I * sinh(y) * cos(x)

    assert cos(I*x).expand(complex=True) == cosh(x)
    assert cos(x+I*y).expand(complex=True) == cos(x)*cosh(y) - \
            I * sinh(y) * sin(x)


def test_pythoncomplex():
    x = Symbol("x")
    assert 4j*x == 4*x*I
    assert 4j*x != 4.0*x*I
    assert 4.1j*x != 4*x*I

def test_rootcomplex():
    assert ((1+I)**Rational(1,2)).expand(complex=True) == 2**Rational(1,4)*cos(Rational(1,2)*atan(1))+2**(Rational(1,4))*sin(Rational(1,2)*atan(1))*I
    assert (sqrt(-10)*I).as_real_imag() == (-sqrt(10), 0)
