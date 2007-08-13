from sympy import *


def test_complex():
    a = Symbol("a", real=True)
    b = Symbol("b", real=True)
    e = (a+I*b)*(a-I*b)
    assert e.expand() == a**2+b**2
    assert sqrt(I) == (-1)**Rational(1,4)
    assert str(abs(a)) == "abs(a)"

def test_conjugate():
    a = Symbol("a", real=True)
    b = Symbol("b", real=True)
    x = Symbol('x')
    z = a + I*b
    zc = a - I*b
    assert conjugate(z) == zc
    assert conjugate(exp(z)) == exp(zc)
    assert conjugate(exp(I*x)) == exp(-I*conjugate(x))
    assert conjugate(z**5) == zc**5
    assert conjugate(abs(x)) == abs(x)
    assert conjugate(sign(x)) == sign(x)
    assert conjugate(sin(z)) == sin(zc)
    assert conjugate(cos(z)) == cos(zc)
    assert conjugate(tan(z)) == tan(zc)
    assert conjugate(cot(z)) == cot(zc)
    assert conjugate(sinh(z)) == sinh(zc)
    assert conjugate(cosh(z)) == cosh(zc)
    assert conjugate(tanh(z)) == tanh(zc)
    assert conjugate(coth(z)) == coth(zc)

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

    assert tan(I*x).expand(complex=True) == tanh(x) * I
    assert tan(x+I*y).expand(complex=True) == \
            ((sin(x)*cos(x) + I*cosh(y)*sinh(y)) / (cos(x)**2 + sinh(y)**2)).expand()

    assert sinh(I*x).expand(complex=True) == I * sin(x)
    assert sinh(x+I*y).expand(complex=True) == sinh(x)*cos(y) + \
            I * sin(y) * cosh(x)

    assert cosh(I*x).expand(complex=True) == cos(x)
    assert cosh(x+I*y).expand(complex=True) == cosh(x)*cos(y) + \
            I * sin(y) * sinh(x)

    assert tanh(I*x).expand(complex=True) == tan(x) * I
    assert tanh(x+I*y).expand(complex=True) == \
            ((sinh(x)*cosh(x) + I*cos(y)*sin(y)) / (sinh(x)**2 + cos(y)**2)).expand()


def test_pythoncomplex():
    x = Symbol("x")
    assert 4j*x == 4*x*I
    assert 4j*x != 4.0*x*I
    assert 4.1j*x != 4*x*I

def test_rootcomplex():
    assert ((1+I)**Rational(1,2)).expand(complex=True) == 2**Rational(1,4)*cos(Rational(1,2)*atan(1))+2**(Rational(1,4))*sin(Rational(1,2)*atan(1))*I
    assert (sqrt(-10)*I).as_real_imag() == (-sqrt(10), 0)
