from sympy import symbols, Symbol, sqrt, oo, re, nan, im, sign, I, E, log, \
        pi, arg, conjugate, expand, exp, sin, cos
from sympy.utilities.pytest import XFAIL


def test_re():

    x, y = symbols('xy')

    r = Symbol('r', real=True)

    assert re(nan) == nan

    assert re(oo) == oo
    assert re(-oo) == -oo

    assert re(0) == 0

    assert re(1) == 1
    assert re(-1) == -1

    assert re(E) == E
    assert re(-E) == -E

    assert re(x) == re(x)
    assert re(x*I) == -im(x)
    assert re(r*I) == 0
    assert re(r) == r

    assert re(x + y) == re(x + y)
    assert re(x + r) == re(x) + r

    assert re(re(x)) == re(x)

    assert re(2 + I) == 2
    assert re(x + I) == re(x)

    assert re(x + y*I) == re(x) - im(y)
    assert re(x + r*I) == re(x)

    assert re(log(2*I)) == log(2)

    assert re((2+I)**2).expand(complex=True) == 3

def test_im():

    x, y = symbols('xy')

    r = Symbol('r', real=True)

    assert im(nan) == nan

    assert im(oo*I) == oo
    assert im(-oo*I) == -oo

    assert im(0) == 0

    assert im(1) == 0
    assert im(-1) == 0

    assert im(E*I) == E
    assert im(-E*I) == -E

    assert im(x) == im(x)
    assert im(x*I) == re(x)
    assert im(r*I) == r
    assert im(r) == 0

    assert im(x + y) == im(x + y)
    assert im(x + r) == im(x)
    assert im(x + r*I) == im(x) + r

    assert im(im(x)*I) == im(x)

    assert im(2 + I) == 1
    assert im(x + I) == im(x) + 1

    assert im(x + y*I) == im(x) + re(y)
    assert im(x + r*I) == im(x) + r

    assert im(log(2*I)) == pi/2

    assert im((2+I)**2).expand(complex=True) == 4

def test_sign():
    assert sign(1.2) == 1
    assert sign(-1.2) == -1
    assert sign(0) == 0

def test_abs():
    x, y = symbols('xy')
    assert abs(0) == 0
    assert abs(1) == 1
    assert abs(-1)== 1
    assert abs(x).diff(x) == sign(x)
    x = Symbol('x',real=True)
    n = Symbol('n',integer=True)
    assert x**(2*n) == abs(x)**(2*n)

def test_abs_real():
    # test some properties of abs that only apply
    # to real numbers
    x = Symbol('x', complex=True)
    assert sqrt(x**2) != abs(x)
    assert abs(x**2) != x**2

    x = Symbol('x', real=True)
    assert sqrt(x**2) == abs(x)
    assert abs(x**2) == x**2

def test_abs_properties():
    x = Symbol('x')
    assert abs(x).is_real == True
    assert abs(x).is_positive == None
    assert abs(x).is_nonnegative == True

    w = Symbol('w', complex=True, zero=False)
    assert abs(w).is_real == True
    assert abs(w).is_positive == True
    assert abs(w).is_zero == False

    q = Symbol('q', positive=True)
    assert abs(q).is_real == True
    assert abs(q).is_positive == True
    assert abs(q).is_zero == False

def test_arg():
    assert arg(0) == nan
    assert arg(1) == 0
    assert arg(-1) == pi
    assert arg(I) == pi/2
    assert arg(-I) == -pi/2
    assert arg(1+I) == pi/4
    assert arg(-1+I) == 3*pi/4
    assert arg(1-I) == -pi/4

    p = Symbol('p', positive=True)
    assert arg(p) == 0

    n = Symbol('n', negative=True)
    assert arg(n) == pi

def test_conjugate():
    a = Symbol('a', real=True)
    assert conjugate(a) == a
    assert conjugate(I*a) == -I*a

    x, y = symbols('xy')
    assert conjugate(conjugate(x)) == x
    assert conjugate(x + y) == conjugate(x) + conjugate(y)
    assert conjugate(x - y) == conjugate(x) - conjugate(y)
    assert conjugate(x * y) == conjugate(x) * conjugate(y)
    assert conjugate(x / y) == conjugate(x) / conjugate(y)
    assert conjugate(-x) == -conjugate(x)

def test_issue936():
    x = Symbol('x')
    assert abs(x).expand(trig=True)     == abs(x)
    assert sign(x).expand(trig=True)    == sign(x)
    assert arg(x).expand(trig=True)     == arg(x)

