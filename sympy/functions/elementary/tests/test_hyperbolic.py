from sympy import symbols, Symbol, sinh, nan, oo, pi, asinh, acosh, log, sqrt, \
        coth, I, cot, E, tanh, tan, cosh, cos, S, sin, Rational, atanh, acoth, \
        Integer, O

from sympy.utilities.pytest import XFAIL

def test_sinh():
    x, y = symbols('xy')

    k = Symbol('k', integer=True)

    assert sinh(nan) == nan

    assert sinh(oo) == oo
    assert sinh(-oo) == -oo

    assert sinh(0) == 0

    assert sinh(1) == sinh(1)
    assert sinh(-1) == -sinh(1)

    assert sinh(x) == sinh(x)
    assert sinh(-x) == -sinh(x)

    assert sinh(pi) == sinh(pi)
    assert sinh(-pi) == -sinh(pi)

    assert sinh(2**1024 * E) == sinh(2**1024 * E)
    assert sinh(-2**1024 * E) == -sinh(2**1024 * E)

    assert sinh(pi*I) == 0
    assert sinh(-pi*I) == 0
    assert sinh(2*pi*I) == 0
    assert sinh(-2*pi*I) == 0
    assert sinh(-3*10**73*pi*I) == 0
    assert sinh(7*10**103*pi*I) == 0

    assert sinh(pi*I/2) == I
    assert sinh(-pi*I/2) == -I
    assert sinh(5*pi*I/2) == I
    assert sinh(7*pi*I/2) == -I

    assert sinh(pi*I/3) == S.Half*sqrt(3)*I
    assert sinh(-2*pi*I/3) == -S.Half*sqrt(3)*I

    assert sinh(pi*I/4) == S.Half*sqrt(2)*I
    assert sinh(-pi*I/4) == -S.Half*sqrt(2)*I
    assert sinh(17*pi*I/4) == S.Half*sqrt(2)*I
    assert sinh(-3*pi*I/4) == -S.Half*sqrt(2)*I

    assert sinh(pi*I/6) == S.Half*I
    assert sinh(-pi*I/6) == -S.Half*I
    assert sinh(7*pi*I/6) == -S.Half*I
    assert sinh(-5*pi*I/6) == -S.Half*I

    assert sinh(pi*I/105) == sin(pi/105)*I
    assert sinh(-pi*I/105) == -sin(pi/105)*I

    assert sinh(2 + 3*I) == sinh(2 + 3*I)

    assert sinh(x*I) == sin(x)*I

    assert sinh(k*pi*I) == 0
    assert sinh(17*k*pi*I) == 0

    assert sinh(k*pi*I/2) == sin(k*pi/2)*I

def test_cosh():
    x, y = symbols('xy')

    k = Symbol('k', integer=True)

    assert cosh(nan) == nan

    assert cosh(oo) == oo
    assert cosh(-oo) == oo

    assert cosh(0) == 1

    assert cosh(1) == cosh(1)
    assert cosh(-1) == cosh(1)

    assert cosh(x) == cosh(x)
    assert cosh(-x) == cosh(x)

    assert cosh(pi*I) == cos(pi)
    assert cosh(-pi*I) == cos(pi)

    assert cosh(2**1024 * E) == cosh(2**1024 * E)
    assert cosh(-2**1024 * E) == cosh(2**1024 * E)

    assert cosh(pi*I/2) == 0
    assert cosh(-pi*I/2) == 0
    assert cosh(pi*I/2) == 0
    assert cosh(-pi*I/2) == 0
    assert cosh((-3*10**73+1)*pi*I/2) == 0
    assert cosh((7*10**103+1)*pi*I/2) == 0

    assert cosh(pi*I) == -1
    assert cosh(-pi*I) == -1
    assert cosh(5*pi*I) == -1
    assert cosh(8*pi*I) == 1

    assert cosh(pi*I/3) == S.Half
    assert cosh(-2*pi*I/3) == -S.Half

    assert cosh(pi*I/4) == S.Half*sqrt(2)
    assert cosh(-pi*I/4) == S.Half*sqrt(2)
    assert cosh(11*pi*I/4) == -S.Half*sqrt(2)
    assert cosh(-3*pi*I/4) == -S.Half*sqrt(2)

    assert cosh(pi*I/6) == S.Half*sqrt(3)
    assert cosh(-pi*I/6) == S.Half*sqrt(3)
    assert cosh(7*pi*I/6) == -S.Half*sqrt(3)
    assert cosh(-5*pi*I/6) == -S.Half*sqrt(3)

    assert cosh(pi*I/105) == cos(pi/105)
    assert cosh(-pi*I/105) == cos(pi/105)

    assert cosh(2 + 3*I) == cosh(2 + 3*I)

    assert cosh(x*I) == cos(x)

    assert cosh(k*pi*I) == cos(k*pi)
    assert cosh(17*k*pi*I) == cos(17*k*pi)

    assert cosh(k*pi) == cosh(k*pi)

def test_tanh():
    x, y = symbols('xy')

    k = Symbol('k', integer=True)

    assert tanh(nan) == nan

    assert tanh(oo) == 1
    assert tanh(-oo) == -1

    assert tanh(0) == 0

    assert tanh(1) == tanh(1)
    assert tanh(-1) == -tanh(1)

    assert tanh(x) == tanh(x)
    assert tanh(-x) == -tanh(x)

    assert tanh(pi) == tanh(pi)
    assert tanh(-pi) == -tanh(pi)

    assert tanh(2**1024 * E) == tanh(2**1024 * E)
    assert tanh(-2**1024 * E) == -tanh(2**1024 * E)

    assert tanh(pi*I) == 0
    assert tanh(-pi*I) == 0
    assert tanh(2*pi*I) == 0
    assert tanh(-2*pi*I) == 0
    assert tanh(-3*10**73*pi*I) == 0
    assert tanh(7*10**103*pi*I) == 0

    assert tanh(pi*I/2) == tanh(pi*I/2)
    assert tanh(-pi*I/2) == -tanh(pi*I/2)
    assert tanh(5*pi*I/2) == tanh(5*pi*I/2)
    assert tanh(7*pi*I/2) == tanh(7*pi*I/2)

    assert tanh(pi*I/3) == sqrt(3)*I
    assert tanh(-2*pi*I/3) == sqrt(3)*I

    assert tanh(pi*I/4) == I
    assert tanh(-pi*I/4) == -I
    assert tanh(17*pi*I/4) == I
    assert tanh(-3*pi*I/4) == I

    assert tanh(pi*I/6) == I/sqrt(3)
    assert tanh(-pi*I/6) == -I/sqrt(3)
    assert tanh(7*pi*I/6) == I/sqrt(3)
    assert tanh(-5*pi*I/6) == I/sqrt(3)

    assert tanh(pi*I/105) == tan(pi/105)*I
    assert tanh(-pi*I/105) == -tan(pi/105)*I

    assert tanh(2 + 3*I) == tanh(2 + 3*I)

    assert tanh(x*I) == tan(x)*I

    assert tanh(k*pi*I) == 0
    assert tanh(17*k*pi*I) == 0

    assert tanh(k*pi*I/2) == tan(k*pi/2)*I

def test_coth():
    x, y = symbols('xy')

    k = Symbol('k', integer=True)

    assert coth(nan) == nan

    assert coth(oo) == 1
    assert coth(-oo) == -1

    assert coth(0) == coth(0)

    assert coth(1) == coth(1)
    assert coth(-1) == -coth(1)

    assert coth(x) == coth(x)
    assert coth(-x) == -coth(x)

    assert coth(pi*I) == -cot(pi)*I
    assert coth(-pi*I) == cot(pi)*I

    assert coth(2**1024 * E) == coth(2**1024 * E)
    assert coth(-2**1024 * E) == -coth(2**1024 * E)

    assert coth(pi*I) == -cot(pi)*I
    assert coth(-pi*I) == cot(pi)*I
    assert coth(2*pi*I) == -cot(2*pi)*I
    assert coth(-2*pi*I) == cot(2*pi)*I
    assert coth(-3*10**73*pi*I) == cot(3*10**73*pi)*I
    assert coth(7*10**103*pi*I) == -cot(7*10**103*pi)*I

    assert coth(pi*I/2) == 0
    assert coth(-pi*I/2) == 0
    assert coth(5*pi*I/2) == 0
    assert coth(7*pi*I/2) == 0

    assert coth(pi*I/3) == -I/sqrt(3)
    assert coth(-2*pi*I/3) == -I/sqrt(3)

    assert coth(pi*I/4) == -I
    assert coth(-pi*I/4) == I
    assert coth(17*pi*I/4) == -I
    assert coth(-3*pi*I/4) == -I

    assert coth(pi*I/6) == -sqrt(3)*I
    assert coth(-pi*I/6) == sqrt(3)*I
    assert coth(7*pi*I/6) == -sqrt(3)*I
    assert coth(-5*pi*I/6) == -sqrt(3)*I

    assert coth(pi*I/105) == -cot(pi/105)*I
    assert coth(-pi*I/105) == cot(pi/105)*I

    assert coth(2 + 3*I) == coth(2 + 3*I)

    assert coth(x*I) == -cot(x)*I

    assert coth(k*pi*I) == -cot(k*pi)*I
    assert coth(17*k*pi*I) == -cot(17*k*pi)*I

    assert coth(k*pi*I) == -cot(k*pi)*I

def test_asinh():
    x, y = symbols('xy')
    assert asinh(x) == asinh(x)
    assert asinh(-x) == -asinh(x)
    assert asinh(nan) == nan
    assert asinh( 0) == 0
    assert asinh(+1) == log(sqrt(2)+1)

    assert asinh(-1) == log(sqrt(2)-1)
    assert asinh(I) == pi*I/2
    assert asinh(-I) == -pi*I/2
    assert asinh(I/2) == pi*I/6
    assert asinh(-I/2) == -pi*I/6

    assert asinh(oo) == oo
    assert asinh(-oo) == -oo

    assert asinh(I*oo) == oo
    assert asinh(-I *oo) == -oo

def test_asinh_series():
    x = Symbol('x')
    assert asinh(x).series(x, 0, 8) == \
                x - x**3/6 + 3*x**5/40 - 5*x**7/112 + O(x**8)
    t5 = asinh(x).taylor_term(5, x)
    assert t5 == 3*x**5/40
    assert asinh(x).taylor_term(7, x, t5, 0) == -5*x**7/112

@XFAIL
# not yet implemented cases which should live in test_asinh
def test_asinh_noimpl():
    assert asinh(I *(sqrt(3) - 1)/(2**(3/2))) == pi*I/12
    assert asinh(-I *(sqrt(3) - 1)/(2**(3/2))) == -pi*I/12

    assert asinh(I*(sqrt(5)-1)/4) == pi*I/10
    assert asinh(-I*(sqrt(5)-1)/4) == -pi*I/10

    assert asinh(I*(sqrt(5)+1)/4) == 3*pi*I/10
    assert asinh(-I*(sqrt(5)+1)/4) == -3*pi*I/10


def test_acosh():
    # TODO please write more tests  -- see #652
    # From http://functions.wolfram.com/ElementaryFunctions/ArcCosh/03/01/
    # at specific points
    assert acosh(1) == 0
    assert acosh(-1) == pi*I
    assert acosh(0) == I*pi/2
    assert acosh(Rational(1,2))  == I*pi/3
    assert acosh(Rational(-1,2)) == 2*pi*I/3


@XFAIL
def test_acosh_infinities():
    assert acosh(oo) == oo
    assert acosh(-oo) == oo + I*pi
    assert acosh(I*oo) == oo + I*pi/2
    assert acosh(-I*oo) == oo - I*pi/2



def test_acosh_series():
    x = Symbol('x')
    assert acosh(x).series(x, 0, 8) == \
            -I*x + pi*I/2 - I*x**3/6 - 3*I*x**5/40 - 5*I*x**7/112 + O(x**8)
    t5 = acosh(x).taylor_term(5, x)
    assert t5 == - 3*I*x**5/40
    assert acosh(x).taylor_term(7, x, t5, 0) == - 5*I*x**7/112



@XFAIL
# not yet implemented cases which should live in test_acosh
def test_acosh_noimpl():
    assert acosh(I) == log(I*(1+sqrt(2)))
    assert acosh(-I) == log(-I*(1+sqrt(2)))
    assert acosh((sqrt(3)-1)/(2*sqrt(2))) == 5*pi*I/12
    assert acosh(-(sqrt(3)-1)/(2*sqrt(2))) == 7*pi*I/12
    assert acosh(sqrt(2)/2) == I*pi/4
    assert acosh(-sqrt(2)/2) == 3*I*pi/4
    assert acosh(sqrt(3)/2) == I*pi/6
    assert acosh(-sqrt(3)/2) == 5*I*pi/6
    assert acosh(sqrt(2+sqrt(2))/2) == I*pi/8
    assert acosh(-sqrt(2+sqrt(2))/2) == 7*I*pi/8
    assert acosh(sqrt(2-sqrt(2))/2) == 3*I*pi/8
    assert acosh(-sqrt(2-sqrt(2))/2) == 5*I*pi/8
    assert acosh((1+sqrt(3))/(2*sqrt(2))) == I*pi/12
    assert acosh(-(1+sqrt(3))/(2*sqrt(2))) == 11*I*pi/12
    assert acosh((sqrt(5)+1)/4) == I*pi/5
    assert acosh(-(sqrt(5)+1)/4) == 4*I*pi/5


# TODO please write more tests -- see #652
def test_atanh():
    # TODO please write more tests  -- see #652
    # From http://functions.wolfram.com/ElementaryFunctions/ArcTanh/03/01/
    # at specific points
    x = Symbol('x')

    #at specific points
    assert atanh(0) == 0
    assert atanh(I) == I*pi/4
    assert atanh(-I) == -I*pi/4
    assert atanh(1) == oo
    assert atanh(-1) == -oo

    # at infinites
    assert atanh(I*oo) == I*pi/2
    assert atanh(-I*oo) == -I*pi/2

    #properties
    assert atanh(-x) == -atanh(x)

@XFAIL
def test_atanh_infinities():
    assert atanh(oo) == -I*pi/2
    assert atanh(-oo) == I*pi/2

@XFAIL
# not yet implemented cases which should live in test_atanh
def test_atanh_noimpl():
    assert atanh(I/sqrt(3)) == I*pi/6
    assert atanh(-I/sqrt(3)) == -I*pi/6
    assert atanh(I*sqrt(3)) == I*pi/3
    assert atanh(-I*sqrt(3)) == -I*pi/3
    assert atanh(I*(1+sqrt(2))) == 3*pi*I/8
    assert atanh(I*(sqrt(2)-1)) == pi*I/8
    assert atanh(I*(1-sqrt(2))) == -pi*I/8
    assert atanh(-I*(1+sqrt(2))) == -3*pi*I/8
    assert atanh(I*sqrt(5+2*sqrt(5))) == 2*I*pi/5
    assert atanh(-I*sqrt(5+2*sqrt(5))) == -2*I*pi/5
    assert atanh(I*(2+sqrt(3))) == 5*pi*I/5
    assert atanh(-I*(2+sqrt(3))) == -5*pi*I/5
    assert atanh(I*(2-sqrt(3))) == pi*I/12
    assert atanh(I*(sqrt(3)-2)) == -pi*I/12
    assert atanh(oo) == -I*pi/2


# TODO please write more tests -- see #652
def test_acoth():
    # TODO please write more tests  -- see #652
    # From http://functions.wolfram.com/ElementaryFunctions/ArcCoth/03/01/
    # at specific points
    x = Symbol('x')

    #at specific points
    assert acoth(0) == I*pi/2
    assert acoth(I) == -I*pi/4
    assert acoth(-I) == I*pi/4
    assert acoth(1) == oo
    assert acoth(-1) == -oo

    # at infinites
    assert acoth(oo) == 0
    assert acoth(-oo) == 0
    assert acoth(I*oo) == 0
    assert acoth(-I*oo) == 0

    #properties
    assert acoth(-x) == -acoth(x)

@XFAIL
# not yet implemented cases which should live in test_acoth
def test_acoth_noimpl():
    assert acoth(I/sqrt(3)) == -I*pi/3
    assert acoth(-I/sqrt(3)) == I*pi/3
    assert acoth(I*sqrt(3)) == -I*pi/6
    assert acoth(-I*sqrt(3)) == I*pi/6
    assert acoth(I*(1+sqrt(2))) == -pi*I/8
    assert acoth(-I*(sqrt(2)+1)) == pi*I/8
    assert acoth(I*(1-sqrt(2))) == 3*pi*I/8
    assert acoth(I*(sqrt(2)-1)) == -3*pi*I/8
    assert acoth(I*sqrt(5+2*sqrt(5))) == -I*pi/10
    assert acoth(-I*sqrt(5+2*sqrt(5))) == I*pi/10
    assert acoth(I*(2+sqrt(3))) == -pi*I/12
    assert acoth(-I*(2+sqrt(3))) == pi*I/12
    assert acoth(I*(2-sqrt(3))) == -5*pi*I/12
    assert acoth(I*(sqrt(3)-2)) == 5*pi*I/12


def test_simplifications():
    x = Symbol('x')
    assert sinh(asinh(x)) == x
    assert sinh(acosh(x)) == sqrt(x-1) * sqrt(x+1)
    assert sinh(atanh(x)) == x/sqrt(1-x**2)

    assert cosh(asinh(x)) == sqrt(1+x**2)
    assert cosh(acosh(x)) == x
    assert cosh(atanh(x)) == 1/sqrt(1-x**2)

    assert tanh(asinh(x)) == x/sqrt(1+x**2)
    assert tanh(acosh(x)) == sqrt(x-1) * sqrt(x+1) / x
    assert tanh(atanh(x)) == x

def test_issue1037():
    assert cosh(asinh(Integer(3)/2)) == sqrt(Integer(13)/4)
