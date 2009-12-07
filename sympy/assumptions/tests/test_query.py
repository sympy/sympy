from sympy.core import Symbol, symbols, S, Rational, Integer
from sympy.functions import exp, log, sin, cos, sign, re, im, sqrt
from sympy.assumptions import (Assume, global_assumptions, Q, ask,
    register_handler, remove_handler)
from sympy.assumptions.handlers import AskHandler
from sympy.utilities.pytest import raises, XFAIL

def test_int_1():
    z = 1
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == True
    assert ask(z, Q.rational)         == True
    assert ask(z, Q.real)             == True
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == False
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == True
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == True
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == True

def test_float_1():
    z = 1.0
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == True
    assert ask(z, Q.rational)         == True
    assert ask(z, Q.real)             == True
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == False
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == True
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == True
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == True

    z = 7.2123
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == False
    assert ask(z, Q.rational)         == True
    assert ask(z, Q.real)             == True
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == False
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == True
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == False
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

def test_zero_0():
    z = Integer(0)
    assert ask(z, Q.nonzero)          == False
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == True
    assert ask(z, Q.rational)         == True
    assert ask(z, Q.real)             == True
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == False
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == True
    assert ask(z, Q.odd)              == False
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == True
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

def test_negativeone():
    z = Integer(-1)
    assert ask(z, Q.nonzero)          == True
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == True
    assert ask(z, Q.rational)         == True
    assert ask(z, Q.real)             == True
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == False
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == False
    assert ask(z, Q.negative)         == True
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == True
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

def test_infinity():
    oo = S.Infinity
    assert ask(oo, Q.commutative)     == True
    assert ask(oo, Q.integer)         == False
    assert ask(oo, Q.rational)        == False
    assert ask(oo, Q.real)            == False
    assert ask(oo, Q.extended_real)   == True
    assert ask(oo, Q.complex)         == False
    assert ask(oo, Q.irrational)      == False
    assert ask(oo, Q.imaginary)       == False
    assert ask(oo, Q.positive)        == True
    assert ask(oo, Q.negative)        == False
    assert ask(oo, Q.even)            == False
    assert ask(oo, Q.odd)             == False
    assert ask(oo, Q.bounded)         == False
    assert ask(oo, Q.infinitesimal)   == False
    assert ask(oo, Q.prime)           == False
    assert ask(oo, Q.composite)       == False

def test_neg_infinity():
    mm = S.NegativeInfinity
    assert ask(mm, Q.commutative)    == True
    assert ask(mm, Q.integer)        == False
    assert ask(mm, Q.rational)       == False
    assert ask(mm, Q.real)           == False
    assert ask(mm, Q.extended_real)  == True
    assert ask(mm, Q.complex)        == False
    assert ask(mm, Q.irrational)     == False
    assert ask(mm, Q.imaginary)      == False
    assert ask(mm, Q.positive)       == False
    assert ask(mm, Q.negative)       == True
    assert ask(mm, Q.even)           == False
    assert ask(mm, Q.odd)            == False
    assert ask(mm, Q.bounded)        == False
    assert ask(mm, Q.infinitesimal)  == False
    assert ask(mm, Q.prime)          == False
    assert ask(mm, Q.composite)      == False

def test_nan():
    nan = S.NaN
    assert ask(nan, Q.commutative)   == True
    assert ask(nan, Q.integer)       == False
    assert ask(nan, Q.rational)      == False
    assert ask(nan, Q.real)          == False
    assert ask(nan, Q.extended_real) == False
    assert ask(nan, Q.complex)       == False
    assert ask(nan, Q.irrational)    == False
    assert ask(nan, Q.imaginary)     == False
    assert ask(nan, Q.positive)      == False
    assert ask(nan, Q.nonzero)       == True
    assert ask(nan, Q.even)          == False
    assert ask(nan, Q.odd)           == False
    assert ask(nan, Q.bounded)       == False
    assert ask(nan, Q.infinitesimal) == False
    assert ask(nan, Q.prime)         == False
    assert ask(nan, Q.composite)     == False

def test_Rational_number():
    r = Rational(3,4)
    assert ask(r, Q.commutative)      == True
    assert ask(r, Q.integer)          == False
    assert ask(r, Q.rational)         == True
    assert ask(r, Q.real)             == True
    assert ask(r, Q.complex)          == True
    assert ask(r, Q.irrational)       == False
    assert ask(r, Q.imaginary)        == False
    assert ask(r, Q.positive)         == True
    assert ask(r, Q.negative)         == False
    assert ask(r, Q.even)             == False
    assert ask(r, Q.odd)              == False
    assert ask(r, Q.bounded)          == True
    assert ask(r, Q.infinitesimal)    == False
    assert ask(r, Q.prime)            == False
    assert ask(r, Q.composite)        == False

    r = Rational(1,4)
    assert ask(r, Q.positive)         == True
    assert ask(r, Q.negative)         == False

    r = Rational(5,4)
    assert ask(r, Q.negative)         == False
    assert ask(r, Q.positive)         == True

    r = Rational(5,3)
    assert ask(r, Q.positive)         == True
    assert ask(r, Q.negative)         == False

    r = Rational(-3,4)
    assert ask(r, Q.positive)         == False
    assert ask(r, Q.negative)         == True

    r = Rational(-1,4)
    assert ask(r, Q.positive)         == False
    assert ask(r, Q.negative)         == True

    r = Rational(-5,4)
    assert ask(r, Q.negative)         == True
    assert ask(r, Q.positive)         == False

    r = Rational(-5,3)
    assert ask(r, Q.positive)         == False
    assert ask(r, Q.negative)         == True

def test_sqrt_2():
    z = sqrt(2)
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == False
    assert ask(z, Q.rational)         == False
    assert ask(z, Q.real)             == True
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == True
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == True
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == False
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

def test_pi():
    z = S.Pi
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == False
    assert ask(z, Q.rational)         == False
    assert ask(z, Q.real)             == True
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == True
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == True
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == False
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

    z = S.Pi + 1
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == False
    assert ask(z, Q.rational)         == False
    assert ask(z, Q.real)             == True
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == True
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == True
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == False
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

    z = 2*S.Pi
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == False
    assert ask(z, Q.rational)         == False
    assert ask(z, Q.real)             == True
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == True
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == True
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == False
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

    z = S.Pi ** 2
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == False
    assert ask(z, Q.rational)         == False
    assert ask(z, Q.real)             == True
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == True
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == True
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == False
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

    z = (1+S.Pi) ** 2
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == False
    assert ask(z, Q.rational)         == False
    assert ask(z, Q.real)             == True
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == True
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == True
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == False
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

def test_E():
    z = S.Exp1
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == False
    assert ask(z, Q.rational)         == False
    assert ask(z, Q.real)             == True
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == True
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == True
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == False
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

def test_I():
    I = S.ImaginaryUnit
    z = I
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == False
    assert ask(z, Q.rational)         == False
    assert ask(z, Q.real)             == False
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == False
    assert ask(z, Q.imaginary)        == True
    assert ask(z, Q.positive)         == False
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == False
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

    z = 1 + I
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == False
    assert ask(z, Q.rational)         == False
    assert ask(z, Q.real)             == False
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == False
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == False
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == False
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

    z = I*(1+I)
    assert ask(z, Q.commutative)      == True
    assert ask(z, Q.integer)          == False
    assert ask(z, Q.rational)         == False
    assert ask(z, Q.real)             == False
    assert ask(z, Q.complex)          == True
    assert ask(z, Q.irrational)       == False
    assert ask(z, Q.imaginary)        == False
    assert ask(z, Q.positive)         == False
    assert ask(z, Q.negative)         == False
    assert ask(z, Q.even)             == False
    assert ask(z, Q.odd)              == False
    assert ask(z, Q.bounded)          == True
    assert ask(z, Q.infinitesimal)    == False
    assert ask(z, Q.prime)            == False
    assert ask(z, Q.composite)        == False

def test_bounded():
    x, y = symbols('xy')
    assert ask(x, Q.bounded) == False
    assert ask(x, Q.bounded, Assume(x, Q.bounded)) == True
    assert ask(x, Q.bounded, Assume(y, Q.bounded)) == False
    assert ask(x, Q.bounded, Assume(x, Q.complex)) == False

    assert ask(x+1, Q.bounded) == False
    assert ask(x+1, Q.bounded, Assume(x, Q.bounded)) == True
    assert ask(x+y, Q.bounded) == None
    assert ask(x+y, Q.bounded, Assume(x, Q.bounded)) == False
    assert ask(x+1, Q.bounded, Assume(x, Q.bounded) & \
                Assume(y, Q.bounded)) == True

    assert ask(2*x, Q.bounded) == False
    assert ask(2*x, Q.bounded, Assume(x, Q.bounded)) == True
    assert ask(x*y, Q.bounded) == None
    assert ask(x*y, Q.bounded, Assume(x, Q.bounded)) == False
    assert ask(x*y, Q.bounded, Assume(x, Q.bounded) & \
                Assume(y, Q.bounded)) == True

    assert ask(x**2, Q.bounded) == False
    assert ask(2**x, Q.bounded) == False
    assert ask(2**x, Q.bounded, Assume(x, Q.bounded)) == True
    assert ask(x**x, Q.bounded) == False
    assert ask(Rational(1,2) ** x, Q.bounded) == True
    assert ask(x ** Rational(1,2), Q.bounded) == False

    # sign function
    assert ask(sign(x), Q.bounded) == True
    assert ask(sign(x), Q.bounded, Assume(x, Q.bounded, False)) == True

    # exponential functions
    assert ask(log(x), Q.bounded) == False
    assert ask(log(x), Q.bounded, Assume(x, Q.bounded)) == True
    assert ask(exp(x), Q.bounded) == False
    assert ask(exp(x), Q.bounded, Assume(x, Q.bounded)) == True
    assert ask(exp(2), Q.bounded) == True

    # trigonometric functions
    assert ask(sin(x), Q.bounded) == True
    assert ask(sin(x), Q.bounded, Assume(x, Q.bounded, False)) == True
    assert ask(cos(x), Q.bounded) == True
    assert ask(cos(x), Q.bounded, Assume(x, Q.bounded, False)) == True
    assert ask(2*sin(x), Q.bounded) == True
    assert ask(sin(x)**2, Q.bounded) == True
    assert ask(cos(x)**2, Q.bounded) == True
    assert ask(cos(x) + sin(x), Q.bounded) == True

@XFAIL
def test_bounded_xfail():
    """We need to support relations in ask for this to work"""
    x = Symbol('x')
    assert ask(sin(x)**x, Q.bounded) == True
    assert ask(cos(x)**x, Q.bounded) == True
    assert ask(sin(x) ** x, Q.bounded) == True

def test_commutative():
    """By default objects are Q.commutative that is why it returns True
    for both key=True and key=False"""
    x, y = symbols('xy')
    assert ask(x, Q.commutative) == True
    assert ask(x, Q.commutative, Assume(x, Q.commutative, False)) == False
    assert ask(x, Q.commutative, Assume(x, Q.complex)) == True
    assert ask(x, Q.commutative, Assume(x, Q.imaginary)) == True
    assert ask(x, Q.commutative, Assume(x, Q.real)) == True
    assert ask(x, Q.commutative, Assume(x, Q.positive)) == True
    assert ask(x, Q.commutative, Assume(y, Q.commutative, False))  == True

    assert ask(2*x, Q.commutative) == True
    assert ask(2*x, Q.commutative, Assume(x, Q.commutative, False)) == False

    assert ask(x + 1, Q.commutative) == True
    assert ask(x + 1, Q.commutative, Assume(x, Q.commutative, False)) == False

    assert ask(x**2, Q.commutative) == True
    assert ask(x**2, Q.commutative, Assume(x, Q.commutative, False)) == False

    assert ask(log(x), Q.commutative) == True

def test_complex():
    x, y = symbols('xy')
    assert ask(x, Q.complex) == None
    assert ask(x, Q.complex, Assume(x, Q.complex)) == True
    assert ask(x, Q.complex, Assume(y, Q.complex)) == None
    assert ask(x, Q.complex, Assume(x, Q.complex, False)) == False
    assert ask(x, Q.complex, Assume(x, Q.real)) == True
    assert ask(x, Q.complex, Assume(x, Q.real, False)) == None
    assert ask(x, Q.complex, Assume(x, Q.rational)) == True
    assert ask(x, Q.complex, Assume(x, Q.irrational)) == True
    assert ask(x, Q.complex, Assume(x, Q.positive)) == True
    assert ask(x, Q.complex, Assume(x, Q.imaginary)) == True

    # a+b
    assert ask(x+1, Q.complex, Assume(x, Q.complex)) == True
    assert ask(x+1, Q.complex, Assume(x, Q.real)) == True
    assert ask(x+1, Q.complex, Assume(x, Q.rational)) == True
    assert ask(x+1, Q.complex, Assume(x, Q.irrational)) == True
    assert ask(x+1, Q.complex, Assume(x, Q.imaginary)) == True
    assert ask(x+1, Q.complex, Assume(x, Q.integer))  == True
    assert ask(x+1, Q.complex, Assume(x, Q.even))  == True
    assert ask(x+1, Q.complex, Assume(x, Q.odd))  == True
    assert ask(x+y, Q.complex, Assume(x, Q.complex) & Assume(y, Q.complex)) == True
    assert ask(x+y, Q.complex, Assume(x, Q.real) & Assume(y, Q.imaginary)) == True

    # a*x +b
    assert ask(2*x+1, Q.complex, Assume(x, Q.complex)) == True
    assert ask(2*x+1, Q.complex, Assume(x, Q.real)) == True
    assert ask(2*x+1, Q.complex, Assume(x, Q.positive)) == True
    assert ask(2*x+1, Q.complex, Assume(x, Q.rational)) == True
    assert ask(2*x+1, Q.complex, Assume(x, Q.irrational)) == True
    assert ask(2*x+1, Q.complex, Assume(x, Q.imaginary)) == True
    assert ask(2*x+1, Q.complex, Assume(x, Q.integer))  == True
    assert ask(2*x+1, Q.complex, Assume(x, Q.even))  == True
    assert ask(2*x+1, Q.complex, Assume(x, Q.odd))  == True

    # x**2
    assert ask(x**2, Q.complex, Assume(x, Q.complex)) == True
    assert ask(x**2, Q.complex, Assume(x, Q.real)) == True
    assert ask(x**2, Q.complex, Assume(x, Q.positive)) == True
    assert ask(x**2, Q.complex, Assume(x, Q.rational)) == True
    assert ask(x**2, Q.complex, Assume(x, Q.irrational)) == True
    assert ask(x**2, Q.complex, Assume(x, Q.imaginary)) == True
    assert ask(x**2, Q.complex, Assume(x, Q.integer))  == True
    assert ask(x**2, Q.complex, Assume(x, Q.even))  == True
    assert ask(x**2, Q.complex, Assume(x, Q.odd))  == True

    # 2**x
    assert ask(2**x, Q.complex, Assume(x, Q.complex)) == True
    assert ask(2**x, Q.complex, Assume(x, Q.real)) == True
    assert ask(2**x, Q.complex, Assume(x, Q.positive)) == True
    assert ask(2**x, Q.complex, Assume(x, Q.rational)) == True
    assert ask(2**x, Q.complex, Assume(x, Q.irrational)) == True
    assert ask(2**x, Q.complex, Assume(x, Q.imaginary)) == True
    assert ask(2**x, Q.complex, Assume(x, Q.integer))  == True
    assert ask(2**x, Q.complex, Assume(x, Q.even))  == True
    assert ask(2**x, Q.complex, Assume(x, Q.odd))  == True
    assert ask(x**y, Q.complex, Assume(x, Q.complex) & \
                     Assume(y, Q.complex)) == True

    # trigonometric expressions
    assert ask(sin(x), Q.complex) == True
    assert ask(sin(2*x + 1), Q.complex) == True
    assert ask(cos(x), Q.complex) == True
    assert ask(cos(2*x+1), Q.complex) == True

    # exponential
    assert ask(exp(x), Q.complex) == True
    assert ask(exp(x), Q.complex) == True

    # Q.complexes
    assert ask(abs(x), Q.complex) == True
    assert ask(re(x),  Q.complex) == True
    assert ask(im(x),  Q.complex) == True

def test_even():
    x, y, z, t = symbols('x y z t')
    assert ask(x, Q.even) == None
    assert ask(x, Q.even, Assume(x, Q.integer)) == None
    assert ask(x, Q.even, Assume(x, Q.integer, False)) == False
    assert ask(x, Q.even, Assume(x, Q.rational)) == None
    assert ask(x, Q.even, Assume(x, Q.positive)) == None

    assert ask(2*x, Q.even) == None
    assert ask(2*x, Q.even, Assume(x, Q.integer)) == True
    assert ask(2*x, Q.even, Assume(x, Q.even)) == True
    assert ask(2*x, Q.even, Assume(x, Q.irrational)) == False
    assert ask(2*x, Q.even, Assume(x, Q.odd)) == True
    assert ask(2*x, Q.even, Assume(x, Q.integer, False)) == None
    assert ask(3*x, Q.even, Assume(x, Q.integer)) == None
    assert ask(3*x, Q.even, Assume(x, Q.even)) == True
    assert ask(3*x, Q.even, Assume(x, Q.odd)) == False

    assert ask(x+1, Q.even, Assume(x, Q.odd)) == True
    assert ask(x+1, Q.even, Assume(x, Q.even)) == False
    assert ask(x+2, Q.even, Assume(x, Q.odd)) == False
    assert ask(x+2, Q.even, Assume(x, Q.even)) == True
    assert ask(7-x, Q.even, Assume(x, Q.odd)) == True
    assert ask(7+x, Q.even, Assume(x, Q.odd)) == True
    assert ask(x+y, Q.even, Assume(x, Q.odd) & Assume(y, Q.odd)) == True
    assert ask(x+y, Q.even, Assume(x, Q.odd) & Assume(y, Q.even)) == False
    assert ask(x+y, Q.even, Assume(x, Q.even) & Assume(y, Q.even)) == True

    assert ask(2*x + 1, Q.even, Assume(x, Q.integer)) == False
    assert ask(2*x*y, Q.even, Assume(x, Q.rational) & Assume(x, Q.rational)) == None
    assert ask(2*x*y, Q.even, Assume(x, Q.irrational) & Assume(x, Q.irrational)) == None

    assert ask(x+y+z, Q.even, Assume(x, Q.odd) & Assume(y, Q.odd) & \
                     Assume(z, Q.even)) == True
    assert ask(x+y+z+t, Q.even, Assume(x, Q.odd) & Assume(y, Q.odd) & \
                     Assume(z, Q.even) & Assume(t, Q.integer)) == None

    assert ask(abs(x), Q.even, Assume(x, Q.even)) == True
    assert ask(abs(x), Q.even, Assume(x, Q.even, False)) == None
    assert ask(re(x),  Q.even, Assume(x, Q.even)) == True
    assert ask(re(x),  Q.even, Assume(x, Q.even, False)) == None
    assert ask(im(x),  Q.even, Assume(x, Q.even)) == True
    assert ask(im(x),  Q.even, Assume(x, Q.real)) == True

def test_extended_real():
    x = symbols('x')
    assert ask(x, Q.extended_real, Assume(x, Q.positive)) == True
    assert ask(-x, Q.extended_real, Assume(x, Q.positive)) == True
    assert ask(-x, Q.extended_real, Assume(x, Q.negative)) == True

    assert ask(x+S.Infinity, Q.extended_real, Assume(x, Q.real)) == True

def test_rational():
    x, y = symbols('xy')
    assert ask(x, Q.rational, Assume(x, Q.integer)) == True
    assert ask(x, Q.rational, Assume(x, Q.irrational)) == False
    assert ask(x, Q.rational, Assume(x, Q.real)) == None
    assert ask(x, Q.rational, Assume(x, Q.positive)) == None
    assert ask(x, Q.rational, Assume(x, Q.negative)) == None
    assert ask(x, Q.rational, Assume(x, Q.nonzero)) == None

    assert ask(2*x, Q.rational, Assume(x, Q.rational)) == True
    assert ask(2*x, Q.rational, Assume(x, Q.integer)) == True
    assert ask(2*x, Q.rational, Assume(x, Q.even)) == True
    assert ask(2*x, Q.rational, Assume(x, Q.odd)) == True
    assert ask(2*x, Q.rational, Assume(x, Q.irrational)) == False

    assert ask(x/2, Q.rational, Assume(x, Q.rational)) == True
    assert ask(x/2, Q.rational, Assume(x, Q.integer)) == True
    assert ask(x/2, Q.rational, Assume(x, Q.even)) == True
    assert ask(x/2, Q.rational, Assume(x, Q.odd)) == True
    assert ask(x/2, Q.rational, Assume(x, Q.irrational)) == False

    assert ask(1/x, Q.rational, Assume(x, Q.rational)) == True
    assert ask(1/x, Q.rational, Assume(x, Q.integer)) == True
    assert ask(1/x, Q.rational, Assume(x, Q.even)) == True
    assert ask(1/x, Q.rational, Assume(x, Q.odd)) == True
    assert ask(1/x, Q.rational, Assume(x, Q.irrational)) == False

    assert ask(2/x, Q.rational, Assume(x, Q.rational)) == True
    assert ask(2/x, Q.rational, Assume(x, Q.integer)) == True
    assert ask(2/x, Q.rational, Assume(x, Q.even)) == True
    assert ask(2/x, Q.rational, Assume(x, Q.odd)) == True
    assert ask(2/x, Q.rational, Assume(x, Q.irrational)) == False

    # with multiple symbols
    assert ask(x*y, Q.rational, Assume(x, Q.irrational) & \
        Assume(y, Q.irrational)) == None
    assert ask(y/x, Q.rational, Assume(x, Q.rational) & \
        Assume(y, Q.rational)) == True
    assert ask(y/x, Q.rational, Assume(x, Q.integer) & \
        Assume(y, Q.rational)) == True
    assert ask(y/x, Q.rational, Assume(x, Q.even) & \
        Assume(y, Q.rational)) == True
    assert ask(y/x, Q.rational, Assume(x, Q.odd) & \
        Assume(y, Q.rational)) == True
    assert ask(y/x, Q.rational, Assume(x, Q.irrational) & \
        Assume(y, Q.rational)) == False

def test_imaginary():
    x, y, z = symbols('x y z')
    I = S.ImaginaryUnit
    assert ask(x, Q.imaginary) == None
    assert ask(x, Q.imaginary, Assume(x, Q.real)) == False
    assert ask(x, Q.imaginary, Assume(x, Q.prime)) == False

    assert ask(x+1, Q.imaginary, Assume(x, Q.real)) == False
    assert ask(x+1, Q.imaginary, Assume(x, Q.imaginary)) == False
    assert ask(x+I, Q.imaginary, Assume(x, Q.real)) == False
    assert ask(x+I, Q.imaginary, Assume(x, Q.imaginary)) == True
    assert ask(x+y, Q.imaginary, Assume(x, Q.imaginary) & \
                     Assume(y, Q.imaginary)) == True
    assert ask(x+y, Q.imaginary, Assume(x, Q.real) & \
                     Assume(y, Q.real)) == False
    assert ask(x+y, Q.imaginary, Assume(x, Q.imaginary) & \
                     Assume(y, Q.real)) == False
    assert ask(x+y, Q.imaginary, Assume(x, Q.complex) & \
                     Assume(y, Q.real)) == None

    assert ask(I*x, Q.imaginary, Assume(x, Q.real)) == True
    assert ask(I*x, Q.imaginary, Assume(x, Q.imaginary)) == False
    assert ask(I*x, Q.imaginary, Assume(x, Q.complex)) == None
    assert ask(x*y, Q.imaginary, Assume(x, Q.imaginary) & \
                 Assume(y, Q.real)) == True

    assert ask(x+y+z, Q.imaginary, Assume(x, Q.real) & \
                     Assume(y, Q.real) & Assume(z, Q.real)) == False
    assert ask(x+y+z, Q.imaginary, Assume(x, Q.real) & \
                     Assume(y, Q.real) & Assume(z, Q.imaginary)) == None
    assert ask(x+y+z, Q.imaginary, Assume(x, Q.real) & \
                     Assume(y, Q.imaginary) & Assume(z, Q.imaginary)) == False

def test_infinitesimal():
    x, y = symbols('x y')
    assert ask(x, Q.infinitesimal) == None
    assert ask(x, Q.infinitesimal, Assume(x, Q.infinitesimal)) == True

    assert ask(2*x, Q.infinitesimal, Assume(x, Q.infinitesimal)) == True
    assert ask(x*y, Q.infinitesimal, Assume(x, Q.infinitesimal)) == None
    assert ask(x*y, Q.infinitesimal, Assume(x, Q.infinitesimal) & \
                     Assume(y, Q.infinitesimal)) == True
    assert ask(x*y, Q.infinitesimal, Assume(x, Q.infinitesimal) & \
                     Assume(y, Q.bounded)) == True

    assert ask(x**2, Q.infinitesimal, Assume(x, Q.infinitesimal)) == True

def test_integer():
    x = symbols('x')
    assert ask(x, Q.integer) == None
    assert ask(x, Q.integer, Assume(x, Q.integer)) == True
    assert ask(x, Q.integer, Assume(x, Q.integer, False)) == False
    assert ask(x, Q.integer, Assume(x, Q.real, False)) == False
    assert ask(x, Q.integer, Assume(x, Q.positive, False)) == None
    assert ask(x, Q.integer, Assume(x, Q.even) | Assume(x, Q.odd)) == True

    assert ask(2*x, Q.integer, Assume(x, Q.integer)) == True
    assert ask(2*x, Q.integer, Assume(x, Q.even)) == True
    assert ask(2*x, Q.integer, Assume(x, Q.prime)) == True
    assert ask(2*x, Q.integer, Assume(x, Q.rational)) == None
    assert ask(2*x, Q.integer, Assume(x, Q.real)) == None
    assert ask(sqrt(2)*x, Q.integer, Assume(x, Q.integer)) == False

    assert ask(x/2, Q.integer, Assume(x, Q.odd)) == False
    assert ask(x/2, Q.integer, Assume(x, Q.even)) == True
    assert ask(x/3, Q.integer, Assume(x, Q.odd)) == None
    assert ask(x/3, Q.integer, Assume(x, Q.even)) == None

def test_negative():
    x, y = symbols('xy')
    assert ask(x, Q.negative, Assume(x, Q.negative)) == True
    assert ask(x, Q.negative, Assume(x, Q.positive)) == False
    assert ask(x, Q.negative, Assume(x, Q.real, False)) == False
    assert ask(x, Q.negative, Assume(x, Q.prime)) == False
    assert ask(x, Q.negative, Assume(x, Q.prime, False)) == None

    assert ask(-x, Q.negative, Assume(x, Q.positive)) == True
    assert ask(-x, Q.negative, Assume(x, Q.positive, False)) == None
    assert ask(-x, Q.negative, Assume(x, Q.negative)) == False
    assert ask(-x, Q.negative, Assume(x, Q.positive)) == True

    assert ask(x-1, Q.negative, Assume(x, Q.negative)) == True
    assert ask(x+y, Q.negative) == None
    assert ask(x+y, Q.negative, Assume(x, Q.negative)) == None
    assert ask(x+y, Q.negative, Assume(x, Q.negative) &\
                     Assume(y, Q.negative)) == True

    assert ask(x**2, Q.negative) == None
    assert ask(x**2, Q.negative, Assume(x, Q.real)) == False
    assert ask(x**1.4, Q.negative, Assume(x, Q.real)) == None

    assert ask(x*y, Q.negative) == None
    assert ask(x*y, Q.negative, Assume(x, Q.positive) & \
                     Assume(y, Q.positive)) == False
    assert ask(x*y, Q.negative, Assume(x, Q.positive) & \
                     Assume(y, Q.negative)) == True
    assert ask(x*y, Q.negative, Assume(x, Q.complex) & \
                     Assume(y, Q.complex)) == None

    assert ask(x**y, Q.negative) == None
    assert ask(x**y, Q.negative, Assume(x, Q.negative) & \
                     Assume(y, Q.even)) == False
    assert ask(x**y, Q.negative, Assume(x, Q.negative) & \
                     Assume(y, Q.odd)) == True
    assert ask(x**y, Q.negative, Assume(x, Q.positive) & \
                     Assume(y, Q.integer)) == False

    assert ask(abs(x), Q.negative) == False

def test_nonzero():
    x, y = symbols('xy')
    assert ask(x, Q.nonzero) == None
    assert ask(x, Q.nonzero, Assume(x, Q.real)) == None
    assert ask(x, Q.nonzero, Assume(x, Q.positive)) == True
    assert ask(x, Q.nonzero, Assume(x, Q.negative)) == True
    assert ask(x, Q.nonzero, Assume(x, Q.negative) | Assume(x, Q.positive)) == True

    assert ask(x+y, Q.nonzero) == None
    assert ask(x+y, Q.nonzero, Assume(x, Q.positive) & Assume(y, Q.positive)) == True
    assert ask(x+y, Q.nonzero, Assume(x, Q.positive) & Assume(y, Q.negative)) == None
    assert ask(x+y, Q.nonzero, Assume(x, Q.negative) & Assume(y, Q.negative)) == True

    assert ask(2*x, Q.nonzero) == None
    assert ask(2*x, Q.nonzero, Assume(x, Q.positive)) == True
    assert ask(2*x, Q.nonzero, Assume(x, Q.negative)) == True
    assert ask(x*y, Q.nonzero, Assume(x, Q.nonzero)) == None
    assert ask(x*y, Q.nonzero, Assume(x, Q.nonzero) & Assume(y, Q.nonzero)) == True

    assert ask(abs(x), Q.nonzero) == None
    assert ask(abs(x), Q.nonzero, Assume(x, Q.nonzero)) == True

def test_odd():
    x, y, z, t = symbols('x y z t')
    assert ask(x, Q.odd) == None
    assert ask(x, Q.odd, Assume(x, Q.odd)) == True
    assert ask(x, Q.odd, Assume(x, Q.integer)) == None
    assert ask(x, Q.odd, Assume(x, Q.integer, False)) == False
    assert ask(x, Q.odd, Assume(x, Q.rational)) == None
    assert ask(x, Q.odd, Assume(x, Q.positive)) == None

    assert ask(-x, Q.odd, Assume(x, Q.odd)) == True

    assert ask(2*x, Q.odd) == None
    assert ask(2*x, Q.odd, Assume(x, Q.integer)) == False
    assert ask(2*x, Q.odd, Assume(x, Q.odd)) == False
    assert ask(2*x, Q.odd, Assume(x, Q.irrational)) == False
    assert ask(2*x, Q.odd, Assume(x, Q.integer, False)) == None
    assert ask(3*x, Q.odd, Assume(x, Q.integer)) == None

    assert ask(x/3, Q.odd, Assume(x, Q.odd)) == None
    assert ask(x/3, Q.odd, Assume(x, Q.even)) == None

    assert ask(x+1, Q.odd, Assume(x, Q.even)) == True
    assert ask(x+2, Q.odd, Assume(x, Q.even)) == False
    assert ask(x+2, Q.odd, Assume(x, Q.odd))  == True
    assert ask(3-x, Q.odd, Assume(x, Q.odd))  == False
    assert ask(3-x, Q.odd, Assume(x, Q.even))  == True
    assert ask(3+x, Q.odd, Assume(x, Q.odd))  == False
    assert ask(3+x, Q.odd, Assume(x, Q.even))  == True
    assert ask(x+y, Q.odd, Assume(x, Q.odd) & Assume(y, Q.odd)) == False
    assert ask(x+y, Q.odd, Assume(x, Q.odd) & Assume(y, Q.even)) == True
    assert ask(x-y, Q.odd, Assume(x, Q.even) & Assume(y, Q.odd)) == True
    assert ask(x-y, Q.odd, Assume(x, Q.odd) & Assume(y, Q.odd)) == False

    assert ask(x+y+z, Q.odd, Assume(x, Q.odd) & Assume(y, Q.odd) & \
                     Assume(z, Q.even)) == False
    assert ask(x+y+z+t, Q.odd, Assume(x, Q.odd) & Assume(y, Q.odd) & \
                     Assume(z, Q.even) & Assume(t, Q.integer)) == None

    assert ask(2*x + 1, Q.odd, Assume(x, Q.integer)) == True
    assert ask(2*x + y, Q.odd, Assume(x, Q.integer) & Assume(y, Q.odd)) == True
    assert ask(2*x + y, Q.odd, Assume(x, Q.integer) & Assume(y, Q.even)) == False
    assert ask(2*x + y, Q.odd, Assume(x, Q.integer) & Assume(y, Q.integer)) == None
    assert ask(x*y,   Q.odd, Assume(x, Q.odd) & Assume(y, Q.even)) == False
    assert ask(x*y,   Q.odd, Assume(x, Q.odd) & Assume(y, Q.odd)) == True
    assert ask(2*x*y, Q.odd, Assume(x, Q.rational) & Assume(x, Q.rational)) == None
    assert ask(2*x*y, Q.odd, Assume(x, Q.irrational) & Assume(x, Q.irrational)) == None

    assert ask(abs(x), Q.odd, Assume(x, Q.odd)) == True

def test_prime():
    x, y = symbols('x y')
    assert ask(x, Q.prime, Assume(x, Q.prime)) == True
    assert ask(x, Q.prime, Assume(x, Q.prime, False)) == False
    assert ask(x, Q.prime, Assume(x, Q.integer)) == None
    assert ask(x, Q.prime, Assume(x, Q.integer, False)) == False

    assert ask(2*x, Q.prime, Assume(x, Q.integer)) == False
    assert ask(x*y, Q.prime) == None
    assert ask(x*y, Q.prime, Assume(x, Q.prime)) == None
    assert ask(x*y, Q.prime, Assume(x, Q.integer) & \
                     Assume(y, Q.integer)) == False

    assert ask(x**2, Q.prime, Assume(x, Q.integer)) == False
    assert ask(x**2, Q.prime, Assume(x, Q.prime)) == False
    assert ask(x**y, Q.prime, Assume(x, Q.integer) & \
                     Assume(y, Q.integer)) == False

def test_positive():
    x, y, z, w = symbols('xyzw')
    assert ask(x, Q.positive, Assume(x, Q.positive)) == True
    assert ask(x, Q.positive, Assume(x, Q.negative)) == False
    assert ask(x, Q.positive, Assume(x, Q.nonzero)) == None

    assert ask(-x, Q.positive, Assume(x, Q.positive)) == False
    assert ask(-x, Q.positive, Assume(x, Q.negative)) == True

    assert ask(x+y, Q.positive, Assume(x, Q.positive) & \
                     Assume(y, Q.positive)) == True
    assert ask(x+y, Q.positive, Assume(x, Q.positive) & \
                     Assume(y, Q.negative)) == None

    assert ask(2*x,  Q.positive, Assume(x, Q.positive)) == True
    assumptions =  Assume(x, Q.positive) & Assume(y, Q.negative) & \
                    Assume(z, Q.negative) & Assume(w, Q.positive)
    assert ask(x*y*z,  Q.positive)  == None
    assert ask(x*y*z,  Q.positive, assumptions) == True
    assert ask(-x*y*z, Q.positive, assumptions) == False

    assert ask(x**2, Q.positive, Assume(x, Q.positive)) == True
    assert ask(x**2, Q.positive, Assume(x, Q.negative)) == True

    #exponential
    assert ask(exp(x),     Q.positive, Assume(x, Q.real)) == True
    assert ask(x + exp(x), Q.positive, Assume(x, Q.real)) == None

    #absolute value
    assert ask(abs(x), Q.positive) == None # abs(0) = 0
    assert ask(abs(x), Q.positive, Assume(x, Q.positive)) == True

@XFAIL
def test_positive_xfail():
    assert ask(1/(1 + x**2), Q.positive, Assume(x, Q.real)) == True

def test_real():
    x, y = symbols('x y')
    assert ask(x, Q.real) == None
    assert ask(x, Q.real, Assume(x, Q.real)) == True
    assert ask(x, Q.real, Assume(x, Q.nonzero)) == True
    assert ask(x, Q.real, Assume(x, Q.positive)) == True
    assert ask(x, Q.real, Assume(x, Q.negative)) == True
    assert ask(x, Q.real, Assume(x, Q.integer)) == True
    assert ask(x, Q.real, Assume(x, Q.even)) == True
    assert ask(x, Q.real, Assume(x, Q.prime)) == True

    assert ask(x/sqrt(2), Q.real, Assume(x, Q.real)) == True
    assert ask(x/sqrt(-2), Q.real, Assume(x, Q.real)) == False

    I = S.ImaginaryUnit
    assert ask(x+1, Q.real, Assume(x, Q.real)) == True
    assert ask(x+I, Q.real, Assume(x, Q.real)) == False
    assert ask(x+I, Q.real, Assume(x, Q.complex)) == None

    assert ask(2*x, Q.real, Assume(x, Q.real)) == True
    assert ask(I*x, Q.real, Assume(x, Q.real)) == False
    assert ask(I*x, Q.real, Assume(x, Q.imaginary)) == True
    assert ask(I*x, Q.real, Assume(x, Q.complex)) == None

    assert ask(x**2, Q.real, Assume(x, Q.real)) == True
    assert ask(sqrt(x), Q.real, Assume(x, Q.negative)) == False
    assert ask(x**y, Q.real, Assume(x, Q.real) & Assume(y, Q.integer)) == True
    assert ask(x**y, Q.real, Assume(x, Q.real) & Assume(y, Q.real)) == None
    assert ask(x**y, Q.real, Assume(x, Q.positive) & \
                     Assume(y, Q.real)) == True

    # trigonometric functions
    assert ask(sin(x), Q.real) == None
    assert ask(cos(x), Q.real) == None
    assert ask(sin(x), Q.real, Assume(x, Q.real)) == True
    assert ask(cos(x), Q.real, Assume(x, Q.real)) == True

    # exponential function
    assert ask(exp(x), Q.real) == None
    assert ask(exp(x), Q.real, Assume(x, Q.real)) == True
    assert ask(x + exp(x), Q.real, Assume(x, Q.real)) == True

    # Q.complexes
    assert ask(re(x), Q.real) == True
    assert ask(im(x), Q.real) == True



def test_global():
    """Test ask with global assumptions"""
    x = symbols('x')
    assert ask(x, Q.integer) == None
    global_assumptions.add(Assume(x, Q.integer))
    assert ask(x, Q.integer) == True
    global_assumptions.clear()
    assert ask(x, Q.integer) == None

def test_incompatible_resolutors():
    x = symbols('x')
    class Prime2AskHandler(AskHandler):
        @staticmethod
        def Number(expr, assumptions):
            return True
    register_handler('prime', Prime2AskHandler)
    raises(ValueError, 'ask(4, Q.prime)')

def test_key_extensibility():
    """test that you can add keys to the ask system at runtime"""
    x = Symbol('x')
    # make sure thie key is not defined
    raises(KeyError, "ask(x, 'my_key')")
    class MyAskHandler(AskHandler):
        @staticmethod
        def Symbol(expr, assumptions):
            return True
    register_handler('my_key', MyAskHandler)
    assert ask(x, 'my_key') == True
    assert ask(x+1, 'my_key') == None
    remove_handler('my_key', MyAskHandler)

def test_type_extensibility():
    """test that new types can be added to the ask system at runtime
    We create a custom type MyType, and override ask Q.prime=True with handler
    MyAskHandler for this type

    TODO: test incompatible resolutors
    """
    from sympy.core import Basic

    class MyType(Basic):
        pass

    class MyAskHandler(AskHandler):
        @staticmethod
        def MyType(expr, assumptions):
            return True

    a = MyType()
    register_handler(Q.prime, MyAskHandler)
    assert ask(a, Q.prime) == True
