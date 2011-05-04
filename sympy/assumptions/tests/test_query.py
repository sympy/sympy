from sympy.utilities.pytest import raises, XFAIL

from sympy.core import Symbol, symbols, S, Rational, Integer, I, pi, oo
from sympy.functions import exp, log, sin, cos, sign, re, im, sqrt, Abs
from sympy.assumptions import (global_assumptions, Q, ask,
    register_handler, remove_handler, AssumptionsContext)
from sympy.assumptions.handlers import AskHandler
from sympy.assumptions.ask import (compute_known_facts,
                                   known_facts_cnf, known_facts_dict)

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
    x, y = symbols('x,y')
    assert ask(x, Q.bounded) == False
    assert ask(x, Q.bounded, Q.bounded(x)) == True
    assert ask(x, Q.bounded, Q.bounded(y)) == False
    assert ask(x, Q.bounded, Q.complex(x)) == False

    assert ask(x+1, Q.bounded) == False
    assert ask(x+1, Q.bounded, Q.bounded(x)) == True
    assert ask(x+y, Q.bounded) == None
    assert ask(x+y, Q.bounded, Q.bounded(x)) == False
    assert ask(x+1, Q.bounded, Q.bounded(x) & Q.bounded(y)) == True

    assert ask(2*x, Q.bounded) == False
    assert ask(2*x, Q.bounded, Q.bounded(x)) == True
    assert ask(x*y, Q.bounded) == None
    assert ask(x*y, Q.bounded, Q.bounded(x)) == False
    assert ask(x*y, Q.bounded, Q.bounded(x) & Q.bounded(y)) == True

    assert ask(x**2, Q.bounded) == False
    assert ask(2**x, Q.bounded) == False
    assert ask(2**x, Q.bounded, Q.bounded(x)) == True
    assert ask(x**x, Q.bounded) == False
    assert ask(Rational(1,2) ** x, Q.bounded) == True
    assert ask(x ** Rational(1,2), Q.bounded) == False

    # sign function
    assert ask(sign(x), Q.bounded) == True
    assert ask(sign(x), Q.bounded, ~Q.bounded(x)) == True

    # exponential functions
    assert ask(log(x), Q.bounded) == False
    assert ask(log(x), Q.bounded, Q.bounded(x)) == True
    assert ask(exp(x), Q.bounded) == False
    assert ask(exp(x), Q.bounded, Q.bounded(x)) == True
    assert ask(exp(2), Q.bounded) == True

    # trigonometric functions
    assert ask(sin(x), Q.bounded) == True
    assert ask(sin(x), Q.bounded, ~Q.bounded(x)) == True
    assert ask(cos(x), Q.bounded) == True
    assert ask(cos(x), Q.bounded, ~Q.bounded(x)) == True
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
    x, y = symbols('x,y')
    assert ask(x, Q.commutative) == True
    assert ask(x, Q.commutative, ~Q.commutative(x)) == False
    assert ask(x, Q.commutative, Q.complex(x)) == True
    assert ask(x, Q.commutative, Q.imaginary(x)) == True
    assert ask(x, Q.commutative, Q.real(x)) == True
    assert ask(x, Q.commutative, Q.positive(x)) == True
    assert ask(x, Q.commutative, ~Q.commutative(y))  == True

    assert ask(2*x, Q.commutative) == True
    assert ask(2*x, Q.commutative, ~Q.commutative(x)) == False

    assert ask(x + 1, Q.commutative) == True
    assert ask(x + 1, Q.commutative, ~Q.commutative(x)) == False

    assert ask(x**2, Q.commutative) == True
    assert ask(x**2, Q.commutative, ~Q.commutative(x)) == False

    assert ask(log(x), Q.commutative) == True

def test_complex():
    x, y = symbols('x,y')
    assert ask(x, Q.complex) == None
    assert ask(x, Q.complex, Q.complex(x)) == True
    assert ask(x, Q.complex, Q.complex(y)) == None
    assert ask(x, Q.complex, ~Q.complex(x)) == False
    assert ask(x, Q.complex, Q.real(x)) == True
    assert ask(x, Q.complex, ~Q.real(x)) == None
    assert ask(x, Q.complex, Q.rational(x)) == True
    assert ask(x, Q.complex, Q.irrational(x)) == True
    assert ask(x, Q.complex, Q.positive(x)) == True
    assert ask(x, Q.complex, Q.imaginary(x)) == True

    # a+b
    assert ask(x+1, Q.complex, Q.complex(x)) == True
    assert ask(x+1, Q.complex, Q.real(x)) == True
    assert ask(x+1, Q.complex, Q.rational(x)) == True
    assert ask(x+1, Q.complex, Q.irrational(x)) == True
    assert ask(x+1, Q.complex, Q.imaginary(x)) == True
    assert ask(x+1, Q.complex, Q.integer(x))  == True
    assert ask(x+1, Q.complex, Q.even(x))  == True
    assert ask(x+1, Q.complex, Q.odd(x))  == True
    assert ask(x+y, Q.complex, Q.complex(x) & Q.complex(y)) == True
    assert ask(x+y, Q.complex, Q.real(x) & Q.imaginary(y)) == True

    # a*x +b
    assert ask(2*x+1, Q.complex, Q.complex(x)) == True
    assert ask(2*x+1, Q.complex, Q.real(x)) == True
    assert ask(2*x+1, Q.complex, Q.positive(x)) == True
    assert ask(2*x+1, Q.complex, Q.rational(x)) == True
    assert ask(2*x+1, Q.complex, Q.irrational(x)) == True
    assert ask(2*x+1, Q.complex, Q.imaginary(x)) == True
    assert ask(2*x+1, Q.complex, Q.integer(x))  == True
    assert ask(2*x+1, Q.complex, Q.even(x))  == True
    assert ask(2*x+1, Q.complex, Q.odd(x))  == True

    # x**2
    assert ask(x**2, Q.complex, Q.complex(x)) == True
    assert ask(x**2, Q.complex, Q.real(x)) == True
    assert ask(x**2, Q.complex, Q.positive(x)) == True
    assert ask(x**2, Q.complex, Q.rational(x)) == True
    assert ask(x**2, Q.complex, Q.irrational(x)) == True
    assert ask(x**2, Q.complex, Q.imaginary(x)) == True
    assert ask(x**2, Q.complex, Q.integer(x))  == True
    assert ask(x**2, Q.complex, Q.even(x))  == True
    assert ask(x**2, Q.complex, Q.odd(x))  == True

    # 2**x
    assert ask(2**x, Q.complex, Q.complex(x)) == True
    assert ask(2**x, Q.complex, Q.real(x)) == True
    assert ask(2**x, Q.complex, Q.positive(x)) == True
    assert ask(2**x, Q.complex, Q.rational(x)) == True
    assert ask(2**x, Q.complex, Q.irrational(x)) == True
    assert ask(2**x, Q.complex, Q.imaginary(x)) == True
    assert ask(2**x, Q.complex, Q.integer(x))  == True
    assert ask(2**x, Q.complex, Q.even(x))  == True
    assert ask(2**x, Q.complex, Q.odd(x))  == True
    assert ask(x**y, Q.complex, Q.complex(x) & Q.complex(y)) == True

    # trigonometric expressions
    assert ask(sin(x), Q.complex) == True
    assert ask(sin(2*x + 1), Q.complex) == True
    assert ask(cos(x), Q.complex) == True
    assert ask(cos(2*x+1), Q.complex) == True

    # exponential
    assert ask(exp(x), Q.complex) == True
    assert ask(exp(x), Q.complex) == True

    # Q.complexes
    assert ask(Abs(x), Q.complex) == True
    assert ask(re(x),  Q.complex) == True
    assert ask(im(x),  Q.complex) == True

def test_even():
    x, y, z, t = symbols('x,y,z,t')
    assert ask(x, Q.even) == None
    assert ask(x, Q.even, Q.integer(x)) == None
    assert ask(x, Q.even, ~Q.integer(x)) == False
    assert ask(x, Q.even, Q.rational(x)) == None
    assert ask(x, Q.even, Q.positive(x)) == None

    assert ask(2*x, Q.even) == None
    assert ask(2*x, Q.even, Q.integer(x)) == True
    assert ask(2*x, Q.even, Q.even(x)) == True
    assert ask(2*x, Q.even, Q.irrational(x)) == False
    assert ask(2*x, Q.even, Q.odd(x)) == True
    assert ask(2*x, Q.even, ~Q.integer(x)) == None
    assert ask(3*x, Q.even, Q.integer(x)) == None
    assert ask(3*x, Q.even, Q.even(x)) == True
    assert ask(3*x, Q.even, Q.odd(x)) == False

    assert ask(x+1, Q.even, Q.odd(x)) == True
    assert ask(x+1, Q.even, Q.even(x)) == False
    assert ask(x+2, Q.even, Q.odd(x)) == False
    assert ask(x+2, Q.even, Q.even(x)) == True
    assert ask(7-x, Q.even, Q.odd(x)) == True
    assert ask(7+x, Q.even, Q.odd(x)) == True
    assert ask(x+y, Q.even, Q.odd(x) & Q.odd(y)) == True
    assert ask(x+y, Q.even, Q.odd(x) & Q.even(y)) == False
    assert ask(x+y, Q.even, Q.even(x) & Q.even(y)) == True

    assert ask(2*x + 1, Q.even, Q.integer(x)) == False
    assert ask(2*x*y, Q.even, Q.rational(x) & Q.rational(x)) == None
    assert ask(2*x*y, Q.even, Q.irrational(x) & Q.irrational(x)) == None

    assert ask(x+y+z, Q.even, Q.odd(x) & Q.odd(y) & Q.even(z)) == True
    assert ask(x+y+z+t, Q.even,
               Q.odd(x) & Q.odd(y) & Q.even(z) & Q.integer(t)) == None

    assert ask(Abs(x), Q.even, Q.even(x)) == True
    assert ask(Abs(x), Q.even, ~Q.even(x)) == None
    assert ask(re(x),  Q.even, Q.even(x)) == True
    assert ask(re(x),  Q.even, ~Q.even(x)) == None
    assert ask(im(x),  Q.even, Q.even(x)) == True
    assert ask(im(x),  Q.even, Q.real(x)) == True

def test_extended_real():
    x = symbols('x')
    assert ask(x, Q.extended_real, Q.positive(x)) == True
    assert ask(-x, Q.extended_real, Q.positive(x)) == True
    assert ask(-x, Q.extended_real, Q.negative(x)) == True

    assert ask(x+S.Infinity, Q.extended_real, Q.real(x)) == True

def test_rational():
    x, y = symbols('x,y')
    assert ask(x, Q.rational, Q.integer(x)) == True
    assert ask(x, Q.rational, Q.irrational(x)) == False
    assert ask(x, Q.rational, Q.real(x)) == None
    assert ask(x, Q.rational, Q.positive(x)) == None
    assert ask(x, Q.rational, Q.negative(x)) == None
    assert ask(x, Q.rational, Q.nonzero(x)) == None

    assert ask(2*x, Q.rational, Q.rational(x)) == True
    assert ask(2*x, Q.rational, Q.integer(x)) == True
    assert ask(2*x, Q.rational, Q.even(x)) == True
    assert ask(2*x, Q.rational, Q.odd(x)) == True
    assert ask(2*x, Q.rational, Q.irrational(x)) == False

    assert ask(x/2, Q.rational, Q.rational(x)) == True
    assert ask(x/2, Q.rational, Q.integer(x)) == True
    assert ask(x/2, Q.rational, Q.even(x)) == True
    assert ask(x/2, Q.rational, Q.odd(x)) == True
    assert ask(x/2, Q.rational, Q.irrational(x)) == False

    assert ask(1/x, Q.rational, Q.rational(x)) == True
    assert ask(1/x, Q.rational, Q.integer(x)) == True
    assert ask(1/x, Q.rational, Q.even(x)) == True
    assert ask(1/x, Q.rational, Q.odd(x)) == True
    assert ask(1/x, Q.rational, Q.irrational(x)) == False

    assert ask(2/x, Q.rational, Q.rational(x)) == True
    assert ask(2/x, Q.rational, Q.integer(x)) == True
    assert ask(2/x, Q.rational, Q.even(x)) == True
    assert ask(2/x, Q.rational, Q.odd(x)) == True
    assert ask(2/x, Q.rational, Q.irrational(x)) == False

    # with multiple symbols
    assert ask(x*y, Q.rational, Q.irrational(x) & Q.irrational(y)) == None
    assert ask(y/x, Q.rational, Q.rational(x) & Q.rational(y)) == True
    assert ask(y/x, Q.rational, Q.integer(x) & Q.rational(y)) == True
    assert ask(y/x, Q.rational, Q.even(x) & Q.rational(y)) == True
    assert ask(y/x, Q.rational, Q.odd(x) & Q.rational(y)) == True
    assert ask(y/x, Q.rational, Q.irrational(x) & Q.rational(y)) == False

def test_imaginary():
    x, y, z = symbols('x,y,z')
    I = S.ImaginaryUnit
    assert ask(x, Q.imaginary) == None
    assert ask(x, Q.imaginary, Q.real(x)) == False
    assert ask(x, Q.imaginary, Q.prime(x)) == False

    assert ask(x+1, Q.imaginary, Q.real(x)) == False
    assert ask(x+1, Q.imaginary, Q.imaginary(x)) == False
    assert ask(x+I, Q.imaginary, Q.real(x)) == False
    assert ask(x+I, Q.imaginary, Q.imaginary(x)) == True
    assert ask(x+y, Q.imaginary, Q.imaginary(x) & Q.imaginary(y)) == True
    assert ask(x+y, Q.imaginary, Q.real(x) & Q.real(y)) == False
    assert ask(x+y, Q.imaginary, Q.imaginary(x) & Q.real(y)) == False
    assert ask(x+y, Q.imaginary, Q.complex(x) & Q.real(y)) == None

    assert ask(I*x, Q.imaginary, Q.real(x)) == True
    assert ask(I*x, Q.imaginary, Q.imaginary(x)) == False
    assert ask(I*x, Q.imaginary, Q.complex(x)) == None
    assert ask(x*y, Q.imaginary, Q.imaginary(x) & Q.real(y)) == True

    assert ask(x+y+z, Q.imaginary, Q.real(x) & Q.real(y) & Q.real(z)) == False
    assert ask(x+y+z, Q.imaginary, Q.real(x) & Q.real(y) & Q.imaginary(z)) == None
    assert ask(x+y+z, Q.imaginary, Q.real(x) & Q.imaginary(y) & Q.imaginary(z)) == False

def test_infinitesimal():
    x, y = symbols('x,y')
    assert ask(x, Q.infinitesimal) == None
    assert ask(x, Q.infinitesimal, Q.infinitesimal(x)) == True

    assert ask(2*x, Q.infinitesimal, Q.infinitesimal(x)) == True
    assert ask(x*y, Q.infinitesimal, Q.infinitesimal(x)) == None
    assert ask(x*y, Q.infinitesimal, Q.infinitesimal(x) & Q.infinitesimal(y)) == True
    assert ask(x*y, Q.infinitesimal, Q.infinitesimal(x) & Q.bounded(y)) == True

    assert ask(x**2, Q.infinitesimal, Q.infinitesimal(x)) == True

def test_integer():
    x = symbols('x')
    assert ask(x, Q.integer) == None
    assert ask(x, Q.integer, Q.integer(x)) == True
    assert ask(x, Q.integer, ~Q.integer(x)) == False
    assert ask(x, Q.integer, ~Q.real(x)) == False
    assert ask(x, Q.integer, ~Q.positive(x)) == None
    assert ask(x, Q.integer, Q.even(x) | Q.odd(x)) == True

    assert ask(2*x, Q.integer, Q.integer(x)) == True
    assert ask(2*x, Q.integer, Q.even(x)) == True
    assert ask(2*x, Q.integer, Q.prime(x)) == True
    assert ask(2*x, Q.integer, Q.rational(x)) == None
    assert ask(2*x, Q.integer, Q.real(x)) == None
    assert ask(sqrt(2)*x, Q.integer, Q.integer(x)) == False

    assert ask(x/2, Q.integer, Q.odd(x)) == False
    assert ask(x/2, Q.integer, Q.even(x)) == True
    assert ask(x/3, Q.integer, Q.odd(x)) == None
    assert ask(x/3, Q.integer, Q.even(x)) == None

def test_negative():
    x, y = symbols('x,y')
    assert ask(x, Q.negative, Q.negative(x)) == True
    assert ask(x, Q.negative, Q.positive(x)) == False
    assert ask(x, Q.negative, ~Q.real(x)) == False
    assert ask(x, Q.negative, Q.prime(x)) == False
    assert ask(x, Q.negative, ~Q.prime(x)) == None

    assert ask(-x, Q.negative, Q.positive(x)) == True
    assert ask(-x, Q.negative, ~Q.positive(x)) == None
    assert ask(-x, Q.negative, Q.negative(x)) == False
    assert ask(-x, Q.negative, Q.positive(x)) == True

    assert ask(x-1, Q.negative, Q.negative(x)) == True
    assert ask(x+y, Q.negative) == None
    assert ask(x+y, Q.negative, Q.negative(x)) == None
    assert ask(x+y, Q.negative, Q.negative(x) & Q.negative(y)) == True

    assert ask(x**2, Q.negative) == None
    assert ask(x**2, Q.negative, Q.real(x)) == False
    assert ask(x**1.4, Q.negative, Q.real(x)) == None

    assert ask(x*y, Q.negative) == None
    assert ask(x*y, Q.negative, Q.positive(x) & Q.positive(y)) == False
    assert ask(x*y, Q.negative, Q.positive(x) & Q.negative(y)) == True
    assert ask(x*y, Q.negative, Q.complex(x) & Q.complex(y)) == None

    assert ask(x**y, Q.negative) == None
    assert ask(x**y, Q.negative, Q.negative(x) & Q.even(y)) == False
    assert ask(x**y, Q.negative, Q.negative(x) & Q.odd(y)) == True
    assert ask(x**y, Q.negative, Q.positive(x) & Q.integer(y)) == False

    assert ask(Abs(x), Q.negative) == False

def test_nonzero():
    x, y = symbols('x,y')
    assert ask(x, Q.nonzero) == None
    assert ask(x, Q.nonzero, Q.real(x)) == None
    assert ask(x, Q.nonzero, Q.positive(x)) == True
    assert ask(x, Q.nonzero, Q.negative(x)) == True
    assert ask(x, Q.nonzero, Q.negative(x) | Q.positive(x)) == True

    assert ask(x+y, Q.nonzero) == None
    assert ask(x+y, Q.nonzero, Q.positive(x) & Q.positive(y)) == True
    assert ask(x+y, Q.nonzero, Q.positive(x) & Q.negative(y)) == None
    assert ask(x+y, Q.nonzero, Q.negative(x) & Q.negative(y)) == True

    assert ask(2*x, Q.nonzero) == None
    assert ask(2*x, Q.nonzero, Q.positive(x)) == True
    assert ask(2*x, Q.nonzero, Q.negative(x)) == True
    assert ask(x*y, Q.nonzero, Q.nonzero(x)) == None
    assert ask(x*y, Q.nonzero, Q.nonzero(x) & Q.nonzero(y)) == True

    assert ask(Abs(x), Q.nonzero) == None
    assert ask(Abs(x), Q.nonzero, Q.nonzero(x)) == True

def test_odd():
    x, y, z, t = symbols('x,y,z,t')
    assert ask(x, Q.odd) == None
    assert ask(x, Q.odd, Q.odd(x)) == True
    assert ask(x, Q.odd, Q.integer(x)) == None
    assert ask(x, Q.odd, ~Q.integer(x)) == False
    assert ask(x, Q.odd, Q.rational(x)) == None
    assert ask(x, Q.odd, Q.positive(x)) == None

    assert ask(-x, Q.odd, Q.odd(x)) == True

    assert ask(2*x, Q.odd) == None
    assert ask(2*x, Q.odd, Q.integer(x)) == False
    assert ask(2*x, Q.odd, Q.odd(x)) == False
    assert ask(2*x, Q.odd, Q.irrational(x)) == False
    assert ask(2*x, Q.odd, ~Q.integer(x)) == None
    assert ask(3*x, Q.odd, Q.integer(x)) == None

    assert ask(x/3, Q.odd, Q.odd(x)) == None
    assert ask(x/3, Q.odd, Q.even(x)) == None

    assert ask(x+1, Q.odd, Q.even(x)) == True
    assert ask(x+2, Q.odd, Q.even(x)) == False
    assert ask(x+2, Q.odd, Q.odd(x))  == True
    assert ask(3-x, Q.odd, Q.odd(x))  == False
    assert ask(3-x, Q.odd, Q.even(x))  == True
    assert ask(3+x, Q.odd, Q.odd(x))  == False
    assert ask(3+x, Q.odd, Q.even(x))  == True
    assert ask(x+y, Q.odd, Q.odd(x) & Q.odd(y)) == False
    assert ask(x+y, Q.odd, Q.odd(x) & Q.even(y)) == True
    assert ask(x-y, Q.odd, Q.even(x) & Q.odd(y)) == True
    assert ask(x-y, Q.odd, Q.odd(x) & Q.odd(y)) == False

    assert ask(x+y+z, Q.odd, Q.odd(x) & Q.odd(y) & Q.even(z)) == False
    assert ask(x+y+z+t, Q.odd,
               Q.odd(x) & Q.odd(y) & Q.even(z) & Q.integer(t)) == None

    assert ask(2*x + 1, Q.odd, Q.integer(x)) == True
    assert ask(2*x + y, Q.odd, Q.integer(x) & Q.odd(y)) == True
    assert ask(2*x + y, Q.odd, Q.integer(x) & Q.even(y)) == False
    assert ask(2*x + y, Q.odd, Q.integer(x) & Q.integer(y)) == None
    assert ask(x*y,   Q.odd, Q.odd(x) & Q.even(y)) == False
    assert ask(x*y,   Q.odd, Q.odd(x) & Q.odd(y)) == True
    assert ask(2*x*y, Q.odd, Q.rational(x) & Q.rational(x)) == None
    assert ask(2*x*y, Q.odd, Q.irrational(x) & Q.irrational(x)) == None

    assert ask(Abs(x), Q.odd, Q.odd(x)) == True

def test_prime():
    x, y = symbols('x,y')
    assert ask(x, Q.prime, Q.prime(x)) == True
    assert ask(x, Q.prime, ~Q.prime(x)) == False
    assert ask(x, Q.prime, Q.integer(x)) == None
    assert ask(x, Q.prime, ~Q.integer(x)) == False

    assert ask(2*x, Q.prime, Q.integer(x)) == False
    assert ask(x*y, Q.prime) == None
    assert ask(x*y, Q.prime, Q.prime(x)) == None
    assert ask(x*y, Q.prime, Q.integer(x) & Q.integer(y)) == False

    assert ask(x**2, Q.prime, Q.integer(x)) == False
    assert ask(x**2, Q.prime, Q.prime(x)) == False
    assert ask(x**y, Q.prime, Q.integer(x) & Q.integer(y)) == False

def test_positive():
    x, y, z, w = symbols('x,y,z,w')
    assert ask(x, Q.positive, Q.positive(x)) == True
    assert ask(x, Q.positive, Q.negative(x)) == False
    assert ask(x, Q.positive, Q.nonzero(x)) == None

    assert ask(-x, Q.positive, Q.positive(x)) == False
    assert ask(-x, Q.positive, Q.negative(x)) == True

    assert ask(x+y, Q.positive, Q.positive(x) & Q.positive(y)) == True
    assert ask(x+y, Q.positive, Q.positive(x) & Q.negative(y)) == None

    assert ask(2*x,  Q.positive, Q.positive(x)) == True
    assumptions =  Q.positive(x) & Q.negative(y) & Q.negative(z) & Q.positive(w)
    assert ask(x*y*z,  Q.positive)  == None
    assert ask(x*y*z,  Q.positive, assumptions) == True
    assert ask(-x*y*z, Q.positive, assumptions) == False

    assert ask(x**2, Q.positive, Q.positive(x)) == True
    assert ask(x**2, Q.positive, Q.negative(x)) == True

    #exponential
    assert ask(exp(x),     Q.positive, Q.real(x)) == True
    assert ask(x + exp(x), Q.positive, Q.real(x)) == None

    #absolute value
    assert ask(Abs(x), Q.positive) == None # Abs(0) = 0
    assert ask(Abs(x), Q.positive, Q.positive(x)) == True

@XFAIL
def test_positive_xfail():
    assert ask(1/(1 + x**2), Q.positive, Q.real(x)) == True

def test_real():
    x, y = symbols('x,y')
    assert ask(x, Q.real) == None
    assert ask(x, Q.real, Q.real(x)) == True
    assert ask(x, Q.real, Q.nonzero(x)) == True
    assert ask(x, Q.real, Q.positive(x)) == True
    assert ask(x, Q.real, Q.negative(x)) == True
    assert ask(x, Q.real, Q.integer(x)) == True
    assert ask(x, Q.real, Q.even(x)) == True
    assert ask(x, Q.real, Q.prime(x)) == True

    assert ask(x/sqrt(2), Q.real, Q.real(x)) == True
    assert ask(x/sqrt(-2), Q.real, Q.real(x)) == False

    I = S.ImaginaryUnit
    assert ask(x+1, Q.real, Q.real(x)) == True
    assert ask(x+I, Q.real, Q.real(x)) == False
    assert ask(x+I, Q.real, Q.complex(x)) == None

    assert ask(2*x, Q.real, Q.real(x)) == True
    assert ask(I*x, Q.real, Q.real(x)) == False
    assert ask(I*x, Q.real, Q.imaginary(x)) == True
    assert ask(I*x, Q.real, Q.complex(x)) == None

    assert ask(x**2, Q.real, Q.real(x)) == True
    assert ask(sqrt(x), Q.real, Q.negative(x)) == False
    assert ask(x**y, Q.real, Q.real(x) & Q.integer(y)) == True
    assert ask(x**y, Q.real, Q.real(x) & Q.real(y)) == None
    assert ask(x**y, Q.real, Q.positive(x) & Q.real(y)) == True

    # trigonometric functions
    assert ask(sin(x), Q.real) == None
    assert ask(cos(x), Q.real) == None
    assert ask(sin(x), Q.real, Q.real(x)) == True
    assert ask(cos(x), Q.real, Q.real(x)) == True

    # exponential function
    assert ask(exp(x), Q.real) == None
    assert ask(exp(x), Q.real, Q.real(x)) == True
    assert ask(x + exp(x), Q.real, Q.real(x)) == True

    # Q.complexes
    assert ask(re(x), Q.real) == True
    assert ask(im(x), Q.real) == True

def test_algebraic():
    x, y = symbols('x,y')

    assert ask(x, 'algebraic') == None

    assert ask(I, 'algebraic') == True
    assert ask(2*I, 'algebraic') == True
    assert ask(I/3, 'algebraic') == True

    assert ask(sqrt(7), 'algebraic') == True
    assert ask(2*sqrt(7), 'algebraic') == True
    assert ask(sqrt(7)/3, 'algebraic') == True

    assert ask(I*sqrt(3), 'algebraic') == True
    assert ask(sqrt(1+I*sqrt(3)), 'algebraic') == True

    assert ask((1+I*sqrt(3)**(S(17)/31)), 'algebraic') == True
    assert ask((1+I*sqrt(3)**(S(17)/pi)), 'algebraic') == False

    assert ask(sin(7), 'algebraic') == None
    assert ask(sqrt(sin(7)), 'algebraic') == None
    assert ask(sqrt(y+I*sqrt(7)), 'algebraic') == None

    assert ask(oo, 'algebraic') == False
    assert ask(-oo, 'algebraic') == False

    assert ask(2.47, 'algebraic') == False

def test_global():
    """Test ask with global assumptions"""
    x = symbols('x')
    assert ask(x, Q.integer) == None
    global_assumptions.add(Q.integer(x))
    assert ask(x, Q.integer) == True
    global_assumptions.clear()
    assert ask(x, Q.integer) == None

def test_custom_context():
    """Test ask with custom assumptions context"""
    x = symbols('x')
    assert ask(x, Q.integer) == None
    local_context = AssumptionsContext()
    local_context.add(Q.integer(x))
    assert ask(x, Q.integer, context = local_context) == True
    assert ask(x, Q.integer) == None

def test_functions_in_assumptions():
    from sympy.logic.boolalg import Equivalent, Xor
    x = symbols('x')
    assert ask(x, Q.negative, Q.real(x) >> Q.positive(x)) is False
    assert ask(x, Q.negative, Equivalent(Q.real(x), Q.positive(x))) is False
    assert ask(x, Q.negative, Xor(Q.real(x), Q.negative(x))) is False

def test_composite_ask_key():
    x = symbols('x')
    assert ask(x, Q.negative & Q.integer, Q.real(x) >> Q.positive(x)) is False

def test_is_true():
    from sympy.logic.boolalg import Equivalent, Implies
    x = symbols('x')
    assert ask(True, Q.is_true) is True
    assert ask(~Q.negative(x), Q.is_true, Q.positive(x)) is True
    assert ask(~Q.real(x), Q.is_true, Q.commutative(x)) is None
    assert ask(Q.negative(x) & Q.integer(x), Q.is_true, Q.positive(x)) is False
    assert ask(Q.negative(x) & Q.integer(x), Q.is_true) is None
    assert ask(Q.real(x) | Q.integer(x), Q.is_true, Q.positive(x)) is True
    assert ask(Q.real(x) | Q.integer(x), Q.is_true) is None
    assert ask(Q.real(x) >> Q.positive(x), Q.is_true, Q.negative(x)) is False
    assert ask(Implies(Q.real(x), Q.positive(x), evaluate=False), Q.is_true,
                    Q.negative(x)) is False
    assert ask(Implies(Q.real(x), Q.positive(x), evaluate=False), Q.is_true) is None
    assert ask(Equivalent(Q.integer(x), Q.even(x)), Q.is_true, Q.even(x)) is True
    assert ask(Equivalent(Q.integer(x), Q.even(x)), Q.is_true) is None
    assert ask(Equivalent(Q.positive(x), Q.integer(x)), Q.is_true, Q.integer(x)) is None

def test_incompatible_resolutors():
    x = symbols('x')
    class Prime2AskHandler(AskHandler):
        @staticmethod
        def Number(expr, assumptions):
            return True
    register_handler('prime', Prime2AskHandler)
    raises(ValueError, 'ask(4, Q.prime)')
    remove_handler('prime', Prime2AskHandler)

    class InconclusiveHandler(AskHandler):
        @staticmethod
        def Number(expr, assumptions):
            return None
    register_handler('prime', InconclusiveHandler)
    assert ask(3, Q.prime) == True


def test_key_extensibility():
    """test that you can add keys to the ask system at runtime"""
    x = Symbol('x')
    # make sure thie key is not defined
    raises(AttributeError, "ask(x, 'my_key')")
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

def test_compute_known_facts():
    ns = {}
    exec 'from sympy.logic.boolalg import And, Or, Not' in globals(), ns
    exec compute_known_facts() in globals(), ns
    assert ns['known_facts_cnf'] == known_facts_cnf
    assert ns['known_facts_dict'] == known_facts_dict
