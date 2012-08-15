from sympy.abc import t, w, x, y, z
from sympy.assumptions import (ask, AssumptionsContext, global_assumptions, Q,
                               register_handler, remove_handler)
from sympy.assumptions.ask import (compute_known_facts, known_facts_cnf,
                                   known_facts_dict)
from sympy.assumptions.handlers import AskHandler
from sympy.core import I, Integer, oo, pi, Rational, S, symbols, Add
from sympy.functions import Abs, cos, exp, im, log, re, sign, sin, sqrt
from sympy.logic import Equivalent, Implies, Xor
from sympy.utilities.pytest import raises, XFAIL, slow

def test_int_1():
    z = 1
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == True
    assert ask(Q.rational(z))         == True
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == False
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == True
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == True
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == True
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

def test_int_11():
    z = 11
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == True
    assert ask(Q.rational(z))         == True
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == False
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == True
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == True
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == True
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

def test_int_12():
    z = 12
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == True
    assert ask(Q.rational(z))         == True
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == False
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == True
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == True
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == True
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

def test_float_1():
    z = 1.0
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == True
    assert ask(Q.rational(z))         == True
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == False
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == True
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == True
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == True
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

    z = 7.2123
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == False
    assert ask(Q.rational(z))         == True
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == False
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == True
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

def test_zero_0():
    z = Integer(0)
    assert ask(Q.nonzero(z))          == False
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == True
    assert ask(Q.rational(z))         == True
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == False
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == True
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == True
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

def test_negativeone():
    z = Integer(-1)
    assert ask(Q.nonzero(z))          == True
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == True
    assert ask(Q.rational(z))         == True
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == False
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == False
    assert ask(Q.negative(z))         == True
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == True
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

def test_infinity():
    assert ask(Q.commutative(oo))     == True
    assert ask(Q.integer(oo))         == False
    assert ask(Q.rational(oo))        == False
    assert ask(Q.real(oo))            == False
    assert ask(Q.extended_real(oo))   == True
    assert ask(Q.complex(oo))         == False
    assert ask(Q.irrational(oo))      == False
    assert ask(Q.imaginary(oo))       == False
    assert ask(Q.positive(oo))        == True
    assert ask(Q.negative(oo))        == False
    assert ask(Q.even(oo))            == False
    assert ask(Q.odd(oo))             == False
    assert ask(Q.bounded(oo))         == False
    assert ask(Q.infinitesimal(oo))   == False
    assert ask(Q.prime(oo))           == False
    assert ask(Q.composite(oo))       == False
    assert ask(Q.hermitian(oo))       == False
    assert ask(Q.antihermitian(oo))   == False

def test_neg_infinity():
    mm = S.NegativeInfinity
    assert ask(Q.commutative(mm))    == True
    assert ask(Q.integer(mm))        == False
    assert ask(Q.rational(mm))       == False
    assert ask(Q.real(mm))           == False
    assert ask(Q.extended_real(mm))  == True
    assert ask(Q.complex(mm))        == False
    assert ask(Q.irrational(mm))     == False
    assert ask(Q.imaginary(mm))      == False
    assert ask(Q.positive(mm))       == False
    assert ask(Q.negative(mm))       == True
    assert ask(Q.even(mm))           == False
    assert ask(Q.odd(mm))            == False
    assert ask(Q.bounded(mm))        == False
    assert ask(Q.infinitesimal(mm))  == False
    assert ask(Q.prime(mm))          == False
    assert ask(Q.composite(mm))      == False
    assert ask(Q.hermitian(mm))      == False
    assert ask(Q.antihermitian(mm))  == False

def test_nan():
    nan = S.NaN
    assert ask(Q.commutative(nan))   == True
    assert ask(Q.integer(nan))       == False
    assert ask(Q.rational(nan))      == False
    assert ask(Q.real(nan))          == False
    assert ask(Q.extended_real(nan)) == False
    assert ask(Q.complex(nan))       == False
    assert ask(Q.irrational(nan))    == False
    assert ask(Q.imaginary(nan))     == False
    assert ask(Q.positive(nan))      == False
    assert ask(Q.nonzero(nan))       == True
    assert ask(Q.even(nan))          == False
    assert ask(Q.odd(nan))           == False
    assert ask(Q.bounded(nan))       == False
    assert ask(Q.infinitesimal(nan)) == False
    assert ask(Q.prime(nan))         == False
    assert ask(Q.composite(nan))     == False
    assert ask(Q.hermitian(nan))     == False
    assert ask(Q.antihermitian(nan)) == False

def test_Rational_number():
    r = Rational(3,4)
    assert ask(Q.commutative(r))      == True
    assert ask(Q.integer(r))          == False
    assert ask(Q.rational(r))         == True
    assert ask(Q.real(r))             == True
    assert ask(Q.complex(r))          == True
    assert ask(Q.irrational(r))       == False
    assert ask(Q.imaginary(r))        == False
    assert ask(Q.positive(r))         == True
    assert ask(Q.negative(r))         == False
    assert ask(Q.even(r))             == False
    assert ask(Q.odd(r))              == False
    assert ask(Q.bounded(r))          == True
    assert ask(Q.infinitesimal(r))    == False
    assert ask(Q.prime(r))            == False
    assert ask(Q.composite(r))        == False
    assert ask(Q.hermitian(r))        == True
    assert ask(Q.antihermitian(r))    == False

    r = Rational(1,4)
    assert ask(Q.positive(r))         == True
    assert ask(Q.negative(r))         == False

    r = Rational(5,4)
    assert ask(Q.negative(r))         == False
    assert ask(Q.positive(r))         == True

    r = Rational(5,3)
    assert ask(Q.positive(r))         == True
    assert ask(Q.negative(r))         == False

    r = Rational(-3,4)
    assert ask(Q.positive(r))         == False
    assert ask(Q.negative(r))         == True

    r = Rational(-1,4)
    assert ask(Q.positive(r))         == False
    assert ask(Q.negative(r))         == True

    r = Rational(-5,4)
    assert ask(Q.negative(r))         == True
    assert ask(Q.positive(r))         == False

    r = Rational(-5,3)
    assert ask(Q.positive(r))         == False
    assert ask(Q.negative(r))         == True

def test_sqrt_2():
    z = sqrt(2)
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == False
    assert ask(Q.rational(z))         == False
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == True
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == True
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

def test_pi():
    z = S.Pi
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == False
    assert ask(Q.rational(z))         == False
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == True
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == True
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

    z = S.Pi + 1
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == False
    assert ask(Q.rational(z))         == False
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == True
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == True
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

    z = 2*S.Pi
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == False
    assert ask(Q.rational(z))         == False
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == True
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == True
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

    z = S.Pi ** 2
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == False
    assert ask(Q.rational(z))         == False
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == True
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == True
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

    z = (1+S.Pi) ** 2
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == False
    assert ask(Q.rational(z))         == False
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == True
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == True
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

def test_E():
    z = S.Exp1
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == False
    assert ask(Q.rational(z))         == False
    assert ask(Q.real(z))             == True
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == True
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == True
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == True
    assert ask(Q.antihermitian(z))    == False

def test_I():
    z = I
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == False
    assert ask(Q.rational(z))         == False
    assert ask(Q.real(z))             == False
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == False
    assert ask(Q.imaginary(z))        == True
    assert ask(Q.positive(z))         == False
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == False
    assert ask(Q.antihermitian(z))    == True

    z = 1 + I
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == False
    assert ask(Q.rational(z))         == False
    assert ask(Q.real(z))             == False
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == False
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == False
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == False
    assert ask(Q.antihermitian(z))    == False

    z = I*(1+I)
    assert ask(Q.commutative(z))      == True
    assert ask(Q.integer(z))          == False
    assert ask(Q.rational(z))         == False
    assert ask(Q.real(z))             == False
    assert ask(Q.complex(z))          == True
    assert ask(Q.irrational(z))       == False
    assert ask(Q.imaginary(z))        == False
    assert ask(Q.positive(z))         == False
    assert ask(Q.negative(z))         == False
    assert ask(Q.even(z))             == False
    assert ask(Q.odd(z))              == False
    assert ask(Q.bounded(z))          == True
    assert ask(Q.infinitesimal(z))    == False
    assert ask(Q.prime(z))            == False
    assert ask(Q.composite(z))        == False
    assert ask(Q.hermitian(z))        == False
    assert ask(Q.antihermitian(z))    == False

@slow
def test_bounded():
    x, y, z = symbols('x,y,z')
    assert ask(Q.bounded(x)) == None
    assert ask(Q.bounded(x), Q.bounded(x)) == True
    assert ask(Q.bounded(x), Q.bounded(y)) == None
    assert ask(Q.bounded(x), Q.complex(x)) == None

    assert ask(Q.bounded(x+1)) == None
    assert ask(Q.bounded(x+1), Q.bounded(x)) == True
    a = x + y
    x, y = a.args
    # B + B
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y)) == True
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & Q.positive(x)) == True
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & Q.positive(y)) == True
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & Q.positive(x) & Q.positive(y)) == True
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & Q.positive(x) & ~Q.positive(y)) == True
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & ~Q.positive(x) & Q.positive(y)) == True
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & ~Q.positive(x) & ~Q.positive(y)) == True
    # B + U
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y) & Q.positive(x)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y) & Q.positive(y)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y) & Q.positive(x) & Q.positive(y)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y) & Q.positive(x) & ~Q.positive(y)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y) & ~Q.positive(x) & Q.positive(y)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y) & ~Q.positive(x) & ~Q.positive(y)) == False
    # B + ?
    assert ask(Q.bounded(a), Q.bounded(x)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(x)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(x) & Q.positive(y)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(x) & ~Q.positive(y)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.positive(x) & Q.positive(y)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.positive(x) & ~Q.positive(y)) == None
    # U + U
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y) & Q.positive(x)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y) & Q.positive(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y) & Q.positive(x) & Q.positive(y)) == False
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y) & Q.positive(x) & ~Q.positive(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y) & ~Q.positive(x) & Q.positive(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y) & ~Q.positive(x) & ~Q.positive(y)) == False
    # U + ?
    assert ask(Q.bounded(a), ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(y) & Q.positive(x)) == None
    assert ask(Q.bounded(a), ~Q.bounded(y) & Q.positive(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(y) & Q.positive(x) & Q.positive(y)) == False
    assert ask(Q.bounded(a), ~Q.bounded(y) & Q.positive(x) & ~Q.positive(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(y) & ~Q.positive(x) & Q.positive(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(y) & ~Q.positive(x) & ~Q.positive(y)) == False
    # ? + ?
    assert ask(Q.bounded(a),) == None
    assert ask(Q.bounded(a),Q.positive(x)) == None
    assert ask(Q.bounded(a),Q.positive(y)) == None
    assert ask(Q.bounded(a),Q.positive(x) & Q.positive(y)) == None
    assert ask(Q.bounded(a),Q.positive(x) & ~Q.positive(y)) == None
    assert ask(Q.bounded(a),~Q.positive(x) & Q.positive(y)) == None
    assert ask(Q.bounded(a),~Q.positive(x) & ~Q.positive(y)) == None
    a = x + y + z
    x, y, z = a.args
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & Q.bounded(y) & Q.negative(z) & Q.bounded(z)) == True
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & Q.bounded(y) & Q.bounded(z)) == True
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & Q.bounded(y) & Q.positive(z) & Q.bounded(z)) == True
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & Q.bounded(y) & Q.negative(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & Q.bounded(y) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.bounded(y) & Q.bounded(z)) == True
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.bounded(y) & Q.positive(z) & Q.bounded(z)) == True
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.bounded(y) & Q.negative(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.bounded(y) & ~Q.bounded(z))== False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.bounded(y) & Q.negative(z))== None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y) & Q.positive(z) & Q.bounded(z)) == True
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y) & Q.negative(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & ~Q.bounded(y) & Q.negative(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & ~Q.bounded(y) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & ~Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & ~Q.bounded(y) & Q.negative(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & ~Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & ~Q.bounded(y) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & ~Q.bounded(y) & Q.positive(z)& ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & ~Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & ~Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(y) & ~Q.bounded(y)& Q.positive(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(y) & ~Q.bounded(y)& Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(y) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(y) & ~Q.bounded(y)& Q.positive(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.negative(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.bounded(x) & Q.positive(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & Q.bounded(z)) == True
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & Q.positive(z) & Q.bounded(z))== True
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & Q.negative(z) & ~Q.bounded(z))== False
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y) & Q.bounded(y) & Q.positive(z)& Q.bounded(z)) == True
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y) & Q.bounded(y) & Q.negative(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y) & Q.bounded(y) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y) & Q.bounded(y) & Q.positive(z)& ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y) & Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y) & Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y) & Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.negative(y) & ~Q.bounded(y) & Q.negative(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & Q.negative(y) & ~Q.bounded(y) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.negative(y) & ~Q.bounded(y) & Q.positive(z)& ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.negative(y) & ~Q.bounded(y) & Q.negative(z)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.negative(y) & ~Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.positive(z)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & Q.negative(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.negative(y)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.negative(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.positive(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y)& Q.positive(z) & Q.bounded(z)) == True
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y)& Q.negative(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y)& ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y)& Q.positive(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y)& Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(y) & Q.bounded(y)& Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)& Q.negative(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)& ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)& Q.positive(z) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)& Q.negative(z)) == False
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)& Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & ~Q.bounded(y) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & ~Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & ~Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & ~Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(y) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.positive(z)) == False
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.negative(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.negative(y)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.negative(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.bounded(x) & Q.positive(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)& Q.negative(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)& ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)& Q.positive(z) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)& Q.negative(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.negative(y) & ~Q.bounded(y)& Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & ~Q.bounded(y) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & ~Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & ~Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & ~Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.positive(y) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.negative(y) & Q.negative(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.negative(y)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.negative(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & ~Q.bounded(x) & Q.positive(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.positive(y) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.negative(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.negative(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.negative(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.positive(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.positive(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & ~Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.positive(z) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.positive(x) & ~Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & ~Q.bounded(x) & Q.positive(y) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.positive(x) & ~Q.bounded(x) & Q.positive(y) & ~Q.bounded(y) & Q.positive(z)) == False
    assert ask(Q.bounded(a), Q.positive(x) & ~Q.bounded(x) & Q.negative(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & ~Q.bounded(x) & Q.negative(y)) == None
    assert ask(Q.bounded(a), Q.positive(x) & ~Q.bounded(x) & Q.negative(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & ~Q.bounded(x)) == None
    assert ask(Q.bounded(a), Q.positive(x) & ~Q.bounded(x) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & ~Q.bounded(x) & Q.positive(y) & Q.positive(z)) == False
    assert ask(Q.bounded(a), Q.negative(x) & Q.negative(y) & Q.negative(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.negative(y)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.negative(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.negative(x) & Q.positive(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a)) == None
    assert ask(Q.bounded(a), Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.positive(y) & Q.positive(z)) == None
    assert ask(Q.bounded(a), Q.positive(x) & Q.positive(y) & Q.positive(z)) == None

    x, y, z = symbols('x,y,z')
    assert ask(Q.bounded(2*x)) == None
    assert ask(Q.bounded(2*x), Q.bounded(x)) == True
    a = x*y
    x, y = a.args
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y)) == True
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y)) == False
    assert ask(Q.bounded(a), Q.bounded(x)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.bounded(y)) == False
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y)) == False
    assert ask(Q.bounded(a), ~Q.bounded(x)) == None
    assert ask(Q.bounded(a), Q.bounded(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(y)) == None
    assert ask(Q.bounded(a)) == None
    a = x*y*z
    x, y, z = a.args
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & Q.bounded(z)) == True
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y) & Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.bounded(x)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.bounded(y) & Q.bounded(z)) == False
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.bounded(y) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.bounded(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y) & Q.bounded(z)) == False
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y) & ~Q.bounded(z)) == False
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & Q.bounded(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(x)) == None
    assert ask(Q.bounded(a), Q.bounded(y) & Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.bounded(y) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), Q.bounded(y)) == None
    assert ask(Q.bounded(a), ~Q.bounded(y) & Q.bounded(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(y) & ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(y)) == None
    assert ask(Q.bounded(a), Q.bounded(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(z) & Q.nonzero(x) & Q.nonzero(y) & Q.nonzero(z)) == None
    assert ask(Q.bounded(a), ~Q.bounded(y) & ~Q.bounded(z) & Q.nonzero(x) & Q.nonzero(y) & Q.nonzero(z)) == False

    x, y, z = symbols('x,y,z')
    assert ask(Q.bounded(x**2)) == None
    assert ask(Q.bounded(2**x)) == None
    assert ask(Q.bounded(2**x), Q.bounded(x)) == True
    assert ask(Q.bounded(x**x)) == None
    assert ask(Q.bounded(Rational(1,2) ** x)) == None
    assert ask(Q.bounded(Rational(1,2) ** x), Q.positive(x)) == True
    assert ask(Q.bounded(Rational(1,2) ** x), Q.negative(x)) == None
    assert ask(Q.bounded(S(2) ** x), Q.negative(x)) == True
    assert ask(Q.bounded(sqrt(x))) == None
    assert ask(Q.bounded(2**x), ~Q.bounded(x))==False
    assert ask(Q.bounded(x**2), ~Q.bounded(x))==False

    # sign function
    assert ask(Q.bounded(sign(x))) == True
    assert ask(Q.bounded(sign(x)), ~Q.bounded(x)) == True

    # exponential functions
    assert ask(Q.bounded(log(x))) == None
    assert ask(Q.bounded(log(x)), Q.bounded(x)) == True
    assert ask(Q.bounded(exp(x))) == None
    assert ask(Q.bounded(exp(x)), Q.bounded(x)) == True
    assert ask(Q.bounded(exp(2))) == True

    # trigonometric functions
    assert ask(Q.bounded(sin(x))) == True
    assert ask(Q.bounded(sin(x)), ~Q.bounded(x)) == True
    assert ask(Q.bounded(cos(x))) == True
    assert ask(Q.bounded(cos(x)), ~Q.bounded(x)) == True
    assert ask(Q.bounded(2*sin(x))) == True
    assert ask(Q.bounded(sin(x)**2)) == True
    assert ask(Q.bounded(cos(x)**2)) == True
    assert ask(Q.bounded(cos(x) + sin(x))) == True

@XFAIL
def test_bounded_xfail():
    """We need to support relations in ask for this to work"""
    assert ask(Q.bounded(sin(x)**x)) == True
    assert ask(Q.bounded(cos(x)**x)) == True
    assert ask(Q.bounded(sin(x) ** x)) == True

def test_commutative():
    """By default objects are Q.commutative that is why it returns True
    for both key=True and key=False"""
    assert ask(Q.commutative(x)) == True
    assert ask(Q.commutative(x), ~Q.commutative(x)) == False
    assert ask(Q.commutative(x), Q.complex(x)) == True
    assert ask(Q.commutative(x), Q.imaginary(x)) == True
    assert ask(Q.commutative(x), Q.real(x)) == True
    assert ask(Q.commutative(x), Q.positive(x)) == True
    assert ask(Q.commutative(x), ~Q.commutative(y))  == True

    assert ask(Q.commutative(2*x)) == True
    assert ask(Q.commutative(2*x), ~Q.commutative(x)) == False

    assert ask(Q.commutative(x + 1)) == True
    assert ask(Q.commutative(x + 1), ~Q.commutative(x)) == False

    assert ask(Q.commutative(x**2)) == True
    assert ask(Q.commutative(x**2), ~Q.commutative(x)) == False

    assert ask(Q.commutative(log(x))) == True

def test_complex():
    assert ask(Q.complex(x)) == None
    assert ask(Q.complex(x), Q.complex(x)) == True
    assert ask(Q.complex(x), Q.complex(y)) == None
    assert ask(Q.complex(x), ~Q.complex(x)) == False
    assert ask(Q.complex(x), Q.real(x)) == True
    assert ask(Q.complex(x), ~Q.real(x)) == None
    assert ask(Q.complex(x), Q.rational(x)) == True
    assert ask(Q.complex(x), Q.irrational(x)) == True
    assert ask(Q.complex(x), Q.positive(x)) == True
    assert ask(Q.complex(x), Q.imaginary(x)) == True

    # a+b
    assert ask(Q.complex(x+1), Q.complex(x)) == True
    assert ask(Q.complex(x+1), Q.real(x)) == True
    assert ask(Q.complex(x+1), Q.rational(x)) == True
    assert ask(Q.complex(x+1), Q.irrational(x)) == True
    assert ask(Q.complex(x+1), Q.imaginary(x)) == True
    assert ask(Q.complex(x+1), Q.integer(x))  == True
    assert ask(Q.complex(x+1), Q.even(x))  == True
    assert ask(Q.complex(x+1), Q.odd(x))  == True
    assert ask(Q.complex(x+y), Q.complex(x) & Q.complex(y)) == True
    assert ask(Q.complex(x+y), Q.real(x) & Q.imaginary(y)) == True

    # a*x +b
    assert ask(Q.complex(2*x+1), Q.complex(x)) == True
    assert ask(Q.complex(2*x+1), Q.real(x)) == True
    assert ask(Q.complex(2*x+1), Q.positive(x)) == True
    assert ask(Q.complex(2*x+1), Q.rational(x)) == True
    assert ask(Q.complex(2*x+1), Q.irrational(x)) == True
    assert ask(Q.complex(2*x+1), Q.imaginary(x)) == True
    assert ask(Q.complex(2*x+1), Q.integer(x))  == True
    assert ask(Q.complex(2*x+1), Q.even(x))  == True
    assert ask(Q.complex(2*x+1), Q.odd(x))  == True

    # x**2
    assert ask(Q.complex(x**2), Q.complex(x)) == True
    assert ask(Q.complex(x**2), Q.real(x)) == True
    assert ask(Q.complex(x**2), Q.positive(x)) == True
    assert ask(Q.complex(x**2), Q.rational(x)) == True
    assert ask(Q.complex(x**2), Q.irrational(x)) == True
    assert ask(Q.complex(x**2), Q.imaginary(x)) == True
    assert ask(Q.complex(x**2), Q.integer(x))  == True
    assert ask(Q.complex(x**2), Q.even(x))  == True
    assert ask(Q.complex(x**2), Q.odd(x))  == True

    # 2**x
    assert ask(Q.complex(2**x), Q.complex(x)) == True
    assert ask(Q.complex(2**x), Q.real(x)) == True
    assert ask(Q.complex(2**x), Q.positive(x)) == True
    assert ask(Q.complex(2**x), Q.rational(x)) == True
    assert ask(Q.complex(2**x), Q.irrational(x)) == True
    assert ask(Q.complex(2**x), Q.imaginary(x)) == True
    assert ask(Q.complex(2**x), Q.integer(x))  == True
    assert ask(Q.complex(2**x), Q.even(x))  == True
    assert ask(Q.complex(2**x), Q.odd(x))  == True
    assert ask(Q.complex(x**y), Q.complex(x) & Q.complex(y)) == True

    # trigonometric expressions
    assert ask(Q.complex(sin(x))) == True
    assert ask(Q.complex(sin(2*x + 1))) == True
    assert ask(Q.complex(cos(x))) == True
    assert ask(Q.complex(cos(2*x+1))) == True

    # exponential
    assert ask(Q.complex(exp(x))) == True
    assert ask(Q.complex(exp(x))) == True

    # Q.complexes
    assert ask(Q.complex(Abs(x))) == True
    assert ask(Q.complex(re(x))) == True
    assert ask(Q.complex(im(x))) == True

def test_even():
    assert ask(Q.even(x)) == None
    assert ask(Q.even(x), Q.integer(x)) == None
    assert ask(Q.even(x), ~Q.integer(x)) == False
    assert ask(Q.even(x), Q.rational(x)) == None
    assert ask(Q.even(x), Q.positive(x)) == None

    assert ask(Q.even(2*x)) == None
    assert ask(Q.even(2*x), Q.integer(x)) == True
    assert ask(Q.even(2*x), Q.even(x)) == True
    assert ask(Q.even(2*x), Q.irrational(x)) == False
    assert ask(Q.even(2*x), Q.odd(x)) == True
    assert ask(Q.even(2*x), ~Q.integer(x)) == None
    assert ask(Q.even(3*x), Q.integer(x)) == None
    assert ask(Q.even(3*x), Q.even(x)) == True
    assert ask(Q.even(3*x), Q.odd(x)) == False

    assert ask(Q.even(x+1), Q.odd(x)) == True
    assert ask(Q.even(x+1), Q.even(x)) == False
    assert ask(Q.even(x+2), Q.odd(x)) == False
    assert ask(Q.even(x+2), Q.even(x)) == True
    assert ask(Q.even(7-x), Q.odd(x)) == True
    assert ask(Q.even(7+x), Q.odd(x)) == True
    assert ask(Q.even(x+y), Q.odd(x) & Q.odd(y)) == True
    assert ask(Q.even(x+y), Q.odd(x) & Q.even(y)) == False
    assert ask(Q.even(x+y), Q.even(x) & Q.even(y)) == True

    assert ask(Q.even(2*x + 1), Q.integer(x)) == False
    assert ask(Q.even(2*x*y), Q.rational(x) & Q.rational(x)) == None
    assert ask(Q.even(2*x*y), Q.irrational(x) & Q.irrational(x)) == None

    assert ask(Q.even(x+y+z), Q.odd(x) & Q.odd(y) & Q.even(z)) == True
    assert ask(Q.even(x+y+z+t), Q.odd(x) & Q.odd(y) & Q.even(z) & Q.integer(t)) == None

    assert ask(Q.even(Abs(x)), Q.even(x)) == True
    assert ask(Q.even(Abs(x)), ~Q.even(x)) == None
    assert ask(Q.even(re(x)), Q.even(x)) == True
    assert ask(Q.even(re(x)), ~Q.even(x)) == None
    assert ask(Q.even(im(x)), Q.even(x)) == True
    assert ask(Q.even(im(x)), Q.real(x)) == True

def test_extended_real():
    assert ask(Q.extended_real(x), Q.positive(x)) == True
    assert ask(Q.extended_real(-x), Q.positive(x)) == True
    assert ask(Q.extended_real(-x), Q.negative(x)) == True

    assert ask(Q.extended_real(x+S.Infinity), Q.real(x)) == True

def test_rational():
    assert ask(Q.rational(x), Q.integer(x)) == True
    assert ask(Q.rational(x), Q.irrational(x)) == False
    assert ask(Q.rational(x), Q.real(x)) == None
    assert ask(Q.rational(x), Q.positive(x)) == None
    assert ask(Q.rational(x), Q.negative(x)) == None
    assert ask(Q.rational(x), Q.nonzero(x)) == None

    assert ask(Q.rational(2*x), Q.rational(x)) == True
    assert ask(Q.rational(2*x), Q.integer(x)) == True
    assert ask(Q.rational(2*x), Q.even(x)) == True
    assert ask(Q.rational(2*x), Q.odd(x)) == True
    assert ask(Q.rational(2*x), Q.irrational(x)) == False

    assert ask(Q.rational(x/2), Q.rational(x)) == True
    assert ask(Q.rational(x/2), Q.integer(x)) == True
    assert ask(Q.rational(x/2), Q.even(x)) == True
    assert ask(Q.rational(x/2), Q.odd(x)) == True
    assert ask(Q.rational(x/2), Q.irrational(x)) == False

    assert ask(Q.rational(1/x), Q.rational(x)) == True
    assert ask(Q.rational(1/x), Q.integer(x)) == True
    assert ask(Q.rational(1/x), Q.even(x)) == True
    assert ask(Q.rational(1/x), Q.odd(x)) == True
    assert ask(Q.rational(1/x), Q.irrational(x)) == False

    assert ask(Q.rational(2/x), Q.rational(x)) == True
    assert ask(Q.rational(2/x), Q.integer(x)) == True
    assert ask(Q.rational(2/x), Q.even(x)) == True
    assert ask(Q.rational(2/x), Q.odd(x)) == True
    assert ask(Q.rational(2/x), Q.irrational(x)) == False

    # with multiple symbols
    assert ask(Q.rational(x*y), Q.irrational(x) & Q.irrational(y)) == None
    assert ask(Q.rational(y/x), Q.rational(x) & Q.rational(y)) == True
    assert ask(Q.rational(y/x), Q.integer(x) & Q.rational(y)) == True
    assert ask(Q.rational(y/x), Q.even(x) & Q.rational(y)) == True
    assert ask(Q.rational(y/x), Q.odd(x) & Q.rational(y)) == True
    assert ask(Q.rational(y/x), Q.irrational(x) & Q.rational(y)) == False

def test_hermitian():
    assert ask(Q.hermitian(x)) == None
    assert ask(Q.hermitian(x), Q.antihermitian(x)) == False
    assert ask(Q.hermitian(x), Q.imaginary(x)) == False
    assert ask(Q.hermitian(x), Q.prime(x)) == True
    assert ask(Q.hermitian(x), Q.real(x)) == True

    assert ask(Q.hermitian(x+1), Q.antihermitian(x)) == False
    assert ask(Q.hermitian(x+1), Q.complex(x)) == None
    assert ask(Q.hermitian(x+1), Q.hermitian(x)) == True
    assert ask(Q.hermitian(x+1), Q.imaginary(x)) == False
    assert ask(Q.hermitian(x+1), Q.real(x)) == True
    assert ask(Q.hermitian(x+I), Q.antihermitian(x)) == None
    assert ask(Q.hermitian(x+I), Q.complex(x)) == None
    assert ask(Q.hermitian(x+I), Q.hermitian(x)) == False
    assert ask(Q.hermitian(x+I), Q.imaginary(x)) == None
    assert ask(Q.hermitian(x+I), Q.real(x)) == False
    assert ask(Q.hermitian(x+y), Q.antihermitian(x) & Q.antihermitian(y)) == None
    assert ask(Q.hermitian(x+y), Q.antihermitian(x) & Q.complex(y)) == None
    assert ask(Q.hermitian(x+y), Q.antihermitian(x) & Q.hermitian(y)) == False
    assert ask(Q.hermitian(x+y), Q.antihermitian(x) & Q.imaginary(y)) == None
    assert ask(Q.hermitian(x+y), Q.antihermitian(x) & Q.real(y)) == False
    assert ask(Q.hermitian(x+y), Q.hermitian(x) & Q.complex(y)) == None
    assert ask(Q.hermitian(x+y), Q.hermitian(x) & Q.hermitian(y)) == True
    assert ask(Q.hermitian(x+y), Q.hermitian(x) & Q.imaginary(y)) == False
    assert ask(Q.hermitian(x+y), Q.hermitian(x) & Q.real(y)) == True
    assert ask(Q.hermitian(x+y), Q.imaginary(x) & Q.complex(y)) == None
    assert ask(Q.hermitian(x+y), Q.imaginary(x) & Q.imaginary(y)) == None
    assert ask(Q.hermitian(x+y), Q.imaginary(x) & Q.real(y)) == False
    assert ask(Q.hermitian(x+y), Q.real(x) & Q.complex(y)) == None
    assert ask(Q.hermitian(x+y), Q.real(x) & Q.real(y)) == True

    assert ask(Q.hermitian(I*x), Q.antihermitian(x)) == True
    assert ask(Q.hermitian(I*x), Q.complex(x)) == None
    assert ask(Q.hermitian(I*x), Q.hermitian(x)) == False
    assert ask(Q.hermitian(I*x), Q.imaginary(x)) == True
    assert ask(Q.hermitian(I*x), Q.real(x)) == False
    assert ask(Q.hermitian(x*y), Q.hermitian(x) & Q.real(y)) == True

    assert ask(Q.hermitian(x+y+z), Q.real(x) & Q.real(y) & Q.real(z)) == True
    assert ask(Q.hermitian(x+y+z), Q.real(x) & Q.real(y) & Q.imaginary(z)) == False
    assert ask(Q.hermitian(x+y+z), Q.real(x) & Q.imaginary(y) & Q.imaginary(z)) == None
    assert ask(Q.hermitian(x+y+z), Q.imaginary(x) & Q.imaginary(y) & Q.imaginary(z)) == None

    assert ask(Q.antihermitian(x)) == None
    assert ask(Q.antihermitian(x), Q.real(x)) == False
    assert ask(Q.antihermitian(x), Q.prime(x)) == False

    assert ask(Q.antihermitian(x+1), Q.antihermitian(x)) == False
    assert ask(Q.antihermitian(x+1), Q.complex(x)) == None
    assert ask(Q.antihermitian(x+1), Q.hermitian(x)) == None
    assert ask(Q.antihermitian(x+1), Q.imaginary(x)) == False
    assert ask(Q.antihermitian(x+1), Q.real(x)) == None
    assert ask(Q.antihermitian(x+I), Q.antihermitian(x)) == True
    assert ask(Q.antihermitian(x+I), Q.complex(x)) == None
    assert ask(Q.antihermitian(x+I), Q.hermitian(x)) == False
    assert ask(Q.antihermitian(x+I), Q.imaginary(x)) == True
    assert ask(Q.antihermitian(x+I), Q.real(x)) == False

    assert ask(Q.antihermitian(x+y), Q.antihermitian(x) & Q.antihermitian(y)) == True
    assert ask(Q.antihermitian(x+y), Q.antihermitian(x) & Q.complex(y)) == None
    assert ask(Q.antihermitian(x+y), Q.antihermitian(x) & Q.hermitian(y)) == False
    assert ask(Q.antihermitian(x+y), Q.antihermitian(x) & Q.imaginary(y)) == True
    assert ask(Q.antihermitian(x+y), Q.antihermitian(x) & Q.real(y)) == False
    assert ask(Q.antihermitian(x+y), Q.hermitian(x) & Q.complex(y)) == None
    assert ask(Q.antihermitian(x+y), Q.hermitian(x) & Q.hermitian(y)) == None
    assert ask(Q.antihermitian(x+y), Q.hermitian(x) & Q.imaginary(y)) == False
    assert ask(Q.antihermitian(x+y), Q.hermitian(x) & Q.real(y)) == None
    assert ask(Q.antihermitian(x+y), Q.imaginary(x) & Q.complex(y)) == None
    assert ask(Q.antihermitian(x+y), Q.imaginary(x) & Q.imaginary(y)) == True
    assert ask(Q.antihermitian(x+y), Q.imaginary(x) & Q.real(y)) == False
    assert ask(Q.antihermitian(x+y), Q.real(x) & Q.complex(y)) == None
    assert ask(Q.antihermitian(x+y), Q.real(x) & Q.real(y)) == None

    assert ask(Q.antihermitian(I*x), Q.real(x)) == True
    assert ask(Q.antihermitian(I*x), Q.antihermitian(x)) == False
    assert ask(Q.antihermitian(I*x), Q.complex(x)) == None
    assert ask(Q.antihermitian(x*y), Q.antihermitian(x) & Q.real(y)) == True

    assert ask(Q.antihermitian(x+y+z), Q.real(x) & Q.real(y) & Q.real(z)) == None
    assert ask(Q.antihermitian(x+y+z), Q.real(x) & Q.real(y) & Q.imaginary(z)) == None
    assert ask(Q.antihermitian(x+y+z), Q.real(x) & Q.imaginary(y) & Q.imaginary(z)) == False
    assert ask(Q.antihermitian(x+y+z), Q.imaginary(x) & Q.imaginary(y) & Q.imaginary(z)) == True

def test_imaginary():
    assert ask(Q.imaginary(x)) == None
    assert ask(Q.imaginary(x), Q.real(x)) == False
    assert ask(Q.imaginary(x), Q.prime(x)) == False

    assert ask(Q.imaginary(x+1), Q.real(x)) == False
    assert ask(Q.imaginary(x+1), Q.imaginary(x)) == False
    assert ask(Q.imaginary(x+I), Q.real(x)) == False
    assert ask(Q.imaginary(x+I), Q.imaginary(x)) == True
    assert ask(Q.imaginary(x+y), Q.imaginary(x) & Q.imaginary(y)) == True
    assert ask(Q.imaginary(x+y), Q.real(x) & Q.real(y)) == False
    assert ask(Q.imaginary(x+y), Q.imaginary(x) & Q.real(y)) == False
    assert ask(Q.imaginary(x+y), Q.complex(x) & Q.real(y)) == None

    assert ask(Q.imaginary(I*x), Q.real(x)) == True
    assert ask(Q.imaginary(I*x), Q.imaginary(x)) == False
    assert ask(Q.imaginary(I*x), Q.complex(x)) == None
    assert ask(Q.imaginary(x*y), Q.imaginary(x) & Q.real(y)) == True

    assert ask(Q.imaginary(x+y+z), Q.real(x) & Q.real(y) & Q.real(z)) == False
    assert ask(Q.imaginary(x+y+z), Q.real(x) & Q.real(y) & Q.imaginary(z)) == None
    assert ask(Q.imaginary(x+y+z), Q.real(x) & Q.imaginary(y) & Q.imaginary(z)) == False

def test_infinitesimal():
    assert ask(Q.infinitesimal(x)) == None
    assert ask(Q.infinitesimal(x), Q.infinitesimal(x)) == True

    assert ask(Q.infinitesimal(2*x), Q.infinitesimal(x)) == True
    assert ask(Q.infinitesimal(x*y), Q.infinitesimal(x)) == None
    assert ask(Q.infinitesimal(x*y), Q.infinitesimal(x) & Q.infinitesimal(y)) == True
    assert ask(Q.infinitesimal(x*y), Q.infinitesimal(x) & Q.bounded(y)) == True

    assert ask(Q.infinitesimal(x**2), Q.infinitesimal(x)) == True

def test_integer():
    assert ask(Q.integer(x)) == None
    assert ask(Q.integer(x), Q.integer(x)) == True
    assert ask(Q.integer(x), ~Q.integer(x)) == False
    assert ask(Q.integer(x), ~Q.real(x)) == False
    assert ask(Q.integer(x), ~Q.positive(x)) == None
    assert ask(Q.integer(x), Q.even(x) | Q.odd(x)) == True

    assert ask(Q.integer(2*x), Q.integer(x)) == True
    assert ask(Q.integer(2*x), Q.even(x)) == True
    assert ask(Q.integer(2*x), Q.prime(x)) == True
    assert ask(Q.integer(2*x), Q.rational(x)) == None
    assert ask(Q.integer(2*x), Q.real(x)) == None
    assert ask(Q.integer(sqrt(2)*x), Q.integer(x)) == False

    assert ask(Q.integer(x/2), Q.odd(x)) == False
    assert ask(Q.integer(x/2), Q.even(x)) == True
    assert ask(Q.integer(x/3), Q.odd(x)) == None
    assert ask(Q.integer(x/3), Q.even(x)) == None

def test_negative():
    assert ask(Q.negative(x), Q.negative(x)) == True
    assert ask(Q.negative(x), Q.positive(x)) == False
    assert ask(Q.negative(x), ~Q.real(x)) == False
    assert ask(Q.negative(x), Q.prime(x)) == False
    assert ask(Q.negative(x), ~Q.prime(x)) == None

    assert ask(Q.negative(-x), Q.positive(x)) == True
    assert ask(Q.negative(-x), ~Q.positive(x)) == None
    assert ask(Q.negative(-x), Q.negative(x)) == False
    assert ask(Q.negative(-x), Q.positive(x)) == True

    assert ask(Q.negative(x-1), Q.negative(x)) == True
    assert ask(Q.negative(x+y)) == None
    assert ask(Q.negative(x+y), Q.negative(x)) == None
    assert ask(Q.negative(x+y), Q.negative(x) & Q.negative(y)) == True

    assert ask(Q.negative(x**2)) == None
    assert ask(Q.negative(x**2), Q.real(x)) == False
    assert ask(Q.negative(x**1.4), Q.real(x)) == None

    assert ask(Q.negative(x*y)) == None
    assert ask(Q.negative(x*y), Q.positive(x) & Q.positive(y)) == False
    assert ask(Q.negative(x*y), Q.positive(x) & Q.negative(y)) == True
    assert ask(Q.negative(x*y), Q.complex(x) & Q.complex(y)) == None

    assert ask(Q.negative(x**y)) == None
    assert ask(Q.negative(x**y), Q.negative(x) & Q.even(y)) == False
    assert ask(Q.negative(x**y), Q.negative(x) & Q.odd(y)) == True
    assert ask(Q.negative(x**y), Q.positive(x) & Q.integer(y)) == False

    assert ask(Q.negative(Abs(x))) == False

def test_nonzero():
    assert ask(Q.nonzero(x)) == None
    assert ask(Q.nonzero(x), Q.real(x)) == None
    assert ask(Q.nonzero(x), Q.positive(x)) == True
    assert ask(Q.nonzero(x), Q.negative(x)) == True
    assert ask(Q.nonzero(x), Q.negative(x) | Q.positive(x)) == True

    assert ask(Q.nonzero(x+y)) == None
    assert ask(Q.nonzero(x+y), Q.positive(x) & Q.positive(y)) == True
    assert ask(Q.nonzero(x+y), Q.positive(x) & Q.negative(y)) == None
    assert ask(Q.nonzero(x+y), Q.negative(x) & Q.negative(y)) == True

    assert ask(Q.nonzero(2*x)) == None
    assert ask(Q.nonzero(2*x), Q.positive(x)) == True
    assert ask(Q.nonzero(2*x), Q.negative(x)) == True
    assert ask(Q.nonzero(x*y), Q.nonzero(x)) == None
    assert ask(Q.nonzero(x*y), Q.nonzero(x) & Q.nonzero(y)) == True

    assert ask(Q.nonzero(Abs(x))) == None
    assert ask(Q.nonzero(Abs(x)), Q.nonzero(x)) == True

def test_odd():
    assert ask(Q.odd(x)) == None
    assert ask(Q.odd(x), Q.odd(x)) == True
    assert ask(Q.odd(x), Q.integer(x)) == None
    assert ask(Q.odd(x), ~Q.integer(x)) == False
    assert ask(Q.odd(x), Q.rational(x)) == None
    assert ask(Q.odd(x), Q.positive(x)) == None

    assert ask(Q.odd(-x), Q.odd(x)) == True

    assert ask(Q.odd(2*x)) == None
    assert ask(Q.odd(2*x), Q.integer(x)) == False
    assert ask(Q.odd(2*x), Q.odd(x)) == False
    assert ask(Q.odd(2*x), Q.irrational(x)) == False
    assert ask(Q.odd(2*x), ~Q.integer(x)) == None
    assert ask(Q.odd(3*x), Q.integer(x)) == None

    assert ask(Q.odd(x/3), Q.odd(x)) == None
    assert ask(Q.odd(x/3), Q.even(x)) == None

    assert ask(Q.odd(x+1), Q.even(x)) == True
    assert ask(Q.odd(x+2), Q.even(x)) == False
    assert ask(Q.odd(x+2), Q.odd(x))  == True
    assert ask(Q.odd(3-x), Q.odd(x))  == False
    assert ask(Q.odd(3-x), Q.even(x))  == True
    assert ask(Q.odd(3+x), Q.odd(x))  == False
    assert ask(Q.odd(3+x), Q.even(x))  == True
    assert ask(Q.odd(x+y), Q.odd(x) & Q.odd(y)) == False
    assert ask(Q.odd(x+y), Q.odd(x) & Q.even(y)) == True
    assert ask(Q.odd(x-y), Q.even(x) & Q.odd(y)) == True
    assert ask(Q.odd(x-y), Q.odd(x) & Q.odd(y)) == False

    assert ask(Q.odd(x+y+z), Q.odd(x) & Q.odd(y) & Q.even(z)) == False
    assert ask(Q.odd(x+y+z+t), Q.odd(x) & Q.odd(y) & Q.even(z) & Q.integer(t)) == None

    assert ask(Q.odd(2*x + 1), Q.integer(x)) == True
    assert ask(Q.odd(2*x + y), Q.integer(x) & Q.odd(y)) == True
    assert ask(Q.odd(2*x + y), Q.integer(x) & Q.even(y)) == False
    assert ask(Q.odd(2*x + y), Q.integer(x) & Q.integer(y)) == None
    assert ask(Q.odd(x*y), Q.odd(x) & Q.even(y)) == False
    assert ask(Q.odd(x*y), Q.odd(x) & Q.odd(y)) == True
    assert ask(Q.odd(2*x*y), Q.rational(x) & Q.rational(x)) == None
    assert ask(Q.odd(2*x*y), Q.irrational(x) & Q.irrational(x)) == None

    assert ask(Q.odd(Abs(x)), Q.odd(x)) == True

def test_prime():
    assert ask(Q.prime(x), Q.prime(x)) == True
    assert ask(Q.prime(x), ~Q.prime(x)) == False
    assert ask(Q.prime(x), Q.integer(x)) == None
    assert ask(Q.prime(x), ~Q.integer(x)) == False

    assert ask(Q.prime(2*x), Q.integer(x)) == False
    assert ask(Q.prime(x*y)) == None
    assert ask(Q.prime(x*y), Q.prime(x)) == None
    assert ask(Q.prime(x*y), Q.integer(x) & Q.integer(y)) == False

    assert ask(Q.prime(x**2), Q.integer(x)) == False
    assert ask(Q.prime(x**2), Q.prime(x)) == False
    assert ask(Q.prime(x**y), Q.integer(x) & Q.integer(y)) == False

def test_positive():
    assert ask(Q.positive(x), Q.positive(x)) == True
    assert ask(Q.positive(x), Q.negative(x)) == False
    assert ask(Q.positive(x), Q.nonzero(x)) == None

    assert ask(Q.positive(-x), Q.positive(x)) == False
    assert ask(Q.positive(-x), Q.negative(x)) == True

    assert ask(Q.positive(x+y), Q.positive(x) & Q.positive(y)) == True
    assert ask(Q.positive(x+y), Q.positive(x) & Q.negative(y)) == None

    assert ask(Q.positive(2*x), Q.positive(x)) == True
    assumptions =  Q.positive(x) & Q.negative(y) & Q.negative(z) & Q.positive(w)
    assert ask(Q.positive(x*y*z))  == None
    assert ask(Q.positive(x*y*z), assumptions) == True
    assert ask(Q.positive(-x*y*z), assumptions) == False

    assert ask(Q.positive(x**2), Q.positive(x)) == True
    assert ask(Q.positive(x**2), Q.negative(x)) == True

    #exponential
    assert ask(Q.positive(exp(x)), Q.real(x)) == True
    assert ask(Q.positive(x + exp(x)), Q.real(x)) == None

    #absolute value
    assert ask(Q.positive(Abs(x))) == None # Abs(0) = 0
    assert ask(Q.positive(Abs(x)), Q.positive(x)) == True

@XFAIL
def test_positive_xfail():
    assert ask(Q.positive(1/(1 + x**2)), Q.real(x)) == True

def test_real():
    assert ask(Q.real(x)) == None
    assert ask(Q.real(x), Q.real(x)) == True
    assert ask(Q.real(x), Q.nonzero(x)) == True
    assert ask(Q.real(x), Q.positive(x)) == True
    assert ask(Q.real(x), Q.negative(x)) == True
    assert ask(Q.real(x), Q.integer(x)) == True
    assert ask(Q.real(x), Q.even(x)) == True
    assert ask(Q.real(x), Q.prime(x)) == True

    assert ask(Q.real(x/sqrt(2)), Q.real(x)) == True
    assert ask(Q.real(x/sqrt(-2)), Q.real(x)) == False

    assert ask(Q.real(x+1), Q.real(x)) == True
    assert ask(Q.real(x+I), Q.real(x)) == False
    assert ask(Q.real(x+I), Q.complex(x)) == None

    assert ask(Q.real(2*x), Q.real(x)) == True
    assert ask(Q.real(I*x), Q.real(x)) == False
    assert ask(Q.real(I*x), Q.imaginary(x)) == True
    assert ask(Q.real(I*x), Q.complex(x)) == None

    assert ask(Q.real(x**2), Q.real(x)) == True
    assert ask(Q.real(sqrt(x)), Q.negative(x)) == False
    assert ask(Q.real(x**y), Q.real(x) & Q.integer(y)) == True
    assert ask(Q.real(x**y), Q.real(x) & Q.real(y)) == None
    assert ask(Q.real(x**y), Q.positive(x) & Q.real(y)) == True

    # trigonometric functions
    assert ask(Q.real(sin(x))) == None
    assert ask(Q.real(cos(x))) == None
    assert ask(Q.real(sin(x)), Q.real(x)) == True
    assert ask(Q.real(cos(x)), Q.real(x)) == True

    # exponential function
    assert ask(Q.real(exp(x))) == None
    assert ask(Q.real(exp(x)), Q.real(x)) == True
    assert ask(Q.real(x + exp(x)), Q.real(x)) == True

    # Q.complexes
    assert ask(Q.real(re(x))) == True
    assert ask(Q.real(im(x))) == True

def test_algebraic():
    assert ask(Q.algebraic(x)) == None

    assert ask(Q.algebraic(I)) == True
    assert ask(Q.algebraic(2*I)) == True
    assert ask(Q.algebraic(I/3)) == True

    assert ask(Q.algebraic(sqrt(7))) == True
    assert ask(Q.algebraic(2*sqrt(7))) == True
    assert ask(Q.algebraic(sqrt(7)/3)) == True

    assert ask(Q.algebraic(I*sqrt(3))) == True
    assert ask(Q.algebraic(sqrt(1+I*sqrt(3)))) == True

    assert ask(Q.algebraic((1+I*sqrt(3)**(S(17)/31)))) == True
    assert ask(Q.algebraic((1+I*sqrt(3)**(S(17)/pi)))) == False

    assert ask(Q.algebraic(sin(7))) == None
    assert ask(Q.algebraic(sqrt(sin(7)))) == None
    assert ask(Q.algebraic(sqrt(y+I*sqrt(7)))) == None

    assert ask(Q.algebraic(oo)) == False
    assert ask(Q.algebraic(-oo)) == False

    assert ask(Q.algebraic(2.47)) == False

def test_global():
    """Test ask with global assumptions"""
    assert ask(Q.integer(x)) == None
    global_assumptions.add(Q.integer(x))
    assert ask(Q.integer(x)) == True
    global_assumptions.clear()
    assert ask(Q.integer(x)) == None

def test_custom_context():
    """Test ask with custom assumptions context"""
    assert ask(Q.integer(x)) == None
    local_context = AssumptionsContext()
    local_context.add(Q.integer(x))
    assert ask(Q.integer(x), context = local_context) == True
    assert ask(Q.integer(x)) == None

def test_functions_in_assumptions():
    assert ask(Q.negative(x), Q.real(x) >> Q.positive(x)) is False
    assert ask(Q.negative(x), Equivalent(Q.real(x), Q.positive(x))) is False
    assert ask(Q.negative(x), Xor(Q.real(x), Q.negative(x))) is False

def test_composite_ask():
    assert ask(Q.negative(x) & Q.integer(x),
           assumptions=Q.real(x) >> Q.positive(x)) is False

def test_composite_proposition():
    assert ask(True) is True
    assert ask(~Q.negative(x), Q.positive(x)) is True
    assert ask(~Q.real(x), Q.commutative(x)) is None
    assert ask(Q.negative(x) & Q.integer(x), Q.positive(x)) is False
    assert ask(Q.negative(x) & Q.integer(x)) is None
    assert ask(Q.real(x) | Q.integer(x), Q.positive(x)) is True
    assert ask(Q.real(x) | Q.integer(x)) is None
    assert ask(Q.real(x) >> Q.positive(x), Q.negative(x)) is False
    assert ask(Implies(Q.real(x), Q.positive(x), evaluate=False), Q.negative(x)) is False
    assert ask(Implies(Q.real(x), Q.positive(x), evaluate=False)) is None
    assert ask(Equivalent(Q.integer(x), Q.even(x)), Q.even(x)) is True
    assert ask(Equivalent(Q.integer(x), Q.even(x))) is None
    assert ask(Equivalent(Q.positive(x), Q.integer(x)), Q.integer(x)) is None

def test_incompatible_resolutors():
    class Prime2AskHandler(AskHandler):
        @staticmethod
        def Number(expr, assumptions):
            return True
    register_handler('prime', Prime2AskHandler)
    raises(ValueError, lambda: ask(Q.prime(4)))
    remove_handler('prime', Prime2AskHandler)

    class InconclusiveHandler(AskHandler):
        @staticmethod
        def Number(expr, assumptions):
            return None
    register_handler('prime', InconclusiveHandler)
    assert ask(Q.prime(3)) == True

def test_key_extensibility():
    """test that you can add keys to the ask system at runtime"""
    # make sure the key is not defined
    raises(AttributeError, lambda: ask(Q.my_key(x)))
    class MyAskHandler(AskHandler):
        @staticmethod
        def Symbol(expr, assumptions):
            return True
    register_handler('my_key', MyAskHandler)
    assert ask(Q.my_key(x)) == True
    assert ask(Q.my_key(x+1)) == None
    remove_handler('my_key', MyAskHandler)
    del Q.my_key
    raises(AttributeError, lambda: ask(Q.my_key(x)))

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
    assert ask(Q.prime(a)) == True

def test_compute_known_facts():
    ns = {}
    exec 'from sympy.logic.boolalg import And, Or, Not' in globals(), ns
    exec compute_known_facts() in globals(), ns
    assert ns['known_facts_cnf'] == known_facts_cnf
    assert ns['known_facts_dict'] == known_facts_dict

def test_Add_queries():
    assert ask(Q.prime(12345678901234567890 + (cos(1)**2 + sin(1)**2))) is True
    assert ask(Q.even(Add(S(2), S(2), evaluate=0))) is True
    assert ask(Q.prime(Add(S(2), S(2), evaluate=0))) is False
    assert ask(Q.integer(Add(S(2), S(2), evaluate=0))) is True
