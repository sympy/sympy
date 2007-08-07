import py
from sympy.core import *


def test_symbol_unset():
    x = Symbol('x',real=True, integer=True)
    assert x.is_real == True
    assert x.is_integer == True
    assert x.is_imaginary == False
    assert x.is_noninteger == False
    x.assume(real=None)
    assert x.is_real == None
    assert x.is_integer == None

def test_zero():
    z = Integer(0)
    assert z.is_commutative == True
    assert z.is_integer == True
    assert z.is_rational == True
    assert z.is_real == True
    assert z.is_complex == True
    assert z.is_noninteger == False
    assert z.is_irrational == False
    assert z.is_imaginary == False
    assert z.is_noncomplex == False
    assert z.is_positive == False
    assert z.is_negative == False
    assert z.is_nonpositive == True
    assert z.is_nonnegative == True
    assert z.is_even == True
    assert z.is_odd == False
    assert z.is_bounded == True
    assert z.is_unbounded == False
    assert z.is_finite == False
    assert z.is_infinitesimal == True
    assert z.is_comparable == True
    assert z.is_prime == False
    assert z.is_composite == True

def test_one():
    z = Integer(1)
    assert z.is_commutative == True
    assert z.is_integer == True
    assert z.is_rational == True
    assert z.is_real == True
    assert z.is_complex == True
    assert z.is_noninteger == False
    assert z.is_irrational == False
    assert z.is_imaginary == False
    assert z.is_noncomplex == False
    assert z.is_positive == True
    assert z.is_negative == False
    assert z.is_nonpositive == False
    assert z.is_nonnegative == True
    assert z.is_even == False
    assert z.is_odd == True
    assert z.is_bounded == True
    assert z.is_unbounded == False
    assert z.is_finite == True
    assert z.is_infinitesimal == False
    assert z.is_comparable == True
    assert z.is_prime == True
    assert z.is_composite == False

def test_negativeone():
    z = Integer(-1)
    assert z.is_commutative == True
    assert z.is_integer == True
    assert z.is_rational == True
    assert z.is_real == True
    assert z.is_complex == True
    assert z.is_noninteger == False
    assert z.is_irrational == False
    assert z.is_imaginary == False
    assert z.is_noncomplex == False
    assert z.is_positive == False
    assert z.is_negative == True
    assert z.is_nonpositive == True
    assert z.is_nonnegative == False
    assert z.is_even == False
    assert z.is_odd == True
    assert z.is_bounded == True
    assert z.is_unbounded == False
    assert z.is_finite == True
    assert z.is_infinitesimal == False
    assert z.is_comparable == True
    assert z.is_prime == False
    assert z.is_composite == True

def test_pos_rational():
    r = Rational(3,4)
    assert r.is_commutative == True
    assert r.is_integer == False
    assert r.is_rational == True
    assert r.is_real == True
    assert r.is_complex == True
    assert r.is_noninteger == True
    assert r.is_irrational == False
    assert r.is_imaginary == False
    assert r.is_noncomplex == False
    assert r.is_positive == True
    assert r.is_negative == False
    assert r.is_nonpositive == False
    assert r.is_nonnegative == True
    assert r.is_even == None
    assert r.is_odd == None
    assert r.is_bounded == True
    assert r.is_unbounded == False
    assert r.is_finite == True
    assert r.is_infinitesimal == False
    assert r.is_comparable == True
    assert r.is_prime == None
    assert r.is_composite == None

    r = Rational(1,4)
    assert r.is_nonpositive == False
    assert r.is_positive == True
    assert r.is_negative == False
    assert r.is_nonnegative == True
    r = Rational(5,4)
    assert r.is_negative == False
    assert r.is_positive == True
    assert r.is_nonpositive == False
    assert r.is_nonnegative == True
    r = Rational(5,3)
    assert r.is_nonnegative == True
    assert r.is_positive == True
    assert r.is_negative == False
    assert r.is_nonpositive == False

def test_neg_rational():
    r = Rational(-3,4)
    assert r.is_positive == False
    assert r.is_nonpositive == True
    assert r.is_negative == True
    assert r.is_nonnegative == False
    r = Rational(-1,4)
    assert r.is_nonpositive == True
    assert r.is_positive == False
    assert r.is_negative == True
    assert r.is_nonnegative == False
    r = Rational(-5,4)
    assert r.is_negative == True
    assert r.is_positive == False
    assert r.is_nonpositive == True
    assert r.is_nonnegative == False
    r = Rational(-5,3)
    assert r.is_nonnegative == False
    assert r.is_positive == False
    assert r.is_negative == True
    assert r.is_nonpositive == True

def test_pi():
    z = Pi
    assert z.is_commutative == True
    assert z.is_integer == None
    assert z.is_rational == False
    assert z.is_real == True
    assert z.is_complex == True
    assert z.is_noninteger == None
    assert z.is_irrational == True
    assert z.is_imaginary == False
    assert z.is_noncomplex == False
    assert z.is_positive == True
    assert z.is_negative == False
    assert z.is_nonpositive == False
    assert z.is_nonnegative == True
    assert z.is_even == None
    assert z.is_odd == None
    assert z.is_bounded == True
    assert z.is_unbounded == False
    assert z.is_finite == True
    assert z.is_infinitesimal == False
    assert z.is_comparable == True
    assert z.is_prime == None
    assert z.is_composite == None

def test_E():
    z = E
    assert z.is_commutative == True
    assert z.is_integer == None
    assert z.is_rational == False
    assert z.is_real == True
    assert z.is_complex == True
    assert z.is_noninteger == None
    assert z.is_irrational == True
    assert z.is_imaginary == False
    assert z.is_noncomplex == False
    assert z.is_positive == True
    assert z.is_negative == False
    assert z.is_nonpositive == False
    assert z.is_nonnegative == True
    assert z.is_even == None
    assert z.is_odd == None
    assert z.is_bounded == True
    assert z.is_unbounded == False
    assert z.is_finite == True
    assert z.is_infinitesimal == False
    assert z.is_comparable == True
    assert z.is_prime == None
    assert z.is_composite == None

def test_I():
    z = I
    assert z.is_commutative == True
    assert z.is_integer == None
    assert z.is_rational == None
    assert z.is_real == False
    assert z.is_complex == True
    assert z.is_noninteger == None
    assert z.is_irrational == None
    assert z.is_imaginary == True
    assert z.is_noncomplex == False
    assert z.is_positive == None
    assert z.is_negative == None
    assert z.is_nonpositive == None
    assert z.is_nonnegative == None
    assert z.is_even == None
    assert z.is_odd == None
    assert z.is_bounded == True
    assert z.is_unbounded == False
    assert z.is_finite == True
    assert z.is_infinitesimal == False
    assert z.is_comparable == None
    assert z.is_prime == None
    assert z.is_composite == None

def test_symbol_positive():
    x = Basic.Symbol('x',positive=True)
    assert x.is_positive == True
    assert x.is_nonpositive == False
    assert x.is_negative == False
    assert x.is_nonnegative == True

def test_neg_symbol_positive():
    x = -Basic.Symbol('x',positive=True)
    assert x.is_positive == False
    assert x.is_nonpositive == True
    assert x.is_negative == True
    assert x.is_nonnegative == False

def test_symbol_nonpositive():
    x = Basic.Symbol('x',nonpositive=True)
    assert x.is_positive == False
    assert x.is_nonpositive == True
    assert x.is_negative == None
    assert x.is_nonnegative == None

def test_neg_symbol_nonpositive():
    x = -Basic.Symbol('x',nonpositive=True)
    assert x.is_positive == None
    assert x.is_nonpositive == None
    assert x.is_negative == False
    assert x.is_nonnegative == True

def test_prime_symbol():
    x = Basic.Symbol('x', prime=True)
    assert x.is_prime == True
    assert x.is_integer == True
    assert x.is_positive == True
    assert x.is_negative == False
    assert x.is_nonpositive == False
    assert x.is_nonnegative == True

    x = Basic.Symbol('x', prime=False)
    assert x.is_prime == False
    assert x.is_integer == None
    assert x.is_positive == None
    assert x.is_negative == None
    assert x.is_nonpositive == None
    assert x.is_nonnegative == None

def _test_other_symbol():
    x = Basic.Symbol('x', integer=True)
    assert x.is_integer == True
    assert x.is_real == True

    x = Basic.Symbol('x', integer=True, nonnegative=True)
    assert x.is_nni == True
    assert x.is_ni == False
    assert x.is_pi == None

    x = Basic.Symbol('x', integer=True, nonpositive=True)
    assert x.is_npi == True
    assert x.is_pi == False
    assert x.is_ni == None

    x = Basic.Symbol('x', odd=True)
    assert x.is_odd == True
    assert x.is_even == False
    assert x.is_integer == True

    # XXX x.is_even currently will be True
    x = Basic.Symbol('x', odd=False)
    assert x.is_odd == False
    assert x.is_even == None
    assert x.is_integer == None

    x = Basic.Symbol('x', even=True)
    assert x.is_even == True
    assert x.is_odd == False
    assert x.is_integer == True

    # XXX x.is_odd currently will be True
    x = Basic.Symbol('x', even=False)
    assert x.is_even == False
    assert x.is_odd == None
    assert x.is_integer == None

    x = Basic.Symbol('x', nni=True)
    assert x.is_integer == True
    assert x.is_nonnegative == True

    x = Basic.Symbol('x', npi=True)
    assert x.is_integer == True
    assert x.is_nonpositive == True

    py.test.raises(AttributeError, "x.is_real = False")