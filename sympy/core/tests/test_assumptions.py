from sympy.core import Symbol, S, Rational, Integer
from sympy.utilities.pytest import XFAIL, raises

def test_symbol_unset():
    x = Symbol('x',real=True, integer=True)
    assert x.is_real == True
    assert x.is_integer == True
    assert x.is_imaginary == False
    assert x.is_noninteger == False

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
    assert z.is_composite == False

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
    assert z.is_composite == False

def test_infinity():
    oo = S.Infinity

    assert oo.is_commutative    == True
    assert oo.is_integer        == None
    assert oo.is_rational       == None
    assert oo.is_real           == True
    assert oo.is_complex        == True
    assert oo.is_noninteger     == None
    assert oo.is_irrational     == None
    assert oo.is_imaginary      == False
    assert oo.is_positive       == True
    assert oo.is_negative       == False
    assert oo.is_nonpositive    == False
    assert oo.is_nonnegative    == True
    assert oo.is_even           == None
    assert oo.is_odd            == None
    assert oo.is_bounded        == False
    assert oo.is_unbounded      == True
    assert oo.is_finite         == False
    assert oo.is_infinitesimal  == False
    assert oo.is_comparable     == True
    assert oo.is_prime          == None
    assert oo.is_composite      == None

def test_neg_infinity():
    mm = S.NegativeInfinity

    assert mm.is_commutative    == True
    assert mm.is_integer        == None
    assert mm.is_rational       == None
    assert mm.is_real           == True
    assert mm.is_complex        == True
    assert mm.is_noninteger     == None
    assert mm.is_irrational     == None
    assert mm.is_imaginary      == False
    assert mm.is_positive       == False
    assert mm.is_negative       == True
    assert mm.is_nonpositive    == True
    assert mm.is_nonnegative    == False
    assert mm.is_even           == None
    assert mm.is_odd            == None
    assert mm.is_bounded        == False
    assert mm.is_unbounded      == True
    assert mm.is_finite         == False
    assert mm.is_infinitesimal  == False
    assert mm.is_comparable     == True
    assert mm.is_prime          == False
    assert mm.is_composite      == False

def test_nan():
    nan = S.NaN

    assert nan.is_commutative   == True
    assert nan.is_integer       == None
    assert nan.is_rational      == None
    assert nan.is_real          == None
    assert nan.is_complex       == None
    assert nan.is_noninteger    == None
    assert nan.is_irrational    == None
    assert nan.is_imaginary     == None
    assert nan.is_positive      == None
    assert nan.is_negative      == None
    assert nan.is_nonpositive   == None
    assert nan.is_nonnegative   == None
    assert nan.is_even          == None
    assert nan.is_odd           == None
    assert nan.is_bounded       == None
    assert nan.is_unbounded     == None
    assert nan.is_finite        == None
    assert nan.is_infinitesimal == None
    assert nan.is_comparable    == False
    assert nan.is_prime         == None
    assert nan.is_composite     == None

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
    assert r.is_positive == True
    assert r.is_negative == False
    assert r.is_nonpositive == False
    assert r.is_nonnegative == True
    assert r.is_even == False
    assert r.is_odd == False
    assert r.is_bounded == True
    assert r.is_unbounded == False
    assert r.is_finite == True
    assert r.is_infinitesimal == False
    assert r.is_comparable == True
    assert r.is_prime == False
    assert r.is_composite == False

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
    z = S.Pi
    assert z.is_commutative == True
    assert z.is_integer == False
    assert z.is_rational == False
    assert z.is_real == True
    assert z.is_complex == True
    assert z.is_noninteger == True
    assert z.is_irrational == True
    assert z.is_imaginary == False
    assert z.is_positive == True
    assert z.is_negative == False
    assert z.is_nonpositive == False
    assert z.is_nonnegative == True
    assert z.is_even == False
    assert z.is_odd == False
    assert z.is_bounded == True
    assert z.is_unbounded == False
    assert z.is_finite == True
    assert z.is_infinitesimal == False
    assert z.is_comparable == True
    assert z.is_prime == False
    assert z.is_composite == False

def test_E():
    z = S.Exp1
    assert z.is_commutative == True
    assert z.is_integer == False
    assert z.is_rational == False
    assert z.is_real == True
    assert z.is_complex == True
    assert z.is_noninteger == True
    assert z.is_irrational == True
    assert z.is_imaginary == False
    assert z.is_positive == True
    assert z.is_negative == False
    assert z.is_nonpositive == False
    assert z.is_nonnegative == True
    assert z.is_even == False
    assert z.is_odd == False
    assert z.is_bounded == True
    assert z.is_unbounded == False
    assert z.is_finite == True
    assert z.is_infinitesimal == False
    assert z.is_comparable == True
    assert z.is_prime == False
    assert z.is_composite == False

def test_I():
    z = S.ImaginaryUnit
    assert z.is_commutative == True
    assert z.is_integer == False
    assert z.is_rational == False
    assert z.is_real == False
    assert z.is_complex == True
    assert z.is_noninteger == False
    assert z.is_irrational == False
    assert z.is_imaginary == True
    assert z.is_positive == False
    assert z.is_negative == False
    assert z.is_nonpositive == False
    assert z.is_nonnegative == False
    assert z.is_even == False
    assert z.is_odd == False
    assert z.is_bounded == True
    assert z.is_unbounded == False
    assert z.is_finite == True
    assert z.is_infinitesimal == False
    assert z.is_comparable == None
    assert z.is_prime == False
    assert z.is_composite == False

def test_symbol_zero():
    x = Symbol('x',zero=True)
    assert x.is_positive == False
    assert x.is_nonpositive == True
    assert x.is_negative == False
    assert x.is_nonnegative == True
    assert x.is_zero == True
    assert x.is_nonzero == False

def test_symbol_positive():
    x = Symbol('x',positive=True)
    assert x.is_positive == True
    assert x.is_nonpositive == False
    assert x.is_negative == False
    assert x.is_nonnegative == True
    assert x.is_zero == False
    assert x.is_nonzero == True

def test_neg_symbol_positive():
    x = -Symbol('x',positive=True)
    assert x.is_positive == False
    assert x.is_nonpositive == True
    assert x.is_negative == True
    assert x.is_nonnegative == False

def test_neg_symbol_positive2():
    x = -Symbol('x',positive=True)
    assert x.is_zero == False
    assert x.is_nonzero == True

def test_symbol_nonpositive():
    x = Symbol('x',nonpositive=True)
    assert x.is_positive == False
    assert x.is_nonpositive == True
    assert x.is_negative == None
    assert x.is_nonnegative == None
    assert x.is_zero == None
    assert x.is_nonzero == None

def test_neg_symbol_nonpositive():
    x = -Symbol('x',nonpositive=True)
    assert x.is_positive == None
    assert x.is_nonpositive == None
    assert x.is_negative == False
    assert x.is_nonnegative == True
    assert x.is_zero == None
    assert x.is_nonzero == None

def test_prime_symbol():
    x = Symbol('x', prime=True)
    assert x.is_prime == True
    assert x.is_integer == True
    assert x.is_positive == True
    assert x.is_negative == False
    assert x.is_nonpositive == False
    assert x.is_nonnegative == True

    x = Symbol('x', prime=False)
    assert x.is_prime == False
    assert x.is_integer == None
    assert x.is_positive == None
    assert x.is_negative == None
    assert x.is_nonpositive == None
    assert x.is_nonnegative == None

def test_other_symbol():
    x = Symbol('x', integer=True)
    assert x.is_integer == True
    assert x.is_real == True

    x = Symbol('x', integer=True, nonnegative=True)
    assert x.is_integer     == True
    assert x.is_nonnegative == True
    assert x.is_negative    == False
    assert x.is_positive    == None

    x = Symbol('x', integer=True, nonpositive=True)
    assert x.is_integer     == True
    assert x.is_nonpositive == True
    assert x.is_positive    == False
    assert x.is_negative    == None

    x = Symbol('x', odd=True)
    assert x.is_odd == True
    assert x.is_even == False
    assert x.is_integer == True

    x = Symbol('x', even=True)
    assert x.is_even == True
    assert x.is_odd == False
    assert x.is_integer == True

    x = Symbol('x', integer=True, nonnegative=True)
    assert x.is_integer == True
    assert x.is_nonnegative == True

    x = Symbol('x', integer=True, nonpositive=True)
    assert x.is_integer == True
    assert x.is_nonpositive == True

    raises(AttributeError, "x.is_real = False")


def test_other_symbol_fail1():
    # XXX x.is_even currently will be True
    x = Symbol('x', odd=False)
    assert x.is_odd == False
    assert x.is_even == None
    assert x.is_integer == None

def test_other_symbol_fail2():
    # XXX x.is_odd currently will be True
    x = Symbol('x', even=False)
    assert x.is_even == False
    assert x.is_odd == None
    assert x.is_integer == None


def test_issue749():
    a = Symbol('a', real=False)

    assert a.is_real        == False
    assert a.is_integer     == False
    assert a.is_negative    == False
    assert a.is_positive    == False
    assert a.is_nonnegative == False
    assert a.is_nonpositive == False
    assert a.is_zero        == False


def test_issue726():
    """catch: hash instability"""
    x = Symbol("x")
    y = Symbol("y")
    a1 = x+y
    a2 = y+x
    a2.is_comparable

    h1 = hash(a1)
    h2 = hash(a2)
    assert h1 == h2


def test_hash_vs_typeinfo():
    """seemingly different typeinfo, but in fact equal"""
    x = Symbol('x')

    # the following two are semantically equal
    x1= Symbol('x', even=True)
    x2= Symbol('x', integer=True, odd=False)

    assert hash(x1) == hash(x2)
    assert x1 == x2


@XFAIL  # our __eq__ just does not look at it (yet)
def test_hash_vs_typeinfo_2():
    """different typeinfo should mean !eq"""
    # the following two are semantically different
    x = Symbol('x')
    x1= Symbol('x', even=True)

    assert x != x1


def test_hash_vs_eq():
    """catch: different hash for equal objects"""
    a = 1 + S.Pi    # important: do not fold it into a Number instance
    ha= hash(a)     #            it should be Add/Mul/... to trigger the bug

    a.is_positive   # this uses .evalf() and deduces it is positive
    assert a.is_positive == True

    # be sure that hash stayed the same
    assert ha == hash(a)

    # now b should be the same expression
    b = a.expand(trig=True)
    hb= hash(b)

    assert a == b
    assert ha== hb

