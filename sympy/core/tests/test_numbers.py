from sympy import Rational, Symbol, Real, I, sqrt, oo, nan, pi, E, Integer, \
        S, factorial, Catalan, EulerGamma, GoldenRatio, cos
from sympy.core.power import integer_nthroot

from sympy.core.numbers import igcd, ilcm, igcdex, seterr, _intcache
from sympy.utilities.pytest import raises
from sympy import mpmath

def test_integers_cache():
    python_int = 2**65 + 3175259

    while python_int in _intcache or hash(python_int) in _intcache:
        value += 1

    sympy_int = Integer(python_int)

    assert python_int in _intcache
    assert hash(python_int) not in _intcache
    assert sympy_int not in _intcache

    sympy_int_int = Integer(sympy_int)

    assert python_int in _intcache
    assert hash(python_int) not in _intcache
    assert sympy_int_int not in _intcache

    sympy_hash_int = Integer(hash(python_int))

    assert python_int in _intcache
    assert hash(python_int) in _intcache
    assert sympy_hash_int not in _intcache

def test_seterr():
    seterr(divide = True)
    raises(ValueError,"S.Zero/S.Zero")
    seterr(divide = False)
    S.Zero / S.Zero == S.NaN

def test_mod():
    x = Rational(1, 2)
    y = Rational(3, 4)
    z = Rational(5, 18043)

    assert x % x == 0
    assert x % y == 1/S(2)
    assert x % z == 3/S(36086)
    assert y % x == 1/S(4)
    assert y % y == 0
    assert y % z == 9/S(72172)
    assert z % x == 5/S(18043)
    assert z % y == 5/S(18043)
    assert z % z == 0

    a = Real('2.6')

    assert round(a % Real('0.2'), 15) == 0.2
    assert round(a % 2, 15) == 0.6
    assert round(a % 0.5, 15) == 0.1
    assert Rational(3,4) % Real(1.1) == 0.75

    a = Integer(7)
    b = Integer(4)

    assert type(a % b) == Integer

def test_divmod():
    assert divmod(S(12), S(8)) == (1, 4)
    assert divmod(-S(12), S(8)) == (-2, 4)
    assert divmod(S(0), S(1)) == (0, 0)
    raises(ZeroDivisionError, "divmod(S(0), S(0))")
    raises(ZeroDivisionError, "divmod(S(1), S(0))")

def test_igcd():
    assert igcd(0, 0) == 0
    assert igcd(0, 1) == 1
    assert igcd(1, 0) == 1
    assert igcd(0, 7) == 7
    assert igcd(7, 0) == 7
    assert igcd(7, 1) == 1
    assert igcd(1, 7) == 1
    assert igcd(-1, 0) == 1
    assert igcd(0, -1) == 1
    assert igcd(-1, -1) == 1
    assert igcd(-1, 7) == 1
    assert igcd(7, -1) == 1
    assert igcd(8, 2) == 2
    assert igcd(4, 8) == 4
    assert igcd(8, 16) == 8
    assert igcd(7, -3) == 1
    assert igcd(-7, 3) == 1
    assert igcd(-7, -3) == 1

def test_ilcm():
    assert ilcm(0, 0) == 0
    assert ilcm(1, 0) == 0
    assert ilcm(0, 1) == 0
    assert ilcm(1, 1) == 1
    assert ilcm(2, 1) == 2
    assert ilcm(8, 2) == 8
    assert ilcm(8, 6) == 24
    assert ilcm(8, 7) == 56

def test_igcdex():
    assert igcdex(2, 3) == (-1, 1, 1)
    assert igcdex(10, 12) == (-1, 1, 2)
    assert igcdex(100, 2004) == (-20, 1, 4)

def _strictly_equal(a, b):
    return (a.p, a.q, type(a.p), type(a.q)) == \
           (b.p, b.q, type(b.p), type(b.q))

def _test_rational_new(cls):
    """
    Tests that are common between Integer and Rational.
    """
    assert cls(0) is S.Zero
    assert cls(1) is S.One
    assert cls(-1) is S.NegativeOne
    # These look odd, but are similar to int():
    assert cls(0.9) is S.Zero
    assert cls('1') is S.One
    assert cls(u'-1') is S.NegativeOne

    i = Integer(10)
    assert _strictly_equal(i, cls('10'))
    assert _strictly_equal(i, cls(u'10'))
    assert _strictly_equal(i, cls(10L))
    assert _strictly_equal(i, cls(10.5))
    assert _strictly_equal(i, cls(i))

    raises(TypeError, "cls(Symbol('x'))")

def test_Integer_new():
    """
    Test for Integer constructor
    """
    _test_rational_new(Integer)

    raises(ValueError, 'Integer("10.5")')

def test_Rational_new():
    """"
    Test for Rational constructor
    """
    _test_rational_new(Rational)

    n1 = Rational(1,2)
    assert n1 == Rational(Integer(1), 2)
    assert n1 == Rational(Integer(1), Integer(2))
    assert n1 == Rational(1, Integer(2))
    assert n1 == Rational(Rational(1,2))
    assert n1 == Rational(1.2,2)

def test_Rational_cmp():
    n1 = Rational(1,4)
    n2 = Rational(1,3)
    n3 = Rational(2,4)
    n4 = Rational(2,-4)
    n5 = Rational(0)
    n6 = Rational(1)
    n7 = Rational(3)
    n8 = Rational(-3)

    assert n8<n5
    assert n5<n6
    assert n6<n7
    assert n8<n7
    assert n7>n8
    assert (n1+1)**n2 < 2
    assert ((n1+n6)/n7) < 1

    assert n4<n3
    assert n2<n3
    assert n1<n2
    assert n3>n1
    assert not n3<n1
    assert not (Rational(-1) > 0)
    assert Rational(-1) < 0


def test_Real():
    def eq(a, b):
        t = Real("1.0E-15")
        return (-t < a-b < t)

    a = Real(2) ** Real(3)
    assert eq(a.evalf(), Real(8))
    assert eq((pi ** -1).evalf(), Real("0.31830988618379067"))
    a = Real(2) ** Real(4)
    assert eq(a.evalf(), Real(16))
    assert (S(.3) == S(.5)) is False
    x_str = Real((0, '13333333333333L', -52, 53))
    x2_str = Real((0, '26666666666666L', -53, 53))
    x_hex = Real((0, 0x13333333333333L, -52, 53))
    x_dec = Real((0, 5404319552844595L, -52, 53))
    x2_hex = Real((0, 0x13333333333333L*2, -53, 53))
    assert x_str == x_hex == x_dec == x2_hex == Real(1.2)
    # x2_str and 1.2 are superficially the same
    assert str(x2_str) == str(Real(1.2))
    # but are different at the mpf level
    assert Real(1.2)._mpf_ == (0, 5404319552844595L, -52, 53)
    assert x2_str._mpf_ == (0, 10808639105689190L, -53, 53)

def test_Real_eval():
    a = Real(3.2)
    assert (a**2).is_Real

def test_Infinity():
    assert oo == oo
    assert oo != 1
    assert 1*oo == oo
    assert 1 != oo
    assert oo != -oo
    assert oo != Symbol("x")**3
    assert oo + 1 == oo + 1
    assert oo + 1 == oo
    assert 2 + oo == oo
    assert 3*oo + 2 == oo
    assert -oo*3 == -oo
    assert oo + oo == oo
    assert -oo + oo*(-5) == -oo
    assert 1/oo  == 0
    assert 1/(-oo)  == 0
    assert 8/oo  == 0
    assert oo % 2 == nan
    assert 2 % oo == nan

def test_Infinity_2():
    x = Symbol('x')
    assert oo*x != oo
    assert oo*(pi-1) == oo
    assert oo*(1-pi) == -oo

    assert (-oo)*x != -oo
    assert (-oo)*(pi-1) == -oo
    assert (-oo)*(1-pi) == oo

def test_NaN():
    assert nan == nan
    assert nan != 1
    assert 1*nan == nan
    assert 1 != nan
    assert nan == -nan
    assert oo != Symbol("x")**3
    assert nan + 1 == nan
    assert 2 + nan == nan
    assert 3*nan + 2 == nan
    assert -nan*3 == nan
    assert nan + nan == nan
    assert -nan + nan*(-5) == nan
    assert 1/nan  == nan
    assert 1/(-nan)  == nan
    assert 8/nan  == nan

def test_powers():
    assert integer_nthroot(1, 2) == (1, True)
    assert integer_nthroot(1, 5) == (1, True)
    assert integer_nthroot(2, 1) == (2, True)
    assert integer_nthroot(2, 2) == (1, False)
    assert integer_nthroot(2, 5) == (1, False)
    assert integer_nthroot(4, 2) == (2, True)
    assert integer_nthroot(123**25, 25) == (123, True)
    assert integer_nthroot(123**25+1, 25) == (123, False)
    assert integer_nthroot(123**25-1, 25) == (122, False)
    assert integer_nthroot(1,1) == (1, True)
    assert integer_nthroot(0,1) == (0, True)
    assert integer_nthroot(0,3) == (0, True)
    assert integer_nthroot(10000, 1) == (10000, True)
    assert integer_nthroot(4,2) == (2, True)
    assert integer_nthroot(16,2) == (4, True)
    assert integer_nthroot(26,2) == (5, False)
    assert integer_nthroot(1234567**7, 7) == (1234567, True)
    assert integer_nthroot(1234567**7+1, 7) == (1234567, False)
    assert integer_nthroot(1234567**7-1, 7) == (1234566, False)
    b = 25**1000
    assert integer_nthroot(b, 1000) == (25, True)
    assert integer_nthroot(b+1, 1000) == (25, False)
    assert integer_nthroot(b-1, 1000) == (24, False)
    c = 10**400
    c2 = c**2
    assert integer_nthroot(c2, 2) == (c, True)
    assert integer_nthroot(c2+1, 2) == (c, False)
    assert integer_nthroot(c2-1, 2) == (c-1, False)
    assert integer_nthroot(2,10**10) == (1, False)

    p, r = integer_nthroot(int(factorial(10000)), 100)
    assert p % (10**10) == 5322420655
    assert not r

    # Test that this is fast
    assert integer_nthroot(2,10**10) == (1, False)

def test_powers_Integer():
    """Test Integer._eval_power"""
    # check infinity
    assert S(1) ** S.Infinity == 1
    assert S(-1)** S.Infinity == S.NaN
    assert S(2) ** S.Infinity == S.Infinity
    assert S(-2)** S.Infinity == S.Infinity + S.Infinity * S.ImaginaryUnit
    assert S(0) ** S.Infinity == 0

    # check Nan
    assert S(1)  ** S.NaN == S.NaN
    assert S(-1) ** S.NaN == S.NaN

    # check for exact roots
    assert S(-1)  ** Rational(6, 5) == - (-1)**(S(1)/5)
    assert S(4)   ** Rational(1, 2) == 2
    assert S(-4)  ** Rational(1, 2) == I * 2
    assert S(16)  ** Rational(1, 4) == 2
    assert S(-16) ** Rational(1, 4) == 2 * (-1)**Rational(1,4)
    assert S(9)   ** Rational(3, 2) == 27
    assert S(-9)  ** Rational(3, 2) == -27*I
    assert S(27)  ** Rational(2, 3) == 9
    assert S(-27) ** Rational(2, 3) == 9 * (S(-1) ** Rational(2, 3))
    assert (-2) ** Rational(-2, 1) == Rational(1, 4)

    # not exact roots
    assert (-3) ** (S(1)/2)  == sqrt(-3)
    assert (3)  ** (S(3)/2)  == 3 * sqrt(3)
    assert (-3) ** (S(3)/2)  == - 3 * sqrt(-3)
    assert (-3) ** (S(5)/2)  ==  9 * I * sqrt(3)
    assert (-3) ** (S(7)/2)  == - I * 27 * sqrt(3)
    assert (2)  ** (S(3)/2)  == 2 * sqrt(2)
    assert (2)  ** (S(-3)/2) == sqrt(2) / 4
    assert (81) ** (S(2)/3)  == 9 * (S(3) ** (S(2)/3))
    assert (-81) ** (S(2)/3)  == 9 * (S(-3) ** (S(2)/3))
    assert (-3) ** Rational(-7, 3) == -(-3) ** Rational(2, 3) / 27
    assert (-3) ** Rational(-2, 3) == -(-3) ** Rational(1, 3) / 3

    # join roots
    assert sqrt(6) + sqrt(24) == 3*sqrt(6)
    assert sqrt(2) * sqrt(3)  == sqrt(6)

    # separate sybols & constansts
    x = Symbol("x")
    assert sqrt(49 * x) == 7 * sqrt(x)
    assert sqrt((3 - sqrt(pi)) ** 2) == 3 - sqrt(pi)

    # check that it is fast for big numbers
    assert (2**64+1) ** Rational(4, 3)
    assert (2**64+1) ** Rational(17,25)

    # negative rational power and negative base
    assert (-3) ** Rational(-7, 3) == -(-3) ** Rational(2, 3) / 27
    assert (-3) ** Rational(-2, 3) == -(-3) ** (S(1) / 3) / 3

def test_powers_Rational():
    """Test Rational._eval_power"""
    # check infinity
    assert Rational(1,2) ** S.Infinity == 0
    assert Rational(3,2) ** S.Infinity == S.Infinity
    assert Rational(-1,2) ** S.Infinity == 0
    assert Rational(-3,2)** S.Infinity == S.Infinity + S.Infinity * S.ImaginaryUnit

    # check Nan
    assert Rational(3,4)  ** S.NaN == S.NaN
    assert Rational(-2,3) ** S.NaN == S.NaN

    # exact roots on numerator
    assert Rational(4,3)  ** Rational(1,2) == 2 * sqrt(3) / 3
    assert Rational(4,3)  ** Rational(3,2) == 8 * sqrt(3) / 9
    assert Rational(-4,3) ** Rational(1,2) == I * 2 * sqrt(3) / 3
    assert Rational(-4,3) ** Rational(3,2) == - I * 8 * sqrt(3) / 9
    assert Rational(27,2) ** Rational(1,3) == 3 * (2 ** Rational(2,3)) / 2
    assert Rational(5**3, 8**3) ** Rational(4,3) == Rational(5**4, 8**4)

    # exact root on denominator
    assert Rational(1,4)  ** Rational(1,2) == Rational(1,2)
    assert Rational(1,-4) ** Rational(1,2) == I * Rational(1,2)
    assert Rational(3,4)  ** Rational(1,2) == sqrt(3) / 2
    assert Rational(3,-4) ** Rational(1,2) == I * sqrt(3) / 2
    assert Rational(5,27) ** Rational(1,3) == (5 ** Rational(1,3)) / 3

    # not exact roots
    assert Rational(1,2)  ** Rational(1,2) == sqrt(2) / 2
    assert Rational(-4,7) ** Rational(1,2) == I * Rational(4,7) ** Rational(1,2)
    assert Rational(-3, 2)**Rational(-7, 3) == \
           -4 * (-3) ** Rational(2, 3)*2 ** Rational(1, 3)/27
    assert Rational(-3, 2)**Rational(-2, 3) == \
           -(-3) ** (S(1) / 3) * 2 ** (S(2) / 3) / 3

    # negative integer power and negative rational base
    assert Rational(-2, 3) ** Rational(-2, 1) == Rational(9, 4)

def test_abs1():
    assert Rational(1,6) != Rational(-1,6)
    assert abs(Rational(1,6)) == abs(Rational(-1,6))

def test_accept_int():
    assert Real(4) == 4

def test_dont_accept_str():
    assert      Real("0.2") != "0.2"
    assert not (Real("0.2") == "0.2")

def test_int():
    a = Rational(5)
    assert int(a)==5

def test_real_bug():
    x = Symbol("x")
    assert str(2.0*x*x) in ["(2.0*x)*x","2.0*x**2","2.00000000000000*x**2"]
    assert str(2.1*x*x)!="(2.0*x)*x"

def test_bug_sqrt():
    assert ((sqrt(Rational(2))+1)*(sqrt(Rational(2))-1)).expand() == 1

def test_pi_Pi():
    "Test, that pi (instance) is imported, but Pi (class) is not"
    from sympy import pi
    raises(ImportError, "from sympy import Pi")

def test_no_len():
    # there should be no len for numbers
    raises(TypeError, "len(Rational(2))")
    raises(TypeError, "len(Rational(2,3))")
    raises(TypeError, "len(Integer(2))")

def test_issue222():
    assert sqrt(Rational(1, 5)) == Rational(1, 5)**S.Half
    assert 5 * Rational(1,5)**Rational(1,2) == 5 * sqrt(Rational(1,5))

def test_issue593():
    assert ((-1)**Rational(1,6)).expand(complex=True) == I/2 + sqrt(3)/2
    assert ((-5)**Rational(1,6)).expand(complex=True) == \
            5**Rational(1,6)*I/2 + 5**Rational(1,6)*sqrt(3)/2
    assert ((-64)**Rational(1,6)).expand(complex=True) == I + sqrt(3)

def test_issue324():
    x = Symbol("x")
    assert sqrt(x-1) == (x-1)**Rational(1,2)
    assert sqrt(x-1) != I*(1-x)**Rational(1,2)

def test_issue350():
    x = Symbol("x", real=True)
    assert sqrt(x**2) == abs(x)
    assert sqrt(x-1).subs(x,5) == 2


def test_Integer_factors():
    def F(i):
        return Integer(i).factors()

    assert F(1)   == { 1:1}
    assert F(2)   == { 2:1}
    assert F(3)   == { 3:1}
    assert F(4)   == { 2:2}
    assert F(5)   == { 5:1}
    assert F(6)   == { 2:1,  3:1}
    assert F(7)   == { 7:1}
    assert F(8)   == { 2:3}
    assert F(9)   == { 3:2}
    assert F(10)  == { 2:1,  5:1}
    assert F(11)  == {11:1}
    assert F(12)  == { 2:2,  3:1}
    assert F(13)  == {13:1}
    assert F(14)  == { 2:1,  7:1}
    assert F(15)  == { 3:1,  5:1}
    assert F(16)  == { 2:4}
    assert F(17)  == {17:1}
    assert F(18)  == { 2:1,  3:2}
    assert F(19)  == {19:1}
    assert F(20)  == { 2:2,  5:1}
    assert F(21)  == { 3:1,  7:1}
    assert F(22)  == { 2:1, 11:1}
    assert F(23)  == {23:1}
    assert F(24)  == { 2:3,  3:1}
    assert F(25)  == { 5:2}
    assert F(26)  == { 2:1, 13:1}
    assert F(27)  == { 3:3}
    assert F(28)  == { 2:2,  7:1}
    assert F(29)  == {29:1}
    assert F(30)  == { 2:1,  3:1,  5:1}
    assert F(31)  == {31:1}
    assert F(32)  == { 2:5}
    assert F(33)  == { 3:1, 11:1}
    assert F(34)  == { 2:1, 17:1}
    assert F(35)  == { 5:1,  7:1}
    assert F(36)  == { 2:2,  3:2}
    assert F(37)  == {37:1}
    assert F(38)  == { 2:1, 19:1}
    assert F(39)  == { 3:1, 13:1}
    assert F(40)  == { 2:3,  5:1}
    assert F(41)  == {41:1}
    assert F(42)  == { 2:1,  3:1,  7:1}
    assert F(43)  == {43:1}
    assert F(44)  == { 2:2, 11:1}
    assert F(45)  == { 3:2,  5:1}
    assert F(46)  == { 2:1, 23:1}
    assert F(47)  == {47:1}
    assert F(48)  == { 2:4,  3:1}
    assert F(49)  == { 7:2}
    assert F(50)  == { 2:1,  5:2}
    assert F(51)  == { 3:1, 17:1}


def test_Rational_factors():
    def F(p,q):
        return Rational(p,q).factors()

    assert F(2,3)   == { 2:1, 3:-1}
    assert F(2,9)   == { 2:1, 3:-2}
    assert F(2,15)  == { 2:1, 3:-1, 5:-1}
    assert F(6,10)  == { 3:1, 5:-1}

    # TODO write more Rational.factor() tests


def test_issue1008():
    assert  pi*(E + 10) + pi*(-E - 10)          == 0
    assert  pi*(E + 10**10) + pi*(-E - 10**10)  == 0
    assert  pi*(E + 10**20) + pi*(-E - 10**20)  == 0
    assert  pi*(E + 10**80) + pi*(-E - 10**80)  == 0


def test_IntegerInteger():
    a = Integer(4)
    b = Integer(a)

    assert a == b

def test_Integer_methods():
    assert Integer(100).half_gcdex(2004) == \
        (Integer(-20), Integer(4))
    assert Integer(100).half_gcdex(Integer(2004)) == \
        (Integer(-20), Integer(4))
    assert Integer(100).gcdex(2004) == \
        (Integer(-20), Integer(1), Integer(4))
    assert Integer(100).gcdex(Integer(2004)) == \
        (Integer(-20), Integer(1), Integer(4))

    raises(ValueError, "Integer(3).half_gcdex(Rational(1,2))")
    raises(ValueError, "Integer(3).gcdex(Rational(1,2))")

    assert Integer(3).invert(7) == Integer(5)
    assert Integer(3).invert(Integer(7)) == Integer(5)

    raises(ValueError, "Integer(3).invert(Rational(1,2))")

    assert Integer(4).cofactors(2) == \
        (Integer(2), Integer(2), Integer(1))
    assert Integer(4).cofactors(Integer(2)) == \
        (Integer(2), Integer(2), Integer(1))

    raises(ValueError, "Integer(4).cofactors(Rational(1,2))")

    assert Integer(4).gcd(2) == Integer(2)
    assert Integer(4).lcm(2) == Integer(4)
    assert Integer(4).gcd(Integer(2)) == Integer(2)
    assert Integer(4).lcm(Integer(2)) == Integer(4)

    raises(ValueError, "Integer(4).gcd(Rational(1,2))")
    raises(ValueError, "Integer(4).lcm(Rational(1,2))")

def test_issue1512():
    assert abs(pi._evalf(50) - 3.14159265358979) < 1e-10
    assert abs(E._evalf(50) - 2.71828182845905) < 1e-10
    assert abs(Catalan._evalf(50) - 0.915965594177219) < 1e-10
    assert abs(EulerGamma._evalf(50) - 0.577215664901533) < 1e-10
    assert abs(GoldenRatio._evalf(50) - 1.61803398874989) < 1e-10
    x = Symbol("x")
    assert (pi+x).evalf() == pi.evalf()+x
    assert (E+x).evalf() == E.evalf()+x
    assert (Catalan+x).evalf() == Catalan.evalf()+x
    assert (EulerGamma+x).evalf() == EulerGamma.evalf()+x
    assert (GoldenRatio+x).evalf() == GoldenRatio.evalf()+x

def test_conversion_to_mpmath():
    assert mpmath.mpmathify(Integer(1)) == mpmath.mpf(1)
    assert mpmath.mpmathify(Rational(1, 2)) == mpmath.mpf(0.5)
    assert mpmath.mpmathify(Real('1.23')) == mpmath.mpf('1.23')

def test_relational():
    # real
    x = S(.1)
    assert (x != cos) is True
    assert (x == cos) is False

    # rational
    x = Rational(1,3)
    assert (x != cos) is True
    assert (x == cos) is False

    # integer defers to rational so these tests are omitted

    # number symbol
    x = pi
    assert (x != cos) is True
    assert (x == cos) is False

