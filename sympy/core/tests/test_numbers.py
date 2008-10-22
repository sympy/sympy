from sympy import Rational, Symbol, Real, I, sqrt, oo, nan, pi, E, Integer, \
        Basic, S
from sympy.core.power import integer_nthroot

from sympy.core.numbers import igcd, ilcm, igcdex
from sympy.utilities.pytest import XFAIL, raises


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

def test_Rational_new():
    """"
    Test for Rational constructor
    """
    n1 = Rational(1,2)
    assert n1 == Rational(Integer(1), 2)
    assert n1 == Rational(Integer(1), Integer(2))
    assert n1 == Rational(1, Integer(2))

    raises(ValueError, 'Rational(1.2)')
    raises(ValueError, "Rational(Symbol('x'))")
    raises(ValueError, 'Rational(Rational(1,2))')

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

    assert 64**(Rational(1)/3)==4
    assert 64**(Rational(2)/3)==16
    assert 24*64**(-Rational(1)/2)==3

    assert Rational(5**3, 8**3)**Rational(4,3) == Rational(5**4, 8**4)
    assert Rational(-4,7)**Rational(1,2) == I*Rational(4,7)**Rational(1,2)

    assert sqrt(6) + sqrt(24) == 3*sqrt(6)
    assert sqrt(2) * sqrt(3) == sqrt(6)
    x = Symbol("x")
    assert sqrt(49*x) == 7*sqrt(x)
    assert sqrt((3-sqrt(pi))**2) == 3 - sqrt(pi)
    assert sqrt(Rational(1,2)) == Rational(1,2) * sqrt(2)

    # Test that this is fast
    assert integer_nthroot(2,10**10) == (1, False)

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
