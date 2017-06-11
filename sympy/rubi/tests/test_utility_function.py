from sympy.rubi.utility_function import *
from sympy.core.symbol import symbols, S
from sympy.functions import log, sin, cos
from sympy.functions.elementary.trigonometric import atan, acsc, asin, asin, acos
from sympy.functions.elementary.hyperbolic import acosh, atanh, asinh

a, b, c, d, x, y, z = symbols('a b c d x y z')

def test_ZeroQ():
    assert ZeroQ(S(0))
    assert not ZeroQ(S(10))
    assert not ZeroQ(S(-2))

def test_NonzeroQ():
    assert NonzeroQ(S(1)) == True

def test_FreeQ():
    l = [a*b, x, a + b]
    assert FreeQ(l, x) == False

    l = [a*b, a + b]
    assert FreeQ(l, x) == True

def test_List():
    assert List(a, b, c) == [a, b, c]

def test_Log():
    assert Log(a) == log(a)

def test_PositiveIntegerQ():
    assert PositiveIntegerQ(S(1))
    assert not PositiveIntegerQ(S(-3))
    assert not PositiveIntegerQ(S(0))

def test_NegativeIntegerQ():
    assert not NegativeIntegerQ(S(1))
    assert NegativeIntegerQ(S(-3))
    assert not NegativeIntegerQ(S(0))

def test_PositiveQ():
    assert PositiveQ(S(1))
    assert not PositiveQ(S(-3))
    assert not PositiveQ(S(0))

def test_IntegerQ():
    assert IntegerQ(S(1))
    assert not IntegerQ(S(-1.9))
    assert not IntegerQ(S(0.0))
    assert IntegerQ(S(-1))

def test_RationalQ():
    assert RationalQ(S(5)/6)
    assert not RationalQ(Sqrt(1.6))

def test_PosQ():
    assert PosQ(S(10))
    assert not PosQ(S(-10))
    assert not PosQ(S(0))

def test_FracPart():
    assert FracPart(S(10)) == 0

def test_IntPart():
    assert IntPart(10) == 10
    assert IntPart(3.6) == 3
    assert IntPart(-3.6) == -4

def test_NegQ():
    assert NegQ(-3)
    assert not NegQ(0)
    assert not NegQ(4)

def test_RationalQ():
    assert RationalQ(S(5)/6)
    assert not RationalQ(Sqrt(1.6))

def test_Sqrt():
    assert Sqrt(S(16)) == 4

def test_ArcCosh():
    assert ArcCosh(x) == acosh(x)

def test_LinearQ():
    assert LinearQ(3*x + y**2, x)
    assert not LinearQ(3*x + y**2, y)

def test_Sqrt():
    assert Sqrt(x) == sqrt(x)
    assert Sqrt(25) == 5

def test_Coefficient():
    assert Coefficient(7 + 2*x + 4*x**3, x, 1) == 2
    assert Coefficient(a + b*x + c*x**3, x, 0) == a
    assert Coefficient(a + b*x + c*x**3, x, 4) == 0
    assert Coefficient(b*x + c*x**3, x, 3) == c

def test_RemoveContent():
    assert RemoveContent(3+6*x**3+8*x+2*x**2, x) == 6*x**3 + 2*x**2 + 8*x
    assert RemoveContent(3+6*x**3+8*x+2, x) == 6*x**3 + 8*x
    assert RemoveContent(3+b*x**3+a+2, x) == b*x**3


def test_Denominator():
    assert Denominator(3/2) == 2
    assert Denominator(x/y) == y
    assert Denominator(S(4)/5) == 5
    assert Denominator(3/6) == 2

def test_Hypergeometric2F1():
    assert Hypergeometric2F1(2, (1,2), 4, 0.75) == 1.303703703703703703703704
    assert Hypergeometric2F1(2, (1,2), 4, 1) == 1.6
    assert Hypergeometric2F1(2, (1,2), 4, 0.25j) == (0.9931169055799728251931672 + 0.06154836525312066938147793j)

def test_ArcTan():
    assert ArcTan(x) == atan(x)

def test_Not():
    a = 10
    assert Not(a == 2)

def test_FractionalPart():
    assert FractionalPart(S(3.0)) == 0.0

def test_IntegerPart():
    assert IntegerPart(3.6) == 3
    assert IntegerPart(-3.6) == -4

def test_SumSimplerQ():
    assert not SumSimplerQ(x**3, 3 + 4*x**2 + 8*x**3)
    assert SumSimplerQ(x**3, -x**3)
    assert SumSimplerQ(1+x**2, 2-x**2+x)

def test_SimplerQ():
    assert SimplerQ(x**3, 3+4*x**2+8*x**3)
    assert SimplerQ(x**3, 3+x**2+x)
    assert SimplerQ(x**3, 3*x**4+3*x**5)
    assert not SimplerQ(x**3, 3)
    assert not SimplerQ(x**3, 3*x)
    assert SimplerQ(x**3, 3*x**2)
    assert not SimplerQ(x**3, x**2)

# utility functions used in tests

def test_AppellF1():
    assert AppellF1(1,0,0.5,1,0.5,0.25) == 1.154700538379251529018298

def test_Simplify():
    assert Simplify(sin(x)**2 + cos(x)**2) == 1
    assert Simplify((x**3 + x**2 - x - 1)/(x**2 + 2*x + 1)) == x - 1

def test_Integrate():
    assert Integrate(x**2, x) == x**3/3
    assert Integrate(x**3, x) == integrate(x**3, x)

def test_EllipticPi():
    assert EllipticPi(0.25, 0.25) == 1.956616279119236207279727
    assert EllipticPi(3, 0) == (0.0 - 1.11072073453959156175397j)

def test_EllipticE():
    assert EllipticE(0) == 1.570796326794896619231322
    assert EllipticE(2) == (0.5990701173677961037199612 + 0.5990701173677961037199612j)
    assert EllipticE(0.5+0.25j) == (1.360868682163129682716687 - 0.1238733442561786843557315j)

def test_EllipticF():
    assert EllipticF(0,1) == 0.0
    assert EllipticF(2+3j,0) == (2.0 + 3.0j)
    assert EllipticF(1,1) == 1.226191170883517070813061

def test_arctanh():
    assert ArcTanh(a) == atanh(a)

def test_arcsin():
    assert ArcSin(a) == asin(a)

def test_arcsinh():
    assert ArcSinh(a) == asinh(a)

def test_arccos():
    assert ArcCos(a) == acos(a)

def test_arccsc():
    assert ArcCsc(a) == acsc(a)

def test_arccsch():
    assert ArcCsch(a) == acsch(a)
