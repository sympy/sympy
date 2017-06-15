from sympy.rubi.utility_function import *
from sympy.core.symbol import symbols, S
from sympy import I

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
    assert RationalQ(S(5)/6, S(4)/5)
    assert not RationalQ(Sqrt(1.6))
    assert not RationalQ(Sqrt(1.6), S(5)/6)

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

def test_Denominator():
    assert Denominator(S(3)/2) == 2
    assert Denominator(x/y) == y
    assert Denominator(S(4)/5) == 5

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

def test_AppellF1():
    assert AppellF1(1,0,0.5,1,0.5,0.25) == 1.154700538379251529018298

def test_Simplify():
    assert Simplify(sin(x)**2 + cos(x)**2) == 1
    assert Simplify((x**3 + x**2 - x - 1)/(x**2 + 2*x + 1)) == x - 1

def test_EllipticPi():
    assert EllipticPi(0.25, 0.25) == 1.956616279119236207279727
    assert EllipticPi(3, 0) == (0.0 - 1.11072073453959156175397j)

def test_EllipticE():
    assert EllipticE(0) == 1.570796326794896619231322
    assert EllipticE(2) == (0.5990701173677961037199612 + 0.5990701173677961037199612j)
    assert EllipticE(0.5 + 0.25j) == (1.360868682163129682716687 - 0.1238733442561786843557315j)

def test_EllipticF():
    assert EllipticF(0,1) == 0.0
    assert EllipticF(2 + 3j,0) == (2.0 + 3.0j)
    assert EllipticF(1,1) == 1.226191170883517070813061

def test_ArcTanh():
    assert ArcTanh(a) == atanh(a)

def test_ArcSin():
    assert ArcSin(a) == asin(a)

def test_ArcSinh():
    assert ArcSinh(a) == asinh(a)

def test_ArcCos():
    assert ArcCos(a) == acos(a)

def test_ArcCsc():
    assert ArcCsc(a) == acsc(a)

def test_ArcCsch():
    assert ArcCsch(a) == acsch(a)

def test_Equal():
    assert Equal(a, a)
    assert not Equal(a, b)

def test_LessEqual():
    assert LessEqual(1, 2, 3)
    assert LessEqual(1, 1)
    assert not LessEqual(3, 2, 1)

def test_With():
    assert With(Set(x, 3), x + y) == 3 + y
    assert With(List(Set(x, 3), Set(y, c)), x + y) == 3 + c

def test_Less():
    assert Less(1, 2, 3)
    assert not Less(1, 1, 3)

def test_Greater():
    assert Greater(3, 2, 1)
    assert not Greater(3, 2, 2)

def test_GreaterEqual():
    assert GreaterEqual(3, 2, 1)
    assert GreaterEqual(3, 2, 2)
    assert not GreaterEqual(2, 3)

def test_Unequal():
    assert Unequal(1, 2)
    assert not Unequal(1, 1)

def test_FractionQ():
    assert FractionQ(S(1), S(2), S(1)/3)
    assert not FractionQ(sqrt(2))

def test_Expand():
    assert Expand((1 + x)**10) == x**10 + 10*x**9 + 45*x**8 + 120*x**7 + 210*x**6 + 252*x**5 + 210*x**4 + 120*x**3 + 45*x**2 + 10*x + 1

def test_Scan():
    assert list(Scan(sin, [a, b])) == [sin(a), sin(b)]

def test_MapAnd():
    assert MapAnd(PositiveQ, [1, 2, 3, 0]) == False
    assert MapAnd(PositiveQ, [1, 2, 3]) == True

def test_FalseQ():
    assert FalseQ(True) == False
    assert FalseQ(False) == True

def test_ComplexNumberQ():
    assert ComplexNumberQ(1 + I*2, I) == True
    assert ComplexNumberQ(a + b, I) == False

def test_RealNumericQ():
    assert RealNumericQ(S(1)) == True

def test_PositiveOrZeroQ():
    assert PositiveOrZeroQ(S(0)) == True
    assert PositiveOrZeroQ(S(1)) == True
    assert PositiveOrZeroQ(-S(1)) == False

def test_RealNumericQ():
    assert RealNumericQ(S(1)) == True
    assert RealNumericQ(-S(1)) == True

def test_NegativeOrZeroQ():
    assert NegativeOrZeroQ(S(0)) == True
    assert NegativeOrZeroQ(-S(1)) == True
    assert NegativeOrZeroQ(S(1)) == False

def test_FractionOrNegativeQ():
    assert FractionOrNegativeQ(S(1)/2) == True
    assert FractionOrNegativeQ(-S(1)) == True

def test_ProductQ():
    assert ProductQ(a*b) == True
    assert ProductQ(a + b) == False

def test_SumQ():
    assert SumQ(a*b) == False
    assert SumQ(a + b) == True

def test_First():
    assert First([1, 2, 3, 4]) == 1

def test_Rest():
    assert Rest([1, 2, 3, 4]) == [2, 3, 4]

def test_SqrtNumberQ():
    assert SqrtNumberQ(sqrt(2)) == True

def test_IntLinearcQ():
    assert IntLinearcQ(1, 2, 3, 4, 5, 6, x) == True
    assert IntLinearcQ(S(1)/100, S(2)/100, S(3)/100, S(4)/100, S(5)/100, S(6)/100, x) == False
