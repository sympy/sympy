from sympy.rubi.utility_function import *
from sympy.core.symbol import symbols, S
from sympy.functions.elementary.trigonometric import atan, acsc, asin, acot, acos, asec
from sympy.functions.elementary.hyperbolic import acosh, asinh, atanh, acsch, cosh, sinh, tanh, coth, sech, csch
from sympy.functions import (log, sin, cos, tan, cot, sec, csc, sqrt)
from sympy import I, E, pi

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
    assert FracPart(S(10)+0.5) == 10.5

def test_IntPart():
    assert IntPart(S(10)) == 10
    #assert IntPart(3*pi)


def test_NegQ():
    assert NegQ(-S(3))
    assert not NegQ(S(0))
    assert not NegQ(S(0))

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

def test_Re():
    assert Re(1 + I) == 1

def test_Im():
    assert Im(1 + 2*I) == 2

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

def test_NonsumQ():
    assert NonsumQ(a*b) == True
    assert NonsumQ(a + b) == False

def test_SqrtNumberQ():
    assert SqrtNumberQ(sqrt(2)) == True

def test_IntLinearcQ():
    assert IntLinearcQ(1, 2, 3, 4, 5, 6, x) == True
    assert IntLinearcQ(S(1)/100, S(2)/100, S(3)/100, S(4)/100, S(5)/100, S(6)/100, x) == False

def test_IndependentQ():
    assert IndependentQ(a + b*x, x) == False
    assert IndependentQ(a + b, x) == True

def test_PowerQ():
    assert PowerQ(a**b) == True
    assert PowerQ(a + b) == False

def test_IntegerPowerQ():
    assert IntegerPowerQ(a**2) == True
    assert IntegerPowerQ(a**0.5) == False

def test_PositiveIntegerPowerQ():
    assert PositiveIntegerPowerQ(a**3) == True
    assert PositiveIntegerPowerQ(a**(-2)) == False

def test_FractionalPowerQ():
    assert FractionalPowerQ(a**2) == True
    assert FractionalPowerQ(a**sqrt(2)) == False

def test_AtomQ():
    assert AtomQ(x)
    assert not AtomQ(x+1)

def test_ExpQ():
    assert ExpQ(E**2)
    assert not ExpQ(2**E)

def test_LogQ():
    assert LogQ(log(x))
    assert not LogQ(sin(x) + log(x))

def test_Head():
    assert Head(sin(x)) == sin
    assert Head(log(x**3 + 3)) == log

def test_MemberQ():
    assert MemberQ([a, b, c], b)
    assert MemberQ([sin, cos, log, tan], Head(sin(x)))

def test_TrigQ():
    assert TrigQ(sin(x))
    assert TrigQ(tan(x**2 + 2))
    assert not TrigQ(sin(x) + tan(x))

def test_SinQ():
    assert SinQ(sin(x))
    assert not SinQ(tan(x))

def test_CosQ():
    assert CosQ(cos(x))
    assert not CosQ(csc(x))

def test_TanQ():
    assert TanQ(tan(x))
    assert not TanQ(cot(x))

def test_CotQ():
    assert not CotQ(tan(x))
    assert CotQ(cot(x))

def test_SecQ():
    assert SecQ(sec(x))
    assert not SecQ(csc(x))

def test_CscQ():
    assert not CscQ(sec(x))
    assert CscQ(csc(x))

def test_HyperbolicQ():
    assert HyperbolicQ(sinh(x))
    assert HyperbolicQ(cosh(x))
    assert HyperbolicQ(tanh(x))
    assert not HyperbolicQ(sinh(x) + cosh(x) + tanh(x))

def test_SinhQ():
    assert SinhQ(sinh(x))
    assert not SinhQ(cosh(x))

def test_CoshQ():
    assert not CoshQ(sinh(x))
    assert CoshQ(cosh(x))

def test_TanhQ():
    assert TanhQ(tanh(x))
    assert not TanhQ(coth(x))

def test_CothQ():
    assert not CothQ(tanh(x))
    assert CothQ(coth(x))

def test_SechQ():
    assert SechQ(sech(x))
    assert not SechQ(csch(x))

def test_CschQ():
    assert not CschQ(sech(x))
    assert CschQ(csch(x))

def test_InverseTrigQ():
    assert InverseTrigQ(acot(x))
    assert InverseTrigQ(asec(x))
    assert not InverseTrigQ(acsc(x) + asec(x))

def test_SinCosQ():
    assert SinCosQ(sin(x))
    assert SinCosQ(cos(x))
    assert SinCosQ(sec(x))
    assert not SinCosQ(acsc(x))

def test_SinhCoshQ():
    assert not SinhCoshQ(sin(x))
    assert SinhCoshQ(cosh(x))
    assert SinhCoshQ(sech(x))
    assert SinhCoshQ(csch(x))

def test_Rt():
    assert Rt(8, 3) == 2
    assert Rt(16807, 5) == 7

def test_LeafCount():
    assert LeafCount(1 + a + x**2) == 6

def test_Numerator():
    assert Numerator(S(3)/2) == 3
    assert Numerator(x/y) == x

def test_Length():
    assert Length(a + b) == 2
    assert Length(sin(a)*cos(a)) == 2

def test_AtomQ():
    assert AtomQ(a)
    assert not AtomQ(a + b)

def test_ListQ():
    assert ListQ([1, 2])
    assert not ListQ(a)

def test_InverseHyperbolicQ():
    assert InverseHyperbolicQ(acosh(a))

def test_InverseFunctionQ():
    assert InverseFunctionQ(log(a))
    assert InverseFunctionQ(acos(a))
    assert not InverseFunctionQ(a)
    assert InverseFunctionQ(acosh(a))
    assert InverseFunctionQ(polylog(a, b))

def test_EqQ():
    assert EqQ(a, a)
    assert not EqQ(a, b)


def test_Rest():
    assert Rest([2, 3, 5, 7]) == [3, 5, 7]
    assert Rest(1/b) == -1
    assert Rest(a + b + c) == b + c
    assert Rest(a*b*c) == b*c

def test_First():
    assert First([2, 3, 5, 7]) == 2
    assert First(y**2) == y
    assert First((1/b)) == b
    assert First(a + b + c) == a
    assert First(a*b*c) == a

def test_ComplexFreeQ():
    assert ComplexFreeQ(x)
    assert not ComplexFreeQ(x+2*I)

def test_FractionalPowerFreeQ():
    assert not FractionalPowerFreeQ(x**(S(2)/3))
    assert FractionalPowerFreeQ(x)

def test_FactorSquareFree():
    assert FactorSquareFree((x**5 - x**3 - x**2 + 1)*(x+1)**2) == (x - 1)**2*(x + 1)**3*(x**2 + x + 1)
    assert FactorSquareFree((x**5 - x**3 - x**2 + 1)*(x-1)) == (x - 1)**3*(x**3 + 2*x**2 + 2*x + 1)
    assert FactorSquareFree(x**5 - x**3 - x**2 + 1) == (x - 1)**2*(x**3 + 2*x**2 + 2*x + 1)

def test_Exponent():
    assert Exponent(x**2+x+1+5, x, List) == [0, 1, 2]
    assert Exponent(x**2+x+1, x, List) == [0, 1, 2]
    assert Exponent(x**2+2*x+1, x, List) == [0, 2, 1]
    assert Exponent(x**3+x+1, x) == 3
    assert Exponent(x**2+2*x+1, x) == 2
    assert Exponent(x**3, x, List) == [3]

def test_QuadraticQ():
    assert not QuadraticQ([x**2+x+1, 5*x**2], x)
    assert QuadraticQ([x**2+x+1, 5*x**2+3*x+6], x)
    assert not QuadraticQ(x**2+1+x**3, x)
    assert QuadraticQ(x**2+1+x, x)
    assert not QuadraticQ(x**2, x)

def test_BinomialParts():
    assert BinomialParts(2 + x*(9*x), x) == [2, 9, 2]
    assert BinomialParts(x**9, x) == [0, 1, 9]
    assert BinomialParts(2*x**3, x) == [0, 2, 3]
    assert BinomialParts(2 + x, x) == [2, 1, 1]

def test_PolynomialQ():
    assert PolynomialQ(x**3, x)
    assert not PolynomialQ(sqrt(x), x)

def test_PolyQ():
    assert PolyQ(x, x, 1)
    assert PolyQ(x**2, x, 2)
    assert not PolyQ(x**3, x, 2)

def test_EvenQ():
    assert EvenQ(S(2))
    assert not EvenQ(S(1))

def test_OddQ():
    assert OddQ(S(1))
    assert not OddQ(S(2))

def test_PerfectSquareQ():
    assert PerfectSquareQ(S(4))
    assert PerfectSquareQ(a**S(2)*b**S(4))
    assert not PerfectSquareQ(S(1)/3)

def test_NiceSqrtQ():
    assert NiceSqrtQ(S(1)/3)
    assert not NiceSqrtQ(-S(1))
    assert NiceSqrtQ(pi**2)
    assert NiceSqrtQ(pi**2*sin(4)**4)
    assert not NiceSqrtQ(pi**2*sin(4)**3)

def test_Together():
    assert Together(1/a + b/2) == (a*b + 2)/(2*a)

def test_PosQ():
    assert not PosQ(S(0))
    assert PosQ(S(1))
    assert PosQ(pi)
    assert PosQ(pi**3)
    assert PosQ((-pi)**4)
    assert PosQ(sin(1)**2*pi**4)

def test_NumericQ():
    assert NumericQ(sin(cos(2)))

def test_NumberQ():
    assert NumberQ(pi)

def test_CoefficientList():
    assert CoefficientList(1 + a*x, x) == [1, a]
    assert CoefficientList(1 + a*x**3, x) == [1, 0, 0, a]
    assert CoefficientList(sqrt(x), x) == []

def test_ReplaceAll():
    assert ReplaceAll(x, x, a) == a
    assert ReplaceAll(a*x, x, a + b) == a*(a + b)

def test_SimplifyTerm():
    assert SimplifyTerm(a/100 + 100/b*x, x) == a/100 + 100/b*x

def test_ExpandLinearProduct():
    assert ExpandLinearProduct(log(x), x**2, a, b, x) == a**2*log(x)/b**2 - 2*a*(a + b*x)*log(x)/b**2 + (a + b*x)**2*log(x)/b**2

def test_ExpandIntegrand():

    assert True

def test_MatchQ():
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x])
    c_ = Wild('c', exclude=[x])
    assert MatchQ(a*b + c, a_*b_ + c_, a_, b_, c_) == (a, b, c)

def test_PolynomialQuotientRemainder():
    assert PolynomialQuotientRemainder(x**2, x+a, x) == [-a + x, a**2]

