from sympy.rubi.utility_function import *
from sympy.core.symbol import symbols, S
from sympy.functions import log, sin, cos
from sympy.functions.elementary.hyperbolic import acosh

a, b, c, d, x, y = symbols('a b c d x y')

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
