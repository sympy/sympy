from sympy import S, symbols, pi, Interval, Rational
from sympy.map import (
    Sin, Cos, Tan, Cot,
    Sec, Csc,
)
from sympy.testing.pytest import raises

sin = Sin(S.Complexes)
cos = Cos(S.Complexes)
tan = Tan(S.Complexes)
cot = Cot(S.Complexes)
sec = Sec(S.Complexes)
csc = Csc(S.Complexes)

def test_Sin():
    x, y = symbols('x y')

    assert sin.nargs == 1

def test_trig_domain():
    assert Sin(S.Complexes).domain == S.Complexes
    assert Sin(S.Complexes).codomain == S.Complexes
    assert Sin(S.Complexes).range == S.Complexes
    assert Sin(S.Reals).domain == S.Reals
    assert Sin(S.Reals).codomain == S.Reals
    assert Sin(S.Reals).range == Interval(-1, 1)

    assert Cos(S.Complexes).domain == S.Complexes
    assert Cos(S.Complexes).codomain == S.Complexes
    assert Cos(S.Complexes).range == S.Complexes
    assert Cos(S.Reals).domain == S.Reals
    assert Cos(S.Reals).codomain == S.Reals
    assert Cos(S.Reals).range == Interval(-1, 1)

    assert Tan(S.Complexes).domain == S.Complexes
    assert Tan(S.Complexes).codomain == S.Complexes
    assert Tan(S.Complexes).range == S.Complexes
    assert Tan(S.Reals).domain == S.Reals
    assert Tan(S.Reals).codomain == S.Reals
    assert Tan(S.Reals).range == S.Reals

    assert Cot(S.Complexes).domain == S.Complexes
    assert Cot(S.Complexes).codomain == S.Complexes
    assert Cot(S.Complexes).range == S.Complexes
    assert Cot(S.Reals).domain == S.Reals
    assert Cot(S.Reals).codomain == S.Reals
    assert Cot(S.Reals).range == S.Reals

    assert Sec(S.Complexes).domain == S.Complexes
    assert Sec(S.Complexes).codomain == S.Complexes
    assert Sec(S.Complexes).range == S.Complexes
    assert Sec(S.Reals).domain == S.Reals
    assert Sec(S.Reals).codomain == S.Reals
    assert Sec(S.Reals).range == S.Reals - Interval(-1, 1)

    assert Csc(S.Complexes).domain == S.Complexes
    assert Csc(S.Complexes).codomain == S.Complexes
    assert Csc(S.Complexes).range == S.Complexes
    assert Csc(S.Reals).domain == S.Reals
    assert Csc(S.Reals).codomain == S.Reals
    assert Csc(S.Reals).range == S.Reals - Interval(-1, 1)

def test_trig_period():
    x, y = symbols('x y')

    assert sin(x).period() == 2*pi
    assert sin(2*x).period() == pi
    assert sin(3*x*y + 2*pi).period(y) == 2*pi/abs(3*x)
    assert cos(x).period() == 2*pi
    assert cos((-3)*x).period() == pi*Rational(2, 3)
    assert cos(x*y).period(x) == 2*pi/abs(y)
    assert tan(x).period() == pi
    assert tan(3*x).period(y) is S.Zero
    assert cot(x).period() == pi
    assert cot(4*x - 6).period() == pi/4
    assert sec(x).period() == 2*pi
    assert csc(x).period() == 2*pi
    raises(NotImplementedError, lambda: sin(x**2).period(x))
