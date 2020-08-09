from sympy import S, symbols, pi, Interval, Rational
from sympy.map import (
    Sine, Cosine, Tangent, Cotangent,
    Secant, Cosecant,
)
from sympy.testing.pytest import raises

sin = Sine(S.Complexes)
cos = Cosine(S.Complexes)
tan = Tangent(S.Complexes)
cot = Cotangent(S.Complexes)
sec = Secant(S.Complexes)
csc = Cosecant(S.Complexes)

def test_Sin():
    x, y = symbols('x y')

    assert sin.nargs == 1

def test_trig_domain():
    assert Sine(S.Complexes).domain == S.Complexes
    assert Sine(S.Complexes).codomain == S.ComplexesField
    assert Sine(S.Complexes).range == S.Complexes
    assert Sine(S.Reals).domain == S.Reals
    assert Sine(S.Reals).codomain == S.RealsField
    assert Sine(S.Reals).range == Interval(-1, 1)

    assert Cosine(S.Complexes).domain == S.Complexes
    assert Cosine(S.Complexes).codomain == S.ComplexesField
    assert Cosine(S.Complexes).range == S.Complexes
    assert Cosine(S.Reals).domain == S.Reals
    assert Cosine(S.Reals).codomain == S.RealsField
    assert Cosine(S.Reals).range == Interval(-1, 1)

    assert Tangent(S.Complexes).domain == S.Complexes
    assert Tangent(S.Complexes).codomain == S.ComplexesField
    assert Tangent(S.Complexes).range == S.Complexes
    assert Tangent(S.Reals).domain == S.Reals
    assert Tangent(S.Reals).codomain == S.RealsField
    assert Tangent(S.Reals).range == S.Reals

    assert Cotangent(S.Complexes).domain == S.Complexes
    assert Cotangent(S.Complexes).codomain == S.ComplexesField
    assert Cotangent(S.Complexes).range == S.Complexes
    assert Cotangent(S.Reals).domain == S.Reals
    assert Cotangent(S.Reals).codomain == S.RealsField
    assert Cotangent(S.Reals).range == S.Reals

    assert Secant(S.Complexes).domain == S.Complexes
    assert Secant(S.Complexes).codomain == S.ComplexesField
    assert Secant(S.Complexes).range == S.Complexes
    assert Secant(S.Reals).domain == S.Reals
    assert Secant(S.Reals).codomain == S.RealsField
    assert Secant(S.Reals).range == S.Reals - Interval(-1, 1)

    assert Cosecant(S.Complexes).domain == S.Complexes
    assert Cosecant(S.Complexes).codomain == S.ComplexesField
    assert Cosecant(S.Complexes).range == S.Complexes
    assert Cosecant(S.Reals).domain == S.Reals
    assert Cosecant(S.Reals).codomain == S.RealsField
    assert Cosecant(S.Reals).range == S.Reals - Interval(-1, 1)

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
