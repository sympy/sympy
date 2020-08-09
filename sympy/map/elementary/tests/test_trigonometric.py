from sympy import S, symbols, pi, Interval, Rational
from sympy.map import Sin, Cos, Tan
from sympy.testing.pytest import raises

sin = Sin(S.Complexes)
cos = Cos(S.Complexes)
tan = Tan(S.Complexes)

def test_Sin():
    x, y = symbols('x y')

    assert sin.nargs == 1

def test_trig_domain():
    assert Sin(S.Complexes).domain == S.Complexes
    assert Sin(S.Complexes).codomain == S.ComplexesField
    assert Sin(S.Complexes).range == S.Complexes
    assert Sin(S.Reals).domain == S.Reals
    assert Sin(S.Reals).codomain == S.RealsField
    assert Sin(S.Reals).range == Interval(-1, 1)

    assert Cos(S.Complexes).domain == S.Complexes
    assert Cos(S.Complexes).codomain == S.ComplexesField
    assert Cos(S.Complexes).range == S.Complexes
    assert Cos(S.Reals).domain == S.Reals
    assert Cos(S.Reals).codomain == S.RealsField
    assert Cos(S.Reals).range == Interval(-1, 1)

    assert Tan(S.Complexes).domain == S.Complexes
    assert Tan(S.Complexes).codomain == S.ComplexesField
    assert Tan(S.Complexes).range == S.Complexes
    assert Tan(S.Reals).domain == S.Reals
    assert Tan(S.Reals).codomain == S.RealsField
    assert Tan(S.Reals).range == S.Reals

def test_trig_period():
    x, y = symbols('x y')

    assert sin(x).period() == 2*pi
    assert sin(2*x).period() == pi
    assert sin(3*x*y + 2*pi).period(y) == 2*pi/abs(3*x)
    assert cos(x).period() == 2*pi
    assert cos((-3)*x).period() == pi*Rational(2, 3)
    assert tan(x).period() == pi
    assert tan(3*x).period(y) is S.Zero
    raises(NotImplementedError, lambda: sin(x**2).period(x))
