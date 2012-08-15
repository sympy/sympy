from sympy import symbols, DiracDelta, Heaviside, nan, oo, sqrt, pi, conjugate

from sympy.utilities.pytest import raises

from sympy.core.function import ArgumentIndexError

from sympy.core import I

x,y = symbols('x y')

def test_DiracDelta():
    assert DiracDelta(1) == 0
    assert DiracDelta(5.1) == 0
    assert DiracDelta(-pi) == 0
    assert DiracDelta(5,7) == 0
    assert DiracDelta(0) == oo
    assert DiracDelta(0,5) == oo
    assert DiracDelta(nan) == nan
    assert DiracDelta(x).func == DiracDelta
    assert conjugate(DiracDelta(x)) == DiracDelta(x)
    assert conjugate(DiracDelta(x - y)) == DiracDelta(x - y)

    assert DiracDelta(x).diff(x) == DiracDelta(x,1)
    assert DiracDelta(x,1).diff(x) == DiracDelta(x,2)

    assert DiracDelta(x).is_simple(x) == True
    assert DiracDelta(3*x).is_simple(x) == True
    assert DiracDelta(x**2).is_simple(x) == False
    assert DiracDelta(sqrt(x)).is_simple(x) == False
    assert DiracDelta(x).is_simple(y) == False

    assert DiracDelta(x*y).simplify(x) == DiracDelta(x)/abs(y)
    assert DiracDelta(x*y).simplify(y) == DiracDelta(y)/abs(x)
    assert DiracDelta(x**2*y).simplify(x) == DiracDelta(x**2*y)
    assert DiracDelta(y).simplify(x) == DiracDelta(y)

    raises(ArgumentIndexError, lambda: DiracDelta(x).fdiff(2))
    raises(ValueError, lambda: DiracDelta(x,-1))

def test_heaviside():
    assert Heaviside(0) == 0.5
    assert Heaviside(-5) == 0
    assert Heaviside(1) == 1
    assert Heaviside(nan) == nan

    assert Heaviside(x).diff(x) == DiracDelta(x)
    assert Heaviside(x+I).is_Function
    assert Heaviside(I*x).is_Function

    raises(ArgumentIndexError, lambda: Heaviside(x).fdiff(2))
    raises(ValueError, lambda: Heaviside(I))
    raises(ValueError, lambda: Heaviside(2+3*I))
