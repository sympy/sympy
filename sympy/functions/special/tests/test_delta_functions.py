from sympy import (
    adjoint, conjugate, DiracDelta, Heaviside, nan, pi, sign, sqrt,
    symbols, transpose, Symbol, Piecewise, I, S, Eq
)

from sympy.utilities.pytest import raises

from sympy.core.function import ArgumentIndexError

x, y = symbols('x y')


def test_DiracDelta():
    assert DiracDelta(1) == 0
    assert DiracDelta(5.1) == 0
    assert DiracDelta(-pi) == 0
    assert DiracDelta(5, 7) == 0
    assert DiracDelta(nan) == nan
    assert DiracDelta(0).func is DiracDelta
    assert DiracDelta(x).func is DiracDelta
    # FIXME: this is generally undefined @ x=0
    #         But then limit(Delta(c)*Heaviside(x),x,-oo)
    #         need's to be implemented.
    #assert 0*DiracDelta(x) == 0

    assert adjoint(DiracDelta(x)) == DiracDelta(x)
    assert adjoint(DiracDelta(x - y)) == DiracDelta(x - y)
    assert conjugate(DiracDelta(x)) == DiracDelta(x)
    assert conjugate(DiracDelta(x - y)) == DiracDelta(x - y)
    assert transpose(DiracDelta(x)) == DiracDelta(x)
    assert transpose(DiracDelta(x - y)) == DiracDelta(x - y)

    assert DiracDelta(x).diff(x) == DiracDelta(x, 1)
    assert DiracDelta(x, 1).diff(x) == DiracDelta(x, 2)

    assert DiracDelta(x).is_simple(x) is True
    assert DiracDelta(3*x).is_simple(x) is True
    assert DiracDelta(x**2).is_simple(x) is False
    assert DiracDelta(sqrt(x)).is_simple(x) is False
    assert DiracDelta(x).is_simple(y) is False

    assert DiracDelta(x*y).simplify(x) == DiracDelta(x)/abs(y)
    assert DiracDelta(x*y).simplify(y) == DiracDelta(y)/abs(x)
    assert DiracDelta(x**2*y).simplify(x) == DiracDelta(x**2*y)
    assert DiracDelta(y).simplify(x) == DiracDelta(y)
    assert DiracDelta((x - 1)*(x - 2)*(x - 3)).simplify(x) == \
        DiracDelta(x - 3)/2 + DiracDelta(x - 2) + DiracDelta(x - 1)/2

    raises(ArgumentIndexError, lambda: DiracDelta(x).fdiff(2))
    raises(ValueError, lambda: DiracDelta(x, -1))


def test_heaviside():
    assert Heaviside(0).func == Heaviside
    assert Heaviside(-5) == 0
    assert Heaviside(1) == 1
    assert Heaviside(nan) == nan

    assert adjoint(Heaviside(x)) == Heaviside(x)
    assert adjoint(Heaviside(x - y)) == Heaviside(x - y)
    assert conjugate(Heaviside(x)) == Heaviside(x)
    assert conjugate(Heaviside(x - y)) == Heaviside(x - y)
    assert transpose(Heaviside(x)) == Heaviside(x)
    assert transpose(Heaviside(x - y)) == Heaviside(x - y)

    assert Heaviside(x).diff(x) == DiracDelta(x)
    assert Heaviside(x + I).is_Function is True
    assert Heaviside(I*x).is_Function is True

    raises(ArgumentIndexError, lambda: Heaviside(x).fdiff(2))
    raises(ValueError, lambda: Heaviside(I))
    raises(ValueError, lambda: Heaviside(2 + 3*I))


def test_rewrite():
    x, y = Symbol('x', real=True), Symbol('y')
    assert Heaviside(x).rewrite(Piecewise) == \
        Piecewise((1, x > 0), (S(1)/2, Eq(x, 0)), (0, True))
    assert Heaviside(y).rewrite(Piecewise) == Heaviside(y)

    assert Heaviside(x).rewrite(sign) == (sign(x)+1)/2
    assert Heaviside(y).rewrite(sign) == Heaviside(y)

    assert DiracDelta(x-5).rewrite(Piecewise) == \
        Piecewise(('undefined', Eq(x-5, 0)), (0, True))
