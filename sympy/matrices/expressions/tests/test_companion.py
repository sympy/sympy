from __future__ import annotations
from sympy.core.expr import unchanged
from sympy.core.numbers import Rational
from sympy.core.symbol import Symbol, symbols
from sympy.matrices.immutable import ImmutableDenseMatrix
from sympy.matrices.expressions.companion import CompanionMatrix
from sympy.matrices.expressions.determinant import Determinant, det
from sympy.matrices.expressions.trace import Trace, trace
from sympy.polys.polytools import Poly
from sympy.testing.pytest import raises


def test_creation():
    x = Symbol('x')
    y = Symbol('y')
    raises(ValueError, lambda: CompanionMatrix(1))
    raises(ValueError, lambda: CompanionMatrix(Poly([1], x)))
    raises(ValueError, lambda: CompanionMatrix(Poly([2, 1], x)))
    raises(ValueError, lambda: CompanionMatrix(Poly(x*y, [x, y])))
    assert unchanged(CompanionMatrix, Poly([1, 2, 3], x))


def test_shape():
    c0, c1, c2 = symbols('c0:3')
    x = Symbol('x')
    assert CompanionMatrix(Poly([1, c0], x)).shape == (1, 1)
    assert CompanionMatrix(Poly([1, c1, c0], x)).shape == (2, 2)
    assert CompanionMatrix(Poly([1, c2, c1, c0], x)).shape == (3, 3)


def test_entry():
    c0, c1, c2 = symbols('c0:3')
    x = Symbol('x')
    A = CompanionMatrix(Poly([1, c2, c1, c0], x))
    assert A[0, 0] == 0
    assert A[1, 0] == 1
    assert A[1, 1] == 0
    assert A[2, 1] == 1
    assert A[0, 2] == -c0
    assert A[1, 2] == -c1
    assert A[2, 2] == -c2


def test_as_explicit():
    c0, c1, c2 = symbols('c0:3')
    x = Symbol('x')
    assert CompanionMatrix(Poly([1, c0], x)).as_explicit() == \
        ImmutableDenseMatrix([-c0])
    assert CompanionMatrix(Poly([1, c1, c0], x)).as_explicit() == \
        ImmutableDenseMatrix([[0, -c0], [1, -c1]])
    assert CompanionMatrix(Poly([1, c2, c1, c0], x)).as_explicit() == \
        ImmutableDenseMatrix([[0, 0, -c0], [1, 0, -c1], [0, 1, -c2]])


def test_trace():
    c0, c1, c2 = symbols('c0:3')
    x = Symbol('x')

    assert Trace(CompanionMatrix(Poly(x + 5, x))).doit() == -5
    assert trace(CompanionMatrix(Poly(x**3 + 2*x**2 + 3*x + 4, x))) == -2
    assert trace(CompanionMatrix(Poly(x**4 + 7*x + 1, x))) == 0
    assert trace(CompanionMatrix(Poly(x**2 + Rational(1, 2)*x + Rational(3, 4), x))) == -Rational(1, 2)
    assert trace(CompanionMatrix(Poly(x**50 + 3*x**49 + 7, x))) == -3

    for p, expected in [
        (Poly([1, c0], x), -c0),
        (Poly([1, c1, c0], x), -c1),
        (Poly([1, c2, c1, c0], x), -c2),
    ]:
        C = CompanionMatrix(p)
        assert Trace(C).doit() == expected
        assert trace(C) == expected
        assert Trace(C).doit() == C.as_explicit().trace()


def test_determinant():
    c0, c1, c2 = symbols('c0:3')
    x = Symbol('x')

    assert Determinant(CompanionMatrix(Poly(x + 5, x))).doit() == -5
    assert det(CompanionMatrix(Poly(x**2 + 2*x + 3, x))) == 3
    assert CompanionMatrix(Poly(x**3 + 2*x**2 + 3*x + 4, x)).det() == -4
    assert det(CompanionMatrix(Poly(x, x))) == 0
    assert det(CompanionMatrix(Poly(x**3 + 2*x**2 + 3*x, x))) == 0
    assert det(CompanionMatrix(Poly(x**2 + Rational(1, 2)*x + Rational(3, 4), x))) == Rational(3, 4)
    assert det(CompanionMatrix(Poly(x**50 + 3*x**49 + 7, x))) == 7
    assert det(CompanionMatrix(Poly(x**51 + 3*x**50 + 7, x))) == -7

    for p, expected in [
        (Poly([1, c0], x), -c0),
        (Poly([1, c1, c0], x), c0),
        (Poly([1, c2, c1, c0], x), -c0),
    ]:
        C = CompanionMatrix(p)
        assert Determinant(C).doit() == expected
        assert det(C) == expected
        assert C.det() == expected
        assert Determinant(C).doit() == C.as_explicit().det()
