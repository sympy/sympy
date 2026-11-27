from __future__ import annotations
from sympy.core.expr import unchanged
from sympy.core.symbol import Symbol, symbols
from sympy.matrices.immutable import ImmutableDenseMatrix
from sympy.matrices.expressions.trace import Trace
from sympy.matrices.expressions.determinant import Determinant
from sympy.matrices.expressions.companion import CompanionMatrix
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


def test_companion_trace_and_determinant():
    x = Symbol('x')
    for coeffs, expected_trace, expected_det in [
        ([1, 2, 3, 4], -2, -4),   # x**3 + 2*x**2 + 3*x + 4
        ([1, 5], -5, -5),         # x + 5
        ([1, 0, 7, 0, 9], 0, 9),  # x**4 + 7*x**2 + 9
        ([1, 0, 0], 0, 0),        # x**2
    ]:
        p = Poly(coeffs, x)
        C = CompanionMatrix(p)
        assert Trace(C).doit() == expected_trace
        assert Determinant(C).doit() == expected_det


def test_companion_trace_and_determinant_match_explicit():
    x = Symbol('x')
    for coeffs in [[1, 3, 2], [1, 0, 0, 7], [1, 4, 2, 0, 9], [1, 1]]:
        p = Poly(coeffs, x)
        C = CompanionMatrix(p)
        explicit = ImmutableDenseMatrix.companion(p)
        assert Trace(C).doit() == explicit.trace()
        assert Determinant(C).doit() == explicit.det()
