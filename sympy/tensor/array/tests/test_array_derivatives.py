from __future__ import annotations
from sympy.core.function import Derivative, Function
from sympy.core.symbol import symbols
from sympy.matrices.dense import Matrix
from sympy.matrices.expressions.matexpr import MatrixSymbol
from sympy.matrices.expressions.trace import Trace
from sympy.tensor.array import Array
from sympy.tensor.array.ndim_array import NDimArray
from sympy.matrices.matrixbase import MatrixBase
from sympy.tensor.array.array_derivatives import ArrayDerivative

x, y, z, t = symbols("x y z t")

m = Matrix([[x, y], [z, t]])

M = MatrixSymbol("M", 3, 2)
N = MatrixSymbol("N", 4, 3)


def test_array_derivative_construction():

    d = ArrayDerivative(x, m, evaluate=False)
    assert d.shape == (2, 2)
    expr = d.doit()
    assert isinstance(expr, MatrixBase)
    assert expr.shape == (2, 2)

    d = ArrayDerivative(m, m, evaluate=False)
    assert d.shape == (2, 2, 2, 2)
    expr = d.doit()
    assert isinstance(expr, NDimArray)
    assert expr.shape == (2, 2, 2, 2)

    d = ArrayDerivative(m, x, evaluate=False)
    assert d.shape == (2, 2)
    expr = d.doit()
    assert isinstance(expr, MatrixBase)
    assert expr.shape == (2, 2)

    d = ArrayDerivative(M, N, evaluate=False)
    assert d.shape == (4, 3, 3, 2)
    expr = d.doit()
    assert isinstance(expr, ArrayDerivative)
    assert expr.shape == (4, 3, 3, 2)

    d = ArrayDerivative(M, (N, 2), evaluate=False)
    assert d.shape == (4, 3, 4, 3, 3, 2)
    expr = d.doit()
    assert isinstance(expr, ArrayDerivative)
    assert expr.shape == (4, 3, 4, 3, 3, 2)

    d = ArrayDerivative(M.as_explicit(), (N.as_explicit(), 2), evaluate=False)
    assert d.doit().shape == (4, 3, 4, 3, 3, 2)
    expr = d.doit()
    assert isinstance(expr, NDimArray)
    assert expr.shape == (4, 3, 4, 3, 3, 2)


def test_derivative_routing_to_array_derivative():
    # A directly constructed Derivative with matrix/array operands must be
    # routed to ArrayDerivative regardless of the construction path, so
    # that Derivative(...) with evaluate=True and Derivative(...).doit()
    # agree with the equivalent .diff() call.
    A = MatrixSymbol("A", 3, 3)
    X = MatrixSymbol("X", 3, 3)
    v = MatrixSymbol("v", 3, 1)
    f = Function("f")
    s = symbols("s")

    # unevaluated construction is re-dispatched to the subclass:
    assert isinstance(Derivative(Trace(A*X), X), ArrayDerivative)
    assert isinstance(Derivative(A*X, X), ArrayDerivative)
    assert isinstance(Derivative(x, m, evaluate=False), ArrayDerivative)

    # scalar-by-matrixsymbol:
    assert Derivative(Trace(A*X), X).doit() == Trace(A*X).diff(X)
    assert Derivative(Trace(A*X), X, evaluate=True) == Trace(A*X).diff(X)
    assert Derivative(Trace(A*X), X).doit() == A.T

    # matrixsymbol-by-matrixsymbol:
    assert Derivative(A*X, X).doit() == (A*X).diff(X)

    # array-by-scalar:
    assert Derivative(Array([f(s), s]), s).doit() == Array([f(s).diff(s), 1])
    assert Derivative(Array([[s, t], [s*t, s**2]]), s).doit() == \
        Array([[1, 0], [t, 2*s]])

    # scalar-by-matrix: the shape of the result is that of the matrix:
    u, w = symbols("u w")
    assert Derivative(x, Matrix([[y, z], [u, w]])).doit() == Matrix.zeros(2, 2)

    # matrix-by-matrix:
    p, q = symbols("p q")
    assert Derivative(Matrix([[y*p, z+q]]), Matrix([[p, q]])).doit() == \
        Matrix([[y*p, z+q]]).diff(Matrix([[p, q]]))

    # n-times (count > 1) form must match repeated .diff():
    assert Derivative(Trace(X*X), X, 2).doit() == Trace(X*X).diff(X).diff(X)
    assert Derivative(v.T*A*v, v, 2).doit() == (v.T*A*v).diff(v).diff(v)
    assert Derivative(v.T*A*v, v, 2).doit() == A + A.T
