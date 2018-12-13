from sympy.matrices.expressions.applyfunc import ElementwiseApplyFunction
from sympy import (Matrix, Lambda, MatrixBase, MatrixSymbol, exp, symbols, MatMul, sin)
from sympy.utilities.pytest import raises
from sympy.matrices.common import ShapeError


X = MatrixSymbol("X", 3, 3)
Y = MatrixSymbol("Y", 3, 3)
Xd = X.as_explicit()

x, y, z, t = symbols("x y z t")


def test_applyfunc_matrix():
    double = Lambda(x, x**2)

    expr = ElementwiseApplyFunction(double, Xd)
    assert isinstance(expr, ElementwiseApplyFunction)
    assert expr.doit() == Xd.applyfunc(lambda x: x**2)
    assert expr.shape == (3, 3)

    expr = ElementwiseApplyFunction(double, X)
    assert isinstance(expr, ElementwiseApplyFunction)
    assert isinstance(expr.doit(), ElementwiseApplyFunction)
    assert expr == X.applyfunc(double)

    expr = ElementwiseApplyFunction(exp, X*Y)
    assert expr.expr == X*Y
    assert expr.function == exp
    assert expr == (X*Y).applyfunc(exp)

    assert isinstance(X*expr, MatMul)
    assert (X*expr).shape == (3, 3)
    Z = MatrixSymbol("Z", 2, 3)
    assert (Z*expr).shape == (2, 3)

    expr = ElementwiseApplyFunction(exp, Z.T)*ElementwiseApplyFunction(exp, Z)
    assert expr.shape == (3, 3)
    expr = ElementwiseApplyFunction(exp, Z)*ElementwiseApplyFunction(exp, Z.T)
    assert expr.shape == (2, 2)

    raises(ShapeError, lambda: ElementwiseApplyFunction(exp, Z)*ElementwiseApplyFunction(exp, Z))

    M = Matrix([[x, y], [z, t]])
    expr = ElementwiseApplyFunction(sin, M)
    assert isinstance(expr, ElementwiseApplyFunction)
    assert expr.function == sin
    assert expr.expr == M
    assert expr.doit() == M.applyfunc(sin)
    assert expr.doit() == Matrix([[sin(x), sin(y)], [sin(z), sin(t)]])
