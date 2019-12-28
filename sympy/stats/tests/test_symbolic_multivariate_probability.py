from sympy.utilities.pytest import raises

from sympy.stats import E, Expectation, Normal, Variance

from sympy import symbols, MatrixSymbol, Matrix, ZeroMatrix, ShapeError

from sympy.stats.rv import RandomSymbol, RandomMatrixSymbol
from sympy.stats.symbolic_multivariate_probability import ExpectationMatrix, VarianceMatrix

i, j, k = symbols("i j k")

A = MatrixSymbol("A", k, k)
B = MatrixSymbol("B", k, k)
C = MatrixSymbol("C", k, k)
D = MatrixSymbol("D", k, k)

a = MatrixSymbol("a", k, 1)
b = MatrixSymbol("b", k, 1)

A2 = MatrixSymbol("A", 2, 2)
B2 = MatrixSymbol("B", 2, 2)
C2 = MatrixSymbol("C", 2, 2)

X = RandomMatrixSymbol("X", k, 1)
Y = RandomMatrixSymbol("Y", k, 1)

R = RandomMatrixSymbol("R", k, k)

X2 = RandomMatrixSymbol("X", 2, 1)
Y2 = RandomMatrixSymbol("Y", 2, 1)

normal = Normal("normal", 0, 1)

m1 = Matrix([
    [1, j],
    [normal, 0]
])


def test_multivariate_expectation():
    expr = Expectation(a)
    assert expr == Expectation(a) == ExpectationMatrix(a)
    assert expr.doit() == a
    assert isinstance(expr, ExpectationMatrix)

    expr = Expectation(X)
    assert expr == Expectation(X) == ExpectationMatrix(X)
    assert expr.shape == (k, 1)
    assert expr.rows == k
    assert expr.cols == 1
    assert isinstance(expr, ExpectationMatrix)

    expr = Expectation(A*X + b)
    assert expr == ExpectationMatrix(A*X + b)
    assert expr.doit() == A*ExpectationMatrix(X) + b
    assert isinstance(expr, ExpectationMatrix)
    assert expr.shape == (k, 1)

    expr: ExpectationMatrix = Expectation(m1*X2)
    assert expr.doit() == expr

    expr = Expectation(A2*m1*B2*X2)
    assert expr.args[0].args == (A2, m1, B2, X2)
    assert expr.doit() == A2*ExpectationMatrix(m1*B2*X2)


def test_multivariate_variance():
    raises(ShapeError, lambda: Variance(A))

    expr = Variance(a)  # type: VarianceMatrix
    assert expr == Variance(a) == VarianceMatrix(a)
    assert expr.doit() == ZeroMatrix(k, k)
    expr = Variance(a.T)
    assert expr == Variance(a.T) == VarianceMatrix(a.T)
    assert expr.doit() == ZeroMatrix(k, k)

    expr = Variance(X)
    assert expr == Variance(X) == VarianceMatrix(X)
    assert expr.shape == (k, k)
    assert expr.rows == k
    assert expr.cols == k
    assert isinstance(expr, VarianceMatrix)

    expr: VarianceMatrix = Variance(A*X)
    assert expr == VarianceMatrix(A*X)
    assert expr.doit() == A*VarianceMatrix(X)*A.T
    assert isinstance(expr, VarianceMatrix)
    assert expr.shape == (k, k)

    expr: VarianceMatrix = Variance(A*B*X)
    assert expr.doit() == A*B*VarianceMatrix(X)*B.T*A.T

    expr: VarianceMatrix = Variance(A*B*X)
    assert expr.doit() == A*B*VarianceMatrix(X)*B.T*A.T

    expr: VarianceMatrix = Variance(m1*X2)
    assert expr.doit() == expr

    expr = Variance(A2*m1*B2*X2)
    assert expr.args[0].args == (A2, m1, B2, X2)
    assert expr.doit() == expr
