"""
Some examples have been taken from:

http://www.math.uwaterloo.ca/~hwolkowi//matrixcookbook.pdf
"""
from sympy import MatrixSymbol, Inverse, symbols, Determinant, Trace, Derivative
from sympy import MatAdd, Identity

k = symbols("k")

X = MatrixSymbol("X", k, k)
x = MatrixSymbol("x", k, 1)

A = MatrixSymbol("A", k, k)
B = MatrixSymbol("B", k, k)
C = MatrixSymbol("C", k, k)
D = MatrixSymbol("D", k, k)

a = MatrixSymbol("a", k, 1)
b = MatrixSymbol("b", k, 1)
c = MatrixSymbol("c", k, 1)
d = MatrixSymbol("d", k, 1)


def test_matrix_derivative_non_matrix_result():
    # This is a 4-dimensional array:
    assert A.diff(A) == Derivative(A, A)
    assert A.T.diff(A) == Derivative(A.T, A)
    assert (2*A).diff(A) == Derivative(2*A, A)
    assert MatAdd(A, A).diff(A) == Derivative(MatAdd(A, A), A)
    assert (A + B).diff(A) == Derivative(A + B, A)  # TODO: `B` can be removed.


def test_matrix_derivative_trivial_cases():
    # Cookbook example 33:
    assert X.diff(A) == 0


def test_matrix_derivative_vectors_and_scalars():

    # Cookbook example 69:
    expr = x.T*a
    assert expr.diff(x) == a
    expr = a.T*x
    assert expr.diff(x) == a

    # Cookbook example 70:
    expr = a.T*X*b
    assert expr.diff(X) == a*b.T

    # Cookbook example 71:
    expr = a.T*X.T*b
    assert expr.diff(X) == b*a.T

    # Cookbook example 77:
    expr = b.T*X.T*X*c
    assert expr.diff(X) == X*b*c.T + X*c*b.T

    # Cookbook example 78:
    expr = (B*x + b).T*C*(D*x + d)
    assert expr.diff(x) == B.T*C*(D*x + d) + D.T*C.T*(B*x + b)

    # Cookbook example 81:
    expr = x.T*B*x
    assert expr.diff(x) == B*x + B.T*x

    # Cookbook example 82:
    expr = b.T*X.T*D*X*c
    assert expr.diff(X) == D.T*X*b*c.T + D*X*c*b.T

    # Cookbook example 83:
    expr = (X*b + c).T*D*(X*b + c)
    assert expr.diff(X) == D*(X*b + c)*b.T + D.T*(X*b + c)*b.T


def test_matrix_derivative_with_inverse():

    # Cookbook example 61:
    expr = a.T*Inverse(X)*b
    assert expr.diff(X) == -Inverse(X).T*a*b.T*Inverse(X).T

    # Cookbook example 63:
    expr = Trace(A*Inverse(X)*B)
    assert expr.diff(X) == -(X**(-1)*B*A*X**(-1)).T

    # Cookbook example 64:
    expr = Trace(Inverse(X + A))
    assert expr.diff(X) == -(Inverse(X + A)*Inverse(X + A)).T


def test_matrix_derivatives_of_traces():

    ## First order:

    # Cookbook example 99:
    expr = Trace(X)
    assert expr.diff(X) == Identity(k)

    # Cookbook example 100:
    expr = Trace(X*A)
    assert expr.diff(X) == A.T

    # Cookbook example 101:
    expr = Trace(A*X*B)
    assert expr.diff(X) == A.T*B.T

    # Cookbook example 102:
    expr = Trace(A*X.T*B)
    assert expr.diff(X) == B*A

    # Cookbook example 103:
    expr = Trace(X.T*A)
    assert expr.diff(X) == A

    # Cookbook example 104:
    expr = Trace(A*X.T)
    assert expr.diff(X) == A

    # Cookbook example 105:
    # TODO: TensorProduct is not supported
    #expr = Trace(TensorProduct(A, X))
    #assert expr.diff(X) == Trace(A)*Identity(k)

    ## Second order:

    # Cookbook example 106:
    expr = Trace(X**2)
    assert expr.diff(X) == 2*X.T

    # Cookbook example 107:
    expr = Trace(X**2*B)
    # TODO: wrong result
    #assert expr.diff(X) == (X*B + B*X).T
    expr = Trace(X*X*B)
    assert expr.diff(X) == (X*B + B*X).T

    # Cookbook example 108:
    expr = Trace(X.T*B*X)
    assert expr.diff(X) == B*X + B.T*X

    # Cookbook example 109:
    expr = Trace(B*X*X.T)
    assert expr.diff(X) == B*X + B.T*X

    # Cookbook example 110:
    expr = Trace(X*X.T*B)
    assert expr.diff(X) == B*X + B.T*X

    # Cookbook example 111:
    expr = Trace(X*B*X.T)
    assert expr.diff(X) == X*B.T + X*B

    # Cookbook example 112:
    expr = Trace(B*X.T*X)
    assert expr.diff(X) == X*B.T + X*B

    # Cookbook example 113:
    expr = Trace(X.T*X*B)
    assert expr.diff(X) == X*B.T + X*B

    # Cookbook example 114:
    expr = Trace(A*X*B*X)
    assert expr.diff(X) == A.T*X.T*B.T + B.T*X.T*A.T

    # Cookbook example 115:
    expr = Trace(X.T*X)
    assert expr.diff(X) == 2*X
    expr = Trace(X*X.T)
    assert expr.diff(X) == 2*X

    # Cookbook example 116:
    expr = Trace(B.T*X.T*C*X*B)
    assert expr.diff(X) == C.T*X*B*B.T + C*X*B*B.T

    # Cookbook example 117:
    expr = Trace(X.T*B*X*C)
    assert expr.diff(X) == B*X*C + B.T*X*C.T

    # Cookbook example 118:
    expr = Trace(A*X*B*X.T*C)
    assert expr.diff(X) == A.T*C.T*X*B.T + C*A*X*B

    # Cookbook example 119:
    expr = Trace((A*X*B + C)*(A*X*B + C).T)
    assert expr.diff(X) == 2*A.T*(A*X*B + C)*B.T

    # Cookbook example 120:
    # TODO: no support for TensorProduct.
    # expr = Trace(TensorProduct(X, X))
    # expr = Trace(X)*Trace(X)
    # expr.diff(X) == 2*Trace(X)*Identity(k)
