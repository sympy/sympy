"""
Some examples have been taken from:

http://www.math.uwaterloo.ca/~hwolkowi//matrixcookbook.pdf
"""
from sympy import MatrixSymbol, Inverse, symbols, Determinant, Trace

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


def test_matrix_derivative_trivial_cases():
    # Cookbook example 33:
    assert X._eval_derivative(A) == 0


def test_matrix_derivative_vectors_and_scalars():

    # Cookbook example 69:
    expr = x.T*a
    assert expr._eval_derivative(x) == a
    expr = a.T*x
    assert expr._eval_derivative(x) == a

    # Cookbook example 70:
    expr = a.T*X*b
    assert expr._eval_derivative(X) == a*b.T

    # Cookbook example 71:
    expr = a.T*X.T*b
    assert expr._eval_derivative(X) == b*a.T

    # Cookbook example 77:
    expr = b.T*X.T*X*c
    assert expr._eval_derivative(X) == X*b*c.T + X*c*b.T

    # Cookbook example 78:
    expr = (B*x + b).T*C*(D*x + d)
    assert expr._eval_derivative(x) == B.T*C*(D*x + d) + D.T*C.T*(B*x + b)

    # Cookbook example 81:
    expr = x.T*B*x
    assert expr._eval_derivative(x) == B*x + B.T*x

    # Cookbook example 82:
    expr = b.T*X.T*D*X*c
    assert expr._eval_derivative(X) == D.T*X*b*c.T + D*X*c*b.T

    # Cookbook example 83:
    expr = (X*b + c).T*D*(X*b + c)
    # TODO: parsing error
    # assert expr._eval_derivative(X) == D*X*b*b.T + D*c*b.T + D.T*X*b*b.T + D.T*c*b.T
    # otherwise: (D + D.T)*(X*b + c)*b.T


def test_matrix_derivative_with_inverse():

    # Cookbook example 61:
    expr = a.T*Inverse(X)*b
    assert expr._eval_derivative(X) == -Inverse(X).T*a*b.T*Inverse(X).T

    # Cookbook example 63:
    expr = Trace(A*Inverse(X)*B)
    assert expr._eval_derivative(X) == -(X**(-1)*B*A*X**(-1)).T

    # Cookbook example 64:
    expr = Trace(Inverse(X + A))
    assert expr._eval_derivative(X) == -(Inverse(X + A)*Inverse(X + A)).T
