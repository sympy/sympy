from sympy.matrices.expressions import MatrixSymbol, MatAdd
from sympy.matrices import eye, ImmutableMatrix
from sympy import Basic

X = MatrixSymbol('X', 2, 2)
Y = MatrixSymbol('Y', 2, 2)

def test_sort_key():
    assert MatAdd(Y, X).doit().args == (X, Y)


def test_matadd_sympify():
    assert isinstance(MatAdd(eye(1), eye(1)).args[0], Basic)


def test_matadd_of_matrices():
    assert MatAdd(eye(2), 4*eye(2), eye(2)).doit() == ImmutableMatrix(6*eye(2))
