from sympy.matrices.expressions.permutation import PermutationMatrix
from sympy.matrices import (Matrix, MatrixSymbol, MatMul, ImmutableMatrix,
        MatrixSymbol)
from sympy import Symbol

def test_basic():
    assert PermutationMatrix(0, 1).is_Identity
    assert PermutationMatrix((0, 1)).is_Identity
    assert Matrix(PermutationMatrix((1, 0))) == Matrix([[0, 1], [1, 0]])
    assert isinstance(PermutationMatrix((1, 0)) * MatrixSymbol('A', 2, 2),
                      MatMul)

def test_inputs():
    P = PermutationMatrix(1, 0)
    inputs = ((1, 0), [1, 0], Matrix([[1, 0]]), Matrix([[1], [0]]),
            ImmutableMatrix([[1, 0]]))
    assert all(PermutationMatrix(inp) == P for inp in inputs)

def test_matrix_interaction():
    P = PermutationMatrix((1, 2, 0))
    assert Matrix(P*Matrix(3, 3, range(9))) == Matrix([[3, 4, 5],
                                                       [6, 7, 8],
                                                       [0, 1, 2]])

def test_matrix_expr_as_input():
    n = Symbol('n')
    A = MatrixSymbol('A', 1, n)
    P = PermutationMatrix(A)
    assert P.shape == (n, n)

def test_identity():
    assert PermutationMatrix(0, 1).is_Identity
    assert not PermutationMatrix(1, 0).is_Identity
    n = Symbol('n')
    A = MatrixSymbol('A', 1, n)
    assert not PermutationMatrix(A).is_Identity
