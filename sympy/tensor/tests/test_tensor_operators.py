from sympy.utilities.pytest import raises

from sympy.tensor.toperators import PartialDerivative
from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorHead, \
    tensor_heads
from sympy import symbols, diag
from sympy import Array


L = TensorIndexType("L")
i, j, k = tensor_indices("i j k", L)
i0 = tensor_indices("i0", L)
L_0 = tensor_indices("L_0", L)

A, B, C, D = tensor_heads("A B C D", [L])

H = TensorHead("H", [L, L])


def test_invalid_partial_derivative_valence():
    raises(ValueError, lambda: PartialDerivative(C(j), D(-j)))
    raises(ValueError, lambda: PartialDerivative(C(-j), D(j)))


def test_tensor_partial_deriv():
    # Test flatten:
    expr = PartialDerivative(PartialDerivative(A(i), A(j)), A(i))
    assert expr.expr == A(L_0)
    assert expr.variables == (A(j), A(L_0))

    expr1 = PartialDerivative(A(i), A(j))
    assert expr1.expr == A(i)
    assert expr1.variables == (A(j),)

    expr2 = A(i)*PartialDerivative(H(k, -i), A(j))
    assert expr2.get_indices() == [L_0, k, -L_0, -j]

    expr2b = A(i)*PartialDerivative(H(k, -i), A(-j))
    assert expr2b.get_indices() == [L_0, k, -L_0, j]

    expr3 = A(i)*PartialDerivative(B(k)*C(-i) + 3*H(k, -i), A(j))
    assert expr3.get_indices() == [L_0, k, -L_0, -j]

    expr4 = (A(i) + B(i))*PartialDerivative(C(j), D(j))
    assert expr4.get_indices() == [i, L_0, -L_0]

    expr4b = (A(i) + B(i))*PartialDerivative(C(-j), D(-j))
    assert expr4b.get_indices() == [i, -L_0, L_0]

    expr5 = (A(i) + B(i))*PartialDerivative(C(-i), D(j))
    assert expr5.get_indices() == [L_0, -L_0, -j]


def test_replace_arrays_partial_derivative():
    x, y, z, t = symbols("x y z t")

    # d(A^i)/d(A_j) = d(g^ik A_k)/d(A_j) = g^ik delta_jk
    expr = PartialDerivative(A(i), A(-j))
    assert expr.get_free_indices() == [i, j]
    assert expr.get_indices() == [i, j]
    assert expr.replace_with_arrays({A(i): [x, y], L: diag(1, 1)}, [i, j]) == Array([[1, 0], [0, 1]])
    assert expr.replace_with_arrays({A(i): [x, y], L: diag(1, -1)}, [i, j]) == Array([[1, 0], [0, -1]])
    assert expr.replace_with_arrays({A(-i): [x, y], L: diag(1, 1)}, [i, j]) == Array([[1, 0], [0, 1]])
    assert expr.replace_with_arrays({A(-i): [x, y], L: diag(1, -1)}, [i, j]) == Array([[1, 0], [0, -1]])

    expr = PartialDerivative(A(i), A(j))
    assert expr.get_free_indices() == [i, -j]
    assert expr.get_indices() == [i, -j]
    assert expr.replace_with_arrays({A(i): [x, y]}, [i, -j]) == Array([[1, 0], [0, 1]])
    assert expr.replace_with_arrays({A(i): [x, y], L: diag(1, 1)}, [i, -j]) == Array([[1, 0], [0, 1]])
    assert expr.replace_with_arrays({A(i): [x, y], L: diag(1, -1)}, [i, -j]) == Array([[1, 0], [0, 1]])
    assert expr.replace_with_arrays({A(-i): [x, y], L: diag(1, 1)}, [i, -j]) == Array([[1, 0], [0, 1]])
    assert expr.replace_with_arrays({A(-i): [x, y], L: diag(1, -1)}, [i, -j]) == Array([[1, 0], [0, 1]])

    expr = PartialDerivative(A(-i), A(-j))
    expr.get_free_indices() == [-i, j]
    expr.get_indices() == [-i, j]
    assert expr.replace_with_arrays({A(-i): [x, y]}, [-i, j]) == Array([[1, 0], [0, 1]])
    assert expr.replace_with_arrays({A(-i): [x, y], L: diag(1, 1)}, [-i, j]) == Array([[1, 0], [0, 1]])
    assert expr.replace_with_arrays({A(-i): [x, y], L: diag(1, -1)}, [-i, j]) == Array([[1, 0], [0, 1]])
    assert expr.replace_with_arrays({A(i): [x, y], L: diag(1, 1)}, [-i, j]) == Array([[1, 0], [0, 1]])
    assert expr.replace_with_arrays({A(i): [x, y], L: diag(1, -1)}, [-i, j]) == Array([[1, 0], [0, 1]])

    expr = PartialDerivative(A(i), A(i))
    assert expr.get_free_indices() == []
    assert expr.get_indices() == [L_0, -L_0]
    assert expr.replace_with_arrays({A(i): [x, y], L: diag(1, 1)}, []) == 2
    assert expr.replace_with_arrays({A(i): [x, y], L: diag(1, -1)}, []) == 2

    expr = PartialDerivative(A(-i), A(-i))
    assert expr.get_free_indices() == []
    assert expr.get_indices() == [-L_0, L_0]
    assert expr.replace_with_arrays({A(i): [x, y], L: diag(1, 1)}, []) == 2
    assert expr.replace_with_arrays({A(i): [x, y], L: diag(1, -1)}, []) == 2
