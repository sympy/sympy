from sympy.tensor.operators import PartialDerivative
from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead

def test_tensor_partial_deriv():
    L = TensorIndexType("L")
    i, j, k = tensor_indices("i j k", L)
    i0 = tensor_indices("i0", L)
    L_0 = tensor_indices("L_0", L)

    A, B, C, D = tensorhead("A B C D", [L], [[1]])

    H = tensorhead("H", [L, L], [[1], [1]])

    # Test flatten:
    expr = PartialDerivative(PartialDerivative(A(i), A(j)), A(-i))
    assert expr.expr == A(L_0)
    assert expr.variables == (A(j), A(-L_0))

    expr1 = PartialDerivative(A(i), A(j))
    assert expr1.expr == A(i)
    assert expr1.variables == (A(j),)

    expr2 = A(i)*PartialDerivative(H(k, -i), A(j))
    assert expr2.get_indices() == [L_0, k, -L_0, j]

    expr3 = A(i)*PartialDerivative(B(k)*C(-i) + 3*H(k, -i), A(j))
    assert expr3.get_indices() == [L_0, k, -L_0, j]

    expr4 = (A(i) + B(i))*PartialDerivative(C(-j), D(j))
    assert expr4.get_indices() == [i, -L_0, L_0]

    expr5 = (A(i) + B(i))*PartialDerivative(C(-i), D(j))
    assert expr5.get_indices() == [L_0, -L_0, j]
