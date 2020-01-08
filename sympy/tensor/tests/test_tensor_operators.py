from sympy.testing.pytest import raises

from sympy.tensor.toperators import PartialDerivative
from sympy.tensor.tensor import (TensorIndexType,
                                 tensor_indices,
                                 TensorHead, tensor_heads,
                                 TensMul)
from sympy import symbols, diag
from sympy import Array, Rational

from sympy import sympify
from random import randint
from sympy import I


L = TensorIndexType("L")
i, j, k, m = tensor_indices("i j k m", L)
i0 = tensor_indices("i0", L)
L_0, L_1 = tensor_indices("L_0 L_1", L)

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

    expr = PartialDerivative(H(i, j) + H(j, i), A(i))
    assert expr.get_indices() == [L_0, j, -L_0]
    assert expr.get_free_indices() == [j]

    expr = PartialDerivative(H(i, j) + H(j, i), A(k))*B(-i)
    assert expr.get_indices() == [L_0, j, -k, -L_0]
    assert expr.get_free_indices() == [j, -k]

    expr = PartialDerivative(A(i)*(H(-i, j) + H(j, -i)), A(j))
    assert expr.get_indices() == [L_0, -L_0, L_1, -L_1]
    assert expr.get_free_indices() == []

    expr = A(j)*A(-j) + expr
    assert expr.get_indices() == [L_0, -L_0, L_1, -L_1]
    assert expr.get_free_indices() == []

    expr = A(i)*(B(j)*PartialDerivative(C(-j), D(i)) + C(j)*PartialDerivative(D(-j), B(i)))
    assert expr.get_indices() == [L_0, L_1, -L_1, -L_0]
    assert expr.get_free_indices() == []

    expr = A(i)*PartialDerivative(C(-j), D(i))
    assert expr.get_indices() == [L_0, -j, -L_0]
    assert expr.get_free_indices() == [-j]


def test_expand_partial_derivative_sum_rule():
    tau = symbols("tau")

    # check sum rule for D(tensor, symbol)
    expr1aa = PartialDerivative(A(i), tau)

    assert expr1aa._expand_partial_derivative() == PartialDerivative(A(i), tau)

    expr1ab = PartialDerivative(A(i) + B(i), tau)

    # This fails in Python2.7
    assert (expr1ab._expand_partial_derivative() ==
            PartialDerivative(A(i), tau) +
            PartialDerivative(B(i), tau))

    expr1ac = PartialDerivative(A(i) + B(i) + C(i), tau)

    # This fails in Python2.7
    assert (expr1ac._expand_partial_derivative() ==
            PartialDerivative(A(i), tau) +
            PartialDerivative(B(i), tau) +
            PartialDerivative(C(i), tau))

    # check sum rule for D(tensor, D(j))
    expr1ba = PartialDerivative(A(i), D(j))

    assert expr1ba._expand_partial_derivative() ==\
        PartialDerivative(A(i), D(j))
    expr1bb = PartialDerivative(A(i) + B(i), D(j))

    # This fails in Python 2.7
    assert (expr1bb._expand_partial_derivative() ==
            PartialDerivative(A(i), D(j)) +
            PartialDerivative(B(i), D(j)))

    expr1bc = PartialDerivative(A(i) + B(i) + C(i), D(j))
    assert expr1bc._expand_partial_derivative() ==\
        PartialDerivative(A(i), D(j))\
        + PartialDerivative(B(i), D(j))\
        + PartialDerivative(C(i), D(j))

    # check sum rule for D(tensor, H(j, k))
    expr1ca = PartialDerivative(A(i), H(j, k))
    assert expr1ca._expand_partial_derivative() ==\
        PartialDerivative(A(i), H(j, k))
    expr1cb = PartialDerivative(A(i) + B(i), H(j, k))
    assert (expr1cb._expand_partial_derivative() ==
            PartialDerivative(A(i), H(j, k))
            + PartialDerivative(B(i), H(j, k)))
    expr1cc = PartialDerivative(A(i) + B(i) + C(i), H(j, k))
    assert (expr1cc._expand_partial_derivative() ==
            PartialDerivative(A(i), H(j, k))
            + PartialDerivative(B(i), H(j, k))
            + PartialDerivative(C(i), H(j, k)))

    # check sum rule for D(D(tensor, D(j)), H(k, m))
    expr1da = PartialDerivative(A(i), (D(j), H(k, m)))
    assert expr1da._expand_partial_derivative() ==\
        PartialDerivative(A(i), (D(j), H(k, m)))
    expr1db = PartialDerivative(A(i) + B(i), (D(j), H(k, m)))
    assert expr1db._expand_partial_derivative() ==\
        PartialDerivative(A(i), (D(j), H(k, m)))\
        + PartialDerivative(B(i), (D(j), H(k, m)))
    expr1dc = PartialDerivative(A(i) + B(i) + C(i), (D(j), H(k, m)))
    assert expr1dc._expand_partial_derivative() ==\
        PartialDerivative(A(i), (D(j), H(k, m)))\
        + PartialDerivative(B(i), (D(j), H(k, m)))\
        + PartialDerivative(C(i), (D(j), H(k, m)))


def test_expand_partial_derivative_constant_factor_rule():
    for count in range(10):
        pos_random_int1 = sympify(randint(0, 1000))
        pos_random_int2 = sympify(randint(0, 1000))
        neg_random_int = sympify(randint(-1000, -1))

        c1 = Rational(pos_random_int1, pos_random_int2)
        c2 = Rational(neg_random_int, pos_random_int2)
        c3 = Rational(pos_random_int1, neg_random_int)

        expr2a = PartialDerivative(pos_random_int1*A(i), D(j))
        assert expr2a._expand_partial_derivative() ==\
            pos_random_int1*PartialDerivative(A(i), D(j))

        expr2b = PartialDerivative(neg_random_int*A(i), D(j))
        assert expr2b._expand_partial_derivative() ==\
            neg_random_int*PartialDerivative(A(i), D(j))

        expr2ca = PartialDerivative(c1*A(i), D(j))
        assert expr2ca._expand_partial_derivative() ==\
            c1*PartialDerivative(A(i), D(j))

        expr2cb = PartialDerivative(c2*A(i), D(j))
        assert expr2cb._expand_partial_derivative() ==\
            c2*PartialDerivative(A(i), D(j))

        expr2cc = PartialDerivative(c3*A(i), D(j))
        assert expr2cc._expand_partial_derivative() ==\
            c3*PartialDerivative(A(i), D(j))


def test_expand_partial_derivative_full_linearity():
    for count in range(10):
        pos_random_int1 = sympify(randint(0, 1000))
        pos_random_int2 = sympify(randint(0, 1000))
        neg_random_int = sympify(randint(-1000, -1))

        c1 = Rational(pos_random_int1, pos_random_int2)
        c2 = Rational(neg_random_int, pos_random_int2)
        c3 = Rational(pos_random_int1, neg_random_int)

        # check full linearity

        expr3a = PartialDerivative(pos_random_int1*A(i) + pos_random_int2*B(i), D(j))
        assert expr3a._expand_partial_derivative() ==\
            pos_random_int1*PartialDerivative(A(i), D(j))\
            + pos_random_int2*PartialDerivative(B(i), D(j))

        expr3b = PartialDerivative(pos_random_int1*A(i) + neg_random_int*B(i), D(j))
        assert expr3b._expand_partial_derivative() ==\
            pos_random_int1*PartialDerivative(A(i), D(j))\
            + neg_random_int*PartialDerivative(B(i), D(j))

        expr3c = PartialDerivative(neg_random_int*A(i) + pos_random_int2*B(i), D(j))
        assert expr3c._expand_partial_derivative() ==\
            neg_random_int*PartialDerivative(A(i), D(j))\
            + pos_random_int2*PartialDerivative(B(i), D(j))

        expr3d = PartialDerivative(c1*A(i) + c2*B(i), D(j))
        assert expr3d._expand_partial_derivative() ==\
            c1*PartialDerivative(A(i), D(j))\
            + c2*PartialDerivative(B(i), D(j))

        expr3e = PartialDerivative(c2*A(i) + c1*B(i), D(j))
        assert expr3e._expand_partial_derivative() ==\
            c2*PartialDerivative(A(i), D(j))\
            + c1*PartialDerivative(B(i), D(j))

        expr3f = PartialDerivative(c2*A(i) + c3*B(i), D(j))
        assert expr3f._expand_partial_derivative() ==\
            c2*PartialDerivative(A(i), D(j))\
            + c3*PartialDerivative(B(i), D(j))

        expr3g = PartialDerivative(c3*A(i) + c2*B(i), D(j))
        assert expr3g._expand_partial_derivative() ==\
            c3*PartialDerivative(A(i), D(j))\
            + c2*PartialDerivative(B(i), D(j))

        expr3h = PartialDerivative(c3*A(i) + c1*B(i), D(j))
        assert expr3h._expand_partial_derivative() ==\
            c3*PartialDerivative(A(i), D(j))\
            + c1*PartialDerivative(B(i), D(j))

        expr3i = PartialDerivative(c1*A(i) + c3*B(i), D(j))
        assert expr3i._expand_partial_derivative() ==\
            c1*PartialDerivative(A(i), D(j))\
            + c3*PartialDerivative(B(i), D(j))


def test_expand_partial_derivative_product_rule():
    # check product rule
    expr4a = PartialDerivative(A(i)*B(j), D(k))
    assert expr4a._expand_partial_derivative() ==\
        A(i)*PartialDerivative(B(j), D(k))\
        + PartialDerivative(A(i), D(k))*B(j)

    expr4b = PartialDerivative(A(i)*B(j)*C(k), D(m))
    assert expr4b._expand_partial_derivative() ==\
        PartialDerivative(A(i), D(m))*B(j)*C(k)\
        + A(i)*PartialDerivative(B(j), D(m))*C(k)\
        + A(i)*B(j)*PartialDerivative(C(k), D(m))

    expr4c = PartialDerivative(A(i)*B(j), C(k), D(m))
    assert expr4c._expand_partial_derivative() ==\
        PartialDerivative(A(i), C(k), D(m))*B(j)\
        + PartialDerivative(A(i), D(m))*PartialDerivative(B(j), C(k))\
        + PartialDerivative(A(i), C(k))*PartialDerivative(B(j), D(m))\
        + A(i)*PartialDerivative(B(j), C(k), D(m))
