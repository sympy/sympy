from __future__ import annotations
from sympy import tanh, expand
from sympy.concrete.summations import Sum
from sympy.core.function import Function, Lambda
from sympy.core.power import Pow
from sympy.core.symbol import symbols, Dummy
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy.matrices.expressions.matexpr import MatrixSymbol
from sympy.matrices.expressions.special import Identity
from sympy.tensor.array.expressions import ArrayElementwiseApplyFunc
from sympy.tensor.indexed import IndexedBase
from sympy.combinatorics import Permutation
from sympy.tensor.array.expressions.array_expressions import ArrayContraction, ArrayTensorProduct, \
    ArrayDiagonal, ArrayAdd, PermuteDims, ArrayElement, _array_tensor_product, _array_contraction, _array_diagonal, \
    _array_add, _permute_dims, ArraySymbol, OneArray
from sympy.tensor.array.expressions.from_array_to_matrix import convert_array_to_matrix
from sympy.tensor.array.expressions.from_indexed_to_array import convert_indexed_to_array, _convert_indexed_to_array
from sympy.testing.pytest import raises


A, B = symbols("A B", cls=IndexedBase)
i, j, k, l, m, n = symbols("i j k l m n")
d0, d1, d2, d3 = symbols("d0:4")

I = Identity(k)

M = MatrixSymbol("M", k, k)
N = MatrixSymbol("N", k, k)
P = MatrixSymbol("P", k, k)
Q = MatrixSymbol("Q", k, k)

a = MatrixSymbol("a", k, 1)
b = MatrixSymbol("b", k, 1)
c = MatrixSymbol("c", k, 1)
d = MatrixSymbol("d", k, 1)


def test_arrayexpr_convert_index_to_array_support_function():
    expr = M[i, j]
    assert _convert_indexed_to_array(expr) == (M, (i, j))
    expr = M[i, j]*N[k, l]
    assert _convert_indexed_to_array(expr) == (ArrayTensorProduct(M, N), (i, j, k, l))
    expr = M[i, j]*N[j, k]
    assert _convert_indexed_to_array(expr) == (ArrayDiagonal(ArrayTensorProduct(M, N), (1, 2)), (i, k, j))
    expr = Sum(M[i, j]*N[j, k], (j, 0, k-1))
    assert _convert_indexed_to_array(expr) == (ArrayContraction(ArrayTensorProduct(M, N), (1, 2)), (i, k))
    expr = M[i, j] + N[i, j]
    assert _convert_indexed_to_array(expr) == (ArrayAdd(M, N), (i, j))
    expr = M[i, j] + N[j, i]
    assert _convert_indexed_to_array(expr) == (ArrayAdd(M, PermuteDims(N, Permutation([1, 0]))), (i, j))
    expr = M[i, j] + M[j, i]
    assert _convert_indexed_to_array(expr) == (ArrayAdd(M, PermuteDims(M, Permutation([1, 0]))), (i, j))
    expr = (M*N*P)[i, j]
    assert _convert_indexed_to_array(expr) == (_array_contraction(ArrayTensorProduct(M, N, P), (1, 2), (3, 4)), (i, j))
    expr = expr.function  # Disregard summation in previous expression
    ret1, ret2 = _convert_indexed_to_array(expr)
    assert ret1 == ArrayDiagonal(ArrayTensorProduct(M, N, P), (1, 2), (3, 4))
    assert str(ret2) == "(i, j, _i_1, _i_2)"
    expr = KroneckerDelta(i, j)*M[i, k]
    assert _convert_indexed_to_array(expr) == (M, ({i, j}, k))
    expr = KroneckerDelta(i, j)*KroneckerDelta(j, k)*M[i, l]
    assert _convert_indexed_to_array(expr) == (M, ({i, j, k}, l))
    expr = KroneckerDelta(j, k)*(M[i, j]*N[k, l] + N[i, j]*M[k, l])
    assert _convert_indexed_to_array(expr) == (_array_diagonal(_array_add(
            ArrayTensorProduct(M, N),
            _permute_dims(ArrayTensorProduct(M, N), Permutation(0, 2)(1, 3))
        ), (1, 2)), (i, l, frozenset({j, k})))
    expr = KroneckerDelta(j, m)*KroneckerDelta(m, k)*(M[i, j]*N[k, l] + N[i, j]*M[k, l])
    assert _convert_indexed_to_array(expr) == (_array_diagonal(_array_add(
            ArrayTensorProduct(M, N),
            _permute_dims(ArrayTensorProduct(M, N), Permutation(0, 2)(1, 3))
        ), (1, 2)), (i, l, frozenset({j, m, k})))
    expr = KroneckerDelta(i, j)*KroneckerDelta(j, k)*KroneckerDelta(k,m)*M[i, 0]*KroneckerDelta(m, n)
    assert _convert_indexed_to_array(expr) == (M, ({i, j, k, m, n}, 0))
    expr = M[i, i]
    assert _convert_indexed_to_array(expr) == (ArrayDiagonal(M, (0, 1)), (i,))


def test_arrayexpr_convert_indexed_to_array_expression():

    s = Sum(A[i]*B[i], (i, 0, 3))
    cg = convert_indexed_to_array(s)
    assert cg == ArrayContraction(ArrayTensorProduct(A, B), (0, 1))

    expr = M*N
    result = ArrayContraction(ArrayTensorProduct(M, N), (1, 2))
    elem = expr[i, j]
    assert convert_indexed_to_array(elem) == result

    expr = M*N*M
    elem = expr[i, j]
    result = _array_contraction(_array_tensor_product(M, M, N), (1, 4), (2, 5))
    cg = convert_indexed_to_array(elem)
    assert cg == result

    cg = convert_indexed_to_array((M * N * P)[i, j])
    assert cg == _array_contraction(ArrayTensorProduct(M, N, P), (1, 2), (3, 4))

    cg = convert_indexed_to_array((M * N.T * P)[i, j])
    assert cg == _array_contraction(ArrayTensorProduct(M, N, P), (1, 3), (2, 4))

    expr = -2*M*N
    elem = expr[i, j]
    cg = convert_indexed_to_array(elem)
    assert cg == ArrayContraction(ArrayTensorProduct(-2, M, N), (1, 2))


def test_arrayexpr_convert_array_element_to_array_expression():
    A = ArraySymbol("A", (k,))
    B = ArraySymbol("B", (k,))

    s = Sum(A[i]*B[i], (i, 0, k-1))
    cg = convert_indexed_to_array(s)
    assert cg == ArrayContraction(ArrayTensorProduct(A, B), (0, 1))

    s = A[i]*B[i]
    cg = convert_indexed_to_array(s)
    assert cg == ArrayDiagonal(ArrayTensorProduct(A, B), (0, 1))

    s = A[i]*B[j]
    cg = convert_indexed_to_array(s, [i, j])
    assert cg == ArrayTensorProduct(A, B)
    cg = convert_indexed_to_array(s, [j, i])
    assert cg == ArrayTensorProduct(B, A)

    s = tanh(A[i]*B[j])
    cg = convert_indexed_to_array(s, [i, j])
    assert cg.dummy_eq(ArrayElementwiseApplyFunc(tanh, ArrayTensorProduct(A, B)))


def test_arrayexpr_convert_indexed_to_array_and_back_to_matrix():

    expr = a.T*b
    elem = expr[0, 0]
    cg = convert_indexed_to_array(elem)
    assert cg == ArrayElement(ArrayContraction(ArrayTensorProduct(a, b), (0, 2)), [0, 0])

    expr = M[i,j] + N[i,j]
    p1, p2 = _convert_indexed_to_array(expr)
    assert convert_array_to_matrix(p1) == M + N

    expr = M[i,j] + N[j,i]
    p1, p2 = _convert_indexed_to_array(expr)
    assert convert_array_to_matrix(p1) == M + N.T

    expr = M[i,j]*N[k,l] + N[i,j]*M[k,l]
    p1, p2 = _convert_indexed_to_array(expr)
    assert convert_array_to_matrix(p1) == ArrayAdd(
        ArrayTensorProduct(M, N),
        ArrayTensorProduct(N, M))

    expr = (M*N*P)[i, j]
    p1, p2 = _convert_indexed_to_array(expr)
    assert convert_array_to_matrix(p1) == M * N * P

    expr = Sum(M[i,j]*(N*P)[j,m], (j, 0, k-1))
    p1, p2 = _convert_indexed_to_array(expr)
    assert convert_array_to_matrix(p1) == M * N * P

    expr = Sum((P[j, m] + P[m, j])*(M[i,j]*N[m,n] + N[i,j]*M[m,n]), (j, 0, k-1), (m, 0, k-1))
    p1, p2 = _convert_indexed_to_array(expr)
    assert convert_array_to_matrix(p1) == M * P * N + M * P.T * N + N * P * M + N * P.T * M


def test_arrayexpr_convert_indexed_to_array_out_of_bounds():

    expr = Sum(M[i, i], (i, 0, 4))
    raises(ValueError, lambda: convert_indexed_to_array(expr))
    expr = Sum(M[i, i], (i, 0, k))
    raises(ValueError, lambda: convert_indexed_to_array(expr))
    expr = Sum(M[i, i], (i, 1, k-1))
    raises(ValueError, lambda: convert_indexed_to_array(expr))

    expr = Sum(M[i, j]*N[j,m], (j, 0, 4))
    raises(ValueError, lambda: convert_indexed_to_array(expr))
    expr = Sum(M[i, j]*N[j,m], (j, 0, k))
    raises(ValueError, lambda: convert_indexed_to_array(expr))
    expr = Sum(M[i, j]*N[j,m], (j, 1, k-1))
    raises(ValueError, lambda: convert_indexed_to_array(expr))


def test_arrayexpr_convert_indexed_to_array_broadcast():
    A = ArraySymbol("A", (3, 3))
    B = ArraySymbol("B", (3, 3))

    expr = A[i, j] + B[k, l]
    O2 = OneArray(3, 3)
    expected = ArrayAdd(ArrayTensorProduct(A, O2), ArrayTensorProduct(O2, B))
    assert convert_indexed_to_array(expr) == expected
    assert convert_indexed_to_array(expr, [i, j, k, l]) == expected
    assert convert_indexed_to_array(expr, [l, k, i, j]) == ArrayAdd(PermuteDims(ArrayTensorProduct(O2, A), [1, 0, 2, 3]), PermuteDims(ArrayTensorProduct(B, O2), [1, 0, 2, 3]))

    expr = A[i, j] + B[j, k]
    O1 = OneArray(3)
    assert convert_indexed_to_array(expr, [i, j, k]) == ArrayAdd(ArrayTensorProduct(A, O1), ArrayTensorProduct(O1, B))

    C = ArraySymbol("C", (d0, d1))
    D = ArraySymbol("D", (d3, d1))

    expr = C[i, j] + D[k, j]
    assert convert_indexed_to_array(expr, [i, j, k]) == ArrayAdd(ArrayTensorProduct(C, OneArray(d3)), PermuteDims(ArrayTensorProduct(OneArray(d0), D), [0, 2, 1]))

    X = ArraySymbol("X", (5, 3))

    expr = X[i, n] - X[j, n]
    assert convert_indexed_to_array(expr, [i, j, n]) == ArrayAdd(ArrayTensorProduct(-1, OneArray(5), X), PermuteDims(ArrayTensorProduct(X, OneArray(5)), [0, 2, 1]))

    raises(ValueError, lambda: convert_indexed_to_array(C[i, j] + D[i, j]))


def test_convert_indexed_to_array_invalid_first_indices():
    # first_indices entries that are not actually free indices of the
    # expression (a summation index, or a symbol not present at all) must be
    # ignored, giving the same result as if they were omitted.  Previously they
    # were filtered by removing from a list while iterating over it, which
    # skipped elements and let invalid indices survive, eventually crashing
    # _af_invert with an IndexError.
    Ash = IndexedBase("Ash", shape=(4, 4))
    Bsh = IndexedBase("Bsh", shape=(4, 4))
    expr = Sum(Ash[i, j]*Bsh[j, k], (j, 0, 3))
    expected = convert_indexed_to_array(expr, first_indices=[k])
    # l is not in the expression, j is the contracted index -- both dropped:
    assert convert_indexed_to_array(expr, first_indices=[l, j, k]) == expected
    assert convert_indexed_to_array(expr, first_indices=[k, l]) == expected


def test_convert_indexed_to_array_sum_over_subset_of_diagonals():
    # Summing over only some of the diagonal groups of the summand used to
    # mix inner and outer coordinates of the ArrayDiagonal, producing
    # out-of-range or overlapping contraction indices (IndexError/ValueError).
    expr = Sum(M[i, j]**2, (j, 0, k-1))
    assert convert_indexed_to_array(expr) == ArrayDiagonal(
        ArrayContraction(ArrayTensorProduct(M, M), (1, 3)), (0, 1))

    expr = Sum(M[i, i]*N[j, j]*M[i, j], (i, 0, k-1))
    assert convert_indexed_to_array(expr) == ArrayDiagonal(
        ArrayContraction(ArrayTensorProduct(M, M, N), (0, 1, 2)), (0, 1, 2))

    # Numeric verification on explicit 2x2 matrices:
    M2 = MatrixSymbol("M2", 2, 2)
    N2 = MatrixSymbol("N2", 2, 2)

    res = convert_indexed_to_array(Sum(M2[i, j]**2, (j, 0, 1))).as_explicit()
    assert res.shape == (2,)
    for i_ in range(2):
        assert expand(res[i_] - sum(M2[i_, j_]**2 for j_ in range(2))) == 0

    res = convert_indexed_to_array(Sum(M2[i, i]*N2[j, j]*M2[i, j], (i, 0, 1))).as_explicit()
    assert res.shape == (2,)
    for j_ in range(2):
        assert expand(res[j_] - sum(M2[i_, i_]*N2[j_, j_]*M2[i_, j_] for i_ in range(2))) == 0


def test_convert_indexed_to_array_pow():
    d = Dummy("d")

    # Negative integer exponent used to silently return 1:
    expr = M[i, j]**(-1)
    assert convert_indexed_to_array(expr).dummy_eq(
        ArrayElementwiseApplyFunc(Lambda(d, 1/d), M))

    # Symbolic exponent used to fall through and drop the indices:
    expr = M[i, j]**n
    assert convert_indexed_to_array(expr).dummy_eq(
        ArrayElementwiseApplyFunc(Lambda(d, d**n), M))

    # Exponent one is the base itself:
    assert _convert_indexed_to_array(Pow(M[i, j], 1, evaluate=False)) == (M, (i, j))

    # Exponent depending on the array indices is not supported:
    raises(NotImplementedError, lambda: convert_indexed_to_array(M[i, j]**i))

    # Numeric verification on explicit 2x2 matrices:
    M2 = MatrixSymbol("M2", 2, 2)
    res = convert_indexed_to_array(M2[i, j]**(-1)).as_explicit()
    assert res[1, 0] == 1/M2[1, 0]
    res = convert_indexed_to_array(M2[i, j]**n).as_explicit()
    assert res[0, 1] == M2[0, 1]**n


def test_convert_indexed_to_array_multiarg_function():
    # Functions of more than one argument used to silently drop all the
    # arguments after the first one:
    f = Function("f")
    raises(NotImplementedError, lambda: convert_indexed_to_array(f(M[i, j], N[i, j])))
