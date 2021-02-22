from sympy import (
    symbols, IndexedBase, Identity, cos, Inverse, tensorcontraction,
    permutedims, tensorproduct, ZeroMatrix, HadamardProduct, HadamardPower, MatPow, sin)
from sympy.tensor.array.expressions.conv_matrix_to_array import convert_matrix_to_array
from sympy.tensor.array.expressions.conv_array_to_matrix import _support_function_tp1_recognize, \
    _array_diag2contr_diagmatrix, convert_array_to_matrix, _remove_trivial_dims, _array2matrix
from sympy.tensor.array.expressions.conv_indexed_to_array import convert_indexed_to_array, _convert_indexed_to_array
from sympy import MatrixSymbol, Sum
from sympy.combinatorics import Permutation
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy.matrices.expressions.diagonal import DiagMatrix
from sympy.matrices import Trace, MatMul, Transpose
from sympy.tensor.array.expressions.array_expressions import ZeroArray, OneArray, ArraySymbol, \
    ArrayTensorProduct, ArrayAdd, PermuteDims, ArrayDiagonal, \
    ArrayElementwiseApplyFunc, ArrayContraction, nest_permutation
from sympy.testing.pytest import raises
import random


A, B = symbols("A B", cls=IndexedBase)
i, j, k, l, m, n = symbols("i j k l m n")

I = Identity(k)

M = MatrixSymbol("M", k, k)
N = MatrixSymbol("N", k, k)
P = MatrixSymbol("P", k, k)
Q = MatrixSymbol("Q", k, k)

a = MatrixSymbol("a", k, 1)
b = MatrixSymbol("b", k, 1)
c = MatrixSymbol("c", k, 1)
d = MatrixSymbol("d", k, 1)


def test_codegen_array_contraction_construction():
    cg = ArrayContraction(A)
    assert cg == A

    s = Sum(A[i]*B[i], (i, 0, 3))
    cg = convert_indexed_to_array(s)
    assert cg == ArrayContraction(ArrayTensorProduct(A, B), (0, 1))

    cg = ArrayContraction(ArrayTensorProduct(A, B), (1, 0))
    assert cg == ArrayContraction(ArrayTensorProduct(A, B), (0, 1))

    expr = M*N
    result = ArrayContraction(ArrayTensorProduct(M, N), (1, 2))
    assert convert_matrix_to_array(expr) == result
    elem = expr[i, j]
    assert convert_indexed_to_array(elem) == result

    expr = M*N*M
    result = ArrayContraction(ArrayTensorProduct(M, N, M), (1, 2), (3, 4))
    assert convert_matrix_to_array(expr) == result
    elem = expr[i, j]
    result = ArrayContraction(ArrayTensorProduct(M, M, N), (1, 4), (2, 5))
    cg = convert_indexed_to_array(elem)
    cg = cg.sort_args_by_name()
    assert cg == result


def test_codegen_array_contraction_indices_types():

    cg = ArrayContraction(ArrayTensorProduct(M, N), (0, 1))
    indtup = cg._get_contraction_tuples()
    assert indtup == [[(0, 0), (0, 1)]]
    assert cg._contraction_tuples_to_contraction_indices(cg.expr, indtup) == [(0, 1)]

    cg = ArrayContraction(ArrayTensorProduct(M, N), (1, 2))
    indtup = cg._get_contraction_tuples()
    assert indtup == [[(0, 1), (1, 0)]]
    assert cg._contraction_tuples_to_contraction_indices(cg.expr, indtup) == [(1, 2)]

    cg = ArrayContraction(ArrayTensorProduct(M, M, N), (1, 4), (2, 5))
    indtup = cg._get_contraction_tuples()
    assert indtup == [[(0, 0), (1, 1)], [(0, 1), (2, 0)]]
    assert cg._contraction_tuples_to_contraction_indices(cg.expr, indtup) == [(0, 3), (1, 4)]


def test_codegen_array_recognize_matrix_mul_lines():

    cg = ArrayContraction(ArrayTensorProduct(M), (0, 1))
    assert convert_array_to_matrix(cg) == Trace(M)

    cg = ArrayContraction(ArrayTensorProduct(M, N), (0, 1), (2, 3))
    assert convert_array_to_matrix(cg) == Trace(M) * Trace(N)

    cg = ArrayContraction(ArrayTensorProduct(M, N), (0, 3), (1, 2))
    assert convert_array_to_matrix(cg) == Trace(M * N)

    cg = ArrayContraction(ArrayTensorProduct(M, N), (0, 2), (1, 3))
    assert convert_array_to_matrix(cg) == Trace(M * N.T)

    cg = convert_indexed_to_array((M * N * P)[i, j])
    assert convert_array_to_matrix(cg) == M * N * P
    cg = convert_matrix_to_array(M * N * P)
    assert convert_array_to_matrix(cg) == M * N * P

    cg = convert_indexed_to_array((M * N.T * P)[i, j])
    assert convert_array_to_matrix(cg) == M * N.T * P
    cg = convert_matrix_to_array(M * N.T * P)
    assert convert_array_to_matrix(cg) == M * N.T * P

    cg = ArrayContraction(ArrayTensorProduct(M,N,P,Q), (1, 2), (5, 6))
    assert convert_array_to_matrix(cg) == ArrayTensorProduct(M * N, P * Q)

    expr = -2*M*N
    elem = expr[i, j]
    cg = convert_indexed_to_array(elem)
    assert convert_array_to_matrix(cg) == -2 * M * N

    a = MatrixSymbol("a", k, 1)
    b = MatrixSymbol("b", k, 1)
    c = MatrixSymbol("c", k, 1)
    cg = PermuteDims(
        ArrayContraction(
            ArrayTensorProduct(
                a,
                ArrayAdd(
                    ArrayTensorProduct(b, c),
                    ArrayTensorProduct(c, b),
                )
            ), (2, 4)), [0, 1, 3, 2])
    assert convert_array_to_matrix(cg) == a * (b.T * c + c.T * b)

    za = ZeroArray(m, n)
    assert convert_array_to_matrix(za) == ZeroMatrix(m, n)

    cg = ArrayTensorProduct(3, M)
    assert convert_array_to_matrix(cg) == 3 * M

    # TODO: not yet supported:

    # cg = ArrayDiagonal(ArrayTensorProduct(M, N, P), (0, 2, 4), (1, 3, 5))
    # assert recognize_matrix_expression(cg) == HadamardProduct(M, N, P)

    # cg = ArrayDiagonal(ArrayTensorProduct(M, N, P), (0, 3, 4), (1, 2, 5))
    # assert recognize_matrix_expression(cg) == HadamardProduct(M, N.T, P)

    x = MatrixSymbol("x", k, 1)
    cg = PermuteDims(
        ArrayContraction(ArrayTensorProduct(OneArray(1), x, OneArray(1), DiagMatrix(Identity(1))),
                                (0, 5)), Permutation(1, 2, 3))
    assert convert_array_to_matrix(cg) == x


def test_codegen_array_flatten():

    # Flatten nested ArrayTensorProduct objects:
    expr1 = ArrayTensorProduct(M, N)
    expr2 = ArrayTensorProduct(P, Q)
    expr = ArrayTensorProduct(expr1, expr2)
    assert expr == ArrayTensorProduct(M, N, P, Q)
    assert expr.args == (M, N, P, Q)

    # Flatten mixed ArrayTensorProduct and ArrayContraction objects:
    cg1 = ArrayContraction(expr1, (1, 2))
    cg2 = ArrayContraction(expr2, (0, 3))

    expr = ArrayTensorProduct(cg1, cg2)
    assert expr == ArrayContraction(ArrayTensorProduct(M, N, P, Q), (1, 2), (4, 7))

    expr = ArrayTensorProduct(M, cg1)
    assert expr == ArrayContraction(ArrayTensorProduct(M, M, N), (3, 4))

    # Flatten nested ArrayContraction objects:
    cgnested = ArrayContraction(cg1, (0, 1))
    assert cgnested == ArrayContraction(ArrayTensorProduct(M, N), (0, 3), (1, 2))

    cgnested = ArrayContraction(ArrayTensorProduct(cg1, cg2), (0, 3))
    assert cgnested == ArrayContraction(ArrayTensorProduct(M, N, P, Q), (0, 6), (1, 2), (4, 7))

    cg3 = ArrayContraction(ArrayTensorProduct(M, N, P, Q), (1, 3), (2, 4))
    cgnested = ArrayContraction(cg3, (0, 1))
    assert cgnested == ArrayContraction(ArrayTensorProduct(M, N, P, Q), (0, 5), (1, 3), (2, 4))

    cgnested = ArrayContraction(cg3, (0, 3), (1, 2))
    assert cgnested == ArrayContraction(ArrayTensorProduct(M, N, P, Q), (0, 7), (1, 3), (2, 4), (5, 6))

    cg4 = ArrayContraction(ArrayTensorProduct(M, N, P, Q), (1, 5), (3, 7))
    cgnested = ArrayContraction(cg4, (0, 1))
    assert cgnested == ArrayContraction(ArrayTensorProduct(M, N, P, Q), (0, 2), (1, 5), (3, 7))

    cgnested = ArrayContraction(cg4, (0, 1), (2, 3))
    assert cgnested == ArrayContraction(ArrayTensorProduct(M, N, P, Q), (0, 2), (1, 5), (3, 7), (4, 6))

    cg = ArrayDiagonal(cg4)
    assert cg == cg4
    assert isinstance(cg, type(cg4))

    # Flatten nested ArrayDiagonal objects:
    cg1 = ArrayDiagonal(expr1, (1, 2))
    cg2 = ArrayDiagonal(expr2, (0, 3))
    cg3 = ArrayDiagonal(ArrayTensorProduct(M, N, P, Q), (1, 3), (2, 4))
    cg4 = ArrayDiagonal(ArrayTensorProduct(M, N, P, Q), (1, 5), (3, 7))

    cgnested = ArrayDiagonal(cg1, (0, 1))
    assert cgnested == ArrayDiagonal(ArrayTensorProduct(M, N), (1, 2), (0, 3))

    cgnested = ArrayDiagonal(cg3, (1, 2))
    assert cgnested == ArrayDiagonal(ArrayTensorProduct(M, N, P, Q), (1, 3), (2, 4), (5, 6))

    cgnested = ArrayDiagonal(cg4, (1, 2))
    assert cgnested == ArrayDiagonal(ArrayTensorProduct(M, N, P, Q), (1, 5), (3, 7), (2, 4))

    cg = ArrayAdd(M, N)
    cg2 = ArrayAdd(cg, P)
    assert isinstance(cg2, ArrayAdd)
    assert cg2.args == (M, N, P)
    assert cg2.shape == (k, k)


def test_codegen_array_parse():
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
    assert _convert_indexed_to_array(expr) == (ArrayContraction(ArrayTensorProduct(M, N, P), (1, 2), (3, 4)), (i, j))
    expr = expr.function  # Disregard summation in previous expression
    ret1, ret2 =  _convert_indexed_to_array(expr)
    assert ret1 == ArrayDiagonal(ArrayTensorProduct(M, N, P), (1, 2), (3, 4))
    assert str(ret2) == "(i, j, _i_1, _i_2)"
    expr = KroneckerDelta(i, j)*M[i, k]
    assert _convert_indexed_to_array(expr) == (M, ({i, j}, k))
    expr = KroneckerDelta(i, j)*KroneckerDelta(j, k)*M[i, l]
    assert _convert_indexed_to_array(expr) == (M, ({i, j, k}, l))
    expr = KroneckerDelta(j, k)*(M[i, j]*N[k, l] + N[i, j]*M[k, l])
    assert _convert_indexed_to_array(expr) == (ArrayDiagonal(ArrayAdd(
            ArrayTensorProduct(M, N),
            PermuteDims(ArrayTensorProduct(M, N), Permutation(0, 2)(1, 3))
        ), (1, 2)), (i, l, frozenset({j, k})))
    expr = KroneckerDelta(j, m)*KroneckerDelta(m, k)*(M[i, j]*N[k, l] + N[i, j]*M[k, l])
    assert _convert_indexed_to_array(expr) == (ArrayDiagonal(ArrayAdd(
            ArrayTensorProduct(M, N),
            PermuteDims(ArrayTensorProduct(M, N), Permutation(0, 2)(1, 3))
        ), (1, 2)), (i, l, frozenset({j, m, k})))
    expr = KroneckerDelta(i, j)*KroneckerDelta(j, k)*KroneckerDelta(k,m)*M[i, 0]*KroneckerDelta(m, n)
    assert _convert_indexed_to_array(expr) == (M, ({i, j, k, m, n}, 0))
    expr = M[i, i]
    assert _convert_indexed_to_array(expr) == (ArrayDiagonal(M, (0, 1)), (i,))


def test_codegen_array_diagonal():
    cg = ArrayDiagonal(M, (1, 0))
    assert cg == ArrayDiagonal(M, (0, 1))

    cg = ArrayDiagonal(ArrayTensorProduct(M, N, P), (4, 1), (2, 0))
    assert cg == ArrayDiagonal(ArrayTensorProduct(M, N, P), (1, 4), (0, 2))


def test_codegen_recognize_matrix_expression():

    expr = ArrayAdd(M, PermuteDims(M, [1, 0]))
    assert convert_array_to_matrix(expr) == M + Transpose(M)

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


def test_codegen_array_shape():
    expr = ArrayTensorProduct(M, N, P, Q)
    assert expr.shape == (k, k, k, k, k, k, k, k)
    Z = MatrixSymbol("Z", m, n)
    expr = ArrayTensorProduct(M, Z)
    assert expr.shape == (k, k, m, n)
    expr2 = ArrayContraction(expr, (0, 1))
    assert expr2.shape == (m, n)
    expr2 = ArrayDiagonal(expr, (0, 1))
    assert expr2.shape == (m, n, k)
    exprp = PermuteDims(expr, [2, 1, 3, 0])
    assert exprp.shape == (m, k, n, k)
    expr3 = ArrayTensorProduct(N, Z)
    expr2 = ArrayAdd(expr, expr3)
    assert expr2.shape == (k, k, m, n)

    # Contraction along axes with discordant dimensions:
    raises(ValueError, lambda: ArrayContraction(expr, (1, 2)))
    # Also diagonal needs the same dimensions:
    raises(ValueError, lambda: ArrayDiagonal(expr, (1, 2)))
    # Diagonal requires at least to axes to compute the diagonal:
    raises(ValueError, lambda: ArrayDiagonal(expr, (1,)))


def test_codegen_array_parse_out_of_bounds():

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


def test_codegen_permutedims_sink():

    cg = PermuteDims(ArrayTensorProduct(M, N), [0, 1, 3, 2], nest_permutation=False)
    sunk = nest_permutation(cg)
    assert sunk == ArrayTensorProduct(M, PermuteDims(N, [1, 0]))
    assert convert_array_to_matrix(sunk) == ArrayTensorProduct(M, N.T)

    cg = PermuteDims(ArrayTensorProduct(M, N), [1, 0, 3, 2], nest_permutation=False)
    sunk = nest_permutation(cg)
    assert sunk == ArrayTensorProduct(PermuteDims(M, [1, 0]), PermuteDims(N, [1, 0]))
    assert convert_array_to_matrix(sunk) == ArrayTensorProduct(M.T, N.T)

    cg = PermuteDims(ArrayTensorProduct(M, N), [3, 2, 1, 0], nest_permutation=False)
    sunk = nest_permutation(cg)
    assert sunk == ArrayTensorProduct(PermuteDims(N, [1, 0]), PermuteDims(M, [1, 0]))
    assert convert_array_to_matrix(sunk) == ArrayTensorProduct(N.T, M.T)

    cg = PermuteDims(ArrayContraction(ArrayTensorProduct(M, N), (1, 2)), [1, 0], nest_permutation=False)
    sunk = nest_permutation(cg)
    assert sunk == ArrayContraction(PermuteDims(ArrayTensorProduct(M, N), [[0, 3]]), (1, 2))

    cg = PermuteDims(ArrayTensorProduct(M, N), [1, 0, 3, 2], nest_permutation=False)
    sunk = nest_permutation(cg)
    assert sunk == ArrayTensorProduct(PermuteDims(M, [1, 0]), PermuteDims(N, [1, 0]))

    cg = PermuteDims(ArrayContraction(ArrayTensorProduct(M, N, P), (1, 2), (3, 4)), [1, 0], nest_permutation=False)
    sunk = nest_permutation(cg)
    assert sunk == ArrayContraction(PermuteDims(ArrayTensorProduct(M, N, P), [[0, 5]]), (1, 2), (3, 4))


def test_parsing_of_matrix_expressions():

    expr = M*N
    assert convert_matrix_to_array(expr) == ArrayContraction(ArrayTensorProduct(M, N), (1, 2))

    expr = Transpose(M)
    assert convert_matrix_to_array(expr) == PermuteDims(M, [1, 0])

    expr = M*Transpose(N)
    assert convert_matrix_to_array(expr) == ArrayContraction(ArrayTensorProduct(M, PermuteDims(N, [1, 0])), (1, 2))

    expr = 3*M*N
    res = convert_matrix_to_array(expr)
    rexpr = convert_array_to_matrix(res)
    assert expr == rexpr

    expr = 3*M + N*M.T*M + 4*k*N
    res = convert_matrix_to_array(expr)
    rexpr = convert_array_to_matrix(res)
    assert expr == rexpr

    expr = Inverse(M)*N
    rexpr = convert_array_to_matrix(convert_matrix_to_array(expr))
    assert expr == rexpr

    expr = M**2
    rexpr = convert_array_to_matrix(convert_matrix_to_array(expr))
    assert expr == rexpr

    expr = M*(2*N + 3*M)
    res = convert_matrix_to_array(expr)
    rexpr = convert_array_to_matrix(res)
    assert expr == rexpr

    expr = Trace(M)
    result = ArrayContraction(M, (0, 1))
    assert convert_matrix_to_array(expr) == result

    expr = 3*Trace(M)
    result = ArrayContraction(ArrayTensorProduct(3, M), (0, 1))
    assert convert_matrix_to_array(expr) == result

    expr = 3*Trace(Trace(M) * M)
    result = ArrayContraction(ArrayTensorProduct(3, M, M), (0, 1), (2, 3))
    assert convert_matrix_to_array(expr) == result

    expr = 3*Trace(M)**2
    result = ArrayContraction(ArrayTensorProduct(3, M, M), (0, 1), (2, 3))
    assert convert_matrix_to_array(expr) == result

    expr = HadamardProduct(M, N)
    result = ArrayDiagonal(ArrayTensorProduct(M, N), (0, 2), (1, 3))
    assert convert_matrix_to_array(expr) == result

    expr = HadamardPower(M, 2)
    result = ArrayDiagonal(ArrayTensorProduct(M, M), (0, 2), (1, 3))
    assert convert_matrix_to_array(expr) == result

    expr = M**2
    assert isinstance(expr, MatPow)
    assert convert_matrix_to_array(expr) == ArrayContraction(ArrayTensorProduct(M, M), (1, 2))


def test_special_matrices():
    a = MatrixSymbol("a", k, 1)
    b = MatrixSymbol("b", k, 1)

    expr = a.T*b
    elem = expr[0, 0]
    cg = convert_indexed_to_array(elem)
    assert cg == ArrayContraction(ArrayTensorProduct(a, b), (0, 2))
    assert convert_array_to_matrix(cg) == a.T * b


def test_push_indices_up_and_down():

    indices = list(range(12))

    contr_diag_indices = [(0, 6), (2, 8)]
    assert ArrayContraction._push_indices_down(contr_diag_indices, indices) == (1, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15)
    assert ArrayContraction._push_indices_up(contr_diag_indices, indices) == (None, 0, None, 1, 2, 3, None, 4, None, 5, 6, 7)

    assert ArrayDiagonal._push_indices_down(contr_diag_indices, indices, 10) == (1, 3, 4, 5, 7, 9, (0, 6), (2, 8), None, None, None, None)
    assert ArrayDiagonal._push_indices_up(contr_diag_indices, indices, 10) == (6, 0, 7, 1, 2, 3, 6, 4, 7, 5, None, None)

    contr_diag_indices = [(1, 2), (7, 8)]
    assert ArrayContraction._push_indices_down(contr_diag_indices, indices) == (0, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15)
    assert ArrayContraction._push_indices_up(contr_diag_indices, indices) == (0, None, None, 1, 2, 3, 4, None, None, 5, 6, 7)

    assert ArrayDiagonal._push_indices_down(contr_diag_indices, indices, 10) == (0, 3, 4, 5, 6, 9, (1, 2), (7, 8), None, None, None, None)
    assert ArrayDiagonal._push_indices_up(contr_diag_indices, indices, 10) == (0, 6, 6, 1, 2, 3, 4, 7, 7, 5, None, None)


def test_recognize_diagonalized_vectors():

    a = MatrixSymbol("a", k, 1)
    b = MatrixSymbol("b", k, 1)
    A = MatrixSymbol("A", k, k)
    B = MatrixSymbol("B", k, k)
    C = MatrixSymbol("C", k, k)
    X = MatrixSymbol("X", k, k)
    x = MatrixSymbol("x", k, 1)
    I1 = Identity(1)
    I = Identity(k)

    # Check matrix recognition over trivial dimensions:

    cg = ArrayTensorProduct(a, b)
    assert convert_array_to_matrix(cg) == a * b.T

    cg = ArrayTensorProduct(I1, a, b)
    assert convert_array_to_matrix(cg) == a * b.T

    # Recognize trace inside a tensor product:

    cg = ArrayContraction(ArrayTensorProduct(A, B, C), (0, 3), (1, 2))
    assert convert_array_to_matrix(cg) == Trace(A * B) * C

    # Transform diagonal operator to contraction:

    cg = ArrayDiagonal(ArrayTensorProduct(A, a), (1, 2))
    assert _array_diag2contr_diagmatrix(cg) == ArrayContraction(ArrayTensorProduct(A, OneArray(1), DiagMatrix(a)), (1, 3))
    assert convert_array_to_matrix(cg) == A * DiagMatrix(a)

    cg = ArrayDiagonal(ArrayTensorProduct(a, b), (0, 2))
    assert _array_diag2contr_diagmatrix(cg) == PermuteDims(
        ArrayContraction(ArrayTensorProduct(DiagMatrix(a), OneArray(1), b), (0, 3)), [1, 2, 0]
    )
    assert convert_array_to_matrix(cg) == b.T * DiagMatrix(a)

    cg = ArrayDiagonal(ArrayTensorProduct(A, a), (0, 2))
    assert _array_diag2contr_diagmatrix(cg) == ArrayContraction(ArrayTensorProduct(A, OneArray(1), DiagMatrix(a)), (0, 3))
    assert convert_array_to_matrix(cg) == A.T * DiagMatrix(a)

    cg = ArrayDiagonal(ArrayTensorProduct(I, x, I1), (0, 2), (3, 5))
    assert _array_diag2contr_diagmatrix(cg) == ArrayContraction(ArrayTensorProduct(I, OneArray(1), I1, DiagMatrix(x)), (0, 5))
    assert convert_array_to_matrix(cg) == DiagMatrix(x)

    cg = ArrayDiagonal(ArrayTensorProduct(I, x, A, B), (1, 2), (5, 6))
    assert _array_diag2contr_diagmatrix(cg) == ArrayDiagonal(ArrayContraction(ArrayTensorProduct(I, OneArray(1), A, B, DiagMatrix(x)), (1, 7)), (5, 6))
    # TODO: not yet working
    #  assert recognize_matrix_expression(cg)

    cg = ArrayDiagonal(ArrayTensorProduct(x, I1), (1, 2))
    assert isinstance(cg, ArrayDiagonal)
    assert cg.diagonal_indices == ((1, 2),)
    assert convert_array_to_matrix(cg) == x

    cg = ArrayDiagonal(ArrayTensorProduct(x, I), (0, 2))
    assert _array_diag2contr_diagmatrix(cg) == ArrayContraction(ArrayTensorProduct(OneArray(1), I, DiagMatrix(x)), (1, 3))
    assert convert_array_to_matrix(cg).doit() == DiagMatrix(x)

    raises(ValueError, lambda: ArrayDiagonal(x, (1,)))

    # Ignore identity matrices with contractions:

    cg = ArrayContraction(ArrayTensorProduct(I, A, I, I), (0, 2), (1, 3), (5, 7))
    assert cg.split_multiple_contractions() == cg
    assert convert_array_to_matrix(cg) == Trace(A) * I

    cg = ArrayContraction(ArrayTensorProduct(Trace(A) * I, I, I), (1, 5), (3, 4))
    assert cg.split_multiple_contractions() == cg
    assert convert_array_to_matrix(cg).doit() == Trace(A) * I

    # Add DiagMatrix when required:

    cg = ArrayContraction(ArrayTensorProduct(A, a), (1, 2))
    assert cg.split_multiple_contractions() == cg
    assert convert_array_to_matrix(cg) == A * a

    cg = ArrayContraction(ArrayTensorProduct(A, a, B), (1, 2, 4))
    assert cg.split_multiple_contractions() == ArrayContraction(ArrayTensorProduct(A, DiagMatrix(a), B), (1, 2), (3, 4))
    assert convert_array_to_matrix(cg) == A * DiagMatrix(a) * B

    cg = ArrayContraction(ArrayTensorProduct(A, a, B), (0, 2, 4))
    assert cg.split_multiple_contractions() == ArrayContraction(ArrayTensorProduct(A, DiagMatrix(a), B), (0, 2), (3, 4))
    assert convert_array_to_matrix(cg) == A.T * DiagMatrix(a) * B

    cg = ArrayContraction(ArrayTensorProduct(A, a, b, a.T, B), (0, 2, 4, 7, 9))
    assert cg.split_multiple_contractions() == ArrayContraction(ArrayTensorProduct(A, DiagMatrix(a), DiagMatrix(b),
                                                                                                 DiagMatrix(a), B),
                                                                       (0, 2), (3, 4), (5, 7), (6, 9))
    assert convert_array_to_matrix(cg).doit() == A.T * DiagMatrix(a) * DiagMatrix(b) * DiagMatrix(a) * B.T

    cg = ArrayContraction(ArrayTensorProduct(I1, I1, I1), (1, 2, 4))
    assert cg.split_multiple_contractions() == ArrayContraction(ArrayTensorProduct(I1, I1, I1), (1, 2), (3, 4))
    assert convert_array_to_matrix(cg).doit() == Identity(1)

    cg = ArrayContraction(ArrayTensorProduct(I, I, I, I, A), (1, 2, 8), (5, 6, 9))
    assert convert_array_to_matrix(cg.split_multiple_contractions()).doit() == A

    cg = ArrayContraction(ArrayTensorProduct(A, a, C, a, B), (1, 2, 4), (5, 6, 8))
    expected = ArrayContraction(ArrayTensorProduct(DiagMatrix(a), DiagMatrix(a), C, A, B), (0, 4), (1, 7), (2, 5), (3, 8))
    assert cg.split_multiple_contractions() == expected
    assert convert_array_to_matrix(cg) == A * DiagMatrix(a) * C * DiagMatrix(a) * B

    cg = ArrayContraction(ArrayTensorProduct(a, I1, b, I1, (a.T*b).applyfunc(cos)), (1, 2, 8), (5, 6, 9))
    assert cg.split_multiple_contractions().dummy_eq(ArrayContraction(ArrayTensorProduct((a.T * b).applyfunc(cos), I1, I1, a, b), (0, 2), (1, 4), (3, 7), (5, 9)))
    assert convert_array_to_matrix(cg).doit().dummy_eq(MatMul(a, (a.T * b).applyfunc(cos), b.T))

    cg = ArrayContraction(ArrayTensorProduct(A.T, a, b, b.T, (A*X*b).applyfunc(cos)), (1, 2, 8), (5, 6, 9))
    assert cg.split_multiple_contractions().dummy_eq(ArrayContraction(ArrayTensorProduct(DiagMatrix(a), (A*X*b).applyfunc(cos), A.T, b, b.T), (0, 2), (1, 5), (3, 7, 8)))
    # assert recognize_matrix_expression(cg)

    # Check no overlap of lines:

    cg = ArrayContraction(ArrayTensorProduct(A, a, C, a, B), (1, 2, 4), (5, 6, 8), (3, 7))
    assert cg.split_multiple_contractions() == cg

    cg = ArrayContraction(ArrayTensorProduct(a, b, A), (0, 2, 4), (1, 3))
    assert cg.split_multiple_contractions() == cg


def test_nested_permutations():

    cg = PermuteDims(PermuteDims(M, (1, 0)), (1, 0))
    assert cg == M

    times = 3
    plist1 = [list(range(6)) for i in range(times)]
    plist2 = [list(range(6)) for i in range(times)]

    for i in range(times):
        random.shuffle(plist1[i])
        random.shuffle(plist2[i])

    plist1.append([2, 5, 4, 1, 0, 3])
    plist2.append([3, 5, 0, 4, 1, 2])

    plist1.append([2, 5, 4, 0, 3, 1])
    plist2.append([3, 0, 5, 1, 2, 4])

    plist1.append([5, 4, 2, 0, 3, 1])
    plist2.append([4, 5, 0, 2, 3, 1])

    Me = M.subs(k, 3).as_explicit()
    Ne = N.subs(k, 3).as_explicit()
    Pe = P.subs(k, 3).as_explicit()
    cge = tensorproduct(Me, Ne, Pe)

    for permutation_array1, permutation_array2 in zip(plist1, plist2):
        p1 = Permutation(permutation_array1)
        p2 = Permutation(permutation_array2)

        cg = PermuteDims(
            PermuteDims(
                ArrayTensorProduct(M, N, P),
                p1),
            p2
        )
        result = PermuteDims(
            ArrayTensorProduct(M, N, P),
            p2*p1
        )
        assert cg == result

        # Check that `permutedims` behaves the same way with explicit-component arrays:
        result1 = permutedims(permutedims(cge, p1), p2)
        result2 = permutedims(cge, p2*p1)
        assert result1 == result2


def test_contraction_permutation_mix():

    Me = M.subs(k, 3).as_explicit()
    Ne = N.subs(k, 3).as_explicit()

    cg1 = ArrayContraction(PermuteDims(ArrayTensorProduct(M, N), Permutation([0, 2, 1, 3])), (2, 3))
    cg2 = ArrayContraction(ArrayTensorProduct(M, N), (1, 3))
    assert cg1 == cg2
    assert convert_array_to_matrix(cg2) == M * N.T
    cge1 = tensorcontraction(permutedims(tensorproduct(Me, Ne), Permutation([0, 2, 1, 3])), (2, 3))
    cge2 = tensorcontraction(tensorproduct(Me, Ne), (1, 3))
    assert cge1 == cge2

    cg1 = PermuteDims(ArrayTensorProduct(M, N), Permutation([0, 1, 3, 2]))
    cg2 = ArrayTensorProduct(M, PermuteDims(N, Permutation([1, 0])))
    assert cg1 == cg2
    assert convert_array_to_matrix(cg1) == ArrayTensorProduct(M, N.T)
    assert convert_array_to_matrix(cg2) == ArrayTensorProduct(M, N.T)

    cg1 = ArrayContraction(
        PermuteDims(
            ArrayTensorProduct(M, N, P, Q), Permutation([0, 2, 3, 1, 4, 5, 7, 6])),
        (1, 2), (3, 5)
    )
    cg2 = ArrayContraction(
        ArrayTensorProduct(M, N, P, PermuteDims(Q, Permutation([1, 0]))),
        (1, 5), (2, 3)
    )
    assert cg1 == cg2
    assert convert_array_to_matrix(cg1) == ArrayTensorProduct(M * P.T * Trace(N), Q.T)
    assert convert_array_to_matrix(cg2) == ArrayTensorProduct(M * P.T * Trace(N), Q.T)

    cg1 = ArrayContraction(
        PermuteDims(
            ArrayTensorProduct(M, N, P, Q), Permutation([1, 0, 4, 6, 2, 7, 5, 3])),
        (0, 1), (2, 6), (3, 7)
    )
    cg2 = PermuteDims(
        ArrayContraction(
            ArrayTensorProduct(M, P, Q, N),
            (0, 1), (2, 3), (4, 7)),
        [1, 0]
    )
    assert cg1 == cg2

    cg1 = ArrayContraction(
        PermuteDims(
            ArrayTensorProduct(M, N, P, Q), Permutation([1, 0, 4, 6, 7, 2, 5, 3])),
        (0, 1), (2, 6), (3, 7)
    )
    cg2 = PermuteDims(
        ArrayContraction(
            ArrayTensorProduct(PermuteDims(M, [1, 0]), N, P, Q),
            (0, 1), (3, 6), (4, 5)
        ),
        Permutation([1, 0])
    )
    assert cg1 == cg2


def test_permute_tensor_product():
    cg1 = PermuteDims(ArrayTensorProduct(M, N, P, Q), Permutation([2, 3, 1, 0, 5, 4, 6, 7]))
    cg2 = ArrayTensorProduct(N, PermuteDims(M, [1, 0]),
                                    PermuteDims(P, [1, 0]), Q)
    assert cg1 == cg2

    # TODO: reverse operation starting with `PermuteDims` and getting down to `bb`...
    cg1 = PermuteDims(ArrayTensorProduct(M, N, P, Q), Permutation([2, 3, 4, 5, 0, 1, 6, 7]))
    cg2 = ArrayTensorProduct(N, P, M, Q)
    assert cg1 == cg2

    cg1 = PermuteDims(ArrayTensorProduct(M, N, P, Q), Permutation([2, 3, 4, 6, 5, 7, 0, 1]))
    assert cg1.expr == ArrayTensorProduct(N, P, Q, M)
    assert cg1.permutation == Permutation([0, 1, 2, 4, 3, 5, 6, 7])

    cg1 = ArrayContraction(
        PermuteDims(
            ArrayTensorProduct(N, Q, Q, M),
            [2, 1, 5, 4, 0, 3, 6, 7]),
        [1, 2, 6])
    cg2 = PermuteDims(ArrayContraction(ArrayTensorProduct(Q, Q, N, M), (3, 5, 6)), [0, 2, 3, 1, 4])
    assert cg1 == cg2

    cg1 = ArrayContraction(
        ArrayContraction(
            ArrayContraction(
                ArrayContraction(
                    PermuteDims(
                        ArrayTensorProduct(N, Q, Q, M),
                        [2, 1, 5, 4, 0, 3, 6, 7]),
                    [1, 2, 6]),
                [1, 3, 4]),
            [1]),
        [0])
    cg2 = ArrayContraction(ArrayTensorProduct(M, N, Q, Q), (0, 3, 5), (1, 4, 7), (2,), (6,))
    assert cg1 == cg2


def test_normalize_diagonal_permutedims():
    tp = ArrayTensorProduct(M, Q, N, P)
    expr = ArrayDiagonal(
        PermuteDims(tp, [0, 1, 2, 4, 7, 6, 3, 5]), (2, 4, 5), (6, 7),
        (0, 3))
    result = ArrayDiagonal(tp, (2, 6, 7), (3, 5), (0, 4))
    assert expr == result

    tp = ArrayTensorProduct(M, N, P, Q)
    expr = ArrayDiagonal(PermuteDims(tp, [0, 5, 2, 4, 1, 6, 3, 7]), (1, 2, 6), (3, 4))
    result = ArrayDiagonal(ArrayTensorProduct(M, P, N, Q), (3, 4, 5), (1, 2))
    assert expr == result


def test_normalize_diagonal_contraction():
    tp = ArrayTensorProduct(M, N, P, Q)
    expr = ArrayContraction(ArrayDiagonal(tp, (1, 3, 4)), (0, 3))
    result = ArrayDiagonal(ArrayContraction(ArrayTensorProduct(M, N, P, Q), (0, 6)), (0, 2, 3))
    assert expr == result

    expr = ArrayContraction(ArrayDiagonal(tp, (0, 1, 2, 3, 7)), (1, 2, 3))
    result = ArrayContraction(ArrayTensorProduct(M, N, P, Q), (0, 1, 2, 3, 5, 6, 7))
    assert expr == result

    expr = ArrayContraction(ArrayDiagonal(tp, (0, 2, 6, 7)), (1, 2, 3))
    result = ArrayDiagonal(ArrayContraction(tp, (3, 4, 5)), (0, 2, 3, 4))
    assert expr == result

    td = ArrayDiagonal(ArrayTensorProduct(M, N, P, Q), (0, 3))
    expr = ArrayContraction(td, (2, 1), (0, 4, 6, 5, 3))
    result = ArrayContraction(ArrayTensorProduct(M, N, P, Q), (0, 1, 3, 5, 6, 7), (2, 4))
    assert expr == result


def test_array_wrong_permutation_size():
    cg = ArrayTensorProduct(M, N)
    raises(ValueError, lambda: PermuteDims(cg, [1, 0]))
    raises(ValueError, lambda: PermuteDims(cg, [1, 0, 2, 3, 5, 4]))


def test_nested_array_elementwise_add():
    cg = ArrayContraction(ArrayAdd(
        ArrayTensorProduct(M, N),
        ArrayTensorProduct(N, M)
    ), (1, 2))
    result = ArrayAdd(
        ArrayContraction(ArrayTensorProduct(M, N), (1, 2)),
        ArrayContraction(ArrayTensorProduct(N, M), (1, 2))
    )
    assert cg == result

    cg = ArrayDiagonal(ArrayAdd(
        ArrayTensorProduct(M, N),
        ArrayTensorProduct(N, M)
    ), (1, 2))
    result = ArrayAdd(
        ArrayDiagonal(ArrayTensorProduct(M, N), (1, 2)),
        ArrayDiagonal(ArrayTensorProduct(N, M), (1, 2))
    )
    assert cg == result


def test_array_expr_zero_array():
    za1 = ZeroArray(k, l, m, n)
    zm1 = ZeroMatrix(m, n)

    za2 = ZeroArray(k, m, m, n)
    zm2 = ZeroMatrix(m, m)
    zm3 = ZeroMatrix(k, k)

    assert ArrayTensorProduct(M, N, za1) == ZeroArray(k, k, k, k, k, l, m, n)
    assert ArrayTensorProduct(M, N, zm1) == ZeroArray(k, k, k, k, m, n)

    assert ArrayContraction(za1, (3,)) == ZeroArray(k, l, m)
    assert ArrayContraction(zm1, (1,)) == ZeroArray(m)
    assert ArrayContraction(za2, (1, 2)) == ZeroArray(k, n)
    assert ArrayContraction(zm2, (0, 1)) == 0

    assert ArrayDiagonal(za2, (1, 2)) == ZeroArray(k, n, m)
    assert ArrayDiagonal(zm2, (0, 1)) == ZeroArray(m)

    assert PermuteDims(za1, [2, 1, 3, 0]) == ZeroArray(m, l, n, k)
    assert PermuteDims(zm1, [1, 0]) == ZeroArray(n, m)

    assert ArrayAdd(za1) == za1
    assert ArrayAdd(zm1) == ZeroArray(m, n)
    tp1 = ArrayTensorProduct(MatrixSymbol("A", k, l), MatrixSymbol("B", m, n))
    assert ArrayAdd(tp1, za1) == tp1
    tp2 = ArrayTensorProduct(MatrixSymbol("C", k, l), MatrixSymbol("D", m, n))
    assert ArrayAdd(tp1, za1, tp2) == ArrayAdd(tp1, tp2)
    assert ArrayAdd(M, zm3) == M
    assert ArrayAdd(M, N, zm3) == ArrayAdd(M, N)


def test_contraction_tp_additions():
    a = ArrayAdd(
        ArrayTensorProduct(M, N),
        ArrayTensorProduct(N, M)
    )
    tp = ArrayTensorProduct(P, a, Q)
    expr = ArrayContraction(tp, (3, 4))
    expected = ArrayTensorProduct(
        P,
        ArrayAdd(
            ArrayContraction(ArrayTensorProduct(M, N), (1, 2)),
            ArrayContraction(ArrayTensorProduct(N, M), (1, 2)),
        ),
        Q
    )
    assert expr == expected
    assert convert_array_to_matrix(expr) == ArrayTensorProduct(P, M * N + N * M, Q)

    expr = ArrayContraction(tp, (1, 2), (3, 4), (5, 6))
    result = ArrayContraction(
        ArrayTensorProduct(
            P,
            ArrayAdd(
                ArrayContraction(ArrayTensorProduct(M, N), (1, 2)),
                ArrayContraction(ArrayTensorProduct(N, M), (1, 2)),
            ),
            Q
        ), (1, 2), (3, 4))
    assert expr == result
    assert convert_array_to_matrix(expr) == P * (M * N + N * M) * Q


def test_recognize_expression_implicit_mul():

    cg = ArrayTensorProduct(a, b)
    assert convert_array_to_matrix(cg) == a * b.T

    cg = ArrayTensorProduct(a, I, b)
    assert convert_array_to_matrix(cg) == a * b.T

    cg = ArrayContraction(ArrayTensorProduct(I, I), (1, 2))
    assert convert_array_to_matrix(cg) == I

    cg = PermuteDims(ArrayTensorProduct(I, Identity(1)), [0, 2, 1, 3])
    assert convert_array_to_matrix(cg) == I


def test_remove_trivial_dims():

    # Tensor Product:
    assert _remove_trivial_dims(ArrayTensorProduct(a, b)) == (a * b.T, [1, 3])
    assert _remove_trivial_dims(ArrayTensorProduct(a.T, b)) == (a * b.T, [0, 3])
    assert _remove_trivial_dims(ArrayTensorProduct(a, b.T)) == (a * b.T, [1, 2])
    assert _remove_trivial_dims(ArrayTensorProduct(a.T, b.T)) == (a * b.T, [0, 2])

    assert _remove_trivial_dims(ArrayTensorProduct(I, a.T, b.T)) == (a * b.T, [0, 1, 2, 4])
    assert _remove_trivial_dims(ArrayTensorProduct(a.T, I, b.T)) == (a * b.T, [0, 2, 3, 4])

    assert _remove_trivial_dims(ArrayTensorProduct(a, I)) == (a, [2, 3])
    assert _remove_trivial_dims(ArrayTensorProduct(I, a)) == (a, [0, 1])

    assert _remove_trivial_dims(ArrayTensorProduct(a.T, b.T, c, d)) == (
        ArrayTensorProduct(a * b.T, c * d.T), [0, 2, 5, 7])
    assert _remove_trivial_dims(ArrayTensorProduct(a.T, I, b.T, c, d, I)) == (
        ArrayTensorProduct(a * b.T, c * d.T, I), [0, 2, 3, 4, 7, 9])

    # Addition:

    cg = ArrayAdd(ArrayTensorProduct(a, b), ArrayTensorProduct(c, d))
    assert _remove_trivial_dims(cg) == (a * b.T + c * d.T, [1, 3])

    # Permute Dims:

    cg = PermuteDims(ArrayTensorProduct(a, b), Permutation(3)(1, 2))
    assert _remove_trivial_dims(cg) == (a * b.T, [2, 3])

    cg = PermuteDims(ArrayTensorProduct(a, I, b), Permutation(5)(1, 2, 3, 4))
    assert _remove_trivial_dims(cg) == (a * b.T, [1, 2, 4, 5])

    cg = PermuteDims(ArrayTensorProduct(I, b, a), Permutation(5)(1, 2, 4, 5, 3))
    assert _remove_trivial_dims(cg) == (b * a.T, [0, 3, 4, 5])

    # Diagonal:

    cg = ArrayDiagonal(ArrayTensorProduct(M, a), (1, 2))
    assert _remove_trivial_dims(cg) == (cg, [])

    # Contraction:

    cg = ArrayContraction(ArrayTensorProduct(M, a), (1, 2))
    assert _remove_trivial_dims(cg) == (cg, [])


def test_diag2contraction_diagmatrix():
    cg = ArrayDiagonal(ArrayTensorProduct(M, a), (1, 2))
    res = _array_diag2contr_diagmatrix(cg)
    assert res.shape == cg.shape
    assert res == ArrayContraction(ArrayTensorProduct(M, OneArray(1), DiagMatrix(a)), (1, 3))

    raises(ValueError, lambda: ArrayDiagonal(ArrayTensorProduct(a, M), (1, 2)))

    cg = ArrayDiagonal(ArrayTensorProduct(a.T, M), (1, 2))
    res = _array_diag2contr_diagmatrix(cg)
    assert res.shape == cg.shape
    assert res == ArrayContraction(ArrayTensorProduct(OneArray(1), M, DiagMatrix(a.T)), (1, 4))

    cg = ArrayDiagonal(ArrayTensorProduct(a.T, M, N, b.T), (1, 2), (4, 7))
    res = _array_diag2contr_diagmatrix(cg)
    assert res.shape == cg.shape
    assert res == ArrayContraction(
        ArrayTensorProduct(OneArray(1), M, N, OneArray(1), DiagMatrix(a.T), DiagMatrix(b.T)), (1, 7), (3, 9))

    cg = ArrayDiagonal(ArrayTensorProduct(a, M, N, b.T), (0, 2), (4, 7))
    res = _array_diag2contr_diagmatrix(cg)
    assert res.shape == cg.shape
    assert res == ArrayContraction(
        ArrayTensorProduct(OneArray(1), M, N, OneArray(1), DiagMatrix(a), DiagMatrix(b.T)), (1, 6), (3, 9))

    cg = ArrayDiagonal(ArrayTensorProduct(a, M, N, b.T), (0, 4), (3, 7))
    res = _array_diag2contr_diagmatrix(cg)
    assert res.shape == cg.shape
    assert res == ArrayContraction(
        ArrayTensorProduct(OneArray(1), M, N, OneArray(1), DiagMatrix(a), DiagMatrix(b.T)), (3, 6), (2, 9))

    I1 = Identity(1)
    x = MatrixSymbol("x", k, 1)
    A = MatrixSymbol("A", k, k)
    cg = ArrayDiagonal(ArrayTensorProduct(x, A.T, I1), (0, 2))
    assert _array_diag2contr_diagmatrix(cg).shape == cg.shape
    assert _array2matrix(cg).shape == cg.shape


def test_recognize_products_support_function():
    A = MatrixSymbol("A", k, k)
    B = MatrixSymbol("B", k, k)
    C = MatrixSymbol("C", k, k)
    D = MatrixSymbol("D", k, k)

    X = MatrixSymbol("X", k, k)
    Y = MatrixSymbol("Y", k, k)

    assert _support_function_tp1_recognize([], [2 * k]) == 2 * k

    assert _support_function_tp1_recognize([(1, 2)], [A, 2 * k, B, 3]) == 6 * k * A * B

    assert _support_function_tp1_recognize([(0, 3), (1, 2)], [A, B]) == Trace(A * B)

    assert _support_function_tp1_recognize([(1, 2)], [A, B]) == A * B
    assert _support_function_tp1_recognize([(0, 2)], [A, B]) == A.T * B
    assert _support_function_tp1_recognize([(1, 3)], [A, B]) == A * B.T
    assert _support_function_tp1_recognize([(0, 3)], [A, B]) == A.T * B.T

    assert _support_function_tp1_recognize([(1, 2), (5, 6)], [A, B, C, D]) == ArrayTensorProduct(A * B, C * D)
    assert _support_function_tp1_recognize([(1, 4), (3, 6)], [A, B, C, D]) == PermuteDims(
        ArrayTensorProduct(A * C, B * D), [0, 2, 1, 3])

    assert _support_function_tp1_recognize([(0, 3), (1, 4)], [A, B, C]) == B * A * C

    assert _support_function_tp1_recognize([(9, 10), (1, 2), (5, 6), (3, 4), (7, 8)],
                                           [X, Y, A, B, C, D]) == X * Y * A * B * C * D

    assert _support_function_tp1_recognize([(9, 10), (1, 2), (5, 6), (3, 4)],
                                           [X, Y, A, B, C, D]) == ArrayTensorProduct(X * Y * A * B, C * D)

    assert _support_function_tp1_recognize([(1, 7), (3, 8), (4, 11)], [X, Y, A, B, C, D]) == PermuteDims(
        ArrayTensorProduct(X * B.T, Y * C, D * A), [0, 2, 5, 1, 3, 4]
    )

    assert _support_function_tp1_recognize([(0, 1), (3, 6), (5, 8)], [X, A, B, C, D]) == PermuteDims(
        ArrayTensorProduct(Trace(X) * A * C, B * D), [0, 2, 1, 3])

    assert _support_function_tp1_recognize([(1, 2), (3, 4), (5, 6), (7, 8)], [A, A, B, C, D]) == A ** 2 * B * C * D
    assert _support_function_tp1_recognize([(1, 2), (3, 4), (5, 6), (7, 8)], [X, A, B, C, D]) == X * A * B * C * D

    assert _support_function_tp1_recognize([(1, 6), (3, 8), (5, 10)], [X, Y, A, B, C, D]) == PermuteDims(
        ArrayTensorProduct(X * B, Y * C, A * D), [0, 2, 4, 1, 3, 5]
    )

    assert _support_function_tp1_recognize([(1, 4), (3, 6)], [A, B, C, D]) == PermuteDims(
        ArrayTensorProduct(A * C, B * D), [0, 2, 1, 3])

    assert _support_function_tp1_recognize([(0, 4), (1, 7), (2, 5), (3, 8)], [X, A, B, C, D]) == C*X.T*B*A*D

    assert _support_function_tp1_recognize([(0, 4), (1, 7), (2, 5), (3, 8)], [X, A, B, C, D]) == C*X.T*B*A*D


def test_array_expr_applyfunc():

    A = ArraySymbol("A", 3, k, 2)
    aaf = ArrayElementwiseApplyFunc(sin, A)
    assert aaf.shape == (3, k, 2)
