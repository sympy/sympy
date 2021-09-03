from sympy import (
    symbols, Identity, cos, ZeroMatrix, OneMatrix, sqrt, HadamardProduct)
from sympy.tensor.array.expressions.conv_matrix_to_array import convert_matrix_to_array
from sympy.tensor.array.expressions.conv_array_to_matrix import _support_function_tp1_recognize, \
    _array_diag2contr_diagmatrix, convert_array_to_matrix, _remove_trivial_dims, _array2matrix
from sympy import MatrixSymbol
from sympy.combinatorics import Permutation
from sympy.matrices.expressions.diagonal import DiagMatrix
from sympy.matrices import Trace, MatMul, Transpose
from sympy.tensor.array.expressions.array_expressions import ZeroArray, OneArray, \
    ArrayTensorProduct, ArrayAdd, PermuteDims, ArrayDiagonal, \
    ArrayContraction
from sympy.testing.pytest import raises


i, j, k, l, m, n = symbols("i j k l m n")

I = Identity(k)
I1 = Identity(1)

M = MatrixSymbol("M", k, k)
N = MatrixSymbol("N", k, k)
P = MatrixSymbol("P", k, k)
Q = MatrixSymbol("Q", k, k)

A = MatrixSymbol("A", k, k)
B = MatrixSymbol("B", k, k)
C = MatrixSymbol("C", k, k)
D = MatrixSymbol("D", k, k)

X = MatrixSymbol("X", k, k)
Y = MatrixSymbol("Y", k, k)

a = MatrixSymbol("a", k, 1)
b = MatrixSymbol("b", k, 1)
c = MatrixSymbol("c", k, 1)
d = MatrixSymbol("d", k, 1)

x = MatrixSymbol("x", k, 1)


def test_arrayexpr_convert_array_to_matrix():

    cg = ArrayContraction(ArrayTensorProduct(M), (0, 1))
    assert convert_array_to_matrix(cg) == Trace(M)

    cg = ArrayContraction(ArrayTensorProduct(M, N), (0, 1), (2, 3))
    assert convert_array_to_matrix(cg) == Trace(M) * Trace(N)

    cg = ArrayContraction(ArrayTensorProduct(M, N), (0, 3), (1, 2))
    assert convert_array_to_matrix(cg) == Trace(M * N)

    cg = ArrayContraction(ArrayTensorProduct(M, N), (0, 2), (1, 3))
    assert convert_array_to_matrix(cg) == Trace(M * N.T)

    cg = convert_matrix_to_array(M * N * P)
    assert convert_array_to_matrix(cg) == M * N * P

    cg = convert_matrix_to_array(M * N.T * P)
    assert convert_array_to_matrix(cg) == M * N.T * P

    cg = ArrayContraction(ArrayTensorProduct(M,N,P,Q), (1, 2), (5, 6))
    assert convert_array_to_matrix(cg) == ArrayTensorProduct(M * N, P * Q)

    cg = ArrayContraction(ArrayTensorProduct(-2, M, N), (1, 2))
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

    # Partial conversion to matrix multiplication:
    expr = ArrayContraction(ArrayTensorProduct(M, N, P, Q), (0, 2), (1, 4, 6))
    assert convert_array_to_matrix(expr) == ArrayContraction(ArrayTensorProduct(M.T*N, P, Q), (0, 2, 4))

    # TODO: not yet supported:

    # cg = ArrayDiagonal(ArrayTensorProduct(M, N, P), (0, 2, 4), (1, 3, 5))
    #  assert recognize_matrix_expression(cg) == HadamardProduct(M, N, P)

    # cg = ArrayDiagonal(ArrayTensorProduct(M, N, P), (0, 3, 4), (1, 2, 5))
    #  assert recognize_matrix_expression(cg) == HadamardProduct(M, N.T, P)

    x = MatrixSymbol("x", k, 1)
    cg = PermuteDims(
        ArrayContraction(ArrayTensorProduct(OneArray(1), x, OneArray(1), DiagMatrix(Identity(1))),
                                (0, 5)), Permutation(1, 2, 3))
    assert convert_array_to_matrix(cg) == x

    expr = ArrayAdd(M, PermuteDims(M, [1, 0]))
    assert convert_array_to_matrix(expr) == M + Transpose(M)


def test_arrayexpr_convert_array_to_matrix2():
    cg = ArrayContraction(ArrayTensorProduct(M, N), (1, 3))
    assert convert_array_to_matrix(cg) == M * N.T

    cg = PermuteDims(ArrayTensorProduct(M, N), Permutation([0, 1, 3, 2]))
    assert convert_array_to_matrix(cg) == ArrayTensorProduct(M, N.T)

    cg = ArrayTensorProduct(M, PermuteDims(N, Permutation([1, 0])))
    assert convert_array_to_matrix(cg) == ArrayTensorProduct(M, N.T)

    cg = ArrayContraction(
        PermuteDims(
            ArrayTensorProduct(M, N, P, Q), Permutation([0, 2, 3, 1, 4, 5, 7, 6])),
        (1, 2), (3, 5)
    )
    assert convert_array_to_matrix(cg) == ArrayTensorProduct(M * P.T * Trace(N), Q.T)

    cg = ArrayContraction(
        ArrayTensorProduct(M, N, P, PermuteDims(Q, Permutation([1, 0]))),
        (1, 5), (2, 3)
    )
    assert convert_array_to_matrix(cg) == ArrayTensorProduct(M * P.T * Trace(N), Q.T)

    cg = ArrayTensorProduct(M, PermuteDims(N, [1, 0]))
    assert convert_array_to_matrix(cg) == ArrayTensorProduct(M, N.T)

    cg = ArrayTensorProduct(PermuteDims(M, [1, 0]), PermuteDims(N, [1, 0]))
    assert convert_array_to_matrix(cg) == ArrayTensorProduct(M.T, N.T)

    cg = ArrayTensorProduct(PermuteDims(N, [1, 0]), PermuteDims(M, [1, 0]))
    assert convert_array_to_matrix(cg) == ArrayTensorProduct(N.T, M.T)


def test_arrayexpr_convert_array_to_diagonalized_vector():

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
    #  assert convert_array_to_matrix(cg)

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


def test_arrayexpr_convert_array_contraction_tp_additions():
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


def test_arrayexpr_convert_array_to_implicit_matmul():
    # Trivial dimensions are suppressed, so the result can be expressed in matrix form:

    cg = ArrayTensorProduct(a, b)
    assert convert_array_to_matrix(cg) == a * b.T

    cg = ArrayTensorProduct(a, I, b)
    assert convert_array_to_matrix(cg) == a * b.T

    cg = ArrayContraction(ArrayTensorProduct(I, I), (1, 2))
    assert convert_array_to_matrix(cg) == I

    cg = PermuteDims(ArrayTensorProduct(I, Identity(1)), [0, 2, 1, 3])
    assert convert_array_to_matrix(cg) == I


def test_arrayexpr_convert_array_to_matrix_remove_trivial_dims():

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

    # A few more cases to test the removal and shift of nested removed axes
    # with array contractions and array diagonals:
    tp = ArrayTensorProduct(
        OneMatrix(1, 1),
        M,
        x,
        OneMatrix(1, 1),
        Identity(1),
    )

    expr = ArrayContraction(tp, (1, 8))
    rexpr, removed = _remove_trivial_dims(expr)
    assert removed == [0, 5, 6, 7]

    expr = ArrayContraction(tp, (1, 8), (3, 4))
    rexpr, removed = _remove_trivial_dims(expr)
    assert removed == [0, 3, 4, 5]

    expr = ArrayDiagonal(tp, (1, 8))
    rexpr, removed = _remove_trivial_dims(expr)
    assert removed == [0, 5, 6, 7, 8]

    expr = ArrayDiagonal(tp, (1, 8), (3, 4))
    rexpr, removed = _remove_trivial_dims(expr)
    assert removed == [0, 3, 4, 5, 6]

    cg = ArrayDiagonal(ArrayTensorProduct(PermuteDims(ArrayTensorProduct(x, I1), Permutation(1, 2, 3)), (x.T*x).applyfunc(sqrt)), (2, 4), (3, 5))
    rexpr, removed = _remove_trivial_dims(cg)
    assert removed == [1, 2, 3]


def test_arrayexpr_convert_array_to_matrix_diag2contraction_diagmatrix():
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


def test_arrayexpr_convert_array_to_matrix_support_function():

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
        ArrayTensorProduct(X * B.T, Y * C, A.T * D.T), [0, 2, 4, 1, 3, 5]
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


def test_convert_array_to_hadamard_products():

    expr = HadamardProduct(M, N)
    cg = convert_matrix_to_array(expr)
    ret = convert_array_to_matrix(cg)
    assert ret == expr

    expr = HadamardProduct(M, N)*P
    cg = convert_matrix_to_array(expr)
    ret = convert_array_to_matrix(cg)
    assert ret == expr

    expr = Q*HadamardProduct(M, N)*P
    cg = convert_matrix_to_array(expr)
    ret = convert_array_to_matrix(cg)
    assert ret == expr

    expr = Q*HadamardProduct(M, N.T)*P
    cg = convert_matrix_to_array(expr)
    ret = convert_array_to_matrix(cg)
    assert ret == expr

    expr = HadamardProduct(M, N)*HadamardProduct(Q, P)
    cg = convert_matrix_to_array(expr)
    ret = convert_array_to_matrix(cg)
    assert expr == ret

    expr = P.T*HadamardProduct(M, N)*HadamardProduct(Q, P)
    cg = convert_matrix_to_array(expr)
    ret = convert_array_to_matrix(cg)
    assert expr == ret

    # ArrayDiagonal should be converted
    cg = ArrayDiagonal(ArrayTensorProduct(M, N, Q), (1, 3), (0, 2, 4))
    ret = convert_array_to_matrix(cg)
    expected = PermuteDims(ArrayDiagonal(ArrayTensorProduct(HadamardProduct(M.T, N.T), Q), (1, 2)), [1, 0, 2])
    assert expected == ret

    # Special case that should return the same expression:
    cg = ArrayDiagonal(ArrayTensorProduct(HadamardProduct(M, N), Q), (0, 2))
    ret = convert_array_to_matrix(cg)
    assert ret == cg

    # Hadamard products with traces:

    expr = Trace(HadamardProduct(M, N))
    cg = convert_matrix_to_array(expr)
    ret = convert_array_to_matrix(cg)
    assert ret == Trace(HadamardProduct(M.T, N.T))

    expr = Trace(A*HadamardProduct(M, N))
    cg = convert_matrix_to_array(expr)
    ret = convert_array_to_matrix(cg)
    assert ret == Trace(HadamardProduct(M, N)*A)

    expr = Trace(HadamardProduct(A, M)*N)
    cg = convert_matrix_to_array(expr)
    ret = convert_array_to_matrix(cg)
    assert ret == Trace(HadamardProduct(M.T, N)*A)

    # These should not be converted into Hadamard products:

    cg = ArrayDiagonal(ArrayTensorProduct(M, N), (0, 1, 2, 3))
    ret = convert_array_to_matrix(cg)
    assert ret == cg

    cg = ArrayDiagonal(ArrayTensorProduct(A), (0, 1))
    ret = convert_array_to_matrix(cg)
    assert ret == cg
