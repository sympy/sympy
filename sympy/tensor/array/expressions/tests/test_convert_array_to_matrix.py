from sympy import (
    symbols, Identity, cos, ZeroMatrix)
from sympy.tensor.array.expressions.conv_matrix_to_array import convert_matrix_to_array
from sympy.tensor.array.expressions.conv_array_to_matrix import _support_function_tp1_recognize, \
    _array_diag2contr_diagmatrix, convert_array_to_matrix, _remove_trivial_dims, _array2matrix
from sympy import MatrixSymbol
from sympy.combinatorics import Permutation
from sympy.matrices.expressions.diagonal import DiagMatrix
from sympy.matrices import Trace, MatMul, Transpose
from sympy.tensor.array.expressions.array_expressions import ZeroArray, OneArray, \
    CodegenArrayTensorProduct, CodegenArrayElementwiseAdd, CodegenArrayPermuteDims, CodegenArrayDiagonal, \
    CodegenArrayContraction
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

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M), (0, 1))
    assert convert_array_to_matrix(cg) == Trace(M)

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (0, 1), (2, 3))
    assert convert_array_to_matrix(cg) == Trace(M) * Trace(N)

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (0, 3), (1, 2))
    assert convert_array_to_matrix(cg) == Trace(M * N)

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (0, 2), (1, 3))
    assert convert_array_to_matrix(cg) == Trace(M * N.T)

    cg = convert_matrix_to_array(M * N * P)
    assert convert_array_to_matrix(cg) == M * N * P

    cg = convert_matrix_to_array(M * N.T * P)
    assert convert_array_to_matrix(cg) == M * N.T * P

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M,N,P,Q), (1, 2), (5, 6))
    assert convert_array_to_matrix(cg) == CodegenArrayTensorProduct(M * N, P * Q)

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(-2, M, N), (1, 2))
    assert convert_array_to_matrix(cg) == -2 * M * N

    a = MatrixSymbol("a", k, 1)
    b = MatrixSymbol("b", k, 1)
    c = MatrixSymbol("c", k, 1)
    cg = CodegenArrayPermuteDims(
        CodegenArrayContraction(
            CodegenArrayTensorProduct(
                a,
                CodegenArrayElementwiseAdd(
                    CodegenArrayTensorProduct(b, c),
                    CodegenArrayTensorProduct(c, b),
                )
            ), (2, 4)), [0, 1, 3, 2])
    assert convert_array_to_matrix(cg) == a * (b.T * c + c.T * b)

    za = ZeroArray(m, n)
    assert convert_array_to_matrix(za) == ZeroMatrix(m, n)

    cg = CodegenArrayTensorProduct(3, M)
    assert convert_array_to_matrix(cg) == 3 * M

    # TODO: not yet supported:

    # cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N, P), (0, 2, 4), (1, 3, 5))
    #  assert recognize_matrix_expression(cg) == HadamardProduct(M, N, P)

    # cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N, P), (0, 3, 4), (1, 2, 5))
    #  assert recognize_matrix_expression(cg) == HadamardProduct(M, N.T, P)

    x = MatrixSymbol("x", k, 1)
    cg = CodegenArrayPermuteDims(
        CodegenArrayContraction(CodegenArrayTensorProduct(OneArray(1), x, OneArray(1), DiagMatrix(Identity(1))),
                                (0, 5)), Permutation(1, 2, 3))
    assert convert_array_to_matrix(cg) == x

    expr = CodegenArrayElementwiseAdd(M, CodegenArrayPermuteDims(M, [1, 0]))
    assert convert_array_to_matrix(expr) == M + Transpose(M)


def test_arrayexpr_convert_array_to_matrix2():
    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (1, 3))
    assert convert_array_to_matrix(cg) == M * N.T

    cg = CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N), Permutation([0, 1, 3, 2]))
    assert convert_array_to_matrix(cg) == CodegenArrayTensorProduct(M, N.T)

    cg = CodegenArrayTensorProduct(M, CodegenArrayPermuteDims(N, Permutation([1, 0])))
    assert convert_array_to_matrix(cg) == CodegenArrayTensorProduct(M, N.T)

    cg = CodegenArrayContraction(
        CodegenArrayPermuteDims(
            CodegenArrayTensorProduct(M, N, P, Q), Permutation([0, 2, 3, 1, 4, 5, 7, 6])),
        (1, 2), (3, 5)
    )
    assert convert_array_to_matrix(cg) == CodegenArrayTensorProduct(M * P.T * Trace(N), Q.T)

    cg = CodegenArrayContraction(
        CodegenArrayTensorProduct(M, N, P, CodegenArrayPermuteDims(Q, Permutation([1, 0]))),
        (1, 5), (2, 3)
    )
    assert convert_array_to_matrix(cg) == CodegenArrayTensorProduct(M * P.T * Trace(N), Q.T)

    cg = CodegenArrayTensorProduct(M, CodegenArrayPermuteDims(N, [1, 0]))
    assert convert_array_to_matrix(cg) == CodegenArrayTensorProduct(M, N.T)

    cg = CodegenArrayTensorProduct(CodegenArrayPermuteDims(M, [1, 0]), CodegenArrayPermuteDims(N, [1, 0]))
    assert convert_array_to_matrix(cg) == CodegenArrayTensorProduct(M.T, N.T)

    cg = CodegenArrayTensorProduct(CodegenArrayPermuteDims(N, [1, 0]), CodegenArrayPermuteDims(M, [1, 0]))
    assert convert_array_to_matrix(cg) == CodegenArrayTensorProduct(N.T, M.T)


def test_arrayexpr_convert_array_to_diagonalized_vector():

    # Check matrix recognition over trivial dimensions:

    cg = CodegenArrayTensorProduct(a, b)
    assert convert_array_to_matrix(cg) == a * b.T

    cg = CodegenArrayTensorProduct(I1, a, b)
    assert convert_array_to_matrix(cg) == a * b.T

    # Recognize trace inside a tensor product:

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, B, C), (0, 3), (1, 2))
    assert convert_array_to_matrix(cg) == Trace(A * B) * C

    # Transform diagonal operator to contraction:

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(A, a), (1, 2))
    assert _array_diag2contr_diagmatrix(cg) == CodegenArrayContraction(CodegenArrayTensorProduct(A, OneArray(1), DiagMatrix(a)), (1, 3))
    assert convert_array_to_matrix(cg) == A * DiagMatrix(a)

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(a, b), (0, 2))
    assert _array_diag2contr_diagmatrix(cg) == CodegenArrayPermuteDims(
        CodegenArrayContraction(CodegenArrayTensorProduct(DiagMatrix(a), OneArray(1), b), (0, 3)), [1, 2, 0]
    )
    assert convert_array_to_matrix(cg) == b.T * DiagMatrix(a)

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(A, a), (0, 2))
    assert _array_diag2contr_diagmatrix(cg) == CodegenArrayContraction(CodegenArrayTensorProduct(A, OneArray(1), DiagMatrix(a)), (0, 3))
    assert convert_array_to_matrix(cg) == A.T * DiagMatrix(a)

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(I, x, I1), (0, 2), (3, 5))
    assert _array_diag2contr_diagmatrix(cg) == CodegenArrayContraction(CodegenArrayTensorProduct(I, OneArray(1), I1, DiagMatrix(x)), (0, 5))
    assert convert_array_to_matrix(cg) == DiagMatrix(x)

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(I, x, A, B), (1, 2), (5, 6))
    assert _array_diag2contr_diagmatrix(cg) == CodegenArrayDiagonal(CodegenArrayContraction(CodegenArrayTensorProduct(I, OneArray(1), A, B, DiagMatrix(x)), (1, 7)), (5, 6))
    # TODO: not yet working
    #  assert recognize_matrix_expression(cg)

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(x, I1), (1, 2))
    assert isinstance(cg, CodegenArrayDiagonal)
    assert cg.diagonal_indices == ((1, 2),)
    assert convert_array_to_matrix(cg) == x

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(x, I), (0, 2))
    assert _array_diag2contr_diagmatrix(cg) == CodegenArrayContraction(CodegenArrayTensorProduct(OneArray(1), I, DiagMatrix(x)), (1, 3))
    assert convert_array_to_matrix(cg).doit() == DiagMatrix(x)

    raises(ValueError, lambda: CodegenArrayDiagonal(x, (1,)))

    # Ignore identity matrices with contractions:

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(I, A, I, I), (0, 2), (1, 3), (5, 7))
    assert cg.split_multiple_contractions() == cg
    assert convert_array_to_matrix(cg) == Trace(A) * I

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(Trace(A) * I, I, I), (1, 5), (3, 4))
    assert cg.split_multiple_contractions() == cg
    assert convert_array_to_matrix(cg).doit() == Trace(A) * I

    # Add DiagMatrix when required:

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, a), (1, 2))
    assert cg.split_multiple_contractions() == cg
    assert convert_array_to_matrix(cg) == A * a

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, a, B), (1, 2, 4))
    assert cg.split_multiple_contractions() == CodegenArrayContraction(CodegenArrayTensorProduct(A, DiagMatrix(a), B), (1, 2), (3, 4))
    assert convert_array_to_matrix(cg) == A * DiagMatrix(a) * B

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, a, B), (0, 2, 4))
    assert cg.split_multiple_contractions() == CodegenArrayContraction(CodegenArrayTensorProduct(A, DiagMatrix(a), B), (0, 2), (3, 4))
    assert convert_array_to_matrix(cg) == A.T * DiagMatrix(a) * B

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, a, b, a.T, B), (0, 2, 4, 7, 9))
    assert cg.split_multiple_contractions() == CodegenArrayContraction(CodegenArrayTensorProduct(A, DiagMatrix(a), DiagMatrix(b),
                                                                                                 DiagMatrix(a), B),
                                                                       (0, 2), (3, 4), (5, 7), (6, 9))
    assert convert_array_to_matrix(cg).doit() == A.T * DiagMatrix(a) * DiagMatrix(b) * DiagMatrix(a) * B.T

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(I1, I1, I1), (1, 2, 4))
    assert cg.split_multiple_contractions() == CodegenArrayContraction(CodegenArrayTensorProduct(I1, I1, I1), (1, 2), (3, 4))
    assert convert_array_to_matrix(cg).doit() == Identity(1)

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(I, I, I, I, A), (1, 2, 8), (5, 6, 9))
    assert convert_array_to_matrix(cg.split_multiple_contractions()).doit() == A

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, a, C, a, B), (1, 2, 4), (5, 6, 8))
    expected = CodegenArrayContraction(CodegenArrayTensorProduct(DiagMatrix(a), DiagMatrix(a), C, A, B), (0, 4), (1, 7), (2, 5), (3, 8))
    assert cg.split_multiple_contractions() == expected
    assert convert_array_to_matrix(cg) == A * DiagMatrix(a) * C * DiagMatrix(a) * B

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(a, I1, b, I1, (a.T*b).applyfunc(cos)), (1, 2, 8), (5, 6, 9))
    assert cg.split_multiple_contractions().dummy_eq(CodegenArrayContraction(CodegenArrayTensorProduct((a.T * b).applyfunc(cos), I1, I1, a, b), (0, 2), (1, 4), (3, 7), (5, 9)))
    assert convert_array_to_matrix(cg).doit().dummy_eq(MatMul(a, (a.T * b).applyfunc(cos), b.T))


def test_arrayexpr_convert_array_contraction_tp_additions():
    a = CodegenArrayElementwiseAdd(
        CodegenArrayTensorProduct(M, N),
        CodegenArrayTensorProduct(N, M)
    )
    tp = CodegenArrayTensorProduct(P, a, Q)
    expr = CodegenArrayContraction(tp, (3, 4))
    expected = CodegenArrayTensorProduct(
        P,
        CodegenArrayElementwiseAdd(
            CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (1, 2)),
            CodegenArrayContraction(CodegenArrayTensorProduct(N, M), (1, 2)),
        ),
        Q
    )
    assert expr == expected
    assert convert_array_to_matrix(expr) == CodegenArrayTensorProduct(P, M * N + N * M, Q)

    expr = CodegenArrayContraction(tp, (1, 2), (3, 4), (5, 6))
    result = CodegenArrayContraction(
        CodegenArrayTensorProduct(
            P,
            CodegenArrayElementwiseAdd(
                CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (1, 2)),
                CodegenArrayContraction(CodegenArrayTensorProduct(N, M), (1, 2)),
            ),
            Q
        ), (1, 2), (3, 4))
    assert expr == result
    assert convert_array_to_matrix(expr) == P * (M * N + N * M) * Q


def test_arrayexpr_convert_array_to_implicit_matmul():
    # Trivial dimensions are suppressed, so the result can be expressed in matrix form:

    cg = CodegenArrayTensorProduct(a, b)
    assert convert_array_to_matrix(cg) == a * b.T

    cg = CodegenArrayTensorProduct(a, I, b)
    assert convert_array_to_matrix(cg) == a * b.T

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(I, I), (1, 2))
    assert convert_array_to_matrix(cg) == I

    cg = CodegenArrayPermuteDims(CodegenArrayTensorProduct(I, Identity(1)), [0, 2, 1, 3])
    assert convert_array_to_matrix(cg) == I


def test_arrayexpr_convert_array_to_matrix_remove_trivial_dims():

    # Tensor Product:
    assert _remove_trivial_dims(CodegenArrayTensorProduct(a, b)) == (a * b.T, [1, 3])
    assert _remove_trivial_dims(CodegenArrayTensorProduct(a.T, b)) == (a * b.T, [0, 3])
    assert _remove_trivial_dims(CodegenArrayTensorProduct(a, b.T)) == (a * b.T, [1, 2])
    assert _remove_trivial_dims(CodegenArrayTensorProduct(a.T, b.T)) == (a * b.T, [0, 2])

    assert _remove_trivial_dims(CodegenArrayTensorProduct(I, a.T, b.T)) == (a * b.T, [0, 1, 2, 4])
    assert _remove_trivial_dims(CodegenArrayTensorProduct(a.T, I, b.T)) == (a * b.T, [0, 2, 3, 4])

    assert _remove_trivial_dims(CodegenArrayTensorProduct(a, I)) == (a, [2, 3])
    assert _remove_trivial_dims(CodegenArrayTensorProduct(I, a)) == (a, [0, 1])

    assert _remove_trivial_dims(CodegenArrayTensorProduct(a.T, b.T, c, d)) == (
        CodegenArrayTensorProduct(a * b.T, c * d.T), [0, 2, 5, 7])
    assert _remove_trivial_dims(CodegenArrayTensorProduct(a.T, I, b.T, c, d, I)) == (
        CodegenArrayTensorProduct(a * b.T, c * d.T, I), [0, 2, 3, 4, 7, 9])

    # Addition:

    cg = CodegenArrayElementwiseAdd(CodegenArrayTensorProduct(a, b), CodegenArrayTensorProduct(c, d))
    assert _remove_trivial_dims(cg) == (a * b.T + c * d.T, [1, 3])

    # Permute Dims:

    cg = CodegenArrayPermuteDims(CodegenArrayTensorProduct(a, b), Permutation(3)(1, 2))
    assert _remove_trivial_dims(cg) == (a * b.T, [2, 3])

    cg = CodegenArrayPermuteDims(CodegenArrayTensorProduct(a, I, b), Permutation(5)(1, 2, 3, 4))
    assert _remove_trivial_dims(cg) == (a * b.T, [1, 2, 4, 5])

    cg = CodegenArrayPermuteDims(CodegenArrayTensorProduct(I, b, a), Permutation(5)(1, 2, 4, 5, 3))
    assert _remove_trivial_dims(cg) == (b * a.T, [0, 3, 4, 5])

    # Diagonal:

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, a), (1, 2))
    assert _remove_trivial_dims(cg) == (cg, [])

    # Contraction:

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, a), (1, 2))
    assert _remove_trivial_dims(cg) == (cg, [])


def test_arrayexpr_convert_array_to_matrix_diag2contraction_diagmatrix():
    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, a), (1, 2))
    res = _array_diag2contr_diagmatrix(cg)
    assert res.shape == cg.shape
    assert res == CodegenArrayContraction(CodegenArrayTensorProduct(M, OneArray(1), DiagMatrix(a)), (1, 3))

    raises(ValueError, lambda: CodegenArrayDiagonal(CodegenArrayTensorProduct(a, M), (1, 2)))

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(a.T, M), (1, 2))
    res = _array_diag2contr_diagmatrix(cg)
    assert res.shape == cg.shape
    assert res == CodegenArrayContraction(CodegenArrayTensorProduct(OneArray(1), M, DiagMatrix(a.T)), (1, 4))

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(a.T, M, N, b.T), (1, 2), (4, 7))
    res = _array_diag2contr_diagmatrix(cg)
    assert res.shape == cg.shape
    assert res == CodegenArrayContraction(
        CodegenArrayTensorProduct(OneArray(1), M, N, OneArray(1), DiagMatrix(a.T), DiagMatrix(b.T)), (1, 7), (3, 9))

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(a, M, N, b.T), (0, 2), (4, 7))
    res = _array_diag2contr_diagmatrix(cg)
    assert res.shape == cg.shape
    assert res == CodegenArrayContraction(
        CodegenArrayTensorProduct(OneArray(1), M, N, OneArray(1), DiagMatrix(a), DiagMatrix(b.T)), (1, 6), (3, 9))

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(a, M, N, b.T), (0, 4), (3, 7))
    res = _array_diag2contr_diagmatrix(cg)
    assert res.shape == cg.shape
    assert res == CodegenArrayContraction(
        CodegenArrayTensorProduct(OneArray(1), M, N, OneArray(1), DiagMatrix(a), DiagMatrix(b.T)), (3, 6), (2, 9))

    I1 = Identity(1)
    x = MatrixSymbol("x", k, 1)
    A = MatrixSymbol("A", k, k)
    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(x, A.T, I1), (0, 2))
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

    assert _support_function_tp1_recognize([(1, 2), (5, 6)], [A, B, C, D]) == CodegenArrayTensorProduct(A * B, C * D)
    assert _support_function_tp1_recognize([(1, 4), (3, 6)], [A, B, C, D]) == CodegenArrayPermuteDims(
        CodegenArrayTensorProduct(A * C, B * D), [0, 2, 1, 3])

    assert _support_function_tp1_recognize([(0, 3), (1, 4)], [A, B, C]) == B * A * C

    assert _support_function_tp1_recognize([(9, 10), (1, 2), (5, 6), (3, 4), (7, 8)],
                                           [X, Y, A, B, C, D]) == X * Y * A * B * C * D

    assert _support_function_tp1_recognize([(9, 10), (1, 2), (5, 6), (3, 4)],
                                           [X, Y, A, B, C, D]) == CodegenArrayTensorProduct(X * Y * A * B, C * D)

    assert _support_function_tp1_recognize([(1, 7), (3, 8), (4, 11)], [X, Y, A, B, C, D]) == CodegenArrayPermuteDims(
        CodegenArrayTensorProduct(X * B.T, Y * C, D * A), [0, 2, 5, 1, 3, 4]
    )

    assert _support_function_tp1_recognize([(0, 1), (3, 6), (5, 8)], [X, A, B, C, D]) == CodegenArrayPermuteDims(
        CodegenArrayTensorProduct(Trace(X) * A * C, B * D), [0, 2, 1, 3])

    assert _support_function_tp1_recognize([(1, 2), (3, 4), (5, 6), (7, 8)], [A, A, B, C, D]) == A ** 2 * B * C * D
    assert _support_function_tp1_recognize([(1, 2), (3, 4), (5, 6), (7, 8)], [X, A, B, C, D]) == X * A * B * C * D

    assert _support_function_tp1_recognize([(1, 6), (3, 8), (5, 10)], [X, Y, A, B, C, D]) == CodegenArrayPermuteDims(
        CodegenArrayTensorProduct(X * B, Y * C, A * D), [0, 2, 4, 1, 3, 5]
    )

    assert _support_function_tp1_recognize([(1, 4), (3, 6)], [A, B, C, D]) == CodegenArrayPermuteDims(
        CodegenArrayTensorProduct(A * C, B * D), [0, 2, 1, 3])

    assert _support_function_tp1_recognize([(0, 4), (1, 7), (2, 5), (3, 8)], [X, A, B, C, D]) == C*X.T*B*A*D

    assert _support_function_tp1_recognize([(0, 4), (1, 7), (2, 5), (3, 8)], [X, A, B, C, D]) == C*X.T*B*A*D
