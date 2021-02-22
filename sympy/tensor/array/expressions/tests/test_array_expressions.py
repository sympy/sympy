import random

from sympy import symbols, ImmutableDenseNDimArray, tensorproduct, tensorcontraction, permutedims, MatrixSymbol, \
    ZeroMatrix, sin, cos, DiagMatrix
from sympy.combinatorics import Permutation
from sympy.tensor.array.expressions.array_expressions import ZeroArray, OneArray, ArraySymbol, ArrayElement, \
    CodegenArrayPermuteDims, CodegenArrayContraction, CodegenArrayTensorProduct, CodegenArrayDiagonal, \
    CodegenArrayElementwiseAdd, nest_permutation, ArrayElementwiseApplyFunc
from sympy.testing.pytest import raises

i, j, k, l, m, n = symbols("i j k l m n")


M = ArraySymbol("M", k, k)
N = ArraySymbol("N", k, k)
P = ArraySymbol("P", k, k)
Q = ArraySymbol("Q", k, k)

A = ArraySymbol("A", k, k)
B = ArraySymbol("B", k, k)
C = ArraySymbol("C", k, k)
D = ArraySymbol("D", k, k)

X = ArraySymbol("X", k, k)
Y = ArraySymbol("Y", k, k)

a = ArraySymbol("a", k, 1)
b = ArraySymbol("b", k, 1)
c = ArraySymbol("c", k, 1)
d = ArraySymbol("d", k, 1)


def test_array_symbol_and_element():
    A = ArraySymbol("A", 2)
    A0 = ArrayElement(A, (0,))
    A1 = ArrayElement(A, (1,))
    assert A.as_explicit() == ImmutableDenseNDimArray([A0, A1])

    A2 = tensorproduct(A, A)
    assert A2.shape == (2, 2)
    # TODO: not yet supported:
    # assert A2.as_explicit() == Array([[A[0]*A[0], A[1]*A[0]], [A[0]*A[1], A[1]*A[1]]])
    A3 = tensorcontraction(A2, (0, 1))
    assert A3.shape == ()
    # TODO: not yet supported:
    # assert A3.as_explicit() == Array([])

    A = ArraySymbol("A", 2, 3, 4)
    Ae = A.as_explicit()
    assert Ae == ImmutableDenseNDimArray(
        [[[ArrayElement(A, (i, j, k)) for k in range(4)] for j in range(3)] for i in range(2)])

    p = permutedims(A, Permutation(0, 2, 1))
    assert isinstance(p, CodegenArrayPermuteDims)


def test_zero_array():
    assert ZeroArray() == 0
    assert ZeroArray().is_Integer

    za = ZeroArray(3, 2, 4)
    assert za.shape == (3, 2, 4)
    za_e = za.as_explicit()
    assert za_e.shape == (3, 2, 4)

    m, n, k = symbols("m n k")
    za = ZeroArray(m, n, k, 2)
    assert za.shape == (m, n, k, 2)
    raises(ValueError, lambda: za.as_explicit())


def test_one_array():
    assert OneArray() == 1
    assert OneArray().is_Integer

    oa = OneArray(3, 2, 4)
    assert oa.shape == (3, 2, 4)
    oa_e = oa.as_explicit()
    assert oa_e.shape == (3, 2, 4)

    m, n, k = symbols("m n k")
    oa = OneArray(m, n, k, 2)
    assert oa.shape == (m, n, k, 2)
    raises(ValueError, lambda: oa.as_explicit())


def test_arrayexpr_contraction_construction():

    cg = CodegenArrayContraction(A)
    assert cg == A

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, B), (1, 0))
    assert cg == CodegenArrayContraction(CodegenArrayTensorProduct(A, B), (0, 1))

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (0, 1))
    indtup = cg._get_contraction_tuples()
    assert indtup == [[(0, 0), (0, 1)]]
    assert cg._contraction_tuples_to_contraction_indices(cg.expr, indtup) == [(0, 1)]

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (1, 2))
    indtup = cg._get_contraction_tuples()
    assert indtup == [[(0, 1), (1, 0)]]
    assert cg._contraction_tuples_to_contraction_indices(cg.expr, indtup) == [(1, 2)]

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, M, N), (1, 4), (2, 5))
    indtup = cg._get_contraction_tuples()
    assert indtup == [[(0, 0), (1, 1)], [(0, 1), (2, 0)]]
    assert cg._contraction_tuples_to_contraction_indices(cg.expr, indtup) == [(0, 3), (1, 4)]


def test_arrayexpr_array_flatten():

    # Flatten nested CodegenArrayTensorProduct objects:
    expr1 = CodegenArrayTensorProduct(M, N)
    expr2 = CodegenArrayTensorProduct(P, Q)
    expr = CodegenArrayTensorProduct(expr1, expr2)
    assert expr == CodegenArrayTensorProduct(M, N, P, Q)
    assert expr.args == (M, N, P, Q)

    # Flatten mixed CodegenArrayTensorProduct and CodegenArrayContraction objects:
    cg1 = CodegenArrayContraction(expr1, (1, 2))
    cg2 = CodegenArrayContraction(expr2, (0, 3))

    expr = CodegenArrayTensorProduct(cg1, cg2)
    assert expr == CodegenArrayContraction(CodegenArrayTensorProduct(M, N, P, Q), (1, 2), (4, 7))

    expr = CodegenArrayTensorProduct(M, cg1)
    assert expr == CodegenArrayContraction(CodegenArrayTensorProduct(M, M, N), (3, 4))

    # Flatten nested CodegenArrayContraction objects:
    cgnested = CodegenArrayContraction(cg1, (0, 1))
    assert cgnested == CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (0, 3), (1, 2))

    cgnested = CodegenArrayContraction(CodegenArrayTensorProduct(cg1, cg2), (0, 3))
    assert cgnested == CodegenArrayContraction(CodegenArrayTensorProduct(M, N, P, Q), (0, 6), (1, 2), (4, 7))

    cg3 = CodegenArrayContraction(CodegenArrayTensorProduct(M, N, P, Q), (1, 3), (2, 4))
    cgnested = CodegenArrayContraction(cg3, (0, 1))
    assert cgnested == CodegenArrayContraction(CodegenArrayTensorProduct(M, N, P, Q), (0, 5), (1, 3), (2, 4))

    cgnested = CodegenArrayContraction(cg3, (0, 3), (1, 2))
    assert cgnested == CodegenArrayContraction(CodegenArrayTensorProduct(M, N, P, Q), (0, 7), (1, 3), (2, 4), (5, 6))

    cg4 = CodegenArrayContraction(CodegenArrayTensorProduct(M, N, P, Q), (1, 5), (3, 7))
    cgnested = CodegenArrayContraction(cg4, (0, 1))
    assert cgnested == CodegenArrayContraction(CodegenArrayTensorProduct(M, N, P, Q), (0, 2), (1, 5), (3, 7))

    cgnested = CodegenArrayContraction(cg4, (0, 1), (2, 3))
    assert cgnested == CodegenArrayContraction(CodegenArrayTensorProduct(M, N, P, Q), (0, 2), (1, 5), (3, 7), (4, 6))

    cg = CodegenArrayDiagonal(cg4)
    assert cg == cg4
    assert isinstance(cg, type(cg4))

    # Flatten nested CodegenArrayDiagonal objects:
    cg1 = CodegenArrayDiagonal(expr1, (1, 2))
    cg2 = CodegenArrayDiagonal(expr2, (0, 3))
    cg3 = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N, P, Q), (1, 3), (2, 4))
    cg4 = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N, P, Q), (1, 5), (3, 7))

    cgnested = CodegenArrayDiagonal(cg1, (0, 1))
    assert cgnested == CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N), (1, 2), (0, 3))

    cgnested = CodegenArrayDiagonal(cg3, (1, 2))
    assert cgnested == CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N, P, Q), (1, 3), (2, 4), (5, 6))

    cgnested = CodegenArrayDiagonal(cg4, (1, 2))
    assert cgnested == CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N, P, Q), (1, 5), (3, 7), (2, 4))

    cg = CodegenArrayElementwiseAdd(M, N)
    cg2 = CodegenArrayElementwiseAdd(cg, P)
    assert isinstance(cg2, CodegenArrayElementwiseAdd)
    assert cg2.args == (M, N, P)
    assert cg2.shape == (k, k)


def test_arrayexpr_array_diagonal():
    cg = CodegenArrayDiagonal(M, (1, 0))
    assert cg == CodegenArrayDiagonal(M, (0, 1))

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N, P), (4, 1), (2, 0))
    assert cg == CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N, P), (1, 4), (0, 2))


def test_arrayexpr_array_shape():
    expr = CodegenArrayTensorProduct(M, N, P, Q)
    assert expr.shape == (k, k, k, k, k, k, k, k)
    Z = MatrixSymbol("Z", m, n)
    expr = CodegenArrayTensorProduct(M, Z)
    assert expr.shape == (k, k, m, n)
    expr2 = CodegenArrayContraction(expr, (0, 1))
    assert expr2.shape == (m, n)
    expr2 = CodegenArrayDiagonal(expr, (0, 1))
    assert expr2.shape == (m, n, k)
    exprp = CodegenArrayPermuteDims(expr, [2, 1, 3, 0])
    assert exprp.shape == (m, k, n, k)
    expr3 = CodegenArrayTensorProduct(N, Z)
    expr2 = CodegenArrayElementwiseAdd(expr, expr3)
    assert expr2.shape == (k, k, m, n)

    # Contraction along axes with discordant dimensions:
    raises(ValueError, lambda: CodegenArrayContraction(expr, (1, 2)))
    # Also diagonal needs the same dimensions:
    raises(ValueError, lambda: CodegenArrayDiagonal(expr, (1, 2)))
    # Diagonal requires at least to axes to compute the diagonal:
    raises(ValueError, lambda: CodegenArrayDiagonal(expr, (1,)))


def test_arrayexpr_permutedims_sink():

    cg = CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N), [0, 1, 3, 2], nest_permutation=False)
    sunk = nest_permutation(cg)
    assert sunk == CodegenArrayTensorProduct(M, CodegenArrayPermuteDims(N, [1, 0]))

    cg = CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N), [1, 0, 3, 2], nest_permutation=False)
    sunk = nest_permutation(cg)
    assert sunk == CodegenArrayTensorProduct(CodegenArrayPermuteDims(M, [1, 0]), CodegenArrayPermuteDims(N, [1, 0]))

    cg = CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N), [3, 2, 1, 0], nest_permutation=False)
    sunk = nest_permutation(cg)
    assert sunk == CodegenArrayTensorProduct(CodegenArrayPermuteDims(N, [1, 0]), CodegenArrayPermuteDims(M, [1, 0]))

    cg = CodegenArrayPermuteDims(CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (1, 2)), [1, 0], nest_permutation=False)
    sunk = nest_permutation(cg)
    assert sunk == CodegenArrayContraction(CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N), [[0, 3]]), (1, 2))

    cg = CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N), [1, 0, 3, 2], nest_permutation=False)
    sunk = nest_permutation(cg)
    assert sunk == CodegenArrayTensorProduct(CodegenArrayPermuteDims(M, [1, 0]), CodegenArrayPermuteDims(N, [1, 0]))

    cg = CodegenArrayPermuteDims(CodegenArrayContraction(CodegenArrayTensorProduct(M, N, P), (1, 2), (3, 4)), [1, 0], nest_permutation=False)
    sunk = nest_permutation(cg)
    assert sunk == CodegenArrayContraction(CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N, P), [[0, 5]]), (1, 2), (3, 4))


def test_arrayexpr_push_indices_up_and_down():

    indices = list(range(12))

    contr_diag_indices = [(0, 6), (2, 8)]
    assert CodegenArrayContraction._push_indices_down(contr_diag_indices, indices) == (1, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15)
    assert CodegenArrayContraction._push_indices_up(contr_diag_indices, indices) == (None, 0, None, 1, 2, 3, None, 4, None, 5, 6, 7)

    assert CodegenArrayDiagonal._push_indices_down(contr_diag_indices, indices, 10) == (1, 3, 4, 5, 7, 9, (0, 6), (2, 8), None, None, None, None)
    assert CodegenArrayDiagonal._push_indices_up(contr_diag_indices, indices, 10) == (6, 0, 7, 1, 2, 3, 6, 4, 7, 5, None, None)

    contr_diag_indices = [(1, 2), (7, 8)]
    assert CodegenArrayContraction._push_indices_down(contr_diag_indices, indices) == (0, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15)
    assert CodegenArrayContraction._push_indices_up(contr_diag_indices, indices) == (0, None, None, 1, 2, 3, 4, None, None, 5, 6, 7)

    assert CodegenArrayDiagonal._push_indices_down(contr_diag_indices, indices, 10) == (0, 3, 4, 5, 6, 9, (1, 2), (7, 8), None, None, None, None)
    assert CodegenArrayDiagonal._push_indices_up(contr_diag_indices, indices, 10) == (0, 6, 6, 1, 2, 3, 4, 7, 7, 5, None, None)


def test_arrayexpr_split_multiple_contractions():
    a = MatrixSymbol("a", k, 1)
    b = MatrixSymbol("b", k, 1)
    A = MatrixSymbol("A", k, k)
    B = MatrixSymbol("B", k, k)
    C = MatrixSymbol("C", k, k)
    X = MatrixSymbol("X", k, k)

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(A.T, a, b, b.T, (A*X*b).applyfunc(cos)), (1, 2, 8), (5, 6, 9))
    assert cg.split_multiple_contractions().dummy_eq(CodegenArrayContraction(CodegenArrayTensorProduct(DiagMatrix(a), (A*X*b).applyfunc(cos), A.T, b, b.T), (0, 2), (1, 5), (3, 7, 8)))
    # assert recognize_matrix_expression(cg)

    # Check no overlap of lines:

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, a, C, a, B), (1, 2, 4), (5, 6, 8), (3, 7))
    assert cg.split_multiple_contractions() == cg

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(a, b, A), (0, 2, 4), (1, 3))
    assert cg.split_multiple_contractions() == cg


def test_arrayexpr_nested_permutations():

    cg = CodegenArrayPermuteDims(CodegenArrayPermuteDims(M, (1, 0)), (1, 0))
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

        cg = CodegenArrayPermuteDims(
            CodegenArrayPermuteDims(
                CodegenArrayTensorProduct(M, N, P),
                p1),
            p2
        )
        result = CodegenArrayPermuteDims(
            CodegenArrayTensorProduct(M, N, P),
            p2*p1
        )
        assert cg == result

        # Check that `permutedims` behaves the same way with explicit-component arrays:
        result1 = permutedims(permutedims(cge, p1), p2)
        result2 = permutedims(cge, p2*p1)
        assert result1 == result2


def test_arrayexpr_contraction_permutation_mix():

    Me = M.subs(k, 3).as_explicit()
    Ne = N.subs(k, 3).as_explicit()

    cg1 = CodegenArrayContraction(CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N), Permutation([0, 2, 1, 3])), (2, 3))
    cg2 = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (1, 3))
    assert cg1 == cg2
    cge1 = tensorcontraction(permutedims(tensorproduct(Me, Ne), Permutation([0, 2, 1, 3])), (2, 3))
    cge2 = tensorcontraction(tensorproduct(Me, Ne), (1, 3))
    assert cge1 == cge2

    cg1 = CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N), Permutation([0, 1, 3, 2]))
    cg2 = CodegenArrayTensorProduct(M, CodegenArrayPermuteDims(N, Permutation([1, 0])))
    assert cg1 == cg2

    cg1 = CodegenArrayContraction(
        CodegenArrayPermuteDims(
            CodegenArrayTensorProduct(M, N, P, Q), Permutation([0, 2, 3, 1, 4, 5, 7, 6])),
        (1, 2), (3, 5)
    )
    cg2 = CodegenArrayContraction(
        CodegenArrayTensorProduct(M, N, P, CodegenArrayPermuteDims(Q, Permutation([1, 0]))),
        (1, 5), (2, 3)
    )
    assert cg1 == cg2

    cg1 = CodegenArrayContraction(
        CodegenArrayPermuteDims(
            CodegenArrayTensorProduct(M, N, P, Q), Permutation([1, 0, 4, 6, 2, 7, 5, 3])),
        (0, 1), (2, 6), (3, 7)
    )
    cg2 = CodegenArrayPermuteDims(
        CodegenArrayContraction(
            CodegenArrayTensorProduct(M, P, Q, N),
            (0, 1), (2, 3), (4, 7)),
        [1, 0]
    )
    assert cg1 == cg2

    cg1 = CodegenArrayContraction(
        CodegenArrayPermuteDims(
            CodegenArrayTensorProduct(M, N, P, Q), Permutation([1, 0, 4, 6, 7, 2, 5, 3])),
        (0, 1), (2, 6), (3, 7)
    )
    cg2 = CodegenArrayPermuteDims(
        CodegenArrayContraction(
            CodegenArrayTensorProduct(CodegenArrayPermuteDims(M, [1, 0]), N, P, Q),
            (0, 1), (3, 6), (4, 5)
        ),
        Permutation([1, 0])
    )
    assert cg1 == cg2


def test_arrayexpr_permute_tensor_product():
    cg1 = CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N, P, Q), Permutation([2, 3, 1, 0, 5, 4, 6, 7]))
    cg2 = CodegenArrayTensorProduct(N, CodegenArrayPermuteDims(M, [1, 0]),
                                    CodegenArrayPermuteDims(P, [1, 0]), Q)
    assert cg1 == cg2

    # TODO: reverse operation starting with `CodegenArrayPermuteDims` and getting down to `bb`...
    cg1 = CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N, P, Q), Permutation([2, 3, 4, 5, 0, 1, 6, 7]))
    cg2 = CodegenArrayTensorProduct(N, P, M, Q)
    assert cg1 == cg2

    cg1 = CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N, P, Q), Permutation([2, 3, 4, 6, 5, 7, 0, 1]))
    assert cg1.expr == CodegenArrayTensorProduct(N, P, Q, M)
    assert cg1.permutation == Permutation([0, 1, 2, 4, 3, 5, 6, 7])

    cg1 = CodegenArrayContraction(
        CodegenArrayPermuteDims(
            CodegenArrayTensorProduct(N, Q, Q, M),
            [2, 1, 5, 4, 0, 3, 6, 7]),
        [1, 2, 6])
    cg2 = CodegenArrayPermuteDims(CodegenArrayContraction(CodegenArrayTensorProduct(Q, Q, N, M), (3, 5, 6)), [0, 2, 3, 1, 4])
    assert cg1 == cg2

    cg1 = CodegenArrayContraction(
        CodegenArrayContraction(
            CodegenArrayContraction(
                CodegenArrayContraction(
                    CodegenArrayPermuteDims(
                        CodegenArrayTensorProduct(N, Q, Q, M),
                        [2, 1, 5, 4, 0, 3, 6, 7]),
                    [1, 2, 6]),
                [1, 3, 4]),
            [1]),
        [0])
    cg2 = CodegenArrayContraction(CodegenArrayTensorProduct(M, N, Q, Q), (0, 3, 5), (1, 4, 7), (2,), (6,))
    assert cg1 == cg2


def test_arrayexpr_normalize_diagonal_permutedims():
    tp = CodegenArrayTensorProduct(M, Q, N, P)
    expr = CodegenArrayDiagonal(
        CodegenArrayPermuteDims(tp, [0, 1, 2, 4, 7, 6, 3, 5]), (2, 4, 5), (6, 7),
        (0, 3))
    result = CodegenArrayDiagonal(tp, (2, 6, 7), (3, 5), (0, 4))
    assert expr == result

    tp = CodegenArrayTensorProduct(M, N, P, Q)
    expr = CodegenArrayDiagonal(CodegenArrayPermuteDims(tp, [0, 5, 2, 4, 1, 6, 3, 7]), (1, 2, 6), (3, 4))
    result = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, P, N, Q), (3, 4, 5), (1, 2))
    assert expr == result


def test_arrayexpr_normalize_diagonal_contraction():
    tp = CodegenArrayTensorProduct(M, N, P, Q)
    expr = CodegenArrayContraction(CodegenArrayDiagonal(tp, (1, 3, 4)), (0, 3))
    result = CodegenArrayDiagonal(CodegenArrayContraction(CodegenArrayTensorProduct(M, N, P, Q), (0, 6)), (0, 2, 3))
    assert expr == result

    expr = CodegenArrayContraction(CodegenArrayDiagonal(tp, (0, 1, 2, 3, 7)), (1, 2, 3))
    result = CodegenArrayContraction(CodegenArrayTensorProduct(M, N, P, Q), (0, 1, 2, 3, 5, 6, 7))
    assert expr == result

    expr = CodegenArrayContraction(CodegenArrayDiagonal(tp, (0, 2, 6, 7)), (1, 2, 3))
    result = CodegenArrayDiagonal(CodegenArrayContraction(tp, (3, 4, 5)), (0, 2, 3, 4))
    assert expr == result

    td = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N, P, Q), (0, 3))
    expr = CodegenArrayContraction(td, (2, 1), (0, 4, 6, 5, 3))
    result = CodegenArrayContraction(CodegenArrayTensorProduct(M, N, P, Q), (0, 1, 3, 5, 6, 7), (2, 4))
    assert expr == result


def test_arrayexpr_array_wrong_permutation_size():
    cg = CodegenArrayTensorProduct(M, N)
    raises(ValueError, lambda: CodegenArrayPermuteDims(cg, [1, 0]))
    raises(ValueError, lambda: CodegenArrayPermuteDims(cg, [1, 0, 2, 3, 5, 4]))


def test_arrayexpr_nested_array_elementwise_add():
    cg = CodegenArrayContraction(CodegenArrayElementwiseAdd(
        CodegenArrayTensorProduct(M, N),
        CodegenArrayTensorProduct(N, M)
    ), (1, 2))
    result = CodegenArrayElementwiseAdd(
        CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (1, 2)),
        CodegenArrayContraction(CodegenArrayTensorProduct(N, M), (1, 2))
    )
    assert cg == result

    cg = CodegenArrayDiagonal(CodegenArrayElementwiseAdd(
        CodegenArrayTensorProduct(M, N),
        CodegenArrayTensorProduct(N, M)
    ), (1, 2))
    result = CodegenArrayElementwiseAdd(
        CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N), (1, 2)),
        CodegenArrayDiagonal(CodegenArrayTensorProduct(N, M), (1, 2))
    )
    assert cg == result


def test_arrayexpr_array_expr_zero_array():
    za1 = ZeroArray(k, l, m, n)
    zm1 = ZeroMatrix(m, n)

    za2 = ZeroArray(k, m, m, n)
    zm2 = ZeroMatrix(m, m)
    zm3 = ZeroMatrix(k, k)

    assert CodegenArrayTensorProduct(M, N, za1) == ZeroArray(k, k, k, k, k, l, m, n)
    assert CodegenArrayTensorProduct(M, N, zm1) == ZeroArray(k, k, k, k, m, n)

    assert CodegenArrayContraction(za1, (3,)) == ZeroArray(k, l, m)
    assert CodegenArrayContraction(zm1, (1,)) == ZeroArray(m)
    assert CodegenArrayContraction(za2, (1, 2)) == ZeroArray(k, n)
    assert CodegenArrayContraction(zm2, (0, 1)) == 0

    assert CodegenArrayDiagonal(za2, (1, 2)) == ZeroArray(k, n, m)
    assert CodegenArrayDiagonal(zm2, (0, 1)) == ZeroArray(m)

    assert CodegenArrayPermuteDims(za1, [2, 1, 3, 0]) == ZeroArray(m, l, n, k)
    assert CodegenArrayPermuteDims(zm1, [1, 0]) == ZeroArray(n, m)

    assert CodegenArrayElementwiseAdd(za1) == za1
    assert CodegenArrayElementwiseAdd(zm1) == ZeroArray(m, n)
    tp1 = CodegenArrayTensorProduct(MatrixSymbol("A", k, l), MatrixSymbol("B", m, n))
    assert CodegenArrayElementwiseAdd(tp1, za1) == tp1
    tp2 = CodegenArrayTensorProduct(MatrixSymbol("C", k, l), MatrixSymbol("D", m, n))
    assert CodegenArrayElementwiseAdd(tp1, za1, tp2) == CodegenArrayElementwiseAdd(tp1, tp2)
    assert CodegenArrayElementwiseAdd(M, zm3) == M
    assert CodegenArrayElementwiseAdd(M, N, zm3) == CodegenArrayElementwiseAdd(M, N)


def test_arrayexpr_array_expr_applyfunc():

    A = ArraySymbol("A", 3, k, 2)
    aaf = ArrayElementwiseApplyFunc(sin, A)
    assert aaf.shape == (3, k, 2)
