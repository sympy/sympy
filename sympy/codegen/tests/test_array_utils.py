from sympy import symbols, IndexedBase
from sympy.codegen.array_utils import (CodegenArrayContraction,
        CodegenArrayTensorProduct, CodegenArrayDiagonal,
        CodegenArrayPermuteDims, CodegenArrayElementwiseAdd,
        _codegen_array_parse, _recognize_matrix_expression, _RecognizeMatAdd)
from sympy import (MatrixSymbol, Sum)
from sympy.combinatorics import Permutation
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy.matrices.expressions.matexpr import MatrixElement
from sympy.matrices import Trace

A, B = symbols("A B", cls=IndexedBase)
i, j, k, l, m, n = symbols("i j k l m n")

M = MatrixSymbol("M", k, k)
N = MatrixSymbol("N", k, k)
P = MatrixSymbol("P", k, k)
Q = MatrixSymbol("Q", k, k)


def test_codegen_array_contraction_construction():
    s = Sum(A[i]*B[i], (i, 0, 3))
    cg = CodegenArrayContraction.from_summation(s)
    assert cg == CodegenArrayContraction(CodegenArrayTensorProduct(A, B), (0, 1))

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(A, B), (1, 0))
    assert cg == CodegenArrayContraction(CodegenArrayTensorProduct(A, B), (0, 1))

    expr = M*N
    result = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (1, 2))
    assert CodegenArrayContraction.from_MatMul(expr) == result
    elem = expr[i, j]
    assert CodegenArrayContraction.from_summation(elem) == result

    expr = M*N*M
    result = CodegenArrayContraction(CodegenArrayTensorProduct(M, N, M), (1, 2), (3, 4))
    assert CodegenArrayContraction.from_MatMul(expr) == result
    elem = expr[i, j]
    result = CodegenArrayContraction(CodegenArrayTensorProduct(M, M, N), (1, 4), (2, 5))
    cg = CodegenArrayContraction.from_summation(elem)
    cg = cg.sort_args_by_name()
    assert cg == result


def test_codegen_array_contraction_indices_types():
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
    assert indtup == [[(0, 1), (2, 0)], [(1, 0), (2, 1)]]
    assert cg._contraction_tuples_to_contraction_indices(cg.expr, indtup) == [(1, 4), (2, 5)]


def test_codegen_array_recognize_matrix_mul_lines():

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M), (0, 1))
    assert cg._recognize_matrix_mul_lines() == [(None, None, [Trace(M)])]

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (0, 1), (2, 3))
    assert cg._recognize_matrix_mul_lines() == [(None, None, [Trace(M)]), (None, None, [Trace(N)])]

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (0, 3), (1, 2))
    assert cg._recognize_matrix_mul_lines() == [(None, None, [Trace(M*N)])]

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (0, 2), (1, 3))
    assert cg._recognize_matrix_mul_lines() == [(None, None, [Trace(M*N.T)])]

    cg = CodegenArrayContraction.from_summation((M*N*P)[i,j])
    assert cg._recognize_matrix_mul_lines() == [(i, j, [M, N, P])]
    cg = CodegenArrayContraction.from_MatMul(M*N*P)
    assert cg._recognize_matrix_mul_lines() == [(0, 5, [M, N, P])]

    cg = CodegenArrayContraction.from_summation((M*N.T*P)[i,j])
    assert cg._recognize_matrix_mul_lines() == [(i, j, [M, N.T, P])]
    cg = CodegenArrayContraction.from_MatMul(M*N.T*P)
    assert cg._recognize_matrix_mul_lines() == [(0, 5, [M, N.T, P])]

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M,N,P,Q), (1, 2), (5, 6))
    assert cg._recognize_matrix_mul_lines() == [(0, 3, [M, N]), (4, 7, [P, Q])]


def test_codegen_array_flatten():

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


def test_codegen_array_parse():
    expr = M[i, j]
    assert _codegen_array_parse(expr) == (M, (i, j))
    expr = M[i, j]*N[k, l]
    assert _codegen_array_parse(expr) == (CodegenArrayTensorProduct(M, N), (i, j, k, l))
    expr = M[i, j]*N[j, k]
    assert _codegen_array_parse(expr) == (CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N), (1, 2)), (i, k, j))
    expr = Sum(M[i, j]*N[j, k], (j, 0, k-1))
    assert _codegen_array_parse(expr) == (CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (1, 2)), (i, k))
    expr = M[i, j] + N[i, j]
    assert _codegen_array_parse(expr) == (CodegenArrayElementwiseAdd(M, N), (i, j))
    expr = M[i, j] + N[j, i]
    assert _codegen_array_parse(expr) == (CodegenArrayElementwiseAdd(M, CodegenArrayPermuteDims(N, Permutation([1,0]))), (i, j))
    expr = M[i, j] + M[j, i]
    assert _codegen_array_parse(expr) == (CodegenArrayElementwiseAdd(M, CodegenArrayPermuteDims(M, Permutation([1,0]))), (i, j))
    expr = KroneckerDelta(i, j)*M[i, k]
    assert _codegen_array_parse(expr) == (M, ({i, j}, k))
    expr = KroneckerDelta(j, k)*(M[i, j]*N[k, l] + N[i, j]*M[k, l])
    assert _codegen_array_parse(expr) == (CodegenArrayDiagonal(CodegenArrayElementwiseAdd(
            CodegenArrayTensorProduct(M, N),
            CodegenArrayPermuteDims(CodegenArrayTensorProduct(M, N), Permutation(0, 2)(1, 3))
        ), (1, 2)), (i, l, frozenset({j, k})))

def test_codegen_array_diagonal():
    cg = CodegenArrayDiagonal(M, (1, 0))
    assert cg == CodegenArrayDiagonal(M, (0, 1))

    cg = CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N, P), (4, 1), (2, 0))
    assert cg == CodegenArrayDiagonal(CodegenArrayTensorProduct(M, N, P), (1, 4), (0, 2))


def test_codegen_recognize_matrix_expression():

    expr = CodegenArrayElementwiseAdd(M, CodegenArrayPermuteDims(M, [1, 0]))
    assert _recognize_matrix_expression(expr) == _RecognizeMatAdd([M, Trace(M)])
