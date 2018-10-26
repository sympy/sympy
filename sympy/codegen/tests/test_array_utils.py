from sympy import symbols, IndexedBase
from sympy.codegen.array_utils import CodegenArrayContraction, CodegenArrayTensorProduct
from sympy import (MatrixSymbol, Sum)
from sympy.matrices.expressions.matexpr import MatrixElement

A, B = symbols("A B", cls=IndexedBase)
i, j, k, l, m, n = symbols("i j k l m n")

M = MatrixSymbol("M", k, k)
N = MatrixSymbol("N", k, k)


def test_codegen_array_contraction_construction():
    s = Sum(A[i]*B[i], (i, 0, 3))
    cg = CodegenArrayContraction.from_summation(s)
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
    result1 = CodegenArrayContraction(CodegenArrayTensorProduct(M, M, N), (1, 4), (2, 5))
    result2 = CodegenArrayContraction(CodegenArrayTensorProduct(M, M, N), (3, 4), (1, 5))
    cg = CodegenArrayContraction.from_summation(elem)
    cg = cg.sort_args_by_name()
    assert cg in (result1, result2)


def test_codegen_array_contraction_indices_types():
    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (0, 1))
    indtup = cg._contraction_indices_to_contraction_tuples()
    assert indtup == [[(0, 0), (0, 1)]]
    assert cg._contraction_tuple_to_contraction_indices(cg.expr, indtup) == [(0, 1)]

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, N), (1, 2))
    indtup = cg._contraction_indices_to_contraction_tuples()
    assert indtup == [[(0, 1), (1, 0)]]
    assert cg._contraction_tuple_to_contraction_indices(cg.expr, indtup) == [(1, 2)]

    cg = CodegenArrayContraction(CodegenArrayTensorProduct(M, M, N), (1, 4), (2, 5))
    indtup = cg._contraction_indices_to_contraction_tuples()
    assert indtup == [[(0, 1), (2, 0)], [(1, 0), (2, 1)]]
    assert cg._contraction_tuple_to_contraction_indices(cg.expr, indtup) == [(1, 4), (2, 5)]
