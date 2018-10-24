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
    result = CodegenArrayContraction(CodegenArrayTensorProduct(M, M, N), (1, 4), (2, 5))
    assert CodegenArrayContraction.from_summation(elem) == result
