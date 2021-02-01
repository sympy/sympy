from sympy import MatrixSymbol, symbols, Identity
from sympy.codegen.array_utils import CodegenArrayTensorProduct, CodegenArrayPermuteDims
from sympy.tensor.array.expressions.arrayexpr_derivatives import array_derive

k = symbols("k")

I = Identity(k)
X = MatrixSymbol("X", k, k)
x = MatrixSymbol("x", k, 1)

A = MatrixSymbol("A", k, k)
B = MatrixSymbol("B", k, k)
C = MatrixSymbol("C", k, k)
D = MatrixSymbol("D", k, k)


def test_arrayexpr_derivatives1():

    res = array_derive(X, X)
    assert res == CodegenArrayPermuteDims(CodegenArrayTensorProduct(I, I), [0, 2, 1, 3])

    cg = CodegenArrayTensorProduct(A, X, B)
    res = array_derive(cg, X)
    assert res == CodegenArrayPermuteDims(
        CodegenArrayTensorProduct(I, A, I, B),
        [0, 4, 2, 3, 1, 5, 6, 7])
