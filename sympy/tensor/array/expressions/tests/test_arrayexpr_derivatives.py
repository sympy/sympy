from sympy import MatrixSymbol, symbols, Identity, sin, cos
from sympy.codegen.array_utils import CodegenArrayTensorProduct, CodegenArrayPermuteDims, CodegenArrayDiagonal, \
    CodegenArrayContraction, ArrayElementwiseApplyFunc
from sympy.matrices.expressions.applyfunc import ElementwiseApplyFunction
from sympy.tensor.array.expressions.array_expressions import ArraySymbol
from sympy.tensor.array.expressions.arrayexpr_derivatives import array_derive

k = symbols("k")

I = Identity(k)
X = MatrixSymbol("X", k, k)
x = MatrixSymbol("x", k, 1)

A = MatrixSymbol("A", k, k)
B = MatrixSymbol("B", k, k)
C = MatrixSymbol("C", k, k)
D = MatrixSymbol("D", k, k)

A1 = ArraySymbol("A", 3, 2, k)


def test_arrayexpr_derivatives1():

    res = array_derive(X, X)
    assert res == CodegenArrayPermuteDims(CodegenArrayTensorProduct(I, I), [0, 2, 1, 3])

    cg = CodegenArrayTensorProduct(A, X, B)
    res = array_derive(cg, X)
    assert res == CodegenArrayPermuteDims(
        CodegenArrayTensorProduct(I, A, I, B),
        [0, 4, 2, 3, 1, 5, 6, 7])

    cg = CodegenArrayContraction(X, (0, 1))
    res = array_derive(cg, X)
    assert res == CodegenArrayContraction(CodegenArrayTensorProduct(I, I), (1, 3))

    cg = CodegenArrayDiagonal(X, (0, 1))
    res = array_derive(cg, X)
    assert res == CodegenArrayDiagonal(CodegenArrayTensorProduct(I, I), (1, 3))

    cg = ElementwiseApplyFunction(sin, X)
    res = array_derive(cg, X)
    assert res.dummy_eq(CodegenArrayDiagonal(
        CodegenArrayTensorProduct(
            ElementwiseApplyFunction(cos, X),
            I,
            I
        ), (0, 3), (1, 5)))

    cg = ArrayElementwiseApplyFunc(sin, X)
    res = array_derive(cg, X)
    assert res.dummy_eq(CodegenArrayDiagonal(
        CodegenArrayTensorProduct(
            I,
            I,
            ArrayElementwiseApplyFunc(cos, X)
        ), (1, 4), (3, 5)))

    res = array_derive(A1, A1)
    assert res == CodegenArrayPermuteDims(
        CodegenArrayTensorProduct(Identity(3), Identity(2), Identity(k)),
        [0, 2, 4, 1, 3, 5]
    )

    cg = ArrayElementwiseApplyFunc(sin, A1)
    res = array_derive(cg, A1)
    assert res.dummy_eq(CodegenArrayDiagonal(
        CodegenArrayTensorProduct(
            Identity(3), Identity(2), Identity(k),
            ArrayElementwiseApplyFunc(cos, A1)
        ), (1, 6), (3, 7), (5, 8)
    ))
