r"""
Array expressions are expressions representing N-dimensional arrays, without evaluating them.
They are *lazy* as they don't evaluate the operations on the arrays until required and
can define operations on unknown array (in a certain way reminiscent of graphs of
operations in TensorFlow).

Every N-dimensional array operator has a corresponding lazy object.

Table of correspondences:

=============================== =============================
         Array operator           Array expression operator
=============================== =============================
        tensorproduct                 ArrayTensorProduct
        tensorcontraction             ArrayContraction
        tensordiagonal                ArrayDiagonal
        permutedims                   PermuteDims
        derive_by_array               array_derive
=============================== =============================

Examples
========

>>> from sympy import MatrixSymbol
>>> from sympy.tensor.array.expressions import ArrayTensorProduct, ArrayContraction
>>> M = MatrixSymbol("M", 3, 3)
>>> N = MatrixSymbol("N", 3, 3)

Express the matrix product in the array expression form:

>>> from sympy.tensor.array.expressions import convert_matrix_to_array
>>> expr = convert_matrix_to_array(M*N)
>>> expr
ArrayContraction(ArrayTensorProduct(M, N), (1, 2))

The expression can be converted back to matrix form:

>>> from sympy.tensor.array.expressions import convert_array_to_matrix
>>> convert_array_to_matrix(expr)
M*N

Add a second contraction on the remaining axes in order to get the trace of `M \cdot N`:

>>> expr_tr = ArrayContraction(expr, (0, 1))
>>> expr_tr
ArrayContraction(ArrayContraction(ArrayTensorProduct(M, N), (1, 2)), (0, 1))

Flatten the expression by calling ``.doit()`` and remove the nested array contraction operations:

>>> expr_tr.doit()
ArrayContraction(ArrayTensorProduct(M, N), (0, 3), (1, 2))

Get the explicit form of the array expression:

>>> expr.as_explicit()
[[M[0, 0]*N[0, 0] + M[0, 1]*N[1, 0] + M[0, 2]*N[2, 0], M[0, 0]*N[0, 1] + M[0, 1]*N[1, 1] + M[0, 2]*N[2, 1], M[0, 0]*N[0, 2] + M[0, 1]*N[1, 2] + M[0, 2]*N[2, 2]], [M[1, 0]*N[0, 0] + M[1, 1]*N[1, 0] + M[1, 2]*N[2, 0], M[1, 0]*N[0, 1] + M[1, 1]*N[1, 1] + M[1, 2]*N[2, 1], M[1, 0]*N[0, 2] + M[1, 1]*N[1, 2] + M[1, 2]*N[2, 2]], [M[2, 0]*N[0, 0] + M[2, 1]*N[1, 0] + M[2, 2]*N[2, 0], M[2, 0]*N[0, 1] + M[2, 1]*N[1, 1] + M[2, 2]*N[2, 1], M[2, 0]*N[0, 2] + M[2, 1]*N[1, 2] + M[2, 2]*N[2, 2]]]

Express the trace of a matrix:

>>> from sympy import Trace
>>> convert_matrix_to_array(Trace(M))
ArrayContraction(M, (0, 1))
>>> convert_matrix_to_array(Trace(M*N))
ArrayContraction(ArrayTensorProduct(M, N), (0, 3), (1, 2))

Express the transposition of a matrix (will be expressed as a permutation of the axes:

>>> convert_matrix_to_array(M.T)
PermuteDims(M, (0 1))

Compute the derivative array expressions:

>>> from sympy.tensor.array.expressions import array_derive
>>> d = array_derive(M, M)
>>> d
PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))

Verify that the derivative corresponds to the form computed with explicit matrices:

>>> d.as_explicit()
[[[[1, 0, 0], [0, 0, 0], [0, 0, 0]], [[0, 1, 0], [0, 0, 0], [0, 0, 0]], [[0, 0, 1], [0, 0, 0], [0, 0, 0]]], [[[0, 0, 0], [1, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 1, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, 1], [0, 0, 0]]], [[[0, 0, 0], [0, 0, 0], [1, 0, 0]], [[0, 0, 0], [0, 0, 0], [0, 1, 0]], [[0, 0, 0], [0, 0, 0], [0, 0, 1]]]]
>>> Me = M.as_explicit()
>>> Me.diff(Me)
[[[[1, 0, 0], [0, 0, 0], [0, 0, 0]], [[0, 1, 0], [0, 0, 0], [0, 0, 0]], [[0, 0, 1], [0, 0, 0], [0, 0, 0]]], [[[0, 0, 0], [1, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 1, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, 1], [0, 0, 0]]], [[[0, 0, 0], [0, 0, 0], [1, 0, 0]], [[0, 0, 0], [0, 0, 0], [0, 1, 0]], [[0, 0, 0], [0, 0, 0], [0, 0, 1]]]]

"""

__all__ = [
    "ArraySymbol", "ArrayElement", "ZeroArray", "OneArray",
    "ArrayTensorProduct",
    "ArrayContraction",
    "ArrayDiagonal",
    "PermuteDims",
    "ArrayAdd",
    "ArrayElementwiseApplyFunc",
    "Reshape",
    "convert_array_to_matrix",
    "convert_matrix_to_array",
    "convert_indexed_to_array",
    "array_derive",
]

from sympy.tensor.array.expressions.array_expressions import ArrayTensorProduct, ArrayAdd, PermuteDims, ArrayDiagonal, \
    ArrayContraction, Reshape, ArraySymbol, ArrayElement, ZeroArray, OneArray, ArrayElementwiseApplyFunc
from sympy.tensor.array.expressions.arrayexpr_derivatives import array_derive
from sympy.tensor.array.expressions.conv_array_to_matrix import convert_array_to_matrix
from sympy.tensor.array.expressions.conv_indexed_to_array import convert_indexed_to_array
from sympy.tensor.array.expressions.conv_matrix_to_array import convert_matrix_to_array
