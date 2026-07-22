"""Algebraic tensor expressions built from matrix factors.

This module provides classes for representing tensor products of matrix-like
objects and linear combinations thereof.  Tensors are built from factors
combined via the tensor product (``⊗``) and addition (``+``).  A third
operation, composition (``*``), performs factor-wise matrix multiplication
and serves as the tensor analogue of matrix multiplication. The multiplication
(``*``) can also be used to multiply tensors with scalars (commutative Symbols).

The module provides three core types:

    - :class:`~sympy.tensor.algebraic.algebraic_pure_tensor.AlgebraicPureTensor`
      A single tensor-product term with internal coefficient storage.

    - :class:`~sympy.tensor.algebraic.algebraic_tensor.AlgebraicTensor`
      A linear combination of same-shape pure tensors.

    - :class:`~sympy.tensor.algebraic.algebraic_zero_tensor.AlgebraicZeroTensor`
      The additive identity (zero tensor) for a given tensor shape.

Examples
========

Create a pure tensor from two matrix symbols:

>>> from sympy.matrices.expressions import MatrixSymbol
>>> from sympy.tensor.algebraic import AlgebraicPureTensor
>>> A = MatrixSymbol("A", 3, 4)
>>> B = MatrixSymbol("B", 4, 5)
>>> T = AlgebraicPureTensor(A, B)
>>> T.shape
((3, 4), (4, 5))

Form a sum of two pure tensors with the same shape:

>>> C = MatrixSymbol("C", 3, 4)
>>> D = MatrixSymbol("D", 4, 5)
>>> from sympy.tensor.algebraic import AlgebraicTensor
>>> S = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
>>> S.shape
((3, 4), (4, 5))

Scale a pure tensor by a numeric or symbolic coefficient:

>>> from sympy.tensor.algebraic import AlgebraicPureTensor
>>> T2 = AlgebraicPureTensor(2, A, B)
>>> T2
2*A ⊗ B
>>> from sympy.abc import x
>>> AlgebraicPureTensor(x, A, B)
x*A ⊗ B

Create a zero tensor for a given shape:

>>> from sympy.tensor.algebraic import algebraic_zero_tensor
>>> Z = algebraic_zero_tensor(((3, 4), (4, 5)))
>>> Z
0_{(3x4), (4x5)}
"""

from sympy.tensor.algebraic.algebraic_tensor import (
    AlgebraicTensor,
    algebraic_tensor_product,
    compose_algebraic_tensors,
)
from sympy.tensor.algebraic.algebraic_pure_tensor import (
    AlgebraicPureTensor,
    ShapeMismatchError,
    compose_algebraic_pure_tensors,
)
from sympy.tensor.algebraic.algebraic_zero_tensor import (
    AlgebraicZeroTensor,
    algebraic_zero_tensor,
)
from sympy.tensor.algebraic.simplify import tensorsimplify

__all__ = [
    "AlgebraicTensor",
    "AlgebraicPureTensor",
    "AlgebraicZeroTensor",
    "ShapeMismatchError",
    "algebraic_tensor_product",
    "algebraic_zero_tensor",
    "compose_algebraic_pure_tensors",
    "compose_algebraic_tensors",
    "tensorsimplify",
]
