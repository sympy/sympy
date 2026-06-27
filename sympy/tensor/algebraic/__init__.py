from sympy.tensor.algebraic.algebraic_tensor import (
    AlgebraicTensor,
    ShapeMismatchError,
)
from sympy.tensor.algebraic.algebraic_pure_tensor import (
    AlgebraicPureTensor,
    algebraic_tensor_product,
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
    "tensorsimplify",
]

# Backward compatibility aliases
PureTensor = AlgebraicPureTensor
ZeroTensor = AlgebraicZeroTensor
tensor_product = algebraic_tensor_product
zero_tensor = algebraic_zero_tensor

__all__.extend([
    "PureTensor",
    "ZeroTensor",
    "tensor_product",
    "zero_tensor",
])
