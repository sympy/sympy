""" A module which handles Matrix Expressions """

__all__ = []

from .slice import MatrixSlice
__all__ += ["MatrixSlice"]

from .blockmatrix import (
    BlockMatrix, BlockDiagMatrix,
    block_collapse, blockcut
)
__all__ += [
    "BlockMatrix", "BlockDiagMatrix",
    "block_collapse", "blockcut"
]

from .funcmatrix import FunctionMatrix
__all__ += ["FunctionMatrix"]

from .inverse import Inverse
__all__ += ["Inverse"]

from .matadd import MatAdd
__all__ += ["MatAdd"]

from .matexpr import (
    Identity, MatrixExpr, MatrixSymbol,
    ZeroMatrix, matrix_symbols
)
__all__ += [
    "Identity", "MatrixExpr", "MatrixSymbol",
    "ZeroMatrix", "matrix_symbols"
]

from .matmul import MatMul
__all__ += ["MatMul"]

from .matpow import MatPow
__all__ += ["MatPow"]

from .trace import Trace, trace
__all__ += ["Trace", "trace"]

from .determinant import Determinant, det
__all__ += ["Determinant", "det"]

from .transpose import Transpose
__all__ += ["Transpose"]

from .adjoint import Adjoint
__all__ += ["Adjoint"]

from .hadamard import hadamard_product, HadamardProduct
__all__ += ["hadamard_product", "HadamardProduct"]

from .diagonal import DiagonalMatrix, DiagonalOf
__all__ += ["DiagonalMatrix", "DiagonalOf"]

from .dotproduct import DotProduct
__all__ += ["DotProduct"]
