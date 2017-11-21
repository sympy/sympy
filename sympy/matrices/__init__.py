"""A module that handles matrices.

Includes functions for fast creating matrices like zero, one/eye, random
matrix, etc.
"""

__all__ = []

from .matrices import (
    DeferredVector, ShapeError,
    NonSquareMatrixError, MatrixBase
)
__all__ += [
    "DeferredVector", "ShapeError",
    "NonSquareMatrixError", "MatrixBase"
]

from .dense import (
    GramSchmidt, Matrix, casoratian, diag, eye,
    hessian, jordan_cell, list2numpy, matrix2numpy,
    matrix_multiply_elementwise, ones, randMatrix,
    rot_axis1, rot_axis2, rot_axis3, symarray,
    wronskian, zeros
)
__all__ += [
    "GramSchmidt", "Matrix", "casoratian", "diag", "eye",
    "hessian", "jordan_cell", "list2numpy", "matrix2numpy",
    "matrix_multiply_elementwise", "ones", "randMatrix",
    "rot_axis1", "rot_axis2", "rot_axis3", "symarray",
    "wronskian", "zeros"
]

from .sparse import MutableSparseMatrix
__all__ += ["MutableSparseMatrix"]

from .immutable import ImmutableMatrix, ImmutableDenseMatrix, ImmutableSparseMatrix
__all__ += ["ImmutableMatrix", "ImmutableDenseMatrix", "ImmutableSparseMatrix"]

from .expressions import (
    MatrixSlice, BlockDiagMatrix, BlockMatrix,
    FunctionMatrix, Identity, Inverse, MatAdd,
    MatMul, MatPow, MatrixExpr, MatrixSymbol,
    Trace, Transpose, ZeroMatrix, blockcut,
    block_collapse, matrix_symbols, Adjoint,
    hadamard_product, HadamardProduct, Determinant,
    det, DiagonalMatrix, DiagonalOf, trace, DotProduct
)
__all__ += [
    "MatrixSlice", "BlockDiagMatrix", "BlockMatrix",
    "FunctionMatrix", "Identity", "Inverse", "MatAdd",
    "MatMul", "MatPow", "MatrixExpr", "MatrixSymbol",
    "Trace", "Transpose", "ZeroMatrix", "blockcut",
    "block_collapse", "matrix_symbols", "Adjoint",
    "hadamard_product", "HadamardProduct", "Determinant",
    "det", "DiagonalMatrix", "DiagonalOf", "trace", "DotProduct"
]

MutableDenseMatrix = MutableMatrix = Matrix
SparseMatrix = MutableSparseMatrix
__all__ += ["MutableDenseMatrix", "MutableMatrix", "SparseMatrix"]
