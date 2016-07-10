"""A module that handles matrices.

Includes functions for fast creating matrices like zero, one/eye, random
matrix, etc.
"""
from .matrices import (DeferredVector, ShapeError, NonSquareMatrixError,
    MatrixBase)

from .dense import (
    GramSchmidt, Matrix, casoratian, diag, eye, hessian, jordan_cell,
    list2numpy, matrix2numpy, matrix_multiply_elementwise, ones,
    randMatrix, rot_axis1, rot_axis2, rot_axis3, symarray, wronskian,
    zeros)

MutableDenseMatrix = MutableMatrix = Matrix

from .sparse import MutableSparseMatrix

SparseMatrix = MutableSparseMatrix

from .immutable import ImmutableMatrix, ImmutableSparseMatrix

MutableSparseMatrix = SparseMatrix
ImmutableDenseMatrix = ImmutableMatrix

from .expressions import (MatrixSlice, BlockDiagMatrix, BlockMatrix,
        FunctionMatrix, Identity, Inverse, MatAdd, MatMul, MatPow, MatrixExpr,
        MatrixSymbol, Trace, Transpose, ZeroMatrix, blockcut, block_collapse,
        matrix_symbols, Adjoint, hadamard_product, HadamardProduct,
        Determinant, det, DiagonalMatrix, DiagonalOf, trace, DotProduct)
