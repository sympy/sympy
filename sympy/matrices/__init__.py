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

from .expressions.slice import MatrixSlice
from .expressions.blockmatrix import BlockDiagMatrix, BlockMatrix
from .expressions.funcmatrix import FunctionMatrix
from .expressions.matexpr import Identity
from .expressions.inverse import Inverse
from .expressions.matadd import MatAdd
from .expressions.matmul import MatMul
from .expressions.matpow import MatPow
from .expressions.matexpr import MatrixExpr, MatrixSymbol
from .expressions.trace import Trace
from .expressions.transpose import Transpose
from .expressions.matexpr import ZeroMatrix
from .expressions.blockmatrix import blockcut, block_collapse
from .expressions.matexpr import matrix_symbols
from .expressions.adjoint import Adjoint
from .expressions.hadamard import hadamard_product, HadamardProduct
from .expressions.determinant import Determinant, det
from .expressions.diagonal import DiagonalMatrix, DiagonalOf
from .expressions.trace import trace
