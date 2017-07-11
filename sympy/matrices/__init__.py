"""A module that handles matrices.

Includes functions for fast creating matrices like zero, one/eye, random
matrix, etc.
"""
from .dense import GramSchmidt, Matrix, casoratian, diag, eye, hessian, \
    jordan_cell, list2numpy, matrix2numpy, matrix_multiply_elementwise, ones, \
    randMatrix, rot_axis1, rot_axis2, rot_axis3, symarray, wronskian, zeros
from .expressions import Adjoint, BlockDiagMatrix, BlockMatrix, Determinant, \
    DiagonalMatrix, DiagonalOf, DotProduct, FunctionMatrix, HadamardProduct, \
    Identity, Inverse, MatAdd, MatMul, MatPow, MatrixExpr, MatrixSlice, \
    MatrixSymbol, Trace, Transpose, ZeroMatrix, block_collapse, blockcut, \
    det, hadamard_product, matrix_symbols, trace
from .immutable import ImmutableDenseMatrix, ImmutableMatrix, \
    ImmutableSparseMatrix
from .matrices import DeferredVector, MatrixBase, NonSquareMatrixError, \
    ShapeError
from .sparse import MutableSparseMatrix

MutableDenseMatrix = MutableMatrix = Matrix


SparseMatrix = MutableSparseMatrix


MutableSparseMatrix = SparseMatrix
