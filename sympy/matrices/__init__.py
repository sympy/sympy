"""A module that handles matrices.

Includes functions for fast creating matrices like zero, one/eye, random
matrix, etc.
"""
from .common import NonSquareMatrixError, ShapeError
from .dense import (GramSchmidt, MutableDenseMatrix, casoratian, diag, eye,
                    hessian, jordan_cell, list2numpy, matrix2numpy,
                    matrix_multiply_elementwise, ones, randMatrix, rot_axis1,
                    rot_axis2, rot_axis3, symarray, wronskian, zeros)
from .expressions import (Adjoint, BlockDiagMatrix, BlockMatrix, Determinant,
                          DiagonalMatrix, DiagonalOf, DotProduct,
                          FunctionMatrix, HadamardProduct, Identity, Inverse,
                          KroneckerProduct, MatAdd, MatMul, MatPow, MatrixExpr,
                          MatrixSlice, MatrixSymbol, Trace, Transpose,
                          ZeroMatrix, block_collapse, blockcut, det,
                          hadamard_product, kronecker_product, matrix_symbols,
                          trace)
from .immutable import ImmutableDenseMatrix, ImmutableSparseMatrix
from .matrices import DeferredVector, MatrixBase
from .sparse import MutableSparseMatrix

Matrix = MutableMatrix = MutableDenseMatrix
ImmutableMatrix = ImmutableDenseMatrix
SparseMatrix = MutableSparseMatrix
