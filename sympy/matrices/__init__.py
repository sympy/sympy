"""A module that handles matrices.

Includes functions for fast creating matrices like zero, one/eye, random
matrix, etc.
"""
from matrices import (SparseMatrix, zeros, ones, eye, diag,
     hessian, randMatrix, GramSchmidt, wronskian, casoratian,
     list2numpy, matrix2numpy, DeferredVector, symarray, ShapeError,
     NonSquareMatrixError, rot_axis1, rot_axis2, rot_axis3)

from matrices import MutableMatrix as Matrix

from immutable_matrix import ImmutableMatrix

from expressions import (BlockDiagMatrix, BlockMatrix, FunctionMatrix,
     Identity, Inverse, MatAdd, MatMul, MatPow, MatrixExpr, MatrixSymbol,
     Trace, Transpose, ZeroMatrix, block_collapse, linear_factors,
     matrix_symbols, matrixify)
