"""A module that handles matrices.

Includes functions for fast creating matrices like zero, one/eye, random
matrix, etc.
"""
from matrices import (DeferredVector, ShapeError, NonSquareMatrixError)

from mutable import (MutableMatrix as Matrix,
    list2numpy, matrix2numpy, symarray,
    rot_axis1, rot_axis2, rot_axis3,
    zeros, ones, eye, diag, jordan_cell,
    hessian, GramSchmidt, wronskian, casoratian,
    randMatrix, matrix_multiply_elementwise)

from immutable import ImmutableMatrix

from sparse import SparseMatrix

from expressions import (BlockDiagMatrix, BlockMatrix, FunctionMatrix,
     Identity, Inverse, MatAdd, MatMul, MatPow, MatrixExpr, MatrixSymbol,
     Trace, Transpose, ZeroMatrix, block_collapse, linear_factors,
     matrix_symbols, matrixify)
