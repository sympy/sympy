"""A module that handles matrices.

Includes functions for fast creating matrices like zero, one/eye, random
matrix, etc.
"""
from .common import ShapeError, NonSquareMatrixError
from .dense import (
    GramSchmidt, casoratian, diag, eye, hessian, jordan_cell,
    list2numpy, matrix2numpy, matrix_multiply_elementwise, ones,
    randMatrix, rot_axis1, rot_axis2, rot_axis3, symarray, wronskian,
    zeros)
from .dense import MutableDenseMatrix
from .matrices import DeferredVector, MatrixBase

Matrix = MutableDenseMatrix

from .sparse import MutableSparseMatrix
from .sparsetools import banded
from .immutable import ImmutableDenseMatrix, ImmutableSparseMatrix
from .numpy import NumPyMatrix

ImmutableMatrix = ImmutableDenseMatrix
MutableMatrix   = MutableDenseMatrix
SparseMatrix    = MutableSparseMatrix
MMatrix         = MutableDenseMatrix
IMatrix         = ImmutableDenseMatrix

from .expressions import (
    MatrixSlice, BlockDiagMatrix, BlockMatrix, FunctionMatrix, Identity,
    Inverse, MatAdd, MatMul, MatPow, MatrixExpr, MatrixSymbol, Trace,
    Transpose, ZeroMatrix, OneMatrix, blockcut, block_collapse, matrix_symbols, Adjoint,
    hadamard_product, HadamardProduct, HadamardPower, Determinant, det,
    diagonalize_vector, DiagMatrix, DiagonalMatrix, DiagonalOf, trace,
    DotProduct, kronecker_product, KroneckerProduct,
    PermutationMatrix, MatrixPermute)

from .determinant import (
    _find_reasonable_pivot, _find_reasonable_pivot_naive,
    _berkowitz_toeplitz_matrix,
    _adjugate, _charpoly, _cofactor, _cofactor_matrix,
    _det, _det_bareiss, _det_berkowitz, _det_LU, _minor, _minor_submatrix,
    adjugate, charpoly, cofactor, cofactor_matrix,
    det, det_bareiss, det_berkowitz, det_LU, minor, minor_submatrix)

from .reductions import (
    _row_reduce_list, _row_reduce,
    _is_echelon, _echelon_form, _rank, _rref,
    is_echelon, echelon_form, rank, rref)


__all__ = [
    'ShapeError', 'NonSquareMatrixError',

    'GramSchmidt', 'casoratian', 'diag', 'eye', 'hessian', 'jordan_cell',
    'list2numpy', 'matrix2numpy', 'matrix_multiply_elementwise', 'ones',
    'randMatrix', 'rot_axis1', 'rot_axis2', 'rot_axis3', 'symarray',
    'wronskian', 'zeros',

    'DeferredVector', 'MatrixBase',

    'banded',

    'Matrix',
    'IMatrix',
    'MMatrix',

    'ImmutableMatrix',
    'MutableMatrix',
    'ImmutableDenseMatrix',
    'MutableDenseMatrix',
    'SparseMatrix',
    'ImmutableSparseMatrix',
    'MutableSparseMatrix',
    'NumPyMatrix',

    'MatrixSlice', 'BlockDiagMatrix', 'BlockMatrix', 'FunctionMatrix',
    'Identity', 'Inverse', 'MatAdd', 'MatMul', 'MatPow', 'MatrixExpr',
    'MatrixSymbol', 'Trace', 'Transpose', 'ZeroMatrix', 'OneMatrix',
    'blockcut', 'block_collapse', 'matrix_symbols', 'Adjoint',
    'hadamard_product', 'HadamardProduct', 'HadamardPower', 'Determinant',
    'det', 'diagonalize_vector', 'DiagMatrix', 'DiagonalMatrix',
    'DiagonalOf', 'trace', 'DotProduct', 'kronecker_product',
    'KroneckerProduct', 'PermutationMatrix', 'MatrixPermute',

    '_find_reasonable_pivot', '_find_reasonable_pivot_naive',
    '_berkowitz_toeplitz_matrix',
    '_adjugate', '_charpoly', '_cofactor', '_cofactor_matrix', '_det',
    '_det_bareiss', '_det_berkowitz', '_det_LU', '_minor', '_minor_submatrix',
    'adjugate', 'charpoly', 'cofactor', 'cofactor_matrix', 'det', 'det_bareiss',
    'det_berkowitz', 'det_LU', 'minor', 'minor_submatrix',

    '_row_reduce_list', '_row_reduce',
    '_is_echelon', '_echelon_form', '_rank', '_rref',
    'is_echelon', 'echelon_form', 'rank', 'rref',
]
