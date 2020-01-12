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

Matrix = MutableMatrix = MutableDenseMatrix

from .sparse import MutableSparseMatrix
from .sparsetools import banded
from .immutable import ImmutableDenseMatrix, ImmutableSparseMatrix

ImmutableMatrix = ImmutableDenseMatrix
SparseMatrix = MutableSparseMatrix

from .expressions import (
    MatrixSlice, BlockDiagMatrix, BlockMatrix, FunctionMatrix, Identity,
    Inverse, MatAdd, MatMul, MatPow, MatrixExpr, MatrixSymbol, Trace,
    Transpose, ZeroMatrix, OneMatrix, blockcut, block_collapse, matrix_symbols, Adjoint,
    hadamard_product, HadamardProduct, HadamardPower, Determinant,
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

from .subspaces import (
    _columnspace, _nullspace, _rowspace, _orthogonalize,
    columnspace, nullspace, rowspace, orthogonalize)

from .eigen import (
    _eigenvals, _eigenvects, _is_diagonalizable, _diagonalize,
    _eval_is_positive_definite,
    _is_positive_definite, _is_positive_semidefinite,
    _is_negative_definite, _is_negative_semidefinite, _is_indefinite,
    _jordan_form, _left_eigenvects, _singular_values,
    eigenvals, eigenvects, is_diagonalizable, diagonalize,
    is_positive_definite, is_positive_semidefinite, is_negative_definite,
    is_negative_semidefinite, is_indefinite,
    jordan_form, left_eigenvects, singular_values)


__all__ = [
    'ShapeError', 'NonSquareMatrixError',

    'GramSchmidt', 'casoratian', 'diag', 'eye', 'hessian', 'jordan_cell',
    'list2numpy', 'matrix2numpy', 'matrix_multiply_elementwise', 'ones',
    'randMatrix', 'rot_axis1', 'rot_axis2', 'rot_axis3', 'symarray',
    'wronskian', 'zeros',

    'MutableDenseMatrix',

    'DeferredVector', 'MatrixBase',

    'Matrix', 'MutableMatrix',

    'MutableSparseMatrix',

    'banded',

    'ImmutableDenseMatrix', 'ImmutableSparseMatrix',

    'ImmutableMatrix', 'SparseMatrix',

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

    '_columnspace', '_nullspace', '_rowspace', '_orthogonalize',
    'columnspace', 'nullspace', 'rowspace', 'orthogonalize',

    '_eigenvals', '_eigenvects', '_is_diagonalizable', '_diagonalize',
    '_eval_is_positive_definite',
    '_is_positive_definite', '_is_positive_semidefinite',
    '_is_negative_definite', '_is_negative_semidefinite', '_is_indefinite',
    '_jordan_form', '_left_eigenvects', '_singular_values',
    'eigenvals', 'eigenvects', 'is_diagonalizable', 'diagonalize',
    'is_positive_definite', 'is_positive_semidefinite', 'is_negative_definite',
    'is_negative_semidefinite', 'is_indefinite',
    'jordan_form', 'left_eigenvects', 'singular_values',
]
