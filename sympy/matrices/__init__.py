"""A module that handles matrices.

Includes functions for fast creating matrices like zero, one/eye, random matrix etc.
"""
from matrices import (Matrix, SparseMatrix, zeros, ones, eye, diag,
     hessian, randMatrix, GramSchmidt, wronskian, casoratian,
     list2numpy, matrix2numpy, DeferredVector, block_diag, symarray, ShapeError,
     NonSquareMatrixError)
