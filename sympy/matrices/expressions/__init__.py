""" A module which handles Matrix Expressions """

from .slice import MatrixSlice
from .blockmatrix import BlockMatrix, BlockDiagMatrix, block_collapse, blockcut
from .funcmatrix import FunctionMatrix
from .inverse import Inverse
from .matadd import MatAdd
from .matexpr import (Identity, MatrixExpr, MatrixSymbol, ZeroMatrix, OneMatrix,
                      matrix_symbols)
from .matmul import MatMul
from .matpow import MatPow
from .trace import Trace, trace
from .determinant import Determinant, det
from .transpose import Transpose
from .adjoint import Adjoint
from .hadamard import hadamard_product, HadamardProduct, hadamard_power, HadamardPower
from .diagonal import DiagonalMatrix, DiagonalOf, DiagMatrix, diagonalize_vector
from .dotproduct import DotProduct
from .kronecker import kronecker_product, KroneckerProduct, combine_kronecker
