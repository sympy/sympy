""" A module which handles Matrix Expressions """

from blockmatrix import BlockMatrix, BlockDiagMatrix, block_collapse
from inverse import Inverse
from matadd import MatAdd
from matexpr import (Identity, MatrixExpr, MatrixSymbol, ZeroMatrix,
     linear_factors, matrix_symbols, matrixify)
from matmul import MatMul
from matpow import MatPow
from transpose import Transpose
