""" A module which handles Matrix Expressions """

from blockmatrix import BlockMatrix, BlockDiagMatrix, block_collapse
from funcmatrix import FunctionMatrix
from inverse import Inverse
from matadd import MatAdd
from matexpr import (Identity, MatrixExpr, MatrixSymbol, ZeroMatrix,
     matrix_symbols)
from matmul import MatMul
from matpow import MatPow
from trace import Trace
from transpose import Transpose
from adjoint import Adjoint
from hadamard import hadamard_product, HadamardProduct
