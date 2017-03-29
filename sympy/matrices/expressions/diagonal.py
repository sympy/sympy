from __future__ import print_function, division

from sympy.matrices.expressions import MatrixExpr
from sympy.core import S
from sympy.functions import KroneckerDelta

class DiagonalMatrix(MatrixExpr):
    arg = property(lambda self: self.args[0])
    shape = property(lambda self: (self.arg.shape[0], self.arg.shape[0]))

    def _entry(self, i, j):
        return self.arg[i, 0]*KroneckerDelta(i,j)

class DiagonalOf(MatrixExpr):
    arg = property(lambda self: self.args[0])
    shape = property(lambda self: (self.arg.shape[0], S.One))

    def _entry(self, i, j):
        return self.arg[i, i]
