from __future__ import print_function, division

from sympy.matrices.expressions import MatrixExpr
from sympy.core import S

class DiagonalMatrix(MatrixExpr):
    arg = property(lambda self: self.args[0])
    shape = property(lambda self: (self.arg.shape[0], self.arg.shape[1]))

    def _entry(self, i, j):
        return S.Zero if i != j else self.arg[i, i]

class DiagonalOf(MatrixExpr):
    arg = property(lambda self: self.args[0])
    shape = property(lambda self: (self.arg.shape[0], self.arg.shape[1]))

    def _entry(self, i, j):
        return S.Zero if i != j else self.arg[i, i]
