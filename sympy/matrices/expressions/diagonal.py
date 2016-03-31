from __future__ import print_function, division

from sympy.matrices.expressions import MatrixExpr
from sympy.core import S

class DiagonalMatrix(MatrixExpr):
    arg = property(lambda self: self.args[0])
    shape = property(lambda self: (self.arg.shape[0], self.arg.shape[0]))

    def _entry(self, i, j):
        return S.Zero if i != j else self.arg[i, 0]

class DiagonalOf(MatrixExpr):
    arg = property(lambda self: self.args[0])
    shape = property(lambda self: (self.arg.shape[0], S.One))

    def _entry(self, i, j):
        return self.arg[i, i]
