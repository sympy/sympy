from __future__ import print_function, division

from sympy.matrices.expressions import MatrixExpr
from sympy.core import S, Eq
from sympy.functions.special.tensor_functions import KroneckerDelta
class DiagonalMatrix(MatrixExpr):
    arg = property(lambda self: self.args[0])
    shape = property(lambda self: (self.arg.shape[0], self.arg.shape[0]))

    def _entry(self, i, j):
        eq = Eq(i, j)
        if eq is S.false:
            return S.Zero
        elif eq is S.true:
            return self.arg[i, i]
        return self.arg[i, j]*KroneckerDelta(i, j)

class DiagonalOf(MatrixExpr):
    arg = property(lambda self: self.args[0])
    shape = property(lambda self: (self.arg.shape[0], S.One))

    def _entry(self, i, j):
        return self.arg[i, i]
