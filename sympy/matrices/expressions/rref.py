from .matexpr import MatrixExpr
from sympy.matrices.matrices import MatrixBase

class ReducedRowEchelonForm(MatrixExpr):

    def __new__(cls, arg, **kwargs):
        return MatrixExpr.__new__(cls, arg, **kwargs)

    def doit(self, deep=True):
        arg = self.args[0]
        if deep == True:
            arg = arg.doit()

        if isinstance(arg, MatrixBase):
            return arg.rref()
        return arg._eval_rref()

    @property
    def rows(self):
        return self.args[0].rows

    @property
    def cols(self):
        return self.args[0].cols
