from __future__ import print_function, division

from sympy.core import Basic
from sympy.core.sympify import _sympify
from sympy.matrices.expressions.transpose import transpose
from sympy.matrices.expressions.matexpr import MatrixExpr


class DotProduct(MatrixExpr):
    """
        Dot Product of vector matrices
    """

    def __new__(cls, arg1, arg2):
        arg1, arg2 = _sympify((arg1, arg2))

        if not arg1.is_Matrix:
            raise TypeError("Argument 1 of DotProduct is not a matrix")
        if not arg2.is_Matrix:
            raise TypeError("Argument 2 of DotProduct is not a matrix")
        if not (1 in arg1.shape):
            raise TypeError("Argument 1 of DotProduct is not a vector")
        if not (1 in arg2.shape):
            raise TypeError("Argument 2 of DotProduct is not a vector")

        if set(arg1.shape) != set(arg2.shape):
            raise TypeError("DotProduct arguments are not the same length" % (str(arg1), str(arg2)))

        return Basic.__new__(cls, arg1, arg2)

    def doit(self, expand=False):
        if self.args[0].shape == self.args[1].shape:
            if self.args[0].shape[0] == 1:
                mul = self.args[0]*transpose(self.args[1])
            else:
                mul = transpose(self.args[0])*self.args[1]
        else:
            if self.args[0].shape[0] == 1:
                mul = self.args[0]*self.args[1]
            else:
                mul = transpose(self.args[0])*transpose(self.args[1])

        return mul[0]
