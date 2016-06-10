from __future__ import print_function, division

from sympy.core import Basic
from sympy.matrices.expressions.transpose import transpose
from sympy.matrices.expressions.matexpr import MatrixExpr


class DotProduct(MatrixExpr):
    """
        Dot Product of vector matrices
    """

    def __new__(cls, arg1, arg2):
        if not arg1.is_Matrix:
            raise TypeError("Input to Dot Product, %s, not a matrix" % str(arg1))
        if not arg2.is_Matrix:
            raise TypeError("Input to Dot Product, %s, not a matrix" % str(arg2))
        if not (1 in arg1.shape):
            raise TypeError("Input to Dot Product, %s, not a vector" % str(arg1))
        if not (1 in arg2.shape):
            raise TypeError("Input to Dot Product, %s, not a vector" % str(arg1))

        if arg1.shape != arg2.shape:
            raise TypeError("Input to Dot Product, %s and %s, are not of same dimensions" % (str(arg1), str(arg2)))

        return Basic.__new__(cls, arg1, arg2)

    def doit(self, expand=False):
        try:
            if self.args[0].shape[0] == 1:
                return (self.args[0]*transpose(self.args[1])).doit()[0]
            else:
                return (transpose(self.args[0])*self.args[1]).doit()[0]
        except (AttributeError, NotImplementedError):
            return self
