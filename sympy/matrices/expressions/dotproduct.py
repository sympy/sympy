from __future__ import print_function, division

from sympy.core import Basic, sympify
from sympy.matrices.expressions.matexpr import (MatrixExpr)


class DotProduct(MatrixExpr):
    """
        Dot Product of vector matrices
    """

    def __new__(cls, *args):
        args = list(map(sympify, args))
        if validate(*args):
            return Basic.__new__(cls, *args)

    @property
    def args(self):
        return self.args

    def doit(self, expand=False):
        try:
            return dot_product(self.args)
        except (AttributeError, NotImplementedError):
            return self


def dot_product(*args):
    """
    Calculates dot product of the arguments.
    :param args: List of arguments for dot product.
    :return: dot product of the arguments provided.
    """
    return True


def validate(*args):
    """
    Validates the matrices for dot product.
    :param args: List of matrices to be validated for being a vector
    :return: True on being valid, else raises the error.
    """
    if len(args) != 2:
        raise
        # Only two vectors could be dot producted
    for i in args:
        if not i.is_Matrix:
            raise
            # Not Matrix
        if not (i.shape[0] == 1 or i.shape[1] == 1):
            raise
            # Not a Vector
    if args[0].shape != args[1].shape:
        raise
        # Vectors shape are not same
    return True
