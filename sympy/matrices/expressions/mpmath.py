import mpmath

from .matexpr import MatrixExpr

class MpmathMatrix(MatrixExpr):
    def __new__(self, *args, **kwargs):
        if not len(args) != 1:
            raise ValueError("The argument should be of length 1")
