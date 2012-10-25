from matexpr import MatrixExpr, ShapeError
from sympy import Mul, Basic
from sympy.rules import (unpack, flatten, sort, condition, exhaust, do_one)

class HadamardProduct(MatrixExpr):
    """Elementwise Product of Matrix Expressions

    >>> from sympy import HadamardProduct, MatrixSymbol
    >>> A = MatrixSymbol('A', 5, 5)
    >>> B = MatrixSymbol('B', 5, 5)
    >>> C = MatrixSymbol('C', 5, 5)
    >>> HadamardProduct(A, B, C)
    A.*B.*C
    """
    is_HadamardProduct = True

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', True)
        check    = kwargs.get('check'   , True)

        obj = Basic.__new__(cls, *args)
        if check:
            validate(*args)
        if evaluate:
            return canonicalize(obj)
        else:
            return obj

    @property
    def shape(self):
        return self.args[0].shape

    def _entry(self, i, j):
        return Mul(*[arg._entry(i, j) for arg in self.args])

    def _eval_transpose(self):
        from transpose import Transpose
        return HadamardProduct(*[Transpose(arg) for arg in self.args])

    def canonicalize(self):
        return canonicalize(self)

def validate(*args):
    if not all(arg.is_Matrix for arg in args):
        raise TypeError("Mix of Matrix and Scalar symbols")
    A = args[0]
    for B in args[1:]:
        if A.shape != B.shape:
            raise ShapeError("Matrices %s and %s are not aligned"%(A,B))

rules = (unpack,
         flatten,
         sort(str))

canonicalize = exhaust(condition(lambda x: isinstance(x, HadamardProduct),
                                 do_one(*rules)))
