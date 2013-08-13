from __future__ import print_function, division

from sympy.core import Add, Basic, sympify
from sympy.functions import adjoint
from sympy.matrices.expressions.transpose import transpose
from sympy.strategies import (rm_id, unpack, flatten, sort, condition, debug,
        exhaust, do_one, glom)
from sympy.matrices.expressions.matexpr import MatrixExpr, ShapeError, ZeroMatrix
from sympy.utilities import default_sort_key

class MatAdd(MatrixExpr):
    """A Sum of Matrix Expressions

    MatAdd inherits from and operates like SymPy Add

    >>> from sympy import MatAdd, MatrixSymbol
    >>> A = MatrixSymbol('A', 5, 5)
    >>> B = MatrixSymbol('B', 5, 5)
    >>> C = MatrixSymbol('C', 5, 5)
    >>> MatAdd(A, B, C)
    A + B + C
    """
    is_MatAdd = True

    def __new__(cls, *args, **kwargs):
        check = kwargs.get('check', True)

        obj = Basic.__new__(cls, *args)
        if check:
            validate(*args)
        return obj

    @property
    def shape(self):
        return self.args[0].shape

    def _entry(self, i, j):
        return Add(*[arg._entry(i, j) for arg in self.args])

    def _eval_transpose(self):
        return MatAdd(*[transpose(arg) for arg in self.args]).doit()

    def _eval_adjoint(self):
        return MatAdd(*[adjoint(arg) for arg in self.args]).doit()

    def _eval_trace(self):
        from trace import Trace
        return MatAdd(*[Trace(arg) for arg in self.args]).doit()

    def doit(self, **ignored):
        return canonicalize(self)

def validate(*args):
    if not all(arg.is_Matrix for arg in args):
        raise TypeError("Mix of Matrix and Scalar symbols")
    A = args[0]
    for B in args[1:]:
        if A.shape != B.shape:
            raise ShapeError("Matrices %s and %s are not aligned"%(A, B))

factor_of = lambda arg: arg.as_coeff_mmul()[0]
matrix_of = lambda arg: unpack(arg.as_coeff_mmul()[1])
def combine(cnt, mat):
    from .matmul import MatMul
    if cnt == 1:
        return mat
    else:
        return cnt * mat

rules = (rm_id(lambda x: x == 0 or isinstance(x, ZeroMatrix)),
         unpack,
         flatten,
         glom(matrix_of, factor_of, combine),
         sort(default_sort_key))

canonicalize = exhaust(condition(lambda x: isinstance(x, MatAdd),
                                 do_one(*rules)))
