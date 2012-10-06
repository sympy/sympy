from matexpr import MatrixExpr, ShapeError, matrixify, ZeroMatrix
from sympy import Add, S, Basic


class MatAdd(MatrixExpr, Add):
    """A Sum of Matrix Expressions

    MatAdd inherits from and operates like SymPy Add

    >>> from sympy import MatAdd, MatrixSymbol
    >>> A = MatrixSymbol('A', 5, 5)
    >>> B = MatrixSymbol('B', 5, 5)
    >>> C = MatrixSymbol('C', 5, 5)
    >>> MatAdd(A, B, C)
    A + B + C
    """

    def __new__(cls, *args, **kwargs):
        simplify = kwargs.get('simplify', True)
        check    = kwargs.get('check'   , True)

        obj = Basic.__new__(cls, *args)
        if check:
            validate(*args)
        if simplify:
            return canonicalize(obj)
        else:
            return obj

    @property
    def shape(self):
        return self.args[0].shape

    def _entry(self, i, j):
        return Add(*[arg._entry(i, j) for arg in self.args])

    def _eval_transpose(self):
        from transpose import Transpose
        return MatAdd(*[Transpose(arg) for arg in self.args])

    def _eval_trace(self):
        from trace import Trace
        return MatAdd(*[Trace(arg) for arg in self.args])

def validate(*args):
    if not all(arg.is_Matrix for arg in args):
        raise ValueError("Mix of Matrix and Scalar symbols")
    A = args[0]
    for B in args[1:]:
        if A.shape != B.shape:
            raise ShapeError("Matrices %s and %s are not aligned"%(A,B))

from sympy.rr import rmid, unpack, flatten, sort, canon, condition, glom

def newadd(*args):
    return Basic.__new__(MatAdd, *args)

def condition_matadd(rule):
    is_matadd = lambda x: x.is_Matrix and x.is_Add
    return condition(is_matadd, rule)

rules = (rmid(lambda x: x == 0 or x.is_Matrix and x.is_ZeroMatrix),
         unpack,
         flatten,
         sort(str),
         glom(lambda num, arg: num*arg))

canonicalize = canon(*map(condition_matadd, rules))

from matmul import MatMul
