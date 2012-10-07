from matexpr import MatrixExpr, ShapeError, ZeroMatrix
from sympy import Add, S, Basic

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
        simplify = kwargs.get('simplify', True)
        check    = kwargs.get('check'   , True)

        # TODO: This is a kludge
        # We still use Matrix + 0 in a few places. This removes it
        # In particular see matrix_multiply
        args = filter(lambda x: not x == 0, args)

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

    def canonicalize(self):
        return canonicalize(self)

def validate(*args):
    if not all(arg.is_Matrix for arg in args):
        raise TypeError("Mix of Matrix and Scalar symbols")
    A = args[0]
    for B in args[1:]:
        if A.shape != B.shape:
            raise ShapeError("Matrices %s and %s are not aligned"%(A,B))

from sympy.rr import rmid, unpack, flatten, sort, canon, condition, glom, debug

def newadd(*args):
    return Basic.__new__(MatAdd, *args)

def condition_matadd(rule):
    is_matadd = lambda x: x.is_Matrix and x.is_MatAdd
    return condition(is_matadd, rule)

def glom_MatAdd(expr):
    def counts(arg):
        if arg.is_MatMul:
            return arg.as_coeff_mmul()
        return 1, arg
    freqs = {}
    for arg in expr.args:
        count, m = counts(arg)
        freqs[m] = freqs.get(m, 0) + count

    # If it is the same expr then return the old one
    if all(v==1 for v in freqs.values()):
        return expr
    args = [m if c == 1 else m*c for m, c in freqs.items()]
    if set(args) == set(expr.args):
        return expr

    return Basic.__new__(MatAdd, *args)


rules = (rmid(lambda x: x == 0 or x.is_Matrix and x.is_ZeroMatrix),
         unpack,
         flatten,
         glom_MatAdd,
         sort(str))

canonicalize = canon(*map(condition_matadd, rules))

from matmul import MatMul
