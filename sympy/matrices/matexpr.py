from sympy import Expr, Symbol, Eq
from sympy.core.basic import Basic
from sympy.core.singleton import S
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.cache import cacheit
from sympy.core.compatibility import any, all, reduce
from sympy.matrices import ShapeError

class MatrixExpr(Expr):

    _op_priority = 11.0

    is_Matrix = True
    is_Identity = False
    is_Inverse = False
    is_Zero = False
    is_BlockMatrix = False

    is_commutative = False

    # The following is adapted from the core Expr object

    def __neg__(self):
        return MatMul(S.NegativeOne, self)
    def __abs__(self):
        raise NotImplementedError

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return MatAdd(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return MatAdd(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return MatAdd(self, -other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return MatAdd(other, -self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return MatMul(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return MatMul(other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        if other == -S.One:
            return Inverse(self)
        return MatPow(self, other)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("Matrix Power not defined")
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        raise NotImplementedError()
        return MatMul(self, Pow(other, S.NegativeOne))
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise NotImplementedError()
        return MatMul(other, Pow(self, S.NegativeOne))

    __truediv__ = __div__
    __rtruediv__ = __rdiv__


    @property
    def n(self):
        return self.shape[0]
    @property
    def m(self):
        return self.shape[1]

    @property
    def is_square(self):
        return Eq(self.n, self.m)

    def subs(self, *args):
        ex = matrixify(Basic.subs(self, *args))
        #if ex.is_Matrix and not ex._check_shape():
#            raise ShapeError("Substituted expression invalid")
        return ex


class MatrixSymbol(MatrixExpr, Symbol):
    def __new__(cls, name, n, m, **assumptions):
        obj = Symbol.__new__(cls, name, False, **assumptions)
        obj.shape = (n,m)
        return obj

    def _hashable_content(self):
        return(self.name, self.shape)

    def _check_shape(self):
        return True

class Identity(MatrixSymbol):
    is_Identity = True
    def __new__(cls, n):
        return MatrixSymbol.__new__(cls, "I", n, n)

class ZeroMatrix(MatrixSymbol):
    is_Zero = True
    def __new__(cls, n, m):
        return MatrixSymbol.__new__(cls, "0", n, m)

def matrixify(expr):
    """
    Mul.flatten is useful but always returns Muls. This is a quick fix to ensure
    that at least the outer class is a Matrix Expression.
    I expect that this will break for complex expression trees.

    Internal use only
    """
    if not (expr.is_Mul or expr.is_Add):
        return expr
    while not expr.is_Matrix and any(arg.is_Matrix for arg in expr.args):
        if expr.is_Mul:
            expr = MatMul(*expr.args)
        if expr.is_Add:
            expr = MatAdd(*expr.args)
    return expr

def matsimp(expr):
    if all(mat.is_Identity for mat in new):
        new = [new[0]]
    else:
        new = [mat for mat in new if not mat.is_Identity] # clear ident


from matmul import MatMul
from matadd import MatAdd
from matpow import MatPow
from transpose import Transpose
from inverse import Inverse
