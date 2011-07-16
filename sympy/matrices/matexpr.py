from sympy import Expr, Symbol, Eq, Mul, Add
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
    is_Transpose = False
    is_ZeroMatrix = False
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

    def eval_transpose(self):
        raise NotImplementedError()

    def eval_inverse(self):
        raise NotImplementedError()


class MatrixSymbol(MatrixExpr, Symbol):
    is_commutative = False
    def __new__(cls, name, n, m):
        obj = Basic.__new__(cls, name, n, m)
        return obj

    def _hashable_content(self):
        return(self.name, self.shape)

    def _check_shape(self):
        return True

    @property
    def shape(self):
        return self.args[1:3]

    @property
    def name(self):
        return self.args[0]

class Identity(MatrixSymbol):
    is_Identity = True
    def __new__(cls, n):
        return MatrixSymbol.__new__(cls, "I", n, n)

class ZeroMatrix(MatrixSymbol):
    is_ZeroMatrix = True
    def __new__(cls, n, m):
        return MatrixSymbol.__new__(cls, "0", n, m)

def matrix_symbols(expr):
    return [sym for sym in expr.free_symbols if sym.is_Matrix]

def matrixify(expr):
    """
    Mul.flatten is useful but always returns Muls. This is a quick fix to ensure
    that at least the outer class is a Matrix Expression.
    I expect that this will break for complex expression trees.

    Internal use only
    """
    if len(matrix_symbols(expr))==0: # No matrix symbols present
        return expr

    class_dict = {Mul:MatMul, Add:MatAdd, MatMul:MatMul, MatAdd:MatAdd}

    if expr.__class__ not in class_dict.keys():
        return expr

    args = map(matrixify, expr.args) # recursively call down the tree

    return Basic.__new__(class_dict[expr.__class__], *args)

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
