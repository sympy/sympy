from sympy import Expr, Symbol, Mul, Add, Pow, expand, sympify, Tuple
from sympy.core.basic import Basic
from sympy.core.singleton import S
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.matrices import ShapeError

class MatrixExpr(Expr):
    """ Matrix Expression Class
    Matrix Expressions subclass SymPy Expr's so that
    MatAdd inherits from Add
    MatMul inherits from Mul
    MatPow inherits from Pow

    They use _op_priority to gain control with binary operations (+, *, -, **)
    are used

    They implement operations specific to Matrix Algebra.
    """

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
        return MatMul(self, other**S.NegativeOne)
    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        raise NotImplementedError()
        #return MatMul(other, Pow(self, S.NegativeOne))

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
        return self.n == self.m

    def eval_transpose(self):
        raise NotImplementedError()

    def eval_inverse(self):
        raise NotImplementedError()

    @property
    def T(self):
        return Transpose(self)

    @property
    def I(self):
        return Inverse(self)


class MatrixSymbol(MatrixExpr, Symbol):
    """Symbolic representation of a Matrix object

    Creates a SymPy Symbol to represent a Matrix. This matrix has a shape and
    can be included in Matrix Expressions

    >>> from sympy import MatrixSymbol, Identity
    >>> A = MatrixSymbol('A', 3, 4) # A 3 by 4 Matrix
    >>> B = MatrixSymbol('B', 4, 3) # A 4 by 3 Matrix
    >>> A.shape
    (3, 4)
    >>> 2*A*B + Identity(3)
    I + 2*A*B
    """
    is_commutative = False

    def __new__(cls, name, n, m):
        obj = Basic.__new__(cls, name, n, m)
        return obj

    def _hashable_content(self):
        return(self.name, self.shape)

    @property
    def shape(self):
        return self.args[1:3]

    @property
    def name(self):
        return self.args[0]

    def _eval_subs(self, old, new):
        if self==old:
            return new
        else:
            shape = Tuple(*self.shape).subs(old, new)
            return MatrixSymbol(self.name, *shape)

    def __call__(self, *args):
        raise TypeError( "%s object is not callable"%self.__class__ )


class Identity(MatrixSymbol):
    """The Matrix Identity I - multiplicative identity
    >>> from sympy.matrices import Identity, MatrixSymbol
    >>> A = MatrixSymbol('A', 3, 5)
    >>> I = Identity(3)
    >>> I*A
    A
    """

    is_Identity = True
    def __new__(cls, n):
        return MatrixSymbol.__new__(cls, "I", n, n)

    def transpose(self):
        return self

class ZeroMatrix(MatrixSymbol):
    """The Matrix Zero 0 - additive identity
    >>> from sympy import MatrixSymbol, ZeroMatrix
    >>> A = MatrixSymbol('A', 3, 5)
    >>> Z = ZeroMatrix(3, 5)
    >>> A+Z
    A
    >>> Z*A.T
    0
    """
    is_ZeroMatrix = True
    def __new__(cls, n, m):
        return MatrixSymbol.__new__(cls, "0", n, m)
    def transpose(self):
        return ZeroMatrix(self.m, self.n)

def matrix_symbols(expr):
    return [sym for sym in expr.free_symbols if sym.is_Matrix]

def matrixify(expr):
    """
    Recursively walks down an expression tree changing Expr's to MatExpr's
    i.e. Add -> MatAdd
         Mul -> MatMul

    Only changes those Exprs which contain MatrixSymbols

    This function is useful when traditional SymPy functions which use Mul and
    Add are called on MatrixExpressions. Examples flatten, expand, simplify...

    Calling matrixify after calling these functions will reset classes back to
    their matrix equivalents

    For internal use
    """
    if len(matrix_symbols(expr))==0: # No matrix symbols present
        return expr

    class_dict = {Mul:MatMul, Add:MatAdd, MatMul:MatMul, MatAdd:MatAdd,
            Pow:MatPow, MatPow:MatPow}

    if expr.__class__ not in class_dict.keys():
        return expr

    args = map(matrixify, expr.args) # Recursively call down the tree

    return Basic.__new__(class_dict[expr.__class__], *args)

def linear_factors(expr, *syms):
    """Reduce a Matrix Expression to a sum of linear factors

    Given symbols and a matrix expression linear in those symbols return a
    dict mapping symbol to the linear factor

    >>> from sympy import MatrixSymbol, linear_factors, symbols
    >>> n, m, l = symbols('n m l')
    >>> A = MatrixSymbol('A', n, m)
    >>> B = MatrixSymbol('B', m, l)
    >>> C = MatrixSymbol('C', n, l)
    >>> linear_factors(2*A*B + C, B, C)
    {B: 2*A, C: I}
    """

    expr = matrixify(expand(expr))
    d = {}
    if expr.is_Matrix and expr.is_Symbol:
        if expr in syms:
            d[expr] = Identity(expr.n)

    if expr.is_Add:
        for sym in syms:
            total_factor = 0
            for arg in expr.args:
                factor = arg.coeff(sym)
                if not factor:
                    # .coeff fails when powers are in the expression
                    if sym in arg.free_symbols:
                        raise ValueError("Expression not linear in symbols")
                    else:
                        factor = 0
                factor = sympify(factor)
                if not factor.is_Matrix:
                    if factor.is_zero:
                        factor = ZeroMatrix(expr.n, sym.n)
                        if not sym.m == expr.m:
                            raise ShapeError(
                            "%s not compatible as factor of %s"%(sym, expr))
                    else:
                        factor = Identity(sym.n)*factor
                total_factor += factor
            d[sym] = total_factor
    elif expr.is_Mul:
        for sym in syms:
            factor = expr.coeff(sym)
            if not factor:
                # .coeff fails when powers are in the expression
                if sym in expr.free_symbols:
                    raise ValueError("Expression not linear in symbols")
                else:
                    factor = 0
            factor = sympify(factor)
            if not factor.is_Matrix:
                if factor.is_zero:
                    factor = ZeroMatrix(expr.n, sym.n)
                    if not sym.m == expr.m:
                        raise ShapeError("%s not compatible as factor of %s"%
                                (sym, expr))
                else:
                    factor = Identity(sym.n)*factor
            d[sym] = factor

    if any(sym in matrix_symbols(Tuple(*d.values())) for sym in syms):
        raise ValueError("Expression not linear in symbols")

    return d

from matmul import MatMul
from matadd import MatAdd
from matpow import MatPow
from transpose import Transpose
from inverse import Inverse
