from sympy import Symbol, sympify, Tuple, Integer
from sympy.core.basic import Basic
from sympy.core.decorators import _sympifyit, call_highest_priority
from sympy.core.singleton import S
from sympy.matrices import ShapeError
from sympy.simplify import simplify


class MatrixExpr(Basic):
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
    is_MatrixExpr = True
    is_Identity = None
    is_Inverse = False
    is_Transpose = False
    is_ZeroMatrix = False
    is_BlockMatrix = False
    is_MatAdd = False
    is_MatMul = False

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
    def rows(self):
        return self.shape[0]

    @property
    def cols(self):
        return self.shape[1]

    @property
    def is_square(self):
        return self.rows == self.cols

    def _eval_transpose(self):
        raise NotImplementedError()

    def _eval_inverse(self):
        raise NotImplementedError()

    def _eval_power(self, exp):
        return MatPow(self, exp)

    def _eval_simplify(self, **kwargs):
        if self.is_Atom:
            return self
        else:
            return self.__class__(*[simplify(x, **kwargs) for x in self.args])

    def _eval_adjoint(self):
        return self.T.conjugate()

    def _entry(self, i, j):
        raise NotImplementedError(
            "Indexing not implemented for %s" % self.__class__.__name__)

    def adjoint(self):
        raise NotImplementedError(
            "adjoint not implemented for %s" % self.__class__.__name__)

    def conjugate(self):
        raise NotImplementedError(
            "conjugate not implemented for %s" % self.__class__.__name__)

    def transpose(self):
        try:
            return self._eval_transpose()
        except (AttributeError, NotImplementedError):
            return Basic.__new__(Transpose, self)

    T = property(transpose, None, None, 'Matrix transposition.')

    @property
    def I(self):
        return Inverse(self)

    def valid_index(self, i, j):
        def is_valid(idx):
            return isinstance(idx, (int, Integer, Symbol))
        return (is_valid(i) and is_valid(j) and
                0 <= i < self.rows and 0 <= j < self.cols)

    def __getitem__(self, key):
        if isinstance(key, tuple) and len(key) == 2:
            i, j = key
            if isinstance(i, slice) or isinstance(j, slice):
                raise NotImplementedError(
                    "Slicing is not implemented for %s" % self.__class__.__name__)
            i, j = sympify(i), sympify(j)
            if self.valid_index(i, j) is not False:
                return self._entry(i, j)
            else:
                raise IndexError("Invalid indices (%s, %s)" % (i, j))
        raise IndexError("Invalid index, wanted %s[i,j]" % self)

    def as_explicit(self):
        """
        Returns a dense Matrix with elements represented explicitly

        Returns an object of type ImmutableMatrix.

        Examples
        ========

        >>> from sympy import Identity
        >>> I = Identity(3)
        >>> I
        I
        >>> I.as_explicit()
        [1, 0, 0]
        [0, 1, 0]
        [0, 0, 1]

        See Also
        ========
        as_mutable: returns mutable Matrix type

        """
        from sympy.matrices.immutable import ImmutableMatrix
        return ImmutableMatrix([[    self[i, j]
                            for j in range(self.cols)]
                            for i in range(self.rows)])

    def as_mutable(self):
        """
        Returns a dense, mutable matrix with elements represented explicitly

        Examples
        ========

        >>> from sympy import Identity
        >>> I = Identity(3)
        >>> I
        I
        >>> I.shape
        (3, 3)
        >>> I.as_mutable()
        [1, 0, 0]
        [0, 1, 0]
        [0, 0, 1]

        See Also
        ========
        as_explicit: returns ImmutableMatrix
        """
        return self.as_explicit().as_mutable()

    def __array__(self):
        from numpy import empty
        a = empty(self.shape, dtype=object)
        for i in range(self.rows):
            for j in range(self.cols):
                a[i, j] = self[i, j]
        return a

    def equals(self, other):
        """
        Test elementwise equality between matrices, potentially of different
        types

        >>> from sympy import Identity, eye
        >>> Identity(3).equals(eye(3))
        True
        """
        return self.as_explicit().equals(other)

    def canonicalize(self):
        return self

    def as_coeff_mmul(self):
        return 1, Basic.__new__(MatMul, self)


class MatrixSymbol(MatrixExpr):
    """Symbolic representation of a Matrix object

    Creates a SymPy Symbol to represent a Matrix. This matrix has a shape and
    can be included in Matrix Expressions

    >>> from sympy import MatrixSymbol, Identity
    >>> A = MatrixSymbol('A', 3, 4) # A 3 by 4 Matrix
    >>> B = MatrixSymbol('B', 4, 3) # A 4 by 3 Matrix
    >>> A.shape
    (3, 4)
    >>> 2*A*B + Identity(3)
    2*A*B + I
    """
    is_commutative = False
    is_Atom = True

    def __new__(cls, name, n, m):
        n, m = sympify(n), sympify(m)
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
        # only do substitutions in shape
        shape = Tuple(*self.shape)._subs(old, new)
        return MatrixSymbol(self.name, *shape)

    def __call__(self, *args):
        raise TypeError( "%s object is not callable" % self.__class__ )

    def _entry(self, i, j):
        # MatMul _entry will pass us a Dummy and ask that we remember it
        # so that it can be summed over later. We'll use the function syntax
        if i.is_Dummy or j.is_Dummy:
            return Symbol(self.name)(i, j)
        # If that isn't the case we'd really rather just make a symbol
        # They are simpler and look much nicer
        else:
            return Symbol('%s_%s%s' % (self.name, str(i), str(j)))

    @property
    def free_symbols(self):
        return set((self,))


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

    def _eval_transpose(self):
        return self

    def _eval_trace(self):
        return self.rows

    def _eval_inverse(self):
        return self

    def conjugate(self):
        return self

    def _entry(self, i, j):
        if i == j:
            return S.One
        else:
            return S.Zero


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

    def _eval_transpose(self):
        return ZeroMatrix(self.cols, self.rows)

    def _eval_trace(self):
        return S.Zero

    def conjugate(self):
        return self

    def _entry(self, i, j):
        return S.Zero


def matrix_symbols(expr):
    return [sym for sym in expr.free_symbols if sym.is_Matrix]

from matmul import MatMul
from matadd import MatAdd
from matpow import MatPow
from transpose import Transpose
from inverse import Inverse
