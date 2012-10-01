import random

from sympy.core.basic import Basic
from sympy.core.compatibility import is_sequence
from sympy.core.expr import Expr
from sympy.core.function import count_ops
from sympy.core.decorators import call_highest_priority
from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.core.sympify import sympify
from sympy.functions.elementary.trigonometric import cos, sin
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.simplify import simplify as _simplify
from sympy.utilities.exceptions import SymPyDeprecationWarning

from sympy.matrices.matrices import MatrixBase, ShapeError, a2idx, classof

# uncomment the import of as_int and delete the function when merged with 0.7.2
#from sympy.core.compatibility import as_int

def as_int(i):
    ii = int(i)
    if i != ii:
        raise TypeError()
    return ii

def _iszero(x):
    """Returns True if x is zero."""
    return x.is_zero

class Matrix(MatrixBase):

    is_MatrixExpr = False

    _op_priority = 12.0
    _class_priority = 10

    @classmethod
    def _new(cls, *args, **kwargs):
        rows, cols, flat_list = MatrixBase._handle_creation_inputs(*args, **kwargs)
        self = object.__new__(cls)
        self.rows = rows
        self.cols = cols
        self._mat = list(flat_list) # create a shallow copy
        return self

    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    @property
    def is_Identity(self):
        if not self.is_square:
            return False
        if not all(self[i, i] == 1 for i in range(self.rows)):
            return False
        for i in range(self.rows):
            for j in range(i + 1, self.cols):
                if self[i, j] or self[j, i]:
                    return False
        return True

    def __getitem__(self, key):
        """
        >>> from sympy import Matrix, I
        >>> m = Matrix([
        ... [1, 2 + I],
        ... [3, 4    ]])
        >>> m[1, 0]
        3

        """
        if type(key) is tuple:
            i, j = key
            if type(i) is slice or type(j) is slice:
                return self.submatrix(key)
            else:
                i, j = self.key2ij(key)
                return self._mat[i*self.cols + j]
        else:
            # row-wise decomposition of matrix
            if type(key) is slice:
                return self._mat[key]
            return self._mat[a2idx(key)]

    def __setitem__(self, key, value):
        """
        >>> from sympy import Matrix, I
        >>> m = Matrix(((1, 2+I), (3, 4)))
        >>> m
        [1, 2 + I]
        [3,     4]
        >>> m[1, 0] = 9
        >>> m
        [1, 2 + I]
        [9,     4]
        >>> m[1, 0] = [[0, 1]]
        >>> m
        [1, 2 + I]
        [0,     1]

        """
        from dense import Matrix

        if type(key) is tuple:
            i, j = key
            if type(i) is slice or type(j) is slice:
                if isinstance(value, MatrixBase):
                    self.copyin_matrix(key, value)
                    return
                if not isinstance(value, Expr) and is_sequence(value):
                    self.copyin_list(key, value)
                    return
                raise ValueError('unexpected value: %s' % value)
            else:
                i, j = self.key2ij(key)
                is_mat = isinstance(value, MatrixBase)
                if not is_mat and \
                    not isinstance(value, Expr) and is_sequence(value):
                    value = Matrix(value)
                    is_mat = True
                if is_mat:
                    self.copyin_matrix((slice(i, i + value.rows),
                                        slice(j, j + value.cols)), value)
                else:
                    self._mat[i*self.cols + j] = sympify(value)
                return
        else:
            # row-wise decomposition of matrix
            if type(key) is slice:
                raise IndexError("Vector slices not implemented yet.")
            else:
                k = a2idx(key)
                self._mat[k] = sympify(value)
                return

    def tolist(self):
        """
        Return the Matrix as a nested Python list.

        Examples
        ========

        >>> from sympy import Matrix, ones
        >>> m = Matrix(3, 3, range(9))
        >>> m
        [0, 1, 2]
        [3, 4, 5]
        [6, 7, 8]
        >>> m.tolist()
        [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
        >>> ones(3, 0).tolist()
        [[], [], []]

        When there are no rows then it will not be possible to tell how
        many columns were in the original matrix:

        >>> ones(0, 3).tolist()
        []

        """
        if not self.rows:
            return []
        if not self.cols:
            return [[] for i in range(self.rows)]
        return [self._mat[i: i + self.cols]
            for i in range(0, len(self), self.cols)]

    def expand(self, **hints):
        """
        Expand each element of the matrix by calling ``expand()``.
        """
        out = self._new(self.rows, self.cols,
                map(lambda i: i.expand(**hints), self._mat))
        return out

    def subs(self, *args, **kwargs):
        """
        Create substituted expressions for each element with ``Expr.subs``.
        """
        out = self._new(self.rows, self.cols,
                map(lambda i: i.subs(*args, **kwargs), self._mat))
        return out

    def copyin_matrix(self, key, value):
        """
        Copy in values from a matrix into the given bounds.

        Parameters
        ==========

        key : slice
            The section of this matrix to replace.
        value : Matrix
            The matrix to copy values from.

        Examples
        ========

        >>> from sympy.matrices import Matrix, eye
        >>> M = Matrix([[0, 1], [2, 3], [4, 5]])
        >>> I = eye(3)
        >>> I[:3, :2] = M
        >>> I
        [0, 1, 0]
        [2, 3, 0]
        [4, 5, 1]
        >>> I[0, 1] = M
        >>> I
        [0, 0, 1]
        [2, 2, 3]
        [4, 4, 5]

        See Also
        ========

        copyin_list
        """
        rlo, rhi, clo, chi = self.key2bounds(key)
        shape = value.shape
        dr, dc = rhi - rlo, chi - clo
        if shape != (dr, dc):
            raise ShapeError("The Matrix `value` doesn't have the same dimensions " +
                "as the in sub-Matrix given by `key`.")

        for i in range(value.rows):
            for j in range(value.cols):
                self[i + rlo, j + clo] = value[i, j]

    def copyin_list(self, key, value):
        """
        Copy in elements from a list.

        Parameters
        ==========

        key : slice
            The section of this matrix to replace.
        value : iterable
            The iterable to copy values from.

        Examples
        ========

        >>> from sympy.matrices import eye
        >>> I = eye(3)
        >>> I[:2, 0] = [1, 2] # col
        >>> I
        [1, 0, 0]
        [2, 1, 0]
        [0, 0, 1]
        >>> I[1, :2] = [[3, 4]]
        >>> I
        [1, 0, 0]
        [3, 4, 0]
        [0, 0, 1]

        See Also
        ========

        copyin_matrix
        """
        if not is_sequence(value):
            raise TypeError("`value` must be an ordered iterable, not %s." % type(value))
        return self.copyin_matrix(key, Matrix(value))

    def row(self, i, f=None):
        """
        Elementary row selector (default) or operation using functor
        which is a function two args interpreted as (self[i, j], j).

        >>> from sympy import ones
        >>> I = ones(3)
        >>> I.row(1, lambda v, i: v*3)
        >>> I
        [1, 1, 1]
        [3, 3, 3]
        [1, 1, 1]
        >>> I.row(1)
        [3, 3, 3]

        See Also
        ========

        col
        row_swap
        row_del
        row_join
        row_insert
        delRowCol
        """
        if f is None:
            return self[i, :]
        for j in range(0, self.cols):
            self[i, j] = f(self[i, j], j)

    def col(self, j, f=None):
        """
        Elementary column selector (default) or operation using functor
        which is a function two args interpreted as (self[i, j], i).

        >>> from sympy import ones
        >>> I = ones(3)
        >>> I.col(0, lambda v, i: v*3)
        >>> I
        [3, 1, 1]
        [3, 1, 1]
        [3, 1, 1]
        >>> I.col(0)
        [3]
        [3]
        [3]

        See Also
        ========

        row
        col_swap
        col_del
        col_join
        col_insert
        delRowCol
        """
        if f is None:
            return self[:, j]
        for i in range(0, self.rows):
            self[i, j] = f(self[i, j], i)

    def row_swap(self, i, j):
        """
        Swap the two given rows of the matrix in-place.

        >>> from sympy.matrices import Matrix
        >>> M = Matrix([[0, 1], [1, 0]])
        >>> M
        [0, 1]
        [1, 0]
        >>> M.row_swap(0, 1)
        >>> M
        [1, 0]
        [0, 1]

        See Also
        ========

        row
        col_swap
        """
        for k in range(0, self.cols):
            self[i, k], self[j, k] = self[j, k], self[i, k]

    def col_swap(self, i, j):
        """
        Swap the two given columns of the matrix in-place.

        >>> from sympy.matrices import Matrix
        >>> M = Matrix([[1, 0], [1, 0]])
        >>> M
        [1, 0]
        [1, 0]
        >>> M.col_swap(0, 1)
        >>> M
        [0, 1]
        [0, 1]

        See Also
        ========

        col
        row_swap
        """
        for k in range(0, self.rows):
            self[k, i], self[k, j] = self[k, j], self[k, i]

    def row_del(self, i):
        """
        Delete the given row.

        >>> from sympy.matrices import eye
        >>> M = eye(3)
        >>> M.row_del(1)
        >>> M
        [1, 0, 0]
        [0, 0, 1]

        See Also
        ========

        row
        col_del
        """
        self._mat = self._mat[:i*self.cols] + self._mat[(i+1)*self.cols:]
        self.rows -= 1

    def col_del(self, i):
        """
        Delete the given column.

        >>> from sympy.matrices import eye
        >>> M = eye(3)
        >>> M.col_del(1)
        >>> M
        [1, 0]
        [0, 0]
        [0, 1]

        See Also
        ========

        col
        row_del
        """
        for j in range(self.rows - 1, -1, -1):
            del self._mat[i + j*self.cols]
        self.cols -= 1

    # Utility functions
    def simplify(self, ratio=1.7, measure=count_ops):
        """Applies simplify to the elements of a matrix in place.

        This is a shortcut for M.applyfunc(lambda x: simplify(x, ratio, measure))

        See Also
        ========

        sympy.simplify.simplify.simplify
        """
        for i in range(len(self._mat)):
            self._mat[i] = _simplify(self._mat[i], ratio=ratio, measure=measure)

    def fill(self, value):
        """Fill the matrix with the scalar value.

        See Also
        ========

        zeros
        ones
        """
        self._mat = [value]*len(self)

    def _eval_trace(self):
        """
        Calculate the trace of a square matrix.

        Examples
        ========

        >>> from sympy.matrices import eye
        >>> eye(3).trace()
        3

        """
        trace = 0
        for i in range(self.cols):
            trace += self._mat[i*self.cols + i]
        return trace

    def _eval_transpose(self):
        """
        Matrix transposition.

        >>> from sympy import Matrix, I
        >>> m=Matrix(((1, 2+I), (3, 4)))
        >>> m
        [1, 2 + I]
        [3,     4]
        >>> m.transpose()
        [    1, 3]
        [2 + I, 4]
        >>> m.T == m.transpose()
        True

        See Also
        ========

        conjugate: By-element conjugation
        """
        a = [0]*len(self)
        for i in range(self.cols):
            a[i*self.rows:(i + 1)*self.rows] = self._mat[i::self.cols]
        return self._new(self.cols, self.rows, a)

    def _eval_conjugate(self):
        """By-element conjugation.

        See Also
        ========

        transpose: Matrix transposition
        H: Hermite conjugation
        D: Dirac conjugation
        """
        out = self._new(self.rows, self.cols,
                lambda i, j: self[i, j].conjugate())
        return out

    def _eval_inverse(self, **kwargs):
        """Return the matrix inverse using the method indicated (default
        is Gauss elimination).

        kwargs
        ======

        method : ('GE', 'LU', or 'ADJ')
        iszerofunc
        try_block_diag

        Notes
        =====

        According to the ``method`` keyword, it calls the appropriate method:

          GE .... inverse_GE(); default
          LU .... inverse_LU()
          ADJ ... inverse_ADJ()

        According to the ``try_block_diag`` keyword, it will try to form block
        diagonal matrices using the method get_diag_blocks(), invert these
        individually, and then reconstruct the full inverse matrix.

        Note, the GE and LU methods may require the matrix to be simplified
        before it is inverted in order to properly detect zeros during
        pivoting. In difficult cases a custom zero detection function can
        be provided by setting the ``iszerosfunc`` argument to a function that
        should return True if its argument is zero. The ADJ routine computes
        the determinant and uses that to detect singular matrices in addition
        to testing for zeros on the diagonal.

        See Also
        ========

        inverse_LU
        inverse_GE
        inverse_ADJ
        """
        from sympy.matrices import diag

        method = kwargs.get('method', 'GE')
        iszerofunc = kwargs.get('iszerofunc', _iszero)
        if kwargs.get('try_block_diag', False):
            blocks = self.get_diag_blocks()
            r = []
            for block in blocks:
                r.append(block.inv(method=method, iszerofunc=iszerofunc))
            return diag(*r)
        if method == "GE":
            return self.inverse_GE(iszerofunc=iszerofunc)
        elif method == "LU":
            return self.inverse_LU(iszerofunc=iszerofunc)
        elif method == "ADJ":
            return self.inverse_ADJ(iszerofunc=iszerofunc)
        else:
            # make sure to add an invertibility check (as in inverse_LU)
            # if a new method is added.
            raise ValueError("Inversion method unrecognized")

    def equals(self, other, failing_expression=False):
        """Applies ``equals`` to corresponding elements of the matrices,
        trying to prove that the elements are equivalent, returning True
        if they are, False if any pair is not, and None (or the first
        failing expression if failing_expression is True) if it cannot
        be decided if the expressions are equivalent or not. This is, in
        general, an expensive operation.

        Examples
        ========

        >>> from sympy.matrices import Matrix
        >>> from sympy.abc import x
        >>> from sympy import cos
        >>> A = Matrix([x*(x - 1), 0])
        >>> B = Matrix([x**2 - x, 0])
        >>> A == B
        False
        >>> A.simplify() == B.simplify()
        True
        >>> A.equals(B)
        True
        >>> A.equals(2)
        False

        See Also
        ========
        sympy.core.expr.equals
        """
        try:
            if self.shape != other.shape:
                return False
            rv = True
            for i in range(self.rows):
                for j in range(self.cols):
                    ans = self[i, j].equals(other[i, j], failing_expression)
                    if ans is False:
                        return False
                    elif ans is not True and rv is True:
                        rv = ans
            return rv
        except AttributeError:
            return False

    def __eq__(self, other):
        try:
            if self.shape != other.shape:
                return False
            if isinstance(other, Matrix):
                return self._mat == other._mat
            elif isinstance(other, MatrixBase):
                return self._mat == other.as_mutable()._mat
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return super(MatrixBase, self).__hash__()

    def _cholesky(self):
        """
        Helper function of cholesky.
        Without the error checks.
        To be used privately. """
        L = zeros(self.rows, self.rows)
        for i in range(self.rows):
            for j in range(i):
                L[i, j] = (1 / L[j, j]) * (self[i, j] -
                    sum(L[i, k] * L[j, k] for k in range(j)))
            L[i, i] = sqrt(self[i, i] -
                    sum(L[i, k]**2 for k in range(i)))
        return self._new(L)

    def _LDLdecomposition(self):
        """
        Helper function of LDLdecomposition.
        Without the error checks.
        To be used privately.
        """
        D = zeros(self.rows, self.rows)
        L = eye(self.rows)
        for i in range(self.rows):
            for j in range(i):
                L[i, j] = (1 / D[j, j]) * (self[i, j] - sum(
                    L[i, k] * L[j, k] * D[k, k] for k in range(j)))
            D[i, i] = self[i, i] - sum(L[i, k]**2 * D[k, k]
                for k in range(i))
        return self._new(L), self._new(D)

    def _lower_triangular_solve(self, rhs):
        """
        Helper function of function lower_triangular_solve.
        Without the error checks.
        To be used privately.
        """
        X = zeros(self.rows, 1)
        for i in range(self.rows):
            if self[i, i] == 0:
                raise TypeError("Matrix must be non-singular.")
            X[i, 0] = (rhs[i, 0] - sum(self[i, k] * X[k, 0]
                for k in range(i))) / self[i, i]
        return self._new(X)

    def _upper_triangular_solve(self, rhs):
        """
        Helper function of function upper_triangular_solve.
        Without the error checks, to be used privately. """
        X = zeros(self.rows, 1)
        for i in reversed(range(self.rows)):
            if self[i, i] == 0:
                raise ValueError("Matrix must be non-singular.")
            X[i, 0] = (rhs[i, 0] - sum(self[i, k] * X[k, 0]
                for k in range(i+1, self.rows))) / self[i, i]
        return self._new(X)

    def _diagonal_solve(self, rhs):
        """
        Helper function of function diagonal_solve,
        without the error checks, to be used privately.
        """
        return self._new(rhs.rows, 1, lambda i, j: rhs[i, 0] / self[i, i])

    ############################
    # Mutable matrix operators #
    ############################

    @call_highest_priority('__radd__')
    def __add__(self, other):
        return MatrixBase.__add__(self, _force_mutable(other))

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return MatrixBase.__radd__(self, _force_mutable(other))

    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return MatrixBase.__sub__(self, _force_mutable(other))

    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return MatrixBase.__rsub__(self, _force_mutable(other))

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return MatrixBase.__mul__(self, _force_mutable(other))

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return MatrixBase.__rmul__(self, _force_mutable(other))

    @call_highest_priority('__div__')
    def __div__(self, other):
        return MatrixBase.__div__(self, _force_mutable(other))

    @call_highest_priority('__truediv__')
    def __truediv__(self, other):
        return MatrixBase.__truediv__(self, _force_mutable(other))

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return MatrixBase.__pow__(self, other)

    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("Matrix Power not defined")

def _force_mutable(x):
    """Return a matrix as a Matrix, otherwise return x."""
    if getattr(x, 'is_Matrix', False):
        return x.as_mutable()
    elif isinstance(x, Basic):
        return x
    elif hasattr(x, '__array__'):
        a = x.__array__()
        if len(a.shape) == 0:
            return sympify(a)
        return Matrix(x)
    return x

###########
# Numpy Utility Functions:
# list2numpy, matrix2numpy, symmarray, rot_axis[123]
###########

def list2numpy(l): # pragma: no cover
    """Converts python list of SymPy expressions to a NumPy array.

    See Also
    ========

    matrix2numpy
    """
    from numpy import empty
    a = empty(len(l), dtype=object)
    for i, s in enumerate(l):
        a[i] = s
    return a

def matrix2numpy(m): # pragma: no cover
    """Converts SymPy's matrix to a NumPy array.

    See Also
    ========

    list2numpy
    """
    from numpy import empty
    a = empty(m.shape, dtype=object)
    for i in range(m.rows):
        for j in range(m.cols):
            a[i, j] = m[i, j]
    return a

def symarray(prefix, shape):
    """Create a numpy ndarray of symbols (as an object array).

    The created symbols are named ``prefix_i1_i2_``...  You should thus provide a
    non-empty prefix if you want your symbols to be unique for different output
    arrays, as SymPy symbols with identical names are the same object.

    Parameters
    ----------

    prefix : string
      A prefix prepended to the name of every symbol.

    shape : int or tuple
      Shape of the created array.  If an int, the array is one-dimensional; for
      more than one dimension the shape must be a tuple.

    Examples
    --------
    These doctests require numpy.

    >>> from sympy import symarray
    >>> symarray('', 3) #doctest: +SKIP
    [_0, _1, _2]

    If you want multiple symarrays to contain distinct symbols, you *must*
    provide unique prefixes:

    >>> a = symarray('', 3) #doctest: +SKIP
    >>> b = symarray('', 3) #doctest: +SKIP
    >>> a[0] is b[0] #doctest: +SKIP
    True
    >>> a = symarray('a', 3) #doctest: +SKIP
    >>> b = symarray('b', 3) #doctest: +SKIP
    >>> a[0] is b[0] #doctest: +SKIP
    False

    Creating symarrays with a prefix:

    >>> symarray('a', 3) #doctest: +SKIP
    [a_0, a_1, a_2]

    For more than one dimension, the shape must be given as a tuple:

    >>> symarray('a', (2, 3)) #doctest: +SKIP
    [[a_0_0, a_0_1, a_0_2],
     [a_1_0, a_1_1, a_1_2]]
    >>> symarray('a', (2, 3, 2)) #doctest: +SKIP
    [[[a_0_0_0, a_0_0_1],
      [a_0_1_0, a_0_1_1],
      [a_0_2_0, a_0_2_1]],
    <BLANKLINE>
     [[a_1_0_0, a_1_0_1],
      [a_1_1_0, a_1_1_1],
      [a_1_2_0, a_1_2_1]]]

    """
    try:
        import numpy as np
    except ImportError:
        raise ImportError("symarray requires numpy to be installed")

    arr = np.empty(shape, dtype=object)
    for index in np.ndindex(shape):
        arr[index] = Symbol('%s_%s' % (prefix, '_'.join(map(str, index))))
    return arr

def rot_axis3(theta):
    """Returns a rotation matrix for a rotation of theta (in radians) about
    the 3-axis.

    Examples
    --------

    >>> from sympy import pi
    >>> from sympy.matrices import rot_axis3

    A rotation of pi/3 (60 degrees):

    >>> theta = pi/3
    >>> rot_axis3(theta)
    [       1/2, sqrt(3)/2, 0]
    [-sqrt(3)/2,       1/2, 0]
    [         0,         0, 1]

    If we rotate by pi/2 (90 degrees):

    >>> rot_axis3(pi/2)
    [ 0, 1, 0]
    [-1, 0, 0]
    [ 0, 0, 1]

    See Also
    ========

    rot_axis1: Returns a rotation matrix for a rotation of theta (in radians)
        about the 1-axis
    rot_axis2: Returns a rotation matrix for a rotation of theta (in radians)
        about the 2-axis
    """
    ct = cos(theta)
    st = sin(theta)
    lil = ((ct, st, 0),
           (-st, ct, 0),
           (0, 0, 1))
    return Matrix(lil)

def rot_axis2(theta):
    """Returns a rotation matrix for a rotation of theta (in radians) about
    the 2-axis.

    Examples
    --------

    >>> from sympy import pi
    >>> from sympy.matrices import rot_axis2

    A rotation of pi/3 (60 degrees):

    >>> theta = pi/3
    >>> rot_axis2(theta)
    [      1/2, 0, -sqrt(3)/2]
    [        0, 1,          0]
    [sqrt(3)/2, 0,        1/2]

    If we rotate by pi/2 (90 degrees):

    >>> rot_axis2(pi/2)
    [0, 0, -1]
    [0, 1,  0]
    [1, 0,  0]

    See Also
    ========

    rot_axis1: Returns a rotation matrix for a rotation of theta (in radians)
        about the 1-axis
    rot_axis3: Returns a rotation matrix for a rotation of theta (in radians)
        about the 3-axis
    """
    ct = cos(theta)
    st = sin(theta)
    lil = ((ct, 0, -st),
           (0, 1, 0),
           (st, 0, ct))
    return Matrix(lil)

def rot_axis1(theta):
    """Returns a rotation matrix for a rotation of theta (in radians) about
    the 1-axis.

    Examples
    --------

    >>> from sympy import pi
    >>> from sympy.matrices import rot_axis1

    A rotation of pi/3 (60 degrees):

    >>> theta = pi/3
    >>> rot_axis1(theta)
    [1,          0,         0]
    [0,        1/2, sqrt(3)/2]
    [0, -sqrt(3)/2,       1/2]

    If we rotate by pi/2 (90 degrees):

    >>> rot_axis1(pi/2)
    [1,  0, 0]
    [0,  0, 1]
    [0, -1, 0]

    See Also
    ========

    rot_axis2: Returns a rotation matrix for a rotation of theta (in radians)
        about the 2-axis
    rot_axis3: Returns a rotation matrix for a rotation of theta (in radians)
        about the 3-axis
    """
    ct = cos(theta)
    st = sin(theta)
    lil = ((1, 0, 0),
           (0, ct, st),
           (0, -st, ct))
    return Matrix(lil)

###############
# Functions
###############
def matrix_add(A, B):
    SymPyDeprecationWarning(
        feature="matrix_add(A, B)",
        useinstead="A + B",
        deprecated_since_version="0.7.2",
    ).warn()
    return A + B

def matrix_multiply(A, B):
    SymPyDeprecationWarning(
        feature="matrix_multiply(A, B)",
        useinstead="A * B",
        deprecated_since_version="0.7.2",
    ).warn()
    return A*B

def matrix_multiply_elementwise(A, B):
    """Return the Hadamard product (elementwise product) of A and B

    >>> from sympy.matrices import matrix_multiply_elementwise
    >>> from sympy.matrices import Matrix
    >>> A = Matrix([[0, 1, 2], [3, 4, 5]])
    >>> B = Matrix([[1, 10, 100], [100, 10, 1]])
    >>> matrix_multiply_elementwise(A, B)
    [  0, 10, 200]
    [300, 40,   5]

    See Also
    ========

    __mul__
    """
    if A.shape != B.shape:
        raise ShapeError()
    shape = A.shape
    return classof(A, B)._new(shape[0], shape[1],
        lambda i, j: A[i, j] * B[i, j])

def zeros(r, c=None, cls=None):
    """Returns a matrix of zeros with ``r`` rows and ``c`` columns;
    if ``c`` is omitted a square matrix will be returned.

    See Also
    ========

    ones
    eye
    diag
    """
    if cls is None:
        from dense import Matrix as cls
    if is_sequence(r):
        SymPyDeprecationWarning(
            feature="The syntax zeros([%i, %i])" % tuple(r),
            useinstead="zeros(%i, %i)." % tuple(r),
            issue=3381, deprecated_since_version="0.7.2",
        ).warn()
        r, c = r
    else:
        c = r if c is None else c
    r, c = [int(i) for i in [r, c]]
    return cls.zeros(r, c)

def ones(r, c=None):
    """Returns a matrix of ones with ``r`` rows and ``c`` columns;
    if ``c`` is omitted a square matrix will be returned.

    See Also
    ========

    zeros
    eye
    diag
    """
    from dense import Matrix

    if is_sequence(r):
        SymPyDeprecationWarning(
                feature="The syntax ones([%i, %i])" % tuple(r),
                useinstead="ones(%i, %i)." % tuple(r),
                issue=3381, deprecated_since_version="0.7.2",
        ).warn()
        r, c = r
    else:
        c = r if c is None else c
    r = as_int(r)
    c = as_int(c)
    return Matrix(r, c, [S.One]*r*c)

def eye(n, cls=None):
    """Create square identity matrix n x n

    See Also
    ========

    diag
    zeros
    ones
    """
    if cls is None:
        from dense import Matrix as cls

    n = as_int(n)
    out = cls.zeros(n)
    for i in range(n):
        out[i, i] = S.One
    return out

def diag(*values):
    """Create a diagonal matrix from a list as a diagonal values.

    When arguments are matrices they are fitted in resultant matrix.

    Examples
    ========

    >>> from sympy.matrices import diag, Matrix, ones
    >>> diag(1, 2, 3)
    [1, 0, 0]
    [0, 2, 0]
    [0, 0, 3]

    >>> from sympy.abc import x, y, z
    >>> a = Matrix([x, y, z])
    >>> b = Matrix([[1, 2], [3, 4]])
    >>> c = Matrix([[5, 6]])
    >>> diag(a, 7, b, c)
    [x, 0, 0, 0, 0, 0]
    [y, 0, 0, 0, 0, 0]
    [z, 0, 0, 0, 0, 0]
    [0, 7, 0, 0, 0, 0]
    [0, 0, 1, 2, 0, 0]
    [0, 0, 3, 4, 0, 0]
    [0, 0, 0, 0, 5, 6]

    A given band off the diagonal can be made by padding with a
    vertical or horizontal "kerning" vector:

    >>> hpad = ones(0, 2)
    >>> vpad = ones(2, 0)
    >>> diag(vpad, 1, 2, 3, hpad) + diag(hpad, 4, 5, 6, vpad)
    [0, 0, 4, 0, 0]
    [0, 0, 0, 5, 0]
    [1, 0, 0, 0, 6]
    [0, 2, 0, 0, 0]
    [0, 0, 3, 0, 0]

    See Also
    ========

    eye
    """
    from sympy.matrices import zeros

    rows = 0
    cols = 0
    for m in values:
        if isinstance(m, MatrixBase):
            rows += m.rows
            cols += m.cols
        else:
            rows += 1
            cols += 1
    res = zeros(rows, cols)
    i_row = 0
    i_col = 0
    for m in values:
        if isinstance(m, MatrixBase):
            res[i_row:i_row + m.rows, i_col:i_col + m.cols] = m
            i_row += m.rows
            i_col += m.cols
        else:
            res[i_row, i_col] = m
            i_row += 1
            i_col += 1
    return res

def jordan_cell(eigenval, n):
    """
    Create matrix of Jordan cell kind:

    Examples
    ========

    >>> from sympy.matrices import jordan_cell
    >>> from sympy.abc import x
    >>> jordan_cell(x, 4)
    [x, 1, 0, 0]
    [0, x, 1, 0]
    [0, 0, x, 1]
    [0, 0, 0, x]
    """
    n = as_int(n)
    out = zeros(n)
    for i in range(n - 1):
        out[i, i] = eigenval
        out[i, i + 1] = S.One
    out[n - 1, n - 1] = eigenval
    return out

def hessian(f, varlist):
    """Compute Hessian matrix for a function f

    see: http://en.wikipedia.org/wiki/Hessian_matrix

    See Also
    ========

    sympy.matrices.mutable.Matrix.jacobian
    wronskian
    """
    # f is the expression representing a function f, return regular matrix
    if is_sequence(varlist):
        m = len(varlist)
        if not m:
            raise ShapeError("`len(varlist)` must not be zero.")
    elif isinstance(varlist, MatrixBase):
        m = varlist.cols
        if not m:
            raise ShapeError("`varlist.cols` must not be zero.")
        if varlist.rows != 1:
            raise ShapeError("`varlist` must be a row vector.")
    else:
        raise ValueError("Improper variable list in hessian function")
    if not getattr(f, 'diff'):
        # check differentiability
        raise ValueError("Function `f` (%s) is not differentiable" % f)
    out = zeros(m)
    for i in range(m):
        for j in range(i, m):
            out[i, j] = f.diff(varlist[i]).diff(varlist[j])
    for i in range(m):
        for j in range(i):
            out[i, j] = out[j, i]
    return out

def GramSchmidt(vlist, orthog=False):
    """
    Apply the Gram-Schmidt process to a set of vectors.

    see: http://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
    """
    out = []
    m = len(vlist)
    for i in range(m):
        tmp = vlist[i]
        for j in range(i):
            tmp -= vlist[i].project(out[j])
        if not tmp.values():
            raise ValueError("GramSchmidt: vector set not linearly independent")
        out.append(tmp)
    if orthog:
        for i in range(len(out)):
            out[i] = out[i].normalized()
    return out

def wronskian(functions, var, method='bareis'):
    """
    Compute Wronskian for [] of functions

    ::

                         | f1       f2        ...   fn      |
                         | f1'      f2'       ...   fn'     |
                         |  .        .        .      .      |
        W(f1, ..., fn) = |  .        .         .     .      |
                         |  .        .          .    .      |
                         |  (n)      (n)            (n)     |
                         | D   (f1) D   (f2)  ...  D   (fn) |

    see: http://en.wikipedia.org/wiki/Wronskian

    See Also
    ========

    sympy.matrices.mutable.Matrix.jacobian
    hessian
    """
    from dense import Matrix

    for index in range(0, len(functions)):
        functions[index] = sympify(functions[index])
    n = len(functions)
    if n == 0:
        return 1
    W = Matrix(n, n, lambda i, j: functions[i].diff(var, j) )
    return W.det(method)

def casoratian(seqs, n, zero=True):
    """Given linear difference operator L of order 'k' and homogeneous
       equation Ly = 0 we want to compute kernel of L, which is a set
       of 'k' sequences: a(n), b(n), ... z(n).

       Solutions of L are linearly independent iff their Casoratian,
       denoted as C(a, b, ..., z), do not vanish for n = 0.

       Casoratian is defined by k x k determinant::

                  +  a(n)     b(n)     . . . z(n)     +
                  |  a(n+1)   b(n+1)   . . . z(n+1)   |
                  |    .         .     .        .     |
                  |    .         .       .      .     |
                  |    .         .         .    .     |
                  +  a(n+k-1) b(n+k-1) . . . z(n+k-1) +

       It proves very useful in rsolve_hyper() where it is applied
       to a generating set of a recurrence to factor out linearly
       dependent solutions and return a basis:

       >>> from sympy import Symbol, casoratian, factorial
       >>> n = Symbol('n', integer=True)

       Exponential and factorial are linearly independent:

       >>> casoratian([2**n, factorial(n)], n) != 0
       True

    """
    from dense import Matrix

    seqs = map(sympify, seqs)

    if not zero:
        f = lambda i, j: seqs[j].subs(n, n+i)
    else:
        f = lambda i, j: seqs[j].subs(n, i)

    k = len(seqs)

    return Matrix(k, k, f).det()

def randMatrix(r, c=None, min=0, max=99, seed=None, symmetric=False, percent=100):
    """Create random matrix with dimensions ``r`` x ``c``. If ``c`` is omitted
    the matrix will be square. If ``symmetric`` is True the matrix must be
    square. If ``percent`` is less than 100 then only approximately the given
    percentage of elements will be non-zero.

    Examples
    ========

    >>> from sympy.matrices import randMatrix
    >>> randMatrix(3) # doctest:+SKIP
    [25, 45, 27]
    [44, 54,  9]
    [23, 96, 46]
    >>> randMatrix(3, 2) # doctest:+SKIP
    [87, 29]
    [23, 37]
    [90, 26]
    >>> randMatrix(3, 3, 0, 2) # doctest:+SKIP
    [0, 2, 0]
    [2, 0, 1]
    [0, 0, 1]
    >>> randMatrix(3, symmetric=True) # doctest:+SKIP
    [85, 26, 29]
    [26, 71, 43]
    [29, 43, 57]
    >>> A = randMatrix(3, seed=1)
    >>> B = randMatrix(3, seed=2)
    >>> A == B # doctest:+SKIP
    False
    >>> A == randMatrix(3, seed=1)
    True
    >>> randMatrix(3, symmetric=True, percent=50) # doctest:+SKIP
    [0, 68, 43]
    [0, 68,  0]
    [0, 91, 34]
    """
    if c is None:
        c = r
    if seed is None:
        prng = random.Random()  # use system time
    else:
        prng = random.Random(seed)
    if symmetric and r != c:
        raise ValueError('For symmetric matrices, r must equal c, but %i != %i' % (r, c))
    if not symmetric:
        m = Matrix._new(r, c, lambda i, j: prng.randint(min, max))
    else:
        m = zeros(r)
        for i in range(r):
            for j in range(i, r):
                m[i, j] = prng.randint(min, max)
        for i in range(r):
            for j in range(i):
                m[i, j] = m[j, i]
    if percent == 100:
        return m
    else:
        z = int(r*c*percent//100)
        m._mat[:z] = [S.Zero]*z
        random.shuffle(m._mat)
    return m

# for compatibility
MutableMatrix = Matrix
