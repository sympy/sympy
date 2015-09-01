from __future__ import print_function, division

import random

from sympy.core.basic import Basic
from sympy.core.compatibility import is_sequence, as_int, range
from sympy.core.function import count_ops
from sympy.core.decorators import call_highest_priority
from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.core.sympify import sympify
from sympy.functions.elementary.trigonometric import cos, sin
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.simplify import simplify as _simplify
from sympy.utilities.misc import filldedent
from sympy.utilities.decorator import doctest_depends_on

from sympy.matrices.matrices import (MatrixBase,
    ShapeError, a2idx, classof)


def _iszero(x):
    """Returns True if x is zero."""
    return x.is_zero


class DenseMatrix(MatrixBase):

    is_MatrixExpr = False

    _op_priority = 10.01
    _class_priority = 4

    def __getitem__(self, key):
        """Return portion of self defined by key. If the key involves a slice
        then a list will be returned (if key is a single slice) or a matrix
        (if key was a tuple involving a slice).

        Examples
        ========

        >>> from sympy import Matrix, I
        >>> m = Matrix([
        ... [1, 2 + I],
        ... [3, 4    ]])

        If the key is a tuple that doesn't involve a slice then that element
        is returned:

        >>> m[1, 0]
        3

        When a tuple key involves a slice, a matrix is returned. Here, the
        first column is selected (all rows, column 0):

        >>> m[:, 0]
        Matrix([
        [1],
        [3]])

        If the slice is not a tuple then it selects from the underlying
        list of elements that are arranged in row order and a list is
        returned if a slice is involved:

        >>> m[0]
        1
        >>> m[::2]
        [1, 3]
        """
        if isinstance(key, tuple):
            i, j = key
            try:
                i, j = self.key2ij(key)
                return self._mat[i*self.cols + j]
            except (TypeError, IndexError):
                if isinstance(i, slice):
                    # XXX remove list() when PY2 support is dropped
                    i = list(range(self.rows))[i]
                elif is_sequence(i):
                    pass
                else:
                    i = [i]
                if isinstance(j, slice):
                    # XXX remove list() when PY2 support is dropped
                    j = list(range(self.cols))[j]
                elif is_sequence(j):
                    pass
                else:
                    j = [j]
                return self.extract(i, j)
        else:
            # row-wise decomposition of matrix
            if isinstance(key, slice):
                return self._mat[key]
            return self._mat[a2idx(key)]

    def __setitem__(self, key, value):
        raise NotImplementedError()

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

    def tolist(self):
        """Return the Matrix as a nested Python list.

        Examples
        ========

        >>> from sympy import Matrix, ones
        >>> m = Matrix(3, 3, range(9))
        >>> m
        Matrix([
        [0, 1, 2],
        [3, 4, 5],
        [6, 7, 8]])
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

    def row(self, i):
        """Elementary row selector.

        Examples
        ========

        >>> from sympy import eye
        >>> eye(2).row(0)
        Matrix([[1, 0]])

        See Also
        ========

        col
        row_op
        row_swap
        row_del
        row_join
        row_insert
        """
        return self[i, :]

    def col(self, j):
        """Elementary column selector.

        Examples
        ========

        >>> from sympy import eye
        >>> eye(2).col(0)
        Matrix([
        [1],
        [0]])

        See Also
        ========

        row
        col_op
        col_swap
        col_del
        col_join
        col_insert
        """
        return self[:, j]

    def _eval_trace(self):
        """Calculate the trace of a square matrix.

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

    def _eval_determinant(self):
        return self.det()

    def _eval_transpose(self):
        """Matrix transposition.

        Examples
        ========

        >>> from sympy import Matrix, I
        >>> m=Matrix(((1, 2+I), (3, 4)))
        >>> m
        Matrix([
        [1, 2 + I],
        [3,     4]])
        >>> m.transpose()
        Matrix([
        [    1, 3],
        [2 + I, 4]])
        >>> m.T == m.transpose()
        True

        See Also
        ========

        conjugate: By-element conjugation
        """
        a = []
        for i in range(self.cols):
            a.extend(self._mat[i::self.cols])
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

    def _eval_adjoint(self):
        return self.T.C

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

        M = self.as_mutable()
        if method == "GE":
            rv = M.inverse_GE(iszerofunc=iszerofunc)
        elif method == "LU":
            rv = M.inverse_LU(iszerofunc=iszerofunc)
        elif method == "ADJ":
            rv = M.inverse_ADJ(iszerofunc=iszerofunc)
        else:
            # make sure to add an invertibility check (as in inverse_LU)
            # if a new method is added.
            raise ValueError("Inversion method unrecognized")
        return self._new(rv)

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
                return self._mat == Matrix(other)._mat
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self == other

    def _cholesky(self):
        """Helper function of cholesky.
        Without the error checks.
        To be used privately. """
        L = zeros(self.rows, self.rows)
        for i in range(self.rows):
            for j in range(i):
                L[i, j] = (1 / L[j, j])*(self[i, j] -
                    sum(L[i, k]*L[j, k] for k in range(j)))
            L[i, i] = sqrt(self[i, i] -
                    sum(L[i, k]**2 for k in range(i)))
        return self._new(L)

    def _LDLdecomposition(self):
        """Helper function of LDLdecomposition.
        Without the error checks.
        To be used privately.
        """
        D = zeros(self.rows, self.rows)
        L = eye(self.rows)
        for i in range(self.rows):
            for j in range(i):
                L[i, j] = (1 / D[j, j])*(self[i, j] - sum(
                    L[i, k]*L[j, k]*D[k, k] for k in range(j)))
            D[i, i] = self[i, i] - sum(L[i, k]**2*D[k, k]
                for k in range(i))
        return self._new(L), self._new(D)

    def _lower_triangular_solve(self, rhs):
        """Helper function of function lower_triangular_solve.
        Without the error checks.
        To be used privately.
        """
        X = zeros(self.rows, rhs.cols)
        for j in range(rhs.cols):
            for i in range(self.rows):
                if self[i, i] == 0:
                    raise TypeError("Matrix must be non-singular.")
                X[i, j] = (rhs[i, j] - sum(self[i, k]*X[k, j]
                    for k in range(i))) / self[i, i]
        return self._new(X)

    def _upper_triangular_solve(self, rhs):
        """Helper function of function upper_triangular_solve.
        Without the error checks, to be used privately. """
        X = zeros(self.rows, rhs.cols)
        for j in range(rhs.cols):
            for i in reversed(range(self.rows)):
                if self[i, i] == 0:
                    raise ValueError("Matrix must be non-singular.")
                X[i, j] = (rhs[i, j] - sum(self[i, k]*X[k, j]
                    for k in range(i + 1, self.rows))) / self[i, i]
        return self._new(X)

    def _diagonal_solve(self, rhs):
        """Helper function of function diagonal_solve,
        without the error checks, to be used privately.
        """
        return self._new(rhs.rows, rhs.cols, lambda i, j: rhs[i, j] / self[i, i])

    def applyfunc(self, f):
        """Apply a function to each element of the matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(2, 2, lambda i, j: i*2+j)
        >>> m
        Matrix([
        [0, 1],
        [2, 3]])
        >>> m.applyfunc(lambda i: 2*i)
        Matrix([
        [0, 2],
        [4, 6]])

        """
        if not callable(f):
            raise TypeError("`f` must be callable.")

        out = self._new(self.rows, self.cols, list(map(f, self._mat)))
        return out

    def reshape(self, rows, cols):
        """Reshape the matrix. Total number of elements must remain the same.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(2, 3, lambda i, j: 1)
        >>> m
        Matrix([
        [1, 1, 1],
        [1, 1, 1]])
        >>> m.reshape(1, 6)
        Matrix([[1, 1, 1, 1, 1, 1]])
        >>> m.reshape(3, 2)
        Matrix([
        [1, 1],
        [1, 1],
        [1, 1]])

        """
        if len(self) != rows*cols:
            raise ValueError("Invalid reshape parameters %d %d" % (rows, cols))
        return self._new(rows, cols, lambda i, j: self._mat[i*cols + j])

    def as_mutable(self):
        """Returns a mutable version of this matrix

        Examples
        ========

        >>> from sympy import ImmutableMatrix
        >>> X = ImmutableMatrix([[1, 2], [3, 4]])
        >>> Y = X.as_mutable()
        >>> Y[1, 1] = 5 # Can set values in Y
        >>> Y
        Matrix([
        [1, 2],
        [3, 5]])
        """
        return Matrix(self)

    def as_immutable(self):
        """Returns an Immutable version of this Matrix
        """
        from .immutable import ImmutableMatrix as cls
        if self.rows and self.cols:
            return cls._new(self.tolist())
        return cls._new(self.rows, self.cols, [])

    @classmethod
    def zeros(cls, r, c=None):
        """Return an r x c matrix of zeros, square if c is omitted."""
        c = r if c is None else c
        r = as_int(r)
        c = as_int(c)
        return cls._new(r, c, [cls._sympify(0)]*r*c)

    @classmethod
    def eye(cls, n):
        """Return an n x n identity matrix."""
        n = as_int(n)
        mat = [cls._sympify(0)]*n*n
        mat[::n + 1] = [cls._sympify(1)]*n
        return cls._new(n, n, mat)

    ############################
    # Mutable matrix operators #
    ############################

    @call_highest_priority('__radd__')
    def __add__(self, other):
        return super(DenseMatrix, self).__add__(_force_mutable(other))

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return super(DenseMatrix, self).__radd__(_force_mutable(other))

    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return super(DenseMatrix, self).__sub__(_force_mutable(other))

    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return super(DenseMatrix, self).__rsub__(_force_mutable(other))

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return super(DenseMatrix, self).__mul__(_force_mutable(other))

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return super(DenseMatrix, self).__rmul__(_force_mutable(other))

    @call_highest_priority('__div__')
    def __div__(self, other):
        return super(DenseMatrix, self).__div__(_force_mutable(other))

    @call_highest_priority('__truediv__')
    def __truediv__(self, other):
        return super(DenseMatrix, self).__truediv__(_force_mutable(other))

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return super(DenseMatrix, self).__pow__(other)

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


class MutableDenseMatrix(DenseMatrix, MatrixBase):
    @classmethod
    def _new(cls, *args, **kwargs):
        rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)
        self = object.__new__(cls)
        self.rows = rows
        self.cols = cols
        self._mat = list(flat_list)  # create a shallow copy
        return self

    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    def as_mutable(self):
        return self.copy()

    def __setitem__(self, key, value):
        """

        Examples
        ========

        >>> from sympy import Matrix, I, zeros, ones
        >>> m = Matrix(((1, 2+I), (3, 4)))
        >>> m
        Matrix([
        [1, 2 + I],
        [3,     4]])
        >>> m[1, 0] = 9
        >>> m
        Matrix([
        [1, 2 + I],
        [9,     4]])
        >>> m[1, 0] = [[0, 1]]

        To replace row r you assign to position r*m where m
        is the number of columns:

        >>> M = zeros(4)
        >>> m = M.cols
        >>> M[3*m] = ones(1, m)*2; M
        Matrix([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [2, 2, 2, 2]])

        And to replace column c you can assign to position c:

        >>> M[2] = ones(m, 1)*4; M
        Matrix([
        [0, 0, 4, 0],
        [0, 0, 4, 0],
        [0, 0, 4, 0],
        [2, 2, 4, 2]])
        """
        rv = self._setitem(key, value)
        if rv is not None:
            i, j, value = rv
            self._mat[i*self.cols + j] = value

    def copyin_matrix(self, key, value):
        """Copy in values from a matrix into the given bounds.

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
        Matrix([
        [0, 1, 0],
        [2, 3, 0],
        [4, 5, 1]])
        >>> I[0, 1] = M
        >>> I
        Matrix([
        [0, 0, 1],
        [2, 2, 3],
        [4, 4, 5]])

        See Also
        ========

        copyin_list
        """
        rlo, rhi, clo, chi = self.key2bounds(key)
        shape = value.shape
        dr, dc = rhi - rlo, chi - clo
        if shape != (dr, dc):
            raise ShapeError(filldedent("The Matrix `value` doesn't have the "
                "same dimensions "
                "as the in sub-Matrix given by `key`."))

        for i in range(value.rows):
            for j in range(value.cols):
                self[i + rlo, j + clo] = value[i, j]

    def copyin_list(self, key, value):
        """Copy in elements from a list.

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
        Matrix([
        [1, 0, 0],
        [2, 1, 0],
        [0, 0, 1]])
        >>> I[1, :2] = [[3, 4]]
        >>> I
        Matrix([
        [1, 0, 0],
        [3, 4, 0],
        [0, 0, 1]])

        See Also
        ========

        copyin_matrix
        """
        if not is_sequence(value):
            raise TypeError("`value` must be an ordered iterable, not %s." % type(value))
        return self.copyin_matrix(key, Matrix(value))

    def zip_row_op(self, i, k, f):
        """In-place operation on row ``i`` using two-arg functor whose args are
        interpreted as ``(self[i, j], self[k, j])``.

        Examples
        ========

        >>> from sympy.matrices import eye
        >>> M = eye(3)
        >>> M.zip_row_op(1, 0, lambda v, u: v + 2*u); M
        Matrix([
        [1, 0, 0],
        [2, 1, 0],
        [0, 0, 1]])

        See Also
        ========
        row
        row_op
        col_op

        """
        i0 = i*self.cols
        k0 = k*self.cols

        ri = self._mat[i0: i0 + self.cols]
        rk = self._mat[k0: k0 + self.cols]

        self._mat[i0: i0 + self.cols] = [ f(x, y) for x, y in zip(ri, rk) ]

    def row_op(self, i, f):
        """In-place operation on row ``i`` using two-arg functor whose args are
        interpreted as ``(self[i, j], j)``.

        Examples
        ========

        >>> from sympy.matrices import eye
        >>> M = eye(3)
        >>> M.row_op(1, lambda v, j: v + 2*M[0, j]); M
        Matrix([
        [1, 0, 0],
        [2, 1, 0],
        [0, 0, 1]])

        See Also
        ========
        row
        zip_row_op
        col_op

        """
        i0 = i*self.cols
        ri = self._mat[i0: i0 + self.cols]
        self._mat[i0: i0 + self.cols] = [ f(x, j) for x, j in zip(ri, list(range(self.cols))) ]

    def col_op(self, j, f):
        """In-place operation on col j using two-arg functor whose args are
        interpreted as (self[i, j], i).

        Examples
        ========

        >>> from sympy.matrices import eye
        >>> M = eye(3)
        >>> M.col_op(1, lambda v, i: v + 2*M[i, 0]); M
        Matrix([
        [1, 2, 0],
        [0, 1, 0],
        [0, 0, 1]])

        See Also
        ========
        col
        row_op
        """
        self._mat[j::self.cols] = [f(*t) for t in list(zip(self._mat[j::self.cols], list(range(self.rows))))]

    def row_swap(self, i, j):
        """Swap the two given rows of the matrix in-place.

        Examples
        ========

        >>> from sympy.matrices import Matrix
        >>> M = Matrix([[0, 1], [1, 0]])
        >>> M
        Matrix([
        [0, 1],
        [1, 0]])
        >>> M.row_swap(0, 1)
        >>> M
        Matrix([
        [1, 0],
        [0, 1]])

        See Also
        ========

        row
        col_swap
        """
        for k in range(0, self.cols):
            self[i, k], self[j, k] = self[j, k], self[i, k]

    def col_swap(self, i, j):
        """Swap the two given columns of the matrix in-place.

        Examples
        ========

        >>> from sympy.matrices import Matrix
        >>> M = Matrix([[1, 0], [1, 0]])
        >>> M
        Matrix([
        [1, 0],
        [1, 0]])
        >>> M.col_swap(0, 1)
        >>> M
        Matrix([
        [0, 1],
        [0, 1]])

        See Also
        ========

        col
        row_swap
        """
        for k in range(0, self.rows):
            self[k, i], self[k, j] = self[k, j], self[k, i]

    def row_del(self, i):
        """Delete the given row.

        Examples
        ========

        >>> from sympy.matrices import eye
        >>> M = eye(3)
        >>> M.row_del(1)
        >>> M
        Matrix([
        [1, 0, 0],
        [0, 0, 1]])

        See Also
        ========

        row
        col_del
        """
        if i < -self.rows or i >= self.rows:
            raise IndexError("Index out of range: 'i = %s', valid -%s <= i"
                             " < %s" % (i, self.rows, self.rows))
        if i < 0:
            i += self.rows
        del self._mat[i*self.cols:(i+1)*self.cols]
        self.rows -= 1

    def col_del(self, i):
        """Delete the given column.

        Examples
        ========

        >>> from sympy.matrices import eye
        >>> M = eye(3)
        >>> M.col_del(1)
        >>> M
        Matrix([
        [1, 0],
        [0, 0],
        [0, 1]])

        See Also
        ========

        col
        row_del
        """
        if i < -self.cols or i >= self.cols:
            raise IndexError("Index out of range: 'i=%s', valid -%s <= i < %s"
                             % (i, self.cols, self.cols))
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
            self._mat[i] = _simplify(self._mat[i], ratio=ratio,
                                     measure=measure)

    def fill(self, value):
        """Fill the matrix with the scalar value.

        See Also
        ========

        zeros
        ones
        """
        self._mat = [value]*len(self)

MutableMatrix = Matrix = MutableDenseMatrix

###########
# Numpy Utility Functions:
# list2numpy, matrix2numpy, symmarray, rot_axis[123]
###########


def list2numpy(l, dtype=object):  # pragma: no cover
    """Converts python list of SymPy expressions to a NumPy array.

    See Also
    ========

    matrix2numpy
    """
    from numpy import empty
    a = empty(len(l), dtype)
    for i, s in enumerate(l):
        a[i] = s
    return a


def matrix2numpy(m, dtype=object):  # pragma: no cover
    """Converts SymPy's matrix to a NumPy array.

    See Also
    ========

    list2numpy
    """
    from numpy import empty
    a = empty(m.shape, dtype)
    for i in range(m.rows):
        for j in range(m.cols):
            a[i, j] = m[i, j]
    return a

@doctest_depends_on(modules=('numpy',))
def symarray(prefix, shape):  # pragma: no cover
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
    ========
    These doctests require numpy.

    >>> from sympy import symarray
    >>> symarray('', 3)
    [_0 _1 _2]

    If you want multiple symarrays to contain distinct symbols, you *must*
    provide unique prefixes:

    >>> a = symarray('', 3)
    >>> b = symarray('', 3)
    >>> a[0] == b[0]
    True
    >>> a = symarray('a', 3)
    >>> b = symarray('b', 3)
    >>> a[0] == b[0]
    False

    Creating symarrays with a prefix:

    >>> symarray('a', 3)
    [a_0 a_1 a_2]

    For more than one dimension, the shape must be given as a tuple:

    >>> symarray('a', (2, 3))
    [[a_0_0 a_0_1 a_0_2]
     [a_1_0 a_1_1 a_1_2]]
    >>> symarray('a', (2, 3, 2))
    [[[a_0_0_0 a_0_0_1]
      [a_0_1_0 a_0_1_1]
      [a_0_2_0 a_0_2_1]]
    <BLANKLINE>
     [[a_1_0_0 a_1_0_1]
      [a_1_1_0 a_1_1_1]
      [a_1_2_0 a_1_2_1]]]

    """
    from numpy import empty, ndindex
    arr = empty(shape, dtype=object)
    for index in ndindex(shape):
        arr[index] = Symbol('%s_%s' % (prefix, '_'.join(map(str, index))))
    return arr


def rot_axis3(theta):
    """Returns a rotation matrix for a rotation of theta (in radians) about
    the 3-axis.

    Examples
    ========

    >>> from sympy import pi
    >>> from sympy.matrices import rot_axis3

    A rotation of pi/3 (60 degrees):

    >>> theta = pi/3
    >>> rot_axis3(theta)
    Matrix([
    [       1/2, sqrt(3)/2, 0],
    [-sqrt(3)/2,       1/2, 0],
    [         0,         0, 1]])

    If we rotate by pi/2 (90 degrees):

    >>> rot_axis3(pi/2)
    Matrix([
    [ 0, 1, 0],
    [-1, 0, 0],
    [ 0, 0, 1]])

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
    ========

    >>> from sympy import pi
    >>> from sympy.matrices import rot_axis2

    A rotation of pi/3 (60 degrees):

    >>> theta = pi/3
    >>> rot_axis2(theta)
    Matrix([
    [      1/2, 0, -sqrt(3)/2],
    [        0, 1,          0],
    [sqrt(3)/2, 0,        1/2]])

    If we rotate by pi/2 (90 degrees):

    >>> rot_axis2(pi/2)
    Matrix([
    [0, 0, -1],
    [0, 1,  0],
    [1, 0,  0]])

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
    ========

    >>> from sympy import pi
    >>> from sympy.matrices import rot_axis1

    A rotation of pi/3 (60 degrees):

    >>> theta = pi/3
    >>> rot_axis1(theta)
    Matrix([
    [1,          0,         0],
    [0,        1/2, sqrt(3)/2],
    [0, -sqrt(3)/2,       1/2]])

    If we rotate by pi/2 (90 degrees):

    >>> rot_axis1(pi/2)
    Matrix([
    [1,  0, 0],
    [0,  0, 1],
    [0, -1, 0]])

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


def matrix_multiply_elementwise(A, B):
    """Return the Hadamard product (elementwise product) of A and B

    >>> from sympy.matrices import matrix_multiply_elementwise
    >>> from sympy.matrices import Matrix
    >>> A = Matrix([[0, 1, 2], [3, 4, 5]])
    >>> B = Matrix([[1, 10, 100], [100, 10, 1]])
    >>> matrix_multiply_elementwise(A, B)
    Matrix([
    [  0, 10, 200],
    [300, 40,   5]])

    See Also
    ========

    __mul__
    """
    if A.shape != B.shape:
        raise ShapeError()
    shape = A.shape
    return classof(A, B)._new(shape[0], shape[1],
        lambda i, j: A[i, j]*B[i, j])


def ones(r, c=None):
    """Returns a matrix of ones with ``r`` rows and ``c`` columns;
    if ``c`` is omitted a square matrix will be returned.

    See Also
    ========

    zeros
    eye
    diag
    """
    from .dense import Matrix

    c = r if c is None else c
    r = as_int(r)
    c = as_int(c)
    return Matrix(r, c, [S.One]*r*c)


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
        from .dense import Matrix as cls
    return cls.zeros(r, c)


def eye(n, cls=None):
    """Create square identity matrix n x n

    See Also
    ========

    diag
    zeros
    ones
    """
    if cls is None:
        from sympy.matrices import Matrix as cls
    return cls.eye(n)

def diag(*values, **kwargs):
    """Create a sparse, diagonal matrix from a list of diagonal values.

    Notes
    =====

    When arguments are matrices they are fitted in resultant matrix.

    The returned matrix is a mutable, dense matrix. To make it a different
    type, send the desired class for keyword ``cls``.

    Examples
    ========

    >>> from sympy.matrices import diag, Matrix, ones
    >>> diag(1, 2, 3)
    Matrix([
    [1, 0, 0],
    [0, 2, 0],
    [0, 0, 3]])
    >>> diag(*[1, 2, 3])
    Matrix([
    [1, 0, 0],
    [0, 2, 0],
    [0, 0, 3]])

    The diagonal elements can be matrices; diagonal filling will
    continue on the diagonal from the last element of the matrix:

    >>> from sympy.abc import x, y, z
    >>> a = Matrix([x, y, z])
    >>> b = Matrix([[1, 2], [3, 4]])
    >>> c = Matrix([[5, 6]])
    >>> diag(a, 7, b, c)
    Matrix([
    [x, 0, 0, 0, 0, 0],
    [y, 0, 0, 0, 0, 0],
    [z, 0, 0, 0, 0, 0],
    [0, 7, 0, 0, 0, 0],
    [0, 0, 1, 2, 0, 0],
    [0, 0, 3, 4, 0, 0],
    [0, 0, 0, 0, 5, 6]])

    When diagonal elements are lists, they will be treated as arguments
    to Matrix:

    >>> diag([1, 2, 3], 4)
    Matrix([
    [1, 0],
    [2, 0],
    [3, 0],
    [0, 4]])
    >>> diag([[1, 2, 3]], 4)
    Matrix([
    [1, 2, 3, 0],
    [0, 0, 0, 4]])

    A given band off the diagonal can be made by padding with a
    vertical or horizontal "kerning" vector:

    >>> hpad = ones(0, 2)
    >>> vpad = ones(2, 0)
    >>> diag(vpad, 1, 2, 3, hpad) + diag(hpad, 4, 5, 6, vpad)
    Matrix([
    [0, 0, 4, 0, 0],
    [0, 0, 0, 5, 0],
    [1, 0, 0, 0, 6],
    [0, 2, 0, 0, 0],
    [0, 0, 3, 0, 0]])



    The type is mutable by default but can be made immutable by setting
    the ``mutable`` flag to False:

    >>> type(diag(1))
    <class 'sympy.matrices.dense.MutableDenseMatrix'>
    >>> from sympy.matrices import ImmutableMatrix
    >>> type(diag(1, cls=ImmutableMatrix))
    <class 'sympy.matrices.immutable.ImmutableMatrix'>

    See Also
    ========

    eye
    """
    from .sparse import MutableSparseMatrix

    cls = kwargs.pop('cls', None)
    if cls is None:
        from .dense import Matrix as cls

    if kwargs:
        raise ValueError('unrecognized keyword%s: %s' % (
            's' if len(kwargs) > 1 else '',
            ', '.join(kwargs.keys())))
    rows = 0
    cols = 0
    values = list(values)
    for i in range(len(values)):
        m = values[i]
        if isinstance(m, MatrixBase):
            rows += m.rows
            cols += m.cols
        elif is_sequence(m):
            m = values[i] = Matrix(m)
            rows += m.rows
            cols += m.cols
        else:
            rows += 1
            cols += 1
    res = MutableSparseMatrix.zeros(rows, cols)
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
    return cls._new(res)


def jordan_cell(eigenval, n):
    """
    Create matrix of Jordan cell kind:

    Examples
    ========

    >>> from sympy.matrices import jordan_cell
    >>> from sympy.abc import x
    >>> jordan_cell(x, 4)
    Matrix([
    [x, 1, 0, 0],
    [0, x, 1, 0],
    [0, 0, x, 1],
    [0, 0, 0, x]])
    """
    n = as_int(n)
    out = zeros(n)
    for i in range(n - 1):
        out[i, i] = eigenval
        out[i, i + 1] = S.One
    out[n - 1, n - 1] = eigenval
    return out


def hessian(f, varlist, constraints=[]):
    """Compute Hessian matrix for a function f wrt parameters in varlist
    which may be given as a sequence or a row/column vector. A list of
    constraints may optionally be given.

    Examples
    ========

    >>> from sympy import Function, hessian, pprint
    >>> from sympy.abc import x, y
    >>> f = Function('f')(x, y)
    >>> g1 = Function('g')(x, y)
    >>> g2 = x**2 + 3*y
    >>> pprint(hessian(f, (x, y), [g1, g2]))
    [                   d               d            ]
    [     0        0    --(g(x, y))     --(g(x, y))  ]
    [                   dx              dy           ]
    [                                                ]
    [     0        0        2*x              3       ]
    [                                                ]
    [                     2               2          ]
    [d                   d               d           ]
    [--(g(x, y))  2*x   ---(f(x, y))   -----(f(x, y))]
    [dx                   2            dy dx         ]
    [                   dx                           ]
    [                                                ]
    [                     2               2          ]
    [d                   d               d           ]
    [--(g(x, y))   3   -----(f(x, y))   ---(f(x, y)) ]
    [dy                dy dx              2          ]
    [                                   dy           ]

    References
    ==========

    http://en.wikipedia.org/wiki/Hessian_matrix

    See Also
    ========

    sympy.matrices.mutable.Matrix.jacobian
    wronskian
    """
    # f is the expression representing a function f, return regular matrix
    if isinstance(varlist, MatrixBase):
        if 1 not in varlist.shape:
            raise ShapeError("`varlist` must be a column or row vector.")
        if varlist.cols == 1:
            varlist = varlist.T
        varlist = varlist.tolist()[0]
    if is_sequence(varlist):
        n = len(varlist)
        if not n:
            raise ShapeError("`len(varlist)` must not be zero.")
    else:
        raise ValueError("Improper variable list in hessian function")
    if not getattr(f, 'diff'):
        # check differentiability
        raise ValueError("Function `f` (%s) is not differentiable" % f)
    m = len(constraints)
    N = m + n
    out = zeros(N)
    for k, g in enumerate(constraints):
        if not getattr(g, 'diff'):
            # check differentiability
            raise ValueError("Function `f` (%s) is not differentiable" % f)
        for i in range(n):
            out[k, i + m] = g.diff(varlist[i])
    for i in range(n):
        for j in range(i, n):
            out[i + m, j + m] = f.diff(varlist[i]).diff(varlist[j])
    for i in range(N):
        for j in range(i + 1, N):
            out[j, i] = out[i, j]
    return out


def GramSchmidt(vlist, orthonormal=False):
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
            raise ValueError(
                "GramSchmidt: vector set not linearly independent")
        out.append(tmp)
    if orthonormal:
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
    from .dense import Matrix

    for index in range(0, len(functions)):
        functions[index] = sympify(functions[index])
    n = len(functions)
    if n == 0:
        return 1
    W = Matrix(n, n, lambda i, j: functions[i].diff(var, j))
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
    from .dense import Matrix

    seqs = list(map(sympify, seqs))

    if not zero:
        f = lambda i, j: seqs[j].subs(n, n + i)
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
        raise ValueError(
            'For symmetric matrices, r must equal c, but %i != %i' % (r, c))
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
        z = int(r*c*percent // 100)
        m._mat[:z] = [S.Zero]*z
        prng.shuffle(m._mat)
    return m
