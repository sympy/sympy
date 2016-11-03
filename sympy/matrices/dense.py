from __future__ import print_function, division

import random
from types import FunctionType
from sympy import Derivative

from sympy.core.basic import Basic
from sympy.core.expr import Expr
from sympy.core.compatibility import is_sequence, as_int, range, reduce
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

from sympy.matrices.matrices import (MatrixBase, MatrixError,
    ShapeError, CommonMatrix, a2idx, classof)


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


def _iszero(x):
    """Returns True if x is zero."""
    return x.is_zero


class _InplaceMatrix(object):
    """A Matrix-like object where all operations happen
    in place.  This is meant to be used inside Matrix
    classes and should not be used directly. It has
    very few safety checks, and is optimized for efficiency."""
    __hash__ = None

    def __init__(self, rows, cols, elements, copy=False):
        if rows*cols != len(elements):
            raise ValueError("elements must have length rows*cols")
        self.rows = rows
        self.cols = cols
        if copy:
            elements = list(elements)
        self._mat = elements

    def __add__(self, other):
        if self.rows != other.rows or \
           self.cols != other.cols or \
           len(self._mat) != len(other._mat):
            raise ValueError("Wrong sized Matrices {} {}".format(self, other))
        # cache self._mat for performance
        mat = self._mat
        for i,x in enumerate(other._mat):
            mat[i] += x
        return self

    def __getitem__(self, key):
        if isinstance(key, tuple):
            i, j = key
            if isinstance(i, slice) or isinstance(j, slice):
                # if the coordinates are not slices, make them so
                # and expand the slices so they don't contain `None`
                row_slice, col_slice = self._normalize_slices(i, j)

                return _InplaceMatrixView(self, row_slice, col_slice)

            # if the key is a tuple of ints, change
            # it to an array index
            key = self._coord_to_index(i,j)
        return self._mat[key]

    def __iter__(self):
        return self._mat.__iter__()

    def __len__(self):
        return self.rows*self.cols

    def __mul__(self, other):
        if self.cols != other.rows:
            raise ValueError("Dimension mismatch when multiplying {} and {}".format(self, other))
        # allocate space for the new matrix
        new_mat_rows = self.rows
        new_mat_cols = other.cols
        new_mat = [S.Zero]*new_mat_rows*new_mat_cols

        def get_coord(i):
            return (i // new_mat_cols, i % new_mat_cols)
        # cache these methods for minor speedups
        row_indices = self._row_indices
        if hasattr(other, '_col_indices'):
            col_indices = other._col_indices
        else:
            # this allows us to multiply by a Matrix,
            cls_col_indices = self.__class__._col_indices
            col_indices = lambda i: cls_col_indices(other, i)

        # if we multiply an n x 0 with a 0 x m, the
        # expected behavior is to produce an n x m matrix of zeros
        if self.cols != 0 and other.rows != 0:
            # cache self._mat and other._mat for performance
            mat = self._mat
            other_mat = other._mat
            for i in range(len(new_mat)):
                row, col = get_coord(i)
                vec = (mat[a]*other_mat[b] for a,b \
                        in zip(row_indices(row), col_indices(col)))
                # `sum` and `Add` cannot be used here because
                # `sum([a])` returns 0+a, which mixes scalars with
                # other types of objects in block matrices, and
                # `Add(a,b)` fails for block matrices
                new_mat[i] = reduce(lambda a,b: a + b, vec)
        # make sure the operation is in-place
        self._mat[:] = new_mat
        self.rows, self.cols = new_mat_rows, new_mat_cols
        return self

    def __neg__(self):
        return self._scalar_mul(-1)

    def __repr__(self):
        return "{}({}, {}, {})".format(self.__class__.__name__, self.rows, self.cols, self._mat)

    def __setitem__(self, key, val):
        """Set item by index or tuple containing a mix
        of coordinates and slices."""
        if isinstance(key, tuple):
            i, j = key
            if isinstance(i, slice) or isinstance(j, slice):
                row_slice, col_slice = self._normalize_slices(i, j)
                inplace_mat = _InplaceMatrixView(self, row_slice, col_slice)
                if len(inplace_mat) != len(val):
                    raise ValueError("Attempting to assign {} values to {} positions".format(len(val), len(inplace_mat)))
                # when we have an appropriate view, we can set items directly
                for i,v in enumerate(val):
                    inplace_mat[i] = v
                return
            else:
                self._set((i, j), val)
                return
        self._set(key, val)

    def _col_indices(self, i):
        """Return an interable that will give the indices of
        the elements in column i"""
        if i < 0 or i >= self.cols:
            raise ValueError("i must be a valid column")
        return range(i, len(self._mat), self.cols)

    def _coord_to_index(self, i, j):
        """Return the index in _mat corresponding
        to the (i,j) position in the matrix. """
        return i*self.cols + j

    def _elementwise_mul(self, other):
        if self.rows != other.rows or \
           self.cols != other.cols or \
           len(self._mat) != len(other._mat):
            raise ValueError("Wrong sized Matrices {} {}".format(self, other))
        # cache self._mat for performance
        mat = self._mat
        for i,x in enumerate(other._mat):
            mat[i] *= x
        return self

    def _index_to_coord(self, i):
        """Return the matrix coordinates of the element
        at position `i` in the flattened matrix."""
        return (i // self.cols, i % self.cols)

    def _normalize_slices(self, row_slice, col_slice):
        """Ensure that row_slice and col_slice don't have
        `None` in their arguments.  Any integers are converted
        to slices of length 1"""
        if not isinstance(row_slice, slice):
            row_slice = slice(row_slice, row_slice+1, None)
        row_slice = slice(*row_slice.indices(self.rows))

        if not isinstance(col_slice, slice):
            col_slice = slice(col_slice, col_slice+1, None)
        col_slice = slice(*col_slice.indices(self.cols))

        return (row_slice, col_slice)

    def _row_indices(self, i):
        """Return an interable that will give the indices of
        the elements in row i"""
        if i < 0 or i >= self.rows:
            raise ValueError("i must be a valid row")
        return range(self.cols*i, self.cols*(i+1))

    def _scalar_mul(self, other):
        # cache self._mat for performance
        mat = self._mat
        for i in range(len(self._mat)):
            mat[i] *= other
        return self

    def _scalar_rmul(self, other):
        # cache self._mat for performance
        mat = self._mat
        for i in range(len(self._mat)):
            mat[i] = other*self._mat[i]
        return self

    def _set(self, coord, val):
        """Set an element in the matrix.  If coord is a tuple,
        those coordinates will be set.  Otherwise, that index will
        be set."""
        if isinstance(coord, tuple):
            i, j = coord
            key = self._coord_to_index(i, j)
            self._mat[key] = val
            return
        self._mat[coord] = val

    def applyfunc(self, f):
        for i,x in enumerate(self._mat):
            self._mat[i] = f(x)

    def col_op(self, i, f):
        self._mat[j::self.cols] = [f(*t) for t in list(zip(self._mat[j::self.cols], list(range(self.rows))))]
        return self

    def col_swap(self, i, j):
        for k in range(0, self.rows):
            self[k, i], self[k, j] = self[k, j], self[k, i]
        return self

    def copy(self):
        return self.__class__(self.rows, self.cols, self._mat, copy=True)

    @classmethod
    def eye(cls, a, b=None):
        b = a if b == None else b
        ret = cls(a, b, [S.Zero]*(a*b))
        for i in range(min(a,b)):
            ret[i,i] = S.One
        return ret

    @classmethod
    def ones(cls, a, b=None):
        b = a if b == None else b
        return cls(a, b, [S.One]*(a*b))

    def row_op(self, i, f):
        i0 = i*self.cols
        ri = self._mat[i0: i0 + self.cols]
        self._mat[i0: i0 + self.cols] = [ f(x, j) for x, j in zip(ri, list(range(self.cols))) ]
        return self

    def row_swap(self, i, j):
        for k in range(0, self.cols):
            self[i, k], self[j, k] = self[j, k], self[i, k]
        return self

    @classmethod
    def zeros(cls, a, b=None):
        b = a if b == None else b
        return cls(a, b, [S.Zero]*(a*b))

    def zip_row_op(self, i, k, f):
        i0 = i*self.cols
        k0 = k*self.cols

        ri = self._mat[i0: i0 + self.cols]
        rk = self._mat[k0: k0 + self.cols]

        self._mat[i0: i0 + self.cols] = [ f(x, y) for x, y in zip(ri, rk) ]


class _InplaceMatrixView(_InplaceMatrix):
    """A view of an InplaceMatrix.  Assigning values
    to elements of the InplaceMatrixView will affect the
    values in the InplaceMatrix.  `row_slice` and `col_slice`
    must be non-empty slices without `None` for any argument."""
    def __init__(self, inplace_mat, row_slice, col_slice):
        if None in (row_slice.start, row_slice.stop, row_slice.step, \
                    col_slice.start, col_slice.stop, col_slice.step):
            raise ValueError("slice objects cannot contain None")

        self.row_slice = row_slice
        self.col_slice = col_slice
        self.inplace_mat = inplace_mat

        self._mat = inplace_mat._mat
        self.rows = (row_slice.stop - row_slice.start) // row_slice.step
        self.cols = (col_slice.stop - col_slice.start) // col_slice.step

    def __getitem__(self, key):
        if isinstance(key, tuple):
            i, j = key
            if isinstance(i, slice) or isinstance(j, slice):
                raise NotImplementedError("Slices of views are not yet implemented")
            key = self._coord_to_index(i,j)

        # adjust the key based on the offsets of the slice
        i, j = self._index_to_coord(key)
        key = self._coord_to_true_index(i, j)

        return self._mat[key]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def _coord_to_true_index(self, i, j):
        """Returns the index in self.inplace_mat._mat corresponding
        to the (i,j) position in the View."""
        true_i = self.row_slice.start + i*self.row_slice.step
        true_j = self.col_slice.start + j*self.col_slice.step
        return self.inplace_mat._coord_to_index(true_i, true_j)

    def _set(self, coord, val):
        # we can easily convert coordinates to indices
        # using _coord_to_true_index, so if we sent an
        # index, convert it to a coord and work from there.
        if not isinstance(coord, tuple):
            coord = self._index_to_coord(coord)
        i, j = coord
        key = self._coord_to_true_index(i, j)
        self._mat[key] = val

    def copy(self):
        raise NotImplementedError("Cannot copy an InplaceMatrixView")


class DenseMatrix(CommonMatrix, MatrixBase):
    _default_inverse_method = "GE"

    _op_priority = 10.01
    _class_priority = 4

    def __eq__(self, other):
        try:
            other = sympify(other)
            if self.shape != other.shape:
                return False
            if isinstance(other, Matrix):
                return self._mat == other._mat
            elif isinstance(other, MatrixBase):
                return self._mat == Matrix(other)._mat
        except AttributeError:
            return False

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
        if is_sequence(key):
            try:
                i, j = key
            except ValueError:
                raise ValueError('Invalid Matrix Index {}'.format(key))
            try:
                i, j = self.key2ij(key)
                return self._mat[i*self.cols + j]
            except (TypeError, IndexError):
                if (isinstance(i, Expr) and not i.is_number) or (isinstance(j, Expr) and not j.is_number):
                    if ((j < 0) is True) or ((j >= self.shape[1]) is True) or\
                       ((i < 0) is True) or ((i >= self.shape[0]) is True):
                        raise ValueError("index out of boundary")
                    from sympy.matrices.expressions.matexpr import MatrixElement
                    return MatrixElement(self, i, j)

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

    def _eval_add(self, other):
        # we assume both arguments are dense matrices since
        # sparse matrices have a higher priority
        a = _InplaceMatrix(self.rows, self.cols, self._mat, copy=True)
        b = _InplaceMatrix(other.rows, other.cols, other._mat)
        a.__add__(b)
        return classof(self, other)._new(a.rows, a.cols, a._mat, copy=False)

    def _eval_col_insert(self, col, other):
        a = _InplaceMatrix(self.rows, self.cols, self._mat, copy=True)
        for i in reversed(range(self.rows)):
            pos = col + i*self.cols
            a._mat[pos:pos] = other[i,:]
        a.cols += other.cols
        return self._new(a.rows, a.cols, a._mat, copy=False)

    def _eval_col_join(self, other):
        a = _InplaceMatrix.zeros(self.rows + other.rows, self.cols)
        a[:self.rows, :] = self
        a[self.rows:, :] = other
        return classof(self, other)._new(a.rows, a.cols, a._mat, copy=False)

    @classmethod
    def _eval_create_with_diag(cls, *entries):
        return diag(*entries)

    def _eval_determinant(self):
        return self.det()

    def _eval_diff(self, *args, **kwargs):
        if kwargs.pop("evaluate", True):
            return self.diff(*args)
        else:
            return Derivative(self, *args, **kwargs)

    def _eval_elementwise_mul(self, other):
        a = _InplaceMatrix(self.rows, self.cols, self._mat, copy=True)
        b = _InplaceMatrix(other.rows, other.cols, other._mat)
        a._elementwise_mul(b)
        return classof(self, other)._new(a.rows, a.cols, a._mat, copy=False)

    def _eval_extract(self, rowsList, colsList):
        mat = self._mat
        cols = self.cols
        indices = (i*cols + j for i in rowsList for j in colsList)
        return self._new(len(rowsList), len(colsList), list(mat[i] for i in indices), copy=False)

    def _eval_has(self, *patterns):
        return any(a.has(*patterns) for a in self._mat)

    def _eval_integral_pow(self, n):
        # n >= 0 is an integer, always
        a = _InplaceMatrix(self.rows, self.cols, self._mat, copy=True)
        b = _InplaceMatrix.eye(self.rows)
        # use iterated squaring to compute the power
        b = MatrixBase._exp_by_squaring(b, a, n)
        return self._new(b.rows, b.cols, b._mat, copy=False)

    def _eval_inverse(self, method=None, **kwargs):
        return self.inv(method, **kwargs)

    def _eval_matrix_mul(self, other):
        a = _InplaceMatrix(self.rows, self.cols, self._mat, copy=True)
        b = _InplaceMatrix(other.rows, other.cols, other._mat)
        a.__mul__(b)
        return classof(self, other)._new(a.rows, a.cols, a._mat, copy=False)

    def _eval_row_insert(self, row, other):
        a = _InplaceMatrix(self.rows, self.cols, self._mat, copy=True)
        pos = row*self.cols
        a._mat[pos:pos] = other
        a.rows += other.rows
        return self._new(a.rows, a.cols, a._mat, copy=False)

    def _eval_row_join(self, other):
        a = _InplaceMatrix.zeros(self.rows, self.cols + other.cols)
        a[:, :self.cols] = self
        a[:, self.cols:] = other
        return classof(self, other)._new(a.rows, a.cols, a._mat, copy=False)

    def _eval_scalar_mul(self, other):
        a = _InplaceMatrix(self.rows, self.cols, self._mat, copy=True)
        a._scalar_mul(other)
        return self._new(a.rows, a.cols, a._mat, copy=False)

    def _eval_scalar_rmul(self, other):
        a = _InplaceMatrix(self.rows, self.cols, self._mat, copy=True)
        a._scalar_rmul(other)
        return self._new(a.rows, a.cols, a._mat, copy=False)

    def _eval_tolist(self):
        cols, mat = self.cols, self._mat
        return [mat[i*cols : (i + 1)*cols] for i in range(self.rows)]

    def _eval_transpose(self):
        a = []
        for i in range(self.cols):
            a.extend(self._mat[i::self.cols])
        return self._new(self.cols, self.rows, a)

    @classmethod
    def _eye(cls, rows, cols):
        mat = [cls._sympify(0)]*rows*cols
        mat[::(cols + 1)] = [cls._sympify(1)]*rows
        return cls._new(rows, cols, mat)

    @classmethod
    def _zeros(cls, rows, cols):
        return cls._new(rows, cols, [cls._sympify(0)]*rows*cols)

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

    def as_immutable(self):
        """Returns an Immutable version of this Matrix
        """
        from .immutable import ImmutableMatrix as cls
        if self.rows and self.cols:
            return cls._new(self.tolist())
        return cls._new(self.rows, self.cols, [])

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


class MutableDenseMatrix(DenseMatrix):
    def __new__(cls, *args, **kwargs):
        # if the copy flag was set to False, the input was rows, cols, [list]
        # and we should not create a copy of the list.
        if not kwargs.get('copy', True):
            if len(args) != 3:
                raise TypeError("'copy' requires a matrix be initialized as rows,cols,[list]")
            rows, cols, flat_list = args
        else:
            rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)
        self = object.__new__(cls)
        self.rows = rows
        self.cols = cols
        self._mat = list(flat_list)  # create a shallow copy
        return self

    @classmethod
    def _new(cls, *args, **kwargs):
        return cls(*args, **kwargs)

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

    def as_mutable(self):
        return self.copy()

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

    def fill(self, value):
        """Fill the matrix with the scalar value.

        See Also
        ========

        zeros
        ones
        """
        self._mat = [value]*len(self)

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


Matrix = MutableMatrix = MutableDenseMatrix

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


@doctest_depends_on(modules=('numpy',))
def symarray(prefix, shape, **kwargs):  # pragma: no cover
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

    \*\*kwargs : dict
      keyword arguments passed on to Symbol

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

    For setting assumptions of the underlying Symbols:

    >>> [s.is_real for s in symarray('a', 2, real=True)]
    [True, True]
    """
    from numpy import empty, ndindex
    arr = empty(shape, dtype=object)
    for index in ndindex(shape):
        arr[index] = Symbol('%s_%s' % (prefix, '_'.join(map(str, index))),
                            **kwargs)
    return arr


###############
# Functions
###############

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


def randMatrix(r, c=None, min=0, max=99, seed=None, symmetric=False,
               percent=100, prng=None):
    """Create random matrix with dimensions ``r`` x ``c``. If ``c`` is omitted
    the matrix will be square. If ``symmetric`` is True the matrix must be
    square. If ``percent`` is less than 100 then only approximately the given
    percentage of elements will be non-zero.

    The pseudo-random number generator used to generate matrix is chosen in the
    following way.

    * If ``prng`` is supplied, it will be used as random number generator.
      It should be an instance of :class:`random.Random`, or at least have
      ``randint`` and ``shuffle`` methods with same signatures.
    * if ``prng`` is not supplied but ``seed`` is supplied, then new
      :class:`random.Random` with given ``seed`` will be created;
    * otherwise, a new :class:`random.Random` with default seed will be used.

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
    # Note that ``Random()`` is equivalent to ``Random(None)``
    prng = prng or random.Random(seed)
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
