# TODO: interpret list as vectors (for multiplication)

from __future__ import division

from mptypes import (convert_lossless, absmax, mpf, mpc, rand, inf)
from functions import nthroot, sqrt

rowsep = '\n'
colsep = '  '

class matrix(object):
    """
    Numerical matrix.

    Specify the dimensions or the data as a nested list.
    Elements default to zero.
    Use a flat list to create a column vector easily.

    By default, only mpf is used to store the data. You can specify another type
    using force_type=type. It's possible to specify None.
    Make sure force_type(force_type()) is fast.
    """

    def __init__(self, *args, **kwargs):
        self.__data = {}
        # LU decompostion cache, this is useful when solving the same system
        # multiple times, when calculating the inverse and when calculating the
        # determinant
        self._LU = None
        if 'force_type' in kwargs:
            self.force_type = kwargs['force_type']
        else:
            self.force_type = convert_lossless
        if isinstance(args[0], (list, tuple)):
            if isinstance(args[0][0], (list, tuple)):
                # interpret nested list as matrix
                A = args[0]
                self.__rows = len(A)
                self.__cols = len(A[0])
                for i, row in enumerate(A):
                    for j, a in enumerate(row):
                        self[i, j] = a
            else:
                # interpret list as row vector
                v = args[0]
                self.__rows = len(v)
                self.__cols = 1
                for i, e in enumerate(v):
                    self[i, 0] = e
        elif isinstance(args[0], int):
            # create empty matrix of given dimensions
            if len(args) == 1:
                self.__rows = self.__cols = args[0]
            else:
                assert isinstance(args[1], int), 'expected int'
                self.__rows = args[0]
                self.__cols = args[1]
        elif isinstance(args[0], matrix):
            A = args[0].copy()
            self.__data = A._matrix__data
            self.__rows = A._matrix__rows
            self.__cols = A._matrix__cols
            # only copy force_type when not specified
            if not 'force_type' in kwargs:
                self.force_type = A.force_type
            elif self.force_type:
                # apply specified force_type
                for i in xrange(A.__rows):
                    for j in xrange(A.__cols):
                        A[i,j] = self.force_type(A[i,j])
        elif hasattr(args[0], 'tolist'):
            A = matrix(args[0].tolist())
            self.__data = A._matrix__data
            self.__rows = A._matrix__rows
            self.__cols = A._matrix__cols
        else:
            raise TypeError('could not interpret given arguments')

    def __nstr__(self, n=None):
        # Build table of string representations of the elements
        res = []
        # Track per-column max lengths for pretty alignment
        maxlen = [0] * self.cols
        for i in range(self.rows):
            res.append([])
            for j in range(self.cols):
                if n:
                    string = nstr(self[i,j], n)
                else:
                    string = str(self[i,j])
                res[-1].append(string)
                maxlen[j] = max(len(string), maxlen[j])
        # Patch strings together
        for i, row in enumerate(res):
            for j, elem in enumerate(row):
                # Pad each element up to maxlen so the columns line up
                row[j] = elem.rjust(maxlen[j])
            res[i] = "[" + colsep.join(row) + "]"
        return rowsep.join(res)

    def __str__(self):
        return self.__nstr__()

    def _toliststr(self, avoid_type=False):
        """
        Create a list string from a matrix.

        If avoid_type: avoid multiple 'mpf's.
        """
        s = '['
        for i in xrange(self.__rows):
            s += '['
            for j in xrange(self.__cols):
                if not avoid_type or not isinstance(self[i,j], mpf):
                    a = repr(self[i,j])
                else:  #TODO: for mpc too
                    a = "'" + str(self[i,j]) + "'"
                s += a + ', '
            s = s[:-2]
            s += '],\n '
        s = s[:-3]
        s += ']'
        return s

    def tolist(self):
        """
        Convert the matrix to a nested list.
        """
        return eval(self._toliststr())

    def __repr__(self):
        s = 'matrix(\n'
        s += self._toliststr(avoid_type=True) + ')'
        return s

    def __getitem__(self, key):
        if isinstance(key, int):
            # only sufficent for vectors
            if self.__rows == 1:
                key = (0, key)
            elif self.__cols == 1:
                key = (key, 0)
            else:
                raise IndexError('insufficient indices for matrix')
        if key in self.__data:
            return self.__data[key]
        else:
            if key[0] >= self.__rows or key[1] >= self.__cols:
                raise IndexError('matrix index out of range')
            if self.force_type:
                return self.force_type(0)
            else:
                return 0

    def __setitem__(self, key, value):
        if isinstance(key, int):
            # only sufficent for vectors
            if self.__rows == 1:
                key = (0, key)
            elif self.__cols == 1:
                key = (key, 0)
            else:
                raise IndexError('insufficient indices for matrix')
        if key[0] >= self.__rows or key[1] >= self.__cols:
            raise IndexError('matrix index out of range')
        if self.force_type: # and not isinstance(value, self.force_type):
            value = self.force_type(value)
        if value: # only store non-zeros
            self.__data[key] = value
        elif key in self.__data:
            del self.__data[key]
        # TODO: maybe do this better, if the performance impact is significant
        if self._LU:
            self._LU = None

    def __iter__(self):
        for i in xrange(self.__rows):
            for j in xrange(self.__cols):
                yield self[i,j]

    def __mul__(self, other):
        if isinstance(other, matrix):
            # dot multiplication  TODO: use Strassen's method?
            if self.__cols != other.__rows:
                raise ValueError('dimensions not compatible for multiplication')
            new = matrix(self.__rows, other.__cols)
            for i in xrange(self.__rows):
                for j in xrange(other.__cols):
                    s = 0
                    for k in xrange(other.__rows):
                        s += self[i,k] * other[k,j]
                    new[i, j] = s
            return new
        else:
            # try scalar multiplication
            new = matrix(self.__rows, self.__cols)
            for i in xrange(self.__rows):
                for j in xrange(self.__cols):
                    new[i, j] = other * self[i, j]
            return new

    def __rmul__(self, other):
        # assume other is scalar and thus commutative
        assert not isinstance(other, matrix)
        return self.__mul__(other)

    def __pow__(self, other):
        # avoid cyclic import problems
        from linalg import inverse
        if not isinstance(other, int):
            raise ValueError('only integer exponents are supported')
        if not self.__rows == self.__cols:
            raise ValueError('only powers of square matrices are defined')
        n = other
        if n == 0:
            return eye(self.__rows)
        if n < 0:
            n = -n
            neg = True
        else:
            neg = False
        i = n
        y = 1
        z = self.copy()
        while i != 0:
            if i % 2 == 1:
                y = y * z
            z = z*z
            i = i // 2
        if neg:
            y = inverse(y)
        return y

    def __div__(self, other):
        # assume other is scalar and do element-wise divison
        assert not isinstance(other, matrix)
        new = matrix(self.__rows, self.__cols)
        for i in xrange(self.__rows):
            for j in xrange(self.__cols):
                new[i,j] = self[i,j] / other
        return new

    __truediv__ = __div__

    def __add__(self, other):
        if isinstance(other, matrix):
            if not (self.__rows == other.__rows and self.__cols == other.__cols):
                raise ValueError('incompatible dimensions for addition')
            new = matrix(self.__rows, self.__cols)
            for i in xrange(self.__rows):
                for j in xrange(self.__cols):
                    new[i,j] = self[i,j] + other[i,j]
            return new
        else:
            # assume other is scalar and add element-wise
            new = matrix(self.__rows, self.__cols)
            for i in xrange(self.__rows):
                for j in xrange(self.__cols):
                    new[i,j] += self[i,j] + other
            return new

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, matrix) and not (self.__rows == other.__rows
                                              and self.__cols == other.__cols):
            raise ValueError('incompatible dimensions for substraction')
        return self.__add__(other * (-1))

    def __neg__(self):
        return (-1) * self

    def __rsub__(self, other):
        return -self + other

    def __eq__(self, other):
        return self.__rows == other.__rows and self.__cols == other.__cols \
               and self.__data == other.__data

    def __len__(self):
        if self.rows == 1:
            return self.cols
        elif self.cols == 1:
            return self.rows
        else:
            return self.rows # do it like numpy

    def __getrows(self):
        return self.__rows

    def __setrows(self, value):
        for key in self.__data.copy().iterkeys():
            if key[0] >= value:
                del self.__data[key]
        self.__rows = value

    rows = property(__getrows, __setrows, doc='number of rows')

    def __getcols(self):
        return self.__cols

    def __setcols(self, value):
        for key in self.__data.copy().iterkeys():
            if key[1] >= value:
                del self.__data[key]
        self.__cols = value

    cols = property(__getcols, __setcols, doc='number of columns')

    def transpose(self):
        new = matrix(self.__cols, self.__rows)
        for i in xrange(self.__rows):
            for j in xrange(self.__cols):
                new[j,i] = self[i,j]
        return new

    T = property(transpose)

    def copy(self):
        new = matrix(self.__rows, self.__cols, force_type=self.force_type)
        new.__data = self.__data.copy()
        return new

    __copy__ = copy

def eye(n, **kwargs):
    """
    Create square identity matrix n x n.
    """
    A = matrix(n, **kwargs)
    for i in xrange(n):
        A[i,i] = 1
    return A

def diag(diagonal, **kwargs):
    """
    Create square diagonal matrix using given list.
    """
    A = matrix(len(diagonal), **kwargs)
    for i in xrange(len(diagonal)):
        A[i,i] = diagonal[i]
    return A

def zeros(*args, **kwargs):
    """
    Create matrix m x n filled with zeros.
    One given dimension will create square matrix n x n.
    """
    if len(args) == 1:
        m = n = args[0]
    elif len(args) == 2:
        m = args[0]
        n = args[1]
    else:
        raise TypeError('zeros expected at most 2 arguments, got %i' % len(args))
    A = matrix(m, n, **kwargs)
    for i in xrange(m):
        for j in xrange(n):
            A[i,j] = 0
    return A

def ones(*args, **kwargs):
    """
    Create matrix m x n filled with ones.
    One given dimension will create square matrix n x n.
    """
    if len(args) == 1:
        m = n = args[0]
    elif len(args) == 2:
        m = args[0]
        n = args[1]
    else:
        raise TypeError('ones expected at most 2 arguments, got %i' % len(args))
    A = matrix(m, n, **kwargs)
    for i in xrange(m):
        for j in xrange(n):
            A[i,j] = 1
    return A

def randmatrix(m, n=None, min=0, max=1, **kwargs):
    """
    Create a random m x n matrix.

    All values are >= min and <max.
    n defaults to m.
    """
    if not n:
        n = m
    A = matrix(m, n, **kwargs)
    for i in xrange(m):
        for j in xrange(n):
            A[i,j] = rand() * (max - min) + min
    return A

def swap_row(A, i, j):
    """
    Swap row i with row j.
    """
    if i == j:
        return
    if isinstance(A, matrix):
        for k in xrange(A.cols):
            A[i,k], A[j,k] = A[j,k], A[i,k]
    elif isinstance(A, list):
        A[i], A[j] = A[j], A[i]
    else:
        raise TypeError('could not interpret type')

def extend(A, b):
    """
    Extend matrix A with column b and return result.
    """
    assert isinstance(A, matrix)
    assert A.rows == len(b)
    A = A.copy()
    A.cols += 1
    for i in xrange(A.rows):
        A[i, A.cols-1] = b[i]
    return A

def mnorm_1(A):
    """
    Calculate the 1-norm (maximal column sum) of a matrix.
    """
    assert isinstance(A, matrix)
    m, n = A.rows, A.cols
    return max((sum((absmax(A[i,j]) for i in xrange(m))) for j in xrange(n)))

def mnorm_oo(A):
    """
    Calculate the oo-norm (maximal row sum) of a matrix.
    """
    assert isinstance(A, matrix)
    m, n = A.rows, A.cols
    return max((sum((absmax(A[i,j]) for j in xrange(n))) for i in xrange(m)))

def mnorm_F(A):
    """
    Calculate the Frobenius norm (root of sum of squares) of a matrix.

    It can be useful for estimating the spectral norm (which is difficult to
    calculate):
    mnorm_F(A) / sqrt(rank(A)) <= mnorm_2(A) <= mnorm_F(A)
    """
    return norm_p(A, 2)

def norm_p(x, p=2):
    """
    Calculate the p-norm of a vector.
    0 < p <= oo

    Note: you may want to use float('inf') or mpmath's equivalent to specify oo.
    """
    if p == inf:
        return max((absmax(i) for i in x))
    elif p > 1:
        return nthroot(sum((abs(i)**p for i in x)), p)
    elif p == 1:
        return sum((abs(i) for i in x))
    else:
        raise ValueError('p has to be an integer greater than 0')

