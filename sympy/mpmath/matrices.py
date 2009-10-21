# TODO: interpret list as vectors (for multiplication)

from __future__ import division

from mptypes import mpmathify, absmax, mpf, mpc, rand, inf, nstr, fsum, fdot
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

    Creating matrices
    -----------------

    Matrices in mpmath are implemented using dictionaries. Only non-zero values
    are stored, so it is cheap to represent sparse matrices.

    The most basic way to create one is to use the ``matrix`` class directly.
    You can create an empty matrix specifying the dimensions:

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> matrix(2)
        matrix(
        [['0.0', '0.0'],
         ['0.0', '0.0']])
        >>> matrix(2, 3)
        matrix(
        [['0.0', '0.0', '0.0'],
         ['0.0', '0.0', '0.0']])

    Calling ``matrix`` with one dimension will create a square matrix.

    To access the dimensions of a matrix, use the ``rows`` or ``cols`` keyword:

        >>> A = matrix(3, 2)
        >>> A
        matrix(
        [['0.0', '0.0'],
         ['0.0', '0.0'],
         ['0.0', '0.0']])
        >>> A.rows
        3
        >>> A.cols
        2

    You can also change the dimension of an existing matrix. This will set the
    new elements to 0. If the new dimension is smaller than before, the
    concerning elements are discarded:

        >>> A.rows = 2
        >>> A
        matrix(
        [['0.0', '0.0'],
         ['0.0', '0.0']])

    Internally ``mpmathify`` is used every time an element is set. This
    is done using the syntax A[row,column], counting from 0:

        >>> A = matrix(2)
        >>> A[1,1] = 1 + 1j
        >>> A
        matrix(
        [['0.0', '0.0'],
         ['0.0', '(1.0 + 1.0j)']])

    You can use the keyword ``force_type`` to change the function which is
    called on every new element:

        >>> matrix(2, 5, force_type=int)
        matrix(
        [[0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0]])

    A more comfortable way to create a matrix lets you use nested lists:

        >>> matrix([[1, 2], [3, 4]])
        matrix(
        [['1.0', '2.0'],
         ['3.0', '4.0']])

    If you want to preserve the type of the elements you can use
    ``force_type=None``:

        >>> matrix([[1, 2.5], [1j, mpf(2)]], force_type=None)
        matrix(
        [[1, 2.5],
         [1j, '2.0']])

    Convenient advanced functions are available for creating various standard
    matrices, see ``zeros``, ``ones``, ``diag``, ``eye``, ``randmatrix`` and
    ``hilbert``.

    Vectors
    .......

    Vectors may also be represented by the ``matrix`` class (with rows = 1 or cols = 1).
    For vectors there are some things which make life easier. A column vector can
    be created using a flat list, a row vectors using an almost flat nested list::

        >>> matrix([1, 2, 3])
        matrix(
        [['1.0'],
         ['2.0'],
         ['3.0']])
        >>> matrix([[1, 2, 3]])
        matrix(
        [['1.0', '2.0', '3.0']])

    Optionally vectors can be accessed like lists, using only a single index::

        >>> x = matrix([1, 2, 3])
        >>> x[1]
        mpf('2.0')
        >>> x[1,0]
        mpf('2.0')

    Other
    .....

    Like you probably expected, matrices can be printed::

        >>> print randmatrix(3) # doctest:+SKIP
        [ 0.782963853573023  0.802057689719883  0.427895717335467]
        [0.0541876859348597  0.708243266653103  0.615134039977379]
        [ 0.856151514955773  0.544759264818486  0.686210904770947]

    Use ``nstr`` or ``nprint`` to specify the number of digits to print::

        >>> nprint(randmatrix(5), 3) # doctest:+SKIP
        [2.07e-1  1.66e-1  5.06e-1  1.89e-1  8.29e-1]
        [6.62e-1  6.55e-1  4.47e-1  4.82e-1  2.06e-2]
        [4.33e-1  7.75e-1  6.93e-2  2.86e-1  5.71e-1]
        [1.01e-1  2.53e-1  6.13e-1  3.32e-1  2.59e-1]
        [1.56e-1  7.27e-2  6.05e-1  6.67e-2  2.79e-1]

    As matrices are mutable, you will need to copy them sometimes::

        >>> A = matrix(2)
        >>> A
        matrix(
        [['0.0', '0.0'],
         ['0.0', '0.0']])
        >>> B = A.copy()
        >>> B[0,0] = 1
        >>> B
        matrix(
        [['1.0', '0.0'],
         ['0.0', '0.0']])
        >>> A
        matrix(
        [['0.0', '0.0'],
         ['0.0', '0.0']])

    Finally, it is possible to convert a matrix to a nested list. This is very useful,
    as most Python libraries involving matrices or arrays (namely NumPy or SymPy)
    support this format::

        >>> B.tolist()
        [[mpf('1.0'), mpf('0.0')], [mpf('0.0'), mpf('0.0')]]


    Matrix operations
    -----------------

    You can add and substract matrices of compatible dimensions::

        >>> A = matrix([[1, 2], [3, 4]])
        >>> B = matrix([[-2, 4], [5, 9]])
        >>> A + B
        matrix(
        [['-1.0', '6.0'],
         ['8.0', '13.0']])
        >>> A - B
        matrix(
        [['3.0', '-2.0'],
         ['-2.0', '-5.0']])
        >>> A + ones(3) # doctest:+ELLIPSIS
        Traceback (most recent call last):
          ...
        ValueError: incompatible dimensions for addition

    It is possible to multiply or add matrices and scalars. In the latter case the
    operation will be done element-wise::

        >>> A * 2
        matrix(
        [['2.0', '4.0'],
         ['6.0', '8.0']])
        >>> A / 4
        matrix(
        [['0.25', '0.5'],
         ['0.75', '1.0']])
        >>> A - 1
        matrix(
        [['0.0', '1.0'],
         ['2.0', '3.0']])

    Of course you can perform matrix multiplication, if the dimensions are
    compatible::

        >>> A * B
        matrix(
        [['8.0', '22.0'],
         ['14.0', '48.0']])
        >>> matrix([[1, 2, 3]]) * matrix([[-6], [7], [-2]])
        matrix(
        [['2.0']])

    You can raise powers of square matrices::

        >>> A**2
        matrix(
        [['7.0', '10.0'],
         ['15.0', '22.0']])

    Negative powers will calculate the inverse::

        >>> A**-1
        matrix(
        [['-2.0', '1.0'],
         ['1.5', '-0.5']])
        >>> A * A**-1
        matrix(
        [['1.0', '1.0842021724855e-19'],
         ['-2.16840434497101e-19', '1.0']])

    Matrix transposition is straightforward::

        >>> A = ones(2, 3)
        >>> A
        matrix(
        [['1.0', '1.0', '1.0'],
         ['1.0', '1.0', '1.0']])
        >>> A.T
        matrix(
        [['1.0', '1.0'],
         ['1.0', '1.0'],
         ['1.0', '1.0']])

    Norms
    .....

    Sometimes you need to know how "large" a matrix or vector is. Due to their
    multidimensional nature it's not possible to compare them, but there are
    several functions to map a matrix or a vector to a positive real number, the
    so called norms.

    For vectors the p-norm is intended, usually the 1-, the 2- and the oo-norm are
    used.

        >>> x = matrix([-10, 2, 100])
        >>> norm(x, 1)
        mpf('112.0')
        >>> norm(x, 2)
        mpf('100.5186549850325')
        >>> norm(x, inf)
        mpf('100.0')

    Please note that the 2-norm is the most used one, though it is more expensive
    to calculate than the 1- or oo-norm.

    It is possible to generalize some vector norms to matrix norm::

        >>> A = matrix([[1, -1000], [100, 50]])
        >>> mnorm(A, 1)
        mpf('1050.0')
        >>> mnorm(A, inf)
        mpf('1001.0')
        >>> mnorm(A, 'F')
        mpf('1006.2310867787777')

    The last norm (the "Frobenius-norm") is an approximation for the 2-norm, which
    is hard to calculate and not available. The Frobenius-norm lacks some
    mathematical properties you might expect from a norm.
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
            self.force_type = mpmathify
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
                if not avoid_type or not isinstance(self[i,j], (mpf, mpc)):
                    a = repr(self[i,j])
                else:
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
        if type(key) is int:
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
        if type(key) is int:
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
                    new[i, j] = fdot((self[i,k], other[k,j])
                                     for k in xrange(other.__rows))
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

    def column(self, n):
        m = matrix(self.rows, 1)
        for i in range(self.rows):
            m[i] = self[i,n]
        return m

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

    Example:
    >>> from mpmath import diag
    >>> diag([1, 2, 3])
    matrix(
    [['1.0', '0.0', '0.0'],
     ['0.0', '2.0', '0.0'],
     ['0.0', '0.0', '3.0']])

    """
    A = matrix(len(diagonal), **kwargs)
    for i in xrange(len(diagonal)):
        A[i,i] = diagonal[i]
    return A

def zeros(*args, **kwargs):
    """
    Create matrix m x n filled with zeros.
    One given dimension will create square matrix n x n.

    Example:
    >>> from mpmath import zeros
    >>> zeros(2)
    matrix(
    [['0.0', '0.0'],
     ['0.0', '0.0']])

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

    Example:
    >>> from mpmath import ones
    >>> ones(2)
    matrix(
    [['1.0', '1.0'],
     ['1.0', '1.0']])

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

def hilbert(m, n=None):
    """
    Create (pseudo) hilbert matrix m x n.
    One given dimension will create hilbert matrix n x n.

    The matrix is very ill-conditioned and symmetric, positive definit if
    square.
    """
    if n is None:
        n = m
    A = matrix(m, n)
    for i in xrange(m):
        for j in xrange(n):
            A[i,j] = 1./ (i + j + 1)
    return A

def randmatrix(m, n=None, min=0, max=1, **kwargs):
    """
    Create a random m x n matrix.

    All values are >= min and <max.
    n defaults to m.

    Example:
    >>> from mpmath import randmatrix
    >>> randmatrix(2) # doctest:+SKIP
    matrix(
    [['0.53491598236191806', '0.57195669543302752'],
     ['0.85589992269513615', '0.82444367501382143']])

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

def norm(x, p=2):
    r"""
    Gives the entrywise `p`-norm of an iterable *x*, i.e. the vector norm
    `\left(\sum_k |x_k|^p\right)^{1/p}`, for any given `1 \le p \le \infty`.

    Special cases:

    If *x* is not iterable, this just returns ``absmax(x)``.

    ``p=1`` gives the sum of absolute values.

    ``p=2`` is the standard Euclidean vector norm.

    ``p=inf`` gives the magnitude of the largest element.

    For *x* a matrix, ``p=2`` is the Frobenius norm.
    For operator matrix norms, use :func:`mnorm` instead.

    You can use the string 'inf' as well as float('inf') or mpf('inf')
    to specify the infinity norm.

    **Examples**

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> x = matrix([-10, 2, 100])
        >>> norm(x, 1)
        mpf('112.0')
        >>> norm(x, 2)
        mpf('100.5186549850325')
        >>> norm(x, inf)
        mpf('100.0')

    """
    try:
        iter(x)
    except TypeError:
        return absmax(x)
    if type(p) is not int:
        p = mpmathify(p)
    if p == inf:
        return max(absmax(i) for i in x)
    elif p == 1:
        return fsum(x, absolute=1)
    elif p == 2:
        return sqrt(fsum(x, absolute=1, squared=1))
    elif p > 1:
        return nthroot(fsum(abs(i)**p for i in x), p)
    else:
        raise ValueError('p has to be >= 1')

def mnorm(A, p=1):
    r"""
    Gives the matrix (operator) `p`-norm of A. Currently ``p=1`` and ``p=inf``
    are supported:

    ``p=1`` gives the 1-norm (maximal column sum)

    ``p=inf`` gives the `\infty`-norm (maximal row sum).
    You can use the string 'inf' as well as float('inf') or mpf('inf')

    ``p=2`` (not implemented) for a square matrix is the usual spectral
    matrix norm, i.e. the largest singular value.

    ``p='f'`` (or 'F', 'fro', 'Frobenius, 'frobenius') gives the
    Frobenius norm, which is the elementwise 2-norm. The Frobenius norm is an
    approximation of the spectral norm and satisfies

    .. math ::

        \frac{1}{\sqrt{\mathrm{rank}(A)}} \|A\|_F \le \|A\|_2 \le \|A\|_F

    The Frobenius norm lacks some mathematical properties that might
    be expected of a norm.

    For general elementwise `p`-norms, use :func:`norm` instead.

    **Examples**

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> A = matrix([[1, -1000], [100, 50]])
        >>> mnorm(A, 1)
        mpf('1050.0')
        >>> mnorm(A, inf)
        mpf('1001.0')
        >>> mnorm(A, 'F')
        mpf('1006.2310867787777')

    """
    A = matrix(A)
    if type(p) is not int:
        if type(p) is str and 'frobenius'.startswith(p.lower()):
            return norm(A, 2)
        p = mpmathify(p)
    m, n = A.rows, A.cols
    if p == 1:
        return max(fsum((A[i,j] for i in xrange(m)), absolute=1) for j in xrange(n))
    elif p == inf:
        return max(fsum((A[i,j] for j in xrange(n)), absolute=1) for i in xrange(m))
    else:
        raise NotImplementedError("matrix p-norm for arbitrary p")

if __name__ == '__main__':
    import doctest
    doctest.testmod()
