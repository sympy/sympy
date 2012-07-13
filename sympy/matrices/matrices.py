from sympy.core.add import Add
from sympy.core.basic import Basic, C
from sympy.core.function import count_ops
from sympy.core.power import Pow
from sympy.core.symbol import Symbol, Dummy
from sympy.core.numbers import Integer, ilcm, Rational, Float
from sympy.core.singleton import S
from sympy.core.sympify import sympify, SympifyError
from sympy.core.compatibility import is_sequence

from sympy.polys import PurePoly, roots, cancel
from sympy.simplify import simplify as _simplify, signsimp, nsimplify
from sympy.utilities.iterables import flatten
from sympy.utilities.misc import filldedent, default_sort_key
from sympy.functions.elementary.miscellaneous import sqrt, Max, Min
from sympy.printing import sstr
from sympy.functions.elementary.trigonometric import cos, sin

from sympy.core.compatibility import callable, reduce
from sympy.utilities.exceptions import SymPyDeprecationWarning
from sympy.core.decorators import call_highest_priority

import random
from types import FunctionType

class MatrixError(Exception):
    pass

class ShapeError(ValueError, MatrixError):
    """Wrong matrix shape"""
    pass

class NonSquareMatrixError(ShapeError):
    pass

def _iszero(x):
    """Returns True if x is zero."""
    return x.is_zero

class DeferredVector(Symbol):
    """A vector whose components are deferred (e.g. for use with lambdify)

    >>> from sympy import DeferredVector, lambdify
    >>> X = DeferredVector( 'X' )
    >>> X
    X
    >>> expr = (X[0] + 2, X[2] + 3)
    >>> func = lambdify( X, expr )
    >>> func( [1, 2, 3] )
    (3, 6)
    """
    def __getitem__(self,i):
        if i == -0:
            i = 0
        if i < 0:
            raise IndexError('DeferredVector index out of range')
        component_name = '%s[%d]'%(self.name,i)
        return Symbol(component_name)

    def __str__(self):
        return sstr(self)

    def __repr__(self):
        return "DeferredVector('%s')"%(self.name)


class MatrixBase(object):

    # Added just for numpy compatibility
    # TODO: investigate about __array_priority__
    __array_priority__ = 10.0

    is_Matrix = True
    _class_priority = 3

    def _sympy_(self):
        #return self.as_immutable()
        raise SympifyError('Matrix cannot be sympified')

    @classmethod
    def _handle_creation_inputs(cls, *args, **kwargs):
        """
        Matrix can be constructed with values or a rule.

        >>> from sympy import Matrix, I
        >>> Matrix( ((1,2+I), (3,4)) ) #doctest:+NORMALIZE_WHITESPACE
        [1, 2 + I]
        [3,     4]
        >>> Matrix(2, 2, lambda i,j: (i+1)*j ) #doctest:+NORMALIZE_WHITESPACE
        [0, 1]
        [0, 2]

        """
        # Matrix(Matrix(...))
        if len(args) == 1 and isinstance(args[0], MatrixBase):
            return args[0].rows, args[0].cols, args[0].mat

        # Matrix(MatrixSymbol('X', 2,2))
        if len(args) == 1 and isinstance(args[0], Basic) and args[0].is_Matrix:
            return args[0].rows, args[0].cols, args[0].as_explicit().mat

        # Matrix(2, 2, lambda i,j: i+j)
        if len(args) == 3 and callable(args[2]):
            operation = args[2]
            rows = args[0]
            cols = args[1]
            mat = []
            for i in range(rows):
                for j in range(cols):
                    mat.append(sympify(operation(i, j)))

        # Matrix(2, 2, [1,2,3,4])
        elif len(args)==3 and is_sequence(args[2]):
            rows = args[0]
            cols = args[1]
            mat  = args[2]
            if len(mat) != rows*cols:
                raise ValueError('List length should be equal to rows*columns')
            mat = map(lambda i: sympify(i), mat)

        # Matrix(numpy.ones((2,2)))
        elif len(args) == 1:
            in_mat = args[0]
            if hasattr(in_mat, "__array__"): #pragma: no cover
                # NumPy array or matrix or some other object that implements
                # __array__. So let's first use this method to get a
                # numpy.array() and then make a python list out of it.
                arr = in_mat.__array__()
                if len(arr.shape) == 2:
                    rows, cols = arr.shape[0], arr.shape[1]
                    mat = map(lambda i: sympify(i), arr.ravel())
                    return rows, cols, mat
                elif len(arr.shape) == 1:
                    rows, cols = 1, arr.shape[0]
                    mat = [0]*cols
                    for i in xrange(len(arr)):
                        mat[i] = sympify(arr[i])
                    return rows, cols, mat
                else:
                    raise NotImplementedError("SymPy supports just 1D and 2D matrices")
            elif not is_sequence(in_mat, include=MatrixBase):
                raise TypeError("Matrix constructor doesn't accept %s as input"
                        % str(type(in_mat)))
            in_mat = []
            for row in args[0]:
                if isinstance(row, MatrixBase):
                    in_mat.extend(row.tolist())
                else:
                    in_mat.append(row)
            rows = len(in_mat)
            if len(in_mat) != 0:
                if not is_sequence(in_mat[0]):
                    cols = 1
                    mat = map(lambda i: sympify(i), in_mat)
                    return rows, cols, mat
                cols = len(in_mat[0])
            else:
                cols = 0
            mat = []
            for j in xrange(rows):
                if len(in_mat[j]) != cols:
                    raise ValueError("Input %s inconsistant to form a Matrix." %
                        args)
                for i in xrange(cols):
                    mat.append(sympify(in_mat[j][i]))

        # Matrix()
        elif len(args) == 0:
            # Empty Matrix
            rows = cols = 0
            mat = []
        else:
            raise TypeError("Data type not understood")

        return rows, cols, mat

    def transpose(self):
        """
        Matrix transposition.

        >>> from sympy import Matrix, I
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        [1, 2 + I]
        [3,     4]
        >>> m.transpose() #doctest: +NORMALIZE_WHITESPACE
        [    1, 3]
        [2 + I, 4]
        >>> m.T == m.transpose()
        True

        See Also
        ========

        conjugate: By-element conjugation
        """
        a = [0]*len(self)
        for i in xrange(self.cols):
            a[i*self.rows:(i+1)*self.rows] = self.mat[i::self.cols]
        return self._new(self.cols,self.rows,a)

    T = property(transpose,None,None,"Matrix transposition.")

    def conjugate(self):
        """By-element conjugation.

        See Also
        ========

        transpose: Matrix transposition
        H: Hermite conjugation
        D: Dirac conjugation
        """
        out = self._new(self.rows,self.cols,
                lambda i,j: self[i,j].conjugate())
        return out

    C = property(conjugate,None,None,"By-element conjugation.")

    @property
    def H(self):
        """
        Hermite conjugation.

        >>> from sympy import Matrix, I
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        [1, 2 + I]
        [3,     4]
        >>> m.H #doctest: +NORMALIZE_WHITESPACE
        [    1, 3]
        [2 - I, 4]

        See Also
        ========

        conjugate: By-element conjugation
        D: Dirac conjugation
        """
        out = self.T.C
        return out

    @property
    def D(self):
        """Dirac conjugation.

        See Also
        ========

        conjugate: By-element conjugation
        H: Hermite conjugation
        """
        from sympy.physics.matrices import mgamma
        try:
            out = self.H * mgamma(0)
            return out
        # In Python 3.2, properties can only return an AttributeError, so we
        # have to catch the ShapeError; see also the commit making this change
        except ShapeError:
            raise AttributeError("Dirac conjugation not possible.")

    def __getitem__(self,key):
        """
        >>> from sympy import Matrix, I
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        [1, 2 + I]
        [3,     4]
        >>> m[1,0]
        3
        >>> m.H[1,0]
        2 - I

        """
        if type(key) is tuple:
            i, j = key
            if type(i) is slice or type(j) is slice:
                return self.submatrix(key)

            else:
                # a2idx inlined
                if type(i) is not int:
                    try:
                        i = i.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % (key,))
                # a2idx inlined
                if type(j) is not int:
                    try:
                        j = j.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % (key,))

                i, j = self.key2ij((i, j))
                if not (i>=0 and i<self.rows and j>=0 and j < self.cols):
                    raise IndexError("Index out of range: a[%s]" % (key,))
                else:
                    return self.mat[i*self.cols + j]


        else:
            # row-wise decomposition of matrix
            if type(key) is slice:
                return self.mat[key]
            else:
                k = a2idx(key)
                if k is not None:
                    return self.mat[k]
        raise IndexError("Invalid index: a[%s]" % repr(key))

    def __setitem__(self, key, value):
        raise NotImplementedError()

    def as_mutable(self):
        """
        Returns a Mutable version of this Matrix

        >>> from sympy import ImmutableMatrix
        >>> X = ImmutableMatrix([[1,2],[3,4]])
        >>> Y = X.as_mutable()
        >>> Y[1,1] = 5 # Can set values in Y
        >>> Y
        [1, 2]
        [3, 5]
        """
        return MutableMatrix(self.rows, self.cols, self.mat)

    def as_immutable(self):
        """
        Returns an Immutable version of this Matrix
        """
        from immutable_matrix import ImmutableMatrix
        return ImmutableMatrix(self.rows, self.cols, self.mat)

    def __array__(self):
        return matrix2numpy(self)

    def __len__(self):
        """
        Return the number of elements of self.

        Implemented mainly so bool(Matrix()) == False.
        """
        return self.rows * self.cols

    def tolist(self):
        """
        Return the Matrix converted in a python list.

        >>> from sympy import Matrix
        >>> m=Matrix(3, 3, range(9))
        >>> m
        [0, 1, 2]
        [3, 4, 5]
        [6, 7, 8]
        >>> m.tolist()
        [[0, 1, 2], [3, 4, 5], [6, 7, 8]]

        """
        ret = [0]*self.rows
        for i in xrange(self.rows):
            ret[i] = self.mat[i*self.cols:(i+1)*self.cols]
        return ret

    def hash(self):
        """Compute a hash every time, because the matrix elements
        could change."""
        return hash(self.__str__() )

    @property
    def shape(self):
        """The shape (dimensions) of the matrix as the 2-tuple (rows, cols)."""
        return (self.rows, self.cols)

    def __rmul__(self,a):
        if hasattr(a, "__array__") and a.shape != ():
            return matrix_multiply(a,self)
        out = self._new(self.rows, self.cols, map(lambda i: a*i,self.mat))
        return out

    def expand(self, **hints):
        """
        Expand each element of the matrix by calling `expand()`.
        """
        out = self._new(self.rows, self.cols,
                map(lambda i: i.expand(**hints), self.mat))
        return out

    def subs(self, *args, **kwargs):
        """
        Create substituted expressions for each element with `Expr.subs`.
        """
        out = self._new(self.rows, self.cols,
                map(lambda i: i.subs(*args, **kwargs), self.mat))
        return out

    def __sub__(self,a):
        return self + (-a)

    def __rsub__(self,a):
        return (-self) + a

    def __mul__(self,a):
        if hasattr(a, "__array__") and a.shape != ():
            return matrix_multiply(self,a)
        out = self._new(self.rows, self.cols, map(lambda i: i*a,self.mat))
        return out

    def __pow__(self, num):
        if not self.is_square:
            raise NonSquareMatrixError()
        if isinstance(num, int) or isinstance(num, Integer):
            n = int(num)
            if n < 0:
                return self.inv() ** -n   # A**-2 = (A**-1)**2
            a = eye(self.cols)
            s = self
            while n:
                if n%2:
                    a *= s
                    n -= 1
                if not n:
                    break
                s *= s
                n //= 2
            return a
        elif isinstance(num, Rational):
            try:
                P, D = self.diagonalize()
            except MatrixError:
                raise NotImplementedError("Implemented only for diagonalizable matrices")
            for i in range(D.rows):
                D[i, i] = D[i, i]**num
            return P * D * P.inv()
        else:
            raise NotImplementedError("Only integer and rational values are supported")

    def __add__(self,a):
        return matrix_add(self,a)

    def __radd__(self,a):
        return matrix_add(a,self)

    def __div__(self,a):
        return self * (S.One/a)

    def __truediv__(self,a):
        return self.__div__(a)

    def multiply(self,b):
        """Returns self*b

        See Also
        ========

        dot
        cross

        multiply_elementwise
        """
        return matrix_multiply(self,b)

    def add(self,b):
        """Return self+b """
        return matrix_add(self,b)

    def __neg__(self):
        return -1*self

    def equals(self, other):
        try:
            return (self.shape == other.shape and
                    all([self[i,j] == other[i,j]
                        for i in xrange(self.rows)
                        for j in xrange(self.cols)]))
        except AttributeError:
            return False

    def __eq__(self, other):
        return self.equals(other)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return super(MatrixBase, self).__hash__()

    def _format_str(self, strfunc, rowsep='\n'):
        # Handle zero dimensions:
        if self.rows == 0 or self.cols == 0:
            return '[]'
        # Build table of string representations of the elements
        res = []
        # Track per-column max lengths for pretty alignment
        maxlen = [0] * self.cols
        for i in range(self.rows):
            res.append([])
            for j in range(self.cols):
                string = strfunc(self[i,j])
                res[-1].append(string)
                maxlen[j] = max(len(string), maxlen[j])
        # Patch strings together
        for i, row in enumerate(res):
            for j, elem in enumerate(row):
                # Pad each element up to maxlen so the columns line up
                row[j] = elem.rjust(maxlen[j])
            res[i] = "[" + ", ".join(row) + "]"
        return rowsep.join(res)

    def __str__(self):
        return sstr(self)

    def __repr__(self):
        return sstr(self)

    def cholesky(self):
        """
        Returns the Cholesky Decomposition L of a Matrix A
        such that L * L.T = A

        A must be a square, symmetric, positive-definite
        and non-singular matrix

        >>> from sympy.matrices import Matrix
        >>> A = Matrix(((25,15,-5),(15,18,0),(-5,0,11)))
        >>> A.cholesky()
        [ 5, 0, 0]
        [ 3, 3, 0]
        [-1, 1, 3]
        >>> A.cholesky() * A.cholesky().T
        [25, 15, -5]
        [15, 18,  0]
        [-5,  0, 11]

        See Also
        ========

        LDLdecomposition
        LUdecomposition
        QRdecomposition
        """

        if not self.is_square:
            raise NonSquareMatrixError("Matrix must be square.")
        if not self.is_symmetric():
            raise ValueError("Matrix must be symmetric.")
        return self._cholesky()

    def _cholesky(self):
        """
        Helper function of cholesky.
        Without the error checks.
        To be used privately. """
        L = zeros(self.rows, self.rows)
        for i in xrange(self.rows):
            for j in xrange(i):
                L[i, j] = (1 / L[j, j]) * (self[i, j] - sum(L[i, k] * L[j, k]
                    for k in xrange(j)))
            L[i, i] = sqrt(self[i, i] - sum(L[i, k] ** 2
                for k in xrange(i)))
        return self._new(L)

    def LDLdecomposition(self):
        """
        Returns the LDL Decomposition (L,D) of matrix A,
        such that L * D * L.T == A
        This method eliminates the use of square root.
        Further this ensures that all the diagonal entries of L are 1.
        A must be a square, symmetric, positive-definite
        and non-singular matrix.

        >>> from sympy.matrices import Matrix, eye
        >>> A = Matrix(((25,15,-5),(15,18,0),(-5,0,11)))
        >>> L, D = A.LDLdecomposition()
        >>> L
        [   1,   0, 0]
        [ 3/5,   1, 0]
        [-1/5, 1/3, 1]
        >>> D
        [25, 0, 0]
        [ 0, 9, 0]
        [ 0, 0, 9]
        >>> L * D * L.T * A.inv() == eye(A.rows)
        True

        See Also
        ========

        cholesky
        LUdecomposition
        QRdecomposition
        """
        if not self.is_square:
            raise NonSquareMatrixError("Matrix must be square.")
        if not self.is_symmetric():
            raise ValueError("Matrix must be symmetric.")
        return self._LDLdecomposition()

    def _LDLdecomposition(self):
        """
        Helper function of LDLdecomposition.
        Without the error checks.
        To be used privately.
        """
        D = zeros(self.rows, self.rows)
        L = eye(self.rows)
        for i in xrange(self.rows):
            for j in xrange(i):
                L[i, j] = (1 / D[j, j]) * (self[i, j] - sum(
                    L[i, k] * L[j, k] * D[k, k] for k in xrange(j)))
            D[i, i] = self[i, i] - sum(L[i, k]**2 * D[k, k]
                for k in xrange(i))
        return self._new(L), self._new(D)

    def lower_triangular_solve(self, rhs):
        """
        Solves Ax = B, where A is a lower triangular matrix.

        See Also
        ========

        upper_triangular_solve
        cholesky_solve
        diagonal_solve
        LDLsolve
        LUsolve
        QRsolve
        """

        if not self.is_square:
            raise NonSquareMatrixError("Matrix must be square.")
        if rhs.rows != self.rows:
            raise ShapeError("Matrices size mismatch.")
        if not self.is_lower():
            raise ValueError("Matrix must be lower triangular.")
        return self._lower_triangular_solve(rhs)

    def _lower_triangular_solve(self, rhs):
        """
        Helper function of function lower_triangular_solve.
        Without the error checks.
        To be used privately.
        """
        X = zeros(self.rows, 1)
        for i in xrange(self.rows):
            if self[i, i] == 0:
                raise TypeError("Matrix must be non-singular.")
            X[i, 0] = (rhs[i, 0] - sum(self[i, k] * X[k, 0]
                for k in xrange(i))) / self[i, i]
        return self._new(X)

    def upper_triangular_solve(self, rhs):
        """
        Solves Ax = B, where A is an upper triangular matrix.

        See Also
        ========

        lower_triangular_solve
        cholesky_solve
        diagonal_solve
        LDLsolve
        LUsolve
        QRsolve
        """
        if not self.is_square:
            raise NonSquareMatrixError("Matrix must be square.")
        if rhs.rows != self.rows:
            raise TypeError("Matrix size mismatch.")
        if not self.is_upper():
            raise TypeError("Matrix is not upper triangular.")
        return self._upper_triangular_solve(rhs)

    def _upper_triangular_solve(self, rhs):
        """
        Helper function of function upper_triangular_solve.
        Without the error checks, to be used privately. """
        X = zeros(self.rows, 1)
        for i in reversed(xrange(self.rows)):
            if self[i, i] == 0:
                raise ValueError("Matrix must be non-singular.")
            X[i, 0] = (rhs[i, 0] - sum(self[i, k] * X[k, 0]
                for k in xrange(i+1, self.rows))) / self[i, i]
        return self._new(X)

    def cholesky_solve(self, rhs):
        """
        Solves Ax = B using Cholesky decomposition,
        for a general square non-singular matrix.
        For a non-square matrix with rows > cols,
        the least squares solution is returned.

        See Also
        ========

        lower_triangular_solve
        upper_triangular_solve
        diagonal_solve
        LDLsolve
        LUsolve
        QRsolve
        """
        if self.is_symmetric():
            L = self._cholesky()
        elif self.rows >= self.cols:
            L = (self.T * self)._cholesky()
            rhs = self.T * rhs
        else:
            raise NotImplementedError("Under-determined System.")
        Y = L._lower_triangular_solve(rhs)
        return (L.T)._upper_triangular_solve(Y)

    def diagonal_solve(self, rhs):
        """
        Solves Ax = B efficiently, where A is a diagonal Matrix,
        with non-zero diagonal entries.

        See Also
        ========

        lower_triangular_solve
        upper_triangular_solve
        cholesky_solve
        LDLsolve
        LUsolve
        QRsolve
        """
        if not self.is_diagonal:
            raise TypeError("Matrix should be diagonal")
        if rhs.rows != self.rows:
            raise TypeError("Size mis-match")
        return self._diagonal_solve(rhs)

    def _diagonal_solve(self, rhs):
        """
        Helper function of function diagonal_solve,
        without the error checks, to be used privately.
        """
        return self._new(rhs.rows, 1, lambda i, j: rhs[i, 0] / self[i, i])

    def LDLsolve(self, rhs):
        """
        Solves Ax = B using LDL decomposition,
        for a general square and non-singular matrix.

        For a non-square matrix with rows > cols,
        the least squares solution is returned.

        See Also
        ========

        LDLdecomposition
        lower_triangular_solve
        upper_triangular_solve
        cholesky_solve
        diagonal_solve
        LUsolve
        QRsolve
        """
        if self.is_symmetric():
            L, D = self.LDLdecomposition()
        elif self.rows >= self.cols:
            L, D = (self.T * self).LDLdecomposition()
            rhs = self.T * rhs
        else:
            raise NotImplementedError("Under-determined System.")
        Y = L._lower_triangular_solve(rhs)
        Z = D._diagonal_solve(Y)
        return (L.T)._upper_triangular_solve(Z)

    def inv(self, method="GE", iszerofunc=_iszero, try_block_diag=False):
        """
        Calculates the matrix inverse.

        According to the "method" parameter, it calls the appropriate method:

          GE .... inverse_GE()
          LU .... inverse_LU()
          ADJ ... inverse_ADJ()

        According to the "try_block_diag" parameter, it will try to form block
        diagonal matrices using the method get_diag_blocks(), invert these
        individually, and then reconstruct the full inverse matrix.

        Note, the GE and LU methods may require the matrix to be simplified
        before it is inverted in order to properly detect zeros during
        pivoting. In difficult cases a custom zero detection function can
        be provided by setting the iszerosfunc argument to a function that
        should return True if its argument is zero. The ADJ routine computes
        the determinant and uses that to detect singular matrices in addition
        to testing for zeros on the diagonal.

        See Also
        ========

        inverse_LU
        inverse_GE
        inverse_ADJ
        """
        if not self.is_square:
            raise NonSquareMatrixError()
        if try_block_diag:
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


    def __mathml__(self):
        mml = ""
        for i in range(self.rows):
            mml += "<matrixrow>"
            for j in range(self.cols):
                mml += self[i,j].__mathml__()
            mml += "</matrixrow>"
        return "<matrix>" + mml + "</matrix>"

    def row_join(self, rhs):
        """
        Concatenates two matrices along self's last and rhs's first column

        >>> from sympy import Matrix
        >>> M = Matrix(3,3,lambda i,j: i+j)
        >>> V = Matrix(3,1,lambda i,j: 3+i+j)
        >>> M.row_join(V)
        [0, 1, 2, 3]
        [1, 2, 3, 4]
        [2, 3, 4, 5]

        See Also
        ========

        row
        col_join
        """
        if self.rows != rhs.rows:
            raise ShapeError("`self` and `rhs` must have the same number of rows.")

        newmat = self.zeros(self.rows, self.cols + rhs.cols)
        newmat[:,:self.cols] = self[:,:]
        newmat[:,self.cols:] = rhs
        return newmat

    def col_join(self, bott):
        """
        Concatenates two matrices along self's last and bott's first row

        >>> from sympy import Matrix, ones
        >>> M = ones(3, 3)
        >>> V = Matrix([[7,7,7]])
        >>> M.col_join(V)
        [1, 1, 1]
        [1, 1, 1]
        [1, 1, 1]
        [7, 7, 7]

        See Also
        ========

        col
        row_join
        """
        if self.cols != bott.cols:
            raise ShapeError("`self` and `bott` must have the same number of columns.")

        newmat = self.zeros(self.rows+bott.rows, self.cols)
        newmat[:self.rows,:] = self[:,:]
        newmat[self.rows:,:] = bott
        return newmat

    def row_insert(self, pos, mti):
        """
        Insert a row at the given position.

        >>> from sympy import Matrix, zeros
        >>> M = Matrix(3,3,lambda i,j: i+j)
        >>> M
        [0, 1, 2]
        [1, 2, 3]
        [2, 3, 4]
        >>> V = zeros(1, 3)
        >>> V
        [0, 0, 0]
        >>> M.row_insert(1,V)
        [0, 1, 2]
        [0, 0, 0]
        [1, 2, 3]
        [2, 3, 4]

        See Also
        ========

        row
        col_insert
        """
        if pos == 0:
            return mti.col_join(self)
        elif pos < 0:
            pos = self.rows + pos
        if pos < 0:
            pos = 0
        elif pos > self.rows:
            pos = self.rows

        if self.cols != mti.cols:
            raise ShapeError("`self` and `mti` must have the same number of columns.")

        newmat = self.zeros(self.rows + mti.rows, self.cols)
        newmat[:pos,:] = self[:pos,:]
        newmat[pos:pos+mti.rows,:] = mti[:,:]
        newmat[pos+mti.rows:,:] = self[pos:,:]
        return newmat

    def col_insert(self, pos, mti):
        """
        Insert a column at the given position.

        >>> from sympy import Matrix, zeros
        >>> M = Matrix(3,3,lambda i,j: i+j)
        >>> M
        [0, 1, 2]
        [1, 2, 3]
        [2, 3, 4]
        >>> V = zeros(3, 1)
        >>> V
        [0]
        [0]
        [0]
        >>> M.col_insert(1,V)
        [0, 0, 1, 2]
        [1, 0, 2, 3]
        [2, 0, 3, 4]

        See Also
        ========

        col
        row_insert
        """
        if pos == 0:
            return mti.row_join(self)
        elif pos < 0:
            pos = self.cols + pos
        if pos < 0:
            pos = 0
        elif pos > self.cols:
            pos = self.cols

        if self.rows != mti.rows:
            raise ShapeError("self and mti must have the same number of rows.")

        newmat = self.zeros(self.rows, self.cols + mti.cols)
        newmat[:,:pos] = self[:,:pos]
        newmat[:,pos:pos+mti.cols] = mti[:,:]
        newmat[:,pos+mti.cols:] = self[:,pos:]
        return newmat

    def trace(self):
        """
        Calculate the trace of a (square) matrix.

        >>> import sympy
        >>> M = sympy.matrices.eye(3)
        >>> M.trace()
        3

        """
        if not self.is_square:
            raise NonSquareMatrixError()

        trace = 0
        for i in range(self.cols):
            trace += self[i,i]
        return trace

    def submatrix(self, keys):
        """
        Get a slice/submatrix of the matrix using the given slice.

        >>> from sympy import Matrix
        >>> m = Matrix(4,4,lambda i,j: i+j)
        >>> m   #doctest: +NORMALIZE_WHITESPACE
        [0, 1, 2, 3]
        [1, 2, 3, 4]
        [2, 3, 4, 5]
        [3, 4, 5, 6]
        >>> m[:1, 1]   #doctest: +NORMALIZE_WHITESPACE
        [1]
        >>> m[:2, :1] #doctest: +NORMALIZE_WHITESPACE
        [0]
        [1]
        >>> m[2:4, 2:4] #doctest: +NORMALIZE_WHITESPACE
        [4, 5]
        [5, 6]

        See Also
        ========

        extract
        """
        rlo, rhi, clo, chi = self.key2bounds(keys)
        outLines, outCols = rhi-rlo, chi-clo
        outMat = [0]*outLines*outCols
        for i in xrange(outLines):
            outMat[i*outCols:(i+1)*outCols] = self.mat[(i+rlo)*self.cols+clo:(i+rlo)*self.cols+chi]
        return self._new(outLines,outCols,outMat)

    def extract(self, rowsList, colsList):
        """
        Extract a submatrix by specifying a list of rows and columns.
        Negative indices can be given. All indices must be in the range
        -n <= i < n where n is the number of rows or columns.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(4, 3, range(12))
        >>> m   #doctest: +NORMALIZE_WHITESPACE
        [0,  1,  2]
        [3,  4,  5]
        [6,  7,  8]
        [9, 10, 11]
        >>> m.extract([0,1,3],[0,1])   #doctest: +NORMALIZE_WHITESPACE
        [0,  1]
        [3,  4]
        [9, 10]

        Rows or columns can be repeated:

        >>> m.extract([0,0,1], [-1])   #doctest: +NORMALIZE_WHITESPACE
        [2]
        [2]
        [5]

        Every other row can be taken by using range to provide the indices:

        >>> m.extract(range(0, m.rows, 2),[-1])   #doctest: +NORMALIZE_WHITESPACE
        [2]
        [8]

        See Also
        ========

        submatrix
        """
        cols = self.cols
        mat = self.mat
        rowsList = [self.key2ij((k,0))[0] for k in rowsList]
        colsList = [self.key2ij((0,k))[1] for k in colsList]
        return self._new(len(rowsList), len(colsList),
                lambda i,j: mat[rowsList[i]*cols + colsList[j]])

    def key2bounds(self, keys):
        """Converts a key with potentially mixed types of keys (integer and slice)
        into a tuple of ranges and raises an error if any index is out of self's
        range.

        See Also
        ========

        key2ij
        slice2bounds
        """

        islice, jslice = [isinstance(k, slice) for k in keys]
        if islice:
            rlo, rhi = self.slice2bounds(keys[0], self.rows)
        else:
            # assuming we don't have a
            rlo, _ = self.slice2bounds((keys[0], 0), self.rows)
            rhi  = rlo + 1
        if jslice:
            clo, chi = self.slice2bounds(keys[1], self.cols)
        else:
            _, clo = self.slice2bounds((0, keys[1]), self.cols)
            chi = clo + 1
        if not ( 0<=rlo<=rhi and 0<=clo<=chi ):
            raise IndexError("Slice indices out of range: a[%s]"%repr(keys))
        return rlo, rhi, clo, chi

    def key2ij(self, key):
        """Converts key=(4,6) to 4,6 and ensures the key is correct. Negative
        indices are also supported and are remapped to positives provided they
        are valid indices.

        See Also
        ========

        key2bounds
        slice2bounds
        """

        if not (is_sequence(key) and len(key) == 2):
            raise TypeError("wrong syntax: a[%s]. Use a[i,j] or a[(i,j)]"
                    %repr(key))
        i, j = [(k + n) if k < 0 else k for k, n in zip(key, (self.rows, self.cols))]
        if not (i>=0 and i<self.rows and j>=0 and j < self.cols):
            raise IndexError("Index out of range: a[%s]"%repr(key))
        return i,j

    def slice2bounds(self, key, defmax):
        """
        Takes slice or number and returns (min,max) for iteration
        Takes a default maxval to deal with the slice ':' which is (none, none)

        See Also
        ========

        key2bounds
        key2ij
        """
        if isinstance(key, slice):
            return key.indices(defmax)[:2]
        elif isinstance(key, tuple):
            key = [defmax - 1 if i == -1 else i for i in key]
            return self.key2ij(key)
        else:
            raise IndexError("Improper index type")

    def applyfunc(self, f):
        """
        Apply a function to each element of the matrix.

        >>> from sympy import Matrix
        >>> m = Matrix(2,2,lambda i,j: i*2+j)
        >>> m   #doctest: +NORMALIZE_WHITESPACE
        [0, 1]
        [2, 3]
        >>> m.applyfunc(lambda i: 2*i)  #doctest: +NORMALIZE_WHITESPACE
        [0, 2]
        [4, 6]

        """
        if not callable(f):
            raise TypeError("`f` must be callable.")

        out = self._new(self.rows, self.cols, map(f,self.mat))
        return out

    def evalf(self, prec=None, **options):
        """
        Evaluate each element of the matrix as a float.
        """
        if prec is None:
            return self.applyfunc(lambda i: i.evalf(**options))
        else:
            return self.applyfunc(lambda i: i.evalf(prec, **options))

    n = evalf

    def reshape(self, _rows, _cols):
        """
        Reshape the matrix. Total number of elements must remain the same.

        >>> from sympy import Matrix
        >>> m = Matrix(2,3,lambda i,j: 1)
        >>> m   #doctest: +NORMALIZE_WHITESPACE
        [1, 1, 1]
        [1, 1, 1]
        >>> m.reshape(1,6)  #doctest: +NORMALIZE_WHITESPACE
        [1, 1, 1, 1, 1, 1]
        >>> m.reshape(3,2)  #doctest: +NORMALIZE_WHITESPACE
        [1, 1]
        [1, 1]
        [1, 1]

        """
        if len(self) != _rows*_cols:
            raise ValueError("Invalid reshape parameters %d %d" % (_rows, _cols))
        return self._new(_rows, _cols, lambda i,j: self.mat[i*_cols + j])

    def print_nonzero (self, symb="X"):
        """
        Shows location of non-zero entries for fast shape lookup.

        Examples
        ========

        >>> from sympy import Matrix, matrices
        >>> m = Matrix(2,3,lambda i,j: i*3+j)
        >>> m           #doctest: +NORMALIZE_WHITESPACE
        [0, 1, 2]
        [3, 4, 5]
        >>> m.print_nonzero()   #doctest: +NORMALIZE_WHITESPACE
        [ XX]
        [XXX]
        >>> m = matrices.eye(4)
        >>> m.print_nonzero("x")    #doctest: +NORMALIZE_WHITESPACE
        [x   ]
        [ x  ]
        [  x ]
        [   x]

        """
        s = []
        for i in range(self.rows):
            line = []
            for j in range(self.cols):
                if self[i,j] == 0:
                    line.append(" ")
                else:
                    line.append(str(symb))
            s.append("[%s]" % ''.join(line))
        print '\n'.join(s)

    def LUsolve(self, rhs, iszerofunc=_iszero):
        """
        Solve the linear system Ax = b for x.
        self is the coefficient matrix A and rhs is the right side b.

        This is for symbolic matrices, for real or complex ones use
        sympy.mpmath.lu_solve or sympy.mpmath.qr_solve.

        See Also
        ========

        lower_triangular_solve
        upper_triangular_solve
        cholesky_solve
        diagonal_solve
        LDLsolve
        QRsolve
        LUdecomposition
        """
        if rhs.rows != self.rows:
            raise ShapeError("`self` and `rhs` must have the same number of rows.")

        A, perm = self.LUdecomposition_Simple(iszerofunc=_iszero)
        n = self.rows
        b = rhs.permuteFwd(perm).as_mutable()
        # forward substitution, all diag entries are scaled to 1
        for i in range(n):
            for j in range(i):
                b.row(i, lambda x,k: x - b[j,k]*A[i,j])
        # backward substitution
        for i in range(n-1,-1,-1):
            for j in range(i+1, n):
                b.row(i, lambda x,k: x - b[j,k]*A[i,j])
            b.row(i, lambda x,k: x / A[i,i])
        return rhs.__class__(b)

    def LUdecomposition(self, iszerofunc=_iszero):
        """
        Returns the decomposition LU and the row swaps p.

        Examples
        ========

        >>> from sympy import Matrix
        >>> a = Matrix([[4, 3], [6, 3]])
        >>> L, U, _ = a.LUdecomposition()
        >>> L
        [  1, 0]
        [3/2, 1]
        >>> U
        [4,    3]
        [0, -3/2]

        See Also
        ========

        cholesky
        LDLdecomposition
        QRdecomposition
        LUdecomposition_Simple
        LUdecompositionFF
        LUsolve
        """
        combined, p = self.LUdecomposition_Simple(iszerofunc=_iszero)
        L = self.zeros(self.rows)
        U = self.zeros(self.rows)
        for i in range(self.rows):
            for j in range(self.rows):
                if i > j:
                    L[i,j] = combined[i,j]
                else:
                    if i == j:
                        L[i,i] = 1
                    U[i,j] = combined[i,j]
        return L, U, p

    def LUdecomposition_Simple(self, iszerofunc=_iszero):
        """
        Returns A comprised of L,U (L's diag entries are 1) and
        p which is the list of the row swaps (in order).

        See Also
        ========

        LUdecomposition
        LUdecompositionFF
        LUsolve
        """
        if not self.is_square:
            raise NonSquareMatrixError()
        n = self.rows
        A = self.as_mutable()
        p = []
        # factorization
        for j in range(n):
            for i in range(j):
                for k in range(i):
                    A[i,j] = A[i,j] - A[i,k]*A[k,j]
            pivot = -1
            for i in range(j,n):
                for k in range(j):
                    A[i,j] = A[i,j] - A[i,k]*A[k,j]
                # find the first non-zero pivot, includes any expression
                if pivot == -1 and not iszerofunc(A[i,j]):
                    pivot = i
            if pivot < 0:
                # this result is based on iszerofunc's analysis of the possible pivots, so even though
                # the element may not be strictly zero, the supplied iszerofunc's evaluation gave True
                raise ValueError("No nonzero pivot found; inversion failed.")
            if pivot != j: # row must be swapped
                A.row_swap(pivot,j)
                p.append([pivot,j])
            scale = 1 / A[j,j]
            for i in range(j+1,n):
                A[i,j] = A[i,j] * scale
        return A, p


    def LUdecompositionFF(self):
        """
        Compute a fraction-free LU decomposition.

        Returns 4 matrices P, L, D, U such that PA = L D**-1 U.
        If the elements of the matrix belong to some integral domain I, then all
        elements of L, D and U are guaranteed to belong to I.

        **Reference**
            - W. Zhou & D.J. Jeffrey, "Fraction-free matrix factors: new forms
              for LU and QR factors". Frontiers in Computer Science in China,
              Vol 2, no. 1, pp. 67-80, 2008.

        See Also
        ========

        LUdecomposition
        LUdecomposition_Simple
        LUsolve
        """
        n, m = self.rows, self.cols
        U, L, P = self.as_mutable(), eye(n), eye(n)
        DD = zeros(n) # store it smarter since it's just diagonal
        oldpivot = 1

        for k in range(n-1):
            if U[k,k] == 0:
                for kpivot in range(k+1, n):
                    if U[kpivot, k] != 0:
                        break
                else:
                    raise ValueError("Matrix is not full rank")
                U[k, k:], U[kpivot, k:] = U[kpivot, k:], U[k, k:]
                L[k, :k], L[kpivot, :k] = L[kpivot, :k], L[k, :k]
                P[k, :], P[kpivot, :] = P[kpivot, :], P[k, :]
            L[k,k] = Ukk = U[k,k]
            DD[k,k] = oldpivot * Ukk
            for i in range(k+1, n):
                L[i,k] = Uik = U[i,k]
                for j in range(k+1, m):
                    U[i,j] = (Ukk * U[i,j] - U[k,j]*Uik) / oldpivot
                U[i,k] = 0
            oldpivot = Ukk
        DD[n-1,n-1] = oldpivot
        return P, L, DD, U

    def cofactorMatrix(self, method="berkowitz"):
        """
        Return a matrix containing the cofactor of each element.

        See Also
        ========

        cofactor
        minorEntry
        minorMatrix
        adjugate
        """
        out = self._new(self.rows, self.cols, lambda i,j:
                self.cofactor(i, j, method))
        return out

    def minorEntry(self, i, j, method="berkowitz"):
        """
        Calculate the minor of an element.

        See Also
        ========

        minorMatrix
        cofactor
        cofactorMatrix
        """
        if not 0 <= i < self.rows or not 0 <= j < self.cols:
            raise ValueError("`i` and `j` must satisfy 0 <= i < `self.rows` " +
                "(%d)" % self.rows + "and 0 <= j < `self.cols` (%d)." % self.cols)
        return self.minorMatrix(i,j).det(method)

    def minorMatrix(self, i, j):
        """
        Creates the minor matrix of a given element.

        See Also
        ========

        minorEntry
        cofactor
        cofactorMatrix
        """
        if not 0 <= i < self.rows or not 0 <= j < self.cols:
            raise ValueError("`i` and `j` must satisfy 0 <= i < `self.rows` " +
                "(%d)" % self.rows + "and 0 <= j < `self.cols` (%d)." % self.cols)
        return self.delRowCol(i,j)

    def cofactor(self, i, j, method="berkowitz"):
        """
        Calculate the cofactor of an element.

        See Also
        ========

        cofactorMatrix
        minorEntry
        minorMatrix
        """
        if (i+j) % 2 == 0:
            return self.minorEntry(i, j, method)
        else:
            return -1 * self.minorEntry(i, j, method)

    def jacobian(self, X):
        """
        Calculates the Jacobian matrix (derivative of a vectorial function).

        *self*
            A vector of expressions representing functions f_i(x_1, ..., x_n).
        *X*
            The set of x_i's in order, it can be a list or a Matrix

        Both self and X can be a row or a column matrix in any order
        (jacobian() should always work).

        Examples
        ========

            >>> from sympy import sin, cos, Matrix
            >>> from sympy.abc import rho, phi
            >>> X = Matrix([rho*cos(phi), rho*sin(phi), rho**2])
            >>> Y = Matrix([rho, phi])
            >>> X.jacobian(Y)
            [cos(phi), -rho*sin(phi)]
            [sin(phi),  rho*cos(phi)]
            [   2*rho,             0]
            >>> X = Matrix([rho*cos(phi), rho*sin(phi)])
            >>> X.jacobian(Y)
            [cos(phi), -rho*sin(phi)]
            [sin(phi),  rho*cos(phi)]

        See Also
        ========

        hessian
        wronskian
        """
        if not isinstance(X, MatrixBase):
            X = self._new(X)
        # Both X and self can be a row or a column matrix, so we need to make
        # sure all valid combinations work, but everything else fails:
        if self.shape[0] == 1:
            m = self.shape[1]
        elif self.shape[1] == 1:
            m = self.shape[0]
        else:
            raise TypeError("self must be a row or a column matrix")
        if X.shape[0] == 1:
            n = X.shape[1]
        elif X.shape[1] == 1:
            n = X.shape[0]
        else:
            raise TypeError("X must be a row or a column matrix")

        # m is the number of functions and n is the number of variables
        # computing the Jacobian is now easy:
        return self._new(m, n, lambda j, i: self[j].diff(X[i]))

    def QRdecomposition(self):
        """
        Return Q,R where A = Q*R, Q is orthogonal and R is upper triangular.

        Examples
        ========

        This is the example from wikipedia:

        >>> from sympy import Matrix, eye
        >>> A = Matrix([[12,-51,4],[6,167,-68],[-4,24,-41]])
        >>> Q, R = A.QRdecomposition()
        >>> Q
        [ 6/7, -69/175, -58/175]
        [ 3/7, 158/175,   6/175]
        [-2/7,    6/35,  -33/35]
        >>> R
        [14,  21, -14]
        [ 0, 175, -70]
        [ 0,   0,  35]
        >>> A == Q*R
        True

        QR factorization of an identity matrix:

        >>> A = Matrix([[1,0,0],[0,1,0],[0,0,1]])
        >>> Q, R = A.QRdecomposition()
        >>> Q
        [1, 0, 0]
        [0, 1, 0]
        [0, 0, 1]
        >>> R
        [1, 0, 0]
        [0, 1, 0]
        [0, 0, 1]

        See Also
        ========

        cholesky
        LDLdecomposition
        LUdecomposition
        QRsolve
        """
        cls = self.__class__
        self = self.as_mutable()

        if not self.rows >= self.cols:
            raise MatrixError("The number of rows must be greater than columns")
        n = self.rows
        m = self.cols
        rank = n
        row_reduced = self.rref()[0]
        for i in range(row_reduced.rows):
            if MutableMatrix(row_reduced[i*m:(i+1)*m]).norm() == 0:
                rank -= 1
        if not rank == self.cols:
            raise MatrixError("The rank of the matrix must match the columns")
        Q, R = self.zeros(n, m), self.zeros(m)
        for j in range(m):      # for each column vector
            tmp = self[:,j]     # take original v
            for i in range(j):
                # subtract the project of self on new vector
                tmp -= Q[:,i] * self[:,j].dot(Q[:,i])
                tmp.expand()
            # normalize it
            R[j,j] = tmp.norm()
            Q[:,j] = tmp / R[j,j]
            if Q[:,j].norm() != 1:
                raise NotImplementedError("Could not normalize the vector %d." % j)
            for i in range(j):
                R[i,j] = Q[:,i].dot(self[:,j])
        return cls(Q), cls(R)

    def QRsolve(self, b):
        """
        Solve the linear system 'Ax = b'.

        'self' is the matrix 'A', the method argument is the vector
        'b'.  The method returns the solution vector 'x'.  If 'b' is a
        matrix, the system is solved for each column of 'b' and the
        return value is a matrix of the same shape as 'b'.

        This method is slower (approximately by a factor of 2) but
        more stable for floating-point arithmetic than the LUsolve method.
        However, LUsolve usually uses an exact arithmetic, so you don't need
        to use QRsolve.

        This is mainly for educational purposes and symbolic matrices, for real
        (or complex) matrices use sympy.mpmath.qr_solve.

        See Also
        ========

        lower_triangular_solve
        upper_triangular_solve
        cholesky_solve
        diagonal_solve
        LDLsolve
        LUsolve
        QRdecomposition
        """

        Q, R = self.as_mutable().QRdecomposition()
        y = Q.T * b

        # back substitution to solve R*x = y:
        # We build up the result "backwards" in the vector 'x' and reverse it
        # only in the end.
        x = []
        n = R.rows
        for j in range(n - 1, -1, -1):
            tmp = y[j,:]
            for k in range(j+1, n):
                tmp -= R[j,k] * x[n-1-k]
            x.append(tmp/R[j,j])
        return self._new([row.mat for row in reversed(x)])

    #def evaluate(self):    # no more eval() so should be removed
    #    for i in range(self.rows):
    #        for j in range(self.cols):
    #            self[i,j] = self[i,j].eval()

    def cross(self, b):
        """
        Calculate the cross product of `self` and `b`.

        See Also
        ========

        dot
        multiply
        multiply_elementwise
        """
        if not is_sequence(b, include=MatrixBase):
            raise TypeError("`b` must be an ordered iterable or Matrix, not %s." %
                type(b))
        if not (self.rows == 1 and self.cols == 3 or \
                self.rows == 3 and self.cols == 1 ) and \
                (b.rows == 1 and b.cols == 3 or \
                b.rows == 3 and b.cols == 1):
            raise ShapeError("Dimensions incorrect for cross product.")
        else:
            return self._new(1,3,((self[1]*b[2] - self[2]*b[1]),
                               (self[2]*b[0] - self[0]*b[2]),
                               (self[0]*b[1] - self[1]*b[0])))

    def dot(self, b):
        """Return the dot product of Matrix self and b relaxing the condition
        of compatible dimensions: if either the number of rows or columns are
        the same as the length of b then the dot product is returned. If self
        is a row or column vector, a scalar is returned. Otherwise, a list
        of results is returned (and in that case the number of columns in self
        must match the length of b).

        >>> from sympy import Matrix
        >>> M = Matrix([[1,2,3], [4,5,6], [7,8,9]])
        >>> v = [1, 1, 1]
        >>> M.row(0).dot(v)
        6
        >>> M.col(0).dot(v)
        12
        >>> M.dot(v)
        [6, 15, 24]

        See Also
        ========

        cross
        multiply
        multiply_elementwise
        """
        if not isinstance(b, MatrixBase):
            if is_sequence(b):
                if len(b) != self.cols and len(b) != self.rows:
                    raise ShapeError("Dimensions incorrect for dot product.")
                return self.dot(MutableMatrix(b))
            else:
                raise TypeError("`b` must be an ordered iterable or Matrix, not %s." %
                type(b))
        if self.cols == b.rows:
            if b.cols != 1:
                self = self.T
                b = b.T
            prod = flatten((self*b).tolist())
            if len(prod) == 1:
                return prod[0]
            return prod
        if self.cols == b.cols:
            return self.dot(b.T)
        elif self.rows == b.rows:
            return self.T.dot(b)
        else:
            raise ShapeError("Dimensions incorrect for dot product.")

    def multiply_elementwise(self, b):
        """Return the Hadamard product (elementwise product) of A and B

        >>> from sympy.matrices.matrices import Matrix
        >>> A = Matrix([[0, 1, 2], [3, 4, 5]])
        >>> B = Matrix([[1, 10, 100], [100, 10, 1]])
        >>> A.multiply_elementwise(B)
        [  0, 10, 200]
        [300, 40,   5]

        See Also
        ========

        cross
        dot
        multiply
        """
        return matrix_multiply_elementwise(self, b)

    def norm(self, ord=None):
        """Return the Norm of a Matrix or Vector.
        In the simplest case this is the geometric size of the vector
        Other norms can be specified by the ord parameter


        =====  ============================  ==========================
        ord    norm for matrices             norm for vectors
        =====  ============================  ==========================
        None   Frobenius norm                2-norm
        'fro'  Frobenius norm                - does not exist
        inf    --                            max(abs(x))
        -inf   --                            min(abs(x))
        1      --                            as below
        -1     --                            as below
        2      2-norm (largest sing. value)  as below
        -2     smallest singular value       as below
        other  - does not exist              sum(abs(x)**ord)**(1./ord)
        =====  ============================  ==========================

        >>> from sympy import Matrix, Symbol, trigsimp, cos, sin
        >>> x = Symbol('x', real=True)
        >>> v = Matrix([cos(x), sin(x)])
        >>> trigsimp( v.norm() )
        1
        >>> v.norm(10)
        (sin(x)**10 + cos(x)**10)**(1/10)
        >>> A = Matrix([[1,1], [1,1]])
        >>> A.norm(2)# Spectral norm (max of |Ax|/|x| under 2-vector-norm)
        2
        >>> A.norm(-2) # Inverse spectral norm (smallest singular value)
        0
        >>> A.norm() # Frobenius Norm
        2

        See Also
        ========

        normalized
        """
        # Row or Column Vector Norms
        if self.rows == 1 or self.cols == 1:
            if ord == 2 or ord == None: # Common case sqrt(<x,x>)
                return sqrt(Add(*(abs(i)**2 for i in self.mat)))

            elif ord == 1: # sum(abs(x))
                return Add(*(abs(i) for i in self.mat))

            elif ord == S.Infinity: # max(abs(x))
                return Max(*self.applyfunc(abs))

            elif ord == S.NegativeInfinity: # min(abs(x))
                return Min(*self.applyfunc(abs))

            # Otherwise generalize the 2-norm, Sum(x_i**ord)**(1/ord)
            # Note that while useful this is not mathematically a norm
            try:
                return Pow( Add(*(abs(i)**ord for i in self.mat)), S(1)/ord )
            except TypeError:
                raise ValueError("Expected order to be Number, Symbol, oo")

        # Matrix Norms
        else:
            if ord == 2: # Spectral Norm
                # Maximum singular value
                return Max(*self.singular_values())

            elif ord == -2:
                # Minimum singular value
                return Min(*self.singular_values())

            elif (ord == None or isinstance(ord,str) and ord.lower() in
                    ['f', 'fro', 'frobenius', 'vector']):
                # Reshape as vector and send back to norm function
                return self.vec().norm(ord=2)

            else:
                raise NotImplementedError("Matrix Norms under development")

    def normalized(self):
        """
        Return the normalized version of `self`.

        See Also
        ========

        norm
        """
        if self.rows != 1 and self.cols != 1:
            raise ShapeError("A Matrix must be a vector to normalize.")
        norm = self.norm()
        out = self.applyfunc(lambda i: i / norm)
        return out

    def project(self, v):
        """Return the projection of ``self`` onto the line containing ``v``.

        >>> from sympy import Matrix, S, sqrt
        >>> V = Matrix([sqrt(3)/2,S.Half])
        >>> x = Matrix([[1, 0]])
        >>> V.project(x)
        [sqrt(3)/2, 0]
        >>> V.project(-x)
        [sqrt(3)/2, 0]
        """
        return v * (self.dot(v) / v.dot(v))

    def permuteBkwd(self, perm):
        """
        Permute the rows of the matrix with the given permutation in reverse.

        >>> import sympy
        >>> M = sympy.matrices.eye(3)
        >>> M.permuteBkwd([[0,1],[0,2]])
        [0, 1, 0]
        [0, 0, 1]
        [1, 0, 0]

        See Also
        ========

        permuteFwd
        """
        copy = self[:,:]
        for i in range(len(perm)-1, -1, -1):
            copy.row_swap(perm[i][0], perm[i][1])
        return copy

    def permuteFwd(self, perm):
        """
        Permute the rows of the matrix with the given permutation.

        >>> import sympy
        >>> M = sympy.matrices.eye(3)
        >>> M.permuteFwd([[0,1],[0,2]])
        [0, 0, 1]
        [1, 0, 0]
        [0, 1, 0]

        See Also
        ========

        permuteBkwd
        """
        copy = self[:,:]
        for i in range(len(perm)):
            copy.row_swap(perm[i][0], perm[i][1])
        return copy

    def delRowCol(self, i, j):
        """
        Creates a copy of the matrix with the given row and column deleted.

        Examples
        ========

        >>> import sympy
        >>> I = sympy.matrices.eye(4)
        >>> I.delRowCol(1,2)
        [1, 0, 0]
        [0, 0, 0]
        [0, 0, 1]

        See Also
        ========

        row_del
        col_del
        """
        # used only for cofactors, makes a copy
        M = self.as_mutable()
        M.row_del(i)
        M.col_del(j)
        return self._new(M)

    def exp(self):
        """ Returns the exponentiation of a matrix

        See Also
        ========

        matrix_multiply
        """
        if not self.is_square:
            raise NonSquareMatrixError("Exponentiation is valid only for square matrices")
        try:
            U, D = self.diagonalize()
        except MatrixError:
            raise NotImplementedError("Exponentiation is implemented only for diagonalizable matrices")
        for i in xrange(0, D.rows):
            D[i, i] = C.exp(D[i, i])
        return U * D * U.inv()

    @classmethod
    def zeros(cls, r, c=None):
        """Returns a matrix of zeros with ``r`` rows and ``c`` columns;
        if ``c`` is omitted a square matrix will be returned.

        See Also
        ========

        ones
        fill
        """
        if is_sequence(r):
            SymPyDeprecationWarning(
            "Pass row and column count as zeros(%i, %i)." % r
            ).warn()
            r, c = r
        else:
            c = r if c is None else c
        r, c = [int(i) for i in [r, c]]
        return MutableMatrix(r, c, [S.Zero]*r*c)

    @classmethod
    def eye(cls, n):
        """Returns the identity matrix of size n."""
        tmp = cls.zeros(n)
        for i in range(tmp.rows):
            tmp[i, i] = S.One
        return tmp

    @property
    def is_square(self):
        """
        Checks if a matrix is square.

        A matrix is square if the number of rows equals the number of columns.
        The empty matrix is square by definition, since the number of rows and
        the number of columns are both zero.

        Examples
        ========

        >>> from sympy import Matrix
        >>> a = Matrix([[1, 2, 3], [4, 5, 6]])
        >>> b = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> c = Matrix([])
        >>> a.is_square
        False
        >>> b.is_square
        True
        >>> c.is_square
        True
        """
        return self.rows == self.cols

    @property
    def is_zero(self):
        """
        Checks if a matrix is a zero matrix.

        A matrix is zero if every element is zero.  A matrix need not be square
        to be considered zero.  The empty matrix is zero by the principle of
        vacuous truth.

        Examples
        ========

        >>> from sympy import Matrix, zeros
        >>> a = Matrix([[0, 0], [0, 0]])
        >>> b = zeros(3, 4)
        >>> c = Matrix([[0, 1], [0, 0]])
        >>> d = Matrix([])
        >>> a.is_zero
        True
        >>> b.is_zero
        True
        >>> c.is_zero
        False
        >>> d.is_zero
        True
        """
        return all(i.is_zero for i in self)

    def is_nilpotent(self):
        """
        Checks if a matrix is nilpotent.

        A matrix B is nilpotent if for some integer k, B**k is
        a zero matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> a = Matrix([[0,0,0],[1,0,0],[1,1,0]])
        >>> a.is_nilpotent()
        True

        >>> a = Matrix([[1,0,1],[1,0,0],[1,1,0]])
        >>> a.is_nilpotent()
        False
        """
        if not self.is_square:
            raise NonSquareMatrixError("Nilpotency is valid only for square matrices")
        x = Dummy('x')
        if self.charpoly(x).args[0] == x**self.rows:
            return True
        return False

    def is_upper(self):
        """
        Check if matrix is an upper triangular matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(2,2,[1, 0, 0, 1])
        >>> m
        [1, 0]
        [0, 1]
        >>> m.is_upper()
        True

        >>> m = Matrix(3,3,[5, 1, 9, 0, 4 , 6, 0, 0, 5])
        >>> m
        [5, 1, 9]
        [0, 4, 6]
        [0, 0, 5]
        >>> m.is_upper()
        True

        >>> m = Matrix(2,3,[4, 2, 5, 6, 1, 1])
        >>> m
        [4, 2, 5]
        [6, 1, 1]
        >>> m.is_upper()
        False

        See Also
        ========

        is_lower
        is_diagonal
        is_upper_hessenberg
        """
        for i in xrange(1, self.rows):
            for j in xrange(0, i):
                if self[i,j] != 0:
                    return False
        return True

    def is_lower(self):
        """
        Check if matrix is a lower triangular matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(2,2,[1, 0, 0, 1])
        >>> m
        [1, 0]
        [0, 1]
        >>> m.is_lower()
        True

        >>> m = Matrix(3,3,[2, 0, 0, 1, 4 , 0, 6, 6, 5])
        >>> m
        [2, 0, 0]
        [1, 4, 0]
        [6, 6, 5]
        >>> m.is_lower()
        True

        >>> from sympy.abc import x, y
        >>> m = Matrix(2,2,[x**2 + y, y**2 + x, 0, x + y])
        >>> m
        [x**2 + y, x + y**2]
        [       0,    x + y]
        >>> m.is_lower()
        False

        See Also
        ========

        is_upper
        is_diagonal
        is_lower_hessenberg
        """
        for i in xrange(0, self.rows):
            for j in xrange(i+1, self.cols):
                if self[i, j] != 0:
                    return False
        return True

    def is_upper_hessenberg(self):
        """
        Checks if the matrix is the upper hessenberg form.

        The upper hessenberg matrix has zero entries
        below the first subdiagonal.

        Examples
        ========

        >>> from sympy.matrices import Matrix
        >>> a = Matrix([[1,4,2,3],[3,4,1,7],[0,2,3,4],[0,0,1,3]])
        >>> a
        [1, 4, 2, 3]
        [3, 4, 1, 7]
        [0, 2, 3, 4]
        [0, 0, 1, 3]
        >>> a.is_upper_hessenberg()
        True

        See Also
        ========

        is_lower_hessenberg
        is_upper
        """
        for i in xrange(2, self.rows):
            for j in xrange(0, i - 1):
                if self[i,j] != 0:
                    return False
        return True

    def is_lower_hessenberg(self):
        r"""
        Checks if the matrix is in the lower hessenberg form.

        The lower hessenberg matrix has zero entries
        above the first superdiagonal.

        Examples
        ========

        >>> from sympy.matrices import Matrix
        >>> a = Matrix([[1,2,0,0],[5,2,3,0],[3,4,3,7],[5,6,1,1]])
        >>> a
        [1, 2, 0, 0]
        [5, 2, 3, 0]
        [3, 4, 3, 7]
        [5, 6, 1, 1]
        >>> a.is_lower_hessenberg()
        True

        See Also
        ========

        is_upper_hessenberg
        is_lower
        """
        for i in xrange(0, self.rows):
            for j in xrange(i + 2, self.cols):
                if self[i, j] != 0:
                    return False
        return True

    def is_symbolic(self):
        """
        Checks if any elements are symbols.

        Examples
        ========

        >>> import sympy
        >>> from sympy.abc import x, y
        >>> M = sympy.matrices.Matrix([[x, y], [1, 0]])
        >>> M.is_symbolic()
        True

        """
        return any(element.has(Symbol) for element in self.mat)

    def is_symmetric(self, simplify=True):
        """
        Check if matrix is symmetric matrix,
        that is square matrix and is equal to its transpose.

        By default, simplifications occur before testing symmetry.
        They can be skipped using 'simplify=False'; while speeding things a bit,
        this may however induce false negatives.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(2,2,[0, 1, 1, 2])
        >>> m
        [0, 1]
        [1, 2]
        >>> m.is_symmetric()
        True

        >>> m = Matrix(2,2,[0, 1, 2, 0])
        >>> m
        [0, 1]
        [2, 0]
        >>> m.is_symmetric()
        False

        >>> m = Matrix(2,3,[0, 0, 0, 0, 0, 0])
        >>> m
        [0, 0, 0]
        [0, 0, 0]
        >>> m.is_symmetric()
        False

        >>> from sympy.abc import x, y
        >>> m = Matrix(3,3,[1, x**2 + 2*x + 1, y, (x + 1)**2 , 2, 0, y, 0, 3])
        >>> m
        [         1, x**2 + 2*x + 1, y]
        [(x + 1)**2,              2, 0]
        [         y,              0, 3]
        >>> m.is_symmetric()
        True

        If the matrix is already simplified, you may speed-up is_symmetric()
        test by using 'simplify=False'.

        >>> m.is_symmetric(simplify=False)
        False
        >>> m1 = m.expand()
        >>> m1.is_symmetric(simplify=False)
        True
        """
        if not self.is_square:
            return False
        if simplify:
            delta = self - self.transpose()
            delta.simplify()
            return delta.equals(self.zeros(self.rows, self.cols))
        else:
            return self == self.transpose()

    def is_anti_symmetric(self, simplify=True):
        """
        Check if matrix M is an antisymmetric matrix,
        that is, M is a square matrix with all M[i, j] == -M[j, i].

        When `simplify=True` (default), the sum M[i, j] + M[j, i] is
        simplified before testing to see if it is zero. By default,
        the SymPy simplify function is used. To use a custom function
        set simplify to a function that accepts a single argument which
        returns a simplified expression. To skip simplification, set
        simplify to False but note that although this will be faster,
        it may induce false negatives.

        Examples
        ========

        >>> from sympy import Matrix, symbols
        >>> m = Matrix(2,2,[0, 1, -1, 0])
        >>> m
        [ 0, 1]
        [-1, 0]
        >>> m.is_anti_symmetric()
        True
        >>> x, y = symbols('x y')
        >>> m = Matrix(2,3,[0, 0, x, -y, 0, 0])
        >>> m
        [ 0, 0, x]
        [-y, 0, 0]
        >>> m.is_anti_symmetric()
        False

        >>> from sympy.abc import x, y
        >>> m = Matrix(3, 3, [0, x**2 + 2*x + 1, y,
        ...                   -(x + 1)**2 , 0, x*y,
        ...                   -y, -x*y, 0])

        Simplification of matrix elements is done by default so even
        though two elements which should be equal and opposite wouldn't
        pass an equality test, the matrix is still reported as
        anti-symmetric:

        >>> m[0, 1] == -m[1, 0]
        False
        >>> m.is_anti_symmetric()
        True

        If 'simplify=False' is used for the case when a Matrix is already
        simplified, this will speed things up. Here, we see that without
        simplification the matrix does not appear anti-symmetric:

        >>> m.is_anti_symmetric(simplify=False)
        False

        But if the matrix were already expanded, then it would appear
        anti-symmetric and simplification in the is_anti_symmetric routine
        is not needed:

        >>> m = m.expand()
        >>> m.is_anti_symmetric(simplify=False)
        True
        """
        # accept custom simplification
        simpfunc = simplify if isinstance(simplify, FunctionType) else \
                   _simplify if simplify else False

        if not self.is_square:
            return False
        n = self.rows
        if simplify:
            for i in xrange(n):
                # diagonal
                if not simpfunc(self[i, i]).is_zero:
                    return False
                # others
                for j in xrange(i + 1, n):
                    diff = self[i, j] + self[j, i]
                    if not simpfunc(diff).is_zero:
                        return False
            return True
        else:
            for i in xrange(n):
                for j in xrange(i, n):
                    if self[i, j] != -self[j, i]:
                        return False
            return True


    def is_diagonal(self):
        """
        Check if matrix is diagonal,
        that is matrix in which the entries outside the main diagonal are all zero.

        Examples
        ========

        >>> from sympy import Matrix, diag
        >>> m = Matrix(2,2,[1, 0, 0, 2])
        >>> m
        [1, 0]
        [0, 2]
        >>> m.is_diagonal()
        True

        >>> m = Matrix(2,2,[1, 1, 0, 2])
        >>> m
        [1, 1]
        [0, 2]
        >>> m.is_diagonal()
        False

        >>> m = diag(1, 2, 3)
        >>> m
        [1, 0, 0]
        [0, 2, 0]
        [0, 0, 3]
        >>> m.is_diagonal()
        True

        See Also
        ========

        is_lower
        is_upper
        is_diagonalizable
        diagonalize
        """
        for i in xrange(self.rows):
            for j in xrange(self.cols):
                if i != j and self[i, j] != 0:
                    return False
        return True

    def clone(self):
        """
        Create a shallow copy of this matrix.
        """
        return self._new(self.rows, self.cols, lambda i, j: self[i, j])

    def det(self, method="bareis"):
        """
        Computes the matrix determinant using the method "method".

        Possible values for "method":
          bareis ... det_bareis
          berkowitz ... berkowitz_det
          lu_decomposition ... det_LU

        See Also
        ========

        det_bareis
        berkowitz_det
        det_LU
        """

        # if methods were made internal and all determinant calculations
        # passed through here, then these lines could be factored out of
        # the method routines
        if not self.is_square:
            raise NonSquareMatrixError()
        if not self:
            return S.One
        if method == "bareis":
            return self.det_bareis()
        elif method == "berkowitz":
            return self.berkowitz_det()
        elif method == "det_LU":
            return self.det_LU_decomposition()
        else:
            raise ValueError("Determinant method unrecognized")

    def det_bareis(self):
        """Compute matrix determinant using Bareis' fraction-free
        algorithm which is an extension of the well known Gaussian
        elimination method. This approach is best suited for dense
        symbolic matrices and will result in a determinant with
        minimal number of fractions. It means that less term
        rewriting is needed on resulting formulae.

        TODO: Implement algorithm for sparse matrices (SFF).

        See Also
        ========

        det
        berkowitz_det
        """
        if not self.is_square:
            raise NonSquareMatrixError()
        if not self:
            return S.One

        M, n = self[:,:], self.rows

        if n == 1:
            det = M[0, 0]
        elif n == 2:
            det = M[0, 0]*M[1, 1] - M[0, 1]*M[1, 0]
        else:
            sign = 1 # track current sign in case of column swap

            for k in range(n-1):
                # look for a pivot in the current column
                # and assume det == 0 if none is found
                if M[k, k] == 0:
                    for i in range(k+1, n):
                        if M[i, k] != 0:
                            M.row_swap(i, k)
                            sign *= -1
                            break
                    else:
                        return S.Zero

                # proceed with Bareis' fraction-free (FF)
                # form of Gaussian elimination algorithm
                for i in range(k+1, n):
                    for j in range(k+1, n):
                        D = M[k, k]*M[i, j] - M[i, k]*M[k, j]

                        if k > 0:
                            D /= M[k-1, k-1]

                        if D.is_Atom:
                            M[i, j] = D
                        else:
                            M[i, j] = cancel(D)

            det = sign * M[n-1, n-1]

        return det.expand()

    def det_LU_decomposition(self):
        """Compute matrix determinant using LU decomposition

        Note that this method fails if the LU decomposition itself
        fails. In particular, if the matrix has no inverse this method
        will fail.

        TODO: Implement algorithm for sparse matrices (SFF).

        See Also
        ========

        det
        det_bareis
        berkowitz_det
        """
        if not self.is_square:
            raise NonSquareMatrixError()
        if not self:
            return S.One

        M, n = self[:,:], self.rows
        p, prod = [] , 1
        l, u, p = M.LUdecomposition()
        if  len(p)%2 != 0:
            prod = -1

        for k in range(n):
            prod = prod*u[k,k]*l[k,k]

        return prod.expand()

    def adjugate(self, method="berkowitz"):
        """
        Returns the adjugate matrix.

        Adjugate matrix is the transpose of the cofactor matrix.

        http://en.wikipedia.org/wiki/Adjugate

        See Also
        ========

        cofactorMatrix
        transpose
        berkowitz
        """

        return self.cofactorMatrix(method).T


    def inverse_LU(self, iszerofunc=_iszero):
        """
        Calculates the inverse using LU decomposition.

        See Also
        ========

        inv
        inverse_GE
        inverse_ADJ
        """
        if not self.is_square:
            raise NonSquareMatrixError()

        ok = self.rref(simplify=True)[0]
        if any(iszerofunc(ok[j, j]) for j in range(ok.rows)):
            raise ValueError("Matrix det == 0; not invertible.")

        return self.LUsolve(self.eye(self.rows), iszerofunc=_iszero)

    def inverse_GE(self, iszerofunc=_iszero):
        """
        Calculates the inverse using Gaussian elimination.

        See Also
        ========

        inv
        inverse_LU
        inverse_ADJ
        """
        if not self.is_square:
            raise NonSquareMatrixError()

        big = self.row_join(self.eye(self.rows))
        red = big.rref(iszerofunc=iszerofunc, simplify=True)[0]
        if any(iszerofunc(red[j, j]) for j in range(red.rows)):
            raise ValueError("Matrix det == 0; not invertible.")

        return red[:,big.rows:]

    def inverse_ADJ(self, iszerofunc=_iszero):
        """
        Calculates the inverse using the adjugate matrix and a determinant.

        See Also
        ========

        inv
        inverse_LU
        inverse_GE
        """
        if not self.is_square:
            raise NonSquareMatrixError()

        d = self.berkowitz_det()
        zero = d.equals(0)
        if zero is None:
            # if equals() can't decide, will rref be able to?
            ok = self.rref(simplify=True)[0]
            zero = any(iszerofunc(ok[j, j]) for j in range(ok.rows))
        if zero:
            raise ValueError("Matrix det == 0; not invertible.")

        return self.adjugate()/d

    def rref(self, simplified=False, iszerofunc=_iszero,
            simplify=False):
        """
        Return reduced row-echelon form of matrix and indices of pivot vars.

        To simplify elements before finding nonzero pivots set simplify=True
        (to use the default SymPy simplify function) or pass a custom
        simplify function.

        >>> from sympy import Matrix
        >>> from sympy.abc import x
        >>> m = Matrix([[1, 2], [x, 1 - 1/x]])
        >>> m.rref()
        ([1, 0]
        [0, 1], [0, 1])
        """
        if simplified is not False:
            SymPyDeprecationWarning(
            feature="'simplified' as a keyword to rref",
            useinstead="simplify=True or set simplify equal to your " +
                       "own custom simplification function"
            ).warn()
            simplify = simplify or True
        simpfunc = simplify if isinstance(simplify, FunctionType) else _simplify
        pivot, r = 0, self[:,:].as_mutable()        # pivot: index of next row to contain a pivot
        pivotlist = []                  # indices of pivot variables (non-free)
        for i in range(r.cols):
            if pivot == r.rows:
                break
            if simplify:
                r[pivot,i] = simpfunc(r[pivot,i])
            if iszerofunc(r[pivot,i]):
                for k in range(pivot, r.rows):
                    if simplify and k > pivot:
                        r[k,i] = simpfunc(r[k,i])
                    if not iszerofunc(r[k,i]):
                        break
                if k == r.rows - 1 and iszerofunc(r[k,i]):
                    continue
                r.row_swap(pivot,k)
            scale = r[pivot,i]
            r.row(pivot, lambda x, _: x/scale)
            for j in range(r.rows):
                if j == pivot:
                    continue
                scale = r[j,i]
                r.row(j, lambda x, k: x - scale*r[pivot,k])
            pivotlist.append(i)
            pivot += 1
        return self._new(r), pivotlist

    def nullspace(self, simplified=False, simplify=False):
        """
        Returns list of vectors (Matrix objects) that span nullspace of self
        """
        if simplified is not False:
            SymPyDeprecationWarning(
            feature="'simplified' as a keyword to rref",
            useinstead="simplify=True or set simplify equal to your " +
                       "own custom simplification function"
            ).warn()
            simplify = simplify or True
        simpfunc = simplify if isinstance(simplify, FunctionType) else _simplify
        reduced, pivots = self.rref(simplify=simpfunc)

        basis = []
        # create a set of vectors for the basis
        for i in range(self.cols - len(pivots)):
            basis.append(zeros(self.cols, 1))
        # contains the variable index to which the vector corresponds
        basiskey, cur = [-1]*len(basis), 0
        for i in range(self.cols):
            if i not in pivots:
                basiskey[cur] = i
                cur += 1
        for i in range(self.cols):
            if i not in pivots: # free var, just set vector's ith place to 1
                basis[basiskey.index(i)][i,0] = 1
            else:               # add negative of nonpivot entry to corr vector
                for j in range(i+1, self.cols):
                    line = pivots.index(i)
                    v = reduced[line, j]
                    if simplify:
                        v = simpfunc(v)
                    if v != 0:
                        if j in pivots:
                            # XXX: Is this the correct error?
                            raise NotImplementedError("Could not compute the nullspace of `self`.")
                        basis[basiskey.index(j)][i,0] = -v
        return [self._new(b) for b in basis]

    def berkowitz(self):
        """The Berkowitz algorithm.

           Given N x N matrix with symbolic content, compute efficiently
           coefficients of characteristic polynomials of 'self' and all
           its square sub-matrices composed by removing both i-th row
           and column, without division in the ground domain.

           This method is particularly useful for computing determinant,
           principal minors and characteristic polynomial, when 'self'
           has complicated coefficients e.g. polynomials. Semi-direct
           usage of this algorithm is also important in computing
           efficiently sub-resultant PRS.

           Assuming that M is a square matrix of dimension N x N and
           I is N x N identity matrix,  then the following following
           definition of characteristic polynomial is begin used:

                          charpoly(M) = det(t*I - M)

           As a consequence, all polynomials generated by Berkowitz
           algorithm are monic.

           >>> from sympy import Matrix
           >>> from sympy.abc import x, y, z

           >>> M = Matrix([[x,y,z], [1,0,0], [y,z,x]])

           >>> p, q, r = M.berkowitz()

           >>> p # 1 x 1 M's sub-matrix
           (1, -x)

           >>> q # 2 x 2 M's sub-matrix
           (1, -x, -y)

           >>> r # 3 x 3 M's sub-matrix
           (1, -2*x, x**2 - y*z - y, x*y - z**2)

           For more information on the implemented algorithm refer to:

           [1] S.J. Berkowitz, On computing the determinant in small
               parallel time using a small number of processors, ACM,
               Information Processing Letters 18, 1984, pp. 147-150

           [2] M. Keber, Division-Free computation of sub-resultants
               using Bezout matrices, Tech. Report MPI-I-2006-1-006,
               Saarbrucken, 2006

        See Also
        ========

        berkowitz_det
        berkowitz_minors
        berkowitz_charpoly
        berkowitz_eigenvals
        """
        if not self.is_square:
            raise NonSquareMatrixError()

        A, N = self, self.rows
        transforms = [0] * (N-1)

        for n in xrange(N, 1, -1):
            T, k = zeros(n+1, n), n - 1

            R, C = -A[k,:k], A[:k,k]
            A, a = A[:k,:k], -A[k,k]

            items = [ C ]

            for i in xrange(0, n-2):
                items.append(A * items[i])

            for i, B in enumerate(items):
                items[i] = (R * B)[0,0]

            items = [ S.One, a ] + items

            for i in xrange(n):
                T[i:,i] = items[:n-i+1]

            transforms[k-1] = T

        polys = [ self._new([S.One, -A[0,0]]) ]

        for i, T in enumerate(transforms):
            polys.append(T * polys[i])

        return tuple(map(tuple, polys))

    def berkowitz_det(self):
        """Computes determinant using Berkowitz method.

        See Also
        ========

        det
        berkowitz
        """
        if not self.is_square:
            raise NonSquareMatrixError()
        if not self:
            return S.One
        poly = self.berkowitz()[-1]
        sign = (-1)**(len(poly)-1)
        return sign * poly[-1]

    def berkowitz_minors(self):
        """Computes principal minors using Berkowitz method.

        See Also
        ========

        berkowitz
        """
        sign, minors = S.NegativeOne, []

        for poly in self.berkowitz():
            minors.append(sign*poly[-1])
            sign = -sign

        return tuple(minors)

    def berkowitz_charpoly(self, x=Dummy('lambda'), simplify=_simplify):
        """Computes characteristic polynomial minors using Berkowitz method.

        A PurePoly is returned so using different variables for ``x`` does
        not affect the comparison or the polynomials:

        >>> from sympy import Matrix
        >>> from sympy.abc import x, y
        >>> A = Matrix([[1, 3], [2, 0]])
        >>> A.berkowitz_charpoly(x) == A.berkowitz_charpoly(y)
        True

        Specifying ``x`` is optional; a Dummy with name ``lambda`` is used by
        default (which looks good when pretty-printed in unicode):

        >>> A.berkowitz_charpoly().as_expr()
        _lambda**2 - _lambda - 6

        No test is done to see that ``x`` doesn't clash with an existing
        symbol, so using the default (``lambda``) or your own Dummy symbol is
        the safest option:

        >>> A = Matrix([[1, 2], [x, 0]])
        >>> A.charpoly().as_expr()
        _lambda**2 - _lambda - 2*x
        >>> A.charpoly(x).as_expr()
        x**2 - 3*x

        See Also
        ========

        berkowitz
        """
        return PurePoly(map(simplify, self.berkowitz()[-1]), x)

    charpoly = berkowitz_charpoly

    def berkowitz_eigenvals(self, **flags):
        """Computes eigenvalues of a Matrix using Berkowitz method.

        See Also
        ========

        berkowitz
        """
        return roots(self.berkowitz_charpoly(Dummy('x')), **flags)

    def eigenvals(self, **flags):
        """Return eigen values using the berkowitz_eigenvals routine.

        Since the roots routine doesn't always work well with Floats,
        they will be replaced with Rationals before calling that
        routine. If this is not desired, set flag ``rational`` to False.
        """
        # roots doesn't like Floats, so replace them with Rationals
        # unless the nsimplify flag indicates that this has already
        # been done, e.g. in eigenvects
        if flags.pop('rational', True):
            float = False
            if any(v.has(Float) for v in self):
                float=True
                self = self._new(self.rows, self.cols, [nsimplify(v, rational=True) for v in self])

        flags.pop('simplify', None) # pop unsupported flag
        return self.berkowitz_eigenvals(**flags)

    def eigenvects(self, **flags):
        """Return list of triples (eigenval, multiplicity, basis).

        The flag ``simplify`` has two effects:
            1) if bool(simplify) is True, as_content_primitive()
            will be used to tidy up normalization artifacts;
            2) if nullspace needs simplification to compute the
            basis, the simplify flag will be passed on to the
            nullspace routine which will interpret it there.

        If the matrix contains any Floats, they will be changed to Rationals
        for computation purposes, but the answers will be returned after being
        evaluated with evalf. If it is desired to removed small imaginary
        portions during the evalf step, pass a value for the ``chop`` flag.
        """

        simplify = flags.get('simplify', True)
        primitive = bool(flags.get('simplify', False))
        chop = flags.pop('chop', False)

        flags.pop('multiple', None) # remove this if it's there

        # roots doesn't like Floats, so replace them with Rationals
        float = False
        if any(v.has(Float) for v in self):
            float = True
            self = self._new(self.rows, self.cols, [nsimplify(v, rational=True) for v in self])
            flags['rational'] = False # to tell eigenvals not to do this

        out, vlist = [], self.eigenvals(**flags)
        vlist = vlist.items()
        vlist.sort(key=default_sort_key)
        flags.pop('rational', None)

        for r, k in vlist:
            tmp = self - eye(self.rows)*r
            basis = tmp.nullspace()
            # whether tmp.is_symbolic() is True or False, it is possible that
            # the basis will come back as [] in which case simplification is
            # necessary.
            if not basis:
                # The nullspace routine failed, try it again with simplification
                basis = tmp.nullspace(simplify=simplify)
                if not basis:
                    raise NotImplementedError("Can't evaluate eigenvector for eigenvalue %s" % r)
            if primitive:
                # the relationship A*e = lambda*e will still hold if we change the
                # eigenvector; so if simplify is True we tidy up any normalization
                # artifacts with as_content_primtive (default) and remove any pure Integer
                # denominators.
                l = 1
                for i, b in enumerate(basis[0]):
                    c, p = signsimp(b).as_content_primitive()
                    if c is not S.One:
                        b = c*p
                        l = ilcm(l, c.q)
                    basis[0][i] = b
                if l != 1:
                    basis[0] *= l
            if float:
                out.append((r.evalf(chop=chop), k, [self._new(b).evalf(chop=chop) for b in basis]))
            else:
                out.append((r, k, [self._new(b) for b in basis]))
        return out

    def singular_values(self):
        """
        Compute the singular values of a Matrix

        >>> from sympy import Matrix, Symbol, eye
        >>> x = Symbol('x', real=True)
        >>> A = Matrix([[0, 1, 0], [0, x, 0], [-1, 0, 0]])
        >>> A.singular_values()
        [sqrt(x**2 + 1), 1, 0]

        See Also
        ========

        condition_number
        """
        self = self.as_mutable()
        # Compute eigenvalues of A.H A
        valmultpairs = (self.H*self).eigenvals()

        # Expands result from eigenvals into a simple list
        vals = []
        for k,v in valmultpairs.items():
            vals += [sqrt(k)]*v # dangerous! same k in several spots!
        # sort them in descending order
        vals.sort(reverse=True, key=default_sort_key)

        return vals

    def condition_number(self):
        """
        Returns the condition number of a matrix.

        This is the maximum singular value divided by the minimum singular value

        >>> from sympy import Matrix, S
        >>> A = Matrix([[1, 0, 0], [0, 10, 0], [0,0,S.One/10]])
        >>> A.condition_number()
        100

        See Also
        ========

        singular_values
        """

        singularvalues = self.singular_values()
        return Max(*singularvalues) / Min(*singularvalues)

    def __getattr__(self, attr):
        if attr in ('diff','integrate','limit'):
            def doit(*args):
                item_doit = lambda item: getattr(item, attr)(*args)
                return self.applyfunc( item_doit )
            return doit
        else:
            raise AttributeError("Matrix has no attribute %s." % attr)

    def integrate(self, *args):
        """
        Integrate each element of the matrix.

        Examples
        ========

        >>> import sympy
        >>> from sympy.abc import x, y
        >>> M = sympy.matrices.Matrix([[x, y], [1, 0]])
        >>> M.integrate((x,))
        [x**2/2, x*y]
        [     x,   0]
        >>> M.integrate((x, 0, 2))
        [2, 2*y]
        [2,   0]

        See Also
        ========

        limit
        diff
        """
        return self._new(self.rows, self.cols,
                lambda i, j: self[i, j].integrate(*args))

    def limit(self, *args):
        """
        Calculate the limit of each element in the matrix.

        Examples
        ========

        >>> import sympy
        >>> from sympy.abc import x, y
        >>> M = sympy.matrices.Matrix([[x, y], [1, 0]])
        >>> M.limit(x, 2)
        [2, y]
        [1, 0]

        See Also
        ========

        integrate
        diff
        """
        return self._new(self.rows, self.cols,
                lambda i, j: self[i, j].limit(*args))

    def diff(self, *args):
        """
        Calculate the derivative of each element in the matrix.

        Examples
        ========

        >>> import sympy
        >>> from sympy.abc import x, y
        >>> M = sympy.matrices.Matrix([[x, y], [1, 0]])
        >>> M.diff(x)
        [1, 0]
        [0, 0]

        See Also
        ========

        integrate
        limit
        """
        return self._new(self.rows, self.cols,
                lambda i, j: self[i, j].diff(*args))

    def vec(self):
        """
        Return the Matrix converted into a one column matrix by stacking columns

        >>> from sympy import Matrix
        >>> m=Matrix([[1,3], [2,4]])
        >>> m
        [1, 3]
        [2, 4]
        >>> m.vec()
        [1]
        [2]
        [3]
        [4]

        See Also
        ========

        vech
        """
        return self.__class__(len(self), 1, self.transpose().mat)

    def vech(self, diagonal=True, check_symmetry=True):
        """
        Return the unique elements of a symmetric Matrix as a one column matrix
        by stacking the elements in the lower triangle.

        Arguments:
        diagonal -- include the diagonal cells of self or not
        check_symmetry -- checks symmetry of self but not completely reliably

        >>> from sympy import Matrix
        >>> m=Matrix([[1,2], [2,3]])
        >>> m
        [1, 2]
        [2, 3]
        >>> m.vech()
        [1]
        [2]
        [3]
        >>> m.vech(diagonal=False)
        [2]

        See Also
        ========

        vec
        """
        c = self.cols
        if c != self.rows:
            raise ShapeError("Matrix must be square")
        if check_symmetry:
            self.simplify()
            if self != self.transpose():
                raise ValueError("Matrix appears to be asymmetric; consider check_symmetry=False")
        count = 0
        if diagonal:
            v = zeros(c * (c + 1) // 2, 1)
            for j in xrange(c):
                for i in xrange(j,c):
                    v[count] = self[i,j]
                    count += 1
        else:
            v = zeros(c * (c - 1) // 2, 1)
            for j in xrange(c):
                for i in xrange(j+1,c):
                    v[count] = self[i,j]
                    count += 1
        return v

    def get_diag_blocks(self):
        """Obtains the square sub-matrices on the main diagonal of a square matrix.

        Useful for inverting symbolic matrices or solving systems of
        linear equations which may be decoupled by having a block diagonal
        structure.

        Examples
        ========

        >>> from sympy import Matrix, symbols
        >>> from sympy.abc import x, y, z
        >>> A = Matrix([[1, 3, 0, 0], [y, z*z, 0, 0], [0, 0, x, 0], [0, 0, 0, 0]])
        >>> a1, a2, a3 = A.get_diag_blocks()
        >>> a1
        [1,    3]
        [y, z**2]
        >>> a2
        [x]
        >>> a3
        [0]
        >>>

        """
        sub_blocks = []
        def recurse_sub_blocks(M):
            i = 1
            while i <= M.shape[0]:
                if i == 1:
                    to_the_right = M[0, i:]
                    to_the_bottom = M[i:, 0]
                else:
                    to_the_right = M[:i, i:]
                    to_the_bottom = M[i:, :i]
                if any(to_the_right) or any(to_the_bottom):
                    i += 1
                    continue
                else:
                    sub_blocks.append(M[:i, :i])
                    if M.shape == M[:i, :i].shape:
                        return
                    else:
                        recurse_sub_blocks(M[i:, i:])
                        return
        recurse_sub_blocks(self)
        return sub_blocks

    def diagonalize(self, reals_only = False):
        """
        Return diagonalized matrix D and transformation P such as

            D = P^-1 * M * P

        where M is current matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(3,3,[1, 2, 0, 0, 3, 0, 2, -4, 2])
        >>> m
        [1,  2, 0]
        [0,  3, 0]
        [2, -4, 2]
        >>> (P, D) = m.diagonalize()
        >>> D
        [1, 0, 0]
        [0, 2, 0]
        [0, 0, 3]
        >>> P
        [-1, 0, -1]
        [ 0, 0, -1]
        [ 2, 1,  2]
        >>> P.inv() * m * P
        [1, 0, 0]
        [0, 2, 0]
        [0, 0, 3]

        See Also
        ========

        is_diagonal
        is_diagonalizable
        """
        if not self.is_square:
            raise NonSquareMatrixError()
        if not self.is_diagonalizable(reals_only, False):
            self._diagonalize_clear_subproducts()
            raise MatrixError("Matrix is not diagonalizable")
        else:
            if self._eigenvects is None:
                self._eigenvects = self.eigenvects(simplify=True)
            diagvals = []
            P = MutableMatrix(self.rows, 0, [])
            for eigenval, multiplicity, vects in self._eigenvects:
                for k in range(multiplicity):
                    diagvals.append(eigenval)
                    vec = vects[k]
                    P = P.col_insert(P.cols, vec)
            D = diag(*diagvals)
            self._diagonalize_clear_subproducts()
            return (P, D)

    def is_diagonalizable(self, reals_only = False, clear_subproducts=True):
        """
        Check if matrix is diagonalizable.

        If reals_only==True then check that diagonalized matrix consists of the only not complex values.

        Some subproducts could be used further in other methods to avoid double calculations,
        By default (if clear_subproducts==True) they will be deleted.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(3,3,[1, 2, 0, 0, 3, 0, 2, -4, 2])
        >>> m
        [1,  2, 0]
        [0,  3, 0]
        [2, -4, 2]
        >>> m.is_diagonalizable()
        True
        >>> m = Matrix(2,2,[0, 1, 0, 0])
        >>> m
        [0, 1]
        [0, 0]
        >>> m.is_diagonalizable()
        False
        >>> m = Matrix(2,2,[0, 1, -1, 0])
        >>> m
        [ 0, 1]
        [-1, 0]
        >>> m.is_diagonalizable()
        True
        >>> m.is_diagonalizable(True)
        False

        See Also
        ========

        is_diagonal
        diagonalize
        """
        if not self.is_square:
            return False
        res = False
        self._is_symbolic = self.is_symbolic()
        self._is_symmetric = self.is_symmetric()
        self._eigenvects = None
        #if self._is_symbolic:
        #    self._diagonalize_clear_subproducts()
        #    raise NotImplementedError("Symbolic matrices are not implemented for diagonalization yet")
        self._eigenvects = self.eigenvects(simplify=True)
        all_iscorrect = True
        for eigenval, multiplicity, vects in self._eigenvects:
            if len(vects) != multiplicity:
                all_iscorrect = False
                break
            elif reals_only and not eigenval.is_real:
                all_iscorrect = False
                break
        res = all_iscorrect
        if clear_subproducts:
            self._diagonalize_clear_subproducts()
        return res

    def _diagonalize_clear_subproducts(self):
        del self._is_symbolic
        del self._is_symmetric
        del self._eigenvects

    def jordan_form(self, calc_transformation = True):
        """
        Return Jordan form J of current matrix.

        If calc_transformation is specified as False, then transformation P such that

              J = P^-1 * M * P

        will not be calculated.

        Notes
        =====

        Calculation of transformation P is not implemented yet.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(4, 4, [6, 5, -2, -3, -3, -1, 3, 3, 2, 1, -2, -3, -1, 1, 5, 5])
        >>> m
        [ 6,  5, -2, -3]
        [-3, -1,  3,  3]
        [ 2,  1, -2, -3]
        [-1,  1,  5,  5]

        >>> (P, J) = m.jordan_form()
        >>> J
        [2, 1, 0, 0]
        [0, 2, 0, 0]
        [0, 0, 2, 1]
        [0, 0, 0, 2]

        See Also
        ========

        jordan_cells
        """
        (P, Jcells) = self.jordan_cells(calc_transformation)
        J = diag(*Jcells)
        return (P, J)

    def jordan_cells(self, calc_transformation = True):
        """
        Return a list of Jordan cells of current matrix.
        This list shape Jordan matrix J.

        If calc_transformation is specified as False, then transformation P such that

              J = P^-1 * M * P

        will not be calculated.

        Notes
        =====

        Calculation of transformation P is not implemented yet.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(4, 4, [6, 5, -2, -3, -3, -1, 3, 3, 2, 1, -2, -3, -1, 1, 5, 5])
        >>> m
        [ 6,  5, -2, -3]
        [-3, -1,  3,  3]
        [ 2,  1, -2, -3]
        [-1,  1,  5,  5]

        >>> (P, Jcells) = m.jordan_cells()
        >>> Jcells[0]
        [2, 1]
        [0, 2]
        >>> Jcells[1]
        [2, 1]
        [0, 2]

        See Also
        ========

        jordan_form
        """
        if not self.is_square:
            raise NonSquareMatrixError()
        _eigenvects = self.eigenvects()
        Jcells = []
        for eigenval, multiplicity, vects in _eigenvects:
            geometrical = len(vects)
            if geometrical == multiplicity:
                Jcell = diag( *([eigenval] * multiplicity))
                Jcells.append(Jcell)
            elif geometrical==0:
                raise MatrixError("Matrix has the eigen vector with geometrical multiplicity equal zero.")
            else:
                sizes = self._jordan_split(multiplicity, geometrical)
                cells = []
                for size in sizes:
                    cell = jordan_cell(eigenval, size)
                    cells.append(cell)
                Jcells += cells
        return (None, Jcells)

    def _jordan_split(self, algebraical, geometrical):
        """return a list of integers with sum equal to 'algebraical'
        and length equal to 'geometrical'"""
        n1 = algebraical // geometrical
        res = [n1] * geometrical
        res[len(res) - 1] += algebraical % geometrical
        assert sum(res) == algebraical
        return res

    def has(self, *patterns):
        """
        Test whether any subexpression matches any of the patterns.

        Examples
        ========

        >>> from sympy import Matrix, Float
        >>> from sympy.abc import x, y
        >>> A = Matrix(((1, x), (0.2, 3)))
        >>> A.has(x)
        True
        >>> A.has(y)
        False
        >>> A.has(Float)
        True
        """
        return any(a.has(*patterns) for a in self.mat)

    @property
    def is_Identity(self):
        if not self.is_square:
            return False
        for i in xrange(self.rows):
            for j in xrange(self.cols):
                if i==j and self[i,j] != 1:
                    return False
                if i!=j and self[i,j] != 0:
                    return False
        return True

    def dual(self):
        """
        Returns the dual of a matrix, which is:

        `(1/2)*levicivita(i,j,k,l)*M(k,l)` summed over indices `k` and `l`

        Since the levicivita method is anti_symmetric for any pairwise
        exchange of indices, the dual of a symmetric matrix is the zero
        matrix. Strictly speaking the dual defined here assumes that the
        'matrix' `M` is a contravariant anti_symmetric second rank tensor,
        so that the dual is a covariant second rank tensor.

        """
        from sympy import LeviCivita

        M, n = self[:,:], self.rows
        work = zeros(n)
        if self.is_symmetric():
            return work

        for i in range(1, n):
            for j in range(1, n):
                acum = 0
                for k in range(1, n):
                    acum += LeviCivita(i, j, 0, k)*M[0, k]
                work[i, j] = acum
                work[j, i] = -acum

        for l in range(1,n):
            acum = 0
            for a in range(1,n):
                for b in range(1,n):
                    acum += LeviCivita(0, l, a, b)*M[a, b]
            acum /= 2
            work[0, l] = -acum
            work[l, 0] = acum

        return work

class MutableMatrix(MatrixBase):

    is_MatrixExpr = False

    _op_priority = 12.0
    _class_priority = 10

    @classmethod
    def _new(cls, *args, **kwargs):
        rows, cols, mat = MatrixBase._handle_creation_inputs(*args, **kwargs)
        self = object.__new__(cls)
        self.rows = rows
        self.cols = cols
        self.mat = list(mat) # create a shallow copy
        return self
    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    def __setitem__(self, key, value):
        """
        >>> from sympy import Matrix, I
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        [1, 2 + I]
        [3,     4]
        >>> m[1,0]=9
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        [1, 2 + I]
        [9,     4]

        """
        if type(key) is tuple:
            i, j = key
            if type(i) is slice or type(j) is slice:
                if isinstance(value, MatrixBase):
                    self.copyin_matrix(key, value)
                    return
                if is_sequence(value):
                    self.copyin_list(key, value)
                    return
            else:
                # a2idx inlined
                if type(i) is not int:
                    try:
                        i = i.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % (key,))

                # a2idx inlined
                if type(j) is not int:
                    try:
                        j = j.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % (key,))


                i, j = self.key2ij((i,j))
                if not (i>=0 and i<self.rows and j>=0 and j < self.cols):
                    raise IndexError("Index out of range: a[%s]" % (key,))
                else:
                    self.mat[i*self.cols + j] = sympify(value)
                    return

        else:
            # row-wise decomposition of matrix
            if type(key) is slice:
                raise IndexError("Vector slices not implemented yet.")
            else:
                k = a2idx(key)
                if k is not None:
                    self.mat[k] = sympify(value)
                    return
        raise IndexError("Invalid index: a[%s]"%repr(key))

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

        >>> import sympy
        >>> M = sympy.matrices.Matrix([[0,1],[2,3]])
        >>> I = sympy.matrices.eye(5)
        >>> I.copyin_matrix((slice(0,2), slice(0,2)),M)
        >>> I
        [0, 1, 0, 0, 0]
        [2, 3, 0, 0, 0]
        [0, 0, 1, 0, 0]
        [0, 0, 0, 1, 0]
        [0, 0, 0, 0, 1]

        See Also
        ========

        copyin_list
        """
        rlo, rhi, clo, chi = self.key2bounds(key)

        if value.rows != rhi - rlo or value.cols != chi - clo:
            raise ShapeError("The Matrix `value` doesn't have the same dimensions " +
                "as the in sub-Matrix given by `key`.")

        for i in range(value.rows):
            for j in range(value.cols):
                self[i+rlo, j+clo] = sympify(value[i,j])

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

        >>> import sympy
        >>> I = sympy.matrices.eye(5)
        >>> I.copyin_list((slice(0,2), slice(0,1)), [1,2])
        >>> I
        [1, 0, 0, 0, 0]
        [2, 1, 0, 0, 0]
        [0, 0, 1, 0, 0]
        [0, 0, 0, 1, 0]
        [0, 0, 0, 0, 1]

        See Also
        ========

        copyin_matrix
        """
        if not is_sequence(value):
            raise TypeError("`value` must be an ordered iterable, not %s." % type(value))
        return self.copyin_matrix(key, MutableMatrix(value))

    def row(self, i, f=None):
        """
        Elementary row selector (default) or operation using functor
        which is a function two args interpreted as (self[i, j], j).

        >>> from sympy import ones
        >>> I = ones(3)
        >>> I.row(1, lambda v,i: v*3)
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
        >>> M = Matrix([[0,1],[1,0]])
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
        >>> M = Matrix([[1,0],[1,0]])
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

        >>> import sympy
        >>> M = sympy.matrices.eye(3)
        >>> M.row_del(1)
        >>> M   #doctest: +NORMALIZE_WHITESPACE
        [1, 0, 0]
        [0, 0, 1]

        See Also
        ========

        row
        col_del
        """
        self.mat = self.mat[:i*self.cols] + self.mat[(i+1)*self.cols:]
        self.rows -= 1

    def col_del(self, i):
        """
        Delete the given column.

        >>> import sympy
        >>> M = sympy.matrices.eye(3)
        >>> M.col_del(1)
        >>> M   #doctest: +NORMALIZE_WHITESPACE
        [1, 0]
        [0, 0]
        [0, 1]

        See Also
        ========

        col
        row_del
        """
        for j in range(self.rows-1, -1, -1):
            del self.mat[i+j*self.cols]
        self.cols -= 1

    # Utility functions
    def simplify(self, ratio=1.7, measure=count_ops):
        """Applies simplify to the elements of a matrix in place.

        This is a shortcut for M.applyfunc(lambda x: simplify(x, ratio, measure))

        See Also
        ========

        sympy.simplify.simplify.simplify
        """
        for i in xrange(len(self.mat)):
            self.mat[i] = _simplify(self.mat[i], ratio=ratio, measure=measure)

    def fill(self, value):
        """Fill the matrix with the scalar value.

        See Also
        ========

        zeros
        ones
        """
        self.mat = [value]*len(self)

    ############################
    # Mutable matrix operators #
    ############################

    @call_highest_priority('__radd__')
    def __add__(self, other):
        return MatrixBase.__add__(self, force_mutable(other))

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return MatrixBase.__radd__(self, force_mutable(other))

    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return MatrixBase.__sub__(self, force_mutable(other))
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return MatrixBase.__rsub__(self, force_mutable(other))

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return MatrixBase.__mul__(self, force_mutable(other))
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return MatrixBase.__rmul__(self, force_mutable(other))

    @call_highest_priority('__div__')
    def __div__(self, other):
        return MatrixBase.__div__(self, force_mutable(other))
    @call_highest_priority('__truediv__')
    def __truediv__(self, other):
        return MatrixBase.__truediv__(self, force_mutable(other))

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return MatrixBase.__pow__(self, other)

    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        raise NotImplementedError("Matrix Power not defined")

def force_mutable(x):
    """
    Safely turn a Matrix object into a Mutable Matrix
    """
    if not isinstance(x, Basic):
        return x
    if not x.is_Matrix:
        return x
    if x.is_MatrixExpr:
        return x.as_mutable()
    return x

def classof(A,B):
    """
    Determines strategy for combining Immutable and Mutable matrices

    Currently the strategy is that Mutability is contagious

    Example

    >>> from sympy import Matrix, ImmutableMatrix
    >>> from sympy.matrices.matrices import classof
    >>> M = Matrix([[1,2],[3,4]]) # a Mutable Matrix
    >>> IM = ImmutableMatrix([[1,2],[3,4]])
    >>> classof(M, IM)
    <class 'sympy.matrices.matrices.Matrix'>
    """
    from immutable_matrix import ImmutableMatrix
    try:
        if A._class_priority > B._class_priority:
            return A.__class__
        else:
            return B.__class__
    except: pass
    try:
        import numpy
        if isinstance(A, numpy.ndarray):
            return B.__class__
        if isinstance(B, numpy.ndarray):
            return A.__class__
    except: pass
    raise TypeError("Incompatible classes %s, %s"%(A.__class__, B.__class__))

def matrix_multiply(A, B):
    """
    Matrix product A*B.

    A and B must be of appropriate dimensions.  If A is an m x k matrix, and B
    is a k x n matrix, the product will be an m x n matrix.

    Examples
    ========

    >>> from sympy import Matrix
    >>> A = Matrix([[1, 2, 3], [4, 5, 6]])
    >>> B = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> A*B
    [30, 36, 42]
    [66, 81, 96]
    >>> B*A
    Traceback (most recent call last):
    ...
    ShapeError: Matrices size mismatch.
    >>>

    See Also
    ========

    matrix_multiply_elementwise
    matrix_add
    """
    # The following implmentation is equivalent, but about 5% slower
    #ma, na = A.shape
    #mb, nb = B.shape
    #
    #if na != mb:
    #    raise ShapeError()
    #product = Matrix(ma, nb, lambda i,j: 0)
    #for i in xrange(ma):
    #    for j in xrange(nb):
    #        s = 0
    #        for k in range(na):
    #            s += A[i, k]*B[k, j]
    #        product[i, j] = s
    #return product
    if A.shape[1] != B.shape[0]:
        raise ShapeError("Matrices size mismatch.")
    blst = B.T.tolist()
    alst = A.tolist()
    return classof(A,B)._new(A.shape[0], B.shape[1], lambda i, j:
                                        reduce(lambda k, l: k+l,
                                        map(lambda n, m: n*m,
                                        alst[i],
                                        blst[j]), 0))

def matrix_multiply_elementwise(A, B):
    """Return the Hadamard product (elementwise product) of A and B

    >>> from sympy.matrices.matrices import Matrix, matrix_multiply_elementwise
    >>> A = Matrix([[0, 1, 2], [3, 4, 5]])
    >>> B = Matrix([[1, 10, 100], [100, 10, 1]])
    >>> matrix_multiply_elementwise(A, B)
    [  0, 10, 200]
    [300, 40,   5]

    See Also
    ========

    matrix_multiply
    matrix_add
    """
    if A.shape != B.shape:
        raise ShapeError()
    shape = A.shape
    return classof(A,B)._new(shape[0], shape[1],
        lambda i, j: A[i,j] * B[i, j])

def matrix_add(A, B):
    """Return A + B

    See Also
    ========

    matrix_multiply
    matrix_multiply_elementwise
    """
    if A.shape != B.shape:
        raise ShapeError()
    alst = A.tolist()
    blst = B.tolist()
    ret = [0]*A.shape[0]
    for i in xrange(A.shape[0]):
        ret[i] = map(lambda j,k: j+k, alst[i], blst[i])
    return classof(A,B)._new(ret)

def zeros(r, c=None, cls=MutableMatrix):
    """Returns a matrix of zeros with ``r`` rows and ``c`` columns;
    if ``c`` is omitted a square matrix will be returned.

    See Also
    ========

    ones
    eye
    diag
    """
    if is_sequence(r):
        SymPyDeprecationWarning(
        "Pass row and column count as zeros(%i, %i)." % r
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

    if is_sequence(r):
        SymPyDeprecationWarning(
        "Pass row and column count as ones(%i, %i)." % r
        ).warn()
        r, c = r
    else:
        c = r if c is None else c
    r, c = [int(i) for i in [r, c]]
    return MutableMatrix(r, c, [S.One]*r*c)

def eye(n, cls=MutableMatrix):
    """Create square identity matrix n x n

    See Also
    ========

    diag
    zeros
    ones
    """
    n = int(n)
    out = cls.zeros(n)
    for i in range(n):
        out[i, i] = S.One
    return out

def diag(*values):
    """Create diagonal matrix from a list as a diagonal values.

    Arguments might be matrices too, in case of it they are fitted in result matrix

    Examples
    ========

    >>> from sympy.matrices import diag, Matrix
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

    See Also
    ========

    eye
    """
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

    >>> from sympy.matrices.matrices import jordan_cell
    >>> from sympy.abc import x
    >>> jordan_cell(x, 4)
    [x, 1, 0, 0]
    [0, x, 1, 0]
    [0, 0, x, 1]
    [0, 0, 0, x]
    """
    n = int(n)
    out = zeros(n)
    for i in range(n-1):
        out[i, i] = eigenval
        out[i, i+1] = S.One
    out[n-1, n-1] = eigenval
    return out

def randMatrix(r, c=None, min=0, max=99, seed=None, symmetric=False):
    """Create random matrix with dimensions ``r`` x ``c``. If ``c`` is omitted
    the matrix will be square. If ``symmetric`` is True the matrix must be
    square.

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
        return MutableMatrix(r, c, lambda i, j: prng.randint(min, max))
    m = zeros(r)
    for i in xrange(r):
        for j in xrange(i, r):
            m[i, j] = prng.randint(min, max)
    for i in xrange(r):
        for j in xrange(i):
            m[i, j] = m[j, i]
    return m

def hessian(f, varlist):
    """Compute Hessian matrix for a function f

    see: http://en.wikipedia.org/wiki/Hessian_matrix

    See Also
    ========

    sympy.matrices.matrices.Matrix.jacobian
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
        for j in range(i,m):
            out[i,j] = f.diff(varlist[i]).diff(varlist[j])
    for i in range(m):
        for j in range(i):
            out[i,j] = out[j,i]
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
        if tmp == MutableMatrix([[0,0,0]]):
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
        W(f1,...,fn) = |  .        .         .     .      |
                       |  .        .          .    .      |
                       |  (n)      (n)            (n)     |
                       | D   (f1) D   (f2)  ...  D   (fn) |

    see: http://en.wikipedia.org/wiki/Wronskian

    See Also
    ========

    sympy.matrices.matrices.Matrix.jacobian
    hessian
    """

    for index in xrange(0, len(functions)):
        functions[index] = sympify(functions[index])
    n = len(functions)
    if n == 0:
        return 1
    W = MutableMatrix(n, n, lambda i,j: functions[i].diff(var, j) )
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
    seqs = map(sympify, seqs)

    if not zero:
        f = lambda i, j: seqs[j].subs(n, n+i)
    else:
        f = lambda i, j: seqs[j].subs(n, i)

    k = len(seqs)

    return MutableMatrix(k, k, f).det()



class SparseMatrix(MatrixBase):
    """
    A sparse matrix (a matrix with a large number of zero elements).

    See Also
    ========

    :class:`Matrix`
    """

    def __init__(self, *args):
        self.mat = {}
        if len(args) == 3 and callable(args[2]):
            op = args[2]
            if not isinstance(args[0], (int, Integer)) or not isinstance(args[1], (int, Integer)):
                raise TypeError("`args[0]` and `args[1]` must both be integers.")
            self.rows = args[0]
            self.cols = args[1]
            for i in range(self.rows):
                for j in range(self.cols):
                    value = sympify(op(i,j))
                    if value != 0:
                        self.mat[(i,j)] = value
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and is_sequence(args[2]):
            self.rows = args[0]
            self.cols = args[1]
            mat = args[2]
            for i in range(self.rows):
                for j in range(self.cols):
                    value = sympify(mat[i*self.cols+j])
                    if value != 0:
                        self.mat[(i,j)] = value
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], dict):
            self.rows = args[0]
            self.cols = args[1]
            # manual copy, copy.deepcopy() doesn't work
            for key in args[2].keys():
                self.mat[key] = args[2][key]
        else:
            # TODO: _handle_creation_inputs creates a temporary dense array
            # and calls sympify on each element. This is a waste of time.
            # Just like above rewrite the rest of the handled cases with
            # "sparse" logic.
            r, c, mat = self._handle_creation_inputs(*args)
            self.rows = r
            self.cols = c
            for i in range(self.rows):
                for j in range(self.cols):
                    value = mat[self.cols*i+j]
                    if value != 0:
                        self.mat[(i,j)] = value

    def __getitem__(self, key):
        if isinstance(key, slice) or isinstance(key, int):
            lo, hi = self.slice2bounds(key, len(self))
            L = []
            for i in range(lo, hi):
                m,n = self.rowdecomp(i)
                if (m,n) in self.mat:
                    L.append(self.mat[(m,n)])
                else:
                    L.append(0)
            if len(L) == 1:
                return L[0]
            else:
                return L
        if len(key) != 2:
            raise ValueError("`key` must be of length 2.")

        if isinstance(key[0], int) and isinstance(key[1], int):
            i,j=self.key2ij(key)
            if (i, j) in self.mat:
                return self.mat[(i,j)]
            else:
                return 0
        elif isinstance(key[0], slice) or isinstance(key[1], slice):
            return self.submatrix(key)
        else:
            raise IndexError("Index out of range: a[%s]"%repr(key))

    def rowdecomp(self, num):
        """
        Perform a row decomposition on the matrix.
        """
        nmax = len(self)
        if not (0 <= num < nmax) or not (0 <= -num < nmax):
            raise ValueError("`num` must satisfy 0 <= `num` < `self.rows*" +
                "*self.cols` (%d) and 0 <= -num < " % nmax +
                "`self.rows*self.cols` (%d) to apply redecomp()." % nmax)
        i, j = 0, num
        while j >= self.cols:
            j -= self.cols
            i += 1
        return i,j

    def __setitem__(self, key, value):
        # almost identical, need to test for 0
        if len(key) != 2:
            raise ValueError("`key` must be of length 2.")

        if isinstance(key[0], slice) or isinstance(key[1], slice):
            if isinstance(value, Matrix):
                self.copyin_matrix(key, value)
            if is_sequence(value):
                self.copyin_list(key, value)
        else:
            i,j=self.key2ij(key)
            testval = sympify(value)
            if testval != 0:
                self.mat[(i,j)] = testval
            elif (i,j) in self.mat:
                del self.mat[(i,j)]

    def row_del(self, k):
        """
        Delete the given row of the matrix.

        Examples
        ========

        >>> import sympy
        >>> M = sympy.matrices.SparseMatrix([[0,0],[0,1]])
        >>> M
        [0, 0]
        [0, 1]
        >>> M.row_del(0)
        >>> M
        [0, 1]

        See Also
        ========

        col_del
        """
        newD = {}
        k = self.key2ij((k,0))[0]
        for (i,j) in self.mat.keys():
            if i==k:
                pass
            elif i > k:
                newD[i-1,j] = self.mat[i,j]
            else:
                newD[i,j] = self.mat[i,j]
        self.mat = newD
        self.rows -= 1

    def col_del(self, k):
        """
        Delete the given column of the matrix.

        Examples
        ========

        >>> import sympy
        >>> M = sympy.matrices.SparseMatrix([[0,0],[0,1]])
        >>> M
        [0, 0]
        [0, 1]
        >>> M.col_del(0)
        >>> M
        [0]
        [1]

        See Also
        ========

        row_del
        """
        newD = {}
        k = self.key2ij((0,k))[1]
        for (i,j) in self.mat.keys():
            if j==k:
                pass
            elif j > k:
                newD[i,j-1] = self.mat[i,j]
            else:
                newD[i,j] = self.mat[i,j]
        self.mat = newD
        self.cols -= 1

    def toMatrix(self):
        """
        Convert this sparse matrix into a matrix.
        """
        l = []
        for i in range(self.rows):
            c = []
            l.append(c)
            for j in range(self.cols):
                if (i, j) in self.mat:
                    c.append(self[i, j])
                else:
                    c.append(0)
        return Matrix(l)

    def row_list(self):
        """
        Returns a Row-sorted list of non-zero elements of the matrix.

        >>> from sympy.matrices import SparseMatrix
        >>> a=SparseMatrix(((1,2),(3,4)))
        >>> a
        [1, 2]
        [3, 4]
        >>> a.RL
        [(0, 0, 1), (0, 1, 2), (1, 0, 3), (1, 1, 4)]

        See Also
        ========

        col_list
        """

        new=[]
        for i in range(self.rows):
            for j in range(self.cols):
                value = self[(i,j)]
                if value!=0:
                    new.append((i,j,value))
        return new

    RL = property(row_list,None,None,"Alternate faster representation")

    def col_list(self):
        """
        Returns a Column-sorted list of non-zero elements of the matrix.
        >>> from sympy.matrices import SparseMatrix
        >>> a=SparseMatrix(((1,2),(3,4)))
        >>> a
        [1, 2]
        [3, 4]
        >>> a.CL
        [(0, 0, 1), (1, 0, 3), (0, 1, 2), (1, 1, 4)]

        See Also
        ========

        row_list
        """
        new=[]
        for j in range(self.cols):
            for i in range(self.rows):
                value = self[(i,j)]
                if value!=0:
                    new.append((i,j,value))
        return new

    CL = property(col_list,None,None,"Alternate faster representation")

    def transpose(self):
        """
        Returns the transposed SparseMatrix of this SparseMatrix
        >>> from sympy.matrices import SparseMatrix
        >>> a = SparseMatrix(((1,2),(3,4)))
        >>> a
        [1, 2]
        [3, 4]
        >>> a.T
        [1, 3]
        [2, 4]
        """
        tran = SparseMatrix(self.cols,self.rows,{})
        for key,value in self.mat.iteritems():
            tran.mat[key[1],key[0]]=value
        return tran

    T = property(transpose,None,None,"Matrix transposition.")


    def __add__(self, other):
        if isinstance(other, SparseMatrix):
            return self.add(other)
        else:
            raise NotImplementedError("Only SparseMatrix + SparseMatrix supported")

    def __radd__(self, other):
        if isinstance(other, SparseMatrix):
            return self.add(other)
        else:
            raise NotImplementedError("Only SparseMatrix + SparseMatrix supported")

    def add(self, other):
        """
        Add two sparse matrices with dictionary representation.

        >>> from sympy.matrices.matrices import SparseMatrix
        >>> A = SparseMatrix(5, 5, lambda i, j: i * j + i)
        >>> A
        [0, 0,  0,  0,  0]
        [1, 2,  3,  4,  5]
        [2, 4,  6,  8, 10]
        [3, 6,  9, 12, 15]
        [4, 8, 12, 16, 20]
        >>> B = SparseMatrix(5, 5, lambda i, j: i + 2 * j)
        >>> B
        [0, 2, 4,  6,  8]
        [1, 3, 5,  7,  9]
        [2, 4, 6,  8, 10]
        [3, 5, 7,  9, 11]
        [4, 6, 8, 10, 12]
        >>> A + B
        [0,  2,  4,  6,  8]
        [2,  5,  8, 11, 14]
        [4,  8, 12, 16, 20]
        [6, 11, 16, 21, 26]
        [8, 14, 20, 26, 32]

        See Also
        ========

        multiply
        """
        if self.shape != other.shape:
            raise ShapeError()
        a, b = self.mat.keys(), other.mat.keys()
        a.sort()
        b.sort()
        i = j = 0
        c = {}
        while i < len(a) or j < len(b):
            if j >= len(b) or (i < len(a) and a[i] < b[j]):
                c[a[i]] = self.mat[a[i]]
                i = i + 1
                continue
            elif i >= len(a) or (j < len(b) and a[i] > b[j]):
                c[b[j]] = other.mat[b[j]]
                j = j + 1
                continue
            else:
                c[a[i]] = self.mat[a[i]] + other.mat[b[j]]
                i = i + 1
                j = j + 1
        return SparseMatrix(self.rows, self.cols, c)



    # from here to end all functions are same as in matrices.py
    # with Matrix replaced with SparseMatrix
    def copyin_list(self, key, value):
        if not is_sequence(value):
            raise TypeError("`value` must be of type list or tuple.")
        self.copyin_matrix(key, SparseMatrix(value))

    def multiply(self,b):
        """Returns self*b

        See Also
        ========

        add
        """

        def dotprod(a,b,i,j):
            if a.cols != b.rows:
                raise ShapeError("`self.cols` must equal `b.rows`.")
            r=0
            for x in range(a.cols):
                r+=a[i,x]*b[x,j]
            return r

        r = SparseMatrix(self.rows, b.cols, lambda i,j: dotprod(self,b,i,j))
        if r.rows == 1 and r.cols ==1:
            return r[0,0]
        return r

    def submatrix(self, keys):
        rlo, rhi, clo, chi = self.key2bounds(keys)
        return SparseMatrix(rhi-rlo, chi-clo, lambda i,j: self[i+rlo, j+clo])

    def reshape(self, _rows, _cols):
        if len(self) != _rows*_cols:
            raise ValueError("Invalid reshape parameters %d %d" % (_rows, _cols))
        newD = {}
        for i in range(_rows):
            for j in range(_cols):
                m,n = self.rowdecomp(i*_cols + j)
                if (m,n) in self.mat:
                    newD[(i,j)] = self.mat[(m,n)]
        return SparseMatrix(_rows, _cols, newD)

    def cross(self, b):
        if not is_sequence(b, include=Matrix):
            raise TypeError("`b` must be an ordered iterable or Matrix, not %s." %
                type(b))
        if not (self.rows == 1 and self.cols == 3 or \
                self.rows == 3 and self.cols == 1 ) and \
                (b.rows == 1 and b.cols == 3 or \
                b.rows == 3 and b.cols == 1):
            raise ShapeError("Dimensions incorrect for cross product")
        else:
            return SparseMatrix(1,3,((self[1]*b[2] - self[2]*b[1]),
                               (self[2]*b[0] - self[0]*b[2]),
                               (self[0]*b[1] - self[1]*b[0])))


    @classmethod
    def zeros(cls, r, c=None):
        """Returns a matrix of zeros with ``r`` rows and ``c`` columns;
        if ``c`` is omitted a square matrix will be returned."""
        if is_sequence(r):
            SymPyDeprecationWarning(
            "Pass row and column count as zeros(%i, %i)." % r
            ).warn()
            r, c = r
        else:
            c = r if c is None else c
        r, c = [int(i) for i in [r, c]]
        return cls(r, c, {})

    @classmethod
    def eye(cls, n):
        n = int(n)
        tmp = cls.zeros(n)
        for i in range(tmp.rows):
            tmp[i,i] = 1
        return tmp

    def __hash__(self):
        return super(Matrix, self).__hash__()

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

def a2idx(a):
    """
    Tries to convert "a" to an index, returns None on failure.

    The result of a2idx() (if not None) can be safely used as an index to
    arrays/matrices.
    """
    if hasattr(a, "__int__"):
        return int(a)
    if hasattr(a, "__index__"):
        return a.__index__()

def symarray(prefix, shape):
    """Create a numpy ndarray of symbols (as an object array).

    The created symbols are named `prefix_i1_i2_`...  You should thus provide a
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
    mat = ((ct,st,0),
           (-st,ct,0),
           (0,0,1))
    return MutableMatrix(mat)

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
    mat = ((ct,0,-st),
           (0,1,0),
           (st,0,ct))
    return MutableMatrix(mat)

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
    mat = ((1,0,0),
           (0,ct,st),
           (0,-st,ct))
    return MutableMatrix(mat)

Matrix = MutableMatrix
Matrix.__name__ = "Matrix"
