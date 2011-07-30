#from __future__ import division

import matrixutils
import densematrix_tools
from datamatrix import DataMatrix
from sympy.core.sympify import sympify as sympy_simplify

class DenseMatrix(DataMatrix):

    def __init__(self, *args, **kwargs):
        """
        DenseMatrix, internal low-level matrix for Matrix
        self.rows ---> number of rows
        self.cols ---> number of cols
        self.mat  ---> data stored in a single array of size rows * cols,
                       with element (i, j) in mat[i*cols + j]
        """
        
        if len(args) == 3:
            rows = args[0]
            cols = args[1]
            mat = args[2]

            if isinstance(mat, dict):
                mat = matrixutils._denserepr_from_dict(rows, cols, mat)
            elif callable(mat):
                mat = matrixutils._denserepr_from_callable(rows, cols, mat)
            elif isinstance(mat, (list, tuple)) and len(mat) == rows * cols:
                mat = matrixutils._denserepr_from_list(rows, cols, mat)
            else:
                raise TypeError('Data type not understood.')
        elif len(args) == 1:
            mat = args[0]
            if isinstance(mat, DenseMatrix):
                rows = mat.rows
                cols = mat.cols
                mat = matrixutils._denserepr_from_list(rows, cols, mat.mat)
            elif isinstance(mat, (list, tuple)):
                rows = len(mat)
                cols = len(mat[0])
                mat = matrixutils._denserepr_from_lil(rows, cols, mat)
            else:
                raise TypeError('Data type not understood.')
        elif len(args) == 0:
            rows = cols = 0
            mat = []

        self.rows = rows
        self.cols = cols
        self.mat = mat

    @property
    def T(self):
        """
        matrix transposition.

        >>> from sympy import DenseMatrix, I
        >>> m=DenseMatrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        [1, 2 + I]
        [3,     4]
        >>> m.transpose() #doctest: +NORMALIZE_WHITESPACE
        [    1, 3]
        [2 + I, 4]
        >>> m.T == m.transpose()
        True

        """
        return DenseMatrix(self.cols, self.rows, lambda i, j: self[j, i])

    def conjugate(self):
        """By-element conjugation."""
        return DenseMatrix(self.rows, self.cols, lambda i, j: self[i,j].conjugate())

    @property
    def H(self):
        """
        Hermite conjugation.

        >>> from sympy import DenseMatrix, I
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        [1, 2 + I]
        [3,     4]
        >>> m.H #doctest: +NORMALIZE_WHITESPACE
        [    1, 3]
        [2 - I, 4]

        """
        return DenseMatrix(self.cols, self.rows, lambda i, j: self[j,i].conjugate())

    def __getitem__(self,key):
        """
        >>> from sympy import DenseMatrix, I
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        [1, 2 + I]
        [3,     4]
        >>> m[1,0]
        3
        >>> m.H[1,0]
        2 - I

        """
        i, j = key
        try:
            return self.mat[i*self.cols + j]
        except:
            return self.submatrix(key)

    def __setitem__(self, key, value):
        """
        >>> from sympy import DenseMatrix, I
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        [1, 2 + I]
        [3,     4]
        >>> m[1,0]=9
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        [1, 2 + I]
        [9,     4]

        """
        i, j = key
        try:
            self.mat[i*self.cols + j] = value
        except:
            if isinstance(value, DenseMatrix):
                self.copyin_matrix(key, value)
            elif isinstance(value, list):
                self.copyin_list(key, value)

    def tolist(self):
        """
        Return the matrix converted in a python list.

        >>> from sympy import MDenseatrix
        >>> m=Matrix(3, 3, range(9))
        >>> m
        [0, 1, 2]
        [3, 4, 5]
        [6, 7, 8]
        >>> m.tolist()
        [[0, 1, 2], [3, 4, 5], [6, 7, 8]]

        """
        li = [0] * self.rows
        for i in xrange(self.rows):
            li[i] = self.mat[i * self.cols : (i + 1) * self.cols]
        return li

    def copyin_matrix(self, key, value):
        rlo, rhi = _slice_to_bounds(key[0], self.rows)
        clo, chi = _slice_to_bounds(key[1], self.cols)

        for i in xrange(value.rows):
            for j in xrange(value.cols):
                self[i + rlo, j + clo] = value[i, j]

    def copyin_list(self, key, value):
        self.copyin_matrix(key, DenseMatrix(value)) # dense list of list repr used here

    def __add__(A, B):
        return DenseMatrix(A.rows, A.cols, lambda i, j: A[i, j] + B[i, j])

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        return -1 * self

    def __mul__(A, B):
        if isinstance(B, DenseMatrix):
            return DenseMatrix(A.rows, B.cols, lambda i, j: 
                sum(A[i, k] * B[k, j] for k in xrange(A.cols)))
        return DenseMatrix(A.rows, A.cols, lambda i, j: A[i, j] * B)

    __rmul__ = __mul__

    def __div__(self, a):
        return self * (1 / a)

    __truediv__ = __div__

    def __pow__(self, num):
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

    def __eq__(self, other):
        if self.shape != other.shape:
            return False
        return all(self[i, j] == other[i, j]
            for i in xrange(self.rows) for j in xrange(self.cols))

    def __ne__(self, other):
        if self.shape != other.shape:
            return True
        return any(self[i, j] != other[i, j]
            for i in xrange(self.rows) for j in xrange(self.cols))

    def to_densematrix(self):
        return self

    def to_dokmatrix(self):
        from sympy.matrices import DOKMatrix
        return DOKMatrix(self.rows, self.cols, lambda i, j: self[i, j])

    def to_lilmatrix(self):
        return LILMatrix(self.rows, self.cols, lambda i, j: self[i, j])

    def inv(self, method=None):
        """
        Calculates the matrix inverse.

        According to the "method" parameter, it calls the appropriate method:

          GE .... inverse_GE()
          LU .... inverse_LU()
          ADJ ... inverse_ADJ()
          LDL ... inverse_solver(method='LDL')

        According to the "try_block_diag" parameter, it will try to form block
        diagonal matrices using the method get_diag_blocks(), invert these
        individually, and then reconstruct the full inverse matrix.

        Note, the GE and LU methods may require the matrix to be simplified
        before it is inverted in order to properly detect zeros during
        pivoting. In difficult cases a custom zero detection function can
        be provided by setting the iszerosfunc argument to a function that
        should return True if its argument is zero.

        """
        if not method:
            method = 'GE'
        if not self.is_square:
            raise NonSquareMatrixError()
        if method == "GE":
            return densematrix_tools.inverse_GE(self)
        elif method == "LU":
            return densematrix_tools.inverse_LU(self)
        elif method == "ADJ":
            return denseatrix_tools.inverse_ADJ(self)
        elif method == "LDL" or method == "CH":
            return densematrix_tools.inverse_solver(self, solver=method)
        else:
            raise ValueError("Inversion method unrecognized")

    def row(self, i, f):
        """
        Elementary row operation using functor

        >>> from sympy import ones
        >>> I = ones(3)
        >>> I.row(1,lambda i,j: i*3)
        >>> I
        [1, 1, 1]
        [3, 3, 3]
        [1, 1, 1]

        """
        for j in range(0, self.cols):
            self[i, j] = f(self[i, j], j)

    def col(self, j, f):
        """
        Elementary column operation using functor

        >>> from sympy import ones
        >>> I = ones(3)
        >>> I.col(0,lambda i,j: i*3)
        >>> I
        [3, 1, 1]
        [3, 1, 1]
        [3, 1, 1]

        """
        for i in range(0, self.rows):
            self[i, j] = f(self[i, j], i)

    def row_insert(self, pos, mti):
        """
        >>> from sympy import DenseMatrix, zeros
        >>> M = DenseMatrix(3,3,lambda i,j: i+j)
        >>> M
        [0, 1, 2]
        [1, 2, 3]
        [2, 3, 4]
        >>> V = zeros((1, 3))
        >>> V
        [0, 0, 0]
        >>> M.row_insert(1,V)
        [0, 1, 2]
        [0, 0, 0]
        [1, 2, 3]
        [2, 3, 4]

        """
        if pos is 0:
            return mti.col_join(self)

        newmat = self.zeros((self.rows + mti.rows, self.cols))
        newmat[:pos,:] = self[:pos,:]
        newmat[pos:pos+mti.rows,:] = mti[:,:]
        newmat[pos+mti.rows:,:] = self[pos:,:]
        return newmat

    def col_insert(self, pos, mti):
        """
        >>> from sympy import DenseMatrix, zeros
        >>> M = DenseMatrix(3,3,lambda i,j: i+j)
        >>> M
        [0, 1, 2]
        [1, 2, 3]
        [2, 3, 4]
        >>> V = zeros((3, 1))
        >>> V
        [0]
        [0]
        [0]
        >>> M.col_insert(1,V)
        [0, 0, 1, 2]
        [1, 0, 2, 3]
        [2, 0, 3, 4]

        """
        if pos is 0:
            return mti.row_join(self)

        newmat = self.zeros((self.rows, self.cols + mti.cols))
        newmat[:,:pos] = self[:,:pos]
        newmat[:,pos:pos+mti.cols] = mti[:,:]
        newmat[:,pos+mti.cols:] = self[:,pos:]
        return newmat

    def trace(self):
        return sum(self[i, i] for i in xrange(self.rows))

    def submatrix(self, *args):
        """
        >>> from sympy import DenseMatrix
        >>> m = DenseMatrix(4,4,lambda i,j: i+j)
        >>> m   #doctest: +NORMALIZE_WHITESPACE
        [0, 1, 2, 3]
        [1, 2, 3, 4]
        [2, 3, 4, 5]
        [3, 4, 5, 6]
        >>> m[0:1, 1]   #doctest: +NORMALIZE_WHITESPACE
        [1]
        >>> m[0:2, 0:1] #doctest: +NORMALIZE_WHITESPACE
        [0]
        [1]
        >>> m[2:4, 2:4] #doctest: +NORMALIZE_WHITESPACE
        [4, 5]
        [5, 6]

        """
        if len(args) == 1:
            from matrixutils import _slice_to_bounds
            keys = args[0]
            rlo, rhi = _slice_to_bounds(keys[0], self.rows)
            clo, chi = _slice_to_bounds(keys[1], self.cols)
        elif len(args) == 4:
            rlo, rhi, clo, chi = args
        else:
            raise TypeError("Data type not understood")
        subrows, subcols = rhi-rlo, chi-clo
        return DenseMatrix(subrows, subcols, lambda i, j: self[i + rlo, j + clo])
    
    def applyfunc(self, f):
        """
        >>> from sympy import DenseMatrix
        >>> m = DenseMatrix(2,2,lambda i,j: i*2+j)
        >>> m   #doctest: +NORMALIZE_WHITESPACE
        [0, 1]
        [2, 3]
        >>> m.applyfunc(lambda i: 2*i)  #doctest: +NORMALIZE_WHITESPACE
        [0, 2]
        [4, 6]

        """
        return DenseMatrix(self.rows, self.cols, lambda i, j: f(self[i, j]))

    def evalf(self, prec=15, **options):
        return self.applyfunc(lambda i: i.evalf(prec, **options))

    def simplify(self, simplify=sympy_simplify, ratio=1.7):
        """Simplify the elements of a matrix in place.

        If (result length)/(input length) > ratio, then input is returned
        unmodified. If 'ratio=oo', then simplify() is applied anyway.

        See also simplify().
        """
        for i in xrange(len(self.mat)):
            self.mat[i] = simplify(self.mat[i], ratio=ratio)

    def copy(self):
        return DenseMatrix(self.rows, self.cols, lambda i, j: self[i, j])

    # This could be removed from the class
    def det(self, method="bareis"):
        """
        Computes the matrix determinant using the method "method".

        Possible values for "method":
          bareis ... det_bareis
          berkowitz ... berkowitz_det
        """

        if method == "bareis":
            return densematrix_tools.det_bareis(self)
        elif method == "berkowitz":
            return densematrix_tools.berkowitz_det(self)
        else:
            raise ValueError("Determinant method unrecognized")

    def fill(self, value):
        """Fill the matrix with the scalar value."""
        self.mat = [value] * len(self)

    def __getattr__(self, attr):
        if hasattr(self[0, 0], attr):
            def doit(*args):
                item_doit = lambda item: getattr(item, attr)(*args)
                return self.applyfunc( item_doit )
            return doit
        else:
            raise AttributeError("DenseMatrix has no attribute %s." % attr)

    def has(self, *patterns):
        """
        Test whether any subexpression matches any of the patterns.

        Examples:
        >>> from sympy import DenseMatrix, Float
        >>> from sympy.abc import x, y
        >>> A = DenseMatrix(((1, x), (0.2, 3)))
        >>> A.has(x)
        True
        >>> A.has(y)
        False
        >>> A.has(Float)
        True
        """
        return any(a.has(*patterns) for a in self.mat)

