import warnings
from sympy.core.basic import Basic
from sympy.core.symbol import Symbol, Dummy
from sympy.core.numbers import Integer
from sympy.core.singleton import S
from sympy.core.sympify import sympify, converter, SympifyError

from sympy.polys import Poly, roots, cancel
from sympy.simplify import simplify
from sympy.utilities import any, all
from sympy.printing import sstr


import random

class NonSquareMatrixException(Exception):
    pass

class ShapeError(ValueError):
    """Wrong matrix shape"""
    pass

class MatrixError(Exception):
    pass

def _dims_to_nm(dims):
    """Converts dimensions tuple (or any object with length 1 or 2) or scalar
    in dims to matrix dimensions n and m."""

    try:
        l = len(dims)
    except TypeError:
        dims = (dims,)
        l = 1

    # This will work for nd-array too when they are added to sympy.
    try:
        for dim in dims:
            assert (dim > 0)
    except AssertionError:
        raise ValueError("Matrix dimensions should be positive integers!")

    if l == 2:
        n, m = map(int, dims)
    elif l == 1:
        n = m = int(dims[0])
    else:
        raise ValueError("Matrix dimensions should be a two-element tuple of ints or a single int!")

    return n, m

def _iszero(x):
    return x == 0


class DeferredVector(object):
    def __init__(self,name):
        self.name=name

    def __getitem__(self,i):
        component_name = '%s[%d]'%(self.name,i)
        return Symbol(component_name)

    def __str__(self):
        return sstr(self)

    def __repr__(self):
        return sstr(self)


class Matrix(object):

    # Added just for numpy compatibility
    # TODO: investigate about __array_priority__
    __array_priority__ = 10.0

    def __init__(self, *args):
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
        if len(args) == 3 and callable(args[2]):
            operation = args[2]
            self.rows = int(args[0])
            self.cols = int(args[1])
            self.mat = []
            for i in range(self.rows):
                for j in range(self.cols):
                    self.mat.append(sympify(operation(i, j)))
        elif len(args)==3 and isinstance(args[2], (list, tuple)):
            self.rows=args[0]
            self.cols=args[1]
            mat = args[2]
            if len(mat) != self.rows*self.cols:
                raise MatrixError('List length should be equal to rows*columns')
            self.mat = map(lambda i: sympify(i), mat)
        elif len(args) == 1:
            mat = args[0]
            if isinstance(mat, Matrix):
                self.rows = mat.rows
                self.cols = mat.cols
                self.mat = mat[:]
                return
            elif hasattr(mat, "__array__"):
                # NumPy array or matrix or some other object that implements
                # __array__. So let's first use this method to get a
                # numpy.array() and then make a python list out of it.
                arr = mat.__array__()
                if len(arr.shape) == 2:
                    self.rows, self.cols = arr.shape[0], arr.shape[1]
                    self.mat = map(lambda i: sympify(i), arr.ravel())
                    return
                elif len(arr.shape) == 1:
                    self.rows, self.cols = 1, arr.shape[0]
                    self.mat = [0]*self.cols
                    for i in xrange(len(arr)):
                        self.mat[i] = sympify(arr[i])
                    return
                else:
                    raise NotImplementedError("Sympy supports just 1D and 2D matrices")
            elif not isinstance(mat, (list, tuple, Matrix)):
                raise TypeError("Matrix constructor doesn't accept %s as input" % str(type(mat)))
            mat = []
            for row in args[0]:
                if isinstance(row, Matrix):
                    mat.extend(row.tolist())
                else:
                    mat.append(row)
            self.rows = len(mat)
            if len(mat) != 0:
                if not isinstance(mat[0], (list, tuple)):
                    self.cols = 1
                    self.mat = map(lambda i: sympify(i), mat)
                    return
                self.cols = len(mat[0])
            else:
                self.cols = 0
            self.mat = []
            for j in xrange(self.rows):
                assert len(mat[j])==self.cols
                for i in xrange(self.cols):
                    self.mat.append(sympify(mat[j][i]))
        elif len(args) == 0:
            # Empty Matrix
            self.rows = self.cols = 0
            self.mat = []
        else:
            raise TypeError("Data type not understood")

    def key2ij(self,key):
        """Converts key=(4,6) to 4,6 and ensures the key is correct."""
        if not (isinstance(key,(list, tuple)) and len(key) == 2):
            raise TypeError("wrong syntax: a[%s]. Use a[i,j] or a[(i,j)]"
                    %repr(key))
        i,j=key
        if not (i>=0 and i<self.rows and j>=0 and j < self.cols):
            print self.rows, " ", self.cols
            raise IndexError("Index out of range: a[%s]"%repr(key))
        return i,j

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

        """
        a = [0]*self.cols*self.rows
        for i in xrange(self.cols):
            a[i*self.rows:(i+1)*self.rows] = self.mat[i::self.cols]
        return Matrix(self.cols,self.rows,a)

    T = property(transpose,None,None,"Matrix transposition.")

    def conjugate(self):
        """By-element conjugation."""
        out = Matrix(self.rows,self.cols,
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

        """
        out = self.T.C
        return out

    @property
    def D(self):
        """Dirac conjugation."""
        from sympy.physics.matrices import mgamma
        out = self.H * mgamma(0)
        return out

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
                try:
                    i = i.__int__()
                except AttributeError:
                    try:
                        i = i.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % (key,))

                # a2idx inlined
                try:
                    j = j.__int__()
                except AttributeError:
                    try:
                        j = j.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % (key,))


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
                if isinstance(value, Matrix):
                    self.copyin_matrix(key, value)
                    return
                if isinstance(value, (list, tuple)):
                    self.copyin_list(key, value)
                    return
            else:
                # a2idx inlined
                try:
                    i = i.__int__()
                except AttributeError:
                    try:
                        i = i.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % (key,))

                # a2idx inlined
                try:
                    j = j.__int__()
                except AttributeError:
                    try:
                        j = j.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % (key,))


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

    def __array__(self):
        return matrix2numpy(self)

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

    def copyin_matrix(self, key, value):
        rlo, rhi = self.slice2bounds(key[0], self.rows)
        clo, chi = self.slice2bounds(key[1], self.cols)
        assert value.rows == rhi - rlo and value.cols == chi - clo
        for i in range(value.rows):
            for j in range(value.cols):
                self[i+rlo, j+clo] = sympify(value[i,j])

    def copyin_list(self, key, value):
        assert isinstance(value, (list, tuple))
        self.copyin_matrix(key, Matrix(value))

    def hash(self):
        """Compute a hash every time, because the matrix elements
        could change."""
        return hash(self.__str__() )

    @property
    def shape(self):
        return (self.rows, self.cols)

    def __rmul__(self,a):
        if hasattr(a, "__array__") and a.shape != ():
            return matrix_multiply(a,self)
        out = Matrix(self.rows,self.cols,map(lambda i: a*i,self.mat))
        return out

    def expand(self):
        out = Matrix(self.rows,self.cols,map(lambda i: i.expand(), self.mat))
        return out

    def combine(self):
        out = Matrix(self.rows,self.cols,map(lambda i: i.combine(),self.mat))
        return out

    def subs(self, *args):
        out = Matrix(self.rows,self.cols,map(lambda i: i.subs(*args),self.mat))
        return out

    def __sub__(self,a):
        return self + (-a)

    def __mul__(self,a):
        if hasattr(a, "__array__") and a.shape != ():
            return matrix_multiply(self,a)
        out = Matrix(self.rows,self.cols,map(lambda i: i*a,self.mat))
        return out

    def __pow__(self, num):
        if not self.is_square:
            raise NonSquareMatrixException()
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
        raise NotImplementedError('Can only raise to the power of an integer for now')

    def __add__(self,a):
        return matrix_add(self,a)

    def __radd__(self,a):
        return matrix_add(a,self)

    def __div__(self,a):
        return self * (S.One/a)

    def __truediv__(self,a):
        return self.__div__(a)

    def multiply(self,b):
        """Returns self*b """
        return matrix_multiply(self,b)

    def add(self,b):
        """Return self+b """
        return matrix_add(self,b)

    def __neg__(self):
        return -1*self

    def __eq__(self, a):
        if not isinstance(a, (Matrix, Basic)):
            a = sympify(a)
        if isinstance(a, Matrix):
            return self.hash() == a.hash()
        else:
            return False

    def __ne__(self,a):
        if not isinstance(a, (Matrix, Basic)):
            a = sympify(a)
        if isinstance(a, Matrix):
            return self.hash() != a.hash()
        else:
            return True

    def __hash__(self):
        return super(Matrix, self).__hash__()

    def _format_str(self, strfunc, rowsep='\n'):
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
        should return True if its argument is zero.

        """
        assert self.cols==self.rows
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
            return self.inverse_ADJ()
        else:
            raise ValueError("Inversion method unrecognized")


    def __mathml__(self):
        mml = ""
        for i in range(self.rows):
            mml += "<matrixrow>"
            for j in range(self.cols):
                mml += self[i,j].__mathml__()
            mml += "</matrixrow>"
        return "<matrix>" + mml + "</matrix>"

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

    def row_swap(self, i, j):
        for k in range(0, self.cols):
            self[i, k], self[j, k] = self[j, k], self[i, k]

    def col_swap(self, i, j):
        for k in range(0, self.rows):
            self[k, i], self[k, j] = self[k, j], self[k, i]

    def row_del(self, i):
        self.mat = self.mat[:i*self.cols] + self.mat[(i+1)*self.cols:]
        self.rows -= 1

    def col_del(self, i):
        """
        >>> import sympy
        >>> M = sympy.matrices.eye(3)
        >>> M.col_del(1)
        >>> M   #doctest: +NORMALIZE_WHITESPACE
        [1, 0]
        [0, 0]
        [0, 1]

        """
        for j in range(self.rows-1, -1, -1):
            del self.mat[i+j*self.cols]
        self.cols -= 1

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

        """
        assert self.rows == rhs.rows
        newmat = self.zeros((self.rows, self.cols + rhs.cols))
        newmat[:,:self.cols] = self[:,:]
        newmat[:,self.cols:] = rhs
        return newmat

    def col_join(self, bott):
        """
        Concatenates two matrices along self's last and bott's first row

        >>> from sympy import Matrix
        >>> M = Matrix(3,3,lambda i,j: i+j)
        >>> V = Matrix(1,3,lambda i,j: 3+i+j)
        >>> M.col_join(V)
        [0, 1, 2]
        [1, 2, 3]
        [2, 3, 4]
        [3, 4, 5]

        """
        assert self.cols == bott.cols
        newmat = self.zeros((self.rows+bott.rows, self.cols))
        newmat[:self.rows,:] = self[:,:]
        newmat[self.rows:,:] = bott
        return newmat

    def row_insert(self, pos, mti):
        """
        >>> from sympy import Matrix, zeros
        >>> M = Matrix(3,3,lambda i,j: i+j)
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
        assert self.cols == mti.cols
        newmat = self.zeros((self.rows + mti.rows, self.cols))
        newmat[:pos,:] = self[:pos,:]
        newmat[pos:pos+mti.rows,:] = mti[:,:]
        newmat[pos+mti.rows:,:] = self[pos:,:]
        return newmat

    def col_insert(self, pos, mti):
        """
        >>> from sympy import Matrix, zeros
        >>> M = Matrix(3,3,lambda i,j: i+j)
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
        assert self.rows == mti.rows
        newmat = self.zeros((self.rows, self.cols + mti.cols))
        newmat[:,:pos] = self[:,:pos]
        newmat[:,pos:pos+mti.cols] = mti[:,:]
        newmat[:,pos+mti.cols:] = self[:,pos:]
        return newmat

    def trace(self):
        assert self.cols == self.rows
        trace = 0
        for i in range(self.cols):
            trace += self[i,i]
        return trace

    def submatrix(self, keys):
        """
        >>> from sympy import Matrix
        >>> m = Matrix(4,4,lambda i,j: i+j)
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
        assert isinstance(keys[0], slice) or isinstance(keys[1], slice)
        rlo, rhi = self.slice2bounds(keys[0], self.rows)
        clo, chi = self.slice2bounds(keys[1], self.cols)
        if not ( 0<=rlo<=rhi and 0<=clo<=chi ):
            raise IndexError("Slice indices out of range: a[%s]"%repr(keys))
        outLines, outCols = rhi-rlo, chi-clo
        outMat = [0]*outLines*outCols
        for i in xrange(outLines):
            outMat[i*outCols:(i+1)*outCols] = self.mat[(i+rlo)*self.cols+clo:(i+rlo)*self.cols+chi]
        return Matrix(outLines,outCols,outMat)

    def extract(self, rowsList, colsList):
        """
        Extract a submatrix by specifying a list of rows and columns

        Examples:

        >>> from sympy import Matrix
        >>> m = Matrix(4, 3, lambda i, j: i*3 + j)
        >>> m   #doctest: +NORMALIZE_WHITESPACE
        [0,  1,  2]
        [3,  4,  5]
        [6,  7,  8]
        [9, 10, 11]
        >>> m.extract([0,1,3],[0,1])   #doctest: +NORMALIZE_WHITESPACE
        [0,  1]
        [3,  4]
        [9, 10]

        See also: .submatrix()
        """
        cols = self.cols
        rows = self.rows
        mat = self.mat
        if not all(i < rows for i in rowsList):
            raise IndexError("Row indices out of range")
        if not all(j < cols for j in colsList):
            raise IndexError("Column indices out of range")
        return Matrix(len(rowsList), len(colsList), lambda i,j: mat[rowsList[i]*cols + colsList[j]])

    def slice2bounds(self, key, defmax):
        """
            Takes slice or number and returns (min,max) for iteration
            Takes a default maxval to deal with the slice ':' which is (none, none)
        """
        if isinstance(key, slice):
            lo, hi = 0, defmax
            if key.start is not None:
                if key.start >= 0:
                    lo = key.start
                else:
                    lo = defmax+key.start
            if key.stop is not None:
                if key.stop >= 0:
                    hi = key.stop
                else:
                    hi = defmax+key.stop
            return lo, hi
        elif isinstance(key, int):
            if key >= 0:
                return key, key+1
            else:
                return defmax+key, defmax+key+1
        else:
            raise IndexError("Improper index type")

    def applyfunc(self, f):
        """
        >>> from sympy import Matrix
        >>> m = Matrix(2,2,lambda i,j: i*2+j)
        >>> m   #doctest: +NORMALIZE_WHITESPACE
        [0, 1]
        [2, 3]
        >>> m.applyfunc(lambda i: 2*i)  #doctest: +NORMALIZE_WHITESPACE
        [0, 2]
        [4, 6]

        """
        assert callable(f)
        out = Matrix(self.rows,self.cols,map(f,self.mat))
        return out

    def evalf(self, prec=None, **options):
        if prec is None:
            return self.applyfunc(lambda i: i.evalf(**options))
        else:
            return self.applyfunc(lambda i: i.evalf(prec, **options))

    def reshape(self, _rows, _cols):
        """
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
        if self.rows*self.cols != _rows*_cols:
            print "Invalid reshape parameters %d %d" % (_rows, _cols)
        return Matrix(_rows, _cols, lambda i,j: self.mat[i*_cols + j])

    def print_nonzero (self, symb="X"):
        """
        Shows location of non-zero entries for fast shape lookup ::

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
        s = ""
        for i in range(self.rows):
            s += "["
            for j in range(self.cols):
                if self[i,j] == 0:
                    s += " "
                else:
                    s += symb + ""
            s += "]\n"
        print s

    def LUsolve(self, rhs, iszerofunc=_iszero):
        """
        Solve the linear system Ax = b for x.
        self is the coefficient matrix A and rhs is the right side b.

        This is for symbolic matrices, for real or complex ones use
        sympy.mpmath.lu_solve or sympy.mpmath.qr_solve.
        """
        assert rhs.rows == self.rows
        A, perm = self.LUdecomposition_Simple(iszerofunc=_iszero)
        n = self.rows
        b = rhs.permuteFwd(perm)
        # forward substitution, all diag entries are scaled to 1
        for i in range(n):
            for j in range(i):
                b.row(i, lambda x,k: x - b[j,k]*A[i,j])
        # backward substitution
        for i in range(n-1,-1,-1):
            for j in range(i+1, n):
                b.row(i, lambda x,k: x - b[j,k]*A[i,j])
            b.row(i, lambda x,k: x / A[i,i])
        return b

    def LUdecomposition(self, iszerofunc=_iszero):
        """
        Returns the decomposition LU and the row swaps p.
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
        """
        assert self.rows == self.cols
        n = self.rows
        A = self[:,:]
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
                raise ValueError("Error: non-invertible matrix passed to LUdecomposition_Simple()")
            if pivot != j: # row must be swapped
                A.row_swap(pivot,j)
                p.append([pivot,j])
            assert not iszerofunc(A[j,j])
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
        """
        n, m = self.rows, self.cols
        U, L, P = self[:,:], eye(n), eye(n)
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
        out = Matrix(self.rows, self.cols, lambda i,j:
                self.cofactor(i, j, method))
        return out

    def minorEntry(self, i, j, method="berkowitz"):
        assert 0 <= i < self.rows and 0 <= j < self.cols
        return self.minorMatrix(i,j).det(method)

    def minorMatrix(self, i, j):
        assert 0 <= i < self.rows and 0 <= j < self.cols
        return self.delRowCol(i,j)

    def cofactor(self, i, j, method="berkowitz"):
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

        Examples::

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

        """
        if not isinstance(X, Matrix):
            X = Matrix(X)
        # Both X and self can be a row or a column matrix, so we need to make
        # sure all valid combinations work, but everything else fails:
        assert len(self.shape) == 2
        assert len(X.shape) == 2
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
        return Matrix(m, n, lambda j, i: self[j].diff(X[i]))

    def QRdecomposition(self):
        """
        Return Q,R where A = Q*R, Q is orthogonal and R is upper triangular.

        Assumes full-rank square (for now).
        """
        assert self.rows == self.cols
        n = self.rows
        Q, R = self.zeros(n), self.zeros(n)
        for j in range(n):      # for each column vector
            tmp = self[:,j]     # take original v
            for i in range(j):
                # subtract the project of self on new vector
                tmp -= Q[:,i] * self[:,j].dot(Q[:,i])
                tmp.expand()
            # normalize it
            R[j,j] = tmp.norm()
            Q[:,j] = tmp / R[j,j]
            assert Q[:,j].norm() == 1
            for i in range(j):
                R[i,j] = Q[:,i].dot(self[:,j])
        return Q,R

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
        """

        Q, R = self.QRdecomposition()
        y = Q.T * b

        # back substitution to solve R*x = y:
        # We build up the result "backwards" in the vector 'x' and reverse it
        # only in the end.
        x = []
        n = R.rows
        for j in range(n-1, -1, -1):
            tmp = y[j,:]
            for k in range(j+1, n):
                tmp -= R[j,k] * x[n-1-k]
            x.append(tmp/R[j,j])
        return Matrix([row.mat for row in reversed(x)])

    # Utility functions
    def simplify(self):
        """Simplify the elements of a matrix in place."""
        for i in xrange(len(self.mat)):
            self.mat[i] = simplify(self.mat[i])

    #def evaluate(self):    # no more eval() so should be removed
    #    for i in range(self.rows):
    #        for j in range(self.cols):
    #            self[i,j] = self[i,j].eval()

    def cross(self, b):
        assert isinstance(b, (list, tuple, Matrix))
        if not (self.rows == 1 and self.cols == 3 or \
                self.rows == 3 and self.cols == 1 ) and \
                (b.rows == 1 and b.cols == 3 or \
                b.rows == 3 and b.cols == 1):
            raise ValueError("Dimensions incorrect for cross product")
        else:
            return Matrix(1,3,((self[1]*b[2] - self[2]*b[1]),
                               (self[2]*b[0] - self[0]*b[2]),
                               (self[0]*b[1] - self[1]*b[0])))

    def dot(self, b):
        assert isinstance(b, (list, tuple, Matrix))
        if isinstance(b, (list, tuple)):
            m = len(b)
        else:
            m = b.rows * b.cols
        assert self.cols*self.rows == m
        prod = 0
        for i in range(m):
            prod += self[i] * b[i]
        return prod

    def multiply_elementwise(self, b):
        """Return the Hadamard product (elementwise product) of A and B

        >>> import sympy
        >>> A = sympy.Matrix([[0, 1, 2], [3, 4, 5]])
        >>> B = sympy.Matrix([[1, 10, 100], [100, 10, 1]])
        >>> print A.multiply_elementwise(B)
        [  0, 10, 200]
        [300, 40,   5]
        """
        return matrix_multiply_elementwise(self, b)

    def norm(self):
        assert self.rows == 1 or self.cols == 1
        out = sympify(0)
        for i in range(self.rows * self.cols):
            out += self[i]*self[i]
        return out**S.Half

    def normalized(self):
        assert self.rows == 1 or self.cols == 1
        norm = self.norm()
        out = self.applyfunc(lambda i: i / norm)
        return out

    def project(self, v):
        """Project onto v."""
        return v * (self.dot(v) / v.dot(v))

    def permuteBkwd(self, perm):
        copy = self[:,:]
        for i in range(len(perm)-1, -1, -1):
            copy.row_swap(perm[i][0], perm[i][1])
        return copy

    def permuteFwd(self, perm):
        copy = self[:,:]
        for i in range(len(perm)):
            copy.row_swap(perm[i][0], perm[i][1])
        return copy

    def delRowCol(self, i, j):
        # used only for cofactors, makes a copy
        M = self[:,:]
        M.row_del(i)
        M.col_del(j)
        return M

    def zeronm(self, n, m):
        # used so that certain functions above can use this
        # then only this func need be overloaded in subclasses
        warnings.warn( 'Deprecated: use zeros() instead.' )
        return Matrix(n,m,[S.Zero]*n*m)

    def zero(self, n):
        """Returns a n x n matrix of zeros."""
        warnings.warn( 'Deprecated: use zeros() instead.' )
        return Matrix(n,n,[S.Zero]*n*n)

    def zeros(self, dims):
        """Returns a dims = (d1,d2) matrix of zeros."""
        n, m = _dims_to_nm( dims )
        return Matrix(n,m,[S.Zero]*n*m)

    def eye(self, n):
        """Returns the identity matrix of size n."""
        tmp = self.zeros(n)
        for i in range(tmp.rows):
            tmp[i,i] = S.One
        return tmp

    @property
    def is_square(self):
        return self.rows == self.cols

    def is_upper(self):
        for i in range(self.rows):
            for j in range(self.cols):
                if i > j and self[i,j] != 0:
                    return False
        return True

    def is_lower(self):
        for i in range(self.rows):
            for j in range(self.cols):
                if i < j and self[i, j] != 0:
                    return False
        return True

    def is_symbolic(self):
        for i in range(self.rows):
            for j in range(self.cols):
                if self[i,j].has(Symbol):
                    return True
        return False

    def is_symmetric(self):
        """
        Check if matrix is symmetric matrix,
        that is square matrix and is equal to its transpose.

        Example:

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
        [         1, 1 + 2*x + x**2, y]
        [(1 + x)**2,              2, 0]
        [         y,              0, 3]
        >>> m.is_symmetric()
        True

        """
        if not self.is_square:
            return False
        m = self.clone()
        m.simplify()
        return (m == m.transpose())

    def is_diagonal(self):
        """
        Check if matrix is diagonal,
        that is matrix in which the entries outside the main diagonal are all zero.

        Example:

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

        See also: .is_lower(), is_upper() .is_diagonalizable()
        """
        if not self.is_square:
            return False
        for i in range(self.cols):
            for j in range(self.rows):
                if i <> j and self[i, j] != 0:
                    return False
        return True

    def clone(self):
        return Matrix(self.rows, self.cols, lambda i, j: self[i, j])

    def det(self, method="bareis"):
        """
        Computes the matrix determinant using the method "method".

        Possible values for "method":
          bareis ... det_bareis
          berkowitz ... berkowitz_det
        """

        if method == "bareis":
            return self.det_bareis()
        elif method == "berkowitz":
            return self.berkowitz_det()
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
        """
        if not self.is_square:
            raise NonSquareMatrixException()

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

    def adjugate(self, method="berkowitz"):
        """
        Returns the adjugate matrix.

        Adjugate matrix is the transpose of the cofactor matrix.

        http://en.wikipedia.org/wiki/Adjugate

        See also: .cofactorMatrix(), .T
        """

        return self.cofactorMatrix(method).T


    def inverse_LU(self, iszerofunc=_iszero):
        """
        Calculates the inverse using LU decomposition.
        """
        return self.LUsolve(self.eye(self.rows), iszerofunc=_iszero)

    def inverse_GE(self, iszerofunc=_iszero):
        """
        Calculates the inverse using Gaussian elimination.
        """
        assert self.rows == self.cols
        assert self.det() != 0
        big = self.row_join(self.eye(self.rows))
        red = big.rref(iszerofunc=iszerofunc)
        return red[0][:,big.rows:]

    def inverse_ADJ(self):
        """
        Calculates the inverse using the adjugate matrix and a determinant.
        """
        assert self.rows == self.cols
        d = self.berkowitz_det()
        assert d != 0
        return self.adjugate()/d

    def rref(self,simplified=False, iszerofunc=_iszero):
        """
        Take any matrix and return reduced row-echelon form and indices of pivot vars

        To simplify elements before finding nonzero pivots set simplified=True
        """
        # TODO: rewrite inverse_GE to use this
        pivots, r = 0, self[:,:]        # pivot: index of next row to contain a pivot
        pivotlist = []                  # indices of pivot variables (non-free)
        for i in range(r.cols):
            if pivots == r.rows:
                break
            if simplified:
                r[pivots,i] = simplify(r[pivots,i])
            if iszerofunc(r[pivots,i]):
                for k in range(pivots, r.rows):
                    if simplified and k>pivots:
                        r[k,i] = simplify(r[k,i])
                    if not iszerofunc(r[k,i]):
                        break
                if k == r.rows - 1 and iszerofunc(r[k,i]):
                    continue
                r.row_swap(pivots,k)
            scale = r[pivots,i]
            r.row(pivots, lambda x, _: x/scale)
            for j in range(r.rows):
                if j == pivots:
                    continue
                scale = r[j,i]
                r.row(j, lambda x, k: x - r[pivots,k]*scale)
            pivotlist.append(i)
            pivots += 1
        return r, pivotlist

    def nullspace(self,simplified=False):
        """
        Returns list of vectors (Matrix objects) that span nullspace of self
        """
        reduced, pivots = self.rref(simplified)
        basis = []
        # create a set of vectors for the basis
        for i in range(self.cols - len(pivots)):
            basis.append(zeros((self.cols, 1)))
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
                    if reduced[line, j] != 0:
                        assert j not in pivots
                        basis[basiskey.index(j)][i,0] = -1 * reduced[line, j]
        return basis

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

           >>> M = Matrix([ [x,y,z], [1,0,0], [y,z,x] ])

           >>> p, q, r = M.berkowitz()

           >>> print p # 1 x 1 M's sub-matrix
           (1, -x)

           >>> print q # 2 x 2 M's sub-matrix
           (1, -x, -y)

           >>> print r # 3 x 3 M's sub-matrix
           (1, -2*x, -y - y*z + x**2, x*y - z**2)

           For more information on the implemented algorithm refer to:

           [1] S.J. Berkowitz, On computing the determinant in small
               parallel time using a small number of processors, ACM,
               Information Processing Letters 18, 1984, pp. 147-150

           [2] M. Keber, Division-Free computation of sub-resultants
               using Bezout matrices, Tech. Report MPI-I-2006-1-006,
               Saarbrucken, 2006

        """
        if not self.is_square:
            raise NonSquareMatrixException()

        A, N = self, self.rows
        transforms = [0] * (N-1)

        for n in xrange(N, 1, -1):
            T, k = zeros((n+1,n)), n - 1

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

        polys = [ Matrix([S.One, -A[0,0]]) ]

        for i, T in enumerate(transforms):
            polys.append(T * polys[i])

        return tuple(map(tuple, polys))

    def berkowitz_det(self):
        """Computes determinant using Berkowitz method."""
        poly = self.berkowitz()[-1]
        sign = (-1)**(len(poly)-1)
        return sign * poly[-1]

    def berkowitz_minors(self):
        """Computes principal minors using Berkowitz method."""
        sign, minors = S.NegativeOne, []

        for poly in self.berkowitz():
            minors.append(sign*poly[-1])
            sign = -sign

        return tuple(minors)

    def berkowitz_charpoly(self, x):
        """Computes characteristic polynomial minors using Berkowitz method."""
        return Poly(map(simplify, self.berkowitz()[-1]), x)

    charpoly = berkowitz_charpoly

    def berkowitz_eigenvals(self, **flags):
        """Computes eigenvalues of a Matrix using Berkowitz method. """
        return roots(self.berkowitz_charpoly(Dummy('x')), **flags)

    eigenvals = berkowitz_eigenvals

    def eigenvects(self, **flags):
        """Return list of triples (eigenval, multiplicity, basis)."""

        if 'multiple' in flags:
            del flags['multiple']

        out, vlist = [], self.eigenvals(**flags)

        for r, k in vlist.iteritems():
            tmp = self - eye(self.rows)*r
            basis = tmp.nullspace()
            # whether tmp.is_symbolic() is True or False, it is possible that
            # the basis will come back as [] in which case simplification is
            # necessary.
            if not basis:
                # The nullspace routine failed, try it again with simplification
                basis = tmp.nullspace(simplified=True)
                if not basis:
                    raise NotImplementedError("Can't evaluate eigenvector for eigenvalue %s" % r)
            out.append((r, k, basis))
        return out

    def fill(self, value):
        """Fill the matrix with the scalar value."""
        self.mat = [value] * self.rows * self.cols

    def __getattr__(self, attr):
        if attr in ('diff','integrate','limit'):
            def doit(*args):
                item_doit = lambda item: getattr(item, attr)(*args)
                return self.applyfunc( item_doit )
            return doit
        else:
            raise AttributeError()

    def vec(self):
        """
        Return the Matrix converted into a one column matrix by stacking columns

        >>> from sympy import Matrix
        >>> m=Matrix([ [1,3], [2,4] ])
        >>> m
        [1, 3]
        [2, 4]
        >>> m.vec()
        [1]
        [2]
        [3]
        [4]

        """
        return Matrix(self.cols*self.rows, 1, self.transpose().mat)

    def vech(self, diagonal=True, check_symmetry=True):
        """
        Return the unique elements of a symmetric Matrix as a one column matrix
        by stacking the elements in the lower triangle.

        Arguments:
        diagonal -- include the diagonal cells of self or not
        check_symmetry -- checks symmetry of self but not completely reliably

        >>> from sympy import Matrix
        >>> m=Matrix([ [1,2], [2,3] ])
        >>> m
        [1, 2]
        [2, 3]
        >>> m.vech()
        [1]
        [2]
        [3]
        >>> m.vech(diagonal=False)
        [2]

        """
        c = self.cols
        if c != self.rows:
            raise TypeError("Matrix must be square")
        if check_symmetry:
            self.simplify()
            if self != self.transpose():
                raise ValueError("Matrix appears to be asymmetric; consider check_symmetry=False")
        count = 0
        if diagonal:
            v = zeros( (c * (c + 1) // 2, 1) )
            for j in xrange(c):
                for i in xrange(j,c):
                    v[count] = self[i,j]
                    count += 1
        else:
            v = zeros( (c * (c - 1) // 2, 1) )
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

        Example:

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
                    to_the_right = M[0:i, i:]
                    to_the_bottom = M[i:, 0:i]
                if any(to_the_right) or any(to_the_bottom):
                    i += 1
                    continue
                else:
                    sub_blocks.append(M[0:i, 0:i])
                    if M.shape == M[0:i, 0:i].shape:
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

        Example:

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
        [-1/2, 0, -1/2]
        [   0, 0, -1/2]
        [   1, 1,    1]
        >>> P.inv() * m * P
        [1, 0, 0]
        [0, 2, 0]
        [0, 0, 3]

        See also: .is_diagonalizable(), .is_diagonal()
        """
        if not self.is_square:
            raise NonSquareMatrixException()
        if not self.is_diagonalizable(reals_only, False):
            self._diagonalize_clear_subproducts()
            raise MatrixError("Matrix is not diagonalizable")
        else:
            if self._eigenvects == None:
                self._eigenvects = self.eigenvects()
            diagvals = []
            P = Matrix(self.rows, 0, [])
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

        Example:

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
        self._eigenvects = self.eigenvects()
        all_iscorrect = True
        for eigenval, multiplicity, vects in self._eigenvects:
            if len(vects) <> multiplicity:
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

        Note:

        Calculation of transformation P is not implemented yet

        Example:

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

        See also: jordan_cells()
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

        Note:

        Calculation of transformation P is not implemented yet

        Example:

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

        See also: jordan_form()
        """
        if not self.is_square:
            raise NonSquareMatrixException()
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
            "return a list which sum is equal to 'algebraical' and length is equal to 'geometrical'"
            n1 = algebraical / geometrical
            res = [n1] * geometrical
            res[len(res)-1] += algebraical % geometrical
            assert sum(res) == algebraical
            return res

def matrix_multiply(A, B):
    """
    Matrix product A*B.

    A and B must be of appropriate dimensions.  If A is a m x k matrix, and B
    is a k x n matrix, the product will be an m x n matrix.

    Example:

    >>> from sympy import Matrix
    >>> A = Matrix([[1, 2, 3], [4, 5, 6]])
    >>> B = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> A*B
    [30, 36, 42]
    [66, 81, 96]
    >>> B*A
    Traceback (most recent call last):
    ...
    ShapeError
    >>>

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
        raise ShapeError()
    blst = B.T.tolist()
    alst = A.tolist()
    return Matrix(A.shape[0], B.shape[1], lambda i, j:
                                        reduce(lambda k, l: k+l,
                                        map(lambda n, m: n*m,
                                        alst[i],
                                        blst[j])))

def matrix_multiply_elementwise(A, B):
    """Return the Hadamard product (elementwise product) of A and B

    >>> import sympy
    >>> A = sympy.Matrix([[0, 1, 2], [3, 4, 5]])
    >>> B = sympy.Matrix([[1, 10, 100], [100, 10, 1]])
    >>> print sympy.matrices.matrix_multiply_elementwise(A, B)
    [  0, 10, 200]
    [300, 40,   5]
    """
    if A.shape != B.shape:
        raise ShapeError()
    shape = A.shape
    return Matrix(shape[0], shape[1],
        lambda i, j: A[i,j] * B[i, j])

def matrix_add(A,B):
    """Return A+B"""
    if A.shape != B.shape:
        raise ShapeError()
    alst = A.tolist()
    blst = B.tolist()
    ret = [0]*A.shape[0]
    for i in xrange(A.shape[0]):
        ret[i] = map(lambda j,k: j+k, alst[i], blst[i])
    return Matrix(ret)

def zero(n):
    """Create square zero matrix n x n"""
    warnings.warn( 'Deprecated: use zeros() instead.' )
    return zeronm(n,n)

def zeronm(n,m):
    """Create zero matrix n x m"""
    warnings.warn( 'Deprecated: use zeros() instead.' )
    assert n>0
    assert m>0
    return Matrix(n,m,[S.Zero]*m*n)

def zeros(dims):
    """Create zero matrix of dimensions dims = (d1,d2)"""
    n, m = _dims_to_nm(dims)
    return Matrix(n, m, [S.Zero]*m*n)

def one(n):
    """Create square all-one matrix n x n"""
    warnings.warn( 'Deprecated: use ones() instead.' )
    return Matrix(n,n,[S.One]*n*n)

def ones(dims):
    """Create all-one matrix of dimensions dims = (d1,d2)"""
    n, m = _dims_to_nm( dims )
    return Matrix(n, m, [S.One]*m*n)

def eye(n):
    """Create square identity matrix n x n

    See also: diag()
    """
    n = int(n)
    out = zeros(n)
    for i in range(n):
        out[i, i] = S.One
    return out

def diag(*values):
    """Create diagonal matrix from a list as a diagonal values.

    Arguments might be matrices too, in case of it they are fitted in result matrix

    Example:

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

    See also: eye()
    """
    rows = 0
    cols = 0
    for m in values:
        if isinstance(m, Matrix):
            rows += m.rows
            cols += m.cols
        else:
            rows += 1
            cols += 1
    res = zeros((rows, cols))
    i_row = 0
    i_col = 0
    for m in values:
        if isinstance(m, Matrix):
            res[i_row:i_row + m.rows, i_col:i_col + m.cols] = m
            i_row += m.rows
            i_col += m.cols
        else:
            res[i_row, i_col] = m
            i_row += 1
            i_col += 1
    return res

def block_diag(matrices):
    """
    Warning: this function is deprecated. See .diag()

    Constructs a block diagonal matrix from a list of square matrices.

    See also: diag(), eye()
    """
    warnings.warn("block_diag() is deprecated, use diag() instead", DeprecationWarning)
    return diag(*matrices)

def jordan_cell(eigenval, n):
    """
    Create matrix of Jordan cell kind:

    Example:

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

def randMatrix(r,c,min=0,max=99,seed=[]):
    """Create random matrix r x c"""
    if seed == []:
        prng = random.Random()  # use system time
    else:
        prng = random.Random(seed)
    return Matrix(r,c,lambda i,j: prng.randint(min,max))

def hessian(f, varlist):
    """Compute Hessian matrix for a function f

    see: http://en.wikipedia.org/wiki/Hessian_matrix
    """
    # f is the expression representing a function f, return regular matrix
    if isinstance(varlist, (list, tuple)):
        m = len(varlist)
    elif isinstance(varlist, Matrix):
        m = varlist.cols
        assert varlist.rows == 1
    else:
        raise ValueError("Improper variable list in hessian function")
    assert m > 0
    try:
        f.diff(varlist[0])   # check differentiability
    except AttributeError:
        raise ValueError("Function %d is not differentiable" % i)
    out = zeros(m)
    for i in range(m):
        for j in range(i,m):
            out[i,j] = f.diff(varlist[i]).diff(varlist[j])
    for i in range(m):
        for j in range(i):
            out[i,j] = out[j,i]
    return out

def GramSchmidt(vlist, orthog=False):
    out = []
    m = len(vlist)
    for i in range(m):
        tmp = vlist[i]
        for j in range(i):
            tmp -= vlist[i].project(out[j])
        if tmp == Matrix([[0,0,0]]):
            raise ValueError("GramSchmidt: vector set not linearly independent")
        out.append(tmp)
    if orthog:
        for i in range(len(out)):
            out[i] = out[i].normalized()
    return out

def wronskian(functions, var, method='bareis'):
    """Compute Wronskian for [] of functions

                   | f1    f2     ...   fn  |
                   | f1'   f2'    ...   fn' |
                   |  .     .     .      .  |
    W(f1,...,fn) = |  .     .      .     .  |
                   |  .     .       .    .  |
                   |  n     n           n   |
                   | D(f1) D(f2)  ...  D(fn)|

    see: http://en.wikipedia.org/wiki/Wronskian
    """

    for index in xrange(0, len(functions)):
        functions[index] = sympify(functions[index])
    n = len(functions)
    if n == 0:
        return 1
    W = Matrix(n, n, lambda i,j: functions[i].diff(var, j) )
    return W.det(method)

def casoratian(seqs, n, zero=True):
    """Given linear difference operator L of order 'k' and homogeneous
       equation Ly = 0 we want to compute kernel of L, which is a set
       of 'k' sequences: a(n), b(n), ... z(n).

       Solutions of L are linearly independent iff their Casoratian,
       denoted as C(a, b, ..., z), do not vanish for n = 0.

       Casoratian is defined by k x k determinant:

                  +  a(n)     b(n)     . . . z(n)     +
                  |  a(n+1)   b(n+1)   . . . z(n+1)   |
                  |    .         .     .        .     |
                  |    .         .       .      .     |
                  |    .         .         .    .     |
                  +  a(n+k-1) b(n+k-1) . . . z(n+k-1) +

       It proves very useful in rsolve_hyper() where it is applied
       to a generating set of a recurrence to factor out linearly
       dependent solutions and return a basis.

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

    return Matrix(k, k, f).det()

# Add sympify converters
def _matrix_sympify(matrix):
    raise SympifyError('Matrix cannot be sympified')
converter[Matrix] = _matrix_sympify
del _matrix_sympify


class SMatrix(Matrix):
    """Sparse matrix"""

    def __init__(self, *args):
        if len(args) == 3 and callable(args[2]):
            op = args[2]
            assert isinstance(args[0], int) and isinstance(args[1], int)
            self.rows = args[0]
            self.cols = args[1]
            self.mat = {}
            for i in range(self.rows):
                for j in range(self.cols):
                    value = sympify(op(i,j))
                    if value != 0:
                        self.mat[(i,j)] = value
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], (list, tuple)):
            self.rows = args[0]
            self.cols = args[1]
            mat = args[2]
            self.mat = {}
            for i in range(self.rows):
                for j in range(self.cols):
                    value = sympify(mat[i*self.cols+j])
                    if value != 0:
                        self.mat[(i,j)] = value
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], dict):
            self.rows = args[0]
            self.cols = args[1]
            self.mat = {}
            # manual copy, copy.deepcopy() doesn't work
            for key in args[2].keys():
                self.mat[key] = args[2][key]
        else:
            if len(args) == 1:
                mat = args[0]
            else:
                mat = args
            if not isinstance(mat[0], (list, tuple)):
                mat = [ [element] for element in mat ]
            self.rows = len(mat)
            self.cols = len(mat[0])
            self.mat = {}
            for i in range(self.rows):
                assert len(mat[i]) == self.cols
                for j in range(self.cols):
                    value = sympify(mat[i][j])
                    if value != 0:
                        self.mat[(i,j)] = value

    def __getitem__(self, key):
        if isinstance(key, slice) or isinstance(key, int):
            lo, hi = self.slice2bounds(key, self.rows*self.cols)
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
        assert len(key) == 2
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
        assert (0 <= num < self.rows * self.cols) or \
               (0 <= -1*num < self.rows * self.cols)
        i, j = 0, num
        while j >= self.cols:
            j -= self.cols
            i += 1
        return i,j

    def __setitem__(self, key, value):
        # almost identical, need to test for 0
        assert len(key) == 2
        if isinstance(key[0], slice) or isinstance(key[1], slice):
            if isinstance(value, Matrix):
                self.copyin_matrix(key, value)
            if isinstance(value, (list, tuple)):
                self.copyin_list(key, value)
        else:
            i,j=self.key2ij(key)
            testval = sympify(value)
            if testval != 0:
                self.mat[(i,j)] = testval
            if (i,j) in self.mat:
                del self.mat[(i,j)]

    def row_del(self, k):
        newD = {}
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
        newD = {}
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

        >>> from sympy.matrices import SMatrix
        >>> a=SMatrix((1,2),(3,4))
        >>> a
        [1, 2]
        [3, 4]
        >>> a.RL
        [(0, 0, 1), (0, 1, 2), (1, 0, 3), (1, 1, 4)]
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
        >>> from sympy.matrices import SMatrix
        >>> a=SMatrix((1,2),(3,4))
        >>> a
        [1, 2]
        [3, 4]
        >>> a.CL
        [(0, 0, 1), (1, 0, 3), (0, 1, 2), (1, 1, 4)]
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
        Returns the transposed SMatrix of this SMatrix
        >>> from sympy.matrices import SMatrix
        >>> a = SMatrix((1,2),(3,4))
        >>> a
        [1, 2]
        [3, 4]
        >>> a.T
        [1, 3]
        [2, 4]
        """
        tran = SMatrix(self.cols,self.rows,{})
        for key,value in self.mat.iteritems():
            tran.mat[key[1],key[0]]=value
        return tran

    T = property(transpose,None,None,"Matrix transposition.")


    def __add__(self, other):
        if isinstance(other, SMatrix):
            return self.add(other)
        else:
            raise NotImplementedError("Only SMatrix + SMatrix supported")

    def __radd__(self, other):
        if isinstance(other, SMatrix):
            return self.add(other)
        else:
            raise NotImplementedError("Only SMatrix + SMatrix supported")

    def add(self, other):
        """
        Add two sparse matrices with dictionary representation.

        >>> from sympy.matrices.matrices import SMatrix
        >>> A = SMatrix(5, 5, lambda i, j : i * j + i)
        >>> A
        [0, 0,  0,  0,  0]
        [1, 2,  3,  4,  5]
        [2, 4,  6,  8, 10]
        [3, 6,  9, 12, 15]
        [4, 8, 12, 16, 20]
        >>> B = SMatrix(5, 5, lambda i, j : i + 2 * j)
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
        return SMatrix(self.rows, self.cols, c)



    # from here to end all functions are same as in matrices.py
    # with Matrix replaced with SMatrix
    def copyin_list(self, key, value):
        assert isinstance(value, (list, tuple))
        self.copyin_matrix(key, SMatrix(value))

    def multiply(self,b):
        """Returns self*b """

        def dotprod(a,b,i,j):
            assert a.cols == b.rows
            r=0
            for x in range(a.cols):
                r+=a[i,x]*b[x,j]
            return r

        r = SMatrix(self.rows, b.cols, lambda i,j: dotprod(self,b,i,j))
        if r.rows == 1 and r.cols ==1:
            return r[0,0]
        return r

    def submatrix(self, keys):
        assert isinstance(keys[0], slice) or isinstance(keys[1], slice)
        rlo, rhi = self.slice2bounds(keys[0], self.rows)
        clo, chi = self.slice2bounds(keys[1], self.cols)
        if not ( 0<=rlo<=rhi and 0<=clo<=chi ):
            raise IndexError("Slice indices out of range: a[%s]"%repr(keys))
        return SMatrix(rhi-rlo, chi-clo, lambda i,j: self[i+rlo, j+clo])

    def reshape(self, _rows, _cols):
        if self.rows*self.cols != _rows*_cols:
            print "Invalid reshape parameters %d %d" % (_rows, _cols)
        newD = {}
        for i in range(_rows):
            for j in range(_cols):
                m,n = self.rowdecomp(i*_cols + j)
                if (m,n) in self.mat:
                    newD[(i,j)] = self.mat[(m,n)]
        return SMatrix(_rows, _cols, newD)

    def cross(self, b):
        assert isinstance(b, (list, tuple, Matrix))
        if not (self.rows == 1 and self.cols == 3 or \
                self.rows == 3 and self.cols == 1 ) and \
                (b.rows == 1 and b.cols == 3 or \
                b.rows == 3 and b.cols == 1):
            raise ValueError("Dimensions incorrect for cross product")
        else:
            return SMatrix(1,3,((self[1]*b[2] - self[2]*b[1]),
                               (self[2]*b[0] - self[0]*b[2]),
                               (self[0]*b[1] - self[1]*b[0])))

    def zeronm(self,n,m):
        warnings.warn( 'Deprecated: use zeros() instead.' )
        return SMatrix(n,m,{})

    def zero(self, n):
        warnings.warn( 'Deprecated: use zeros() instead.' )
        return SMatrix(n,n,{})

    def zeros(self, dims):
        """Returns a dims = (d1,d2) matrix of zeros."""
        n, m = _dims_to_nm( dims )
        return SMatrix(n,m,{})

    def eye(self, n):
        tmp = SMatrix(n,n,lambda i,j:0)
        for i in range(tmp.rows):
            tmp[i,i] = 1
        return tmp


def list2numpy(l):
    """Converts python list of SymPy expressions to a NumPy array."""
    from numpy import empty
    a = empty(len(l), dtype=object)
    for i, s in enumerate(l):
        a[i] = s
    return a

def matrix2numpy(m):
    """Converts SymPy's matrix to a NumPy array."""
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

    The created symbols are named prefix_i1_i2_...  You should thus provide a
    non-empty prefix if you want your symbols to be unique for different output
    arrays, as Sympy symbols with identical names are the same object.

    Parameters
    ----------

    prefix : string
      A prefix prepended to the name of every symbol.

    shape : int or tuple
      Shape of the created array.  If an int, the array is one-dimensional; for
      more than one dimension the shape must be a tuple.

    Examples
    --------

    >> from sympy import symarray
    >> symarray('', 3)
    [_0 _1 _2]

    If you want multiple symarrays to contain distinct symbols, you *must*
    provide unique prefixes:

    >> a = symarray('', 3)
    >> b = symarray('', 3)
    >> a[0] is b[0]
    True
    >> a = symarray('a', 3)
    >> b = symarray('b', 3)
    >> a[0] is b[0]
    False

    Creating symarrays with a prefix:
    >> symarray('a', 3)
    [a_0 a_1 a_2]

    For more than one dimension, the shape must be given as a tuple:
    >> symarray('a', (2,3))
    [[a_0_0 a_0_1 a_0_2]
    [a_1_0 a_1_1 a_1_2]]
    >> symarray('a', (2,3,2))
    [[[a_0_0_0 a_0_0_1]
      [a_0_1_0 a_0_1_1]
      [a_0_2_0 a_0_2_1]]
    <BLANKLINE>
     [[a_1_0_0 a_1_0_1]
      [a_1_1_0 a_1_1_1]
      [a_1_2_0 a_1_2_1]]]

    """
    try:
        import numpy as np
    except ImportError:
        raise ImportError("symarray requires numpy to be installed")

    arr = np.empty(shape, dtype=object)
    for index in np.ndindex(shape):
        arr[index] = Symbol('%s_%s' % (prefix, '_'.join(map(str, index))))
    return arr

