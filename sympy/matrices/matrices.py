import warnings
from sympy import Basic, Symbol, Integer
from sympy.core import sympify

from sympy.core.basic import S, C
from sympy.polys import Poly, roots
from sympy.simplify import simplify
from sympy.utilities import any

# from sympy.printing import StrPrinter /cyclic/

import random

class NonSquareMatrixException(Exception):
    pass

class ShapeError(ValueError):
    """Wrong matrix shape"""
    pass

class MatrixError(Exception):
    pass

def _dims_to_nm( dims ):
    """Converts dimensions tuple (or any object with length 1 or 2) or scalar
    in dims to matrix dimensions n and m."""

    try:
        l = len( dims )
    except TypeError:
        dims = (dims,)
        l = 1

    # This will work for nd-array too when they are added to sympy.
    try:
        for dim in dims:
            assert (dim > 0) and isinstance( dim, int )
    except AssertionError:
        raise ValueError("Matrix dimensions should positive integers!")

    if l == 2:
        n, m = dims
    elif l == 1:
        n, m = dims[0], dims[0]
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
        return StrPrinter.doprint(self)
    def __repr__(self):
        return StrPrinter.doprint(self)

class Matrix(object):

    # Added just for numpy compatibility
    # TODO: investigate about __array_priority__
    __array_priority__ = 10.0

    def __init__(self, *args):
        """
        Matrix can be constructed with values or a rule.

        >>> from sympy import *
        >>> Matrix( ((1,2+I), (3,4)) ) #doctest:+NORMALIZE_WHITESPACE
        [1, 2 + I]
        [3,     4]
        >>> Matrix(2, 2, lambda i,j: (i+1)*j ) #doctest:+NORMALIZE_WHITESPACE
        [0, 1]
        [0, 2]

        """
        if len(args) == 3 and callable(args[2]):
            operation = args[2]
            assert isinstance(args[0], int) and isinstance(args[1], int)
            self.rows = args[0]
            self.cols = args[1]
            self.mat = []
            for i in range(self.rows):
                for j in range(self.cols):
                    self.mat.append(sympify(operation(i, j)))
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], (list, tuple)):
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
            elif not isinstance(mat, (list, tuple)):
                raise TypeError("Matrix constructor doesn't accept %s as input" % str(type(mat)))
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
            # TODO: on 0.7.0 delete this and uncomment the last line
            mat = args
            if not isinstance(mat[0], (list, tuple)):
                # make each element a singleton
                mat = [ [element] for element in mat ]
            warnings.warn("Deprecated constructor, use brackets: Matrix(%s)" % str(mat))
            self.rows=len(mat)
            self.cols=len(mat[0])
            self.mat=[]
            for j in xrange(self.rows):
                assert len(mat[j])==self.cols
                for i in xrange(self.cols):
                    self.mat.append(sympify(mat[j][i]))
            #raise TypeError("Data type not understood")

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

        >>> from sympy import *
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

        >>> from sympy import *
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
        >>> from sympy import *
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
        >>> from sympy import *
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

        >>> from sympy import *
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
        if hasattr(a, "__array__"):
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
        if hasattr(a, "__array__"):
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
            while n:
                if n % 2:
                    a = a * self
                    n -= 1
                self = self * self
                n = n // 2
            return a
        raise NotImplementedError('Can only rise to the power of an integer for now')

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
        return StrPrinter.doprint(self)

    def __repr__(self):
        return StrPrinter.doprint(self)

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
            return block_diag(r)
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
        """Elementary row operation using functor"""
        for j in range(0, self.cols):
            self[i, j] = f(self[i, j], j)

    def col(self, j, f):
        """Elementary column operation using functor"""
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

        >>> from sympy import *
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

        >>> from sympy import *
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
        >>> from sympy import *
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
        >>> from sympy import *
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
        >>> from sympy import *
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

    def slice2bounds(self, key, defmax):
        """
            Takes slice or number and returns (min,max) for iteration
            Takes a default maxval to deal with the slice ':' which is (none, none)
        """
        if isinstance(key, slice):
            lo, hi = 0, defmax
            if key.start != None:
                if key.start >= 0:
                    lo = key.start
                else:
                    lo = defmax+key.start
            if key.stop != None:
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
        >>> from sympy import *
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
        >>> from sympy import *
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
        Shows location of non-zero entries for fast shape lookup
        >>> from sympy import *
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
        s="";
        for i in range(self.rows):
            s+="["
            for j in range(self.cols):
                if self[i,j] == 0:
                    s+=" "
                else:
                    s+= symb+""
            s+="]\n"
        print s

    def LUsolve(self, rhs, iszerofunc=_iszero):
        """
        Solve the linear system Ax = b.
        self is the coefficient matrix A and rhs is the right side b.
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
        Returns the decompositon LU and the row swaps p.
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
        Returns A compused of L,U (L's diag entries are 1) and
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
        Returns 4 matrices P, L, D, U such that PA = L D**-1 U.

        From the paper "fraction-free matrix factors..." by Zhou and Jeffrey
        """
        n, m = self.rows, self.cols
        U, L, P = self[:,:], eye(n), eye(n)
        DD = zeros(n) # store it smarter since it's just diagonal
        oldpivot = 1

        for k in range(n-1):
            if U[k,k] == 0:
                kpivot = k+1
                Notfound = True
                while kpivot < n and Notfound:
                    if U[kpivot, k] != 0:
                        Notfound = False
                    else:
                        kpivot = kpivot + 1
                if kpivot == n+1:
                    raise ValueError("Matrix is not full rank")
                else:
                    swap = U[k, k:]
                    U[k,k:] = U[kpivot,k:]
                    U[kpivot, k:] = swap
                    swap = P[k, k:]
                    P[k, k:] = P[kpivot, k:]
                    P[kpivot, k:] = swap
            assert U[k, k] != 0
            L[k,k] = U[k,k]
            DD[k,k] = oldpivot * U[k,k]
            assert DD[k,k] != 0
            Ukk = U[k,k]
            for i in range(k+1, n):
                L[i,k] = U[i,k]
                Uik = U[i,k]
                for j in range(k+1, m):
                    U[i,j] = (Ukk * U[i,j] - U[k,j]*Uik) / oldpivot
                U[i,k] = 0
            oldpivot = U[k,k]
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

        self ... a vector of expressions representing functions f_i(x_1, ...,
                        x_n).
        X ...... is the set of x_i's in order, it can be a list or a Matrix

        Both self and X can be a row or a column matrix in any order
        (jacobian() should always work).

        Examples::
        >>> from sympy import symbols, sin, cos
        >>> rho, phi = symbols("rho phi")
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
        Return Q*R where Q is orthogonal and R is upper triangular.

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

    # TODO: QRsolve

    # Utility functions
    def simplify(self):
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
        for i in range(self.cols):
            for j in range(self.rows):
                if i > j and self[i,j] != 0:
                    return False
        return True

    def is_lower(self):
        for i in range(self.cols):
            for j in range(self.rows):
                if i < j and self[i, j] != 0:
                    return False
        return True

    def is_symbolic(self):
        for i in range(self.cols):
            for j in range(self.rows):
                if self[i,j].atoms(Symbol):
                    return True
        return False

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
           minimal numer of fractions. It means that less term
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
                            M[i, j] = Poly.cancel(D)

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
        assert self.cols >= self.rows
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

           This method is particulary useful for computing determinant,
           principal minors and characteristic polynomial, when 'self'
           has complicated coefficients eg. polynomials. Semi-direct
           usage of this algorithm is also important in computing
           efficiently subresultant PRS.

           Assuming that M is a square matrix of dimension N x N and
           I is N x N identity matrix,  then the following following
           definition of characteristic polynomial is begin used:

                          charpoly(M) = det(t*I - M)

           As a consequence, all polynomials generated by Berkowitz
           algorithm are monic.

           >>> from sympy import *
           >>> x,y,z = symbols('xyz')

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

           [2] M. Keber, Division-Free computation of subresultants
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
        coeffs, monoms = self.berkowitz()[-1], range(self.rows+1)
        return Poly(list(zip(coeffs, reversed(monoms))), x)

    charpoly = berkowitz_charpoly

    def berkowitz_eigenvals(self, **flags):
        """Computes eigenvalues of a Matrix using Berkowitz method. """
        return roots(self.berkowitz_charpoly(Symbol('x', dummy=True)), **flags)

    eigenvals = berkowitz_eigenvals

    def eigenvects(self, **flags):
        """Return list of triples (eigenval, multiplicty, basis)."""

        if 'multiple' in flags:
            del flags['multiple']

        out, vlist = [], self.eigenvals(**flags)

        for r, k in vlist.iteritems():
            tmp = self - eye(self.rows)*r
            basis = tmp.nullspace()
            # check if basis is right size, don't do it if symbolic - too many solutions
            if not tmp.is_symbolic():
                assert len(basis) == k
            elif len(basis) != k:
                # The nullspace routine failed, try it again with simplification
                basis = tmp.nullspace(simplified=True)
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

    def vech(self, diagonal=True):
        """
        Return the unique elements of a symmetric Matrix as a one column matrix by stacking
        the elements in the lower triangle

        diagonal ...... include the diagonal cells of self or not

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
        if self != self.transpose():
            raise TypeError("Matrix must be symmetric")
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
        >>> x, y, z = symbols('x y z')
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
    n, m = _dims_to_nm( dims )
    return Matrix(n,m,[S.Zero]*m*n)

def one(n):
    """Create square all-one matrix n x n"""
    warnings.warn( 'Deprecated: use ones() instead.' )
    return Matrix(n,n,[S.One]*n*n)

def ones(dims):
    """Create all-one matrix of dimensions dims = (d1,d2)"""
    n, m = _dims_to_nm( dims )
    return Matrix(n,m,[S.One]*m*n)

def eye(n):
    """Create square identity matrix n x n"""
    assert n>0
    out = zeros(n)
    for i in range(n):
        out[i,i]=S.One
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
    """Compute wronskian for [] of functions

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

       Solutions of L are lineary independent iff their Casoratian,
       denoted as C(a, b, ..., z), do not vanish for n = 0.

       Casoratian is defined by k x k determinant:

                  +  a(n)     b(n)     . . . z(n)     +
                  |  a(n+1)   b(n+1)   . . . z(n+1)   |
                  |    .         .     .        .     |
                  |    .         .       .      .     |
                  |    .         .         .    .     |
                  +  a(n+k-1) b(n+k-1) . . . z(n+k-1) +

       It proves very useful in rsolve_hyper() where it is applied
       to a generating set of a recurrence to factor out lineary
       dependent solutions and return a basis.

       >>> from sympy import *
       >>> n = Symbol('n', integer=True)

       Exponential and factorial are lineary independent:

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

def block_diag(matrices):
    """
    Constructs a block diagonal matrix from a list of square matrices.

    Example:
    >>> from sympy import block_diag, symbols
    >>> x, y, z = symbols("x y z")
    >>> a = Matrix([[1, 2], [2, 3]])
    >>> b = Matrix([[3, x], [y, 3]])
    >>> block_diag([a, b, b])
    [1, 2, 0, 0, 0, 0]
    [2, 3, 0, 0, 0, 0]
    [0, 0, 3, x, 0, 0]
    [0, 0, y, 3, 0, 0]
    [0, 0, 0, 0, 3, x]
    [0, 0, 0, 0, y, 3]

    """
    rows = 0
    for m in matrices:
        assert m.rows == m.cols, "All matrices must be square."
        rows += m.rows
    A = zeros((rows, rows))
    i = 0
    for m in matrices:
        A[i+0:i+m.rows, i+0:i+m.cols] = m
        i += m.rows
    return A


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
                if self.mat.has_key((m,n)):
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
            elif self.mat.has_key((i,j)):
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
                if self.mat.has_key((m,n)):
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
