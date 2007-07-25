
from sympy import *

import sympy.modules.polynomials
import random

class NonSquareMatrixException(Exception):
    pass

class Matrix(object):
    
    def __init__(self, *args):
        """
        Matrix can be constructed with values or a rule.
        
        >>> from sympy import *
        >>> Matrix( (1,2+I), (3,4) ) #doctest:+NORMALIZE_WHITESPACE
        1 2 + I
        3 4
        >>> Matrix(2, 2, lambda i,j: (i+1)*j ) #doctest:+NORMALIZE_WHITESPACE
        0 1
        0 2
        
        Note: in SymPy we count indices from 0. The rule however counts from 1.
        """
        if len(args) == 3 and callable(args[2]):
            operation = args[2]
            assert isinstance(args[0], int) and isinstance(args[1], int)
            self.lines = args[0]
            self.cols = args[1]
            self.mat = []
            for i in range(self.lines):
                for j in range(self.cols):
                    self.mat.append(Basic.sympify(operation(i, j)))
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], (list, tuple)):
            self.lines=args[0]
            self.cols=args[1]
            mat = args[2]
            self.mat=[]
            for j in range(self.lines):
                for i in range(self.cols):
                    self.mat.append(Basic.sympify(mat[j*self.cols+i]))
        else:
            if len(args) == 1:
                mat = args[0]
            else:
                mat = args
            if not isinstance(mat[0], (list, tuple)):
                # make each element a singleton
                mat = [ [element] for element in mat ]
            self.lines=len(mat)
            self.cols=len(mat[0])
            self.mat=[]
            for j in range(self.lines):
                assert len(mat[j])==self.cols
                for i in range(self.cols):
                    self.mat.append(Basic.sympify(mat[j][i]))
    
    def key2ij(self,key):
        """Converts key=(4,6) to 4,6 and ensures the key is correct."""
        if not (isinstance(key,(list, tuple)) and len(key) == 2):
            raise TypeError("wrong syntax: a[%s]. Use a[i,j] or a[(i,j)]"
                    %repr(key))
        i,j=key
        if not (i>=0 and i<self.lines and j>=0 and j < self.cols):
            print self.lines, " ", self.cols
            raise IndexError("Index out of range: a[%s]"%repr(key))
        return i,j
    
    def __getattr__(self,name):
        """
        >>> from sympy import *
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        1 2 + I
        3 4
        >>> m.T #doctest: +NORMALIZE_WHITESPACE
        1 3
        2 + I 4
        >>> m.H #doctest: +NORMALIZE_WHITESPACE
        1 3
        2 - I 4
        
        """
        out = self[:,:]
        if name == "T":
            #transposition
            out = out.reshape(self.cols, self.lines)
            out[:,:] = Matrix(self.cols,self.lines, lambda i,j: self[j,i])
            return out
        if name == "C":
            #conjugation
            out[:,:] = Matrix(self.lines,self.cols,
                    lambda i,j: self[i,j].conjugate())
            return out
        if name == "H":
            #hermite conjugation
            out.reshape(self.cols, self.lines)
            out[:,:] = self.T.C
            return out
        if name == "D":
            #dirac conjugation
            out.reshape(self.cols, self.lines)
            out = self.H * gamma(0)
            return out
        raise AttributeError("'%s' object has no attribute '%s'"%
                (self.__class__.__name__, name))
    
    def __getitem__(self,key):
        """
        >>> from sympy import *
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        1 2 + I
        3 4
        >>> m[1,0]
        3
        >>> m.H[1,0]
        2 - I
        
        """
        # row-wise decomposition of matrix
        if isinstance(key, slice) or isinstance(key, int):
            return self.mat[key]
        # proper 2-index access
        assert len(key) == 2
        if isinstance(key[0], int) and isinstance(key[1], int):
            i,j=self.key2ij(key)
            return self.mat[i*self.cols+j]
        elif isinstance(key[0], slice) or isinstance(key[1], slice):
            return self.submatrix(key)
        else:
            raise IndexError("Index out of range: a[%s]"%repr(key))
    
    def __setitem__(self,key,value):
        """
        >>> from sympy import *
        >>> m=Matrix(((1,2+I),(3,4)))
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        1 2 + I
        3 4
        >>> m[1,0]=9
        >>> m  #doctest: +NORMALIZE_WHITESPACE
        1 2 + I
        9 4
        
        """
        assert len(key) == 2
        if isinstance(key[0], slice) or isinstance(key[1], slice):
            if isinstance(value, Matrix):
                self.copyin_matrix(key, value)
            if isinstance(value, (list, tuple)):
                self.copyin_list(key, value)
        else:
            i,j=self.key2ij(key)
            self.mat[i*self.cols + j] = Basic.sympify(value)
    
    def copyin_matrix(self, key, value):
        rlo, rhi = self.slice2bounds(key[0], self.lines)
        clo, chi = self.slice2bounds(key[1], self.cols)
        assert value.lines == rhi - rlo and value.cols == chi - clo
        for i in range(value.lines):
            for j in range(value.cols):
                self[i+rlo, j+clo] = Basic.sympify(value[i,j])
    
    def copyin_list(self, key, value):
        assert isinstance(value, (list, tuple))
        self.copyin_matrix(key, Matrix(value))
    
    def hash(self):
        """Compute a hash every time, because the matrix elements
        could change."""
        return hash(self.__str__() )
    
    def __rmul__(self,a):
        assert not isinstance(a,Matrix)
        r=self.zeronm(self.lines,self.cols)
        for i in range(self.lines):
            for j in range(self.cols):
                r[i,j]=a*self[i,j]
        return r
    
    def expand(self):
        out = self[:,:]
        out[:,:] = Matrix(self.lines, self.cols, lambda i,j: self[i,j].expand())
        return out
    
    def combine(self):
        out = self[:,:]
        out[:,:] = Matrix(self.lines, self.cols, lambda i,j: self[i,j].combine())
        return out
    
    def subs(self,a,b):
        out = self[:,:]
        out[:,:] = Matrix(self.lines, self.cols, lambda i,j: self[i,j].subs())
        return out
    
    def __sub__(self,a):
        return self + (-a)
    
    def __mul__(self,a):
        if isinstance(a,Matrix):
            return self.multiply(a)
        out = self[:,:]
        out[:,:] = Matrix(self.lines, self.cols, lambda i,j: self[i,j]*a)
        return out
    
    def __pow__(self, num):
        if isinstance(num, int) or (isinstance(num, Rational) and num.isinteger()):
            if num < 0:
                a = self.inv() # A**-2 = (A**-1)**2
                num = -num
            else:
                a = self
            for i in range(1, num):
                a *= self
            return a
        raise NotImplementedError('Can only rise to the power of an integer for now')
    
    def __add__(self,a):
        return self.add(a)
    
    def __div__(self,a):
        return self * (Rational(1)/a)
    
    def multiply(self,b):
        """Returns self*b """
        
        def dotprod(a,b,i,j):
            assert a.cols == b.lines
            r=0
            for x in range(a.cols):
                r+=a[i,x]*b[x,j]
            return r
        
        r = Matrix(self.lines,b.cols, lambda i,j: dotprod(self,b,i,j))
        if r.lines == 1 and r.cols ==1:
            return r[0,0]
        return r
    
    def add(self,b):
        """Returns self+b """
        
        assert self.lines == b.lines
        assert self.cols == b.cols
        out = self[:,:]
        out[:,:] = Matrix(self.lines, self.cols, lambda i,j: self[i,j]+b[i,j])
        return out
    
    def __neg__(self):
        return -1*self
    
    def __eq__(self,a):
        if not isinstance(a, (Matrix, Basic)):
            a = Basic.sympify(a)
        return self.hash() == a.hash()
    
    def __ne__(self,a):
        if not isinstance(a, (Matrix, Basic)):
            a = Basic.sympify(a)
        return self.hash() != a.hash()
    
    def __repr__(self):
        return str(self)
    
    def inv(self, method="GE"):
        assert self.cols==self.lines
        # gaussian by default, also can be done by LU
        if method == "GE":
            return self.inverse_GE()
        elif method == "LU":
            return self.inverse_LU()
        else:
            raise "Inversion method unrecognized"
    
    def __str__(self):
        s="";
        for i in range(self.lines):
            for j in range(self.cols):
                s+="%s "%repr(self[i,j]);
            if i != self.lines - 1:
                s+="\n"
        return s
    
    def __mathml__(self):
        mml = ""
        for i in range(self.lines):
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
        for i in range(0, self.lines):
            self[i, j] = f(self[i, j], j)
    
    def row_swap(self, i, j):
        for k in range(0, self.cols):
            self[i, k], self[j, k] = self[j, k], self[i, k]
    
    def col_swap(self, i, j):
        for k in range(0, self.lines):
            self[k, i], self[k, j] = self[k, j], self[k, i]
    
    def row_del(self, i):
        self.mat = self.mat[:i*self.cols] + self.mat[(i+1)*self.cols:]
        self.lines -= 1
    
    def col_del(self, i):
        """
        >>> import sympy
        >>> M = sympy.modules.matrices.eye(3)
        >>> M.col_del(1)
        >>> M   #doctest: +NORMALIZE_WHITESPACE
        1 0
        0 0
        0 1
        """
        for j in range(self.lines-1, -1, -1):
            del self.mat[i+j*self.cols]
        self.cols -= 1
    
    def row_join(self, rhs):
        # concatenates two matrices along self's last and rhs's first col
        assert self.lines == rhs.lines
        newmat = self.zeronm(self.lines, self.cols + rhs.cols)
        newmat[:,:self.cols] = self[:,:]
        newmat[:,self.cols:] = rhs
        return newmat
    
    def col_join(self, bott):
        assert self.cols == bott.cols
        newmat = self.zeronm(self.lines+bott.lines, self.cols)
        newmat[:self.lines,:] = self[:,:]
        newmat[self.lines:,:] = bott
        return newmat
    
    def trace(self):
        assert self.cols == self.lines
        trace = 0
        for i in range(self.cols):
            trace += self[i,i]
        return self
    
    def submatrix(self, keys):
        """
        >>> from sympy import *
        >>> m = Matrix(4,4,lambda i,j: i+j)
        >>> m   #doctest: +NORMALIZE_WHITESPACE
        0 1 2 3
        1 2 3 4
        2 3 4 5
        3 4 5 6
        >>> m[0:1, 1]   #doctest: +NORMALIZE_WHITESPACE
        1
        >>> m[0:2, 0:1] #doctest: +NORMALIZE_WHITESPACE
        0
        1
        >>> m[2:4, 2:4] #doctest: +NORMALIZE_WHITESPACE
        4 5
        5 6
        """
        assert isinstance(keys[0], slice) or isinstance(keys[1], slice)
        rlo, rhi = self.slice2bounds(keys[0], self.lines)
        clo, chi = self.slice2bounds(keys[1], self.cols)
        if not ( 0<=rlo<=rhi and 0<=clo<=chi ):
            raise IndexError("Slice indices out of range: a[%s]"%repr(keys))
        return Matrix(rhi-rlo, chi-clo, lambda i,j: self[i+rlo, j+clo])
    
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
        0 1
        2 3
        >>> m.applyfunc(lambda i: 2*i)  #doctest: +NORMALIZE_WHITESPACE
        0 2
        4 6
        """
        assert callable(f)
        out = self[:,:]
        for i in range(self.lines):
            for j in range(self.cols):
                out[i,j] = f(self[i,j])
        return out
    
    def reshape(self, _rows, _cols):
        """
        >>> from sympy import *
        >>> m = Matrix(2,3,lambda i,j: 1)
        >>> m   #doctest: +NORMALIZE_WHITESPACE
        1 1 1
        1 1 1
        >>> m.reshape(1,6)  #doctest: +NORMALIZE_WHITESPACE
        1 1 1 1 1 1
        >>> m.reshape(3,2)  #doctest: +NORMALIZE_WHITESPACE
        1 1
        1 1
        1 1
        """
        if self.lines*self.cols != _rows*_cols:
            print "Invalid reshape parameters %d %d" % (_rows, _cols)
        return Matrix(_rows, _cols, lambda i,j: self.mat[i*_cols + j])
    
    def print_nonzero (self, symb="X"):
        """
        Shows location of non-zero entries for fast shape lookup
        >>> from sympy import *
        >>> m = Matrix(2,3,lambda i,j: i*3+j)
        >>> m           #doctest: +NORMALIZE_WHITESPACE
        0 1 2
        3 4 5
        >>> m.print_nonzero()   #doctest: +NORMALIZE_WHITESPACE
        [ XX]
        [XXX]
        >>> m = modules.matrices.eye(4)
        >>> m.print_nonzero("x")    #doctest: +NORMALIZE_WHITESPACE
        [x   ]
        [ x  ]
        [  x ]
        [   x]
        """
        s="";
        for i in range(self.lines):
            s+="["
            for j in range(self.cols):
                if self[i,j] == 0:
                    s+=" "
                else:
                    s+= symb+""
            s+="]\n"
        print s
    
    def LUsolve(self, rhs):
        assert rhs.lines == self.lines
        A, perm = self.LUdecomposition_Simple()
        n = self.lines
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
    
    def LUdecomposition(self):
        combined, p = self.LUdecomposition_Simple()
        L = self.zero(self.lines)
        U = self.zero(self.lines)
        for i in range(self.lines):
            for j in range(self.lines):
                if i > j:
                    L[i,j] = combined[i,j]
                else:
                    if i == j:
                        L[i,i] = 1
                    U[i,j] = combined[i,j]
        return L, U, p
    
    def LUdecomposition_Simple(self):
        # returns A compused of L,U (L's diag entries are 1) and
        # p which is the list of the row swaps (in order)
        assert self.lines == self.cols
        n = self.lines
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
                if pivot == -1 and A[i,j] != 0:
                    pivot = i
            if pivot < 0:
                raise "Error: non-invertible matrix passed to LUdecomposition_Simple()"
            if pivot != j: # row must be swapped
                A.row_swap(pivot,j)
                p.append([pivot,j])
            assert A[j,j] != 0
            scale = 1 / A[j,j]
            for i in range(j+1,n):
                A[i,j] = A[i,j] * scale
        return A, p
    
    def LUdecomposition_Block(self):
        raise NotImplementedError("Not yet implemented")
    
    def cofactorMatrix(self):
        out = self[:,:]
        out[:,:] = Matrix(self.lines, self.cols, lambda i,j: self.cofactor(i,j))
        return out
    
    def minorEntry(self, i, j):
        assert 0 <= i < self.lines and 0 <= j < self.cols
        return self.minorMatrix(i,j).det()
    
    def minorMatrix(self, i, j):
        assert 0 <= i < self.lines and 0 <= j < self.cols
        return self.delRowCol(i,j)
    
    def cofactor(self, i, j):
        if (i+j) % 2 == 0:
            return self.minorEntry(i,j)
        else:
            return -1 * self.minorEntry(i,j)
    
    def jacobian(self, varlist):
        # self is a vector of expression representing functions f_i(x_1, ..., x_n)
        # varlist is the set of x_i's in order
        assert self.lines == 1
        m = self.cols
        if isinstance(varlist, Matrix):
            assert varlist.lines == 1
            n = varlist.cols
        elif isinstance(varlist, (list, tuple)):
            n = len(varlist)
        assert n > 0 # need to diff by something
        J = self.zeronm(m,n) # maintain subclass type
        for i in range(m):
            if isinstance(self[i], (float, int)):
                continue    # constant function, jacobian row is zero
            try:
                tmp = self[i].diff(varlist[0])   # check differentiability
                J[i,0] = tmp
            except AttributeError:
                raise "Function %d is not differentiable" % i
            for j in range(1,n):
                J[i,j] = self[i].diff(varlist[j])
        return J
    
    def QRdecomposition(self):
        # TODO: still doesn't work for large expressions, there's a bug in an eval somewhere
        # return Q*R where Q is orthogonal and R is upper triangular
        # assume full-rank square, for now
        assert self.lines == self.cols
        n = self.lines
        Q, R = self.zero(n), self.zero(n)
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
    
    # Utility functions
    def simplify(self):
        for i in range(self.lines):
            for j in range(self.cols):
                try:
                    self[i,j] = self[i,j].simplify()
                except:
                    pass
    
    def expand(self):
        for i in range(self.lines):
            for j in range(self.cols):
                self[i,j] = self[i,j].expand()
    
    def evaluate(self):
        for i in range(self.lines):
            for j in range(self.cols):
                self[i,j] = self[i,j].eval()
    
    def cross(self, b):
        assert isinstance(b, (list, tuple, Matrix))
        if not (self.lines == 1 and self.cols == 3 or \
                self.lines == 3 and self.cols == 1 ) and \
                (b.lines == 1 and b.cols == 3 or \
                b.lines == 3 and b.cols == 1):
            raise "Dimensions incorrect for cross product"
        else:
            return Matrix(1,3,((self[1]*b[2] - self[2]*b[1]),
                               (self[2]*b[0] - self[0]*b[2]),
                               (self[0]*b[1] - self[1]*b[0])))
    
    def dot(self, b):
        assert isinstance(b, (list, tuple, Matrix))
        if isinstance(b, (list, tuple)):
            m = len(b)
        else:
            m = b.lines * b.cols
        assert self.cols*self.lines == m
        prod = 0
        for i in range(m):
            prod += self[i] * b[i]
        return prod
    
    def norm(self):
        assert self.lines == 1 or self.cols == 1
        out = Basic.sympify(0)
        for i in range(self.lines * self.cols):
            out += self[i]*self[i]
        return out**Basic.Half()
    
    def project(self, v):
        # project ONTO v
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
        #used only for cofactors, makes a copy
        M = self[:,:]
        M.row_del(i)
        M.col_del(j)
        return M
    
    def zeronm(self, n, m):
        # used so that certain functions above can use this
        # then only this func need be overloaded in subclasses
        return Matrix(n,m,lambda i,j:0)
    
    def zero(self, n):
        return Matrix(n,n,lambda i,j:0)
    
    def eye(self, n):
        tmp = self.zero(n)
        for i in range(tmp.lines):
            tmp[i,i] = 1
        return tmp
    
    @property
    def is_square(self):
        return self.lines == self.cols
    
    def is_upper(self):
        for i in range(self.cols):
            for j in range(self.lines):
                if i > j and self[i,j] != 0:
                    return False
        return True
    
    def is_lower(self):
        for i in range(self.cols):
            for j in range(self.lines):
                if i < j and self[i, j] != 0:
                    return False
        return True
    
    def is_symbolic(self):
        for i in range(self.cols):
            for j in range(self.lines):
                if not self[i,j].atoms(type=Symbol):
                    return True
        return False
    
    def clone(self):
        return Matrix(self.lines, self.cols, lambda i, j: self[i, j])
    
    def det(self):
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
        
        M, n = self[:,:], self.lines
        
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
                        return Rational(0)
                
                # proceed with Bareis' fraction-free (FF)
                # form of Gaussian elimination algorithm
                for i in range(k+1, n):
                    for j in range(k+1, n):
                        D = M[k, k]*M[i, j] - M[i, k]*M[k, j]
                        
                        if k > 0:
                            M[i, j] = D / M[k-1, k-1]
                        else:
                            M[i, j] = D
            
            det = sign * M[n-1, n-1]
        
        return simplify(det)
    
    def inverse_LU(self):
        return self.LUsolve(self.eye(self.lines))
    
    def inverse_GE(self):
        assert self.lines == self.cols
        assert self.det() != 0
        big = self.row_join(self.eye(self.lines))
        red = big.rref()
        return red[0][:,big.lines:]
    
    def rref(self):
        # take any matrix and return reduced row-ech form and indices of pivot vars
        # TODO: rewrite inverse_GE to use this
        pivots, r = 0, self[:,:]        # pivot: index of next row to contain a pivot
        pivotlist = []                  # indices of pivot variables (non-free)
        for i in range(r.cols):
            if pivots == r.lines:
                break
            if r[pivots,i] == 0:
                for k in range(pivots, r.lines):
                    if r[k,i] != 0:
                        break
                if k == r.lines - 1:
                    continue
                r.row_swap(pivots,k)
            scale = r[pivots,i]
            r.row(pivots, lambda x, _: x/scale)
            for j in range(0, r.lines):
                if j == pivots:
                    continue
                scale = r[j,i]
                r.row(j, lambda x, k: x - r[pivots,k]*scale)
            pivotlist.append(i)
            pivots += 1
        return r, pivotlist
    
    def nullspace(self):
        # Returns list of vectors (Matrix objects) that span nullspace of self
        assert self.cols >= self.lines
        reduced, pivots = self.rref()
        basis = []
        # create a set of vectors for the basis
        for i in range(self.cols - len(pivots)):
            basis.append(zeronm(1,self.cols))
        # contains the variable index to which the vector corresponds
        basiskey, cur = [-1]*len(basis), 0
        for i in range(self.cols):
            if i not in pivots:
                basiskey[cur] = i
                cur += 1
        for i in range(self.cols):
            if i not in pivots: # free var, just set vector's ith place to 1
                basis[basiskey.index(i)][0,i] = 1
            else:               # add negative of nonpivot entry to corr vector
                for j in range(i+1, self.cols):
                    line = pivots.index(i)
                    if reduced[line, j] != 0:
                        assert j not in pivots
                        basis[basiskey.index(j)][0,i] = -1 * reduced[line, j]
        return basis
    
    def charpoly(self, var):
        assert self.lines == self.cols
        if isinstance(var, Symbol):
            x = var
        else:
            raise "Input variable to charpoly not valid"
        copy = self[:,:]
        for i in range(self.lines):
            copy[i,i] -= x
        return copy.det()
    
    def eigenvals(self, var=None):
        """ Calls wrapper.py's roots(), doesn't support coeff type right now """
        # returns list of pairs (eigenval, multiplicty)
        if var == None:
            var = Symbol('x')
        p = self.charpoly(var)
        num, den = p.as_numer_denom()
        if den != 1:
            divop = sympy.modules.polynomials.div(num, den)
            assert divop[1] == 0
            p = divop[0]
        rl = sympy.modules.polynomials.roots(p, var)
        assert len(rl) == self.lines
        outlist = []
        def f(num):
            def g(n):
                if num == n:
                    return True
                else:
                    return False
            return g
        while rl != []:
            n = len(filter(f(rl[0]), rl))
            outlist.append([rl[0], n])
            for i in range(n):
                rl.remove(rl[0])
        return outlist
    
    def eigenvects(self):
        # return list of triples (eigenval, multiplicty, basis)
        out, vlist = [], self.eigenvals()
        for i in range(len(vlist)):
            tmp = self - eye(self.lines)*vlist[i][0]
            basis = tmp.nullspace()
            # check if basis is right size, don't do it if symbolic - too many solutions
            if not tmp.is_symbolic():
                assert len(basis) == vlist[i][1]
            vlist[i].append(basis)
            out.append(vlist[i])
        return out

def zero(n):
    return zeronm(n,n)

def zeronm(n,m):
    assert n>0
    assert m>0
    return Matrix(n,m, lambda i,j: 0)

def one(n):
    m = zero(n)
    for i in range(n):
        m[i,i]=1
    return m

def eye(n):
    assert n>0
    out = zeronm(n,n)
    for i in range(n):
        out[i,i]=1
    return out

def randMatrix(r,c,min=0,max=99,seed=[]):
    if seed == []:
        random.seed()
    else:
        random.seed(seed)       # use system time
    return Matrix(r,c,lambda i,j: random.randint(min,max))

def hessian(f, varlist):
    # f is the expression representing a function f, return regular matrix
    if isinstance(varlist, (list, tuple)):
        m = len(varlist)
    elif isinstance(varlist, Matrix):
        m = varlist.cols
        assert varlist.lines == 1
    else:
        raise "Improper variable list in hessian function"
    assert m > 0
    try:
        f.diff(varlist[0])   # check differentiability
    except AttributeError:
        raise "Function %d is not differentiable" % i
    out = zero(m)
    for i in range(m):
        for j in range(i,m):
            out[i,j] = f.diff(varlist[i]).diff(varlist[j])
    for i in range(m):
        for j in range(i):
            out[i,j] = out[j,i]
    return out

def GramSchmidt(vlist):
    out = []
    m = len(vlist)
    for i in range(m):
        tmp = vlist[i]
        for j in range(i):
            tmp -= vlist[i].project(out[j])
        if tmp == Matrix([[0,0,0]]):
            raise "GramSchmidt: vector set not linearly independent"
        out.append(tmp)
    return out

class SMatrix(Matrix):
    def __init__(self, *args):
        if len(args) == 3 and callable(args[2]):
            op = args[2]
            assert isinstance(args[0], int) and isinstance(args[1], int)
            self.lines = args[0]
            self.cols = args[1]
            self.mat = {}
            for i in range(self.lines):
                for j in range(self.cols):
                    value = Basic.sympify(op(i,j))
                    if value != 0:
                        self.mat[(i,j)] = value
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], (list, tuple)):
            self.lines = args[0]
            self.cols = args[1]
            mat = args[2]
            self.mat = {}
            for i in range(self.lines):
                for j in range(self.cols):
                    value = Basic.sympify(mat[i*self.cols+j])
                    if value != 0:
                        self.mat[(i,j)] = value
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], dict):
            self.lines = args[0]
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
            self.lines = len(mat)
            self.cols = len(mat[0])
            self.mat = {}
            for i in range(self.lines):
                assert len(mat[i]) == self.cols
                for j in range(self.cols):
                    value = Basic.sympify(mat[i][j])
                    if value != 0:
                        self.mat[(i,j)] = value
    
    def __getitem__(self, key):
        if isinstance(key, slice) or isinstance(key, int):
            lo, hi = self.slice2bounds(key, self.lines*self.cols)
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
            if self.mat.has_key((i,j)):
                return self.mat[(i,j)]
            else:
                return 0
        elif isinstance(key[0], slice) or isinstance(key[1], slice):
            return self.submatrix(key)
        else:
            raise IndexError("Index out of range: a[%s]"%repr(key))
    
    def rowdecomp(self, num):
        assert (0 <= num < self.lines * self.cols) or \
               (0 <= -1*num < self.lines * self.cols)
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
            testval = Basic.sympify(value)
            if testval != 0:
                self.mat[(i,j)] = testval
            elif self.mat.has_key((i,j)):
                del self.mat[(i,j)]
    
    def __str__(self):
        s = ""
        for i in range(self.lines):
            for j in range(self.cols):
                if self.mat.has_key((i,j)):
                    s += "%s " % repr(self[i,j])
                else:
                    s += "0 "
            s += "\n"
        return s
    
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
        self.lines -= 1
    
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
    
    # from here to end all functions are same as in matrices.py
    # with Matrix replaced with SMatrix
    def copyin_list(self, key, value):
        assert isinstance(value, (list, tuple))
        self.copyin_matrix(key, SMatrix(value))
    
    def multiply(self,b):
        """Returns self*b """
        
        def dotprod(a,b,i,j):
            assert a.cols == b.lines
            r=0
            for x in range(a.cols):
                r+=a[i,x]*b[x,j]
            return r
        
        r = SMatrix(self.lines, b.cols, lambda i,j: dotprod(self,b,i,j))
        if r.lines == 1 and r.cols ==1:
            return r[0,0]
        return r
    
    def submatrix(self, keys):
        assert isinstance(keys[0], slice) or isinstance(keys[1], slice)
        rlo, rhi = self.slice2bounds(keys[0], self.lines)
        clo, chi = self.slice2bounds(keys[1], self.cols)
        if not ( 0<=rlo<=rhi and 0<=clo<=chi ):
            raise IndexError("Slice indices out of range: a[%s]"%repr(keys))
        return SMatrix(rhi-rlo, chi-clo, lambda i,j: self[i+rlo, j+clo])
    
    def reshape(self, _rows, _cols):
        if self.lines*self.cols != _rows*_cols:
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
        if not (self.lines == 1 and self.cols == 3 or \
                self.lines == 3 and self.cols == 1 ) and \
                (b.lines == 1 and b.cols == 3 or \
                b.lines == 3 and b.cols == 1):
            raise "Dimensions incorrect for cross product"
        else:
            return SMatrix(1,3,((self[1]*b[2] - self[2]*b[1]),
                               (self[2]*b[0] - self[0]*b[2]),
                               (self[0]*b[1] - self[1]*b[0])))
    
    def zeronm(self,n,m):
        return SMatrix(n,m,{})
    
    def zero(self, n):
        return SMatrix(n,n,{})
    
    def eye(self, n):
        tmp = SMatrix(n,n,lambda i,j:0)
        for i in range(tmp.lines):
            tmp[i,i] = 1
        return tmp

