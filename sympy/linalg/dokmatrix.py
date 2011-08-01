from __future__ import division
from sympy import Basic, Symbol, Integer, C, S, Dummy, Rational, Add
from sympy.core.sympify import sympify, converter, SympifyError
from sympy.functions.elementary.miscellaneous import sqrt


from sympy.core.numbers import nan, oo

from sympy.polys import Poly, roots, cancel
from sympy.simplify import simplify as sympy_simplify
from sympy.printing import sstr

from sympy.linalg import DenseMatrix
from exceptions import *
from datamatrix import DataMatrix

import random



class DOKMatrix(DataMatrix):
    """Sparse matrix"""
    def __init__(self, *args, **kwargs):
        self.type = lambda i: i
        self._cached = kwargs.get('cached', True)
        if len(args) == 3 and callable(args[2]):
            op = args[2]
            if not isinstance(args[0], int) or not isinstance(args[1], int):
                raise TypeError("`args[0]` and `args[1]` must both be integers.")
            self.rows = args[0]
            self.cols = args[1]
            self.mat = {}
            for i in range(self.rows):
                for j in range(self.cols):
                    value = self.type(op(i,j))
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
                    value = self.type(mat[i*self.cols+j])
                    if value != 0:
                        self.mat[(i,j)] = value
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], dict):
            self.rows = args[0]
            self.cols = args[1]
            self.mat = {}
            # manual copy, copy.deepcopy() doesn't work
            for key in args[2].keys():
                val = args[2][key]
                if val != 0:
                    self.mat[key] = val
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
                if len(mat[i]) != self.cols:
                    raise ValueError("All arguments must have the same length.")
                for j in range(self.cols):
                    # value = sympify(mat[i][j])
                    value = self.type(mat[i][j])
                    if value != 0:
                        self.mat[(i,j)] = value

    def __getitem__(self, key):
        return self.mat.get(key, 0)
   
    def rowdecomp(self, num):
        nmax = len(self)
        if not (0 <= num < nmax) and not (0 <= -num < nmax):
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
            if isinstance(value, (DOKMatrix, Matrix)):
                self.copyin_matrix(key, value)
            if isinstance(value, (list, tuple)):
                self.copyin_list(key, value)
        else:
            i,j=self.key2ij(key)
            testval = self.type(value)
            if testval != 0:
                self.mat[(i, j)] = testval
            elif (i,j) in self.mat:
                del self.mat[(i, j)]

    def copyin_matrix(self, key, value):
        rlo, rhi = self.slice2bounds(key[0], self.rows)
        clo, chi = self.slice2bounds(key[1], self.cols)
        if value.rows != rhi - rlo or value.cols != chi - clo:
            raise ShapeError("The Matrix `value` doesn't have the same dimensions " +
                "as the in sub-Matrix given by `key`.")

        for i in range(value.rows):
            for j in range(value.cols):
                self[i+rlo, j+clo] = sympify(value[i,j])

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

    def _nonzero_elements_in_row(self, m):
        keys = sorted(self.mat.keys())
        start = 0
        end = len(keys)
        for i in xrange(len(keys)):
            if keys[i][0] == m - 1:
                start = i + 1
            if keys[i][0] == m + 1:
                end = i
                return keys[start:end]
        return keys[start:end]

    def _nonzero_elements_in_col(self, n):
        keys = [(i, j) for (j, i) in sorted((j,i) for (i, j) in self.mat.keys())]
        start = 0
        end = len(keys)
        for i in xrange(len(keys)):
            if keys[i][1] == n - 1:
                start = i + 1
            if keys[i][1] == n + 1:
                end = i
                return keys[start:end]
        return keys[start:end]

    def _lil_row_major(self):
        n = self.rows
        keys = sorted(self.mat.keys())
        lil = [[]] * n
        if not keys:
            return lil
        k = 0
        start = 0
        for i in xrange(len(keys)):
            if keys[i][0] > k:
                lil[k] = keys[start:i]
                start = i
                k = keys[i][0]
        lil[keys[-1][0]] = keys[start:]
        return lil

    _lil = _lil_row_major

    def _lil_col_major(self):
        n = self.cols
        keys = [(i, j) for (j, i) in sorted((j,i) for (i, j) in self.mat.keys())]
        lil = [[]] * n
        if not keys:
            return lil
        k = 0
        start = 0
        for i in xrange(len(keys)):
            if keys[i][1] > k:
                lil[k] = keys[start:i]
                start = i
                k = keys[i][0]
        lil[keys[-1][1]] = keys[start:]
        return lil
            
    def elem_ops_row(self, row1, row2, alpha):
        ''' row1 = row1 + alpha * row2 '''

        self._lil_row_major()
        rowlist1 = lil[row1]
        rowlist2 = lil[row2]
        for key in rowlist2:
            key1 = (row1, key[1])
            if key1 not in rowlist1:
                self.mat[key1] = alpha * self.mat[key]
            else:
                self.mat[key1] += alpha * self.mat[key]

    def elem_ops_col(self, col1, col2, alpha):
        "col1 = col1 + alpha * col2"

        self._lil_col_major()
        collist1 = lil[col1]
        collist2 = lil[col2]
        for key in collist2:
            key1 = (key[0], col1)
            if key1 not in collist1:
                self.mat[key1] = alpha * self.mat[key]
            else:
                self.mat[key1] += alpha * self.mat[key]

    def transpose(self):
        """
        Returns the transposed DOKMatrix of this DOKMatrix
        >>> from sympy.matrices import DOKMatrix
        >>> a = DOKMatrix((1,2),(3,4))
        >>> a
        [1, 2]
        [3, 4]
        >>> a.T
        [1, 3]
        [2, 4]
        """
        tran = DOKMatrix(self.cols,self.rows,{})
        for key,value in self.mat.iteritems():
            tran.mat[key[1],key[0]]=value
        return tran

    T = property(transpose,None,None,"Matrix transposition.")

    def __add__(self, other):
        if isinstance(other, DOKMatrix):
            return self.DOKMatrix_add(other) # DOKMatrix + DOKMatrix --> DOKMatrix 
        elif isinstance(other, Matrix):
            return other + self.toMatrix() # DOKMatrix + Matrix --> Matrix
        else:
            raise NotImplementedError("Addition/Subtraction of %s and %s not supported" % (type(self), type(other)))

    def __mul__(self, other):
        if isinstance(other, DOKMatrix):
            return DOK_matrix_multiply_row_major(self, other)
        elif isinstance(other, Matrix):
            return self * DOKMatrix.fromMatrix(other)
        else:
            return DOK_scalar_multiply(self, other)

    def __neg__(self):
        return -1 * self

    def __sub__(self,a):
        return self + (-a)
    
    def __rmul__(self, other):
        return self.__mul__(other)  

    def __eq__(self, other):
        if not isinstance(other, (Matrix, DOKMatrix)):
            return False
        if self.rows != other.rows or self.cols != other.cols:
            return False
        for i in xrange(self.rows):
            for j in xrange(self.cols):
                if self[i, j] != other[i, j]:
                    return False
        return True

    def DOKMatrix_add(self, other):
        if self.shape != other.shape:
            raise ShapeError()
        M = DOKMatrix(self.rows, self.cols, self.mat) # self[:,:] should be as fast as this
        for i in other.mat:
            if i in M.mat:
                M[i] += other[i]
            else:
                M[i] = other[i]
        return M

    # from here to end all functions are same as in matrices.py
    # with Matrix replaced with DOKMatrix
    def copyin_list(self, key, value):
        if not isinstance(value, (list, tuple)):
            raise TypeError("`value` must be of type list or tuple.")
        self.copyin_matrix(key, DOKMatrix(value))

    @classmethod
    def fromMatrix(cls, matrix):
        return DOKMatrix(matrix.rows, matrix.cols, lambda i, j: matrix[i, j])

    def LUdecom(self):
        cols = self._lil_col_major()
        A = self[:,:]
        for k in xrange(self.rows):
            for i, _ in cols[k]:
                if i > k:
                    A[i, k] = A[i, k] / A[k, k]
        pass
        return A

    def factor_lower(self):
        L = [None] * self.rows
        for i in xrange(self.rows):
            l = DOKMatrix.eye(self.rows)
            l[:, i] = self[:, i]
            L[i] = l
        return L 
                

    def submatrix(self, keys):
        if not isinstance(keys[0], slice) and not isinstance(keys[1], slice):
            raise TypeError("Both elements of `keys` must be slice objects.")
        rlo, rhi = self.slice2bounds(keys[0], self.rows)
        clo, chi = self.slice2bounds(keys[1], self.cols)
        if not ( 0<=rlo<=rhi and 0<=clo<=chi ):
            raise IndexError("Slice indices out of range: a[%s]"%repr(keys))
        return DOKMatrix(rhi-rlo, chi-clo, lambda i,j: self[i+rlo, j+clo])

    def reshape(self, _rows, _cols):
        if len(self) != _rows*_cols:
            print "Invalid reshape parameters %d %d" % (_rows, _cols)
        newD = {}
        for i in range(_rows):
            for j in range(_cols):
                m,n = self.rowdecomp(i*_cols + j)
                if (m,n) in self.mat:
                    newD[(i,j)] = self.mat[(m,n)]
        return DOKMatrix(_rows, _cols, newD)

    def cross(self, b):
        if not isinstance(b, (list, tuple, Matrix)):
            raise TypeError("`b` must be of type list, tuple, or Matrix, not %s." %
                type(b))
        if not (self.rows == 1 and self.cols == 3 or \
                self.rows == 3 and self.cols == 1 ) and \
                (b.rows == 1 and b.cols == 3 or \
                b.rows == 3 and b.cols == 1):
            raise ShapeError("Dimensions incorrect for cross product")
        else:
            return DOKMatrix(1,3,((self[1]*b[2] - self[2]*b[1]),
                               (self[2]*b[0] - self[0]*b[2]),
                               (self[0]*b[1] - self[1]*b[0])))


    def zeros(self, dims):
        """Returns a dims = (d1,d2) matrix of zeros."""
        n, m = _dims_to_nm( dims )
        return DOKMatrix(n,m,{})

    def __str__(self):
        return sstr(self.toMatrix())

    def __repr__(self):
        return sstr(self.toMatrix())

    def key2ij(self,key):
        """Converts key=(4,6) to 4,6 and ensures the key is correct."""
        if not (isinstance(key,(list, tuple)) and len(key) == 2):
            raise TypeError("wrong syntax: a[%s]. Use a[i,j] or a[(i,j)]"
                    %repr(key))
        i,j=key
        if not (i>=0 and i<self.rows and j>=0 and j < self.cols):
            raise IndexError("Index out of range: a[%s]"%repr(key))
        return i,j

    @property
    def shape(self):
        return (self.rows, self.cols)

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

    def __len__(self):
        return self.rows*self.cols

    def sparsity(self):
        return len(self.mat)/len(self)

    def _lower_columnar_nonzero_structure(self):
        n = self.cols
        NZlist = [[]] * n
        keys = [(i, j) for (j, i) in sorted((j,i) for (i, j) in self.mat.keys())]
        k = 0
        startset = False
        for i in xrange(len(keys)):
            if startset and keys[i][1] > k:
                NZlist[k] = keys[start:i]
                startset = False
            k = keys[i][1]
            if not startset and keys[i][1] == k and keys[i][1] <= keys[i][0]:
                start = i
                startset = True
        if (n-1, n-1) == keys[-1]:
            NZlist[-1].append((n-1,n-1))
        return NZlist

    def _lower_row_nonzero_structure(self):
        n = self.rows
        NZlist = [[]] * n
        keys = sorted(self.mat.keys())
        k = 0
        startset = False
        for i in xrange(len(keys)):
            if startset and (keys[i][0] > k or keys[i][1] > keys[i][0]):
                NZlist[k] = keys[start:i]
                startset = False
            k = keys[i][0]
            if not startset and keys[i][1] <= keys[i][0] and keys[i][0] == k and not NZlist[k]:
                start = i
                startset = True
        NZlist[keys[start][0]] = keys[start:]
        return NZlist

    def cholesky_structure_wiki(self):
        A = self._lower_columnar_nonzero_structure()
        n = self.cols
        L = [[]] * n
        parent = [0] * n
        for i in range(n):
            L[i] = A[i]
            for j in range(n):
                if parent[j] == i:
                    L[i].extend(L[j])
                    if (j,i) in L[i]:
                        L[i].remove((j,i))
                    L[i].sort()
            parent[i] = min(elem[0] for elem in L[i])
        return L

            
    def cholesky_structure(self): # :(
        def add(E, ij):
            if ij not in E[ij[1]]:
                E[ij[1]].append(ij)
                E[ij[1]].sort()
        def setmin(a):
            try:
                return min(a)
            except:
                return 0

        n = self.cols
        parent = [0] * n
        E = self._lower_columnar_nonzero_structure()
        for k in range(n):
            parent[k] = setmin(i[0] for i in E[k])
            for i, j in E[k]:
                add(E, (i, parent[i]))
        return E

    def test_cholesky_structure(self):
        A = self.cholesky_structure()
        L = []
        for i in A:
            L.extend(i)
        C = DOKMatrix(self.rows, self.cols, lambda i, j: 1 if (i,j) in L else 0)
        return self, C, C - self

    def test_cholesky_structure_wiki(self):
        A = self.cholesky_structure_wiki()
        L = []
        for i in A:
            L.extend(i)
        C = DOKMatrix(self.rows, self.cols, lambda i, j: 1 if (i,j) in L else 0)
        return self, C, C - self

    def liu(self):
        R = self._lower_row_nonzero_structure()
        parent = [None] * self.rows
        for j in xrange(self.rows):
            parent[j] = None
            for i in xrange(len(R[j])):
                r = R[j][i][1]
                while parent[r] != None:
                    r = parent[r]
                if r != j:
                    parent[r] = j
        return parent

    def liupc(self):
        R = self._lower_row_nonzero_structure()
        parent = [None] * self.rows
        virtual = [None] * self.rows
        for j in xrange(self.rows):
            for i, (_,r) in enumerate(R[j][:-1]):
                while virtual[r] and virtual[r] < j:
                    t = virtual[r]
                    virtual[r] = j
                    r = t
                if not virtual[r]:
                    virtual[r] = j
                    parent[r] = j
        return R, parent

    def row_structure_symbolic_cholesky(self):
        from copy import deepcopy
        R, parent = self.liupc()
        Lrow = deepcopy(R)
        for k in xrange(self.rows):
            for _, j in R[k]:
                while j != None and j != k:
                    if (k, j) not in Lrow[k]:
                        Lrow[k].append((k, j))
                        Lrow[k].sort()
                    j = parent[j]
        return Lrow

    def elementary_symbolic_cholesky(self):
        C = self._lower_columnar_nonzero_structure()
        for col in C:
            for iup, up in enumerate(col[1:]):
                for down in col[iup+1:]:
                    if (down[0], up[0]) not in C[up[0]]:
                        C[up[0]].append((down[0], up[0]))
                        C[up[0]].sort()
        return C

    def fast_symbolic_cholesky(self): 
        """ implement algorithm 1.3 from 10.1.1.163.7506 """
        C = self._lower_columnar_nonzero_structure()
        p=[]
        for k, col in enumerate(C):
            if len(C[k]) > 1:
                parent = C[k][1][0]
                p.append(parent)
            else:
                continue
            for i, _ in col[1:]: # add (parent[j], j) to C
                if (i, parent) not in C[parent]:
                    C[parent].append((i, parent))
                    C[parent].sort()
        return C

    def _cholesky(self):

        cached_result = self._get_cache('CH')
        if cached_result and self._cached:
            return cached_result

        CHstruc = self.fast_symbolic_cholesky()
        C = DOKMatrix(self.rows, self.cols, {})
        for col in CHstruc:
            for i, j in col:
                if i != j:
                    C[i, j] = (1 / C[j, j]) * (self[i, j] - sum(C[i, k] * C[j, k] 
                        for k in xrange(j)))
                elif i==j:
                    C[j, j] = (self[j, j] - sum(C[j, k] ** 2
                     for k in xrange(j))) ** (0.5)

        self._set_cache('CH', C)
        return C

    def _cholesky_sparse(self):

        cached_result = self._get_cache('CH')
        if cached_result and self._cached:
            return cached_result

        Crowstruc = self.row_structure_symbolic_cholesky()
        C = DOKMatrix.zeros(self.rows)
        for row in Crowstruc:
            for i, j in row:
                if i != j:
                    C[i, j] = self[i, j]
                    summ = 0
                    for p1 in Crowstruc[i]:
                        if p1[1] < j:
                            for p2 in Crowstruc[j]:
                                if p2[1] < j:
                                    if p1[1] == p2[1]:
                                        summ += C[p1] * C[p2]
                                else:
                                    break
                            else:
                                break
                    C[i, j] -= summ
                    C[i, j] /= C[j, j]
                else:
                    C[j, j] = self[j, j]
                    summ = 0
                    for _, k in Crowstruc[j]:
                        if k < j:
                            summ += C[j, k] ** 2
                        else:
                            break
                    C[j, j] -= summ
                    C[j, j] = C[j, j] ** 0.5

        self._set_cache('CH', C)
        return C
                    
            
    def _LDL(self):

        cached_result = self._get_cache('LDL')
        if cached_result and self._cached:
            return cached_result

        CHstruc = self.fast_symbolic_cholesky()
        L = DOKMatrix.eye(self.rows)
        D = DOKMatrix(self.rows, self.cols, {})
        
        for col in CHstruc:
            for i, j in col:
                if i != j:
                    L[i, j] = (self[i, j] - sum(L[i, k] * L[j, k] * D[k, k] 
                        for k in range(j))) / D[j, j]
                elif i == j:
                    D[i, i] = self[i, i] - sum(L[i, k]**2 * D[k, k] 
                        for k in range(i))

        self._set_cache('LDL', (L, D))
        return L, D

    
    def _LDL_sparse(self):

        cached_result = self._get_cache('LDL')
        if cached_result and self._cached:
            return cached_result

        Lrowstruc = self.row_structure_symbolic_cholesky()
        L = DOKMatrix.eye(self.rows)
        D = DOKMatrix(self.rows, self.cols, {})
        
        for row in Lrowstruc:
            for i, j in row:
                if i != j:
                    L[i, j] = self[i, j]
                    summ = 0
                    for p1 in Lrowstruc[i]:
                        if p1[1] < j:
                            for p2 in Lrowstruc[j]:
                                if p2[1] < j:
                                    if p1[1] == p2[1]:
                                        summ += L[p1] * L[p2] * D[p1[1], p1[1]]
                                else:
                                    break
                        else:
                            break
                    L[i, j] -= summ
                    L[i, j] /= D[j, j]
                elif i == j:
                    D[i, i] = self[i, i]
                    summ = 0
                    for _, k in Lrowstruc[i]:
                        if k < i:
                            summ += L[i, k]**2 * D[k, k]
                        else:
                            break
                    D[i, i] -= summ

        self._set_cache('LDL', (L, D)) 
        return L, D


    def _lower_triangular_solve(self, rhs):
        rows = self._lil_row_major()
        X = DOKMatrix(rhs.rows, 1, rhs.mat)
        for i in xrange(self.rows):
            for key in rows[i]:
                if key[1] < i:
                    X[i, 0] -= self[key] * X[key[1], 0]
                else:
                    break
            X[i, 0] /= self[i, i]
        return X

    def _upper_triangular_solve(self, rhs):
        """ Helper function of function upper_triangular_solve, without the error checks, to be used privately. """
        rows = self._lil_row_major()
        X = DOKMatrix(rhs.rows, 1, rhs.mat)
        for i in reversed(xrange(self.rows)):
            for key in reversed(rows[i]):
                if key[1] > i:
                    X[i, 0] -= self[key] * X[key[1], 0]
                else:
                    break
            X[i, 0] /= self[i, i]
        return X

    def _diagonal_solve(self, rhs):
        return DOKMatrix(self.rows, 1, lambda i, j: rhs[i, 0] / self[i, i])
    
    def _cholesky_solve(self, rhs):
        L = self._cholesky()
        Y = L._lower_triangular_solve(rhs)
        X = L.T._upper_triangular_solve(Y)
        return X

    def LDL_solve(self, rhs):
        if not self.is_square():
            raise Exception("Make a square matrix exception")
        if not self.rows == rhs.rows:
            raise Exception
        if not self.is_symmetric():
            rhs = self.T * rhs
            self = self.T * self
        L, D = self._LDL_sparse()
        z = L._lower_triangular_solve(rhs)
        y = D._diagonal_solve(z)
        x = L.T._upper_triangular_solve(y)
        if x.has(nan) or x.has(oo): # import this
            raise Exception
        return x

    def _LDLsolve(self, rhs):
        L, D = self._LDL_sparse()
        z = L._lower_triangular_solve(rhs)
        y = D._diagonal_solve(z)
        x = L.T._upper_triangular_solve(y)
        return x

    def is_symmetric(self):
        return all(self[i, j] == self[j, i] for i, j in self.mat.keys())

    def is_square(self):
        return self.rows == self.cols

    def has(self, expr):
        any(self[i, j].has(expr) for i, j in self.mat.keys())

    def det(self, method="LDL"):
        if method == "CH":
            L = (self.T * self)._cholesky_sparse()
            det = 1
            for i in xrange(L.rows):
                det *= L[i,i]
            return det
        elif method =="LDL":
            _, D = (self.T * self)._LDL_sparse()
            det = 1
            for i in xrange(D.rows):
                det *= D[i, i]
            return sqrt(det) 
        else:
            raise Exception('No such method')

    def conjgrad(A, b, x = None):
        if x == None:
            x = DOKMatrix.zeros((A.rows, 1))
        r = b - A * x
        p = r[:,:]
        rsold = r.T * r
        i=0
        # for i in xrange(A.rows):
        while(True):
            i += 1
            Ap = A * p
            alpha = rsold / (p.T * Ap)
            x = x + alpha * p
            r = r - alpha * Ap
            rsnew = r.T * r
            print rsnew
            if rsnew < 10 ** -15:
                break
            p = r + (rsnew/rsold ) * p
            rsold = rsnew
        print i
        return x 

    def conjgrad_pre(A, b, T = None, x = None, maxiter = None, 
        epsilon = 10**-15):
        if not maxiter:
            maxiter = 10 * A.rows
        if not T:
            T = DOKMatrix(A.rows, A.cols, lambda i, j:
                1/A[i, j] if i == j else 0)
        if not x:
            x = DOKMatrix.zeros((A.rows, 1))
        i = 0
        r = b - A * x
        d = T * r
        delnew = r.T * d
        del0 = delnew
        while(i < maxiter and delnew > (epsilon**2) * del0):
            q = A * d
            alpha = delnew/(d.T * q)
            x += alpha * d
            if not i % 50:
                r = b - A * x
            else:
                r = r - alpha * q
            s = T * r
            delold = delnew
            delnew = r.T * s
            beta = delnew / delold
            d = s + beta * d
            i = i + 1
        return x  
        
    
    def test_sparse_dense(self):
        A = self.toMatrix().cholesky().applyfunc(lambda i: 1 if i!=0 else 0)
        B = self.test_elementary_symbolic_cholesky()[1].toMatrix()
        return A, B, A - B

    def applyfunc(self, op):
        for i in self.mat:
            self.mat[i] = op(self.mat[i])

    @classmethod
    def eye(cls, n):
        A = cls(n, n, {})
        for i in xrange(n):
            A[i, i] = 1
        return A

    @classmethod
    def zeros(cls, dims):
        if isinstance(dims, tuple):
            return cls(dims[0], dims[1], {})
        else:
            return cls(dims, dims, {})

    @classmethod
    def _from_dict(cls, rows, cols, dok):
        return cls(rows, cols, dok)

    def scalar_multiply(self, scalar):
        mat = DOKMatrix(self.rows, self.cols, {})
        for key, value in self.mat.iteritems():
            mat[key] = scalar * value
        return mat

    def to_dokmatrix(self):
        return self

    def to_lilmatrix(self):
        from sympy.matrices.lilmatrix import LILMatrix
        return LILMatrix._from_dict(self.rows, self.cols, self.mat)

    def to_densematrix(self):
        return Matrix(self.rows, self.cols, lambda i, j: self[i, j])

    def copy(self):
        return DOKMatrix(self.rows, self.cols, self.mat)

    def inverse_solver(self, solver=None):
        from matrixutils import vecs2matrix
        if not solver or solver == 'CH':
            solver = self._cholesky_solve
        elif solver == 'LDL':
            solver = self.LDL_solve
        else:
            raise Exception('solver method not recognized')
        I = Matrix.eye(self.rows)
        return vecs2matrix([solver(I[:, i])
            for i in xrange(self.cols)], repr='dok')

    def _lil_row_lower(A):
        return A._lil_lower(0)

    def _lil_col_lower(A):
        return A._lil_lower(1)

    def _lil_lower(A, j):
        lower = [k for k in A.mat.iterkeys() if k[0] >= k[1]]
        return _lil(lower, A.cols, j)

    def _lil_row(A):
        return _lil(A.mat.keys(), A.cols, 0)

    def _lil_col(A):
        return _lil(A.mat.keys(), A.cols, 1) 


def list2lil(keys, n, j):
    lil = [[]] * n

def _lil(keys, n, j):
    keys.sort()
    zeros = set(xrange(n)) - set(k[j] for k in keys)
    lil = [ [] if i in zeros 
        else [k for k in keys if k[j] == i ]
        for i in xrange(n) ]
    return lil

lil = _lil

def DOK_matrix_multiply(self, other):
    C = DOKMatrix(self.rows, other.cols, {})
    assert self.cols == other.rows
    a = self.mat.keys()
    b = other.mat.keys()
    for i in a:
        for j in b:
            if i[1] == j[0]: 
                if (i[0], j[1]) in C.mat:
                    C[(i[0], j[1])] += self[i] * other[j]
                else:
                    C[(i[0], j[1])] = self[i] * other[j]
    if C.shape == (1, 1):
        return C[0, 0]
    return C

def DOK_matrix_multiply_row_major(A, B):
    rows1 = A._lil_row_major()
    rows2 = B._lil_row_major()
    Cdict = {}
    for k in xrange(A.rows):
        for _, j in rows1[k]:
            for _, n  in rows2[j]:
                temp = A[k, j] * B[j, n]
                if (k, n) in Cdict:
                    Cdict[k, n] += temp
                else:
                    Cdict[k, n] = temp
    C = DOKMatrix(A.rows, B.cols, Cdict)
    if C.shape == (1, 1):
        return C[0, 0]
    return C

def DOK_matrix_multiply_elem_major(A, B):
    rows = A._lil_row_major()
    cols = B._lil_col_major()
    C = DOKMatrix(A.rows, B.cols, {})
    for i in xrange(A.rows):
        for j in xrange(B.cols):
            temp = 0
            for _, q in rows[i]:
                for m, _ in cols[j]:
                    if q == m:
                        temp += A[i, q] * B[m, j]
            C[i, j] = temp
    return C
                        
                    
                

def DOK_matrix_multiply_col_major(A, B):
    cols = A._lil_col_major()
    rows = B._lil_row_major()
    

def DOK_scalar_multiply(matrix, scalar):
    C = DOKMatrix(matrix.rows, matrix.cols, {})
    for i in matrix.mat:
        C.mat[i] = scalar * matrix.mat[i]
    return C

def randDOKMatrix(i, j, min=1, max=1, sparsity=0.5):
    return DOKMatrix(i, j, lambda i, j: random.randint(min, max) if random.random() < sparsity else 0)

def randInvDOKMatrix(n, d, min=-5, max=10):
    A = DOKMatrix(n, n, lambda i, j: random.randint(min, max) if abs(i - j) <= d-1 else 0)
    return A

def randSymMatrix(n, d):
    A = randInvDOKMatrix(n, d)
    return A.T * A        

def test(n, d):
    A = randInvDOKMatrix(n, d)
    A = A.T * A
    L, D = A._LDL()
    Ls = DOKMatrix(L.rows, L.cols, lambda i, j: 1 if L[i,j] else 0)
    ld, dd = A.toMatrix().LDLdecomposition()
    Lsd = DOKMatrix(ld.rows, ld.cols, lambda i, j: 1 if ld[i,j]!=0 else 0)
    return A, Ls, Lsd, Ls - Lsd 
 
def test2(n, d):
    A = randInvDOKMatrix(n, d)
    A = A.T * A
    L, D = A._LDL()
    Lstruc = A.elementary_symbolic_cholesky()
    L = []
    for i in Lstruc: 
        L.extend(i)
    S = Matrix(n,n,lambda i,j:1 if (i,j) in L else 0)
    return ('MLDL, SLDL, MC, SC, struc',A.toMatrix().LDLdecomposition()[0].nonzero(), 
    A._LDL()[0].toMatrix().nonzero(), 
    A.toMatrix().cholesky().nonzero(),
    A._cholesky().toMatrix().nonzero(),
    S.nonzero())
           
def vecs2matrix(vecs):
    """Join column vectors to a matrix."""
    m = len(vecs[0])
    n = len(vecs)
    A = DOKMatrix.zeros((m, n))
    for i in xrange(m):
        for j in xrange(n):
            A[i,j] = vecs[j][i,0]
    return A

def transformation_matrix(f, A, B):
    """Return the transformation matrix representing the linear transformation f.

    By definition this is a_{ij} in

        f(v_j) == \sum a_{ij} w_i

    where v_j are the base vectors in A and w_i the base vectors in B.

    f is a linear mapping K^n -> K^m where K is a field. It can be a function
    accepting a vector as argument or a matrix.
    The transformation matrix is calculated for transformations from the vector
    space represented by base A to the space represented by base B (each given
    as a Matrix.)
    """
    if isinstance(f, (Matrix, DOKMatrix)):
        C = f
        f = lambda x: C*x
    return vecs2matrix([B.LDL_solve(f(A[:,i])) for i in xrange(A.cols)])

def mat(n, d):
    A = randInvDOKMatrix(n, d)
    return A.T * A
