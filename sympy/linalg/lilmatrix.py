from __future__ import division
import random
from sympy.matrices.matrices import Matrix
from sympy.printing import sstr, pretty
from sympy.simplify.simplify import simplify as sympy_simplify
from sympy.core.singleton import S
from matrixbase import MatrixBase
def _iszero(x):
    return x == 0

class LILMatrix(MatrixBase):
    def __init__(self, *args, **kwargs):
        if len(args) == 3 and callable(args[2]):
            "LILMatrix from lambda function"
            op = args[2]
            if not isinstance(args[0], int) or not isinstance(args[1], int):
                raise TypeError("`args[0]` and `args[1]` must both be integers.")
            self.rows = args[0]
            self.cols = args[1]
            self.mat = [[] for i in xrange(self.rows)]
            for i in xrange(self.rows):
                for j in xrange(self.cols):
                    value = op(i,j)
                    if value != 0:
                        self.mat[i].append((j, value))
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], (list, tuple)):
            "LILmatrix from a single list"
            self.rows = args[0]
            self.cols = args[1]
            mat = args[2]
            self.mat = [[] for i in xrange(self.rows)]
            for i in range(self.rows):
                for j in range(self.cols):
                    value = mat[i*self.cols+j]
                    if value != 0:
                        self.mat[i].append((j, value))
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], dict):
            "LILMatrix from a dictionary"
            self.rows = args[0]
            self.cols = args[1]
            self.mat = [[] for i in xrange(self.rows)]
            # manual copy, copy.deepcopy() doesn't work
            for key in args[2].keys():
                val = args[2][key]
                if val != 0:
                    self.mat[i].append((j, value))
        else:
            if len(args) == 1:
                mat = args[0]
            else:
                mat = args
            "LILMatri from a list of list"
            if not isinstance(mat[0], (list, tuple)):
                mat = [ [element] for element in mat ]
            self.rows = len(mat)
            self.cols = len(mat[0])
            self.mat = [[] for i in xrange(self.rows)]
            for i in range(self.rows):
                if len(mat[i]) != self.cols:
                    raise ValueError("All arguments must have the same length.")
                for j in range(self.cols):
                    value = mat[i][j]
                    if value != 0:
                        self.mat[i].append((j, value))

    def __str__(self):
        return sstr(self.to_densematrix())

    def __repr__(self):
        return sstr(self.to_densematrix())

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

    def __getitem__(self, key):
        i, j = key
        if type(i) is slice or type(j) is slice:
                return self.submatrix2(key)
        
        for ind, val in self.mat[i]:
            if ind >= j:
                if ind == j:
                    return val
                else:
                    return 0
        return 0

    def __setitem__(self, key, value):
        i, j = key
        for ind, (j2, val) in enumerate(self.mat[i]):
            if j2 >= j:
                if j2 == j:
                    if value == 0:
                        self.mat[i].pop(ind)
                    else:
                        self.mat[i][ind] = (j, value)
                    return
                else:
                    if value != 0:
                        self.mat[i].insert(ind, (j, value))
                    return
        if value != 0:
            self.mat[i].append((j, value))

    def submatrix2(self, keys):
        if not isinstance(keys[0], slice) and not isinstance(keys[1], slice):
            raise TypeError("At least one element of `keys` must be a slice object.")

        rlo, rhi = self.slice2bounds(keys[0], self.rows)
        clo, chi = self.slice2bounds(keys[1], self.cols)
        if not ( 0<=rlo<=rhi<=self.rows and 0<=clo<=chi<=self.cols ):
            raise IndexError("Slice indices out of range: a[%s]"%repr(keys))
        outRows, outCols = rhi-rlo, chi-clo
        outMat = []
        for i in xrange(rlo, rhi):
            startset = False
            start = end = len(self.mat[i])
            for ind, (j, val) in enumerate(self.mat[i]):
                if j >= clo and not startset:
                    start = ind
                    startset = True
                if j >= chi:
                    end = ind
                    break
            outMat.append(self.mat[i][start:end])
        for i in xrange(len(outMat)):
            for ind, (j, val) in enumerate(outMat[i]):
                 outMat[i][ind] = (j - clo, val)
        A = LILMatrix.zeros((outRows, outCols))
        A.mat = outMat
        return A

    def submatrix(self, keys):
        if not isinstance(keys[0], slice) and not isinstance(keys[1], slice):
            raise TypeError("At least one element of `keys` must be a slice object.")

        rlo, rhi = self.slice2bounds(keys[0], self.rows)
        clo, chi = self.slice2bounds(keys[1], self.cols)
        if not ( 0<=rlo<=rhi and 0<=clo<=chi ):
            raise IndexError("Slice indices out of range: a[%s]"%repr(keys))
        outRows, outCols = rhi-rlo, chi-clo
        outMat = []
        for i in xrange(rlo, rhi):
            startset = False
            start = 0
            end = len(self.mat[i])
            for ind, (j, val) in enumerate(self.mat[i]):
                if not startset and j >= clo:
                    start = ind
                    startset = True
                if j >= chi:
                    end = ind
                    break
            outMat.append(self.mat[i][start:end])
        A = LILMatrix.zeros((outRows, outCols))
        A.mat = outMat
        return A

    def copyin_matrix(self, key, value):
        rlo, rhi = self.slice2bounds(key[0], self.rows)
        clo, chi = self.slice2bounds(key[1], self.cols)
        if value.rows != rhi - rlo or value.cols != chi - clo:
            raise ShapeError("The Matrix `value` doesn't have the same dimensions " +
                "as the in sub-Matrix given by `key`.")
        raise NotImplemented

    def row_add(self, r1, r2, alpha):
        if r1 == r2:
            return
        row1 = self.mat[r1]
        row2 = self.mat[r2]
        self.mat[r1] = _row_add(row1, row2, alpha)

    def row_scale(self, r, alpha):
        for ind, (j, val) in enumerate(self.mat[r]):
            self.mat[r][ind] = (j, alpha * val)

    def row_add_bad(self, r1, r2, alpha):
        row2 = self.mat[r2]
        for elem in row2:
            self[r1,elem[0]] += alpha * elem[1]

    def row_functor(self, r, f):
        for i in xrange(len(self.mat[r])):
            self.mat[r][i] = self.mat[r][i][0],f(self.mat[r][i][1],self.mat[r][i][0])

    row = row_functor

    def row_swap(self, r1, r2):
        self.mat[r1], self.mat[r2] = self.mat[r2], self.mat[r1]

    @classmethod
    def _from_dict(cls, rows, cols, dok):
        mat = cls.zeros((rows, cols))
        for i in dok:
            mat[i] = dok[i]
        return mat    

    @classmethod
    def zeros(cls, shape):
        if isinstance(shape, tuple):
            return cls(shape[0], shape[1], lambda i,j:0)
        else:
            return cls(shape, shape, lambda i, j: 0)

    def to_densematrix(self):
        return Matrix(self.rows, self.cols, lambda i, j: self[i, j])

    def to_dokmatrix(self):
        from sympy import DOKMatrix
        mat = DOKMatrix(self.rows, self.cols, {})
        for i in xrange(self.rows):
            for j, value in self.mat[i]:
                mat[i, j] = value
        return mat

    def to_lilmatrix(self):
        return self

    def gauss_sparse(self):
        A = self.copy()
        for i in xrange(A.rows):
            if A[i, i] == 0:
                print 'bad pivot, exchanging', i
                for k in xrange(i + 1, A.rows):
                    print '\tconsidering', k, i
                    if A[k, i] != 0:
                        print '\t\tfound', k, i
                        print pretty(A)
                        print

                        A.row_swap(k, i)

                        print pretty(A)
                        break
                if A[i, i] == 0:
                    print 'bad bad pivot', i
                    print pretty(A)
                    raise Exception
            ind = 0
            l = len(A.mat[i])
            while ind < l:
                j, v = A.mat[i][ind]
                if j < i:
                    print 'zeroing out', i, j, A.mat[i]
                    A.row_add(i, j, - v / A[j, j])
                    print 'zeroed out', i, j, A.mat[i]
                    l = len(A.mat[i])
                else:
                    break
        return A

    def gauss_col(self):
        "Gaussian elimnation, currently tested only on square matrices"
        A = self.copy()
        row_swaps = []
        for j in xrange(A.rows):
            rlist = A.nz_col_lower(j)
            if A[j, j] == 0:
                if rlist:
                    A.row_swap(j, rlist[0])
                    row_swaps.append((j, rlist[0]))
                    rlist.pop(0)
                else:
                    continue
            for i in rlist:
                A.row_add(i, j, - A[i, j] / A[j, j])
        return A, row_swaps

    def _upper_triangular_solve(self, rhs):
        X = LILMatrix.zeros((rhs.rows, 1))
        for i in reversed(xrange(self.rows)):
            X[i, 0] = (rhs[i, 0] - sum(value * X[j, 0] for j, value in self.mat[i])) / self[i, i]
        return X

    def _lower_triangular_solve(self, rhs):
        X = LILMatrix.zeros((rhs.rows, 1))
        for i in xrange(self.rows):
            X[i, 0] = (rhs[i, 0] - sum(value * X[j, 0] for j, value in self.mat[i])) / self[i, i]
        return X

    def LUsolve(self, rhs):
        L, U, p = self.LU_sparse()
        b = rhs.copy()
        b.permute(p)
        Y = L._lower_triangular_solve(b)
        return U._upper_triangular_solve(Y)

    def solve_gauss(self, rhs):
        U, p = self.gauss_col()
        b = rhs.clone()
        b.permute(p)
        return U._upper_triangular_solve(b)

    def solve_rref(self, rhs):
        big = self.join_rows(rhs)
        rref = big.rref()
        return rref[:, self.cols]

    def solve(self, rhs, method="GE"):
        if method == "GE":
            return self.solve_gauss(rhs)
        elif method == "RREF":
            return self.solve_rref(rhs)
        else:
            raise ValueError('Unrecognised method')

    def __mul__(self, other):
        if isinstance(other, MatrixBase):
            prod = self.to_dokmatrix() * other.to_dokmatrix()
            return prod.to_lilmatrix() 
        else:
            # Scalar multiplication
            prod = self.clone()
            for i in xrange(self.rows):
                for ind, (j, value) in enumerate(self.mat[i]):
                    prod.mat[i][ind] = (j, other * value)
            return prod

    def __rmul__(self, other):
        ## assume other is scalar
        return self.__mul__(other)

    def __sub__(self, other):
        return self + (-1 * other)

    def det_gauss(self):
        ref, p = self.gauss_col()
        det = 1
        for i in xrange(ref.rows):
            det *= ref[i, i]
        if len(p) % 2 == 1:
            det *= -1
        return det

    def det_LU(self):
        L, U, p = self.LU_sparse()
        det = 1
        for i in xrange(U.rows):
            det *= U[i, i]
        if len(p) % 2 == 1:
            det *= -1
        return det

    def det(self, method='GE'):
        if method == "GE":
            return self.det_gauss()
        elif method == 'LU':
            return self.det_LU()
        else:
            raise Exception

    def rref2(self):
        pivot, r = 0, self.copy()
        pivotlist = []
        for i in range(r.cols):
            if pivot == r.rows:
                break
            if r[pivot,i] == 0:
                for k in xrange(pivot + 1, r.rows):
                    if r[k, i] != 0:
                        r.row_swap(pivot, k)
                        break
                if r[pivot, i] == 0:
                    continue
            r.row_scale(pivot, 1 / r[pivot, i])
            for j in r.nz_col(i):
                if j == pivot:
                    continue
                scale = r[j,i]
                r.row_add(j, pivot, - scale)
            pivotlist.append(i)
            pivot += 1
        return r, pivotlist

    def rref(self):
        "rref"
        A = self.copy()
        for j in xrange(A.rows):
            rlist = A.nz_col_lower(j)
            if A[j, j] == 0:
                if rlist:
                    A.row_swap(j, rlist[0])
                    rlist.pop(0)
                else:
                    continue
            A.row_scale(j, 1 /  A[j, j])
            for i in A.nz_col(j):
                if i != j:
                    A.row_add(i, j, - A[i, j])
        return A
                    
    def nz_col(self, j):
        li = []
        for i in xrange(self.rows):
            if self[i, j] != 0:
                li.append(i)
        return li

    def nz_col_lower(self, j):
        " Returns the row indices of non-zero elements in column j, below the diagonal"
        li = []
        for i in xrange(j + 1, self.rows):
            if self[i, j] != 0:
                li.append(i)
        return li

    def is_upper(self):
        return all(j >= i for i in xrange(self.rows) for j, _ in self.mat[i])
            
    def applyfunc(self, f):
        for i in xrange(self.rows):
            for ind, (j, value) in enumerate(self.mat[i]):
                self.mat[i][ind] = (self.mat[i][ind][0], f(self.mat[i][ind][1]))

    def sparsity(self):
        return float(self.nnz()) / (self.rows * self.cols)

    def nnz(self):
        return sum(len(self.mat[i]) for i in xrange(self.rows))

    def join_rows(A, B):
        A = A.copy()
        B = B.copy()
        assert A.rows == B.rows
        for i in xrange(B.rows):
            for ind, (j, val) in enumerate(B.mat[i]):
                B.mat[i][ind] = (j + A.cols, val)
            A.mat[i].extend(B.mat[i])
        A.cols += B.cols
        return A

    @classmethod
    def eye(cls, n, one = 1, zero = 0):
        return cls(n, n, lambda i, j: one if i==j else zero)

    def inv_rref(self):
        aug = self.join_rows(LILMatrix.eye(self.rows, one = 1, zero = 0))
        reduced = aug.rref()
        return reduced[:,self.rows:]

    def copy(self):
        # Could be better
        return LILMatrix(self.rows, self.cols, lambda i, j: self[i, j])
        
    def __add__(self, other):
        if not isinstance(other, LILMatrix):
            if self.is_square():
                other = other * LILMatrix.eye(self.rows)
            else:
                raise Exception
        add = self.copy()
        for i in xrange(self.rows):
            add.mat[i] = _row_add(self.mat[i], other.mat[i], 1)
        return add

    def is_square(self):
        return self.rows == self.cols

    def transpose(self):
        T = LILMatrix.zeros(self.rows)
        for i in xrange(self.rows):
            for j, value in self.mat[i]:
                T.mat[j].append((i, value))
        return T
    
    T = property(transpose)

    @property
    def shape(self):
        return (self.rows, self.cols)
    
    def __eq__(self, other):
        if not self.shape == other.shape:
            return False
        return all(self.mat[i][ind] == other.mat[i][ind] for i in xrange(self.rows) for ind in xrange(len(self.mat[i]))) 

    def __ne__(self, other):
        if not self.shape == other.shape:
            return True
        return any(self.mat[i][ind] != other.mat[i][ind] for i in xrange(self.rows) for ind in xrange(len(self.mat[i]))) 

    def doolittle(self):
        A = self.clone()
        n = A.rows
        for i in xrange(n):
            for j in xrange(i):
                a = A[i, j]
                for p in xrange(j):
                    a -= A[i, p] * A[p, j]
                A[i, j] = a / A[j, j]
            for j in xrange(i, n):
                a = A[i, j]
                for p in xrange(i):
                    a -= A[i, p] * A[p, j]
                A[i, j] = a
        return A

    def LU_sparse(self):

        cached_result = self._get_cache('LUP')
        if cached_result and self._cached:
            return cached_result

        row_swaps = []
        A = self.clone()
        n = self.rows
        for k in xrange(n):

            rlist = A.nz_col_lower(k)

            # Pivoting
            if A[k, k] == 0:
                if not rlist:
                    raise Exception('Singular')
                A.row_swap(k, rlist[0])
                row_swaps.append((k, rlist[0]))
                rlist.pop(0)
            assert A[k, k] != 0

            # Algorithm
            for i in rlist:
                A[i, k] /= A[k, k]

            for j, val in A.mat[k]:
                if j <= k:
                    continue
                # j > k
                for i in rlist:
                    A[i, j] -= A[i, k] * val # i > k, j > k
    
        L = LILMatrix.eye(self.rows)
        for i in xrange(L.rows):
            for j in xrange(i):
                L[i, j] = A[i, j]

        U = LILMatrix.zeros(self.rows)
        for i in xrange(U.rows):
            for j in xrange(i, U.rows):
                U[i, j] = A[i, j]

        self._set_cache('LUP', (L, U, row_swaps))

        return L, U, row_swaps

    def permute(self, row_swaps):
        for r1, r2 in row_swaps:
            self.row_swap(r1, r2)

    def scalar_multiply(self, scalar):
        prod = self.clone()
        for i in xrange(self.rows):
            for ind, (j, value) in enumerate(self.mat[i]):
                prod.mat[i][ind] = (j, scalar * value)
        return prod
        
            
            
                
def _row_add(row1, row2, alpha):
    "li = row1 + alpha * row2 "
    li = []
    i1 = i2 = 0
    n1 = len(row1)
    n2 = len(row2)
    while i1 < n1 or i2 < n2:
        # print i1, i2, len(row1), len(row2)
        if i1 < n1 and (i2 >= n2 or row1[i1][0] < row2[i2][0]):
            li.append(row1[i1])
            i1 += 1
        elif i1 >= n1 or row1[i1][0] > row2[i2][0]:
            li.append((row2[i2][0], alpha * row2[i2][1]))
            i2 += 1
        else:
            val = row1[i1][1] + alpha * row2[i2][1]
            if val != 0:
                li.append((row1[i1][0], val))
            i1 += 1
            i2 += 1
    return li   

def randInvLILMatrix(n, d, min=-5, max=10):
    A = LILMatrix(n, n, lambda i, j: random.randint(min, max) if abs(i - j) <= d-1 else 0)
    return A
    
def iszero(A):
    return all(len(i) == 0 for i in A.mat)

def randLILMatrix(i, j, min=1, max=1, sparsity=0.5):
    return LILMatrix(i, j, lambda i, j: random.randint(min, max) if random.random() < sparsity else 0)
    
