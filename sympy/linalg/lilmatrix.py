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

    def __getitem__(self, key):
        i, j = key
        if type(i) is slice or type(j) is slice:
                return self.submatrix(key)
        
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

    def submatrix(self, keys):
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

    def copyin_matrix(self, key, value):
        from matrixutils import _slice_to_bounds
        rlo, rhi = _slice_to_bounds(self, key[0], self.rows)
        clo, chi = _slice_to_bounds(self, key[1], self.cols)
        lilrepr = value.mat
        for i in xrange(len(value.mat)):
            self.mat[rlo + i].extend([(j + clo, val) for j, val in value.mat[i]])
            self.mat[rlo + i].sort()

    def row(self, r, f):
        for i, (j, val) in enumerate(self.mat[r]):
            self.mat[r][i] = j, f(val, j)

    def row_swap(self, r1, r2):
        self.mat[r1], self.mat[r2] = self.mat[r2], self.mat[r1]

    def to_densematrix(self):
        return Matrix(self.rows, self.cols, lambda i, j: self[i, j])

    def to_dokmatrix(self):
        from sympy.linalg import DOKMatrix
        mat = DOKMatrix(self.rows, self.cols, {})
        for i in xrange(self.rows):
            for j, value in self.mat[i]:
                mat[i, j] = value
        return mat

    def to_lilmatrix(self):
        return self

    def solve(self, rhs, method="GE"):
        if method == "GE":
            from lilmatrix_tools import solve_gauss
            return solve_gauss(self, rhs)
        elif method == "RREF":
            from lilmatrix_tools import solve_rref
            return solve_rref(self, rhs)
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

    def det(self, method='GE'):
        if method == "GE":
            from lilmatrix_tools import det_gauss
            return det_gauss(self)
        elif method == 'LU':
            from lilmatrix_tools import det_LU
            return det_LU(self)
        else:
            raise Exception

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
    
    def __eq__(self, other):
        if not self.shape == other.shape:
            return False
        return all(self.mat[i][ind] == other.mat[i][ind] for i in xrange(self.rows) for ind in xrange(len(self.mat[i]))) 

    def __ne__(self, other):
        if not self.shape == other.shape:
            return True
        return any(self.mat[i][ind] != other.mat[i][ind] for i in xrange(self.rows) for ind in xrange(len(self.mat[i])))            
                
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
    
