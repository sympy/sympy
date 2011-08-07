from __future__ import division
import random
from sympy.printing import sstr, pretty
from sympy.simplify.simplify import simplify as sympy_simplify
from sympy.core.singleton import S
from datamatrix import DataMatrix
import matrixutils

def _iszero(x):
    return x == 0

class LILMatrix(DataMatrix):
    """
        LILMatrix, internal low-level sparse list of list representation matrix
        self.rows ---> number of rows
        self.cols ---> number of cols
        self.mat  ---> data stored in a list of rows, where each row i is a
                       list of non-zero elements in the form of tuple (j, value)
        """

    def __init__(self, *args, **kwargs):
        if len(args) == 3:
            rows = args[0]
            cols = args[1]
            mat = args[2]

            if isinstance(mat, dict):
                mat = matrixutils._lilrepr_from_dict(rows, cols, mat)
            elif callable(mat):
                mat = matrixutils._lilrepr_from_callable(rows, cols, mat)
            elif isinstance(mat, (list, tuple)) and len(mat) == rows * cols:
                mat = matrixutils._lilrepr_from_list(rows, cols, mat)
            else:
                raise TypeError('Data type not understood.')
        elif len(args) == 1:
            mat = args[0]
            if isinstance(mat, DataMatrix):
                matrix = mat.to_lilmatrix()
                rows = matrix.rows
                col = matrix.cols
                mat = matrix.mat
            elif isinstance(mat, (list, tuple)):
                rows = len(mat)
                cols = len(mat[0])
                mat = matrixutils._lilrepr_from_lil(rows, cols, mat)
            else:
                raise TypeError('Data type not understood.')
        elif len(args) == 0:
            rows = cols = 0
            mat = []

        self.rows = rows
        self.cols = cols
        self.mat = mat

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

    def __eq__(self, other):
        if not self.shape == other.shape:
            return False
        try:
            return all(self.mat[i][ind] == other.mat[i][ind]
                for i in xrange(self.rows)
                for ind in xrange(len(self.mat[i])))
        except:
            return False

    def __ne__(self, other):
        if not self.shape == other.shape:
            return True
        return any(self.mat[i][ind] != other.mat[i][ind]
            for i in xrange(self.rows)
            for ind in xrange(len(self.mat[i])))

    def __add__(self, other):
        from lilmatrix_tools import _row_add
        if not isinstance(other, LILMatrix):
            if self.is_square():
                other = other * LILMatrix.eye(self.rows)
            else:
                raise Exception
        add = self.copy()
        for i in xrange(self.rows):
            add.mat[i] = _row_add(self.mat[i], other.mat[i], 1)
        return add

    def __sub__(self, other):
        return self + (-1 * other)
    
    def __mul__(self, other):
        if isinstance(other, DataMatrix):
            prod = self.to_dokmatrix() * other.to_dokmatrix()
            return prod.to_lilmatrix() 
        else:
            # Scalar multiplication
            if other == 0:
                from lilmatrix_tools import zeros
                return zeros(self.rows, self.cols)
            prod = self.copy()
            for i in xrange(self.rows):
                for ind, (j, value) in enumerate(self.mat[i]):
                    prod.mat[i][ind] = (j, other * value)
            return prod

    def __rmul__(self, other):
        ## assume other is scalar
        return self.__mul__(other)

    def submatrix(self, *args):
        from matrixutils import _slice_to_bounds
        from lilmatrix_tools import zeros
        if len(args) == 1:
            keys = args[0]
            rlo, rhi = _slice_to_bounds(keys[0], self.rows)
            clo, chi = _slice_to_bounds(keys[1], self.cols)
        elif len(args) == 4:
            rlo, rhi, clo, chi = args
        else:
            raise TypeError("Data type not understood")

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
        A = zeros(outRows, outCols)
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

    def to_densematrix(self):
        from sympy.linalg.densematrix import DenseMatrix
        return DenseMatrix(self.rows, self.cols, lambda i, j: self[i, j])

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

    def is_square(self):
        return self.rows == self.cols

    def transpose(self):
        from lilmatrix_tools import zeros
        T = zeros(self.cols, self.rows)
        for i in xrange(self.rows):
            for j, value in self.mat[i]:
                T.mat[j].append((i, value))
        return T

    T = property(transpose)

    def is_symmetric(self):
        for i in xrange(self.rows):
            for j in xrange(i):
                if self[i, j] != self[j, i]:
                    return False
        return True

    def is_lower(self):
        return all (i >= j for i in xrange(len(self.mat)) for j, val in self.mat[i])

    def permute(self, row_swaps):
        for r1, r2 in row_swaps:
            self.row_swap(r1, r2)


