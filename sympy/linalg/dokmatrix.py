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

import matrixutils

import random



class DOKMatrix(DataMatrix):
    """
        DOKMatrix, internal low-level sparse list of list representation matrix
        self.rows ---> number of rows
        self.cols ---> number of cols
        self.mat  ---> data stored in a dictionary of keys representation matrix
                       self.mat[i, j] denotes the (i, j) element if it is non-zero.
        """

    def __init__(self, *args, **kwargs):
        if len(args) == 3:
            rows = args[0]
            cols = args[1]
            mat = args[2]

            if isinstance(mat, dict):
                mat = matrixutils._dokrepr_from_dict(rows, cols, mat)
            elif callable(mat):
                mat = matrixutils._dokrepr_from_callable(rows, cols, mat)
            elif isinstance(mat, (list, tuple)) and len(mat) == rows * cols:
                mat = matrixutils._dokrepr_from_list(rows, cols, mat)
            else:
                raise TypeError('Data type not understood.')
        elif len(args) == 1:
            mat = args[0]
            if isinstance(mat, DataMatrix):
                matrix = mat.to_dokmatrix()
                rows = matrix.rows
                col = matrix.cols
                mat = matrix.mat
            elif isinstance(mat, (list, tuple)):
                rows = len(mat)
                cols = len(mat[0])
                mat = matrixutils._dokrepr_from_lil(rows, cols, mat)
            else:
                raise TypeError('Data type not understood.')
        elif len(args) == 0:
            rows = cols = 0
            mat = []

        self.rows = rows
        self.cols = cols
        self.mat = mat

    def __getitem__(self, key):
        return self.mat.get(key, 0)

    def __setitem__(self, key, value):
        if value != 0:
            self.mat[key] == value
        else:
            del self.mat[key]

    def copyin_matrix(self, key, value):
        rlo, rhi = slice_to_bounds(key[0], self.rows)
        clo, chi = slice_to_bounds(key[1], self.cols)

        for i in range(value.rows):
            for j in range(value.cols):
                self[i+rlo, j+clo] = sympify(value[i,j])

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
        M = self.copy()
        for i in other.mat:
            if i in M.mat:
                M[i] += other[i]
            else:
                M[i] = other[i]
        return M

    def __mul__(self, other):
        if isinstance(other, DOKMatrix):
            rows1 = _lil_row_major(A)
            rows2 = _lil_row_major(B)
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
        else:
            C = DOKMatrix(matrix.rows, matrix.cols, {})
            if scalar != 0:
                for i in matrix.mat:
                    C.mat[i] = scalar * matrix.mat[i]
            return C

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

    def submatrix(self, keys):
        if not isinstance(keys[0], slice) and not isinstance(keys[1], slice):
            raise TypeError("Both elements of `keys` must be slice objects.")
        rlo, rhi = self.slice2bounds(keys[0], self.rows)
        clo, chi = self.slice2bounds(keys[1], self.cols)
        if not ( 0<=rlo<=rhi and 0<=clo<=chi ):
            raise IndexError("Slice indices out of range: a[%s]"%repr(keys))
        return DOKMatrix(rhi-rlo, chi-clo, lambda i,j: self[i+rlo, j+clo])

    @property
    def shape(self):
        return (self.rows, self.cols)

    def __len__(self):
        return self.rows*self.cols

    def sparsity(self):
        return len(self.mat)/len(self)

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

    def applyfunc(self, op):
        for i in self.mat:
            self.mat[i] = op(self.mat[i])

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

def _lil(keys, n, j):
    keys.sort()
    zeros = set(xrange(n)) - set(k[j] for k in keys)
    lil = [ [] if i in zeros 
        else [k for k in keys if k[j] == i ]
        for i in xrange(n) ]
    return lil

lil = _lil


