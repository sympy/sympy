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

        for i in mat:
            mat[i] = S(mat[i])
        self.rows = rows
        self.cols = cols
        self.mat = mat

    def __getitem__(self, key):
        "Only element accessing. Slicing disabled for performance."
        return self.mat.get(key, 0)

    def __setitem__(self, key, value):
        "Key cannot be a slice"
        if value != 0:
            self.mat[key] = value
        else:
            del self.mat[key]

    def copyin_matrix(self, key, value):
        from matrixutils import _slice_to_bounds
        rlo, rhi = _slice_to_bounds(key[0], self.rows)
        clo, chi = _slice_to_bounds(key[1], self.cols)

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
        """
        Addition of two matrices.
        """
        M = self.copy()
        for i in other.mat:
            if i in M.mat:
                M[i] += other[i]
            else:
                M[i] = other[i]
        return M

    def __mul__(self, other):
        """
        Multiplication a matrix with another matrix or a scalar.
        """
        if isinstance(other, DOKMatrix):
            from dokmatrix_tools import _lil_row_major
            rows1 = _lil_row_major(self)
            rows2 = _lil_row_major(other)
            Cdict = {}
            for k in xrange(self.rows):
                for _, j in rows1[k]:
                    for _, n  in rows2[j]:
                        temp = self[k, j] * other[j, n]
                        if (k, n) in Cdict:
                            Cdict[k, n] += temp
                        else:
                            Cdict[k, n] = temp
            C = DOKMatrix(self.rows, other.cols, Cdict)
            if C.shape == (1, 1):
                return C[0, 0]
            return C
        else:
            "other is a scalar"
            C = DOKMatrix(self.rows, self.cols, {})
            if B != 0:
                for i in self.mat:
                    C.mat[i] = B * self.mat[i]
            return C

    def __neg__(self):
        return -1 * self

    def __sub__(self,a):
        return self + (-a)
    
    def __rmul__(self, other):
        return self.__mul__(other)  

    def __eq__(self, other):
        if not isinstance(other, DOKMatrix):
            return False
        if self.rows != other.rows or self.cols != other.cols:
            return False
        for i in xrange(self.rows):
            for j in xrange(self.cols):
                if self[i, j] != other[i, j]:
                    return False
        return True

    def submatrix(self, keys):
        from matrixutils import _slice_to_bounds
        rlo, rhi = _slice_to_bounds(keys[0], self.rows)
        clo, chi = _slice_to_bounds(keys[1], self.cols)
        return DOKMatrix(rhi-rlo, chi-clo, lambda i,j: self[i+rlo, j+clo])

    @property
    def shape(self):
        return (self.rows, self.cols)

    def __len__(self):
        return self.rows * self.cols

    def sparsity(self):
        return float(len(self.mat)) / len(self)

    def is_symmetric(self):
        return all(self[i, j] == self[j, i] for i, j in self.mat.keys())

    def is_square(self):
        return self.rows == self.cols

    def has(self, expr):
        any(self[i, j].has(expr) for i, j in self.mat.keys())

    def solve(self, rhs, method="LDL"):
        """
        Solves self * X = rhs using method specified.
        
        LDL     ---     LDL decomposition
        CH      ---     Cholesky decomposition
        """
        if method == "LDL":
            from dokmatrix_tools import _LDLsolve
            return _LDLsolve(self, rhs)
        elif method == "CH":
            from dokmatrix_tools import _cholesky_solve
            return _cholesky_solve(self, rhs)
        else:
            raise ValueError("Unrecognized method")

    def applyfunc(self, op):
        "Applies the function op to all non-zero elements of the matrix"
        for i in self.mat:
            self.mat[i] = op(self.mat[i])

    def to_dokmatrix(self):
        return self

    def to_lilmatrix(self):
        from sympy.linalg.lilmatrix import LILMatrix
        return LILMatrix._from_dict(self.rows, self.cols, self.mat)

    def to_densematrix(self):
        from sympy.linalg.densematrix import DenseMatrix
        return DenseMatrix(self.rows, self.cols, lambda i, j: self[i, j])

    def copy(self):
        return DOKMatrix(self.rows, self.cols, self.mat)

