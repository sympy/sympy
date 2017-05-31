from __future__ import print_function

from sympy.core import S
from sympy.matrices.dense import DenseMatrix
from sympy.matrices.matrices import classof


class PolyMatrix(DenseMatrix):
    """
    Matrix of polynomials.

    >>> from sympy.polys.polymatric import PolyMatrix
    >>> from sympy.matrices.matrices import Matrix
    >>> from sympy import Symbol
    >>> pm1 = PolyMatrix([[Poly(x**2, x), Poly(-x, x)], [Poly(x**3, x), Poly(-1 + x, x)]])
    >>> v1 = Matrix([[1, 0], [-1, 0]])
    >>> x = Symbol('x')
    >>> pm1*v1
    Matrix([
    [    Poly(x**2 + x, x, domain='ZZ'), Poly(0, x, domain='ZZ')],
    [Poly(x**3 - x + 1, x, domain='ZZ'), Poly(0, x, domain='ZZ')]])

    >>> pm2 = PolyMatrix([[Poly(x**2, x, domain='QQ'), Poly(0, x, domain='QQ'), Poly(-x**2, x, domain='QQ'),
            Poly(x**3, x, domain='QQ'), Poly(0, x, domain='QQ'), Poly(-x**3, x, domain='QQ')]])
    >>> v2 = Matrix([1, 0, 0, 0, 0, 0])
    >>> pm2*v2
    Matrix([[Poly(x**2, x, domain='QQ')]])

    """
    _class_priority = 5

    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    @classmethod
    def _new(cls, *args, **kwargs):
        if kwargs.get('copy', True) is False:
            if len(args) != 3:
                 raise TypeError("'copy=False' requires a matrix be initialized as rows,cols,[list]")
            rows, cols, flat_list = args
        else:
            rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)
            flat_list = list(flat_list) # create a shallow copy
        self = object.__new__(cls)
        self.rows = rows
        self.cols = cols
        self._mat = flat_list
        return self

    def _eval_matrix_mul(self, other):
        self_rows, self_cols = self.rows, self.cols
        other_rows, other_cols = other.rows, other.cols
        other_len = other_rows * other_cols
        new_mat_rows = self.rows
        new_mat_cols = other.cols

        new_mat = [S.Zero]*new_mat_rows*new_mat_cols

        if self.cols != 0 and other.rows != 0:
            mat = self._mat
            other_mat = other._mat
            for i in range(len(new_mat)):
                row, col = i // new_mat_cols, i % new_mat_cols
                row_indices = range(self_cols*row, self_cols*(row+1))
                col_indices = range(col, other_len, other_cols)
                vec = (mat[a]*other_mat[b] for a,b in zip(row_indices, col_indices))
                new_mat[i] = sum(vec)

        return classof(self, other)._new(new_mat_rows, new_mat_cols, new_mat, copy=False)
