from __future__ import annotations
from sympy.core.numbers import Integer
from sympy.matrices.dense import eye, zeros, Matrix

i3 = Integer(3)
M = eye(100)


def timeit_Matrix__getitem_ii():
    M[3, 3]


def timeit_Matrix__getitem_II():
    M[i3, i3]


def timeit_Matrix__getslice():
    M[:, :]


def timeit_Matrix_zeronm():
    zeros(100, 100)


class TimeMatrixOps:
    """Benchmarks for common dense matrix operations."""

    def setup(self):
        self.A = (Matrix(5, 5, lambda i, j: 1 if i == j else 0)
                  + Matrix(5, 5, lambda i, j: i + j))
        self.B = Matrix(5, 5, lambda i, j: i - j + 1)

    def time_matrix_multiply(self):
        """Matrix multiplication"""
        self.A * self.B

    def time_matrix_det(self):
        """Determinant computation"""
        self.A.det()

    def time_matrix_inv(self):
        """Matrix inverse"""
        self.A.inv()

    def time_eigenvals(self):
        """Eigenvalue computation"""
        self.A.eigenvals()