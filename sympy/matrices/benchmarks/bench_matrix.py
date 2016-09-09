from __future__ import print_function, division

from sympy import eye, zeros, Integer, Matrix, SparseMatrix

i3 = Integer(3)
M = eye(100)
A = Matrix(10,10,range(100))
MS = SparseMatrix(200, 200, {(3, 15): 1, (6, 9): 1, (10, 11): 1, (8, 0): 1, (16, 2): 1, (5, 19): 1, (18, 19): 1, (2, 9): 1, (16, 14): 1, (16, 13): 1, (15, 19): 1, (10, 3): 1, (7, 2): 1, (1, 2): 1, (8, 14): 1, (9, 0): 1, (2, 19): 1, (0, 15): 1, (10, 7): 1, (8, 5): 1, (7, 6): 1, (14, 18): 1, (8, 10): 1, (19, 16): 1, (3, 14): 1, (18, 9): 1, (2, 2): 1, (15, 4): 1, (19, 10): 1, (7, 9): 1, (15, 15): 1, (6, 4): 1, (16, 11): 1, (10, 6): 1, (8, 2): 1, (17, 10): 1, (17, 9): 1, (10, 19): 1, (19, 18): 1, (14, 6): 1, (14, 15): 1, (16, 0): 1, (0, 5): 1, (3, 19): 1, (1, 0): 1, (0, 8): 1, (9, 6): 1, (18, 5): 1, (3, 5): 1, (14, 11): 1, (11, 7): 1, (16, 5): 1, (4, 6): 1, (11, 8): 1, (15, 16): 1, (18, 10): 1, (3, 1): 1, (2, 11): 1, (9, 9): 1, (0, 6): 1, (19, 12): 1, (18, 11): 1, (7, 17): 1, (13, 9): 1, (3, 4): 1, (17, 17): 1, (9, 12): 1, (15, 7): 1})

def timeit_Matrix__getitem_ii():
    M[3, 3]


def timeit_Matrix__getitem_II():
    M[i3, i3]


def timeit_Matrix__getslice():
    M[:, :]


def timeit_Matrix_zeronm():
    zeros(100, 100)

def timeit_Matrix_extract():
    A.extract(range(0,10,2), range(0,10,3))

def timeit_Matrix_mul():
    A*A

def timeit_Matrix_scalar_mul():
    i3*A

def timeit_Matrix_add():
    A+A

def timeit_SpareMatrix_mul():
    MS*MS
