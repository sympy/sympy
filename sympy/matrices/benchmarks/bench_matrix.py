from __future__ import print_function, division

from sympy import eye, zeros, Integer, Matrix

i3 = Integer(3)
M = eye(100)
A = Matrix(10,10,range(100))

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
