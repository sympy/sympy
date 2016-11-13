from __future__ import print_function, division
from sympy import eye, zeros, Integer, randMatrix
import time
import itertools
import pandas

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

output_times = pandas.DataFrame(columns=["n", "l", "time_recursion", "time_jordan"])

def timeit_Matrix_pow_I():
    # n is the size of the matrix
    nlist = [2, 3, 4, 5]
    l_list = [25*2**i for i in range(20)]
    repetitions = range(7)
    for i, (n, l) in enumerate(itertools.product(nlist, l_list)):
        tr = 0.0
        tj = 0.0
        for j in repetitions:
            A = randMatrix(n)
            t0 = time.time()
            A._matrix_pow_by_recursion(l)
            t1 = time.time()
            A._matrix_pow_by_jordan_blocks(l)
            t2 = time.time()
            tr = t1 - t0
            tj = t2 - t1
        tr /= len(repetitions)
        tj /= len(repetitions)
        output_times.loc[i] = [n, l, tr, tj]
