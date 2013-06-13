from sympy import SparseMatrix
from collections import defaultdict
from sympy import oo
def doktocsr(sparse):
    """
    Converts a sparse matrix to Compressed Sparse Row(CSR) format.
    """
    A = [] #array of non-zero element values
    IA = [] #array of index of first non-zero elements of row i
    JA = [] #array of column index of each A element
    row, JA, A = [list(i) for i in zip(*sparse.row_list())]
    d = row[0] - 0
    while d >= 0:
        IA.append(0)
        d = d - 1
    i = 1
    while i < len(row):
        if row[i] != row[i -1]:
            d = row[i] - row[i - 1]
            while d > 0:
                IA.append(i)
                d = d - 1
        i = i + 1
    IA.append(len(row))
    return [A, JA, IA]
