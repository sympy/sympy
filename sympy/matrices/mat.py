from sympy import SparseMatrix
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
    IA = [0]*(row[0] + 1)
    i = 1
    while i < len(row):
        if row[i] != row[i -1]:
            IA.extend([i]*(row[i] - row[i - 1]))
        i = i + 1
    IA.extend([len(A)]*(sparse.rows - len(IA) + 1))
    print IA
