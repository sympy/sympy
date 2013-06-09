from sympy import SparseMatrix
def doktocsr(sparse):
    """
    Converts a sparse matrix to Compressed Sparse Row(CSR) format.
    """
    A = [] #array of non-zero element values
    IA = [] #array of index of first non-zero elements of row i
    JA = [] #array of column index of each A element
    [A.append(e[-1]) for e in sparse.row_list()]
    a = sparse.row_list()
    j = 0
    check = - 1
    elements = []
    while j < sparse.nnz():
        if check !=sparse.row_list()[j : ][0][0]:
            elements.append(sparse.row_list()[j :][0][-1])
            check = sparse.row_list()[j : ][0][0]
        j = j + 1
    i = 0
    j = 0
    while i < len(elements):
        if A[j] == elements[i]:
            IA.append(j)
            i = i + 1
        j = j + 1
    IA.append(len(A))
    [JA.append(e[1]) for e in sparse.row_list()]
    return [A, JA, IA]
