from sympy import SparseMatrix


def _doktocsr(dok):
    """Converts a sparse matrix to Compressed Sparse Row (CSR) format.

    Parameters
    ==========

    A : contains non-zero elements sorted by key (row, column)
    JA : JA[i] is the column corresponding to A[i]
    IA : IA[i] contains the index in A for the first non-zero element
        of row[i]. Thus IA[i+1] - IA[i] gives number of non-zero
        elements row[i]. The length of IA is always 1 more than the
        number of rows in the matrix.
    """
    row, JA, A = [list(i) for i in zip(*dok.row_list())]
    IA = [0]*((row[0] if row else 0) + 1)
    for i, r in enumerate(row):
        IA.extend([i]*(r - row[i - 1]))  # if i = 0 nothing is extended
    IA.extend([len(A)]*(dok.rows - len(IA) + 1))
    shape = [dok.rows, dok.cols]
    return [A, JA, IA, shape]


def _csrtodok(csr):
    """Converts a CSR representation to DOK representation"""
    smat = {}
    A, JA, IA, shape = csr
    for i in range(len(IA) - 1):
        indices = slice(IA[i], IA[i + 1])
        for l, m in zip(A[indices], JA[indices]):
            smat[i, m] = l
    return SparseMatrix(*(shape + [smat]))


def _add(sp1, sp2, K):
    if sp1.nnz() > sp2.nnz():
        temp = sp1
        sp1 = sp2
        sp2 = temp
    keys1 = sorted(sp1._smat.keys())
    keys2 = sorted(sp2._smat.keys())
    result = {}
    i = 0
    for k in keys1:
        if _binsearch(k, keys2, 0, len(keys2)) or _binsearch(k, keys2, 0, len(keys2)) == 0:
            l = keys2[_binsearch(k, keys2, 0, len(keys2))]
            result[k] = K(sp1[k]) + K(sp2[l])
            keys2.remove(l)
        else:
            result[k] = K(sp1[k]) 
    for k in keys2:
        result[k] = K(sp2[k])
    return SparseMatrix(sp1.rows, sp1.cols, result)


def _sub(sp1, sp2, K):
    return add(sp1, (- sp2), K)
    

def _mulspsp(sp1, sp2, K):
    pass
    
        

def _mulspvec(csr, vec, K):
    a, ja, ia, shape = csr
    smat = {}
    for i in range(len(ia) - 1):
        stripe = slice(ia[i], ia[i + 1])
        for m, n in zip(a[stripe], ja[stripe]):
            try:
                smat[i, 0] = smat[i, 0] +  m*vec[n]
            except KeyError:
                smat[i, 0] = K.zero +  m*vec[n]
    return SparseMatrix(shape[0], 1, smat)


def _binsearch(i, v, beg, end):
    if beg > end or i > v[end - 1]:
        return False
    else:
        mid = (beg + end)/2
        if i == v[mid]:
            return mid
        elif i < v[mid]:
            return _binsearch(i, v, beg, mid - 1)
        elif i > v[mid]:
            return _binsearch(i, v, mid + 1, end)


def _invtuple(t):
    t = list(t)
    temp = t[0]
    t[0] = t[1]
    t[1] = temp
    return tuple(t)
