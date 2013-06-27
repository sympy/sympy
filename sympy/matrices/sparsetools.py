from sympy import SparseMatrix
from collections import defaultdict


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


def _mulscsp(v, csr, K):
    a, ja, ia, shape = csr
    return [[K(v)*i for i in a], ja, ia, shape]


def _mulspsp(sp1, sp2, K):
    csr1, csr2 = _doktocsr(sp1), _doktocsr(sp2)
    a1, ja1, ia1, shape1 = csr1
    a2, ja2, ia2, shape2 = csr2
    sp2cols = list(set(ja2))
    smat = {}
    i = 0
    while i < len(ia1) - 1:
        stripe = slice(ia1[i], ia1[i + 1])
        csr_row = [a1[stripe], ja1[stripe], [ia1[i] - ia1[i], ia1[i + 1] - ia1[i]], [1, sp1.cols]]
        j = 0
        while j < sp2.cols:
            if _binsearch(j, sp2cols, 0, len(sp2cols)) or _binsearch(j, sp2cols, 0, len(sp2cols)) == 0:
                smat[i, j] = _mulspvec(csr_row, sp2.col(j), K)[0]
            j = j + 1
        i = i + 1
    return SparseMatrix(sp1.rows, sp2.cols, smat)


def _mulspsp2(sp1, sp2, K):
    i = 0
    smat = {}
    while i < sp1.rows:
        if len(sp1.row(i)._smat) != 0:
            crnt_row = _doktocsr(sp1.row(i))
            j = 0
            while j < sp2.cols:
                if len(sp2.col(j)._smat) != 0:
                    crnt_col = sp2.col(j)
                    smat[i, j] = _mulspvec(crnt_row, crnt_col, K)[0]
                j = j + 1
        i = i + 1
    return SparseMatrix(sp1.rows, sp2.cols, smat)


def _mulspsp3(self, other, K):
    A = self
    B = other
    # sort B's row_list into list of rows
    Blist = [[] for i in range(B.rows)]
    for i, j, v in B.row_list():
        Blist[i].append((j, v))
    Cdict = defaultdict(int)
    for k, j, Akj in A.row_list():
        for n, Bjn in Blist[j]:
            temp = K(Akj)*K(Bjn)
            Cdict[k, n] += temp
    rv = self.zeros(A.rows, B.cols)
    rv._smat = dict([(k, v) for k, v in Cdict.iteritems() if v])
    return rv


def _mulspvec(csr, vec, K):
    a, ja, ia, shape = csr
    smat = {}
    for i in range(len(ia) - 1):
        stripe = slice(ia[i], ia[i + 1])
        for m, n in zip(a[stripe], ja[stripe]):
            try:
                smat[i, 0] = K(smat[i, 0]) +  K(m)*K(vec[n])
            except KeyError:
                smat[i, 0] = K.zero +  m*vec[n]
    return SparseMatrix(shape[0], 1, smat)


def _applyfunc(csr, f, K):
    if not callable(f):
        return TypeError("f must be callable")
    else:
        a, ja, ia, shape = csr
        return [[f(i) for i in a], ja, ia, shape]


def _binsearch(i, v, beg, end):
    if beg > end:
        return False
    else:
        mid = (beg + end)/2
        if v[mid] == i:
            return mid
        elif v[mid] > i:
            return _binsearch(i, v, beg, mid - 1)
        elif v[mid] < i:
            return _binsearch(i, v, mid + 1, end)
