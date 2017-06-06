from __future__ import print_function, division

from sympy.core.compatibility import range
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


def add(csr1, csr2):
    a1, ja1, ia1, shape = csr1
    a2, ja2, ia2, shape = csr2
    length = len(ia1)
    a_res = []
    ja_res = []
    ia_res = [0]
    for i in range(length - 1):
        stripe1 = slice(ia1[i], ia1[i + 1])
        stripe2 = slice(ia2[i], ia2[i + 1])
        a1_row, ja1_row = a1[stripe1], ja1[stripe1]
        a2_row, ja2_row = a2[stripe2], ja2[stripe2]
        ia_res.append(ia_res[-1])
        for m in range(len(ja1_row)):
            if ja1_row[m] in ja2_row and ja1_row:
                n = ja2_row.index(ja1_row[m])
                a_res.append(a1_row[m] + a2_row[n])
                ja_res.append(ja1_row[m])
                ia_res[-1] = ia_res[-1] + 1
            else:
                a_res.append(a1_row[m])
                ja_res.append(ja1_row[m])
                ia_res[-1] = ia_res[-1] + 1
        for k in range(len(ja2_row)):
            if ja2_row[k] not in ja1_row:
                a_res.append(a2_row[k])
                ja_res.append(ja2_row[k])
                ia_res[-1] = ia_res[-1] + 1
    return [a_res, ja_res, ia_res , shape]


def sub(csr1, csr2):
    return sparseadd(csr1, -csr2)
