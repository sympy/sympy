"""
Fundamental Arithmetic of Sparse Matrices. The Data structure for Sparse
Matrix is Dictionary of Keys(DOK) and Compressed Sparse Row(CSR).

This is part of the level 0 of the architecture described in [1].

[1] http://www.saurabhjha.me/proposal.html

"""
from sympy import SparseMatrix
from collections import defaultdict


def _doktocsr(dok):
    """Converts a Sparse Matrix to Compressed Sparse Row (CSR) format.

    Parameters
    ==========

    A : contains non-zero elements sorted by key (row, column)
    JA : JA[i] is the column corresponding to A[i]
    IA : IA[i] contains the index in A for the first non-zero element
        of row[i]. Thus IA[i+1] - IA[i] gives number of non-zero
        elements row[i]. The length of IA is always 1 more than the
        number of rows in the matrix.
    """

    try:
        row, JA, A = [list(i) for i in zip(*dok.row_list())]
    except ValueError:
        row, JA, A = [], [], []
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


def add(csr1, csr2, K):
    """
    Adds two Sparse Matrices which are represeneted in CSR format
    using _merge_rowadd. The result is returned in CSR format.

    Examples
    ========
    >>> from sympy.matrices.sparsearith import _doktocsr
    >>> from sympy.matrices.sparsearith import add
    >>> from sympy import ZZ
    >>> from sympy import SparseMatrix

    >>> dok1 = SparseMatrix(2, 2, {(1, 0): ZZ(32), (1, 1): ZZ(19)})
    >>> dok2 = SparseMatrix(2, 2, {(0, 0): 12, (1, 1): ZZ(14)})
    >>> csr1 = _doktocsr(dok1)
    >>> csr2 = _doktocsr(dok2)

    >>> add(csr1, csr2, ZZ)
    [[12, 32, 33], [0, 0, 1], [0, 1, 3], [2, 2]]

    """
    a1, ja1, ia1, shape = csr1
    a2, ja2, ia2, shape = csr2
    a, ia, ja, shape = [], [0], [], shape
    nrow = len(ia1) - 1
    for i in xrange(nrow):
        index1 = slice(ia1[i], ia1[i + 1])
        index2 = slice(ia2[i], ia2[i + 1])
        a_r, ja_r =_merge_rowadd([a1[index1], ja1[index1]], [a2[index2], ja2[index2]], K)
        ia.append(ia[-1] + len(a_r))
        for i in xrange(len(a_r)):
            a.append(a_r[i])
            ja.append(ja_r[i])
    return [a, ja, ia, shape]


def sub(csr1, csr2, K):
    """
    Subtracts two Sparse Matrices which are represeneted in CSR
    format using _merge_rowadd and neg. The result is
    returned in CSR format.

    Examples
    ========
    >>> from sympy.matrices.sparsearith import _doktocsr
    >>> from sympy.matrices.sparsearith import sub
    >>> from sympy import ZZ
    >>> from sympy import SparseMatrix

    >>> dok1 = SparseMatrix(2, 2, {(1, 0): ZZ(32), (1, 1): ZZ(19)})
    >>> dok2 = SparseMatrix(2, 2, {(0, 0): ZZ(12), (1, 1): ZZ(14)})
    >>> csr1 = _doktocsr(dok1)
    >>> csr2 = _doktocsr(dok2)

    >>> sub(csr1, csr2, ZZ)
    [[-12, 32, 5], [0, 0, 1], [0, 1, 3], [2, 2]]

    """
    return add(csr1, neg(csr2, K), K)


def _merge_rowadd(row1, row2, K):
    """
    Helper function of add and sub. The "merge" in the function name
    emphasis the fact that it's working is similar to merge sort.

    """
    a1, ja1 = row1
    a2, ja2 = row2
    m1, m2 = len(row1[0]), len(row2[0])
    k1, k2 = 0, 0
    a_r, ja_r = [], []
    while k1 < m1 and k2 < m2:
        if ja1[k1] < ja2[k2]:
            a_r.append(a1[k1]), ja_r.append(ja1[k1])
            k1 += 1
        elif ja1[k1] > ja2[k2]:
            a_r.append(a2[k2]), ja_r.append(ja2[k2])
            k2 += 1
        else:
            a_r.append(a1[k1] + a2[k2]), ja_r.append(ja1[k1])
            k1 += 1
            k2 += 1
    if k1 == m1:
        for a2, ja2 in zip(a2[k2:], ja2[k2:]):
            a_r.append(a2), ja_r.append(ja2)
    else:
        for a1, ja1 in zip(a1[k1:], ja1[k1:]):
            a_r.append(a1), ja_r.append(ja1)
    return a_r, ja_r


def neg(csr, K):
    """
    Negates the elements of a Sparse Matrix given in CSR format. The
    result is returned in CSR format.

    Examples
    ========

    >>> from sympy.matrices.sparsearith import neg, _doktocsr
    >>> from sympy import ZZ
    >>> from sympy import SparseMatrix

    >>> dok1 = SparseMatrix(2, 2, {(0, 0): ZZ(12), (1, 1): ZZ(14)})
    >>> csr1 = _doktocsr(dok1)

    >>> neg(csr1, ZZ)
    [[-12, -14], [0, 1], [0, 1, 2], [2, 2]]

    """
    a, ja, ia, shape = csr
    a_new = []
    for i in xrange(len(a)):
        a_new.append(-a[i])
    return [a_new, ja, ia, shape]


def mulspvec(csr, vec, K):
    """
    Performs the Sparse Matrix-Vector calculation. The vector is assumed to be
    dense matrix in list of lists format. The sparse matrix is represented in
    CSR format.

    The result is returned in DOK
    format.

    Examples
    ========

    >>> from sympy.matrices.sparsearith import mulspvec
    >>> from sympy.matrices.sparsearith import _doktocsr
    >>> from sympy import ZZ
    >>> from sympy import SparseMatrix

    >>> dok1 = SparseMatrix(2, 2, {(0, 0): ZZ(12), (1, 1): ZZ(14)})
    >>> csr1 = _doktocsr(dok1)
    >>> vec = [[ZZ(1)], [ZZ(2)]]

    >>> mulspvec(csr1, vec, ZZ)
    [12]
    [28]

    """
    a, ja, ia, shape = csr
    smat = {}
    for i in xrange(len(ia) - 1):
        stripe = slice(ia[i], ia[i + 1])
        for m, n in zip(a[stripe], ja[stripe]):
            try:
                smat[i, 0] = smat[i, 0] + m*vec[n][0]
            except KeyError:
                smat[i, 0] = K.zero + m*vec[n][0]
    return SparseMatrix(shape[0], 1, smat)


def mulspsp(dok1, dok2, K):
    """
    Performs the Sparse Matrix-Sparse Matrix calculation using _mulrowcol.
    Both matrices are represented in DOK format.

    The result is returned in DOK format.

    Examples
    ========

    >>> from sympy.matrices.sparsearith import mulspsp
    >>> from sympy import SparseMatrix
    >>> from sympy import ZZ

    >>> dok1 = SparseMatrix(2, 2, {(0, 0): ZZ(12), (1, 1): ZZ(14)})
    >>> dok2 = SparseMatrix(2, 2, {(1, 0): 32, (1, 1): 19})

    >>> mulspsp(dok1, dok2, ZZ)
    [  0,   0]
    [448, 266]

    """
    nrow, ncol = dok1.rows, dok2.cols
    smat = {}
    for i in xrange(nrow):
        row = dok1.row(i)
        row_smat = row._smat
        for j in xrange(ncol):
            col = dok2.col(j)
            col_smat = col._smat
            smat[i, j] = _mulrowcol(row, col, K)
    return SparseMatrix(nrow, ncol, smat)


def _mulrowcol(row, col, K):
    """
    Helper function of mulspsp. Multiplies a row and a column.
    The row and column is represented in DOK format.

    """
    smat1, smat2 = row._smat, col._smat
    keys1, keys2 = sorted(smat1.keys()), sorted(smat2.keys())
    result = K.zero
    for i, j in keys1:
        k = _binsearch((j, 0), keys2, 0, len(keys2) - 1)
        if k != None:
            result += row[i, j]*col[keys2[k]]
    return result


def _binsearch(element, array, beg, end):
    """
    An implementation of binary search. It is expected to be
    used as follows for an element i and list v--
    _binsearch(i, v, 0, len(v) - 1).

    """
    if beg > end:
        return None
    mid = (beg + end)/2
    if element < array[mid]:
        return _binsearch(element, array, beg, mid - 1)
    elif element > array[mid]:
        return _binsearch(element, array, mid + 1, end)
    else:
        return mid
