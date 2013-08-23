"""
Fundamental Operations of Sparse Matrices. The Data structure for Sparse
Matrix is Dictionary of Keys(DOK) and Compressed Row(CSR).

This is part of the level 0 of the architecture described in [1].

[1] http://www.saurabhjha.me/proposal.html

"""
from sympy.matrices.sparsearith import _binsearch
from sympy import SparseMatrix

def trace(dok, K):
    """
    Returns the trace of a Sparse Matrix. The input and output matrix
    is given in DOK format.

    Examples
    ========

    >>> from sympy.matrices.sparsetools import trace
    >>> from sympy.matrices.sparsearith import _doktocsr
    >>> from sympy import SparseMatrix
    >>> from sympy import ZZ

    >>> d = SparseMatrix(10, 10, {(1, 1): ZZ(12), (3, 5): ZZ(7), (7, 8): ZZ(12)})

    >>> trace(d, ZZ)
    12

    """
    smat = dok._smat
    keys = sorted(smat.keys())
    nrows, ncols = dok.rows, dok.cols
    result = K.zero
    for i in xrange(nrows):
        if _binsearch((i, i), keys, 0, len(keys) - 1) != None:
            result += smat[i, i]
    return result


def conjugate(dok, K):#not working for complex entries
    """
    Returns the conjugate of a Sparse Matrix. The input and output matrix
    is given in DOK format.

    Examples
    ========

    >>> from sympy.matrices.sparsetools import conjugate
    >>> from sympy import SparseMatrix
    >>> from sympy import ZZ

    >>> a = SparseMatrix(2, 3, {(0, 0): 3, (1, 2): 12})

    >>> conjugate(a, ZZ)
    [3, 0,  0]
    [0, 0, 12]

    """
    smat = dok._smat
    for s in smat:
        smat[s] = smat[s].conjugate()
    return dok


def transpose(dok, K):
    """
    Returns the transpose of a Sparse Matrix. The input and output
    matrix is given in DOK format.

    Examples
    ========

    >>> from sympy.matrices.sparsetools import transpose
    >>> from sympy import SparseMatrix
    >>> from sympy import ZZ

    >>> a = SparseMatrix(2, 3, {(0, 0): ZZ(3), (1, 2): ZZ(12)})

    >>> transpose(a, ZZ)
    [3,  0]
    [0,  0]
    [0, 12]

    """
    smat_r = {}
    smat = dok._smat
    for k in smat:
        i, j = k
        smat_r[j, i] = smat[i, j]
    return SparseMatrix(dok.cols, dok.rows, smat_r)


def conjugate_transpose(dok, K):
    """
    Returns the conjugate-transpose of a Sparse Matrix. The input
    and output matrix is given in DOK format.

    Examples
    ========

    >>> from sympy.matrices.sparsetools import conjugate_transpose
    >>> from sympy import SparseMatrix
    >>> from sympy import ZZ

    >>> a = SparseMatrix(2, 3, {(0, 0): 3, (1, 2): 12})

    >>> conjugate_transpose(a, ZZ)
    [3,  0]
    [0,  0]
    [0, 12]

    """
    return conjugate(transpose(dok, K), K)


def csr_row(csr, i, K):
    """
    Returns the i'th row of a matrix. The input and output matrices
    are in CSR format.

    """
    a, ja, ia, shape = csr
    index = slice(ia[i], ia[i + 1])
    return [a[index], ja[index], [0, len(a[index])], [1, shape[1]]]
