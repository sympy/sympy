"""
Fundamental Operations of Sparse Matrices. The Data structure for Sparse
Matrix is Dictionary of Keys(DOK) and Compressed Sparse Row(CSR).

"""
from __future__ import print_function, division
from sympy import SparseMatrix
from sympy.core.compatibility import xrange


def trace(dok, K):
    """
    Returns the trace of a Sparse Matrix. The input and output matrices
    are in Dictionary of Keys (DOK) format.

    Examples
    ========

    >>> from sympy.matrices.sparsetools import trace
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


def conjugate(dok, K):
    """
    Note: Not working for Complex Entries
    Returns the conjugate of a Sparse Matrix. The input and output matrices
    are in Dicionary of Keys (DOK) format.

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
    Returns the transpose of a Sparse Matrix.The input and output matrices
    are in Dicionary of Keys (DOK) format.

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
    Returns the conjugate-transpose of a Sparse Matrix. The input and output
    matrices are in Dicionary of Keys (DOK) format.

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


def _binsearch(element, array, beg, end):
    """
    An implementation of binary search. It is expected to be
    used as follows for an element i and list v--
    _binsearch(i, v, 0, len(v) - 1).

    """
    if beg > end:
        return None
    mid = (beg + end)//2
    if element < array[mid]:
        return _binsearch(element, array, beg, mid - 1)
    elif element > array[mid]:
        return _binsearch(element, array, mid + 1, end)
    else:
        return mid
