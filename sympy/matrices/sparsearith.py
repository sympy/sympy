"""
Fundamental Arithmetic of Sparse Matrices. The Data structure for Sparse
Matrix is Dictionary of Keys(DOK) and Compressed Sparse Row(CSR).

"""
from sympy import SparseMatrix
from collections import defaultdict
from sympy.core.compatibility import xrange


def _doktocsr(dok):
    """
    Converts a Sparse Matrix in Dictionary of Keys (DOK) format
    to Compressed Sparse Row (CSR) format.

    CONSTITUENTS
    ==========

    A : contains non-zero elements sorted by key (row, column).
    JA : JA[i] is the column corresponding to the element A[i].
    IA : IA[i] contains the index in A for the first non-zero element
        of row[i]. Thus IA[i+1] - IA[i] gives number of non-zero
        elements row[i]. The length of IA is always 1 more than the
        number of rows in the matrix.

    Examples
    =========
    >>> from sympy import SparseMatrix, ZZ
    >>> from sympy.matrices.sparsearith import _doktocsr
    >>> test_matrix = SparseMatrix([[ZZ(1), ZZ(2), ZZ(0), ZZ(0)], [ZZ(0), ZZ(3), ZZ(9), ZZ(0)], [ZZ(0), ZZ(1), ZZ(4), ZZ(0)]])
    >>> _doktocsr(test_matrix)
    [[1, 2, 3, 9, 1, 4], [0, 1, 1, 2, 1, 2], [0, 2, 4, 6], [3, 4]]

    See Also
    ========

    _csrtodok

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
    """
    Converts a Sparse Matrix in Compressed Sparse Row (CSR) format
    to Dictionary of Keys (DOK) format.

    Examples
    =========
    >>> from sympy import SparseMatrix, ZZ
    >>> from sympy.matrices.sparsearith import _csrtodok
    >>> test_matrix_csr = [[ZZ(5), ZZ(7), ZZ(5)], [2, 1, 3], [0, 1, 1, 3], [3, 4]]
    >>> _csrtodok(test_matrix_csr)
    [[0, 0, 5, 0], [0, 0, 0, 0], [0, 7, 0, 5]]

    See Also
    ========

    _doktocsr

    """
    smat = {}
    A, JA, IA, shape = csr
    for i in range(len(IA) - 1):
        indices = slice(IA[i], IA[i + 1])
        for l, m in zip(A[indices], JA[indices]):
            smat[i, m] = l
    return SparseMatrix(*(shape + [smat]))


def add(dok1, dok2, K):
    """
    Adds two matrices which are in Dictionary of Keys (DOK)
    format. The result returned is in DOK representation.

    Examples
    =========
    >>> from sympy.matrices.sparsearith import add
    >>> from sympy import SparseMatrix, ZZ
    >>> a = SparseMatrix(3, 4, {(0, 0): ZZ(1), (0, 1): ZZ(2), (2, 2): ZZ(3)})
    >>> b = SparseMatrix(3, 4, {(0, 0): ZZ(1), (0, 1): ZZ(2), (1, 1): ZZ(3), (1, 2): ZZ(9), (2, 1): ZZ(1), (2, 2): ZZ(4)})
    >>> add(a, b, ZZ)
    [2, 4, 0, 0]
    [0, 3, 9, 0]
    [0, 1, 7, 0]

    See Also
    ========

    sub

    """
    dok_result = SparseMatrix(dok1.rows, dok2.cols, {})
    non_zero_positions1, non_zero_positions2 = dok1._smat.keys(), dok2._smat.keys()
    common_positions = list(set(non_zero_positions1).intersection(set(non_zero_positions2)))
    all_positions = iter(set(non_zero_positions1).union(set(non_zero_positions2)))
    for position in all_positions:
        if position in common_positions:
            dok_result._smat[position] = dok1._smat[position] + dok2._smat[position]
        else:
            dok_result._smat[position] = dok1._smat[position] if position in non_zero_positions1 else dok2._smat[position]
    return dok_result


def sub(dok1, dok2, K):
    """
    Subtracts two matrices which are in Dictionary of Keys (DOK)
    format. The result returned is in DOK format.

    Examples
    =========
    >>> from sympy.matrices.sparsearith import sub
    >>> from sympy import SparseMatrix, ZZ
    >>> a = SparseMatrix(3, 4, {(0, 0): ZZ(1), (0, 1): ZZ(2), (2, 2): ZZ(3)})
    >>> b = SparseMatrix(3, 4, {(0, 0): ZZ(1), (0, 1): ZZ(2), (1, 1): ZZ(3), (1, 2): ZZ(9), (2, 1): ZZ(1), (2, 2): ZZ(4)})
    >>> sub(a, b, ZZ)
    [0, 0, 0, 0]
    [0, -3, -9, 0]
    [0, -1, -1, 0]

    See Also
    ========

    add
    neg

    """
    return add(dok1, neg(dok2, K), K)


def mulspvec(csr, vector, K):
    """
    Multiplies a Sparse Matrix and a vector. The Sparse Matrix is
    represented in Compressed Sparse Row (CSR) format. The vector
    maybe dense or sparse. The returned vector is a Dense Matrix.

    Examples
    =========
    >>> from sympy.matrices.sparsearith import mulspvec
    >>> from sympy import SparseMatrix, ZZ
    >>> csr = [[ZZ(1), ZZ(2), ZZ(3)], [0, 1, 2], [0, 2, 2, 3], [3, 4]]
    >>> vector = SparseMatrix([ZZ(1), ZZ(2), ZZ(3), ZZ(4)])
    >>> mulspvec(csr, vector, ZZ)
    [5, 0, 9]

    See Also
    ========

    mulrowvec
    get_row

    """
    a, ja, ia, shape = csr
    number_of_rows = len(ia) - 1
    result_vector = []
    for i in xrange(number_of_rows):
        current_row_indices = slice(ia[i], ia[i + 1])
        current_row = get_row(csr, current_row_indices)
        result_vector.append(mulrowvec(current_row, vector, K))
    return SparseMatrix(result_vector)
    

def get_row(csr, indices):
    """
    Returns a row of a Sparse Matrix in compressed Sparse row (CSR)
    format given a Sparse Matrix in CSR format and indices.

    Examples
    =========
    >>> from sympy.matrices.sparsearith import get_row
    >>> from sympy import ZZ
    >>> csr = [[ZZ(1), ZZ(2), ZZ(3)], [0, 1, 2], [0, 2, 2, 3], [3, 4]]
    >>> get_row(csr, slice(1, 2))
    [[2], [1], [0, 1], [1, 4]]

    See Also
    ========

    mulspvec

    """
    a, ja, ia, shape = csr
    a_row, ja_row = a[indices], ja[indices]
    ia_row = [0, len(a_row)]
    return [a_row, ja_row, ia_row, [1, shape[1]]]


def mulrowvec(row, vector, K):
    """
    Multiplies the row of a Sparse Matrix represented in Compressed
    Sparse Row (CSR) format and a vector.

    Examples
    =========
    >>> from sympy.matrices.sparsearith import mulrowvec
    >>> from sympy import SparseMatrix
    >>> row = [[ZZ(2)], [1], [0, 1], [1, 4]]
    >>> vector = SparseMatrix([ZZ(1), ZZ(2), ZZ(3), ZZ(4)])
    >>> mulrowvec(row, vector, ZZ)
    4

    See Also
    ========

    mulspvec

    """
    result = 0
    a_row, ja_row, ia_row, shape = row
    for index, column_index in enumerate(ja_row):
        result += vector[column_index] * a_row[index]
    return result


def neg(dok, K):
    """
    Negates each element of a Sparse Matrix represeneted in
    Dictionary of Keys (DOK) format.

    Examples
    ========
    >>> from sympy import SparseMatrix
    >>> from sympy.matrices.sparsearith import neg
    >>> matrix = SparseMatrix(3, 4, {(0, 0): ZZ(1), (0, 1): ZZ(2), (2, 2): ZZ(3)})
    >>> neg(matrix, ZZ)
    [-1, -2, 0, 0]
    [0, 0, 0, 0]
    [0, 0, -3, 0]

    See Also
    =========

    sub

    """
    result = SparseMatrix(dok.rows, dok.cols, {})
    for key in dok._smat.keys():
        result._smat[key] = -(dok._smat[key])
    return result
