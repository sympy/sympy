from __future__ import print_function, division

from sympy import SparseMatrix, Matrix


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


def _doktodia(dok):
    """Return a tuple of the lower and upper bandwidths and converts the
       DOK to diagonal ordered form:

       ab[u + i - j, j] == a[i, j], where 'u' is the upper bandwidth of 'a'.
    """
    offsets = sorted(set(map(lambda x: x[1]-x[0], dok.CL)))[::-1]
    offmap = {i:o for i,o in enumerate(offsets)}
    def todia(i, j):
        im = j-offmap[i]
        if 0 <= im < dok.shape[0]:
            return dok[im,j]
        return 0
    return ((abs(offsets[-1]), offsets[0]),
            Matrix(len(offsets), dok.shape[0], todia))

def _banded_LUdecomposition(A, l_and_u=None, overwrite=False):
    """
    Returns the LU decomposition of a matrix A in diagonal ordered form.

    Note that a single matrix utilizing diagonal ordered form will be
    returned with LU[0:u, :] corresponding to U and LU[u:, :] corresponds
    to L.

    Parameters
    ==========

    A : Either a `SparseMatrix` or a `Matrix` in diagonal ordered form.
    l_and_u : Tuple of lower and upper bandwidths.  Must be specified if
    A is in diagonal ordered form.
    overwrite : If True, overwrite input matrix A.  Applies if A is in diagonal
    form.

    References
    ==========

    Matrix Computations, 4th Ed.
    Gene H. Golub, Charles F. Van Loan  (2013)
    """
    if isinstance(A, SparseMatrix):
        (l, u), A = _doktodia(A)
    else:
        if l_and_u is None:
            raise TypeError("Expected tuple of lower and upper bandwidths.")
        l, u = l_and_u
        if not overwrite:
            A = A.copy()

    n = A.shape[-1]
    for k in range(n-1):
        A[u+1:u+l+1, k] /= A[u, k]

        i, ii = k+1, min(k+l, n)+1
        for j in range(k+1, min(k+u+1, n)):
            A[u+i-j:u+ii-j, j] -= A[u+i-k:u+ii-k, k]*A[u+k-j,j]

    return (l, u), m

def banded_LUsolve(A, rhs, l_and_u=None, have_LU=False, overwrite=False):
    """
    Solves the linear system A * x = rhs via LU decomposition where A is banded

    Parameters
    ==========

    A : Either a `SparseMatrix` or a `Matrix` in diagonal ordered form.  If the
    same system will be solved multiple times, A should be passed in as an LU
    factorized matrix in diagonal form (i.e. the result of
    _banded_LUdecomposition).  This should significantly increase performance.
    rhs : solution vector/matrix.
    l_and_u : Tuple of lower and upper bandwidths.  Must be specified if
    A is in diagonal ordered form.
    have_LU : Set this to True if A is already an LU factorized matrix in
    diagonal form.  If set to False (default), The LU decomposition
    will be computed.
    overwrite : If True, overwrite solution matrix rhs.

    References
    ==========

    Matrix Computations, 4th Ed.
    Gene H. Golub, Charles F. Van Loan  (2013)
    """
    if isinstance(A, SparseMatrix):
        (l, u), A = _banded_LUdecomposition(A)
    else:
        if l_and_u is None:
            raise TypeError("Expected tuple of lower, upper bandwidths.")
        if not have_LU:
            (l, u), A = _banded_LUdecomposition(A, (l, u), overwrite=overwrite)
    rhs = rhs.copy()
    l, u = l_and_u
    n = lu.shape[-1]
    # forward subs
    for j in range(n):
        i, ii = j+1, min(j+l+1, n)
        rhs[i:ii, :] -= lu[u+i-j:u+ii-j, j]*rhs[j, :]
    # backward subs
    for j in range(n-1, -1, -1):
        rhs[j, :] /= lu[u, j]
        i, ii = max(0, j-u), j
        rhs[i:ii, :] -= lu[u+i-j:u+ii-j, j]*rhs[j, :]
    return rhs
