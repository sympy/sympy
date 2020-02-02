'''Functions returning normal forms of matrices'''
from __future__ import division, print_function

from sympy.matrices.dense import diag, zeros


def smith_normal_form(m, domain = None):
    '''
    Return the Smith Normal Form of a matrix `m` over the ring `domain`.
    This will only work if the ring is a principal ideal domain.

    Examples
    ========

    >>> from sympy.polys.solvers import RawMatrix as Matrix
    >>> from sympy.polys.domains import ZZ
    >>> from sympy.matrices.normalforms import smith_normal_form
    >>> m = Matrix([[12, 6, 4], [3, 9, 6], [2, 16, 14]])
    >>> setattr(m, "ring", ZZ)
    >>> print(smith_normal_form(m))
    Matrix([[1, 0, 0], [0, 10, 0], [0, 0, -30]])

    '''
    invs = invariant_factors(m, domain=domain)
    smf = diag(*invs)
    n = len(invs)
    if m.rows > n:
        smf = smf.row_insert(m.rows, zeros(m.rows-n, m.cols))
    elif m.cols > n:
        smf = smf.col_insert(m.cols, zeros(m.rows, m.cols-n))
    return smf

def invariant_factors(m, domain = None):
    '''
    Return the tuple of abelian invariants for a matrix `m`
    (as in the Smith-Normal form)

    References
    ==========

    [1] https://en.wikipedia.org/wiki/Smith_normal_form#Algorithm
    [2] http://sierra.nmsu.edu/morandi/notes/SmithNormalForm.pdf

    '''
    if not domain:
        if not (hasattr(m, "ring") and m.ring.is_PID):
            raise ValueError(
                "The matrix entries must be over a principal ideal domain")
        else:
            domain = m.ring

    if len(m) == 0:
        return ()

    m = m[:, :]

    def add_rows(m, i, j, a, b, c, d):
        # replace m[i, :] by a*m[i, :] + b*m[j, :]
        # and m[j, :] by c*m[i, :] + d*m[j, :]
        for k in range(m.cols):
            e = m[i, k]
            m[i, k] = a*e + b*m[j, k]
            m[j, k] = c*e + d*m[j, k]

    def add_columns(m, i, j, a, b, c, d):
        # replace m[:, i] by a*m[:, i] + b*m[:, j]
        # and m[:, j] by c*m[:, i] + d*m[:, j]
        for k in range(m.rows):
            e = m[k, i]
            m[k, i] = a*e + b*m[k, j]
            m[k, j] = c*e + d*m[k, j]

    def clear_column(m):
        # make m[1:, 0] zero by row and column operations
        if m[0,0] == 0:
            return m
        pivot = m[0, 0]
        for j in range(1, m.rows):
            if m[j, 0] == 0:
                continue
            d, r = domain.div(m[j,0], pivot)
            if r == 0:
                add_rows(m, 0, j, 1, 0, -d, 1)
            else:
                a, b, g = domain.gcdex(pivot, m[j,0])
                d_0 = domain.div(m[j, 0], g)[0]
                d_j = domain.div(pivot, g)[0]
                add_rows(m, 0, j, a, b, d_0, -d_j)
                pivot = g
        return m

    def clear_row(m):
        # make m[0, 1:] zero by row and column operations
        if m[0] == 0:
            return m
        pivot = m[0, 0]
        for j in range(1, m.cols):
            if m[0, j] == 0:
                continue
            d, r = domain.div(m[0, j], pivot)
            if r == 0:
                add_columns(m, 0, j, 1, 0, -d, 1)
            else:
                a, b, g = domain.gcdex(pivot, m[0, j])
                d_0 = domain.div(m[0, j], g)[0]
                d_j = domain.div(pivot, g)[0]
                add_columns(m, 0, j, a, b, d_0, -d_j)
                pivot = g
        return m

    # permute the rows and columns until m[0,0] is non-zero if possible
    ind = [i for i in range(m.rows) if m[i,0] != 0]
    if ind and ind[0] != 0:
        m = m.permute_rows([[0, ind[0]]])
    else:
        ind = [j for j in range(m.cols) if m[0,j] != 0]
        if ind and ind[0] != 0:
            m = m.permute_cols([[0, ind[0]]])

    # make the first row and column except m[0,0] zero
    while (any([m[0,i] != 0 for i in range(1,m.cols)]) or
           any([m[i,0] != 0 for i in range(1,m.rows)])):
        m = clear_column(m)
        m = clear_row(m)

    if 1 in m.shape:
        invs = ()
    else:
        invs = invariant_factors(m[1:,1:], domain=domain)

    if m[0,0]:
        result = [m[0,0]]
        result.extend(invs)
        # in case m[0] doesn't divide the invariants of the rest of the matrix
        for i in range(len(result)-1):
            if result[i] and domain.div(result[i+1], result[i])[1] != 0:
                g = domain.gcd(result[i+1], result[i])
                result[i+1] = domain.div(result[i], g)[0]*result[i+1]
                result[i] = g
            else:
                break
    else:
        result = invs + (m[0,0],)
    return tuple(result)



def _nonzero(m):
    """
    returns a list of containing (i, j) representing the indices of
    non-zero elements in m.
    Help function for hermite_norlal_form
    """

    idx = []
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            if m[i, j] != 0:
                idx.append((i, j))
    return idx

def _nonzero_negative(m):
    """
    returns 0 if the first nonzero column j of A contains more than one nonzero
    entry, or contains only one nonzero entry and which is positive and  returns 1
    if the first nonzero column j of A contains only one nonzero entry, which
    is negative.
    """
    columns = list(zip(*_nonzero(m)))[1]
    nonzero_col = m[:, min(columns)]
    nums = [i for i in nonzero_col if i !=0]
    return len(nums) == 1 and nums[0] < 0


def hermite_normal_form(m, domain = None):
    """
    Returns the hermite_normal_form of integers matices.

    Examples
    ========

    >>> from sympy.polys.solvers import RawMatrix as Matrix
    >>> from sympy.polys.domains import ZZ
    >>> from sympy.matrices.normalforms import hermite_normal_form
    >>> m = Matrix([[2, 3, 6, 2], [5, 6, 1, 6], [8, 3, 1, 1]])
    >>> setattr(m, "ring", ZZ)
    >>> print(hermite_normal_form(m))
    (Matrix([[1, 0, 50, -11], [0, 3, 28, -2], [0, 0, 61, -13]]),
    Matrix([[ 9, -5, 1], [ 5, -2, 0], [11, -6, 1]]))

    """

    from sympy import eye, ones

    if domain != 'ZZ':
        raise ValueError("Hermite normal for is defined only for integer domain")
    row = m.shape[0]
    col = m.shape[1]

    A = m
    B = eye(row)
    L =zeros(row, row)
    D = ones(row + 1, 1)

    if _nonzero_negative(m):
        B[row - 1, row - 1] = -1
        A[row - 1, :] *= -1
    k = 1
    while(k < row):
        col1, col2 = _reduce_matrix(A, B, L, k, k - 1, D)
        u = (int(D[k - 1]) * int(D[k + 1]) +
                  int(L[k, k - 1]) * int(L[k, k - 1]))
        v = int(D[k]) * int(D[k])
        if col1 <= min(col2, col - 1) or (col1 == col and col2 == col and u < v):
            _swap_rows(k, A, B, L, D)
            if k > 1:
                k = k - 1
        else:
            for i in reversed(range(k - 1)):
                _reduce_matrix(A, B, L, k, i, D)
            k = k + 1
    return A[::-1, :], B[::-1, :]



def lnearint(a, b):
    if b < 0:
        a = -a
        b = -b
    y = a // b
    x = b * y
    z = a - x
    z = 2 * z
    if z > b:
        y = y + 1
    return y


def _reduce_matrix(A, B, L, k, i, D):
    """
    Helper function usd by hermite_normal_form function
    """

    non_zero_elems = list(zip(*_nonzero(A[i, :])))

    if len(non_zero_elems):
        col1 = non_zero_elems[1][0]
        if A[i, col1] < 0:
            L[i, :] = -L[i, :]
            L[:, i] = -L[:, i]
            A[i, :] *= -1
            B[i, :] *= -1
    else:
        col1 = A.shape[1]
    non_zero_elems = list(zip(*_nonzero(A[k, :])))
    if len(non_zero_elems):
        col2 = non_zero_elems[1][0]
    else:
        col2 = A.shape[1]
    if col1 < A.shape[1]:
        q = A[k, col1] // A[i, col1]
    else:
        if 2 * abs(L[k, i]) > D[i + 1]:
            q = lnearint(L[k, i], D[i + 1])
        else:
            return col1, col2
    A[k, :] -= q * A[i, :]
    B[k, :] -= q * B[i, :]
    L[k, i] -= q * D[i + 1]
    L[k, :i] -= q * L[i, :i]
    return col1, col2


def _swap_rows(k, A, B, L, D):
    """
    Helper function usd by hermite_normal_form function
    """

    # To avoid the interpretation of -1 as the last index of the matrix create
    # a reverse stop that ends past the negative of the length of the matrix
    reverse_stop = k - 2 if k > 1 else -(A.shape[0] + 1)
    # Swap rows of the matrices
    A[(k - 1):(k + 1), :] = A[k:reverse_stop:-1, :]
    B[(k - 1):(k + 1), :] = B[k:reverse_stop:-1, :]
    L[(k - 1):(k + 1), :(k - 1)] = L[k:reverse_stop:-1, :(k - 1)]
    t = (L[(k + 1):, k - 1] * D[k + 1] / D[k] -
         L[(k + 1):, k] * L[k, k - 1] / D[k])
    L[(k + 1):, k - 1] = (L[(k + 1):, k - 1] * L[k, k - 1] +
                          L[(k + 1):, k] * D[k - 1]) / D[k]
    L[(k + 1):, k] = t
    t = int(D[k - 1]) * int(D[k + 1]) + int(L[k, k - 1]) * int(L[k, k - 1])
    D[k] = t / D[k]
