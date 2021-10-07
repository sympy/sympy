"""

Tools for working with non-invertible (or not-necessarily invertible)
DomainMatrices.

"""


def invertible_supplement(M):
    """
    Given an n x r matrix M of rank r (so r <= n), find an invertible n x n
    matrix B such that the first r columns of B equal M.

    Parameters
    ==========

    M: the DomainMatrix to be supplemented

    Returns
    =======

    DomainMatrix
        the supplemented Matrix B

    Raises
    ======

    ValueError if M was not of rank r

    """
    n, r = M.shape
    M = M.hstack(M.eye(n, M.domain))

    R, pivots = M.rref()
    if pivots[:r] != tuple(range(r)):
        raise ValueError('M was not of maximal rank')

    A = R[:, r:]
    B = A.inv()
    return B
