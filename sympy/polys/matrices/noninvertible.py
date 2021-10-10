"""

Tools for working with non-invertible (or not-necessarily invertible)
DomainMatrices.

"""


from .exceptions import DMRankError


def invertible_supplement(M):
    """
    Append carefully chosen columns to M to make a square invertible matrix.

    Explanation
    ===========

    Given an n x r matrix M of rank r (so r <= n), this function computes an
    invertible n x n matrix B such that the first r columns of B equal M.

    This operation can be interpreted as a way of extending a basis for a
    subspace, to give a basis for the whole space.

    To be precise, suppose you have an n-dimensional vector space V, with basis
    vectors v1, v2, ..., vn, and an r-dimensional subspace W of V, spanned by a
    basis w1, w2, ..., wr, where the w_j are given as linear combinations of
    the v_i. If the columns of M represent the w_j as such lin. combs., then
    the columns of the matrix B computed by this function give a new basis
    u1, u2, ..., un for V, again relative to the {v_i} basis, and such that
    u_j == w_j for 1 <= j <= r.

    Parameters
    ==========

    M: the DomainMatrix to be supplemented

    Returns
    =======

    DomainMatrix
        the supplemented Matrix B

    Raises
    ======

    DMRankError if M was not of maximal rank

    References
    ==========

    Cohen, H. *A Course in Computational Algebraic Number Theory*
    (See Sec. 2.3.2.)

    """
    n, r = M.shape
    # Let In be the n x n identity matrix.
    # Form the augmented matrix [M | In] and compute RREF.
    Maug = M.hstack(M.eye(n, M.domain))
    R, pivots = Maug.rref()
    if pivots[:r] != tuple(range(r)):
        raise DMRankError('M was not of maximal rank')
    # Let J be the n x r matrix equal to the first r columns of In.
    # Since M is of rank r, RREF reduces [M | In] to [J | A], where A is the product of
    # elementary matrices Ei corresp. to the row ops performed by RREF. Since the Ei are
    # invertible, so is A. Let B = A^(-1).
    A = R[:, r:]
    B = A.inv()
    # Then B is the desired matrix. It is invertible, since B^(-1) == A.
    # And A * [M | In] == [J | A]
    #  => A * M == J
    #  => M == B * J == the first r columns of B.
    return B
