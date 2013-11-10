from __future__ import print_function, division

from sympy.core import Dummy
from sympy.matrices import Matrix, zeros


def solve_general_linear(M, v):
    r"""This function solves the general linear equation

    .. math:: A x = b

    for any matrix :math:`A` having :math:`m` rows and :math:`n`
    columns and arbitrary rank :math:`r`. Returned are a matrix
    of shape :math:`(n,1)` containing the solutions and a second
    matrix of shape :math:`(m-r,1)` with all free parameters.
    """
    R, C = M.shape
    U = M.hstack(M.copy(), v.copy())

    # Gauss Jordan elimination to produce a row echelon form
    U, pivots = U.rref()
    U, v = U[:,:-1], U[:,-1]
    pivots = filter(lambda p: p < C, pivots)
    rank = len(pivots)

    # Bring to block form
    permutation = Matrix(range(C)).T
    U = U.vstack(U, permutation)

    for i, c in enumerate(pivots):
        U.col_swap(i,c)

    U, permutation = U[:-1,:], U[-1,:]

    # Check if there are solutions
    vzero = v[rank:,0]
    if not vzero.is_zero:
        raise ValueError("Linear system has no solution")

    # Free parameters
    # The double List and T operation are a hack to get the shapes correct
    tau = Matrix([[ Dummy("t_" + str(k)) for k in xrange(C - rank) ]]).T

    # Full parametric solution
    V = U[:rank,rank:]
    vt = v[:rank,0]
    xsigma = tau.vstack(vt - V*tau, tau)

    # Undo permutation
    sol = zeros(C, 1)
    for k, xksigma in enumerate(xsigma):
        sol[permutation[k],0] = xksigma

    return sol, tau
