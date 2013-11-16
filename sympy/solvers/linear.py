from __future__ import print_function, division

from sympy.core import Dummy
from sympy.utilities.iterables import numbered_symbols
from sympy.matrices import Matrix, zeros


def solve_general_linear(M, v, dummygen=None):
    r"""This function solves the general linear equation

    .. math:: A x = b

    for any matrix :math:`A` having :math:`m` rows and :math:`n`
    columns and arbitrary rank :math:`r`. Returned are a matrix
    of shape :math:`(n,1)` containing the solutions and a second
    matrix of shape :math:`(m-r,1)` with all free parameters.

    Example
    =======

    >>> from sympy import Matrix
    >>> from sympy.solvers.linear import solve_general_linear
    >>> M = Matrix([[1,2,3],[4,5,6],[7,8,9]])

    >>> M.det()
    0

    >>> v = Matrix([3,6,9])

    >>> sol, params = solve_general_linear(M, v)
    >>> sol
    Matrix([
    [   _tau0 - 1],
    [-2*_tau0 + 2],
    [       _tau0]])
    >>> params
    Matrix([[_tau0]])

    It is possible to provide in the third argument a generator
    from which the free variables if any will be taken.

    >>> from sympy import Symbol
    >>> x = Symbol("x")
    >>> y = Symbol("y")
    >>> z = Symbol("z")

    >>> M = Matrix([[1,2,3],[2,4,6],[3,6,9]])
    >>> v = Matrix([0,0,0])

    >>> sol, params = solve_general_linear(M, v, (s for s in (x,y,z)))
    >>> sol
    Matrix([
    [-2*x - 3*y],
    [         x],
    [         y]])
    >>> params
    Matrix([
    [x],
    [y]])

    The solver will raise a ValueError exception in case the
    linear system has no solution.
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
    if dummygen is None:
        dummygen = numbered_symbols("tau", Dummy)
    tau = Matrix([[ dummygen.next() for k in xrange(C - rank) ]]).T

    # Full parametric solution
    V = U[:rank,rank:]
    vt = v[:rank,0]
    xsigma = tau.vstack(vt - V*tau, tau)

    # Undo permutation
    sol = zeros(C, 1)
    for k, xksigma in enumerate(xsigma):
        sol[permutation[k],0] = xksigma

    return sol, tau
