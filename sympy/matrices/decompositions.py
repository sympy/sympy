from __future__ import division, print_function

from sympy.core.function import expand_mul
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.simplify.simplify import dotprodsimp as _dotprodsimp

from .common import NonSquareMatrixError, NonPositiveDefiniteMatrixError


def _cholesky(M, hermitian=True, dotprodsimp=None):
    """Returns the Cholesky-type decomposition L of a matrix A
    such that L * L.H == A if hermitian flag is True,
    or L * L.T == A if hermitian is False.

    A must be a Hermitian positive-definite matrix if hermitian is True,
    or a symmetric matrix if it is False.

    Parameters
    ==========

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    Examples
    ========

    >>> from sympy.matrices import Matrix
    >>> A = Matrix(((25, 15, -5), (15, 18, 0), (-5, 0, 11)))
    >>> A.cholesky()
    Matrix([
    [ 5, 0, 0],
    [ 3, 3, 0],
    [-1, 1, 3]])
    >>> A.cholesky() * A.cholesky().T
    Matrix([
    [25, 15, -5],
    [15, 18,  0],
    [-5,  0, 11]])

    The matrix can have complex entries:

    >>> from sympy import I
    >>> A = Matrix(((9, 3*I), (-3*I, 5)))
    >>> A.cholesky()
    Matrix([
    [ 3, 0],
    [-I, 2]])
    >>> A.cholesky() * A.cholesky().H
    Matrix([
    [   9, 3*I],
    [-3*I,   5]])

    Non-hermitian Cholesky-type decomposition may be useful when the
    matrix is not positive-definite.

    >>> A = Matrix([[1, 2], [2, 1]])
    >>> L = A.cholesky(hermitian=False)
    >>> L
    Matrix([
    [1,         0],
    [2, sqrt(3)*I]])
    >>> L*L.T == A
    True

    See Also
    ========

    sympy.matrices.matrices.MatrixBase.LDLdecomposition
    LUdecomposition
    QRdecomposition
    """

    from .dense import MutableDenseMatrix

    if not M.is_square:
        raise NonSquareMatrixError("Matrix must be square.")
    if hermitian and not M.is_hermitian:
        raise ValueError("Matrix must be Hermitian.")
    if not hermitian and not M.is_symmetric():
        raise ValueError("Matrix must be symmetric.")

    dps = _dotprodsimp if dotprodsimp else expand_mul
    L   = MutableDenseMatrix.zeros(M.rows, M.rows)

    if hermitian:
        for i in range(M.rows):
            for j in range(i):
                L[i, j] = dps((1 / L[j, j])*(M[i, j] -
                    sum(L[i, k]*L[j, k].conjugate() for k in range(j))))

            Lii2 = dps(M[i, i] -
                sum(L[i, k]*L[i, k].conjugate() for k in range(i)))

            if Lii2.is_positive is False:
                raise NonPositiveDefiniteMatrixError(
                    "Matrix must be positive-definite")

            L[i, i] = sqrt(Lii2)

    else:
        for i in range(M.rows):
            for j in range(i):
                L[i, j] = dps((1 / L[j, j])*(M[i, j] -
                    sum(L[i, k]*L[j, k] for k in range(j))))

            L[i, i] = sqrt(dps(M[i, i] -
                sum(L[i, k]**2 for k in range(i))))

    return M._new(L)


def _cholesky_sparse(M, hermitian=True, dotprodsimp=None):
    """
    Returns the Cholesky decomposition L of a matrix A
    such that L * L.T = A

    A must be a square, symmetric, positive-definite
    and non-singular matrix

    Examples
    ========

    >>> from sympy.matrices import SparseMatrix
    >>> A = SparseMatrix(((25,15,-5),(15,18,0),(-5,0,11)))
    >>> A.cholesky()
    Matrix([
    [ 5, 0, 0],
    [ 3, 3, 0],
    [-1, 1, 3]])
    >>> A.cholesky() * A.cholesky().T == A
    True
    """

    from sympy.core.numbers import nan, oo

    def cholesky_sparse(M, dotprodsimp=None):
        """Algorithm for numeric Cholesky factorization of a sparse matrix."""

        dps       = _dotprodsimp if dotprodsimp else expand_mul
        Crowstruc = M.row_structure_symbolic_cholesky()
        C         = M.zeros(M.rows)

        for i in range(len(Crowstruc)):
            for j in Crowstruc[i]:
                if i != j:
                    C[i, j] = M[i, j]
                    summ    = 0

                    for p1 in Crowstruc[i]:
                        if p1 < j:
                            for p2 in Crowstruc[j]:
                                if p2 < j:
                                    if p1 == p2:
                                        summ += C[i, p1]*C[j, p1]
                                else:
                                    break
                            else:
                                break

                    C[i, j] = dps((C[i, j] - summ) / C[j, j])

                else:
                    C[j, j] = M[j, j]
                    summ    = 0

                    for k in Crowstruc[j]:
                        if k < j:
                            summ += C[j, k]**2
                        else:
                            break

                    C[j, j] = dps(C[j, j] - summ)
                    C[j, j] = sqrt(C[j, j])

        return C

    if not M.is_square:
        raise NonSquareMatrixError("Matrix must be square.")
    if hermitian and not M.is_hermitian:
        raise ValueError("Matrix must be Hermitian.")
    if not hermitian and not M.is_symmetric():
        raise ValueError("Matrix must be symmetric.")

    ret = cholesky_sparse(M.as_mutable(), dotprodsimp=dotprodsimp)

    if ret.has(nan) or ret.has(oo):
        raise ValueError('Cholesky decomposition applies only to '
            'positive-definite matrices')

    return M._new(ret)
