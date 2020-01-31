from __future__ import division, print_function

import copy

from sympy.core.function import expand_mul
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.simplify.simplify import dotprodsimp as _dotprodsimp

from .common import NonSquareMatrixError, NonPositiveDefiniteMatrixError


def _liupc(M):
    """Liu's algorithm, for pre-determination of the Elimination Tree of
    the given matrix, used in row-based symbolic Cholesky factorization.

    Examples
    ========

    >>> from sympy.matrices import SparseMatrix
    >>> S = SparseMatrix([
    ... [1, 0, 3, 2],
    ... [0, 0, 1, 0],
    ... [4, 0, 0, 5],
    ... [0, 6, 7, 0]])
    >>> S.liupc()
    ([[0], [], [0], [1, 2]], [4, 3, 4, 4])

    References
    ==========

    Symbolic Sparse Cholesky Factorization using Elimination Trees,
    Jeroen Van Grondelle (1999)
    http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.39.7582
    """
    # Algorithm 2.4, p 17 of reference

    # get the indices of the elements that are non-zero on or below diag
    R = [[] for r in range(M.rows)]

    for r, c, _ in M.row_list():
        if c <= r:
            R[r].append(c)

    inf     = len(R)  # nothing will be this large
    parent  = [inf]*M.rows
    virtual = [inf]*M.rows

    for r in range(M.rows):
        for c in R[r][:-1]:
            while virtual[c] < r:
                t          = virtual[c]
                virtual[c] = r
                c          = t

            if virtual[c] == inf:
                parent[c] = virtual[c] = r

    return R, parent

def _row_structure_symbolic_cholesky(M):
    """Symbolic cholesky factorization, for pre-determination of the
    non-zero structure of the Cholesky factororization.

    Examples
    ========

    >>> from sympy.matrices import SparseMatrix
    >>> S = SparseMatrix([
    ... [1, 0, 3, 2],
    ... [0, 0, 1, 0],
    ... [4, 0, 0, 5],
    ... [0, 6, 7, 0]])
    >>> S.row_structure_symbolic_cholesky()
    [[0], [], [0], [1, 2]]

    References
    ==========

    Symbolic Sparse Cholesky Factorization using Elimination Trees,
    Jeroen Van Grondelle (1999)
    http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.39.7582
    """

    R, parent = M.liupc()
    inf       = len(R)  # this acts as infinity
    Lrow      = copy.deepcopy(R)

    for k in range(M.rows):
        for j in R[k]:
            while j != inf and j != k:
                Lrow[k].append(j)
                j = parent[j]

        Lrow[k] = list(sorted(set(Lrow[k])))

    return Lrow


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

    sympy.matrices.dense.DenseMatrix.LDLdecomposition
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

    The matrix can have complex entries:

    >>> from sympy import I
    >>> A = SparseMatrix(((9, 3*I), (-3*I, 5)))
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

    >>> A = SparseMatrix([[1, 2], [2, 1]])
    >>> L = A.cholesky(hermitian=False)
    >>> L
    Matrix([
    [1,         0],
    [2, sqrt(3)*I]])
    >>> L*L.T == A
    True

    See Also
    ========

    sympy.matrices.sparse.SparseMatrix.LDLdecomposition
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

    dps       = _dotprodsimp if dotprodsimp else expand_mul
    Crowstruc = M.row_structure_symbolic_cholesky()
    C         = MutableDenseMatrix.zeros(M.rows)

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
                                    if hermitian:
                                        summ += C[i, p1]*C[j, p1].conjugate()
                                    else:
                                        summ += C[i, p1]*C[j, p1]
                            else:
                                break
                        else:
                            break

                C[i, j] = dps((C[i, j] - summ) / C[j, j])

            else: # i == j
                C[j, j] = M[j, j]
                summ    = 0

                for k in Crowstruc[j]:
                    if k < j:
                        if hermitian:
                            summ += C[j, k]*C[j, k].conjugate()
                        else:
                            summ += C[j, k]**2
                    else:
                        break

                Cjj2 = dps(C[j, j] - summ)

                if hermitian and Cjj2.is_positive is False:
                    raise NonPositiveDefiniteMatrixError(
                        "Matrix must be positive-definite")

                C[j, j] = sqrt(Cjj2)

    return M._new(C)


def _LDLdecomposition(M, hermitian=True, dotprodsimp=None):
    """Returns the LDL Decomposition (L, D) of matrix A,
    such that L * D * L.H == A if hermitian flag is True, or
    L * D * L.T == A if hermitian is False.
    This method eliminates the use of square root.
    Further this ensures that all the diagonal entries of L are 1.
    A must be a Hermitian positive-definite matrix if hermitian is True,
    or a symmetric matrix otherwise.

    Parameters
    ==========

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    Examples
    ========

    >>> from sympy.matrices import Matrix, eye
    >>> A = Matrix(((25, 15, -5), (15, 18, 0), (-5, 0, 11)))
    >>> L, D = A.LDLdecomposition()
    >>> L
    Matrix([
    [   1,   0, 0],
    [ 3/5,   1, 0],
    [-1/5, 1/3, 1]])
    >>> D
    Matrix([
    [25, 0, 0],
    [ 0, 9, 0],
    [ 0, 0, 9]])
    >>> L * D * L.T * A.inv() == eye(A.rows)
    True

    The matrix can have complex entries:

    >>> from sympy import I
    >>> A = Matrix(((9, 3*I), (-3*I, 5)))
    >>> L, D = A.LDLdecomposition()
    >>> L
    Matrix([
    [   1, 0],
    [-I/3, 1]])
    >>> D
    Matrix([
    [9, 0],
    [0, 4]])
    >>> L*D*L.H == A
    True

    See Also
    ========

    sympy.matrices.dense.DenseMatrix.cholesky
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
    D   = MutableDenseMatrix.zeros(M.rows, M.rows)
    L   = MutableDenseMatrix.eye(M.rows)

    if hermitian:
        for i in range(M.rows):
            for j in range(i):
                L[i, j] = dps((1 / D[j, j])*(M[i, j] - sum(
                    L[i, k]*L[j, k].conjugate()*D[k, k] for k in range(j))))

            D[i, i] = dps(M[i, i] -
                sum(L[i, k]*L[i, k].conjugate()*D[k, k] for k in range(i)))

            if D[i, i].is_positive is False:
                raise NonPositiveDefiniteMatrixError(
                    "Matrix must be positive-definite")

    else:
        for i in range(M.rows):
            for j in range(i):
                L[i, j] = dps((1 / D[j, j])*(M[i, j] - sum(
                    L[i, k]*L[j, k]*D[k, k] for k in range(j))))

            D[i, i] = dps(M[i, i] - sum(L[i, k]**2*D[k, k] for k in range(i)))

    return M._new(L), M._new(D)

def _LDLdecomposition_sparse(M, hermitian=True, dotprodsimp=None):
    """
    Returns the LDL Decomposition (matrices ``L`` and ``D``) of matrix
    ``A``, such that ``L * D * L.T == A``. ``A`` must be a square,
    symmetric, positive-definite and non-singular.

    This method eliminates the use of square root and ensures that all
    the diagonal entries of L are 1.

    Examples
    ========

    >>> from sympy.matrices import SparseMatrix
    >>> A = SparseMatrix(((25, 15, -5), (15, 18, 0), (-5, 0, 11)))
    >>> L, D = A.LDLdecomposition()
    >>> L
    Matrix([
    [   1,   0, 0],
    [ 3/5,   1, 0],
    [-1/5, 1/3, 1]])
    >>> D
    Matrix([
    [25, 0, 0],
    [ 0, 9, 0],
    [ 0, 0, 9]])
    >>> L * D * L.T == A
    True

    """

    from .dense import MutableDenseMatrix

    if not M.is_square:
        raise NonSquareMatrixError("Matrix must be square.")
    if hermitian and not M.is_hermitian:
        raise ValueError("Matrix must be Hermitian.")
    if not hermitian and not M.is_symmetric():
        raise ValueError("Matrix must be symmetric.")

    dps       = _dotprodsimp if dotprodsimp else expand_mul
    Lrowstruc = M.row_structure_symbolic_cholesky()
    L         = MutableDenseMatrix.eye(M.rows)
    D         = MutableDenseMatrix.zeros(M.rows, M.cols)

    for i in range(len(Lrowstruc)):
        for j in Lrowstruc[i]:
            if i != j:
                L[i, j] = M[i, j]
                summ    = 0

                for p1 in Lrowstruc[i]:
                    if p1 < j:
                        for p2 in Lrowstruc[j]:
                            if p2 < j:
                                if p1 == p2:
                                    if hermitian:
                                        summ += L[i, p1]*L[j, p1].conjugate()*D[p1, p1]
                                    else:
                                        summ += L[i, p1]*L[j, p1]*D[p1, p1]
                            else:
                                break
                    else:
                        break

                L[i, j] = dps((L[i, j] - summ) / D[j, j])

            else: # i == j
                D[i, i] = M[i, i]
                summ    = 0

                for k in Lrowstruc[i]:
                    if k < i:
                        if hermitian:
                            summ += L[i, k]*L[i, k].conjugate()*D[k, k]
                        else:
                            summ += L[i, k]**2*D[k, k]
                    else:
                        break

                D[i, i] = dps(D[i, i] - summ)

                if hermitian and D[i, i].is_positive is False:
                    raise NonPositiveDefiniteMatrixError(
                        "Matrix must be positive-definite")

    return M._new(L), M._new(D)
