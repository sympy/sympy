from __future__ import division, print_function

from sympy.simplify.simplify import dotprodsimp as _dotprodsimp

from .common import ShapeError, NonSquareMatrixError


def _lower_triangular_solve(M, rhs, dotprodsimp=None):
    """Solves ``Ax = B``, where A is a lower triangular matrix.

    Parameters
    ==========

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    See Also
    ========

    upper_triangular_solve
    gauss_jordan_solve
    cholesky_solve
    diagonal_solve
    LDLsolve
    LUsolve
    QRsolve
    pinv_solve
    """

    from .dense import MutableDenseMatrix

    if not M.is_square:
        raise NonSquareMatrixError("Matrix must be square.")
    if rhs.rows != M.rows:
        raise ShapeError("Matrices size mismatch.")
    if not M.is_lower:
        raise ValueError("Matrix must be lower triangular.")

    dps = _dotprodsimp if dotprodsimp else lambda x: x
    X   = MutableDenseMatrix.zeros(M.rows, rhs.cols)

    for j in range(rhs.cols):
        for i in range(M.rows):
            if M[i, i] == 0:
                raise TypeError("Matrix must be non-singular.")

            X[i, j] = dps((rhs[i, j] - sum(M[i, k]*X[k, j]
                                        for k in range(i))) / M[i, i])

    return M._new(X)


def _lower_triangular_solve_sparse(M, rhs, dotprodsimp=None):
    """Solves ``Ax = B``, where A is a lower triangular matrix.

    Parameters
    ==========

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    See Also
    ========

    upper_triangular_solve
    gauss_jordan_solve
    cholesky_solve
    diagonal_solve
    LDLsolve
    LUsolve
    QRsolve
    pinv_solve
    """

    if not M.is_square:
        raise NonSquareMatrixError("Matrix must be square.")
    if rhs.rows != M.rows:
        raise ShapeError("Matrices size mismatch.")
    if not M.is_lower:
        raise ValueError("Matrix must be lower triangular.")

    dps  = _dotprodsimp if dotprodsimp else lambda x: x
    rows = [[] for i in range(M.rows)]

    for i, j, v in M.row_list():
        if i > j:
            rows[i].append((j, v))

    X = rhs.as_mutable()

    for j in range(rhs.cols):
        for i in range(rhs.rows):
            for u, v in rows[i]:
                X[i, j] -= v*X[u, j]

            X[i, j] = dps(X[i, j] / M[i, i])

    return M._new(X)


def _upper_triangular_solve(M, rhs, dotprodsimp=None):
    """Solves ``Ax = B``, where A is an upper triangular matrix.

    Parameters
    ==========

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    See Also
    ========

    lower_triangular_solve
    gauss_jordan_solve
    cholesky_solve
    diagonal_solve
    LDLsolve
    LUsolve
    QRsolve
    pinv_solve
    """

    from .dense import MutableDenseMatrix

    if not M.is_square:
        raise NonSquareMatrixError("Matrix must be square.")
    if rhs.rows != M.rows:
        raise ShapeError("Matrix size mismatch.")
    if not M.is_upper:
        raise TypeError("Matrix is not upper triangular.")

    dps = _dotprodsimp if dotprodsimp else lambda x: x
    X   = MutableDenseMatrix.zeros(M.rows, rhs.cols)

    for j in range(rhs.cols):
        for i in reversed(range(M.rows)):
            if M[i, i] == 0:
                raise ValueError("Matrix must be non-singular.")

            X[i, j] = dps((rhs[i, j] - sum(M[i, k]*X[k, j]
                                        for k in range(i + 1, M.rows))) / M[i, i])

    return M._new(X)


def _upper_triangular_solve_sparse(M, rhs, dotprodsimp=None):
    """Solves ``Ax = B``, where A is an upper triangular matrix.

    Parameters
    ==========

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    See Also
    ========

    lower_triangular_solve
    gauss_jordan_solve
    cholesky_solve
    diagonal_solve
    LDLsolve
    LUsolve
    QRsolve
    pinv_solve
    """

    if not M.is_square:
        raise NonSquareMatrixError("Matrix must be square.")
    if rhs.rows != M.rows:
        raise ShapeError("Matrix size mismatch.")
    if not M.is_upper:
        raise TypeError("Matrix is not upper triangular.")

    dps  = _dotprodsimp if dotprodsimp else lambda x: x
    rows = [[] for i in range(M.rows)]

    for i, j, v in M.row_list():
        if i < j:
            rows[i].append((j, v))

    X = rhs.as_mutable()

    for j in range(rhs.cols):
        for i in reversed(range(rhs.rows)):
            for u, v in reversed(rows[i]):
                X[i, j] -= v*X[u, j]

            X[i, j] = dps(X[i, j] / M[i, i])

    return M._new(X)


def _cholesky_solve(M, rhs, dotprodsimp=None):
    """Solves ``Ax = B`` using Cholesky decomposition,
    for a general square non-singular matrix.
    For a non-square matrix with rows > cols,
    the least squares solution is returned.

    Parameters
    ==========

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    See Also
    ========

    sympy.matrices.dense.DenseMatrix.lower_triangular_solve
    sympy.matrices.dense.DenseMatrix.upper_triangular_solve
    gauss_jordan_solve
    diagonal_solve
    LDLsolve
    LUsolve
    QRsolve
    pinv_solve
    """

    hermitian = True

    if M.is_symmetric():
        hermitian = False
        L = M.cholesky(hermitian=hermitian, dotprodsimp=dotprodsimp)
    elif M.is_hermitian:
        L = M.cholesky(hermitian=hermitian, dotprodsimp=dotprodsimp)
    elif M.rows >= M.cols:
        L = M.H.multiply(M, dotprodsimp=dotprodsimp).cholesky(
                hermitian=hermitian, dotprodsimp=dotprodsimp)
        rhs = M.H.multiply(rhs, dotprodsimp=dotprodsimp)
    else:
        raise NotImplementedError('Under-determined System. '
                                    'Try M.gauss_jordan_solve(rhs)')

    Y = L.lower_triangular_solve(rhs, dotprodsimp=dotprodsimp)

    if hermitian:
        return (L.H).upper_triangular_solve(Y, dotprodsimp=dotprodsimp)
    else:
        return (L.T).upper_triangular_solve(Y, dotprodsimp=dotprodsimp)
