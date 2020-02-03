from __future__ import division, print_function

from sympy.core.function import expand_mul
from sympy.core.symbol import _uniquely_named_symbol
from sympy.simplify.simplify import dotprodsimp as _dotprodsimp
from sympy.utilities.iterables import numbered_symbols

from .common import ShapeError, NonSquareMatrixError
from .utilities import _iszero


def _diagonal_solve(M, rhs):
    """Solves ``Ax = B`` efficiently, where A is a diagonal Matrix,
    with non-zero diagonal entries.

    Examples
    ========

    >>> from sympy.matrices import Matrix, eye
    >>> A = eye(2)*2
    >>> B = Matrix([[1, 2], [3, 4]])
    >>> A.diagonal_solve(B) == B/2
    True

    See Also
    ========

    sympy.matrices.dense.DenseMatrix.lower_triangular_solve
    sympy.matrices.dense.DenseMatrix.upper_triangular_solve
    gauss_jordan_solve
    cholesky_solve
    LDLsolve
    LUsolve
    QRsolve
    pinv_solve
    """

    if not M.is_diagonal():
        raise TypeError("Matrix should be diagonal")
    if rhs.rows != M.rows:
        raise TypeError("Size mis-match")

    return M._new(
        rhs.rows, rhs.cols, lambda i, j: rhs[i, j] / M[i, i])


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


def _LDLsolve(M, rhs, dotprodsimp=None):
    """Solves ``Ax = B`` using LDL decomposition,
    for a general square and non-singular matrix.

    For a non-square matrix with rows > cols,
    the least squares solution is returned.

    Parameters
    ==========

    dotprodsimp : bool, optional
        Specifies whether intermediate term algebraic simplification is used
        during matrix multiplications to control expression blowup and thus
        speed up calculation.

    Examples
    ========

    >>> from sympy.matrices import Matrix, eye
    >>> A = eye(2)*2
    >>> B = Matrix([[1, 2], [3, 4]])
    >>> A.LDLsolve(B) == B/2
    True

    See Also
    ========

    sympy.matrices.dense.DenseMatrix.LDLdecomposition
    sympy.matrices.dense.DenseMatrix.lower_triangular_solve
    sympy.matrices.dense.DenseMatrix.upper_triangular_solve
    gauss_jordan_solve
    cholesky_solve
    diagonal_solve
    LUsolve
    QRsolve
    pinv_solve
    """

    hermitian = True

    if M.is_symmetric():
        hermitian = False
        L, D = M.LDLdecomposition(hermitian=hermitian, dotprodsimp=dotprodsimp)
    elif M.is_hermitian:
        L, D = M.LDLdecomposition(hermitian=hermitian, dotprodsimp=dotprodsimp)
    elif M.rows >= M.cols:
        L, D = M.H.multiply(M, dotprodsimp=dotprodsimp) \
                .LDLdecomposition(hermitian=hermitian, dotprodsimp=dotprodsimp)
        rhs = M.H.multiply(rhs, dotprodsimp=dotprodsimp)
    else:
        raise NotImplementedError('Under-determined System. '
                                    'Try M.gauss_jordan_solve(rhs)')

    Y = L.lower_triangular_solve(rhs, dotprodsimp=dotprodsimp)
    Z = D.diagonal_solve(Y)

    if hermitian:
        return (L.H).upper_triangular_solve(Z, dotprodsimp=dotprodsimp)
    else:
        return (L.T).upper_triangular_solve(Z, dotprodsimp=dotprodsimp)


def _LUsolve(M, rhs, iszerofunc=_iszero, dotprodsimp=None):
    """Solve the linear system ``Ax = rhs`` for ``x`` where ``A = M``.

    This is for symbolic matrices, for real or complex ones use
    mpmath.lu_solve or mpmath.qr_solve.

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
    cholesky_solve
    diagonal_solve
    LDLsolve
    QRsolve
    pinv_solve
    LUdecomposition
    """

    if rhs.rows != M.rows:
        raise ShapeError(
            "``M`` and ``rhs`` must have the same number of rows.")

    m = M.rows
    n = M.cols

    if m < n:
        raise NotImplementedError("Underdetermined systems not supported.")

    try:
        A, perm = M.LUdecomposition_Simple(
            iszerofunc=_iszero, rankcheck=True, dotprodsimp=dotprodsimp)
    except ValueError:
        raise NotImplementedError("Underdetermined systems not supported.")

    dps = _dotprodsimp if dotprodsimp else lambda e: e
    b   = rhs.permute_rows(perm).as_mutable()

    # forward substitution, all diag entries are scaled to 1
    for i in range(m):
        for j in range(min(i, n)):
            scale = A[i, j]
            b.zip_row_op(i, j, lambda x, y: dps(x - y * scale))

    # consistency check for overdetermined systems
    if m > n:
        for i in range(n, m):
            for j in range(b.cols):
                if not iszerofunc(b[i, j]):
                    raise ValueError("The system is inconsistent.")

        b = b[0:n, :]   # truncate zero rows if consistent

    # backward substitution
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            scale = A[i, j]
            b.zip_row_op(i, j, lambda x, y: dps(x - y * scale))

        scale = A[i, i]
        b.row_op(i, lambda x, _: dps(x / scale))

    return rhs.__class__(b)


def _QRsolve(M, b, dotprodsimp=None):
    """Solve the linear system ``Ax = b``.

    ``M`` is the matrix ``A``, the method argument is the vector
    ``b``.  The method returns the solution vector ``x``.  If ``b`` is a
    matrix, the system is solved for each column of ``b`` and the
    return value is a matrix of the same shape as ``b``.

    This method is slower (approximately by a factor of 2) but
    more stable for floating-point arithmetic than the LUsolve method.
    However, LUsolve usually uses an exact arithmetic, so you don't need
    to use QRsolve.

    This is mainly for educational purposes and symbolic matrices, for real
    (or complex) matrices use mpmath.qr_solve.

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
    cholesky_solve
    diagonal_solve
    LDLsolve
    LUsolve
    pinv_solve
    QRdecomposition
    """

    dps  = _dotprodsimp if dotprodsimp else expand_mul
    Q, R = M.QRdecomposition(dotprodsimp=dotprodsimp)
    y    = Q.T * b

    # back substitution to solve R*x = y:
    # We build up the result "backwards" in the vector 'x' and reverse it
    # only in the end.
    x = []
    n = R.rows

    for j in range(n - 1, -1, -1):
        tmp = y[j, :]

        for k in range(j + 1, n):
            tmp -= R[j, k] * x[n - 1 - k]

        tmp = dps(tmp)

        x.append(tmp / R[j, j])

    return M._new([row._mat for row in reversed(x)])


def _gauss_jordan_solve(M, B, freevar=False):
    """
    Solves ``Ax = B`` using Gauss Jordan elimination.

    There may be zero, one, or infinite solutions.  If one solution
    exists, it will be returned. If infinite solutions exist, it will
    be returned parametrically. If no solutions exist, It will throw
    ValueError.

    Parameters
    ==========

    B : Matrix
        The right hand side of the equation to be solved for.  Must have
        the same number of rows as matrix A.

    freevar : List
        If the system is underdetermined (e.g. A has more columns than
        rows), infinite solutions are possible, in terms of arbitrary
        values of free variables. Then the index of the free variables
        in the solutions (column Matrix) will be returned by freevar, if
        the flag `freevar` is set to `True`.

    Returns
    =======

    x : Matrix
        The matrix that will satisfy ``Ax = B``.  Will have as many rows as
        matrix A has columns, and as many columns as matrix B.

    params : Matrix
        If the system is underdetermined (e.g. A has more columns than
        rows), infinite solutions are possible, in terms of arbitrary
        parameters. These arbitrary parameters are returned as params
        Matrix.

    Examples
    ========

    >>> from sympy import Matrix
    >>> A = Matrix([[1, 2, 1, 1], [1, 2, 2, -1], [2, 4, 0, 6]])
    >>> B = Matrix([7, 12, 4])
    >>> sol, params = A.gauss_jordan_solve(B)
    >>> sol
    Matrix([
    [-2*tau0 - 3*tau1 + 2],
    [                 tau0],
    [           2*tau1 + 5],
    [                 tau1]])
    >>> params
    Matrix([
    [tau0],
    [tau1]])
    >>> taus_zeroes = { tau:0 for tau in params }
    >>> sol_unique = sol.xreplace(taus_zeroes)
    >>> sol_unique
        Matrix([
    [2],
    [0],
    [5],
    [0]])


    >>> A = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
    >>> B = Matrix([3, 6, 9])
    >>> sol, params = A.gauss_jordan_solve(B)
    >>> sol
    Matrix([
    [-1],
    [ 2],
    [ 0]])
    >>> params
    Matrix(0, 1, [])

    >>> A = Matrix([[2, -7], [-1, 4]])
    >>> B = Matrix([[-21, 3], [12, -2]])
    >>> sol, params = A.gauss_jordan_solve(B)
    >>> sol
    Matrix([
    [0, -2],
    [3, -1]])
    >>> params
    Matrix(0, 2, [])

    See Also
    ========

    sympy.matrices.dense.DenseMatrix.lower_triangular_solve
    sympy.matrices.dense.DenseMatrix.upper_triangular_solve
    cholesky_solve
    diagonal_solve
    LDLsolve
    LUsolve
    QRsolve
    pinv

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Gaussian_elimination

    """

    from sympy.matrices import Matrix, zeros

    cls      = M.__class__
    aug      = M.hstack(M.copy(), B.copy())
    B_cols   = B.cols
    row, col = aug[:, :-B_cols].shape

    # solve by reduced row echelon form
    A, pivots = aug.rref(simplify=True)
    A, v      = A[:, :-B_cols], A[:, -B_cols:]
    pivots    = list(filter(lambda p: p < col, pivots))
    rank      = len(pivots)

    # Bring to block form
    permutation = Matrix(range(col)).T

    for i, c in enumerate(pivots):
        permutation.col_swap(i, c)

    # check for existence of solutions
    # rank of aug Matrix should be equal to rank of coefficient matrix
    if not v[rank:, :].is_zero_matrix:
        raise ValueError("Linear system has no solution")

    # Get index of free symbols (free parameters)
    # non-pivots columns are free variables
    free_var_index = permutation[len(pivots):]

    # Free parameters
    # what are current unnumbered free symbol names?
    name = _uniquely_named_symbol('tau', aug,
            compare=lambda i: str(i).rstrip('1234567890')).name
    gen  = numbered_symbols(name)
    tau  = Matrix([next(gen) for k in range((col - rank)*B_cols)]).reshape(
            col - rank, B_cols)

    # Full parametric solution
    V        = A[:rank, [c for c in range(A.cols) if c not in pivots]]
    vt       = v[:rank, :]
    free_sol = tau.vstack(vt - V * tau, tau)

    # Undo permutation
    sol = zeros(col, B_cols)

    for k in range(col):
        sol[permutation[k], :] = free_sol[k,:]

    sol, tau = cls(sol), cls(tau)

    if freevar:
        return sol, tau, free_var_index
    else:
        return sol, tau
