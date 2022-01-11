from sympy.core.numbers import mod_inverse

from .common import MatrixError, NonSquareMatrixError, NonInvertibleMatrixError
from .utilities import _iszero


def _pinv_full_rank(M):
    """Subroutine for full row or column rank matrices.

    For full row rank matrices, inverse of ``A * A.H`` Exists.
    For full column rank matrices, inverse of ``A.H * A`` Exists.

    This routine can apply for both cases by checking the shape
    and have small decision.
    """

    if M.is_zero_matrix:
        return M.H

    if M.rows >= M.cols:
        return M.H.multiply(M).inv().multiply(M.H)
    else:
        return M.H.multiply(M.multiply(M.H).inv())

def _pinv_rank_decomposition(M):
    """Subroutine for rank decomposition

    With rank decompositions, `A` can be decomposed into two full-
    rank matrices, and each matrix can take pseudoinverse
    individually.
    """

    if M.is_zero_matrix:
        return M.H

    B, C = M.rank_decomposition()

    Bp = _pinv_full_rank(B)
    Cp = _pinv_full_rank(C)

    return Cp.multiply(Bp)

def _pinv_diagonalization(M):
    """Subroutine using diagonalization

    This routine can sometimes fail if SymPy's eigenvalue
    computation is not reliable.
    """

    if M.is_zero_matrix:
        return M.H

    A  = M
    AH = M.H

    try:
        if M.rows >= M.cols:
            P, D   = AH.multiply(A).diagonalize(normalize=True)
            D_pinv = D.applyfunc(lambda x: 0 if _iszero(x) else 1 / x)

            return P.multiply(D_pinv).multiply(P.H).multiply(AH)

        else:
            P, D   = A.multiply(AH).diagonalize(
                        normalize=True)
            D_pinv = D.applyfunc(lambda x: 0 if _iszero(x) else 1 / x)

            return AH.multiply(P).multiply(D_pinv).multiply(P.H)

    except MatrixError:
        raise NotImplementedError(
            'pinv for rank-deficient matrices where '
            'diagonalization of A.H*A fails is not supported yet.')

def _pinv(M, method='RD'):
    """Calculate the Moore-Penrose pseudoinverse of the matrix.

    The Moore-Penrose pseudoinverse exists and is unique for any matrix.
    If the matrix is invertible, the pseudoinverse is the same as the
    inverse.

    Parameters
    ==========

    method : String, optional
        Specifies the method for computing the pseudoinverse.

        If ``'RD'``, Rank-Decomposition will be used.

        If ``'ED'``, Diagonalization will be used.

    Examples
    ========

    Computing pseudoinverse by rank decomposition :

    >>> from sympy import Matrix
    >>> A = Matrix([[1, 2, 3], [4, 5, 6]])
    >>> A.pinv()
    Matrix([
    [-17/18,  4/9],
    [  -1/9,  1/9],
    [ 13/18, -2/9]])

    Computing pseudoinverse by diagonalization :

    >>> B = A.pinv(method='ED')
    >>> B.simplify()
    >>> B
    Matrix([
    [-17/18,  4/9],
    [  -1/9,  1/9],
    [ 13/18, -2/9]])

    See Also
    ========

    inv
    pinv_solve

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Moore-Penrose_pseudoinverse

    """

    # Trivial case: pseudoinverse of all-zero matrix is its transpose.
    if M.is_zero_matrix:
        return M.H

    if method == 'RD':
        return _pinv_rank_decomposition(M)
    elif method == 'ED':
        return _pinv_diagonalization(M)
    else:
        raise ValueError('invalid pinv method %s' % repr(method))


def _inv_mod(M, m):
    r"""
    Returns the inverse of the matrix `K` (mod `m`), if it exists.

    Method to find the matrix inverse of `K` (mod `m`) implemented in this function:

    * Compute `\mathrm{adj}(K) = \mathrm{cof}(K)^t`, the adjoint matrix of `K`.

    * Compute `r = 1/\mathrm{det}(K) \pmod m`.

    * `K^{-1} = r\cdot \mathrm{adj}(K) \pmod m`.

    Examples
    ========

    >>> from sympy import Matrix
    >>> A = Matrix(2, 2, [1, 2, 3, 4])
    >>> A.inv_mod(5)
    Matrix([
    [3, 1],
    [4, 2]])
    >>> A.inv_mod(3)
    Matrix([
    [1, 1],
    [0, 1]])

    """

    if not M.is_square:
        raise NonSquareMatrixError()

    N       = M.cols
    det_K   = M.det()
    det_inv = None

    try:
        det_inv = mod_inverse(det_K, m)
    except ValueError:
        raise NonInvertibleMatrixError('Matrix is not invertible (mod %d)' % m)

    K_adj = M.adjugate()
    K_inv = M.__class__(N, N,
            [det_inv * K_adj[i, j] % m for i in range(N) for j in range(N)])

    return K_inv


def _verify_invertible(M, iszerofunc=_iszero):
    """Initial check to see if a matrix is invertible. Raises or returns
    determinant for use in _inv_ADJ."""

    if not M.is_square:
        raise NonSquareMatrixError("A Matrix must be square to invert.")

    d    = M.det(method='berkowitz')
    zero = d.equals(0)

    if zero is None: # if equals() can't decide, will rref be able to?
        ok   = M.rref(simplify=True)[0]
        zero = any(iszerofunc(ok[j, j]) for j in range(ok.rows))

    if zero:
        raise NonInvertibleMatrixError("Matrix det == 0; not invertible.")

    return d

def _inv_ADJ(M, iszerofunc=_iszero):
    """Calculates the inverse using the adjugate matrix and a determinant.

    See Also
    ========

    inv
    inverse_GE
    inverse_LU
    inverse_CH
    inverse_LDL
    """

    d = _verify_invertible(M, iszerofunc=iszerofunc)

    return M.adjugate() / d

def _inv_GE(M, iszerofunc=_iszero):
    """Calculates the inverse using Gaussian elimination.

    See Also
    ========

    inv
    inverse_ADJ
    inverse_LU
    inverse_CH
    inverse_LDL
    """

    from .dense import Matrix

    if not M.is_square:
        raise NonSquareMatrixError("A Matrix must be square to invert.")

    big = Matrix.hstack(M.as_mutable(), Matrix.eye(M.rows))
    red = big.rref(iszerofunc=iszerofunc, simplify=True)[0]

    if any(iszerofunc(red[j, j]) for j in range(red.rows)):
        raise NonInvertibleMatrixError("Matrix det == 0; not invertible.")

    return M._new(red[:, big.rows:])

def _inv_LU(M, iszerofunc=_iszero):
    """Calculates the inverse using LU decomposition.

    See Also
    ========

    inv
    inverse_ADJ
    inverse_GE
    inverse_CH
    inverse_LDL
    """

    if not M.is_square:
        raise NonSquareMatrixError("A Matrix must be square to invert.")
    if M.free_symbols:
        _verify_invertible(M, iszerofunc=iszerofunc)

    return M.LUsolve(M.eye(M.rows), iszerofunc=_iszero)

def _inv_CH(M, iszerofunc=_iszero):
    """Calculates the inverse using cholesky decomposition.

    See Also
    ========

    inv
    inverse_ADJ
    inverse_GE
    inverse_LU
    inverse_LDL
    """

    _verify_invertible(M, iszerofunc=iszerofunc)

    return M.cholesky_solve(M.eye(M.rows))

def _inv_LDL(M, iszerofunc=_iszero):
    """Calculates the inverse using LDL decomposition.

    See Also
    ========

    inv
    inverse_ADJ
    inverse_GE
    inverse_LU
    inverse_CH
    """

    _verify_invertible(M, iszerofunc=iszerofunc)

    return M.LDLsolve(M.eye(M.rows))

def _inv_QR(M, iszerofunc=_iszero):
    """Calculates the inverse using QR decomposition.

    See Also
    ========

    inv
    inverse_ADJ
    inverse_GE
    inverse_CH
    inverse_LDL
    """

    _verify_invertible(M, iszerofunc=iszerofunc)

    return M.QRsolve(M.eye(M.rows))

def _inv_block(M, iszerofunc=_iszero):
    """Calculates the inverse using BLOCKWISE inversion.

    See Also
    ========

    inv
    inverse_ADJ
    inverse_GE
    inverse_CH
    inverse_LDL
    """
    from sympy.matrices.expressions.blockmatrix import BlockMatrix
    i = M.shape[0]
    if i <= 20 :
        return M.inv(method="LU", iszerofunc=_iszero)
    A = M[:i // 2, :i //2]
    B = M[:i // 2, i // 2:]
    C = M[i // 2:, :i // 2]
    D = M[i // 2:, i // 2:]
    try:
        D_inv = _inv_block(D)
    except NonInvertibleMatrixError:
        return M.inv(method="LU", iszerofunc=_iszero)
    B_D_i = B*D_inv
    BDC = B_D_i*C
    A_n = A - BDC
    try:
        A_n = _inv_block(A_n)
    except NonInvertibleMatrixError:
        return M.inv(method="LU", iszerofunc=_iszero)
    B_n = -A_n*B_D_i
    dc = D_inv*C
    C_n = -dc*A_n
    D_n = D_inv + dc*-B_n
    nn = BlockMatrix([[A_n, B_n], [C_n, D_n]]).as_explicit()
    return nn

def _inv(M, method=None, iszerofunc=_iszero, try_block_diag=False):
    """
    Return the inverse of a matrix using the method indicated. Default for
    dense matrices is is Gauss elimination, default for sparse matrices is LDL.

    Parameters
    ==========

    method : ('GE', 'LU', 'ADJ', 'CH', 'LDL')

    iszerofunc : function, optional
        Zero-testing function to use.

    try_block_diag : bool, optional
        If True then will try to form block diagonal matrices using the
        method get_diag_blocks(), invert these individually, and then
        reconstruct the full inverse matrix.

    Examples
    ========

    >>> from sympy import SparseMatrix, Matrix
    >>> A = SparseMatrix([
    ... [ 2, -1,  0],
    ... [-1,  2, -1],
    ... [ 0,  0,  2]])
    >>> A.inv('CH')
    Matrix([
    [2/3, 1/3, 1/6],
    [1/3, 2/3, 1/3],
    [  0,   0, 1/2]])
    >>> A.inv(method='LDL') # use of 'method=' is optional
    Matrix([
    [2/3, 1/3, 1/6],
    [1/3, 2/3, 1/3],
    [  0,   0, 1/2]])
    >>> A * _
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])
    >>> A = Matrix(A)
    >>> A.inv('CH')
    Matrix([
    [2/3, 1/3, 1/6],
    [1/3, 2/3, 1/3],
    [  0,   0, 1/2]])
    >>> A.inv('ADJ') == A.inv('GE') == A.inv('LU') == A.inv('CH') == A.inv('LDL') == A.inv('QR')
    True

    Notes
    =====

    According to the ``method`` keyword, it calls the appropriate method:

        GE .... inverse_GE(); default for dense matrices
        LU .... inverse_LU()
        ADJ ... inverse_ADJ()
        CH ... inverse_CH()
        LDL ... inverse_LDL(); default for sparse matrices
        QR ... inverse_QR()

    Note, the GE and LU methods may require the matrix to be simplified
    before it is inverted in order to properly detect zeros during
    pivoting. In difficult cases a custom zero detection function can
    be provided by setting the ``iszerofunc`` argument to a function that
    should return True if its argument is zero. The ADJ routine computes
    the determinant and uses that to detect singular matrices in addition
    to testing for zeros on the diagonal.

    See Also
    ========

    inverse_ADJ
    inverse_GE
    inverse_LU
    inverse_CH
    inverse_LDL

    Raises
    ======

    ValueError
        If the determinant of the matrix is zero.
    """

    from sympy.matrices import diag, SparseMatrix

    if method is None:
        method = 'LDL' if isinstance(M, SparseMatrix) else 'GE'

    if try_block_diag:
        blocks = M.get_diag_blocks()
        r      = []

        for block in blocks:
            r.append(block.inv(method=method, iszerofunc=iszerofunc))

        return diag(*r)

    if method == "GE":
        rv = M.inverse_GE(iszerofunc=iszerofunc)
    elif method == "LU":
        rv = M.inverse_LU(iszerofunc=iszerofunc)
    elif method == "ADJ":
        rv = M.inverse_ADJ(iszerofunc=iszerofunc)
    elif method == "CH":
        rv = M.inverse_CH(iszerofunc=iszerofunc)
    elif method == "LDL":
        rv = M.inverse_LDL(iszerofunc=iszerofunc)
    elif method == "QR":
        rv = M.inverse_QR(iszerofunc=iszerofunc)
    elif method == "BLOCK":
        rv = M.inverse_BLOCK(iszerofunc=iszerofunc)
    else:
        raise ValueError("Inversion method unrecognized")

    return M._new(rv)
