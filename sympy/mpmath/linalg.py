# TODO:
# *implement high-level qr()
# *test unitvector
# *iterative solving
# *iterative improving of solution

from __future__ import division

from mptypes import extraprec, absmin, mp, eps
from functions import sqrt, sign
from matrices import matrix, eye, swap_row, extend, mnorm_1, norm_p
from copy import copy

def LU_decomp(A, overwrite=False):
    """
    LU-factorization of a n*n matrix using the Gauss algorithm.
    Returns L and U in one matrix and the pivot indices.

    Use overwrite to specify whether A will be overwritten with L and U.
    """
    if not A.rows == A.cols:
        raise ValueError('need n*n matrix')
    # get from cache if possilbe
    if isinstance(A, matrix) and A._LU:
        return A._LU
    if not overwrite:
        orig = A
        A = A.copy()
    tol = absmin(mnorm_1(A) * eps) # each pivot element has to be bigger
    n = A.rows
    p = [None]*(n - 1)
    for j in xrange(n - 1):
        # pivoting, choose max(abs(reciprocal row sum)*abs(pivot element))
        biggest = 0
        for k in xrange(j, n):
            current = 1/sum([absmin(A[k,l]) for l in xrange(j, n)]) \
                      * absmin(A[k,j])
            if current > biggest: # TODO: what if equal?
                biggest = current
                p[j] = k
        # swap rows according to p
        swap_row(A, j, p[j])
        if absmin(A[j,j]) < tol:
            raise ZeroDivisionError('matrix is numerically singular')
        # calculate elimination factors and add rows
        for i in xrange(j + 1, n):
            A[i,j] /= A[j,j]
            for k in xrange(j + 1, n):
                A[i,k] -= A[i,j]*A[j,k]
    # cache decomposition
    if not overwrite and isinstance(orig, matrix):
        orig._LU = (A, p)
    return A, p

def L_solve(L, b, p=None):
    """
    Solve the lower part of a LU factorized matrix for y.
    """
    L.rows == L.cols, 'need n*n matrix'
    n = L.rows
    assert len(b) == n
    b = copy(b)
    if p: # swap b according to p
        for k in xrange(0, len(p)):
            swap_row(b, k, p[k])
    # solve
    for i in xrange(1, n):
        for j in xrange(i):
            b[i] -= L[i,j] * b[j]
    return b

def U_solve(U, y):
    """
    Solve the upper part of a LU factorized matrix for x.
    """
    assert U.rows == U.cols, 'need n*n matrix'
    n = U.rows
    assert len(y) == n
    x = copy(y)
    for i in xrange(n - 1, -1, -1):
        for j in xrange(i + 1, n):
            x[i] -= U[i,j] * x[j]
        x[i] /= U[i,i]
    return x

@extraprec(10)
def lu_solve(A, b, **kwargs):
    """
    Ax = b => x

    Solve a determined or overdetermined linear equations system.
    Fast LU decomposition is used, which is less accurate than QR decomposition
    (especially for overdetermined systems), but it's twice as efficient.
    Use qr_solve if you want more precison or have to solve a very ill-
    conditioned system.
    """
    # do not overwrite A nor b
    A, b = matrix(A, **kwargs).copy(), matrix(b, **kwargs).copy()
    if A.rows < A.cols:
        raise ValueError('cannot solve underdetermined system')
    if A.rows > A.cols:
        # use least-squares method if overdetermined
        # (this increases errors)
        AT = A.T
        A = AT * A
        b = AT * b
        return cholesky_solve(A, b)
    else:
        # LU factorization
        A, p = LU_decomp(A)
        b = L_solve(A, b, p)
        x = U_solve(A, b)
        return x

def lu(A):
    """
    A -> P, L, U

    LU factorisation of a square matrix A. L is the lower, U the upper part.
    P is the permutation matrix indicating the row swaps.

    P*A = L*U

    If you need efficiency, use the low-level method LU_decomp instead, it's
    much more memory efficient.
    """
    # get factorization
    A, p = LU_decomp(A.copy())
    n = A.rows
    L = matrix(n)
    U = matrix(n)
    for i in xrange(n):
        for j in xrange(n):
            if i > j:
                L[i,j] = A[i,j]
            elif i == j:
                L[i,j] = 1
                U[i,j] = A[i,j]
            else:
                U[i,j] = A[i,j]
    # calculate permutation matrix
    P = eye(n)
    for k in xrange(len(p)):
        swap_row(P, k, p[k])
    return P, L, U

def unitvector(n, i):
    """
    Return the i-th n-dimensional unit vector.
    """
    assert 0 < i <= n, 'this unit vector does not exist'
    return [0]*(i-1) + [1] + [0]*(n-i)

@extraprec(10)
def inverse(A, **kwargs):
    """
    Calculate the inverse of a matrix.

    If you want to solve an equation system Ax = b, it's recommended to use
    solve(A, b) instead, it's about 3 times more efficient.
    """
    # do not overwrite A
    A = matrix(A, **kwargs).copy()
    n = A.rows
    # get LU factorisation
    A, p = LU_decomp(A)
    cols = []
    # calculate unit vectors and solve corresponding system to get columns
    for i in xrange(1, n + 1):
        e = unitvector(n, i)
        y = L_solve(A, e, p)
        cols.append(U_solve(A, y))
    # convert columns to matrix
    inv = []
    for i in xrange(n):
        row = []
        for j in xrange(n):
            row.append(cols[j][i])
        inv.append(row)
    return matrix(inv, **kwargs)

def householder(A):
    """
    (A|b) -> H, p, x, res

    (A|b) is the coefficient matrix with left hand side of an optionally
    overdetermined linear equation system.
    H and p contain all information about the transformation matrices.
    x is the solution, res the residual.
    """
    assert isinstance(A, matrix)
    m = A.rows
    n = A.cols
    assert m >= n - 1
    # calculate Householder matrix
    p = []
    for j in xrange(0, n - 1):
        s = 0.
        for i in xrange(j, m):
            s += (A[i,j])**2
        if not abs(s) > eps:
            raise ValueError('matrix is numerically singular')
        p.append(-sign(A[j,j]) * sqrt(s))
        kappa = s - p[j] * A[j,j]
        A[j,j] -= p[j]
        for k in xrange(j+1, n):
            y = 0.
            for i in xrange(j, m):
                y += A[i,j] * A[i,k]
            y /= kappa
            for i in xrange(j, m):
                A[i,k] -= A[i,j] * y
    # solve Rx = c1
    x = []
    for i in xrange(n - 1):
        x.append(A[i,n - 1])
    for i in xrange(n - 2, -1, -1):
        for j in xrange(i + 1, n - 1):
            x[i] -= A[i,j] * x[j]
        x[i] /= p[i]
    # calculate residual
    if not m == n - 1:
        r = []
        for i in xrange(m - n + 1):
            r.append(A[m-1-i, n-1])
    else:
        # determined system, residual should be 0
        r = [0]*m
    return A, p, x, r

#def qr(A):
#    """
#    A -> Q, R
#
#    QR factorisation of a square matrix A using Householder decomposition.
#    Q is orthogonal, this leads to very few numerical errors.
#
#    A = Q*R
#    """
#    H, p, x, res = householder(A)
# TODO: implement this

def residual(A, x, b, **kwargs):
    """
    Calculate the residual of a solution to a linear equation system.

    r = A*x - b for A*x = b
    """
    oldprec = mp.prec
    try:
        mp.prec *= 2
        A, x, b = matrix(A, **kwargs), matrix(x, **kwargs), matrix(b, **kwargs)
        return A*x - b
    finally:
        mp.prec = oldprec

@extraprec(10)
def qr_solve(A, b, norm=lambda x: norm_p(x, 2), **kwargs):
    """
    Ax = b => x, ||Ax - b||

    Solve a determined or overdetermined linear equations system and
    calculate the norm of the residual (error).
    QR decompostion using Householder factorization is applied, which gives very
    accurate results even for ill-conditioned matrices. qr_solve is twice as
    efficient.
    """
    # do not overwrite A nor b
    A, b = matrix(A, **kwargs).copy(), matrix(b, **kwargs).copy()
    if A.rows < A.cols:
        raise ValueError('cannot solve underdetermined system')
    H, p, x, r = householder(extend(A, b))
    res = norm(r)
    # calculate residual "manually" for determined systems
    if res == 0:
        res = norm(residual(A, x, b))
    return matrix(x, **kwargs), res

def cholesky(A):
    """
    Cholesky decompositon of a symmetric positive-definite matrix.

    Can be used to solve linear equation systems twice as efficient compared
    to LU decomposition or to test whether A is positive-definite.

    A = L * L.T
    Only L (the lower part) is returned.
    """
    assert isinstance(A, matrix)
    if not A.rows == A.cols:
        raise ValueError('need n*n matrix')
    n = A.rows
    L = matrix(n)
    for j in xrange(n):
        s = A[j,j] - sum((L[j,k]**2 for k in xrange(j)))
        if s < eps:
            raise ValueError('matrix not positive-definite')
        L[j,j] = sqrt(s)
        for i in xrange(j, n):
            L[i,j] = (A[i,j] - sum((L[i,k] * L[j,k] for k in xrange(j)))) \
                     / L[j,j]
    return L

@extraprec(10)
def cholesky_solve(A, b, **kwargs):
    """
    Ax = b => x

    Solve a symmetric positive-definite linear equation system.
    This is twice as efficient as lu_solve.

    Typical use cases:
    * A.T*A
    * Hessian matrix
    * differential equations
    """
    # do not overwrite A nor b
    A, b = matrix(A, **kwargs).copy(), matrix(b, **kwargs).copy()
    if A.rows !=  A.cols:
        raise ValueError('can only solve determined system')
    # Cholesky factorization
    L = cholesky(A)
    # solve
    n = L.rows
    assert len(b) == n
    for i in xrange(n):
        b[i] -= sum((L[i,j] * b[j] for j in xrange(i)))
        b[i] /= L[i,i]
    x = U_solve(L.T, b)
    return x

@extraprec(10)
def det(A):
    """
    Calculate the determinant of a matrix.
    """
    # do not overwrite A
    A = matrix(A).copy()
    # use LU factorization to calculate determinant
    try:
        R, p = LU_decomp(A)
    except ZeroDivisionError:
        return 0
    z = 1
    for i, e in enumerate(p):
        if i != e:
            z *= -1
    for i in xrange(A.rows):
        z *= R[i,i]
    return z

def cond(A, norm=mnorm_1):
    """
    Calculate the condition number of a matrix using a specified matrix norm.

    The condition number estimates the sensitivity of a matrix to errors.
    Example: small input errors for ill-conditionded coefficient matrices
    alter the solution of the system dramatically.

    For ill-conditioned matrices it's recommended to use qr_solve() instead
    of lu_solve(). This does not help with input errors however, it just avoids
    to add additional errors.

    Definition:    cond(A) = ||A|| * ||A**-1||
    """
    return norm(A) * norm(inverse(A))


