"""
Linear algebra
--------------

Linear equations
................

Basic linear algebra is implemented; you can for example solve the linear
equation system::

      x + 2*y = -10
    3*x + 4*y =  10

using ``lu_solve``::

    >>> A = matrix([[1, 2], [3, 4]])
    >>> b = matrix([-10, 10])
    >>> x = lu_solve(A, b)
    >>> x
    matrix(
    [['30.0'],
     ['-20.0']])

If you don't trust the result, use ``residual`` to calculate the residual ||A*x-b||::

    >>> residual(A, x, b)
    matrix(
    [['3.46944695195361e-18'],
     ['3.46944695195361e-18']])
    >>> str(eps)
    '2.22044604925031e-16'

As you can see, the solution is quite accurate. The error is caused by the
inaccuracy of the internal floating point arithmetic. Though, it's even smaller
than the current machine epsilon, which basically means you can trust the
result.

If you need more speed, use NumPy. Or choose a faster data type using the
keyword ``force_type``::

    >>> lu_solve(A, b, force_type=float)
    matrix(
    [[29.999999999999996],
     [-19.999999999999996]])

``lu_solve`` accepts overdetermined systems. It is usually not possible to solve
such systems, so the residual is minimized instead. Internally this is done
using Cholesky decomposition to compute a least squares approximation. This means
that that ``lu_solve`` will square the errors. If you can't afford this, use
``qr_solve`` instead. It is twice as slow but more accurate, and it calculates
the residual automatically.


Matrix factorization
....................

The function ``lu`` computes an explicit LU factorization of a matrix::

    >>> P, L, U = lu(matrix([[0,2,3],[4,5,6],[7,8,9]]))
    >>> print P
    [0.0  0.0  1.0]
    [1.0  0.0  0.0]
    [0.0  1.0  0.0]
    >>> print L
    [              1.0                0.0  0.0]
    [              0.0                1.0  0.0]
    [0.571428571428571  0.214285714285714  1.0]
    >>> print U
    [7.0  8.0                9.0]
    [0.0  2.0                3.0]
    [0.0  0.0  0.214285714285714]
    >>> print P.T*L*U
    [0.0  2.0  3.0]
    [4.0  5.0  6.0]
    [7.0  8.0  9.0]

Interval matrices
-----------------

Matrices may contain interval elements. This allows one to perform
basic linear algebra operations such as matrix multiplication
and equation solving with rigorous error bounds::

    >>> a = matrix([['0.1','0.3','1.0'],
    ...             ['7.1','5.5','4.8'],
    ...             ['3.2','4.4','5.6']], force_type=mpi)
    >>>
    >>> b = matrix(['4','0.6','0.5'], force_type=mpi)
    >>> c = lu_solve(a, b)
    >>> c
    matrix(
    [[[5.2582327113062393041, 5.2582327113062749951]],
     [[-13.155049396267856583, -13.155049396267821167]],
     [[7.4206915477497212555, 7.4206915477497310922]]])
    >>> print a*c
    [  [3.9999999999999866773, 4.0000000000000133227]]
    [[0.59999999999972430942, 0.60000000000027142733]]
    [[0.49999999999982236432, 0.50000000000018474111]]

"""

# TODO:
# *implement high-level qr()
# *test unitvector
# *iterative solving

from __future__ import division

from mptypes import extraprec, absmin, mp, eps, mpf, fsum
from functions import sqrt, sign, log, factorial
from matrices import matrix, eye, swap_row, extend, norm, mnorm
from copy import copy

def LU_decomp(A, overwrite=False, use_cache=True):
    """
    LU-factorization of a n*n matrix using the Gauss algorithm.
    Returns L and U in one matrix and the pivot indices.

    Use overwrite to specify whether A will be overwritten with L and U.
    """
    if not A.rows == A.cols:
        raise ValueError('need n*n matrix')
    # get from cache if possible
    if use_cache and isinstance(A, matrix) and A._LU:
        return A._LU
    if not overwrite:
        orig = A
        A = A.copy()
    tol = absmin(mnorm(A,1) * eps) # each pivot element has to be bigger
    n = A.rows
    p = [None]*(n - 1)
    for j in xrange(n - 1):
        # pivoting, choose max(abs(reciprocal row sum)*abs(pivot element))
        biggest = 0
        for k in xrange(j, n):
            s = fsum([absmin(A[k,l]) for l in xrange(j, n)])
            if absmin(s) <= tol:
                raise ZeroDivisionError('matrix is numerically singular')
            current = 1/s * absmin(A[k,j])
            if current > biggest: # TODO: what if equal?
                biggest = current
                p[j] = k
        # swap rows according to p
        swap_row(A, j, p[j])
        if absmin(A[j,j]) <= tol:
            raise ZeroDivisionError('matrix is numerically singular')
        # calculate elimination factors and add rows
        for i in xrange(j + 1, n):
            A[i,j] /= A[j,j]
            for k in xrange(j + 1, n):
                A[i,k] -= A[i,j]*A[j,k]
    if absmin(A[n - 1,n - 1]) <= tol:
        raise ZeroDivisionError('matrix is numerically singular')
    # cache decomposition
    if not overwrite and isinstance(orig, matrix):
        orig._LU = (A, p)
    return A, p

def L_solve(L, b, p=None):
    """
    Solve the lower part of a LU factorized matrix for y.
    """
    assert L.rows == L.cols, 'need n*n matrix'
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

def improve_solution(A, x, b, maxsteps=1):
    """
    Improve a solution to a linear equation system iteratively.

    This re-uses the LU decomposition and is thus cheap.
    Usually 3 up to 4 iterations are giving the maximal improvement.
    """
    assert A.rows == A.cols, 'need n*n matrix' # TODO: really?
    for _ in xrange(maxsteps):
        r = residual(A, x, b)
        if norm(r, 2) < 10*eps:
            break
        # this uses cached LU decomposition and is thus cheap
        dx = lu_solve(A, -r)
        x += dx
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
    A, p = LU_decomp(A)
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
        s = fsum((A[i,j])**2 for i in xrange(j, m))
        if not abs(s) > eps:
            raise ValueError('matrix is numerically singular')
        p.append(-sign(A[j,j]) * sqrt(s))
        kappa = s - p[j] * A[j,j]
        A[j,j] -= p[j]
        for k in xrange(j+1, n):
            y = fsum(A[i,j] * A[i,k] for i in xrange(j, m)) /  kappa
            for i in xrange(j, m):
                A[i,k] -= A[i,j] * y
    # solve Rx = c1
    x = [A[i,n - 1] for i in xrange(n - 1)]
    for i in xrange(n - 2, -1, -1):
        x[i] -= fsum(A[i,j] * x[j] for j in xrange(i + 1, n - 1))
        x[i] /= p[i]
    # calculate residual
    if not m == n - 1:
        r = [A[m-1-i, n-1] for i in xrange(m - n + 1)]
    else:
        # determined system, residual should be 0
        r = [0]*m # maybe a bad idea, changing r[i] will change all elements
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
def qr_solve(A, b, norm=norm, **kwargs):
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
        s = A[j,j] - fsum(L[j,k]**2 for k in xrange(j))
        if s < eps:
            raise ValueError('matrix not positive-definite')
        L[j,j] = sqrt(s)
        for i in xrange(j, n):
            L[i,j] = (A[i,j] - fsum(L[i,k] * L[j,k] for k in xrange(j))) \
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
        b[i] -= fsum(L[i,j] * b[j] for j in xrange(i))
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

def cond(A, norm=lambda x: mnorm(x,1)):
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


def lu_solve_mat(a, b):
    """Solve a * x = b  where a and b are matrices."""
    r = matrix(a.rows, b.cols)
    for i in range(b.cols):
        c = lu_solve(a, b.column(i))
        for j in range(len(c)):
            r[j, i] = c[j]
    return r

def exp_pade(a):
    """Exponential of a matrix using Pade approximants.

       See G. H. Golub, C. F. van Loan 'Matrix Computations',
         third Ed., page 572

       TODO:
         - find a good estimate for q
         - reduce the number of matrix multiplications to improve
           performance
    """
    def eps_pade(p):
        return mpf(2)**(3-2*p) * factorial(p)**2/(factorial(2*p)**2 * (2*p + 1))
    q = 4
    extraq = 8
    while 1:
        if eps_pade(q) < eps:
            break
        q += 1
    q += extraq
    j = max(1, int(log(mnorm(a,'inf'),2)))
    extra = q
    mp.dps += extra
    try:
        a = a/2**j
        na = a.rows
        den = eye(na)
        num = eye(na)
        x = eye(na)
        c = mpf(1)
        for k in range(1, q+1):
            c *= mpf(q - k + 1)/((2*q - k + 1) * k)
            x = a*x
            cx = c*x
            num += cx
            den += (-1)**k * cx
        f = lu_solve_mat(den, num)
        for k in range(j):
            f = f*f
    finally:
        mp.dps -= extra
    return f

