from sympy.core.sympify import sympify as sympy_simplify

def _iszero(x):
    return x == 0

#############
## inverse ##
#############

def inverse_LU(self, iszerofunc=_iszero):
    """
    Calculates the inverse using LU decomposition.
    """
    return self.LUsolve(self.eye(self.rows), iszerofunc=_iszero)

def inverse_GE(self, iszerofunc=_iszero):
    """
    Calculates the inverse using Gaussian elimination.
    """
    if not self.is_square:
        raise NonSquareMatrixError()

    if self.det() == 0:
        raise ValueError("A Matrix must have non-zero determinant to invert.")

    big = self.row_join(self.eye(self.rows))
    red = big.rref(iszerofunc=iszerofunc)
    return red[0][:,big.rows:]

def inverse_ADJ(self):
    """
    Calculates the inverse using the adjugate matrix and a determinant.
    """
    if not self.is_square:
        raise NonSquareMatrixError()

    d = berkowitz_det(self)
    if d == 0:
        raise ValueError("A Matrix must have non-zero determinant to invert.")

    return self.adjugate()/d

#############
#### LU #####
#############

def LUsolve(self, rhs, iszerofunc=_iszero):
    """
    Solve the linear system Ax = b for x.
    self is the coefficient matrix A and rhs is the right side b.

    This is for symbolic matrices, for real or complex ones use
    sympy.mpmath.lu_solve or sympy.mpmath.qr_solve.

    """
    if rhs.rows != self.rows:
        raise ShapeError("`self` and `rhs` must have the same number of rows.")

    A, perm = self.LUdecomposition_Simple(iszerofunc=_iszero)
    n = self.rows
    b = rhs.permuteFwd(perm)
    # forward substitution, all diag entries are scaled to 1
    for i in range(n):
        for j in range(i):
            b.row(i, lambda x,k: x - b[j,k]*A[i,j])
    # backward substitution
    for i in range(n-1,-1,-1):
        for j in range(i+1, n):
            b.row(i, lambda x,k: x - b[j,k]*A[i,j])
        b.row(i, lambda x,k: x / A[i,i])
    return b

def LUdecomposition(self, iszerofunc=_iszero):
    """
    Returns the decomposition LU and the row swaps p.

    Example:
    >>> from sympy import Matrix
    >>> a = Matrix([[4, 3], [6, 3]])
    >>> L, U, _ = a.LUdecomposition()
    >>> L
    [  1, 0]
    [3/2, 1]
    >>> U
    [4,    3]
    [0, -3/2]

    """
    combined, p = self.LUdecomposition_Simple(iszerofunc=_iszero)
    L = self.zeros(self.rows)
    U = self.zeros(self.rows)
    for i in range(self.rows):
        for j in range(self.rows):
            if i > j:
                L[i,j] = combined[i,j]
            else:
                if i == j:
                    L[i,i] = 1
                U[i,j] = combined[i,j]
    return L, U, p

def LUdecomposition_Simple(self, iszerofunc=_iszero):
    """
    Returns A comprised of L,U (L's diag entries are 1) and
    p which is the list of the row swaps (in order).
    """
    if not self.is_square:
        raise NonSquareMatrixError()
    n = self.rows
    A = self.copy()
    p = []
    # factorization
    for j in range(n):
        for i in range(j):
            for k in range(i):
                A[i,j] = A[i,j] - A[i,k]*A[k,j]
        pivot = -1
        for i in range(j,n):
            for k in range(j):
                A[i,j] = A[i,j] - A[i,k]*A[k,j]
            # find the first non-zero pivot, includes any expression
            if pivot == -1 and not iszerofunc(A[i,j]):
                pivot = i
        if pivot < 0:
            # this result is based on iszerofunc's analysis of the possible pivots, so even though
            # the element may not be strictly zero, the supplied iszerofunc's evaluation gave True
            raise ValueError("No nonzero pivot found; inversion failed.")
        if pivot != j: # row must be swapped
            A.row_swap(pivot,j)
            p.append([pivot,j])
        scale = 1 / A[j,j]
        for i in range(j+1,n):
            A[i,j] = A[i,j] * scale
    return A, p

def LUdecompositionFF(self):
    """
    Compute a fraction-free LU decomposition.

    Returns 4 matrices P, L, D, U such that PA = L D**-1 U.
    If the elements of the matrix belong to some integral domain I, then all
    elements of L, D and U are guaranteed to belong to I.

    **Reference**
        - W. Zhou & D.J. Jeffrey, "Fraction-free matrix factors: new forms
            for LU and QR factors". Frontiers in Computer Science in China,
            Vol 2, no. 1, pp. 67-80, 2008.
    """
    n, m = self.rows, self.cols
    U, L, P = self.copy(), eye(n), eye(n)
    DD = zeros(n) # store it smarter since it's just diagonal
    oldpivot = 1

    for k in range(n-1):
        if U[k,k] == 0:
            for kpivot in range(k+1, n):
                if U[kpivot, k] != 0:
                    break
            else:
                raise ValueError("Matrix is not full rank")
            U[k, k:], U[kpivot, k:] = U[kpivot, k:], U[k, k:]
            L[k, :k], L[kpivot, :k] = L[kpivot, :k], L[k, :k]
            P[k, :], P[kpivot, :] = P[kpivot, :], P[k, :]
        L[k,k] = Ukk = U[k,k]
        DD[k,k] = oldpivot * Ukk
        for i in range(k+1, n):
            L[i,k] = Uik = U[i,k]
            for j in range(k+1, m):
                U[i,j] = (Ukk * U[i,j] - U[k,j]*Uik) / oldpivot
            U[i,k] = 0
        oldpivot = Ukk
    DD[n-1,n-1] = oldpivot
    return P, L, DD, U

#############
#### det ####
#############

def det_bareis(self):
    """Compute matrix determinant using Bareis' fraction-free
        algorithm which is an extension of the well known Gaussian
        elimination method. This approach is best suited for dense
        symbolic matrices and will result in a determinant with
        minimal number of fractions. It means that less term
        rewriting is needed on resulting formulae.

        TODO: Implement algorithm for sparse matrices (SFF).
    """
    if not self.is_square:
        raise NonSquareMatrixError()

    M, n = self[:,:], self.rows

    if n == 1:
        det = M[0, 0]
    elif n == 2:
        det = M[0, 0]*M[1, 1] - M[0, 1]*M[1, 0]
    else:
        sign = 1 # track current sign in case of column swap

        for k in range(n-1):
            # look for a pivot in the current column
            # and assume det == 0 if none is found
            if M[k, k] == 0:
                for i in range(k+1, n):
                    if M[i, k] != 0:
                        M.row_swap(i, k)
                        sign *= -1
                        break
                else:
                    return S.Zero

            # proceed with Bareis' fraction-free (FF)
            # form of Gaussian elimination algorithm
            for i in range(k+1, n):
                for j in range(k+1, n):
                    D = M[k, k]*M[i, j] - M[i, k]*M[k, j]

                    if k > 0:
                        D /= M[k-1, k-1]

                    if D.is_Atom:
                        M[i, j] = D
                    else:
                        M[i, j] = cancel(D)

        det = sign * M[n-1, n-1]

    return det.expand()

#############
# berkowitz #
#############

def berkowitz(self):
    """The Berkowitz algorithm.

        Given N x N matrix with symbolic content, compute efficiently
        coefficients of characteristic polynomials of 'self' and all
        its square sub-matrices composed by removing both i-th row
        and column, without division in the ground domain.

        This method is particularly useful for computing determinant,
        principal minors and characteristic polynomial, when 'self'
        has complicated coefficients e.g. polynomials. Semi-direct
        usage of this algorithm is also important in computing
        efficiently sub-resultant PRS.

        Assuming that M is a square matrix of dimension N x N and
        I is N x N identity matrix,  then the following following
        definition of characteristic polynomial is begin used:

                        charpoly(M) = det(t*I - M)

        As a consequence, all polynomials generated by Berkowitz
        algorithm are monic.

        >>> from sympy import Matrix
        >>> from sympy.abc import x, y, z

        >>> M = Matrix([ [x,y,z], [1,0,0], [y,z,x] ])

        >>> p, q, r = berkowitz(M)

        >>> print p # 1 x 1 M's sub-matrix
        (1, -x)

        >>> print q # 2 x 2 M's sub-matrix
        (1, -x, -y)

        >>> print r # 3 x 3 M's sub-matrix
        (1, -2*x, x**2 - y*z - y, x*y - z**2)

        For more information on the implemented algorithm refer to:

        [1] S.J. Berkowitz, On computing the determinant in small
            parallel time using a small number of processors, ACM,
            Information Processing Letters 18, 1984, pp. 147-150

        [2] M. Keber, Division-Free computation of sub-resultants
            using Bezout matrices, Tech. Report MPI-I-2006-1-006,
            Saarbrucken, 2006

    """
    if not self.is_square:
        raise NonSquareMatrixError()

    A, N = self, self.rows
    transforms = [0] * (N-1)

    for n in xrange(N, 1, -1):
        T, k = zeros((n+1,n)), n - 1

        R, C = -A[k,:k], A[:k,k]
        A, a = A[:k,:k], -A[k,k]

        items = [ C ]

        for i in xrange(0, n-2):
            items.append(A * items[i])

        for i, B in enumerate(items):
            items[i] = (R * B)[0,0]

        items = [ S.One, a ] + items

        for i in xrange(n):
            T[i:,i] = items[:n-i+1]

        transforms[k-1] = T

    polys = [ Matrix([S.One, -A[0,0]]) ]

    for i, T in enumerate(transforms):
        polys.append(T * polys[i])

    return tuple(map(tuple, polys))

def berkowitz_det(self):
    """Computes determinant using Berkowitz method."""
    poly = berkowitz(self)[-1]
    sign = (-1)**(len(poly)-1)
    return sign * poly[-1]

def berkowitz_minors(self):
    """Computes principal minors using Berkowitz method."""
    sign, minors = S.NegativeOne, []

    for poly in berkowitz(self):
        minors.append(sign*poly[-1])
        sign = -sign

    return tuple(minors)

def berkowitz_charpoly(self, x, simplify=sympy_simplify):
    """Computes characteristic polynomial minors using Berkowitz method."""
    return Poly(map(simplify, berkowitz(self)[-1]), x)

def berkowitz_eigenvals(self, **flags):
    """Computes eigenvalues of a Matrix using Berkowitz method. """
    return roots(berkowitz_charpoly(self, Dummy('x')), **flags)



