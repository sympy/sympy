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

    d = self.berkowitz_det()
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

def eye():
def zeros():

