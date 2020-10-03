from operator import mul

from sympy.core.sympify import _sympify
from sympy.matrices.common import (NonInvertibleMatrixError,
    NonSquareMatrixError, ShapeError)
from sympy.polys.constructor import construct_domain


class DDMError(Exception):
    """Base class for errors raised by DDM"""
    pass


class DDMBadInputError(DDMError):
    """list of lists is inconsistent with shape"""
    pass


class DDMDomainError(DDMError):
    """domains do not match"""
    pass


class DDMShapeError(DDMError):
    """shapes are inconsistent"""
    pass


class DDM(list):
    """Dense matrix based on polys domain elements

    This is a list subclass and is a wrapper for a list of lists that supports
    basic matrix arithmetic +, -, *, **.
    """
    def __init__(self, rowslist, shape, domain):
        super().__init__(rowslist)
        self.shape = self.rows, self.cols = m, n = shape
        self.domain = domain

        if not (len(self) == m and all(len(row) == n for row in self)):
            raise DDMBadInputError("Inconsistent row-list/shape")

    def __str__(self):
        cls = type(self).__name__
        rows = list.__str__(self)
        return '%s(%s, %s, %s)' % (cls, rows, self.shape, self.domain)

    def __eq__(self, other):
        if not isinstance(other, DDM):
            return False
        return (super().__eq__(other) and self.domain == other.domain)

    def __ne__(self, other):
        return not self.__eq__(other)

    @classmethod
    def zeros(cls, shape, domain):
        z = domain.zero
        m, n = shape
        rowslist = ([z] * n for _ in range(m))
        return DDM(rowslist, shape, domain)

    @classmethod
    def eye(cls, size, domain):
        one = domain.one
        ddm = cls.zeros((size, size), domain)
        for i in range(size):
            ddm[i][i] = one
        return ddm

    def copy(self):
        copyrows = (row[:] for row in self)
        return DDM(copyrows, self.shape, self.domain)

    def __add__(a, b):
        if not isinstance(b, DDM):
            return NotImplemented
        return a.add(b)

    def __sub__(a, b):
        if not isinstance(b, DDM):
            return NotImplemented
        return a.sub(b)

    def __neg__(a):
        return a.neg()

    def __mul__(a, b):
        if b in a.domain:
            return a.mul(b)
        else:
            return NotImplemented

    def __matmul__(a, b):
        if isinstance(b, DDM):
            return a.matmul(b)
        else:
            return NotImplemented

    @classmethod
    def _check(cls, a, op, b, ashape, bshape):
        if a.domain != b.domain:
            msg = "Domain mismatch: %s %s %s" % (a.domain, op, b.domain)
            raise DDMDomainError(msg)
        if ashape != bshape:
            msg = "Shape mismatch: %s %s %s" % (a.shape, op, b.shape)
            raise DDMShapeError(msg)

    def add(a, b):
        """a + b"""
        a._check(a, '+', b, a.shape, b.shape)
        c = a.copy()
        ddm_iadd(c, b)
        return c

    def sub(a, b):
        """a - b"""
        a._check(a, '-', b, a.shape, b.shape)
        c = a.copy()
        ddm_isub(c, b)
        return c

    def neg(a):
        """-a"""
        b = a.copy()
        ddm_ineg(b)
        return b

    def matmul(a, b):
        """a @ b (matrix product)"""
        m, o = a.shape
        o2, n = b.shape
        a._check(a, '*', b, o, o2)
        c = a.zeros((m, n), a.domain)
        ddm_imatmul(c, a, b)
        return c

    def rref(a):
        """Reduced-row echelon form of a and list of pivots"""
        b = a.copy()
        pivots = ddm_irref(b)
        return b, pivots

    def det(a):
        """Determinant of a"""
        m, n = a.shape
        if m != n:
            raise DDMShapeError("Determinant of non-square matrix")
        b = a.copy()
        K = b.domain
        deta = ddm_idet(b, K)
        return deta

    def inv(a):
        """Inverse of a"""
        m, n = a.shape
        if m != n:
            raise DDMShapeError("Determinant of non-square matrix")
        ainv = a.copy()
        K = a.domain
        ddm_iinv(ainv, a, K)
        return ainv

    def lu(a):
        """L, U decomposition of a"""
        m, n = a.shape
        K = a.domain

        U = a.copy()
        L = a.eye(m, K)
        swaps = ddm_ilu_split(L, U, K)

        return L, U, swaps

    def lu_solve(a, b):
        """x where a*x = b"""
        m, n = a.shape
        m2, o = b.shape
        a._check(a, 'lu_solve', b, m, m2)

        L, U, swaps = a.lu()
        x = a.zeros((n, o), a.domain)
        ddm_ilu_solve(x, L, U, swaps, b)
        return x

    def charpoly(a):
        """Coefficients of characteristic polynomial of a"""
        K = a.domain
        m, n = a.shape
        if m != n:
            raise DDMShapeError("Charpoly of non-square matrix")
        vec = ddm_berk(a, K)
        coeffs = [vec[i][0] for i in range(n+1)]
        return coeffs


def ddm_iadd(a, b):
    """a += b"""
    for ai, bi in zip(a, b):
        for j, bij in enumerate(bi):
            ai[j] += bij


def ddm_isub(a, b):
    """a -= b"""
    for ai, bi in zip(a, b):
        for j, bij in enumerate(bi):
            ai[j] -= bij


def ddm_ineg(a):
    """a  <--  -a"""
    for ai in a:
        for j, aij in enumerate(ai):
            ai[j] = -aij


def ddm_imatmul(a, b, c):
    """a += b @ c"""
    cT = list(zip(*c))

    for bi, ai in zip(b, a):
        for j, cTj in enumerate(cT):
            ai[j] = sum(map(mul, bi, cTj), ai[j])


def ddm_irref(a):
    """a  <--  rref(a)"""
    # a is (m x n)
    m = len(a)
    if not m:
        return []
    n = len(a[0])

    i = 0
    pivots = []

    for j in range(n):
        # pivot
        aij = a[i][j]

        # zero-pivot
        if not aij:
            for ip in range(i+1, m):
                aij = a[ip][j]
                # row-swap
                if aij:
                    a[i], a[ip] = a[ip], a[i]
                    break
            else:
                # next column
                continue

        # normalise row
        ai = a[i]
        for l in range(j, n):
            ai[l] /= aij # ai[j] = one

        # eliminate above and below to the right
        for k, ak in enumerate(a):
            if k == i or not ak[j]:
                continue
            akj = ak[j]
            ak[j] -= akj # ak[j] = zero
            for l in range(j+1, n):
                ak[l] -= akj * ai[l]

        # next row
        pivots.append(j)
        i += 1

        # no more rows?
        if i >= m:
            break

    return pivots


def ddm_idet(a, K):
    """a  <--  echelon(a); return det"""
    # Fraction-free Gaussian elimination
    # https://www.math.usm.edu/perry/Research/Thesis_DRL.pdf

    # a is (m x n)
    m = len(a)
    if not m:
        return K.one
    n = len(a[0])

    is_field = K.is_Field
    # uf keeps track of the effect of row swaps and multiplies
    uf = K.one
    for j in range(n-1):
        # if zero on the diagonal need to swap
        if not a[j][j]:
            for l in range(j+1, n):
                if a[l][j]:
                    a[j], a[l] = a[l], a[j]
                    uf = -uf
                    break
            else:
                # unable to swap: det = 0
                return K.zero
        for i in range(j+1, n):
            if a[i][j]:
                if not is_field:
                    d = K.gcd(a[j][j], a[i][j])
                    b = a[j][j] // d
                    c = a[i][j] // d
                else:
                    b = a[j][j]
                    c = a[i][j]
                # account for multiplying row i by b
                uf = b * uf
                for k in range(j+1, n):
                    a[i][k] = b*a[i][k] - c*a[j][k]

    # triangular det is product of diagonal
    prod = K.one
    for i in range(n):
        prod = prod * a[i][i]
    # incorporate swaps and multiplies
    if not is_field:
        D = prod // uf
    else:
        D = prod / uf
    return D


def ddm_iinv(ainv, a, K):
    if not K.is_Field:
        raise ValueError('Not a field')

    # a is (m x n)
    m = len(a)
    if not m:
        return
    n = len(a[0])
    if m != n:
        raise NonSquareMatrixError

    eye = [[K.one if i==j else K.zero for j in range(n)] for i in range(n)]
    Aaug = [row + eyerow for row, eyerow in zip(a, eye)]
    pivots = ddm_irref(Aaug)
    if pivots != list(range(n)):
        raise NonInvertibleMatrixError('Matrix det == 0; not invertible.')
    ainv[:] = [row[n:] for row in Aaug]


def ddm_ilu_split(L, U, K):
    """L, U  <--  LU(U)"""
    m = len(U)
    if not m:
        return []
    n = len(U[0])

    swaps = ddm_ilu(U)

    zeros = [K.zero] * min(m, n)
    for i in range(1, m):
        j = min(i, n)
        L[i][:j] = U[i][:j]
        U[i][:j] = zeros[:j]

    return swaps


def ddm_ilu(a):
    """a  <--  LU(a)"""
    m = len(a)
    if not m:
        return []
    n = len(a[0])

    swaps = []

    for i in range(min(m, n)):
        if not a[i][i]:
            for ip in range(i+1, m):
                if a[ip][i]:
                    swaps.append((i, ip))
                    a[i], a[ip] = a[ip], a[i]
                    break
            else:
                # M = Matrix([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 1, 2]])
                continue
        for j in range(i+1, m):
            l_ji = a[j][i] / a[i][i]
            a[j][i] = l_ji
            for k in range(i+1, n):
                a[j][k] -= l_ji * a[i][k]

    return swaps


def ddm_ilu_solve(x, L, U, swaps, b):
    """x  <--  solve(L*U*x = swaps(b))"""
    m = len(U)
    if not m:
        return
    n = len(U[0])

    m2 = len(b)
    if not m2:
        raise DDMShapeError("Shape mismtch")
    o = len(b[0])

    if m != m2:
        raise DDMShapeError("Shape mismtch")
    if m < n:
        raise NotImplementedError("Underdetermined")

    if swaps:
        b = [row[:] for row in b]
        for i1, i2 in swaps:
            b[i1], b[i2] = b[i2], b[i1]

    # solve Ly = b
    y = [[None] * o for _ in range(m)]
    for k in range(o):
        for i in range(m):
            rhs = b[i][k]
            for j in range(i):
                rhs -= L[i][j] * y[j][k]
            y[i][k] = rhs

    if m > n:
        for i in range(n, m):
            for j in range(o):
                if y[i][j]:
                    raise NonInvertibleMatrixError

    # Solve Ux = y
    for k in range(o):
        for i in reversed(range(n)):
            if not U[i][i]:
                raise NonInvertibleMatrixError
            rhs = y[i][k]
            for j in range(i+1, n):
                rhs -= U[i][j] * x[j][k]
            x[i][k] = rhs / U[i][i]


def ddm_berk(M, K):
    m = len(M)
    if not m:
        return [[K.one]]
    n = len(M[0])

    if m != n:
        raise DDMShapeError("Not square")

    if n == 1:
        return [[K.one], [-M[0][0]]]

    a = M[0][0]
    R = [M[0][1:]]
    C = [[row[0]] for row in M[1:]]
    A = [row[1:] for row in M[1:]]

    q = ddm_berk(A, K)

    T = [[K.zero] * n for _ in range(n+1)]
    for i in range(n):
        T[i][i] = K.one
        T[i+1][i] = -a
    for i in range(2, n+1):
        if i == 2:
            AnC = C
        else:
            C = AnC
            AnC = [[K.zero] for row in C]
            ddm_imatmul(AnC, A, C)
        RAnC = [[K.zero]]
        ddm_imatmul(RAnC, R, AnC)
        for j in range(0, n+1-i):
            T[i+j][j] = -RAnC[0][0]

    qout = [[K.zero] for _ in range(n+1)]
    ddm_imatmul(qout, T, q)
    return qout


class DomainMatrix:

    def __init__(self, rows, shape, domain):
        self.rep = DDM(rows, shape, domain)
        self.shape = shape
        self.domain = domain

    @classmethod
    def from_ddm(cls, ddm):
        return cls(ddm, ddm.shape, ddm.domain)

    @classmethod
    def from_list_sympy(cls, nrows, ncols, rows):
        assert len(rows) == nrows
        assert all(len(row) == ncols for row in rows)

        items_sympy = [_sympify(item) for row in rows for item in row]

        domain, items_domain = cls.get_domain(items_sympy)

        domain_rows = [[items_domain[ncols*r + c] for c in range(ncols)] for r in range(nrows)]

        return DomainMatrix(domain_rows, (nrows, ncols), domain)

    @classmethod
    def get_domain(cls, items_sympy, **kwargs):
        K, items_K = construct_domain(items_sympy, **kwargs)
        return K, items_K

    def convert_to(self, K):
        Kold = self.domain
        new_rows = [[K.convert_from(e, Kold) for e in row] for row in self.rep]
        return DomainMatrix(new_rows, self.shape, K)

    def to_field(self):
        K = self.domain.get_field()
        return self.convert_to(K)

    def unify(self, other):
        K1 = self.domain
        K2 = other.domain
        if K1 == K2:
            return self, other
        K = K1.unify(K2)
        if K1 != K:
            self = self.convert_to(K)
        if K2 != K:
            other = other.convert_to(K)
        return self, other

    def to_Matrix(self):
        from sympy.matrices.dense import MutableDenseMatrix
        rows_sympy = [[self.domain.to_sympy(e) for e in row] for row in self.rep]
        return MutableDenseMatrix(rows_sympy)

    def __repr__(self):
        rows_str = ['[%s]' % (', '.join(map(str, row))) for row in self.rep]
        rowstr = '[%s]' % ', '.join(rows_str)
        return 'DomainMatrix(%s, %r, %r)' % (rowstr, self.shape, self.domain)

    def __add__(A, B):
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        return A.add(B)

    def __sub__(A, B):
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        return A.sub(B)

    def __neg__(A):
        return A.neg()

    def __mul__(A, B):
        """A * B"""
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        return A.matmul(B)

    def __pow__(A, n):
        """A ** n"""
        if not isinstance(n, int):
            return NotImplemented
        return A.pow(n)

    def add(A, B):
        if A.shape != B.shape:
            raise ShapeError("shape")
        if A.domain != B.domain:
            raise ValueError("domain")
        return A.from_ddm(A.rep.add(B.rep))

    def sub(A, B):
        if A.shape != B.shape:
            raise ShapeError("shape")
        if A.domain != B.domain:
            raise ValueError("domain")
        return A.from_ddm(A.rep.sub(B.rep))

    def neg(A):
        return A.from_ddm(A.rep.neg())

    def matmul(A, B):
        return A.from_ddm(A.rep.matmul(B.rep))

    def pow(A, n):
        if n < 0:
            raise NotImplementedError('Negative powers')
        elif n == 0:
            m, n = A.shape
            rows = [[A.domain.zero] * m for _ in range(m)]
            for i in range(m):
                rows[i][i] = A.domain.one
            return type(A)(rows, A.shape, A.domain)
        elif n == 1:
            return A
        elif n % 2 == 1:
            return A * A**(n - 1)
        else:
            sqrtAn = A ** (n // 2)
            return sqrtAn * sqrtAn

    def rref(self):
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        rref_ddm, pivots = self.rep.rref()
        return self.from_ddm(rref_ddm), tuple(pivots)

    def inv(self):
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        m, n = self.shape
        if m != n:
            raise NonSquareMatrixError
        inv = self.rep.inv()
        return self.from_ddm(inv)

    def det(self):
        m, n = self.shape
        if m != n:
            raise NonSquareMatrixError
        return self.rep.det()

    def lu(self):
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        L, U, swaps = self.rep.lu()
        return self.from_ddm(L), self.from_ddm(U), swaps

    def lu_solve(self, rhs):
        if self.shape[0] != rhs.shape[0]:
            raise ShapeError("Shape")
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        sol = self.rep.lu_solve(rhs.rep)
        return self.from_ddm(sol)

    def charpoly(self):
        m, n = self.shape
        if m != n:
            raise NonSquareMatrixError("not square")
        return self.rep.charpoly()

    def __eq__(A, B):
        """A == B"""
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        return A.rep == B.rep
