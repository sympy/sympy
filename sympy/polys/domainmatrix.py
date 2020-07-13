from __future__ import print_function

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
        if A.shape != B.shape:
            raise ShapeError("shape")
        if A.domain != B.domain:
            raise ValueError("domain")
        return A.from_ddm(A.rep + B.rep)

    def __sub__(A, B):
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        if A.shape != B.shape:
            raise ShapeError("shape")
        if A.domain != B.domain:
            raise ValueError("domain")
        return A.from_ddm(A.rep - B.rep)

    def __neg__(A):
        return A.from_ddm(-A.rep)

    def __mul__(A, B):
        """A * B"""
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        return A.from_ddm(A.rep @ B.rep)

    def __pow__(A, n):
        """A ** n"""
        if not isinstance(n, int):
            return NotImplemented
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
        dom = self.domain
        rows = self.rep
        eye = [[dom.one if i==j else dom.zero for j in range(n)] for i in range(n)]
        rows = [row + eyerow for row, eyerow in zip(rows, eye)]
        Aaug = DomainMatrix(rows, (n, 2*n), dom)
        Aaug_rref, pivots = Aaug.rref()
        if pivots != tuple(range(n)):
            raise NonInvertibleMatrixError('Matrix det == 0; not invertible.')
        Ainv_rows = [row[n:] for row in Aaug_rref.rep]
        Ainv = DomainMatrix(Ainv_rows, (n, n), dom)
        return Ainv

    def det(self):
        m, n = self.shape
        if m != n:
            raise NonSquareMatrixError
        return self.rep.det()

    def __eq__(A, B):
        """A == B"""
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        return A.rep == B.rep


def lu_decomp(M):

    dom = M.domain

    if not dom.is_Field:
        raise ValueError("Domain")

    N, O = M.shape

    lu = [row[:] for row in M.rep]
    swaps = []

    for n in range(min(N, O)):
        if not lu[n][n]:
            for m in range(n+1, N):
                if lu[m][n]:
                    swaps.append((n, m))
                    lu[n], lu[m] = lu[m], lu[n]
                    break
            else:
                # M = Matrix([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 1, 2]])
                continue
        for i in range(n+1, N):
            l_in = lu[i][n] / lu[n][n]
            lu[i][n] = l_in
            for j in range(n+1, O):
                lu[i][j] -= l_in * lu[n][j]

    L = []
    U = []
    for i in range(N):
        if i < O:
            L.append(lu[i][:i] + [dom.one] + [dom.zero] * (N-i-1))
            U.append([dom.zero] * i + lu[i][i:])
        else:
            L.append(lu[i] + [dom.zero] * (i-O) + [dom.one] + [dom.zero] * (N-i-1))
            U.append([dom.zero] * O)

    L = DomainMatrix(L, (N, N), dom)
    U = DomainMatrix(U, (N, O), dom)

    return L, U, swaps


def lu_solve(M, b):
    if M.domain != b.domain:
        raise ValueError("Domain")

    dom = M.domain
    m, n = M.shape
    m2, o = b.shape

    if m != m2:
        raise ShapeError("Shape")
    if m < n:
        raise NotImplementedError("Underdetermined")

    L, U, swaps = lu_decomp(M)

    b_rows = [row[:] for row in b.rep]
    if swaps:
        for i1, i2 in swaps:
            b_rows[i1], b_rows[i2] = b_rows[i2], b_rows[i1]

    # solve Ly = b
    y = [[None] * o for _ in range(m)]
    for k in range(o):
        for i in range(m):
            rhs = b_rows[i][k]
            for j in range(i):
                rhs -= L.rep[i][j] * y[j][k]
            y[i][k] = rhs

    if m > n:
        for i in range(n, m):
            for j in range(o):
                if y[i][j]:
                    raise NonInvertibleMatrixError

    # Solve Ux = y
    x = [[None] * o for _ in range(n)]
    for k in range(o):
        for i in reversed(range(n)):
            if not U.rep[i][i]:
                raise NonInvertibleMatrixError
            rhs = y[i][k]
            for j in range(i+1, n):
                rhs -= U.rep[i][j] * x[j][k]
            x[i][k] = rhs / U.rep[i][i]

    x = DomainMatrix(x, (n, o), dom)
    return x

def berk(M):
    dom = M.domain
    m, n = M.shape
    assert m == n

    if n == 1:
        v = [[1], [-M.rep[0][0]]]
        return DomainMatrix(v, (2, 1), dom)

    a = M.rep[0][0]
    R = DomainMatrix([M.rep[0][1:]], (1, n-1), dom)
    C = DomainMatrix([[row[0]] for row in M.rep[1:]], (n-1, 1), dom)
    A = DomainMatrix([row[1:] for row in M.rep[1:]], (n-1, n-1), dom)

    q = berk(A)

    T = [[dom.zero] * n for _ in range(n+1)]
    for i in range(n):
        T[i][i] = dom.one
        T[i+1][i] = -a
    for i in range(2, n+1):
        if i == 2:
            AnC = C
        else:
            AnC = A * AnC
        RAnC = R * AnC
        assert RAnC.shape == (1, 1)
        for j in range(0, n+1-i):
            T[i+j][j] = -RAnC.rep[0][0]
    T = DomainMatrix(T, (n+1, n), dom)

    return T * q
