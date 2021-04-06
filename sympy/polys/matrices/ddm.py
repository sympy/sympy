"""

Module for the DDM class.

The DDM class is an internal representation used by DomainMatrix. The letters
DDM stand for Dense Domain Matrix. A DDM instance represents a matrix using
elements from a polynomial Domain (e.g. ZZ, QQ, ...) in a dense-matrix
representation.

Basic usage:

    >>> from sympy import ZZ, QQ
    >>> from sympy.polys.matrices.ddm import DDM
    >>> A = DDM([[ZZ(0), ZZ(1)], [ZZ(-1), ZZ(0)]], (2, 2), ZZ)
    >>> A.shape
    (2, 2)
    >>> A
    [[0, 1], [-1, 0]]
    >>> type(A)
    <class 'sympy.polys.matrices.ddm.DDM'>
    >>> A @ A
    [[-1, 0], [0, -1]]

The ddm_* functions are designed to operate on DDM as well as on an ordinary
list of lists:

    >>> from sympy.polys.matrices.dense import ddm_idet
    >>> ddm_idet(A, QQ)
    1
    >>> ddm_idet([[0, 1], [-1, 0]], QQ)
    1
    >>> A
    [[-1, 0], [0, 1]]

Note that ddm_idet modifies the input matrix in-place. It is recommended to
use the DDM.det method as a friendlier interface to this instead which takes
care of copying the matrix:

    >>> B = DDM([[ZZ(0), ZZ(1)], [ZZ(-1), ZZ(0)]], (2, 2), ZZ)
    >>> B.det()
    1

Normally DDM would not be used directly and is just part of the internal
representation of DomainMatrix which adds further functionality including e.g.
unifying domains.

The dense format used by DDM is a list of lists of elements e.g. the 2x2
identity matrix is like [[1, 0], [0, 1]]. The DDM class itself is a subclass
of list and its list items are plain lists. Elements are accessed as e.g.
ddm[i][j] where ddm[i] gives the ith row and ddm[i][j] gets the element in the
jth column of that row. Subclassing list makes e.g. iteration and indexing
very efficient. We do not override __getitem__ because it would lose that
benefit.

The core routines are implemented by the ddm_* functions defined in dense.py.
Those functions are intended to be able to operate on a raw list-of-lists
representation of matrices with most functions operating in-place. The DDM
class takes care of copying etc and also stores a Domain object associated
with its elements. This makes it possible to implement things like A + B with
domain checking and also shape checking so that the list of lists
representation is friendlier.

"""
from .exceptions import DDMBadInputError, DDMShapeError, DDMDomainError

from .dense import (
        ddm_iadd,
        ddm_isub,
        ddm_ineg,
        ddm_imul,
        ddm_imatmul,
        ddm_irref,
        ddm_idet,
        ddm_iinv,
        ddm_ilu_split,
        ddm_ilu_solve,
        ddm_berk,
        )


class DDM(list):
    """Dense matrix based on polys domain elements

    This is a list subclass and is a wrapper for a list of lists that supports
    basic matrix arithmetic +, -, *, **.
    """

    fmt = 'dense'

    def __init__(self, rowslist, shape, domain):
        super().__init__(rowslist)
        self.shape = self.rows, self.cols = m, n = shape
        self.domain = domain

        if not (len(self) == m and all(len(row) == n for row in self)):
            raise DDMBadInputError("Inconsistent row-list/shape")

    def to_list(self):
        return list(self)

    def to_ddm(self):
        return self

    def to_sdm(self):
        return SDM.from_list(self, self.shape, self.domain)

    def convert_to(self, K):
        Kold = self.domain
        if K == Kold:
            return self.copy()
        rows = ([K.convert_from(e, Kold) for e in row] for row in self)
        return DDM(rows, self.shape, K)

    def __str__(self):
        rowsstr = ['[%s]' % ', '.join(map(str, row)) for row in self]
        return '[%s]' % ', '.join(rowsstr)

    def __repr__(self):
        cls = type(self).__name__
        rows = list.__repr__(self)
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

    def __rmul__(a, b):
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

    def mul(a, b):
        c = a.copy()
        ddm_imul(c, b)
        return c

    def matmul(a, b):
        """a @ b (matrix product)"""
        m, o = a.shape
        o2, n = b.shape
        a._check(a, '*', b, o, o2)
        c = a.zeros((m, n), a.domain)
        ddm_imatmul(c, a, b)
        return c

    def hstack(A, B):
        Anew = list(A.copy())
        rows, cols = A.shape
        domain = A.domain

        Brows, Bcols = B.shape
        assert Brows == rows
        assert B.domain == domain

        cols += Bcols

        for i, Bi in enumerate(B):
            Anew[i].extend(Bi)

        return DDM(Anew, (rows, cols), A.domain)

    def rref(a):
        """Reduced-row echelon form of a and list of pivots"""
        b = a.copy()
        pivots = ddm_irref(b)
        return b, pivots

    def nullspace(a):
        rref, pivots = a.rref()
        rows, cols = a.shape
        domain = a.domain

        basis = []
        nonpivots = []
        for i in range(cols):
            if i in pivots:
                continue
            nonpivots.append(i)
            vec = [domain.one if i == j else domain.zero for j in range(cols)]
            for ii, jj in enumerate(pivots):
                vec[jj] -= rref[ii][i]
            basis.append(vec)

        return DDM(basis, (len(basis), cols), domain), nonpivots

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


from .sdm import SDM
