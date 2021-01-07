from sympy.core.sympify import _sympify

from ..constructor import construct_domain

from .exceptions import NonSquareMatrixError, ShapeError

from .ddm import DDM


class DomainMatrix:

    def __init__(self, rows, shape, domain):
        self.rep = DDM(rows, shape, domain)
        self.shape = shape
        self.domain = domain

    @classmethod
    def from_ddm(cls, ddm):
        return cls(ddm, ddm.shape, ddm.domain)

    @classmethod
    def from_list_sympy(cls, nrows, ncols, rows, **kwargs):
        assert len(rows) == nrows
        assert all(len(row) == ncols for row in rows)

        items_sympy = [_sympify(item) for row in rows for item in row]

        domain, items_domain = cls.get_domain(items_sympy, **kwargs)

        domain_rows = [[items_domain[ncols*r + c] for c in range(ncols)] for r in range(nrows)]

        return DomainMatrix(domain_rows, (nrows, ncols), domain)

    @classmethod
    def from_Matrix(cls, M, **kwargs):
        return cls.from_list_sympy(*M.shape, M.tolist(), **kwargs)

    @classmethod
    def get_domain(cls, items_sympy, **kwargs):
        K, items_K = construct_domain(items_sympy, **kwargs)
        return K, items_K

    def convert_to(self, K):
        Kold = self.domain
        if K == Kold:
            return self.from_ddm(self.rep.copy())
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
        if isinstance(B, DomainMatrix):
            return A.matmul(B)
        elif B in A.domain:
            return A.from_ddm(A.rep * B)
        else:
            return NotImplemented

    def __rmul__(A, B):
        if B in A.domain:
            return A.from_ddm(A.rep * B)
        else:
            return NotImplemented

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

    def mul(A, b):
        return A.from_ddm(A.rep.mul(b))

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

    def nullspace(self):
        return self.from_ddm(self.rep.nullspace())

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

    @classmethod
    def eye(cls, n, domain):
        return cls.from_ddm(DDM.eye(n, domain))

    def __eq__(A, B):
        """A == B"""
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        return A.rep == B.rep
