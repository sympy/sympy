"""

Module for the DomainMatrix class.

A DomainMatrix represents a matrix with elements that are in a particular
Domain. Each DomainMatrix internally wraps a DDM which is used for the
lower-level operations. The idea is that the DomainMatrix class provides the
convenience routines for converting between Expr and the poly domains as well
as unifying matrices with different domains.

"""
from sympy.core.sympify import _sympify

from ..constructor import construct_domain

from .exceptions import NonSquareMatrixError, ShapeError

from .ddm import DDM

from .sdm import SDM


class DomainMatrix:
    r"""

        Module for the DomainMatrix class.

        The idea behind DomainMatrix is to associate "domain" with matrix,
        which could be found analogous to numpy arrays with "dtype".
        In the DomainMatrix, each matrix has a domain such as ZZ(Integers) or QQ<sqrt(2)>[Algebraic Number Field].

        DomainMatrix uses "domain elements" for its internal representation which makes it more faster for many common operations
        than current sympy Matrix class, but this advantage makes it not entirely compatible with Matrix.


        Creating a DomainMatrix from the existing Matrix class:

        >>> from sympy import Matrix
        >>> from sympy.polys.matrices import DomainMatrix
        >>> matrix = Matrix([[99, 16, 30, 62], [65, 79, 94, 15], [38, 13, 31, 47], [17, 2, 53, 78]])
        >>> domainMatrix = DomainMatrix.from_Matrix(matrix)
        >>> domainMatrix
        DomainMatrix([[99, 16, 30, 62], [65, 79, 94, 15], [38, 13, 31, 47], [17, 2, 53, 78]], (4, 4), ZZ)


        Methods:
        ========

        >>> from sympy import RR, CC

        @convert_to: Returns a DomainMatrix with the desired domain or field
        >>> domainMatrix.convert_to(CC)
        DomainMatrix([[(99.0 + 0.0j), (16.0 + 0.0j), (30.0 + 0.0j), (62.0 + 0.0j)], [(65.0 + 0.0j), (79.0 + 0.0j), (94.0 + 0.0j), (15.0 + 0.0j)],
        [(38.0 + 0.0j), (13.0 + 0.0j), (31.0 + 0.0j), (47.0 + 0.0j)], [(17.0 + 0.0j), (2.0 + 0.0j), (53.0 + 0.0j), (78.0 + 0.0j)]], (4, 4), CC)

        @to_field: Returns a DomainMatrix with the appropriate field
        >>> domainMatrix.to_field()
        DomainMatrix([[99, 16, 30, 62], [65, 79, 94, 15], [38, 13, 31, 47], [17, 2, 53, 78]], (4, 4), QQ)

        @to_Matrix: Returns a MutableDenseMatrix for a DomainMatrix. This method can be used to make changes to a DomainMatrix.
        >>> mutableDomainMatrix = domainMatrix.to_Matrix()
        >>> mutableDomainMatrix
        Matrix([
        [99, 16, 30, 62],
        [65, 79, 94, 15],
        [38, 13, 31, 47],
        [17,  2, 53, 78]])
        >>> mutableDomainMatrix[0]
        99
        >>> mutableDomainMatrix[0] = 56
        >>> mutableDomainMatrix
        Matrix([
        [56, 16, 30, 62],
        [65, 79, 94, 15],
        [38, 13, 31, 47],
        [17,  2, 53, 78]])

        @det: Returns the determinant value for a Square DomainMatrix
        >>> domainMatrix.det()
        -3026882

        @charpoly: Returns the coefficients of the characteristic polynomial of the DomainMatrix
        >>> domainMatrix.charpoly()
        [1, -287, 22664, -473672, -3026882]


        >>> Matrix1 = Matrix([[97, 36, 48, 61], [7, 38, 9, 69], [4, 29, 63, 21], [44, 84, 40, 19]])
        >>> domainMatrix1 = DomainMatrix.from_Matrix(Matrix1)
        >>> domainMatrix1
        DomainMatrix([[97, 36, 48, 61], [7, 38, 9, 69], [4, 29, 63, 21], [44, 84, 40, 19]], (4, 4), ZZ)

        @unify: Unifies any two DomainMatrix with different domains
        >>> domainMatrix.unify(domainMatrix1)
        (DomainMatrix([[99, 16, 30, 62], [65, 79, 94, 15], [38, 13, 31, 47], [17, 2, 53, 78]], (4, 4), ZZ),
        DomainMatrix([[97, 36, 48, 61], [7, 38, 9, 69], [4, 29, 63, 21], [44, 84, 40, 19]], (4, 4), ZZ))

        @__repr__: Represents the DomainMatrix
        >>> domainMatrix.__repr__()
        'DomainMatrix([[99, 16, 30, 62], [65, 79, 94, 15], [38, 13, 31, 47], [17, 2, 53, 78]], (4, 4), ZZ)'

        @add: Adds two DomainMatrix matrices
        >>> domainMatrix.add(domainMatrix1)
        DomainMatrix([[196, 52, 78, 123], [72, 117, 103, 84], [42, 42, 94, 68], [61, 86, 93, 97]], (4, 4), ZZ)

        @sub: Subtracts two DomainMatrix matrices
        >>> domainMatrix.sub(domainMatrix1)
        DomainMatrix([[2, -20, -18, 1], [58, 41, 85, -54], [34, -16, -32, 26], [-27, -82, 13, 59]], (4, 4), ZZ)

        @neg: Returns the negative of DomainMatrix
        >>> domainMatrix.neg()
        DomainMatrix([[-99, -16, -30, -62], [-65, -79, -94, -15], [-38, -13, -31, -47], [-17, -2, -53, -78]], (4, 4), ZZ)

        @__eq__: Checks for two DomainMatrix matrices two be equal or not
        >>> domainMatrix.__eq__(domainMatrix)
        True
        >>> domainMatrix.__eq__(domainMatrix1)
        False

        @rref: Returns reduced-row echelon form and list of pivots for the DomainMatrix
        >>> dMF = domainMatrix.convert_to(RR)
        >>> dMF
        DomainMatrix([[99.0, 16.0, 30.0, 62.0], [65.0, 79.0, 94.0, 15.0], [38.0, 13.0, 31.0, 47.0], [17.0, 2.0, 53.0, 78.0]], (4, 4), RR)
        >>> dMF.rref()
        (DomainMatrix([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]], (4, 4), RR),
        (0, 1, 2, 3))

        @nullspace: Returns the Null Space for the DomainMatrix
        >>> dMF.nullspace()
        DomainMatrix([], (0, 4), RR)

        @lu: Returns Lower and Upper decomposition(L, U triangular matrices) of the DomainMatrix
        >>> dMF.lu()
        (DomainMatrix([[1.0, 0.0, 0.0, 0.0], [0.656565656565657, 1.0, 0.0, 0.0], [0.383838383838384, 0.100132723787052, 1.0, 0.0], [0.171717171717172, -0.0109128447131691, 4.03990205081114, 1.0]], (4, 4), RR),
        DomainMatrix([[99.0, 16.0, 30.0, 62.0], [0.0, 68.4949494949495, 74.3030303030303, -25.7070707070707], [0.0, 0.0, 12.0446836749742, 25.7761392125055], [0.0, 0.0, 0.0, -37.060079583716]], (4, 4), RR), [])


        @lu_solve: Solver for "DomainMatrix" x in the A*x = B, where A, B are DomainMatrix matrices.
        >>> A = Matrix([[77, 97], [76, 31]])
        >>> B = Matrix([[54, 30], [72, 77]])
        >>> A1 = DomainMatrix.from_Matrix(A)
        >>> B1 = DomainMatrix.from_Matrix(B)
        >>> A1
        DomainMatrix([[77, 97], [76, 31]], (2, 2), ZZ)
        >>> B1
        DomainMatrix([[54, 30], [72, 77]], (2, 2), ZZ)
        >>> A1f = A1.convert_to(RR)
        >>> B1f = B1.convert_to(RR)
        >>> A1f.lu_solve(B1f)
        DomainMatrix([[1.06519558676028, 1.31173520561685], [-0.288866599799398, -0.731995987963892]], (2, 2), RR)


        >>> Matrix1 = Matrix([[2, 73], [27, 86]])
        >>> domainMatrix = DomainMatrix.from_Matrix(Matrix1)
        >>> Matrix2 = Matrix([[19, 30], [56, 14]])
        >>> domainMatrix1 = DomainMatrix.from_Matrix(Matrix2)
        >>> domainMatrix
        DomainMatrix([[2, 73], [27, 86]], (2, 2), ZZ)
        >>> domainMatrix1
        DomainMatrix([[19, 30], [56, 14]], (2, 2), ZZ)

        @mul: Performs term by term multiplication for the second DomainMatrix w.r.t first DomainMatrix.
        Returns a DomainMatrix whose rows are list of DomainMatrix matrices created after term by term multiplication.
        >>> domainMatrix.mul(domainMatrix1)
        DomainMatrix([[DomainMatrix([[38, 60], [112, 28]], (2, 2), ZZ), DomainMatrix([[1387, 2190], [4088, 1022]], (2, 2), ZZ)],
        [DomainMatrix([[513, 810], [1512, 378]], (2, 2), ZZ), DomainMatrix([[1634, 2580], [4816, 1204]], (2, 2), ZZ)]], (2, 2), ZZ)

        @matmul: Performs matrix multiplication of two DomainMatrix matrices
        >>> domainMatrix.matmul(domainMatrix1)
        DomainMatrix([[4126, 1082], [5329, 2014]], (2, 2), ZZ)

        @pow: Perfoms multiplication on the given DomainMatrix for the given number of times, or else computes A^n
        >>> domainMatrix.pow(2)
        DomainMatrix([[1975, 6424], [2376, 9367]], (2, 2), ZZ)

        @inv: Finds the inverse of the DomainMatrix if exists
        >>> dMF = domainMatrix.convert_to(RR)
        >>> dMF
        DomainMatrix([[2.0, 73.0], [27.0, 86.0]], (2, 2), RR)
        >>> dMF.inv()
        DomainMatrix([[-0.0478043357420789, 0.0405780989438577], [0.0150083379655364, -0.00111172873818788]], (2, 2), RR)


        References
        ==========

        ..[1] https://docs.sympy.org/dev/modules/polys/domainsintro.html

    """

    def __init__(self, rows, shape, domain):
        if isinstance(rows, list):
            self.rep = SDM.from_list(rows, shape, domain)
        else:
            self.rep = SDM(rows, shape, domain)
        self.shape = shape
        self.domain = domain

    @classmethod
    def from_rep(cls, ddm):
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
        return self.from_rep(self.rep.convert_to(K))

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
        elemlist = self.rep.to_list()
        rows_sympy = [[self.domain.to_sympy(e) for e in row] for row in elemlist]
        return MutableDenseMatrix(rows_sympy)

    def __repr__(self):
        elemlist = self.rep.to_list()
        rows_str = ['[%s]' % (', '.join(map(str, row))) for row in elemlist]
        rowstr = '[%s]' % ', '.join(rows_str)
        return 'DomainMatrix(%s, %r, %r)' % (rowstr, self.shape, self.domain)

    def hstack(A, B):
        A, B = A.unify(B)
        return A.from_rep(A.rep.hstack(B.rep))

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
            return A.from_rep(A.rep * B)
        else:
            return NotImplemented

    def __rmul__(A, B):
        if B in A.domain:
            return A.from_rep(A.rep * B)
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
        return A.from_rep(A.rep.add(B.rep))

    def sub(A, B):
        if A.shape != B.shape:
            raise ShapeError("shape")
        if A.domain != B.domain:
            raise ValueError("domain")
        return A.from_rep(A.rep.sub(B.rep))

    def neg(A):
        return A.from_rep(A.rep.neg())

    def mul(A, b):
        return A.from_rep(A.rep.mul(b))

    def matmul(A, B):
        return A.from_rep(A.rep.matmul(B.rep))

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
        return self.from_rep(rref_ddm), tuple(pivots)

    def nullspace(self):
        return self.from_rep(self.rep.nullspace()[0])

    def inv(self):
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        m, n = self.shape
        if m != n:
            raise NonSquareMatrixError
        inv = self.rep.inv()
        return self.from_rep(inv)

    def det(self):
        m, n = self.shape
        if m != n:
            raise NonSquareMatrixError
        return self.rep.det()

    def lu(self):
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        L, U, swaps = self.rep.lu()
        return self.from_rep(L), self.from_rep(U), swaps

    def lu_solve(self, rhs):
        if self.shape[0] != rhs.shape[0]:
            raise ShapeError("Shape")
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        sol = self.rep.lu_solve(rhs.rep)
        return self.from_rep(sol)

    def charpoly(self):
        m, n = self.shape
        if m != n:
            raise NonSquareMatrixError("not square")
        return self.rep.charpoly()

    @classmethod
    def eye(cls, n, domain):
        return cls.from_rep(DDM.eye(n, domain))

    def __eq__(A, B):
        """A == B"""
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        return A.rep == B.rep
