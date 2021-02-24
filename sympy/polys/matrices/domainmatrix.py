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
    Associate Matrix with :py:class:`~.Domain`

    Explanation
    ===========

    DomainMatrix uses "domain elements" for its internal representation
    which makes it more faster for many common operations
    than current sympy Matrix class, but this advantage makes it not
    entirely compatible with Matrix.
    DomainMatrix could be found analogous to numpy arrays with "dtype".
    In the DomainMatrix, each matrix has a domain such as ZZ(Integers)
    or QQ<sqrt(2)>[Algebraic Number Field].


    Examples
    ========

    Creating a DomainMatrix from the existing Matrix class:

    >>> from sympy import Matrix
    >>> from sympy.polys.matrices import DomainMatrix
    >>> Matrix1 = Matrix([
    ...    [1, 2],
    ...    [3, 4]])
    >>> domainMatrix = DomainMatrix.from_Matrix(Matrix1)
    >>> domainMatrix
    DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ)

    Driectly forming a DomainMatrix:

    >>> from sympy import ZZ
    >>> from sympy.polys.matrices import DomainMatrix
    >>> domainMatrix = DomainMatrix([
    ...    [ZZ(1), ZZ(2)],
    ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    >>> domainMatrix
    DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ)

    See Also
    ========

    :py:class:`~.ddm`
    :py:class:`~.Domain`
    :py:class:`~.Poly`

    """

    def __init__(self, rows, shape, domain):
        """
        Creates a :py:class:~.DomainMatrix.

        Parameters
        ==========

        rows : Represents elements of DomainMatrix as list of lists
        shape : Represents dimension of DomainMatrix
        domain : Represents :py:class:`~.Domain` of DomainMatrix

        Raises
        ======

        TypeError
            If any of rows, shape and domain are not provided

        """
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
        r"""
        Returns a DomainMatrix with the desired domain or field

        Parameters
        ==========

        K : Represents the desired domain or field

        Examples
        =====

        >>> from sympy import ZZ, CC
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)

        >>> domainMatrix.convert_to(CC)
        DomainMatrix([[(1.0 + 0.0j), (2.0 + 0.0j)],
        [(3.0 + 0.0j), (4.0 + 0.0j)]], (2, 2), CC)

        Returns
        =======

        DomainMatrix

        """
        return self.from_rep(self.rep.convert_to(K))

    def to_field(self):
        r"""
        Returns a DomainMatrix with the appropriate field

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)

        >>> domainMatrix.to_field()
        DomainMatrix([[1, 2], [3, 4]], (2, 2), QQ)

        Returns
        =======

        DomainMatrix

        """
        K = self.domain.get_field()
        return self.convert_to(K)

    def unify(self, other):
        r"""
        Unifies any two DomainMatrix with different domains

        Parameters
        ==========

        other : Represents the second(other) DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ, QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix1 = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> domainMatrix2 = domainMatrix1.convert_to(QQ)

        >>> dM1_new, dM2_new = domainMatrix1.unify(domainMatrix2)
        >>> dM1_new
        DomainMatrix([[1, 2], [3, 4]], (2, 2), QQ)
        >>> dM2_new
        DomainMatrix([[1, 2], [3, 4]], (2, 2), QQ)

        See Also
        ========

        :py:meth:`~.convert_to`

        Returns
        =======

        (dM1, dM2)
            dM1, dM2 are now two DomainMatrix with unified Domain

        """
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
        r"""
        Returns a MutableDenseMatrix for a DomainMatrix.

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)

        >>> mutableDomainMatrix = domainMatrix.to_Matrix()
        >>> mutableDomainMatrix
        Matrix([
        [1, 2],
        [3, 4]])
        >>> mutableDomainMatrix[0] = 4
        >>> mutableDomainMatrix
        Matrix([
        [4, 2],
        [3, 4]])

        See Also
        ========

        :py:meth:`~.from_Matrix`

        Returns
        =======

        Matrix

        """
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
        r"""
        Returns a DomainMatrix formed by stacking the rows horizontally (i.e. column wise)
        of the two DomainMatrix.
        Analogous to "numpy.hstack()"

        Parameters
        ==========

        A : Represents first DomainMatrix
        B : Represents second DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix1 = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> domainMatrix2 = DomainMatrix([
        ...    [ZZ(4), ZZ(3)],
        ...    [ZZ(2), ZZ(1)]], (2, 2), ZZ)

        >>> domainMatrix1.hstack(domainMatrix2)
        DomainMatrix([[1, 2, 4, 3], [3, 4, 2, 1]], (2, 4), ZZ)

        See Also
        ========

        :py:meth:`~.unify`

        Returns
        =======

        DomainMatrix

        """
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
        r"""
        Adds two DomainMatrix matrices of the same Domain

        Parameters
        ==========

        A : Represents first DomainMatrix
        B : Represents second DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix1 = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> domainMatrix2 = DomainMatrix([
        ...    [ZZ(4), ZZ(3)],
        ...    [ZZ(2), ZZ(1)]], (2, 2), ZZ)

        >>> domainMatrix1.add(domainMatrix2)
        DomainMatrix([[5, 5], [5, 5]], (2, 2), ZZ)

        See Also
        ========

        :py:meth:`~.sub`

        Returns
        =======

        DomainMatrix
            Addition of the DomainMatrix

        Raises
        ======

        ShapeError
            If the dimensions of the two DomainMatrix are not equal

        ValueError
            If the domain of the two DomainMatrix are not same

        """
        if A.shape != B.shape:
            raise ShapeError("shape")
        if A.domain != B.domain:
            raise ValueError("domain")
        return A.from_rep(A.rep.add(B.rep))

    def sub(A, B):
        r"""
        Subtracts two DomainMatrix matrices of the same Domain

        Parameters
        ==========

        A : Represents first DomainMatrix
        B : Represents second DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix1 = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> domainMatrix2 = DomainMatrix([
        ...    [ZZ(4), ZZ(3)],
        ...    [ZZ(2), ZZ(1)]], (2, 2), ZZ)

        >>> domainMatrix1.sub(domainMatrix2)
        DomainMatrix([[-3, -1], [1, 3]], (2, 2), ZZ)

        See Also
        ========

        :py:meth:`~.add`

        Returns
        =======

        DomainMatrix
            Substraction of the DomainMatrix

        Raises
        ======

        ShapeError
            If the dimensions of the two DomainMatrix are not equal

        ValueError
            If the domain of the two DomainMatrix are not same

        """
        if A.shape != B.shape:
            raise ShapeError("shape")
        if A.domain != B.domain:
            raise ValueError("domain")
        return A.from_rep(A.rep.sub(B.rep))

    def neg(A):
        r"""
        Returns the negative of DomainMatrix

        Parameters
        ==========

        A : Represents a DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)

        >>> domainMatrix.neg()
        DomainMatrix([[-1, -2], [-3, -4]], (2, 2), ZZ)

        See Also
        ========

        :py:meth:`~.inv`

        Returns
        =======

        DomainMatrix
            Negation of the DomainMatrix

        """
        return A.from_rep(A.rep.neg())

    def mul(A, b):
        r"""
        Performs term by term multiplication for the second DomainMatrix
        w.r.t first DomainMatrix. Returns a DomainMatrix whose rows are
        list of DomainMatrix matrices created after term by term multiplication.

        Parameters
        ==========

        A : Represents first DomainMatrix

        B : Represents second DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix1 = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> domainMatrix2 = DomainMatrix([
        ...    [ZZ(1), ZZ(1)],
        ...    [ZZ(0), ZZ(1)]], (2, 2), ZZ)

        >>> domainMatrix1.mul(domainMatrix2)
        DomainMatrix([[DomainMatrix([[1, 1], [0, 1]], (2, 2), ZZ),
        DomainMatrix([[2, 2], [0, 2]], (2, 2), ZZ)],
        [DomainMatrix([[3, 3], [0, 3]], (2, 2), ZZ),
        DomainMatrix([[4, 4], [0, 4]], (2, 2), ZZ)]], (2, 2), ZZ)

        See Also
        ========

        :py:meth:`~.matmul`

        Returns
        =======

        DomainMatrix
            term by term multiplication of two DomainMatrix matrices

        """
        return A.from_rep(A.rep.mul(b))

    def matmul(A, B):
        r"""
        Performs matrix multiplication of two DomainMatrix matrices

        Parameters
        ==========

        A : Represents first DomainMatrix

        B : Represents second DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix1 = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> domainMatrix2 = DomainMatrix([
        ...    [ZZ(1), ZZ(1)],
        ...    [ZZ(0), ZZ(1)]], (2, 2), ZZ)

        >>> domainMatrix1.matmul(domainMatrix2)
        DomainMatrix([[1, 3], [3, 7]], (2, 2), ZZ)

        See Also
        ========

        :py:meth:`~.mul`

        :py:meth:`~.pow`

        Returns
        =======

        DomainMatrix
            matrix multiplication of two DomainMatrix matrices

        """
        return A.from_rep(A.rep.matmul(B.rep))

    def pow(A, n):
        r"""
        Perfoms multiplication on the given DomainMatrix for the given number of times.
        Or else computes A^n

        Parameters
        ==========

        A : Represents a DomainMatrix

        n : Represents exponent for A

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix = DomainMatrix([
        ...    [ZZ(1), ZZ(1)],
        ...    [ZZ(0), ZZ(1)]], (2, 2), ZZ)

        >>> domainMatrix.pow(2)
        DomainMatrix([[1, 2], [0, 1]], (2, 2), ZZ)

        See Also
        ========

        :py:meth:`~.matmul`

        Returns
        =======

        DomainMatrix

        Raises
        ======

        NotImplementedError
            if n is negative.

        """
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
        r"""
        Returns reduced-row echelon form and list of pivots for the DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ, QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix = DomainMatrix([
        ...    [ZZ(1), ZZ(1)],
        ...    [ZZ(0), ZZ(1)]], (2, 2), ZZ)
        >>> dMF = domainMatrix.convert_to(QQ)
        >>> dMF
        DomainMatrix([[1, 1], [0, 1]], (2, 2), QQ)

        >>> dMF.rref()
        (DomainMatrix([[1, 0], [0, 1]], (2, 2), QQ), (0, 1))

        See Also
        ========

        :py:meth:`~.convert_to`

        :py:meth:`~.lu`

        Returns
        =======

        (DomainMatrix, list)
            reduced-row echelon form and list of pivots for the DomainMatrix

        Raises
        ======

        ValueError
            If the domain of DomainMatrix not a Field

        """
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        rref_ddm, pivots = self.rep.rref()
        return self.from_rep(rref_ddm), tuple(pivots)

    def nullspace(self):
        r"""
        Returns the Null Space for the DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ, QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix = DomainMatrix([
        ...    [ZZ(1), ZZ(1)],
        ...    [ZZ(0), ZZ(1)]], (2, 2), ZZ)
        >>> dMF = domainMatrix.convert_to(QQ)
        >>> dMF
        DomainMatrix([[1, 1], [0, 1]], (2, 2), QQ)

        >>> dMF.nullspace()
        DomainMatrix([], (0, 2), QQ)

        Returns
        =======

        DomainMatrix
            Null Space of the DomainMatrix

        """
        return self.from_rep(self.rep.nullspace()[0])

    def inv(self):
        r"""
        Finds the inverse of the DomainMatrix if exists

        Examples
        ========

        >>> from sympy import ZZ, QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix = DomainMatrix([
        ...    [ZZ(1), ZZ(1)],
        ...    [ZZ(0), ZZ(1)]], (2, 2), ZZ)
        >>> dMF = domainMatrix.convert_to(QQ)
        >>> dMF
        DomainMatrix([[1, 1], [0, 1]], (2, 2), QQ)

        >>> dMF.inv()
        DomainMatrix([[1, -1], [0, 1]], (2, 2), QQ)

        See Also
        ========

        :py:meth:`~.neg`

        Returns
        =======

        DomainMatrix
            Inverse of the DomainMatrix

        Raises
        ======

        ValueError
            If the domain of DomainMatrix not a Field

        NonSquareMatrixError
            If the DomainMatrix is not a not Square DomainMatrix

        """
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        m, n = self.shape
        if m != n:
            raise NonSquareMatrixError
        inv = self.rep.inv()
        return self.from_rep(inv)

    def det(self):
        r"""
        Returns the determinant value for a Square DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)

        >>> domainMatrix.det()
        -2

        Returns
        =======

        CC
            determinant value for a Square DomainMatrix

        Raises
        ======

        ValueError
            If the domain of DomainMatrix not a Field

        """
        m, n = self.shape
        if m != n:
            raise NonSquareMatrixError
        return self.rep.det()

    def lu(self):
        r"""
        Returns Lower and Upper decomposition(L, U triangular matrices)
        of the DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ, QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> dMF = domainMatrix.convert_to(QQ)
        >>> dMF
        DomainMatrix([[1, 2], [3, 4]], (2, 2), QQ)

        >>> dMF.lu()
        (DomainMatrix([[1, 0], [3, 1]], (2, 2), QQ),
        DomainMatrix([[1, 2], [0, -2]], (2, 2), QQ), [])

        See Also
        ========

        :py:meth:`~.lu_solve`

        Returns
        =======

        (L, U, exchange)
            L, U are Lower and Upper decomposition of the DomainMatrix,
            exchange is the list of indices of rows exchanged in the decomposition.

        Raises
        ======

        ValueError
            If the domain of DomainMatrix not a Field

        """
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        L, U, swaps = self.rep.lu()
        return self.from_rep(L), self.from_rep(U), swaps

    def lu_solve(self, rhs):
        r"""
        Solver for DomainMatrix x in the A*x = B, where A, B
        are DomainMatrix matrices.

        Parameters
        ==========

        rhs : Represents the DomainMatrix B

        Examples
        ========

        >>> from sympy import RR
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix1 = DomainMatrix([
        ...    [RR(1), RR(2)],
        ...    [RR(3), RR(4)]], (2, 2), RR)
        >>> domainMatrix2 = DomainMatrix([
        ...    [RR(1), RR(1)],
        ...    [RR(0), RR(1)]], (2, 2), RR)

        >>> domainMatrix1.lu_solve(domainMatrix2)
        DomainMatrix([[-2.0, -1.0], [1.5, 1.0]], (2, 2), RR)

        See Also
        ========

        :py:meth:`~.lu`

        Returns
        =======

        DomainMatrix
            Solution x of the equation, A*x = B

        Raises
        ======

        ShapeError
            If the DomainMatrix A and rhs have different number of rows

        ValueError
            If the domain of DomainMatrix A not a Field

        """
        if self.shape[0] != rhs.shape[0]:
            raise ShapeError("Shape")
        if not self.domain.is_Field:
            raise ValueError('Not a field')
        sol = self.rep.lu_solve(rhs.rep)
        return self.from_rep(sol)

    def charpoly(self):
        r"""
        Returns the coefficients of the characteristic polynomial of the DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)

        >>> domainMatrix.charpoly()
        [1, -5, -2]

        Returns
        =======

        list

        Raises
        ======

        NonSquareMatrixError
            If the DomainMatrix is not a not Square DomainMatrix

        """
        m, n = self.shape
        if m != n:
            raise NonSquareMatrixError("not square")
        return self.rep.charpoly()

    @classmethod
    def eye(cls, n, domain):
        return cls.from_rep(DDM.eye(n, domain))

    def __eq__(A, B):
        r"""
        Checks for two DomainMatrix matrices two be equal or not

        Parameters
        ==========

        A : Represents first DomainMatrix

        B : Represents second DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> domainMatrix1 = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> domainMatrix2 = DomainMatrix([
        ...    [ZZ(1), ZZ(1)],
        ...    [ZZ(0), ZZ(1)]], (2, 2), ZZ)

        >>> domainMatrix1.__eq__(domainMatrix1)
        True
        >>> domainMatrix1.__eq__(domainMatrix2)
        False

        Returns
        =======

        Boolean

        Raises
        ======

        NotImplementedError
            If B is not a DomainMatrix

        """
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        return A.rep == B.rep
