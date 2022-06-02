"""

Module for the DomainMatrix class.

A DomainMatrix represents a matrix with elements that are in a particular
Domain. Each DomainMatrix internally wraps a DDM which is used for the
lower-level operations. The idea is that the DomainMatrix class provides the
convenience routines for converting between Expr and the poly domains as well
as unifying matrices with different domains.

"""
from functools import reduce
from typing import Union as tUnion, Tuple as tTuple

from sympy.core.sympify import _sympify

from ..domains import Domain

from ..constructor import construct_domain

from .exceptions import (DMNonSquareMatrixError, DMShapeError,
                         DMDomainError, DMFormatError, DMBadInputError,
                         DMNotAField)

from .ddm import DDM

from .sdm import SDM

from .domainscalar import DomainScalar

from sympy.polys.domains import ZZ, EXRAW


def DM(rows, domain):
    """Convenient alias for DomainMatrix.from_list

    Examples
    =======

    >>> from sympy import ZZ
    >>> from sympy.polys.matrices import DM
    >>> DM([[1, 2], [3, 4]], ZZ)
    DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ)

    See also
    =======

    DomainMatrix.from_list
    """
    return DomainMatrix.from_list(rows, domain)


class DomainMatrix:
    r"""
    Associate Matrix with :py:class:`~.Domain`

    Explanation
    ===========

    DomainMatrix uses :py:class:`~.Domain` for its internal representation
    which makes it more faster for many common operations
    than current SymPy Matrix class, but this advantage makes it not
    entirely compatible with Matrix.
    DomainMatrix could be found analogous to numpy arrays with "dtype".
    In the DomainMatrix, each matrix has a domain such as :ref:`ZZ`
    or  :ref:`QQ(a)`.


    Examples
    ========

    Creating a DomainMatrix from the existing Matrix class:

    >>> from sympy import Matrix
    >>> from sympy.polys.matrices import DomainMatrix
    >>> Matrix1 = Matrix([
    ...    [1, 2],
    ...    [3, 4]])
    >>> A = DomainMatrix.from_Matrix(Matrix1)
    >>> A
    DomainMatrix({0: {0: 1, 1: 2}, 1: {0: 3, 1: 4}}, (2, 2), ZZ)

    Driectly forming a DomainMatrix:

    >>> from sympy import ZZ
    >>> from sympy.polys.matrices import DomainMatrix
    >>> A = DomainMatrix([
    ...    [ZZ(1), ZZ(2)],
    ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
    >>> A
    DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ)

    See Also
    ========

    DDM
    SDM
    Domain
    Poly

    """
    rep: tUnion[SDM, DDM]
    shape: tTuple[int, int]
    domain: Domain

    def __new__(cls, rows, shape, domain, *, fmt=None):
        """
        Creates a :py:class:`~.DomainMatrix`.

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
        if isinstance(rows, (DDM, SDM)):
            raise TypeError("Use from_rep to initialise from SDM/DDM")
        elif isinstance(rows, list):
            rep = DDM(rows, shape, domain)
        elif isinstance(rows, dict):
            rep = SDM(rows, shape, domain)
        else:
            msg = "Input should be list-of-lists or dict-of-dicts"
            raise TypeError(msg)

        if fmt is not None:
            if fmt == 'sparse':
                rep = rep.to_sdm()
            elif fmt == 'dense':
                rep = rep.to_ddm()
            else:
                raise ValueError("fmt should be 'sparse' or 'dense'")

        return cls.from_rep(rep)

    def __getnewargs__(self):
        rep = self.rep
        if isinstance(rep, DDM):
            arg = list(rep)
        elif isinstance(rep, SDM):
            arg = dict(rep)
        else:
            raise RuntimeError # pragma: no cover
        return arg, self.shape, self.domain

    def __getitem__(self, key):
        i, j = key
        m, n = self.shape
        if not (isinstance(i, slice) or isinstance(j, slice)):
            return DomainScalar(self.rep.getitem(i, j), self.domain)

        if not isinstance(i, slice):
            if not -m <= i < m:
                raise IndexError("Row index out of range")
            i = i % m
            i = slice(i, i+1)
        if not isinstance(j, slice):
            if not -n <= j < n:
                raise IndexError("Column index out of range")
            j = j % n
            j = slice(j, j+1)

        return self.from_rep(self.rep.extract_slice(i, j))

    def getitem_sympy(self, i, j):
        return self.domain.to_sympy(self.rep.getitem(i, j))

    def extract(self, rowslist, colslist):
        return self.from_rep(self.rep.extract(rowslist, colslist))

    def __setitem__(self, key, value):
        i, j = key
        if not self.domain.of_type(value):
            raise TypeError
        if isinstance(i, int) and isinstance(j, int):
            self.rep.setitem(i, j, value)
        else:
            raise NotImplementedError

    @classmethod
    def from_rep(cls, rep):
        """Create a new DomainMatrix efficiently from DDM/SDM.

        Examples
        ========

        Create a :py:class:`~.DomainMatrix` with an dense internal
        representation as :py:class:`~.DDM`:

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> from sympy.polys.matrices.ddm import DDM
        >>> drep = DDM([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> dM = DomainMatrix.from_rep(drep)
        >>> dM
        DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ)

        Create a :py:class:`~.DomainMatrix` with a sparse internal
        representation as :py:class:`~.SDM`:

        >>> from sympy.polys.matrices import DomainMatrix
        >>> from sympy.polys.matrices.sdm import SDM
        >>> from sympy import ZZ
        >>> drep = SDM({0:{1:ZZ(1)},1:{0:ZZ(2)}}, (2, 2), ZZ)
        >>> dM = DomainMatrix.from_rep(drep)
        >>> dM
        DomainMatrix({0: {1: 1}, 1: {0: 2}}, (2, 2), ZZ)

        Parameters
        ==========

        rep: SDM or DDM
            The internal sparse or dense representation of the matrix.

        Returns
        =======

        DomainMatrix
            A :py:class:`~.DomainMatrix` wrapping *rep*.

        Notes
        =====

        This takes ownership of rep as its internal representation. If rep is
        being mutated elsewhere then a copy should be provided to
        ``from_rep``. Only minimal verification or checking is done on *rep*
        as this is supposed to be an efficient internal routine.

        """
        if not isinstance(rep, (DDM, SDM)):
            raise TypeError("rep should be of type DDM or SDM")
        self = super().__new__(cls)
        self.rep = rep
        self.shape = rep.shape
        self.domain = rep.domain
        return self


    @classmethod
    def from_list(cls, rows, domain):
        r"""
        Convert a list of lists into a DomainMatrix

        Parameters
        ==========

        rows: list of lists
            Each element of the inner lists should be either the single arg,
            or tuple of args, that would be passed to the domain constructor
            in order to form an element of the domain. See examples.

        Returns
        =======

        DomainMatrix containing elements defined in rows

        Examples
        ========

        >>> from sympy.polys.matrices import DomainMatrix
        >>> from sympy import FF, QQ, ZZ
        >>> A = DomainMatrix.from_list([[1, 0, 1], [0, 0, 1]], ZZ)
        >>> A
        DomainMatrix([[1, 0, 1], [0, 0, 1]], (2, 3), ZZ)
        >>> B = DomainMatrix.from_list([[1, 0, 1], [0, 0, 1]], FF(7))
        >>> B
        DomainMatrix([[1 mod 7, 0 mod 7, 1 mod 7], [0 mod 7, 0 mod 7, 1 mod 7]], (2, 3), GF(7))
        >>> C = DomainMatrix.from_list([[(1, 2), (3, 1)], [(1, 4), (5, 1)]], QQ)
        >>> C
        DomainMatrix([[1/2, 3], [1/4, 5]], (2, 2), QQ)

        See Also
        ========

        from_list_sympy

        """
        nrows = len(rows)
        ncols = 0 if not nrows else len(rows[0])
        conv = lambda e: domain(*e) if isinstance(e, tuple) else domain(e)
        domain_rows = [[conv(e) for e in row] for row in rows]
        return DomainMatrix(domain_rows, (nrows, ncols), domain)


    @classmethod
    def from_list_sympy(cls, nrows, ncols, rows, **kwargs):
        r"""
        Convert a list of lists of Expr into a DomainMatrix using construct_domain

        Parameters
        ==========

        nrows: number of rows
        ncols: number of columns
        rows: list of lists

        Returns
        =======

        DomainMatrix containing elements of rows

        Examples
        ========

        >>> from sympy.polys.matrices import DomainMatrix
        >>> from sympy.abc import x, y, z
        >>> A = DomainMatrix.from_list_sympy(1, 3, [[x, y, z]])
        >>> A
        DomainMatrix([[x, y, z]], (1, 3), ZZ[x,y,z])

        See Also
        ========

        sympy.polys.constructor.construct_domain, from_dict_sympy

        """
        assert len(rows) == nrows
        assert all(len(row) == ncols for row in rows)

        items_sympy = [_sympify(item) for row in rows for item in row]

        domain, items_domain = cls.get_domain(items_sympy, **kwargs)

        domain_rows = [[items_domain[ncols*r + c] for c in range(ncols)] for r in range(nrows)]

        return DomainMatrix(domain_rows, (nrows, ncols), domain)

    @classmethod
    def from_dict_sympy(cls, nrows, ncols, elemsdict, **kwargs):
        """

        Parameters
        ==========

        nrows: number of rows
        ncols: number of cols
        elemsdict: dict of dicts containing non-zero elements of the DomainMatrix

        Returns
        =======

        DomainMatrix containing elements of elemsdict

        Examples
        ========

        >>> from sympy.polys.matrices import DomainMatrix
        >>> from sympy.abc import x,y,z
        >>> elemsdict = {0: {0:x}, 1:{1: y}, 2: {2: z}}
        >>> A = DomainMatrix.from_dict_sympy(3, 3, elemsdict)
        >>> A
        DomainMatrix({0: {0: x}, 1: {1: y}, 2: {2: z}}, (3, 3), ZZ[x,y,z])

        See Also
        ========

        from_list_sympy

        """
        if not all(0 <= r < nrows for r in elemsdict):
            raise DMBadInputError("Row out of range")
        if not all(0 <= c < ncols for row in elemsdict.values() for c in row):
            raise DMBadInputError("Column out of range")

        items_sympy = [_sympify(item) for row in elemsdict.values() for item in row.values()]
        domain, items_domain = cls.get_domain(items_sympy, **kwargs)

        idx = 0
        items_dict = {}
        for i, row in elemsdict.items():
                items_dict[i] = {}
                for j in row:
                    items_dict[i][j] = items_domain[idx]
                    idx += 1

        return DomainMatrix(items_dict, (nrows, ncols), domain)

    @classmethod
    def from_Matrix(cls, M, fmt='sparse',**kwargs):
        r"""
        Convert Matrix to DomainMatrix

        Parameters
        ==========

        M: Matrix

        Returns
        =======

        Returns DomainMatrix with identical elements as M

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.polys.matrices import DomainMatrix
        >>> M = Matrix([
        ...    [1.0, 3.4],
        ...    [2.4, 1]])
        >>> A = DomainMatrix.from_Matrix(M)
        >>> A
        DomainMatrix({0: {0: 1.0, 1: 3.4}, 1: {0: 2.4, 1: 1.0}}, (2, 2), RR)

        We can keep internal representation as ddm using fmt='dense'
        >>> from sympy import Matrix, QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix.from_Matrix(Matrix([[QQ(1, 2), QQ(3, 4)], [QQ(0, 1), QQ(0, 1)]]), fmt='dense')
        >>> A.rep
        [[1/2, 3/4], [0, 0]]

        See Also
        ========

        Matrix

        """
        if fmt == 'dense':
            return cls.from_list_sympy(*M.shape, M.tolist(), **kwargs)

        return cls.from_dict_sympy(*M.shape, M.todod(), **kwargs)

    @classmethod
    def get_domain(cls, items_sympy, **kwargs):
        K, items_K = construct_domain(items_sympy, **kwargs)
        return K, items_K

    def copy(self):
        return self.from_rep(self.rep.copy())

    def convert_to(self, K):
        r"""
        Change the domain of DomainMatrix to desired domain or field

        Parameters
        ==========

        K : Represents the desired domain or field.
            Alternatively, ``None`` may be passed, in which case this method
            just returns a copy of this DomainMatrix.

        Returns
        =======

        DomainMatrix
            DomainMatrix with the desired domain or field

        Examples
        ========

        >>> from sympy import ZZ, ZZ_I
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)

        >>> A.convert_to(ZZ_I)
        DomainMatrix([[1, 2], [3, 4]], (2, 2), ZZ_I)

        """
        if K is None:
            return self.copy()
        return self.from_rep(self.rep.convert_to(K))

    def to_sympy(self):
        return self.convert_to(EXRAW)

    def to_field(self):
        r"""
        Returns a DomainMatrix with the appropriate field

        Returns
        =======

        DomainMatrix
            DomainMatrix with the appropriate field

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)

        >>> A.to_field()
        DomainMatrix([[1, 2], [3, 4]], (2, 2), QQ)

        """
        K = self.domain.get_field()
        return self.convert_to(K)

    def to_sparse(self):
        """
        Return a sparse DomainMatrix representation of *self*.

        Examples
        ========

        >>> from sympy.polys.matrices import DomainMatrix
        >>> from sympy import QQ
        >>> A = DomainMatrix([[1, 0],[0, 2]], (2, 2), QQ)
        >>> A.rep
        [[1, 0], [0, 2]]
        >>> B = A.to_sparse()
        >>> B.rep
        {0: {0: 1}, 1: {1: 2}}
        """
        if self.rep.fmt == 'sparse':
            return self

        return self.from_rep(SDM.from_ddm(self.rep))

    def to_dense(self):
        """
        Return a dense DomainMatrix representation of *self*.

        Examples
        ========

        >>> from sympy.polys.matrices import DomainMatrix
        >>> from sympy import QQ
        >>> A = DomainMatrix({0: {0: 1}, 1: {1: 2}}, (2, 2), QQ)
        >>> A.rep
        {0: {0: 1}, 1: {1: 2}}
        >>> B = A.to_dense()
        >>> B.rep
        [[1, 0], [0, 2]]

        """
        if self.rep.fmt == 'dense':
            return self

        return self.from_rep(SDM.to_ddm(self.rep))

    @classmethod
    def _unify_domain(cls, *matrices):
        """Convert matrices to a common domain"""
        domains = {matrix.domain for matrix in matrices}
        if len(domains) == 1:
            return matrices
        domain = reduce(lambda x, y: x.unify(y), domains)
        return tuple(matrix.convert_to(domain) for matrix in matrices)

    @classmethod
    def _unify_fmt(cls, *matrices, fmt=None):
        """Convert matrices to the same format.

        If all matrices have the same format, then return unmodified.
        Otherwise convert both to the preferred format given as *fmt* which
        should be 'dense' or 'sparse'.
        """
        formats = {matrix.rep.fmt for matrix in matrices}
        if len(formats) == 1:
            return matrices
        if fmt == 'sparse':
            return tuple(matrix.to_sparse() for matrix in matrices)
        elif fmt == 'dense':
            return tuple(matrix.to_dense() for matrix in matrices)
        else:
            raise ValueError("fmt should be 'sparse' or 'dense'")

    def unify(self, *others, fmt=None):
        """
        Unifies the domains and the format of self and other
        matrices.

        Parameters
        ==========

        others : DomainMatrix

        fmt: string 'dense', 'sparse' or `None` (default)
            The preferred format to convert to if self and other are not
            already in the same format. If `None` or not specified then no
            conversion if performed.

        Returns
        =======

        Tuple[DomainMatrix]
            Matrices with unified domain and format

        Examples
        ========

        Unify the domain of DomainMatrix that have different domains:

        >>> from sympy import ZZ, QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([[ZZ(1), ZZ(2)]], (1, 2), ZZ)
        >>> B = DomainMatrix([[QQ(1, 2), QQ(2)]], (1, 2), QQ)
        >>> Aq, Bq = A.unify(B)
        >>> Aq
        DomainMatrix([[1, 2]], (1, 2), QQ)
        >>> Bq
        DomainMatrix([[1/2, 2]], (1, 2), QQ)

        Unify the format (dense or sparse):

        >>> A = DomainMatrix([[ZZ(1), ZZ(2)]], (1, 2), ZZ)
        >>> B = DomainMatrix({0:{0: ZZ(1)}}, (2, 2), ZZ)
        >>> B.rep
        {0: {0: 1}}

        >>> A2, B2 = A.unify(B, fmt='dense')
        >>> B2.rep
        [[1, 0], [0, 0]]

        See Also
        ========

        convert_to, to_dense, to_sparse

        """
        matrices = (self,) + others
        matrices = DomainMatrix._unify_domain(*matrices)
        if fmt is not None:
            matrices = DomainMatrix._unify_fmt(*matrices, fmt=fmt)
        return matrices

    def to_Matrix(self):
        r"""
        Convert DomainMatrix to Matrix

        Returns
        =======

        Matrix
            MutableDenseMatrix for the DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)

        >>> A.to_Matrix()
        Matrix([
            [1, 2],
            [3, 4]])

        See Also
        ========

        from_Matrix

        """
        from sympy.matrices.dense import MutableDenseMatrix
        elemlist = self.rep.to_list()
        elements_sympy = [self.domain.to_sympy(e) for row in elemlist for e in row]
        return MutableDenseMatrix(*self.shape, elements_sympy)

    def to_list(self):
        return self.rep.to_list()

    def to_list_flat(self):
        return self.rep.to_list_flat()

    def to_dok(self):
        return self.rep.to_dok()

    def __repr__(self):
        return 'DomainMatrix(%s, %r, %r)' % (str(self.rep), self.shape, self.domain)

    def transpose(self):
        """Matrix transpose of ``self``"""
        return self.from_rep(self.rep.transpose())

    def flat(self):
        rows, cols = self.shape
        return [self[i,j].element for i in range(rows) for j in range(cols)]

    @property
    def is_zero_matrix(self):
        return self.rep.is_zero_matrix()

    @property
    def is_upper(self):
        """
        Says whether this matrix is upper-triangular. True can be returned
        even if the matrix is not square.
        """
        return self.rep.is_upper()

    @property
    def is_lower(self):
        """
        Says whether this matrix is lower-triangular. True can be returned
        even if the matrix is not square.
        """
        return self.rep.is_lower()

    @property
    def is_square(self):
        return self.shape[0] == self.shape[1]

    def rank(self):
        rref, pivots = self.rref()
        return len(pivots)

    def hstack(A, *B):
        r"""Horizontally stack the given matrices.

        Parameters
        ==========

        B: DomainMatrix
            Matrices to stack horizontally.

        Returns
        =======

        DomainMatrix
            DomainMatrix by stacking horizontally.

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix

        >>> A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> B = DomainMatrix([[ZZ(5), ZZ(6)], [ZZ(7), ZZ(8)]], (2, 2), ZZ)
        >>> A.hstack(B)
        DomainMatrix([[1, 2, 5, 6], [3, 4, 7, 8]], (2, 4), ZZ)

        >>> C = DomainMatrix([[ZZ(9), ZZ(10)], [ZZ(11), ZZ(12)]], (2, 2), ZZ)
        >>> A.hstack(B, C)
        DomainMatrix([[1, 2, 5, 6, 9, 10], [3, 4, 7, 8, 11, 12]], (2, 6), ZZ)

        See Also
        ========

        unify
        """
        A, *B = A.unify(*B, fmt='dense')
        return DomainMatrix.from_rep(A.rep.hstack(*(Bk.rep for Bk in B)))

    def vstack(A, *B):
        r"""Vertically stack the given matrices.

        Parameters
        ==========

        B: DomainMatrix
            Matrices to stack vertically.

        Returns
        =======

        DomainMatrix
            DomainMatrix by stacking vertically.

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix

        >>> A = DomainMatrix([[ZZ(1), ZZ(2)], [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> B = DomainMatrix([[ZZ(5), ZZ(6)], [ZZ(7), ZZ(8)]], (2, 2), ZZ)
        >>> A.vstack(B)
        DomainMatrix([[1, 2], [3, 4], [5, 6], [7, 8]], (4, 2), ZZ)

        >>> C = DomainMatrix([[ZZ(9), ZZ(10)], [ZZ(11), ZZ(12)]], (2, 2), ZZ)
        >>> A.vstack(B, C)
        DomainMatrix([[1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [11, 12]], (6, 2), ZZ)

        See Also
        ========

        unify
        """
        A, *B = A.unify(*B, fmt='dense')
        return DomainMatrix.from_rep(A.rep.vstack(*(Bk.rep for Bk in B)))

    def applyfunc(self, func, domain=None):
        if domain is None:
            domain = self.domain
        return self.from_rep(self.rep.applyfunc(func, domain))

    def __add__(A, B):
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        A, B = A.unify(B, fmt='dense')
        return A.add(B)

    def __sub__(A, B):
        if not isinstance(B, DomainMatrix):
            return NotImplemented
        A, B = A.unify(B, fmt='dense')
        return A.sub(B)

    def __neg__(A):
        return A.neg()

    def __mul__(A, B):
        """A * B"""
        if isinstance(B, DomainMatrix):
            A, B = A.unify(B, fmt='dense')
            return A.matmul(B)
        elif B in A.domain:
            return A.scalarmul(B)
        elif isinstance(B, DomainScalar):
            A, B = A.unify(B)
            return A.scalarmul(B.element)
        else:
            return NotImplemented

    def __rmul__(A, B):
        if B in A.domain:
            return A.rscalarmul(B)
        elif isinstance(B, DomainScalar):
            A, B = A.unify(B)
            return A.rscalarmul(B.element)
        else:
            return NotImplemented

    def __pow__(A, n):
        """A ** n"""
        if not isinstance(n, int):
            return NotImplemented
        return A.pow(n)

    def _check(a, op, b, ashape, bshape):
        if a.domain != b.domain:
            msg = "Domain mismatch: %s %s %s" % (a.domain, op, b.domain)
            raise DMDomainError(msg)
        if ashape != bshape:
            msg = "Shape mismatch: %s %s %s" % (a.shape, op, b.shape)
            raise DMShapeError(msg)
        if a.rep.fmt != b.rep.fmt:
            msg = "Format mismatch: %s %s %s" % (a.rep.fmt, op, b.rep.fmt)
            raise DMFormatError(msg)

    def add(A, B):
        r"""
        Adds two DomainMatrix matrices of the same Domain

        Parameters
        ==========

        A, B: DomainMatrix
            matrices to add

        Returns
        =======

        DomainMatrix
            DomainMatrix after Addition

        Raises
        ======

        DMShapeError
            If the dimensions of the two DomainMatrix are not equal

        ValueError
            If the domain of the two DomainMatrix are not same

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> B = DomainMatrix([
        ...    [ZZ(4), ZZ(3)],
        ...    [ZZ(2), ZZ(1)]], (2, 2), ZZ)

        >>> A.add(B)
        DomainMatrix([[5, 5], [5, 5]], (2, 2), ZZ)

        See Also
        ========

        sub, matmul

        """
        A._check('+', B, A.shape, B.shape)
        return A.from_rep(A.rep.add(B.rep))


    def sub(A, B):
        r"""
        Subtracts two DomainMatrix matrices of the same Domain

        Parameters
        ==========

        A, B: DomainMatrix
            matrices to substract

        Returns
        =======

        DomainMatrix
            DomainMatrix after Substraction

        Raises
        ======

        DMShapeError
            If the dimensions of the two DomainMatrix are not equal

        ValueError
            If the domain of the two DomainMatrix are not same

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> B = DomainMatrix([
        ...    [ZZ(4), ZZ(3)],
        ...    [ZZ(2), ZZ(1)]], (2, 2), ZZ)

        >>> A.sub(B)
        DomainMatrix([[-3, -1], [1, 3]], (2, 2), ZZ)

        See Also
        ========

        add, matmul

        """
        A._check('-', B, A.shape, B.shape)
        return A.from_rep(A.rep.sub(B.rep))

    def neg(A):
        r"""
        Returns the negative of DomainMatrix

        Parameters
        ==========

        A : Represents a DomainMatrix

        Returns
        =======

        DomainMatrix
            DomainMatrix after Negation

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)

        >>> A.neg()
        DomainMatrix([[-1, -2], [-3, -4]], (2, 2), ZZ)

        """
        return A.from_rep(A.rep.neg())

    def mul(A, b):
        r"""
        Performs term by term multiplication for the second DomainMatrix
        w.r.t first DomainMatrix. Returns a DomainMatrix whose rows are
        list of DomainMatrix matrices created after term by term multiplication.

        Parameters
        ==========

        A, B: DomainMatrix
            matrices to multiply term-wise

        Returns
        =======

        DomainMatrix
            DomainMatrix after term by term multiplication

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> B = DomainMatrix([
        ...    [ZZ(1), ZZ(1)],
        ...    [ZZ(0), ZZ(1)]], (2, 2), ZZ)

        >>> A.mul(B)
        DomainMatrix([[DomainMatrix([[1, 1], [0, 1]], (2, 2), ZZ),
        DomainMatrix([[2, 2], [0, 2]], (2, 2), ZZ)],
        [DomainMatrix([[3, 3], [0, 3]], (2, 2), ZZ),
        DomainMatrix([[4, 4], [0, 4]], (2, 2), ZZ)]], (2, 2), ZZ)

        See Also
        ========

        matmul

        """
        return A.from_rep(A.rep.mul(b))

    def rmul(A, b):
        return A.from_rep(A.rep.rmul(b))

    def matmul(A, B):
        r"""
        Performs matrix multiplication of two DomainMatrix matrices

        Parameters
        ==========

        A, B: DomainMatrix
            to multiply

        Returns
        =======

        DomainMatrix
            DomainMatrix after multiplication

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> B = DomainMatrix([
        ...    [ZZ(1), ZZ(1)],
        ...    [ZZ(0), ZZ(1)]], (2, 2), ZZ)

        >>> A.matmul(B)
        DomainMatrix([[1, 3], [3, 7]], (2, 2), ZZ)

        See Also
        ========

        mul, pow, add, sub

        """

        A._check('*', B, A.shape[1], B.shape[0])
        return A.from_rep(A.rep.matmul(B.rep))

    def _scalarmul(A, lamda, reverse):
        if lamda == A.domain.zero:
            return DomainMatrix.zeros(A.shape, A.domain)
        elif lamda == A.domain.one:
            return A.copy()
        elif reverse:
            return A.rmul(lamda)
        else:
            return A.mul(lamda)

    def scalarmul(A, lamda):
        return A._scalarmul(lamda, reverse=False)

    def rscalarmul(A, lamda):
        return A._scalarmul(lamda, reverse=True)

    def mul_elementwise(A, B):
        assert A.domain == B.domain
        return A.from_rep(A.rep.mul_elementwise(B.rep))

    def __truediv__(A, lamda):
        """ Method for Scalar Divison"""
        if isinstance(lamda, int) or ZZ.of_type(lamda):
            lamda = DomainScalar(ZZ(lamda), ZZ)

        if not isinstance(lamda, DomainScalar):
            return NotImplemented

        A, lamda = A.to_field().unify(lamda)
        if lamda.element == lamda.domain.zero:
            raise ZeroDivisionError
        if lamda.element == lamda.domain.one:
            return A.to_field()

        return A.mul(1 / lamda.element)

    def pow(A, n):
        r"""
        Computes A**n

        Parameters
        ==========

        A : DomainMatrix

        n : exponent for A

        Returns
        =======

        DomainMatrix
            DomainMatrix on computing A**n

        Raises
        ======

        NotImplementedError
            if n is negative.

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [ZZ(1), ZZ(1)],
        ...    [ZZ(0), ZZ(1)]], (2, 2), ZZ)

        >>> A.pow(2)
        DomainMatrix([[1, 2], [0, 1]], (2, 2), ZZ)

        See Also
        ========

        matmul

        """
        nrows, ncols = A.shape
        if nrows != ncols:
            raise DMNonSquareMatrixError('Power of a nonsquare matrix')
        if n < 0:
            raise NotImplementedError('Negative powers')
        elif n == 0:
            return A.eye(nrows, A.domain)
        elif n == 1:
            return A
        elif n % 2 == 1:
            return A * A**(n - 1)
        else:
            sqrtAn = A ** (n // 2)
            return sqrtAn * sqrtAn

    def scc(self):
        """Compute the strongly connected components of a DomainMatrix

        Explanation
        ===========

        A square matrix can be considered as the adjacency matrix for a
        directed graph where the row and column indices are the vertices. In
        this graph if there is an edge from vertex ``i`` to vertex ``j`` if
        ``M[i, j]`` is nonzero. This routine computes the strongly connected
        components of that graph which are subsets of the rows and columns that
        are connected by some nonzero element of the matrix. The strongly
        connected components are useful because many operations such as the
        determinant can be computed by working with the submatrices
        corresponding to each component.

        Examples
        ========

        Find the strongly connected components of a matrix:

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> M = DomainMatrix([[ZZ(1), ZZ(0), ZZ(2)],
        ...                   [ZZ(0), ZZ(3), ZZ(0)],
        ...                   [ZZ(4), ZZ(6), ZZ(5)]], (3, 3), ZZ)
        >>> M.scc()
        [[1], [0, 2]]

        Compute the determinant from the components:

        >>> MM = M.to_Matrix()
        >>> MM
        Matrix([
        [1, 0, 2],
        [0, 3, 0],
        [4, 6, 5]])
        >>> MM[[1], [1]]
        Matrix([[3]])
        >>> MM[[0, 2], [0, 2]]
        Matrix([
        [1, 2],
        [4, 5]])
        >>> MM.det()
        -9
        >>> MM[[1], [1]].det() * MM[[0, 2], [0, 2]].det()
        -9

        The components are given in reverse topological order and represent a
        permutation of the rows and columns that will bring the matrix into
        block lower-triangular form:

        >>> MM[[1, 0, 2], [1, 0, 2]]
        Matrix([
        [3, 0, 0],
        [0, 1, 2],
        [6, 4, 5]])

        Returns
        =======

        List of lists of integers
            Each list represents a strongly connected component.

        See also
        ========

        sympy.matrices.matrices.MatrixBase.strongly_connected_components
        sympy.utilities.iterables.strongly_connected_components

        """
        rows, cols = self.shape
        assert rows == cols
        return self.rep.scc()

    def rref(self):
        r"""
        Returns reduced-row echelon form and list of pivots for the DomainMatrix

        Returns
        =======

        (DomainMatrix, list)
            reduced-row echelon form and list of pivots for the DomainMatrix

        Raises
        ======

        ValueError
            If the domain of DomainMatrix not a Field

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...     [QQ(2), QQ(-1), QQ(0)],
        ...     [QQ(-1), QQ(2), QQ(-1)],
        ...     [QQ(0), QQ(0), QQ(2)]], (3, 3), QQ)

        >>> rref_matrix, rref_pivots = A.rref()
        >>> rref_matrix
        DomainMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]], (3, 3), QQ)
        >>> rref_pivots
        (0, 1, 2)

        See Also
        ========

        convert_to, lu

        """
        if not self.domain.is_Field:
            raise DMNotAField('Not a field')
        rref_ddm, pivots = self.rep.rref()
        return self.from_rep(rref_ddm), tuple(pivots)

    def columnspace(self):
        r"""
        Returns the columnspace for the DomainMatrix

        Returns
        =======

        DomainMatrix
            The columns of this matrix form a basis for the columnspace.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [QQ(1), QQ(-1)],
        ...    [QQ(2), QQ(-2)]], (2, 2), QQ)
        >>> A.columnspace()
        DomainMatrix([[1], [2]], (2, 1), QQ)

        """
        if not self.domain.is_Field:
            raise DMNotAField('Not a field')
        rref, pivots = self.rref()
        rows, cols = self.shape
        return self.extract(range(rows), pivots)

    def rowspace(self):
        r"""
        Returns the rowspace for the DomainMatrix

        Returns
        =======

        DomainMatrix
            The rows of this matrix form a basis for the rowspace.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [QQ(1), QQ(-1)],
        ...    [QQ(2), QQ(-2)]], (2, 2), QQ)
        >>> A.rowspace()
        DomainMatrix([[1, -1]], (1, 2), QQ)

        """
        if not self.domain.is_Field:
            raise DMNotAField('Not a field')
        rref, pivots = self.rref()
        rows, cols = self.shape
        return self.extract(range(len(pivots)), range(cols))

    def nullspace(self):
        r"""
        Returns the nullspace for the DomainMatrix

        Returns
        =======

        DomainMatrix
            The rows of this matrix form a basis for the nullspace.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [QQ(1), QQ(-1)],
        ...    [QQ(2), QQ(-2)]], (2, 2), QQ)
        >>> A.nullspace()
        DomainMatrix([[1, 1]], (1, 2), QQ)

        """
        if not self.domain.is_Field:
            raise DMNotAField('Not a field')
        return self.from_rep(self.rep.nullspace()[0])

    def inv(self):
        r"""
        Finds the inverse of the DomainMatrix if exists

        Returns
        =======

        DomainMatrix
            DomainMatrix after inverse

        Raises
        ======

        ValueError
            If the domain of DomainMatrix not a Field

        DMNonSquareMatrixError
            If the DomainMatrix is not a not Square DomainMatrix

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...     [QQ(2), QQ(-1), QQ(0)],
        ...     [QQ(-1), QQ(2), QQ(-1)],
        ...     [QQ(0), QQ(0), QQ(2)]], (3, 3), QQ)
        >>> A.inv()
        DomainMatrix([[2/3, 1/3, 1/6], [1/3, 2/3, 1/3], [0, 0, 1/2]], (3, 3), QQ)

        See Also
        ========

        neg

        """
        if not self.domain.is_Field:
            raise DMNotAField('Not a field')
        m, n = self.shape
        if m != n:
            raise DMNonSquareMatrixError
        inv = self.rep.inv()
        return self.from_rep(inv)

    def det(self):
        r"""
        Returns the determinant of a Square DomainMatrix

        Returns
        =======

        S.Complexes
            determinant of Square DomainMatrix

        Raises
        ======

        ValueError
            If the domain of DomainMatrix not a Field

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)

        >>> A.det()
        -2

        """
        m, n = self.shape
        if m != n:
            raise DMNonSquareMatrixError
        return self.rep.det()

    def lu(self):
        r"""
        Returns Lower and Upper decomposition of the DomainMatrix

        Returns
        =======

        (L, U, exchange)
            L, U are Lower and Upper decomposition of the DomainMatrix,
            exchange is the list of indices of rows exchanged in the decomposition.

        Raises
        ======

        ValueError
            If the domain of DomainMatrix not a Field

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [QQ(1), QQ(-1)],
        ...    [QQ(2), QQ(-2)]], (2, 2), QQ)
        >>> A.lu()
        (DomainMatrix([[1, 0], [2, 1]], (2, 2), QQ), DomainMatrix([[1, -1], [0, 0]], (2, 2), QQ), [])

        See Also
        ========

        lu_solve

        """
        if not self.domain.is_Field:
            raise DMNotAField('Not a field')
        L, U, swaps = self.rep.lu()
        return self.from_rep(L), self.from_rep(U), swaps

    def lu_solve(self, rhs):
        r"""
        Solver for DomainMatrix x in the A*x = B

        Parameters
        ==========

        rhs : DomainMatrix B

        Returns
        =======

        DomainMatrix
            x in A*x = B

        Raises
        ======

        DMShapeError
            If the DomainMatrix A and rhs have different number of rows

        ValueError
            If the domain of DomainMatrix A not a Field

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [QQ(1), QQ(2)],
        ...    [QQ(3), QQ(4)]], (2, 2), QQ)
        >>> B = DomainMatrix([
        ...    [QQ(1), QQ(1)],
        ...    [QQ(0), QQ(1)]], (2, 2), QQ)

        >>> A.lu_solve(B)
        DomainMatrix([[-2, -1], [3/2, 1]], (2, 2), QQ)

        See Also
        ========

        lu

        """
        if self.shape[0] != rhs.shape[0]:
            raise DMShapeError("Shape")
        if not self.domain.is_Field:
            raise DMNotAField('Not a field')
        sol = self.rep.lu_solve(rhs.rep)
        return self.from_rep(sol)

    def _solve(A, b):
        # XXX: Not sure about this method or its signature. It is just created
        # because it is needed by the holonomic module.
        if A.shape[0] != b.shape[0]:
            raise DMShapeError("Shape")
        if A.domain != b.domain or not A.domain.is_Field:
            raise DMNotAField('Not a field')
        Aaug = A.hstack(b)
        Arref, pivots = Aaug.rref()
        particular = Arref.from_rep(Arref.rep.particular())
        nullspace_rep, nonpivots = Arref[:,:-1].rep.nullspace()
        nullspace = Arref.from_rep(nullspace_rep)
        return particular, nullspace

    def charpoly(self):
        r"""
        Returns the coefficients of the characteristic polynomial
        of the DomainMatrix. These elements will be domain elements.
        The domain of the elements will be same as domain of the DomainMatrix.

        Returns
        =======

        list
            coefficients of the characteristic polynomial

        Raises
        ======

        DMNonSquareMatrixError
            If the DomainMatrix is not a not Square DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)

        >>> A.charpoly()
        [1, -5, -2]

        """
        m, n = self.shape
        if m != n:
            raise DMNonSquareMatrixError("not square")
        return self.rep.charpoly()

    @classmethod
    def eye(cls, shape, domain):
        r"""
        Return identity matrix of size n

        Examples
        ========

        >>> from sympy.polys.matrices import DomainMatrix
        >>> from sympy import QQ
        >>> DomainMatrix.eye(3, QQ)
        DomainMatrix({0: {0: 1}, 1: {1: 1}, 2: {2: 1}}, (3, 3), QQ)

        """
        if isinstance(shape, int):
            shape = (shape, shape)
        return cls.from_rep(SDM.eye(shape, domain))

    @classmethod
    def diag(cls, diagonal, domain, shape=None):
        r"""
        Return diagonal matrix with entries from ``diagonal``.

        Examples
        ========

        >>> from sympy.polys.matrices import DomainMatrix
        >>> from sympy import ZZ
        >>> DomainMatrix.diag([ZZ(5), ZZ(6)], ZZ)
        DomainMatrix({0: {0: 5}, 1: {1: 6}}, (2, 2), ZZ)

        """
        if shape is None:
            N = len(diagonal)
            shape = (N, N)
        return cls.from_rep(SDM.diag(diagonal, domain, shape))

    @classmethod
    def zeros(cls, shape, domain, *, fmt='sparse'):
        """Returns a zero DomainMatrix of size shape, belonging to the specified domain

        Examples
        ========

        >>> from sympy.polys.matrices import DomainMatrix
        >>> from sympy import QQ
        >>> DomainMatrix.zeros((2, 3), QQ)
        DomainMatrix({}, (2, 3), QQ)

        """
        return cls.from_rep(SDM.zeros(shape, domain))

    @classmethod
    def ones(cls, shape, domain):
        """Returns a DomainMatrix of 1s, of size shape, belonging to the specified domain

        Examples
        ========

        >>> from sympy.polys.matrices import DomainMatrix
        >>> from sympy import QQ
        >>> DomainMatrix.ones((2,3), QQ)
        DomainMatrix([[1, 1, 1], [1, 1, 1]], (2, 3), QQ)

        """
        return cls.from_rep(DDM.ones(shape, domain))

    def __eq__(A, B):
        r"""
        Checks for two DomainMatrix matrices to be equal or not

        Parameters
        ==========

        A, B: DomainMatrix
            to check equality

        Returns
        =======

        Boolean
            True for equal, else False

        Raises
        ======

        NotImplementedError
            If B is not a DomainMatrix

        Examples
        ========

        >>> from sympy import ZZ
        >>> from sympy.polys.matrices import DomainMatrix
        >>> A = DomainMatrix([
        ...    [ZZ(1), ZZ(2)],
        ...    [ZZ(3), ZZ(4)]], (2, 2), ZZ)
        >>> B = DomainMatrix([
        ...    [ZZ(1), ZZ(1)],
        ...    [ZZ(0), ZZ(1)]], (2, 2), ZZ)
        >>> A.__eq__(A)
        True
        >>> A.__eq__(B)
        False

        """
        if not isinstance(A, type(B)):
            return NotImplemented
        return A.domain == B.domain and A.rep == B.rep

    def unify_eq(A, B):
        if A.shape != B.shape:
            return False
        if A.domain != B.domain:
            A, B = A.unify(B)
        return A == B
