"""

Module for the SDM class.

"""

from operator import add, neg, pos, sub
from collections import defaultdict

from .exceptions import DDMBadInputError, DDMDomainError, DDMShapeError

from .ddm import DDM


class SDM(dict):
    """Sparse matrix based on polys domain elements

    This is a dict subclass and is a wrapper for a dict of dicts that supports
    basic matrix arithmetic +, -, *, **.
    """

    fmt = 'sparse'

    def __init__(self, elemsdict, shape, domain):
        super().__init__(elemsdict)
        self.shape = self.rows, self.cols = m, n = shape
        self.domain = domain

        if not all(0 <= r < m for r in self):
            raise DDMBadInputError("Row out of range")
        if not all(0 <= c < n for row in self.values() for c in row):
            raise DDMBadInputError("Column out of range")

    def __str__(self):
        rowsstr = []
        for i, row in self.items():
            elemsstr = ', '.join('%s: %s' % (j, elem) for j, elem in row.items())
            rowsstr.append('%s: {%s}' % (i, elemsstr))
        return '{%s}' % ', '.join(rowsstr)

    def __repr__(self):
        cls = type(self).__name__
        rows = dict.__repr__(self)
        return '%s(%s, %s, %s)' % (cls, rows, self.shape, self.domain)

    @classmethod
    def new(cls, sdm, shape, domain):
        return cls(sdm, shape, domain)

    def copy(A):
        Ac = {i: Ai.copy() for i, Ai in A.items()}
        return A.new(Ac, A.shape, A.domain)

    @classmethod
    def from_list(cls, ddm, shape, domain):
        m, n = shape
        if not (len(ddm) == m and all(len(row) == n for row in ddm)):
            raise DDMBadInputError("Inconsistent row-list/shape")
        getrow = lambda i: {j:ddm[i][j] for j in range(n) if ddm[i][j]}
        irows = ((i, getrow(i)) for i in range(m))
        sdm = {i: row for i, row in irows if row}
        return cls(sdm, shape, domain)

    @classmethod
    def from_ddm(cls, ddm):
        return cls.from_list(ddm, ddm.shape, ddm.domain)

    def to_list(M):
        m, n = M.shape
        zero = M.domain.zero
        ddm = [[zero] * n for _ in range(m)]
        for i, row in M.items():
            for j, e in row.items():
                ddm[i][j] = e
        return ddm

    def to_ddm(M):
        return DDM(M.to_list(), M.shape, M.domain)

    def to_sdm(M):
        return M

    @classmethod
    def zeros(cls, shape, domain):
        return cls({}, shape, domain)

    @classmethod
    def eye(cls, size, domain):
        one = domain.one
        sdm = {i: {i: one} for i in range(size)}
        return cls(sdm, (size, size), domain)

    def transpose(M):
        MT = sdm_transpose(M)
        return M.new(MT, M.shape[::-1], M.domain)

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

    def matmul(A, B):
        if A.domain != B.domain:
            raise DDMDomainError
        m, n = A.shape
        n2, o = B.shape
        if n != n2:
            raise DDMShapeError
        C = sdm_matmul(A, B)
        return A.new(C, (m, o), A.domain)

    def mul(A, b):
        Csdm = unop_dict(A, lambda aij: aij*b)
        return A.new(Csdm, A.shape, A.domain)

    def add(A, B):
        Csdm = binop_dict(A, B, add, pos, pos)
        return A.new(Csdm, A.shape, A.domain)

    def sub(A, B):
        Csdm = binop_dict(A, B, sub, pos, neg)
        return A.new(Csdm, A.shape, A.domain)

    def neg(A):
        Csdm = unop_dict(A, neg)
        return A.new(Csdm, A.shape, A.domain)

    def convert_to(A, K):
        Kold = A.domain
        if K == Kold:
            return A.copy()
        Ak = unop_dict(A, lambda e: K.convert_from(e, Kold))
        return A.new(Ak, A.shape, K)

    def rref(A):
        B, pivots, _ = sdm_irref(A)
        return A.new(B, A.shape, A.domain), pivots

    def inv(A):
        return A.from_ddm(A.to_ddm().inv())

    def det(A):
        return A.to_ddm().det()

    def lu(A):
        L, U, swaps = A.to_ddm().lu()
        return A.from_ddm(L), A.from_ddm(U), swaps

    def lu_solve(A, b):
        return A.from_ddm(A.to_ddm().lu_solve(b.to_ddm()))

    def nullspace(A):
        ncols = A.shape[1]
        one = A.domain.one
        B, pivots, nzcols = sdm_irref(A)
        K, nonpivots = sdm_nullspace_from_rref(B, one, ncols, pivots, nzcols)
        K = dict(enumerate(K))
        shape = (len(K), ncols)
        return A.new(K, shape, A.domain), nonpivots

    def particular(A):
        ncols = A.shape[1]
        B, pivots, nzcols = sdm_irref(A)
        P = sdm_particular_from_rref(B, ncols, pivots)
        rep = {0:P} if P else {}
        return A.new(rep, (1, A.shape[1]), A.domain)

    def hstack(A, *B):
        Anew = dict(A.copy())
        rows, cols = A.shape
        domain = A.domain

        for Bk in B:
            Bkrows, Bkcols = Bk.shape
            assert Bkrows == rows
            assert Bk.domain == domain

            for i, Bki in Bk.items():
                Ai = Anew.get(i, None)
                if Ai is None:
                    Anew[i] = Ai = {}
                for j, Bkij in Bki.items():
                    Ai[j + cols] = Bkij
            cols += Bkcols

        return A.new(Anew, (rows, cols), A.domain)

    def charpoly(A):
        return A.to_ddm().charpoly()


def binop_dict(A, B, fab, fa, fb):
    Anz, Bnz = set(A), set(B)
    C = {}
    for i in Anz & Bnz:
        Ai, Bi = A[i], B[i]
        Ci = {}
        Anzi, Bnzi = set(Ai), set(Bi)
        for j in Anzi & Bnzi:
            elem = fab(Ai[j], Bi[j])
            if elem:
                Ci[j] = elem
        for j in Anzi - Bnzi:
            Ci[j] = fa(Ai[j])
        for j in Bnzi - Anzi:
            Ci[j] = fb(Bi[j])
        if Ci:
            C[i] = Ci
    for i in Anz - Bnz:
        Ai = A[i]
        C[i] = {j: fa(Aij) for j, Aij in Ai.items()}
    for i in Bnz - Anz:
        Bi = B[i]
        C[i] = {j: fb(Bij) for j, Bij in Bi.items()}
    return C


def unop_dict(A, f):
    B = {}
    for i, Ai in A.items():
        Bi = {}
        for j, Aij in Ai.items():
            Bij = f(Aij)
            if Bij:
                Bi[j] = Bij
        if Bi:
            B[i] = Bi
    return B


def sdm_transpose(M):
    MT = {}
    for i, Mi in M.items():
        for j, Mij in Mi.items():
            try:
                MT[j][i] = Mij
            except KeyError:
                MT[j] = {i: Mij}
    return MT


def sdm_matmul(A, B):
    #
    # Should be fast if A and B are very sparse.
    # Consider e.g. A = B = eye(1000).
    #
    # The idea here is that we compute C = A*B in terms of the rows of C and
    # B since the dict of dicts representation naturally stores the matrix as
    # rows. The ith row of C (Ci) is equal to the sum of Aik * Bk where Bk is
    # the kth row of B. The algorithm below loops over each nonzero element
    # Aik of A and if the corresponding row Bj is nonzero then we do
    #    Ci += Aik * Bk.
    # To make this more efficient we don't need to loop over all elements Aik.
    # Instead for each row Ai we compute the intersection of the nonzero
    # columns in Ai with the nonzero rows in B. That gives the k such that
    # Aik and Bk are both nonzero. In Python the intersection of two sets
    # of int can be computed very efficiently.
    #
    C = {}
    B_knz = set(B)
    for i, Ai in A.items():
        Ci = {}
        Ai_knz = set(Ai)
        for k in Ai_knz & B_knz:
            Aik = Ai[k]
            for j, Bkj in B[k].items():
                Cij = Ci.get(j, None)
                if Cij is not None:
                    Cij = Cij + Aik * Bkj
                    if Cij:
                        Ci[j] = Cij
                    else:
                        Ci.pop(j)
                else:
                    Cij = Aik * Bkj
                    if Cij:
                        Ci[j] = Cij
        if Ci:
            C[i] = Ci
    return C


def sdm_irref(A):
    """RREF and pivots of a sparse matrix *A*.

    Compute the reduced row echelon form (RREF) of the matrix *A* and return a
    list of the pivot columns. This routine does not work in place and leaves
    the original matrix *A* unmodified.

    Examples
    ========

    This routine works with a dict of dicts sparse representation of a matrix:

    >>> from sympy import QQ
    >>> from sympy.polys.matrices.sdm import sdm_irref
    >>> A = {0: {0: QQ(1), 1: QQ(2)}, 1: {0: QQ(3), 1: QQ(4)}}
    >>> Arref, pivots, _ = sdm_irref(A)
    >>> Arref
    {0: {0: 1}, 1: {1: 1}}
    >>> pivots
    [0, 1]

    The analogous calculation with :py:class:`~.Matrix` would be

    >>> from sympy import Matrix
    >>> M = Matrix([[1, 2], [3, 4]])
    >>> Mrref, pivots = M.rref()
    >>> Mrref
    Matrix([
    [1, 0],
    [0, 1]])
    >>> pivots
    (0, 1)

    Notes
    =====

    The cost of this algorithm is determined purely by the nonzero elements of
    the matrix. No part of the cost of any step in this algorithm depends on
    the number of rows or columns in the matrix. No step depends even on the
    number of nonzero rows apart from the primary loop over those rows. The
    implementation is much faster than ddm_rref for sparse matrices. In fact
    at the time of writing it is also (slightly) faster than the dense
    implementation even if the input is a fully dense matrix so it seems to be
    faster in all cases.

    The elements of the matrix should support exact division with ``/``. For
    example elements of any domain that is a field (e.g. ``QQ``) should be
    fine. No attempt is made to handle inexact arithmetic.

    """
    #
    # Any zeros in the matrix are not stored at all so an element is zero if
    # its row dict has no index at that key. A row is entirely zero if its
    # row index is not in the outer dict. Since rref reorders the rows and
    # removes zero rows we can completely discard the row indices. The first
    # step then copies the row dicts into a list sorted by the index of the
    # first nonzero column in each row.
    #
    # The algorithm then processes each row Ai one at a time. Previously seen
    # rows are used to cancel their pivot columns from Ai. Then a pivot from
    # Ai is chosen and is cancelled from all previously seen rows. At this
    # point Ai joins the previously seen rows. Once all rows are seen all
    # elimination has occurred and the rows are sorted by pivot column index.
    #
    # The previously seen rows are stored in two separate groups. The reduced
    # group consists of all rows that have been reduced to a single nonzero
    # element (the pivot). There is no need to attempt any further reduction
    # with these. Rows that still have other nonzeros need to be considered
    # when Ai is cancelled from the previously seen rows.
    #
    # A dict nonzerocolumns is used to map from a column index to a set of
    # previously seen rows that still have a nonzero element in that column.
    # This means that we can cancel the pivot from Ai into the previously seen
    # rows without needing to loop over each row that might have a zero in
    # that column.
    #

    # Row dicts sorted by index of first nonzero column
    # (Maybe sorting is not needed/useful.)
    Arows = sorted((Ai.copy() for Ai in A.values()), key=min)

    # Each processed row has an associated pivot column.
    # pivot_row_map maps from the pivot column index to the row dict.
    # This means that we can represent a set of rows purely as a set of their
    # pivot indices.
    pivot_row_map = {}

    # Set of pivot indices for rows that are fully reduced to a single nonzero.
    reduced_pivots = set()

    # Set of pivot indices for rows not fully reduced
    nonreduced_pivots = set()

    # Map from column index to a set of pivot indices representing the rows
    # that have a nonzero at that column.
    nonzero_columns = defaultdict(set)

    while Arows:
        # Select pivot element and row
        Ai = Arows.pop()

        # Nonzero columns from fully reduced pivot rows can be removed
        Ai = {j: Aij for j, Aij in Ai.items() if j not in reduced_pivots}

        # Others require full row cancellation
        for j in nonreduced_pivots & set(Ai):
            Aj = pivot_row_map[j]
            Aij = Ai[j]
            Ainz = set(Ai)
            Ajnz = set(Aj)
            for k in Ajnz - Ainz:
                Ai[k] = - Aij * Aj[k]
            for k in Ajnz & Ainz:
                Aik = Ai[k] - Aij * Aj[k]
                if Aik:
                    Ai[k] = Aik
                else:
                    Ai.pop(k)

        # We have now cancelled previously seen pivots from Ai.
        # If it is zero then discard it.
        if not Ai:
            continue

        # Choose a pivot from Ai:
        j = min(Ai)
        Aij = Ai[j]
        pivot_row_map[j] = Ai
        Ainz = set(Ai)

        # Normalise the pivot row to make the pivot 1.
        #
        # This approach is slow for some domains. Cross cancellation might be
        # better for e.g. QQ(x) with division delayed to the final steps.
        Aijinv = Aij**-1
        for l in Ai:
            Ai[l] *= Aijinv

        # Use Aij to cancel column j from all previously seen rows
        for k in nonzero_columns.pop(j, ()):
            Ak = pivot_row_map[k]
            Akj = Ak[j]
            Aknz = set(Ak)
            for l in Ainz - Aknz:
                Ak[l] = - Akj * Ai[l]
                nonzero_columns[l].add(k)
            for l in Ainz & Aknz:
                Akl = Ak[l] - Akj * Ai[l]
                if Akl:
                    Ak[l] = Akl
                else:
                    # Drop nonzero elements
                    Ak.pop(l)
                    if l != j:
                        nonzero_columns[l].remove(k)
            if len(Ak) == 1:
                reduced_pivots.add(k)
                nonreduced_pivots.remove(k)

        if len(Ai) == 1:
            reduced_pivots.add(j)
        else:
            nonreduced_pivots.add(j)
            for l in Ai:
                if l != j:
                    nonzero_columns[l].add(j)

    # All done!
    pivots = sorted(reduced_pivots | nonreduced_pivots)
    pivot2row = {p: n for n, p in enumerate(pivots)}
    nonzero_columns = {c: set(pivot2row[p] for p in s) for c, s in nonzero_columns.items()}
    rows = [pivot_row_map[i] for i in pivots]
    rref = dict(enumerate(rows))
    return rref, pivots, nonzero_columns


def sdm_nullspace_from_rref(A, one, ncols, pivots, nonzero_cols):
    """Get nullspace from A which is in RREF"""
    nonpivots = sorted(set(range(ncols)) - set(pivots))

    K = []
    for j in nonpivots:
        Kj = {j:one}
        for i in nonzero_cols.get(j, ()):
            Kj[pivots[i]] = -A[i][j]
        K.append(Kj)

    return K, nonpivots


def sdm_particular_from_rref(A, ncols, pivots):
    """Get a particular solution from A which is in RREF"""
    P = {}
    for i, j in enumerate(pivots):
        Ain = A[i].get(ncols-1, None)
        if Ain is not None:
            P[j] = Ain / A[i][j]
    return P
