from __future__ import print_function, division

import copy
from collections import defaultdict

from sympy.core.containers import Dict
from sympy.core.compatibility import is_sequence, as_int, range
from sympy.core.logic import fuzzy_and
from sympy.core.singleton import S
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.polys import cancel
from sympy.utilities.iterables import uniq

from .matrices import MatrixBase, ShapeError, a2idx, NonSquareMatrixError
from .dense import Matrix
import collections


class SparseMatrix(MatrixBase):
    """
    A sparse matrix (a matrix with a large number of zero elements).

    Examples
    ========

    >>> from sympy.matrices import SparseMatrix
    >>> SparseMatrix(2, 2, range(4))
    Matrix([
    [0, 1],
    [2, 3]])
    >>> SparseMatrix(2, 2, {(1, 1): 2})
    Matrix([
    [0, 0],
    [0, 2]])

    See Also
    ========
    sympy.matrices.dense.Matrix
    """

    def __init__(self, *args):

        if len(args) == 1 and isinstance(args[0], SparseMatrix):
            self.rows = args[0].rows
            self.cols = args[0].cols
            self._smat = dict(args[0]._smat)
            return

        self._smat = {}

        if len(args) == 3:
            self.rows = as_int(args[0])
            self.cols = as_int(args[1])

            if isinstance(args[2], collections.Callable):
                op = args[2]
                for i in range(self.rows):
                    for j in range(self.cols):
                        value = self._sympify(
                            op(self._sympify(i), self._sympify(j)))
                        if value:
                            self._smat[(i, j)] = value
            elif isinstance(args[2], (dict, Dict)):
                # manual copy, copy.deepcopy() doesn't work
                for key in args[2].keys():
                    v = args[2][key]
                    if v:
                        self._smat[key] = self._sympify(v)
            elif is_sequence(args[2]):
                if len(args[2]) != self.rows*self.cols:
                    raise ValueError(
                        'List length (%s) != rows*columns (%s)' %
                        (len(args[2]), self.rows*self.cols))
                flat_list = args[2]
                for i in range(self.rows):
                    for j in range(self.cols):
                        value = self._sympify(flat_list[i*self.cols + j])
                        if value:
                            self._smat[(i, j)] = value
        else:
            # handle full matrix forms with _handle_creation_inputs
            r, c, _list = Matrix._handle_creation_inputs(*args)
            self.rows = r
            self.cols = c
            for i in range(self.rows):
                for j in range(self.cols):
                    value = _list[self.cols*i + j]
                    if value:
                        self._smat[(i, j)] = value

    def __getitem__(self, key):

        if isinstance(key, tuple):
            i, j = key
            try:
                i, j = self.key2ij(key)
                return self._smat.get((i, j), S.Zero)
            except (TypeError, IndexError):
                if isinstance(i, slice):
                    # XXX remove list() when PY2 support is dropped
                    i = list(range(self.rows))[i]
                elif is_sequence(i):
                    pass
                else:
                    if i >= self.rows:
                        raise IndexError('Row index out of bounds')
                    i = [i]
                if isinstance(j, slice):
                    # XXX remove list() when PY2 support is dropped
                    j = list(range(self.cols))[j]
                elif is_sequence(j):
                    pass
                else:
                    if j >= self.cols:
                        raise IndexError('Col index out of bounds')
                    j = [j]
                return self.extract(i, j)

        # check for single arg, like M[:] or M[3]
        if isinstance(key, slice):
            lo, hi = key.indices(len(self))[:2]
            L = []
            for i in range(lo, hi):
                m, n = divmod(i, self.cols)
                L.append(self._smat.get((m, n), S.Zero))
            return L

        i, j = divmod(a2idx(key, len(self)), self.cols)
        return self._smat.get((i, j), S.Zero)

    def __setitem__(self, key, value):
        raise NotImplementedError()

    def copy(self):
        return self._new(self.rows, self.cols, self._smat)

    @property
    def is_Identity(self):
        if not self.is_square:
            return False
        if not all(self[i, i] == 1 for i in range(self.rows)):
            return False
        return len(self._smat) == self.rows

    def tolist(self):
        """Convert this sparse matrix into a list of nested Python lists.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix, ones
        >>> a = SparseMatrix(((1, 2), (3, 4)))
        >>> a.tolist()
        [[1, 2], [3, 4]]

        When there are no rows then it will not be possible to tell how
        many columns were in the original matrix:

        >>> SparseMatrix(ones(0, 3)).tolist()
        []

        """
        if not self.rows:
            return []
        if not self.cols:
            return [[] for i in range(self.rows)]
        I, J = self.shape
        return [[self[i, j] for j in range(J)] for i in range(I)]

    def row(self, i):
        """Returns column i from self as a row vector.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> a = SparseMatrix(((1, 2), (3, 4)))
        >>> a.row(0)
        Matrix([[1, 2]])

        See Also
        ========
        col
        row_list
        """
        return self[i,:]

    def col(self, j):
        """Returns column j from self as a column vector.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> a = SparseMatrix(((1, 2), (3, 4)))
        >>> a.col(0)
        Matrix([
        [1],
        [3]])

        See Also
        ========
        row
        col_list
        """
        return self[:, j]

    def row_list(self):
        """Returns a row-sorted list of non-zero elements of the matrix.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> a = SparseMatrix(((1, 2), (3, 4)))
        >>> a
        Matrix([
        [1, 2],
        [3, 4]])
        >>> a.RL
        [(0, 0, 1), (0, 1, 2), (1, 0, 3), (1, 1, 4)]

        See Also
        ========
        row_op
        col_list
        """
        return [tuple(k + (self[k],)) for k in
            sorted(list(self._smat.keys()), key=lambda k: list(k))]

    RL = property(row_list, None, None, "Alternate faster representation")

    def col_list(self):
        """Returns a column-sorted list of non-zero elements of the matrix.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> a=SparseMatrix(((1, 2), (3, 4)))
        >>> a
        Matrix([
        [1, 2],
        [3, 4]])
        >>> a.CL
        [(0, 0, 1), (1, 0, 3), (0, 1, 2), (1, 1, 4)]

        See Also
        ========
        col_op
        row_list
        """
        return [tuple(k + (self[k],)) for k in sorted(list(self._smat.keys()), key=lambda k: list(reversed(k)))]

    CL = property(col_list, None, None, "Alternate faster representation")

    def _eval_trace(self):
        """Calculate the trace of a square matrix.

        Examples
        ========

        >>> from sympy.matrices import eye
        >>> eye(3).trace()
        3

        """
        trace = S.Zero
        for i in range(self.cols):
            trace += self._smat.get((i, i), 0)
        return trace

    def _eval_transpose(self):
        """Returns the transposed SparseMatrix of this SparseMatrix.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> a = SparseMatrix(((1, 2), (3, 4)))
        >>> a
        Matrix([
        [1, 2],
        [3, 4]])
        >>> a.T
        Matrix([
        [1, 3],
        [2, 4]])
        """
        tran = self.zeros(self.cols, self.rows)
        for key, value in self._smat.items():
            key = key[1], key[0]  # reverse
            tran._smat[key] = value
        return tran

    def _eval_conjugate(self):
        """Return the by-element conjugation.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> from sympy import I
        >>> a = SparseMatrix(((1, 2 + I), (3, 4), (I, -I)))
        >>> a
        Matrix([
        [1, 2 + I],
        [3,     4],
        [I,    -I]])
        >>> a.C
        Matrix([
        [ 1, 2 - I],
        [ 3,     4],
        [-I,     I]])

        See Also
        ========

        transpose: Matrix transposition
        H: Hermite conjugation
        D: Dirac conjugation
        """
        conj = self.copy()
        for key, value in self._smat.items():
            conj._smat[key] = value.conjugate()
        return conj

    def multiply(self, other):
        """Fast multiplication exploiting the sparsity of the matrix.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix, ones
        >>> A, B = SparseMatrix(ones(4, 3)), SparseMatrix(ones(3, 4))
        >>> A.multiply(B) == 3*ones(4)
        True

        See Also
        ========

        add
        """
        A = self
        B = other
        # sort B's row_list into list of rows
        Blist = [[] for i in range(B.rows)]
        for i, j, v in B.row_list():
            Blist[i].append((j, v))
        Cdict = defaultdict(int)
        for k, j, Akj in A.row_list():
            for n, Bjn in Blist[j]:
                temp = Akj*Bjn
                Cdict[k, n] += temp
        rv = self.zeros(A.rows, B.cols)
        rv._smat = dict([(k, v) for k, v in Cdict.items() if v])
        return rv

    def scalar_multiply(self, scalar):
        "Scalar element-wise multiplication"
        M = self.zeros(*self.shape)
        if scalar:
            for i in self._smat:
                v = scalar*self._smat[i]
                if v:
                    M._smat[i] = v
                else:
                    M._smat.pop(i, None)
        return M

    def __mul__(self, other):
        """Multiply self and other, watching for non-matrix entities.

        When multiplying be a non-sparse matrix, the result is no longer
        sparse.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix, eye, zeros
        >>> I = SparseMatrix(eye(3))
        >>> I*I == I
        True
        >>> Z = zeros(3)
        >>> I*Z
        Matrix([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]])
        >>> I*2 == 2*I
        True
        """
        if isinstance(other, SparseMatrix):
            return self.multiply(other)
        if isinstance(other, MatrixBase):
            return other._new(self*self._new(other))
        return self.scalar_multiply(other)

    def __rmul__(self, other):
        """Return product the same type as other (if a Matrix).

        When multiplying be a non-sparse matrix, the result is no longer
        sparse.

        Examples
        ========

        >>> from sympy.matrices import Matrix, SparseMatrix
        >>> A = Matrix(2, 2, range(1, 5))
        >>> S = SparseMatrix(2, 2, range(2, 6))
        >>> A*S == S*A
        False
        >>> (isinstance(A*S, SparseMatrix) ==
        ...  isinstance(S*A, SparseMatrix) == False)
        True
        """
        if isinstance(other, MatrixBase):
            return other*other._new(self)
        return self.scalar_multiply(other)

    def __add__(self, other):
        """Add other to self, efficiently if possible.

        When adding a non-sparse matrix, the result is no longer
        sparse.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix, eye
        >>> A = SparseMatrix(eye(3)) + SparseMatrix(eye(3))
        >>> B = SparseMatrix(eye(3)) + eye(3)
        >>> A
        Matrix([
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 2]])
        >>> A == B
        True
        >>> isinstance(A, SparseMatrix) and isinstance(B, SparseMatrix)
        False

        """
        if isinstance(other, SparseMatrix):
            return self.add(other)
        elif isinstance(other, MatrixBase):
            return other._new(other + self)
        else:
            raise NotImplementedError(
                "Cannot add %s to %s" %
                tuple([c.__class__.__name__ for c in (other, self)]))

    def __neg__(self):
        """Negate all elements of self.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix, eye
        >>> -SparseMatrix(eye(3))
        Matrix([
        [-1,  0,  0],
        [ 0, -1,  0],
        [ 0,  0, -1]])

        """

        rv = self.copy()
        for k, v in rv._smat.items():
            rv._smat[k] = -v
        return rv

    def add(self, other):
        """Add two sparse matrices with dictionary representation.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix, eye, ones
        >>> SparseMatrix(eye(3)).add(SparseMatrix(ones(3)))
        Matrix([
        [2, 1, 1],
        [1, 2, 1],
        [1, 1, 2]])
        >>> SparseMatrix(eye(3)).add(-SparseMatrix(eye(3)))
        Matrix([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]])

        Only the non-zero elements are stored, so the resulting dictionary
        that is used to represent the sparse matrix is empty:
        >>> _._smat
        {}

        See Also
        ========

        multiply
        """
        if not isinstance(other, SparseMatrix):
            raise ValueError('only use add with %s, not %s' %
                tuple([c.__class__.__name__ for c in (self, other)]))
        if self.shape != other.shape:
            raise ShapeError()
        M = self.copy()
        for i, v in other._smat.items():
            v = M[i] + v
            if v:
                M._smat[i] = v
            else:
                M._smat.pop(i, None)
        return M

    def extract(self, rowsList, colsList):
        urow = list(uniq(rowsList))
        ucol = list(uniq(colsList))
        smat = {}
        if len(urow)*len(ucol) < len(self._smat):
            # there are fewer elements requested than there are elements in the matrix
            for i, r in enumerate(urow):
                for j, c in enumerate(ucol):
                    smat[i, j] = self._smat.get((r, c), 0)
        else:
            # most of the request will be zeros so check all of self's entries,
            # keeping only the ones that are desired
            for rk, ck in self._smat:
                if rk in urow and ck in ucol:
                    smat[(urow.index(rk), ucol.index(ck))] = self._smat[(rk, ck)]

        rv = self._new(len(urow), len(ucol), smat)
        # rv is nominally correct but there might be rows/cols
        # which require duplication
        if len(rowsList) != len(urow):
            for i, r in enumerate(rowsList):
                i_previous = rowsList.index(r)
                if i_previous != i:
                    rv = rv.row_insert(i, rv.row(i_previous))
        if len(colsList) != len(ucol):
            for i, c in enumerate(colsList):
                i_previous = colsList.index(c)
                if i_previous != i:
                    rv = rv.col_insert(i, rv.col(i_previous))
        return rv
    extract.__doc__ = MatrixBase.extract.__doc__

    @property
    def is_hermitian(self):
        """Checks if the matrix is Hermitian.

        In a Hermitian matrix element i,j is the complex conjugate of
        element j,i.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> from sympy import I
        >>> from sympy.abc import x
        >>> a = SparseMatrix([[1, I], [-I, 1]])
        >>> a
        Matrix([
        [ 1, I],
        [-I, 1]])
        >>> a.is_hermitian
        True
        >>> a[0, 0] = 2*I
        >>> a.is_hermitian
        False
        >>> a[0, 0] = x
        >>> a.is_hermitian
        >>> a[0, 1] = a[1, 0]*I
        >>> a.is_hermitian
        False
        """
        def cond():
            d = self._smat
            yield self.is_square
            if len(d) <= self.rows:
                yield fuzzy_and(
                    d[i, i].is_real for i, j in d if i == j)
            else:
                yield fuzzy_and(
                    d[i, i].is_real for i in range(self.rows) if (i, i) in d)
            yield fuzzy_and(
                    ((self[i, j] - self[j, i].conjugate()).is_zero
                    if (j, i) in d else False) for (i, j) in d)
        return fuzzy_and(i for i in cond())

    def is_symmetric(self, simplify=True):
        """Return True if self is symmetric.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix, eye
        >>> M = SparseMatrix(eye(3))
        >>> M.is_symmetric()
        True
        >>> M[0, 2] = 1
        >>> M.is_symmetric()
        False
        """
        if simplify:
            return all((k[1], k[0]) in self._smat and
                not (self[k] - self[(k[1], k[0])]).simplify()
                for k in self._smat)
        else:
            return all((k[1], k[0]) in self._smat and
                self[k] == self[(k[1], k[0])] for k in self._smat)

    def has(self, *patterns):
        """Test whether any subexpression matches any of the patterns.

        Examples
        ========

        >>> from sympy import SparseMatrix, Float
        >>> from sympy.abc import x, y
        >>> A = SparseMatrix(((1, x), (0.2, 3)))
        >>> A.has(x)
        True
        >>> A.has(y)
        False
        >>> A.has(Float)
        True
        """
        return any(self[key].has(*patterns) for key in self._smat)

    def applyfunc(self, f):
        """Apply a function to each element of the matrix.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> m = SparseMatrix(2, 2, lambda i, j: i*2+j)
        >>> m
        Matrix([
        [0, 1],
        [2, 3]])
        >>> m.applyfunc(lambda i: 2*i)
        Matrix([
        [0, 2],
        [4, 6]])

        """
        if not callable(f):
            raise TypeError("`f` must be callable.")

        out = self.copy()
        for k, v in self._smat.items():
            fv = f(v)
            if fv:
                out._smat[k] = fv
            else:
                out._smat.pop(k, None)
        return out

    def reshape(self, rows, cols):
        """Reshape matrix while retaining original size.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> S = SparseMatrix(4, 2, range(8))
        >>> S.reshape(2, 4)
        Matrix([
        [0, 1, 2, 3],
        [4, 5, 6, 7]])

        """
        if len(self) != rows*cols:
            raise ValueError("Invalid reshape parameters %d %d" % (rows, cols))
        smat = {}
        for k, v in self._smat.items():
            i, j = k
            n = i*self.cols + j
            ii, jj = divmod(n, cols)
            smat[(ii, jj)] = self._smat[(i, j)]
        return self._new(rows, cols, smat)

    def liupc(self):
        """Liu's algorithm, for pre-determination of the Elimination Tree of
        the given matrix, used in row-based symbolic Cholesky factorization.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> S = SparseMatrix([
        ... [1, 0, 3, 2],
        ... [0, 0, 1, 0],
        ... [4, 0, 0, 5],
        ... [0, 6, 7, 0]])
        >>> S.liupc()
        ([[0], [], [0], [1, 2]], [4, 3, 4, 4])

        References
        ==========

        Symbolic Sparse Cholesky Factorization using Elimination Trees,
        Jeroen Van Grondelle (1999)
        http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.39.7582
        """
        # Algorithm 2.4, p 17 of reference

        # get the indices of the elements that are non-zero on or below diag
        R = [[] for r in range(self.rows)]
        for r, c, _ in self.row_list():
            if c <= r:
                R[r].append(c)

        inf = len(R)  # nothing will be this large
        parent = [inf]*self.rows
        virtual = [inf]*self.rows
        for r in range(self.rows):
            for c in R[r][:-1]:
                while virtual[c] < r:
                    t = virtual[c]
                    virtual[c] = r
                    c = t
                if virtual[c] == inf:
                    parent[c] = virtual[c] = r
        return R, parent

    def row_structure_symbolic_cholesky(self):
        """Symbolic cholesky factorization, for pre-determination of the
        non-zero structure of the Cholesky factororization.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> S = SparseMatrix([
        ... [1, 0, 3, 2],
        ... [0, 0, 1, 0],
        ... [4, 0, 0, 5],
        ... [0, 6, 7, 0]])
        >>> S.row_structure_symbolic_cholesky()
        [[0], [], [0], [1, 2]]

        References
        ==========

        Symbolic Sparse Cholesky Factorization using Elimination Trees,
        Jeroen Van Grondelle (1999)
        http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.39.7582
        """

        R, parent = self.liupc()
        inf = len(R)  # this acts as infinity
        Lrow = copy.deepcopy(R)
        for k in range(self.rows):
            for j in R[k]:
                while j != inf and j != k:
                    Lrow[k].append(j)
                    j = parent[j]
            Lrow[k] = list(sorted(set(Lrow[k])))
        return Lrow

    def _cholesky_sparse(self):
        """Algorithm for numeric Cholesky factorization of a sparse matrix."""
        Crowstruc = self.row_structure_symbolic_cholesky()
        C = self.zeros(self.rows)
        for i in range(len(Crowstruc)):
            for j in Crowstruc[i]:
                if i != j:
                    C[i, j] = self[i, j]
                    summ = 0
                    for p1 in Crowstruc[i]:
                        if p1 < j:
                            for p2 in Crowstruc[j]:
                                if p2 < j:
                                    if p1 == p2:
                                        summ += C[i, p1]*C[j, p1]
                                else:
                                    break
                            else:
                                break
                    C[i, j] -= summ
                    C[i, j] /= C[j, j]
                else:
                    C[j, j] = self[j, j]
                    summ = 0
                    for k in Crowstruc[j]:
                        if k < j:
                            summ += C[j, k]**2
                        else:
                            break
                    C[j, j] -= summ
                    C[j, j] = sqrt(C[j, j])

        return C

    def _LDL_sparse(self):
        """Algorithm for numeric LDL factization, exploiting sparse structure.
        """
        Lrowstruc = self.row_structure_symbolic_cholesky()
        L = self.eye(self.rows)
        D = self.zeros(self.rows, self.cols)

        for i in range(len(Lrowstruc)):
            for j in Lrowstruc[i]:
                if i != j:
                    L[i, j] = self[i, j]
                    summ = 0
                    for p1 in Lrowstruc[i]:
                        if p1 < j:
                            for p2 in Lrowstruc[j]:
                                if p2 < j:
                                    if p1 == p2:
                                        summ += L[i, p1]*L[j, p1]*D[p1, p1]
                                else:
                                    break
                        else:
                            break
                    L[i, j] -= summ
                    L[i, j] /= D[j, j]
                elif i == j:
                    D[i, i] = self[i, i]
                    summ = 0
                    for k in Lrowstruc[i]:
                        if k < i:
                            summ += L[i, k]**2*D[k, k]
                        else:
                            break
                    D[i, i] -= summ

        return L, D

    def _lower_triangular_solve(self, rhs):
        """Fast algorithm for solving a lower-triangular system,
        exploiting the sparsity of the given matrix.
        """
        rows = [[] for i in range(self.rows)]
        for i, j, v in self.row_list():
            if i > j:
                rows[i].append((j, v))
        X = rhs.copy()
        for i in range(self.rows):
            for j, v in rows[i]:
                X[i, 0] -= v*X[j, 0]
            X[i, 0] /= self[i, i]
        return self._new(X)

    def _upper_triangular_solve(self, rhs):
        """Fast algorithm for solving an upper-triangular system,
        exploiting the sparsity of the given matrix.
        """
        rows = [[] for i in range(self.rows)]
        for i, j, v in self.row_list():
            if i < j:
                rows[i].append((j, v))
        X = rhs.copy()
        for i in range(self.rows - 1, -1, -1):
            rows[i].reverse()
            for j, v in rows[i]:
                X[i, 0] -= v*X[j, 0]
            X[i, 0] /= self[i, i]
        return self._new(X)

    def _diagonal_solve(self, rhs):
        "Diagonal solve."
        return self._new(self.rows, 1, lambda i, j: rhs[i, 0] / self[i, i])

    def _cholesky_solve(self, rhs):
        # for speed reasons, this is not uncommented, but if you are
        # having difficulties, try uncommenting to make sure that the
        # input matrix is symmetric

        #assert self.is_symmetric()
        L = self._cholesky_sparse()
        Y = L._lower_triangular_solve(rhs)
        rv = L.T._upper_triangular_solve(Y)
        return rv

    def _LDL_solve(self, rhs):
        # for speed reasons, this is not uncommented, but if you are
        # having difficulties, try uncommenting to make sure that the
        # input matrix is symmetric

        #assert self.is_symmetric()
        L, D = self._LDL_sparse()
        Z = L._lower_triangular_solve(rhs)
        Y = D._diagonal_solve(Z)
        return L.T._upper_triangular_solve(Y)

    def cholesky(self):
        """
        Returns the Cholesky decomposition L of a matrix A
        such that L * L.T = A

        A must be a square, symmetric, positive-definite
        and non-singular matrix

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> A = SparseMatrix(((25,15,-5),(15,18,0),(-5,0,11)))
        >>> A.cholesky()
        Matrix([
        [ 5, 0, 0],
        [ 3, 3, 0],
        [-1, 1, 3]])
        >>> A.cholesky() * A.cholesky().T == A
        True
        """

        from sympy.core.numbers import nan, oo
        if not self.is_symmetric():
            raise ValueError('Cholesky decomposition applies only to '
                'symmetric matrices.')
        M = self.as_mutable()._cholesky_sparse()
        if M.has(nan) or M.has(oo):
            raise ValueError('Cholesky decomposition applies only to '
                'positive-definite matrices')
        return self._new(M)

    def LDLdecomposition(self):
        """
        Returns the LDL Decomposition (matrices ``L`` and ``D``) of matrix
        ``A``, such that ``L * D * L.T == A``. ``A`` must be a square,
        symmetric, positive-definite and non-singular.

        This method eliminates the use of square root and ensures that all
        the diagonal entries of L are 1.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> A = SparseMatrix(((25, 15, -5), (15, 18, 0), (-5, 0, 11)))
        >>> L, D = A.LDLdecomposition()
        >>> L
        Matrix([
        [   1,   0, 0],
        [ 3/5,   1, 0],
        [-1/5, 1/3, 1]])
        >>> D
        Matrix([
        [25, 0, 0],
        [ 0, 9, 0],
        [ 0, 0, 9]])
        >>> L * D * L.T == A
        True

        """
        from sympy.core.numbers import nan, oo
        if not self.is_symmetric():
            raise ValueError('LDL decomposition applies only to '
                'symmetric matrices.')
        L, D = self.as_mutable()._LDL_sparse()
        if L.has(nan) or L.has(oo) or D.has(nan) or D.has(oo):
            raise ValueError('LDL decomposition applies only to '
                'positive-definite matrices')

        return self._new(L), self._new(D)

    def solve_least_squares(self, rhs, method='LDL'):
        """Return the least-square fit to the data.

        By default the cholesky_solve routine is used (method='CH'); other
        methods of matrix inversion can be used. To find out which are
        available, see the docstring of the .inv() method.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix, Matrix, ones
        >>> A = Matrix([1, 2, 3])
        >>> B = Matrix([2, 3, 4])
        >>> S = SparseMatrix(A.row_join(B))
        >>> S
        Matrix([
        [1, 2],
        [2, 3],
        [3, 4]])

        If each line of S represent coefficients of Ax + By
        and x and y are [2, 3] then S*xy is:

        >>> r = S*Matrix([2, 3]); r
        Matrix([
        [ 8],
        [13],
        [18]])

        But let's add 1 to the middle value and then solve for the
        least-squares value of xy:

        >>> xy = S.solve_least_squares(Matrix([8, 14, 18])); xy
        Matrix([
        [ 5/3],
        [10/3]])

        The error is given by S*xy - r:

        >>> S*xy - r
        Matrix([
        [1/3],
        [1/3],
        [1/3]])
        >>> _.norm().n(2)
        0.58

        If a different xy is used, the norm will be higher:

        >>> xy += ones(2, 1)/10
        >>> (S*xy - r).norm().n(2)
        1.5

        """
        t = self.T
        return (t*self).inv(method=method)*t*rhs

    def solve(self, rhs, method='LDL'):
        """Return solution to self*soln = rhs using given inversion method.

        For a list of possible inversion methods, see the .inv() docstring.
        """
        if not self.is_square:
            if self.rows < self.cols:
                raise ValueError('Under-determined system.')
            elif self.rows > self.cols:
                raise ValueError('For over-determined system, M, having '
                    'more rows than columns, try M.solve_least_squares(rhs).')
        else:
            return self.inv(method=method)*rhs

    def _eval_inverse(self, **kwargs):
        """Return the matrix inverse using Cholesky or LDL (default)
        decomposition as selected with the ``method`` keyword: 'CH' or 'LDL',
        respectively.

        Examples
        ========

        >>> from sympy import SparseMatrix, Matrix
        >>> A = SparseMatrix([
        ... [ 2, -1,  0],
        ... [-1,  2, -1],
        ... [ 0,  0,  2]])
        >>> A.inv('CH')
        Matrix([
        [2/3, 1/3, 1/6],
        [1/3, 2/3, 1/3],
        [  0,   0, 1/2]])
        >>> A.inv(method='LDL') # use of 'method=' is optional
        Matrix([
        [2/3, 1/3, 1/6],
        [1/3, 2/3, 1/3],
        [  0,   0, 1/2]])
        >>> A * _
        Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]])

        """
        sym = self.is_symmetric()
        M = self.as_mutable()
        I = M.eye(M.rows)
        if not sym:
            t = M.T
            r1 = M[0, :]
            M = t*M
            I = t*I
        method = kwargs.get('method', 'LDL')
        if method in "LDL":
            solve = M._LDL_solve
        elif method == "CH":
            solve = M._cholesky_solve
        else:
            raise NotImplementedError(
                'Method may be "CH" or "LDL", not %s.' % method)
        rv = M.hstack(*[solve(I[:, i]) for i in range(I.cols)])
        if not sym:
            scale = (r1*rv[:, 0])[0, 0]
            rv /= scale
        return self._new(rv)

    def __eq__(self, other):
        try:
            if self.shape != other.shape:
                return False
            if isinstance(other, SparseMatrix):
                return self._smat == other._smat
            elif isinstance(other, MatrixBase):
                return self._smat == MutableSparseMatrix(other)._smat
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self == other

    def as_mutable(self):
        """Returns a mutable version of this matrix.

        Examples
        ========

        >>> from sympy import ImmutableMatrix
        >>> X = ImmutableMatrix([[1, 2], [3, 4]])
        >>> Y = X.as_mutable()
        >>> Y[1, 1] = 5 # Can set values in Y
        >>> Y
        Matrix([
        [1, 2],
        [3, 5]])
        """
        return MutableSparseMatrix(self)

    def as_immutable(self):
        """Returns an Immutable version of this Matrix."""
        from .immutable import ImmutableSparseMatrix
        return ImmutableSparseMatrix(self)

    def nnz(self):
        """Returns the number of non-zero elements in Matrix."""
        return len(self._smat)

    @classmethod
    def zeros(cls, r, c=None):
        """Return an r x c matrix of zeros, square if c is omitted."""
        c = r if c is None else c
        r = as_int(r)
        c = as_int(c)
        return cls(r, c, {})

    @classmethod
    def eye(cls, n):
        """Return an n x n identity matrix."""
        n = as_int(n)
        return cls(n, n, dict([((i, i), S.One) for i in range(n)]))

    def det_bareis(self):
        """
        Compute matrix determinant using a variant of Bareiss's fraction-free
        algorithm. This variant is best suited for sparse matrices.
        (SFF = Sparse Fraction Free Algorithm)
        See http://www.eecis.udel.edu/~saunders/papers/sffge/it5.ps.
        Maybe implement the variant "division first" (see appendix) ?

        See Also
        ========

        det
        berkowitz_det
        """
        if not self.is_square:
            raise NonSquareMatrixError()
        if not self:
            return S.One

        M, n = self.copy().as_mutable(), self.rows

        if n == 1:
            det = M[0, 0]
        elif n == 2:
            det = M[0, 0]*M[1, 1] - M[0, 1]*M[1, 0]
        elif n == 3:
            det = (M[0, 0]*M[1, 1]*M[2, 2] + M[0, 1]*M[1, 2]*M[2, 0] + M[0, 2]*M[1, 0]*M[2, 1]) - \
                  (M[0, 2]*M[1, 1]*M[2, 0] + M[0, 0]*M[1, 2]*M[2, 1] + M[0, 1]*M[1, 0]*M[2, 2])
        else:
            sign = 1  # track current sign in case of column swap

            B = self.copy().as_mutable()  # data
            h = Matrix.zeros(n, n)  # history

            def P(m):  # See paper : P_m = M^(m)_(m+1, m+1) and P_(-1) = 1
                if m == -1:
                    return 1
                else:
                    return M[m, m]

            for k in range(n - 1):
                # look for a pivot in the current column
                # and assume det == 0 if none is found
                if M[k, k] == 0:
                    for i in range(k + 1, n):
                        if M[i, k]:
                            M.row_swap(i, k)
                            B.row_swap(i, k)
                            h.row_swap(i, k)
                            sign *= -1
                            break
                    else:
                        return S.Zero

                # proceed with SFF
                for i in range(k + 1, n):
                    for j in range(k + 1, n):
                        if B[i, k] * B[k, j] != 0:
                            B[i, j] = (P(k-1) * B[k, k] * B[i, j]) / (P(h[k, k] - 1) * P(h[i, j] - 1)) - \
                                      (P(k-1) * B[i, k] * B[k, j]) / (P(h[i, k] - 1) * P(h[k, j] - 1))
                            h[i, j] = k+1

                        D = B[i, j] * P(k) / P(h[i, j] - 1)
                        if D.is_Atom:
                            M[i, j] = D
                        else:
                            M[i, j] = cancel(D)

            det = sign*M[n - 1, n - 1]
        return det.expand()


class MutableSparseMatrix(SparseMatrix, MatrixBase):
    @classmethod
    def _new(cls, *args, **kwargs):
        return cls(*args)

    def as_mutable(self):
        return self.copy()

    def __setitem__(self, key, value):
        """Assign value to position designated by key.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix, ones
        >>> M = SparseMatrix(2, 2, {})
        >>> M[1] = 1; M
        Matrix([
        [0, 1],
        [0, 0]])
        >>> M[1, 1] = 2; M
        Matrix([
        [0, 1],
        [0, 2]])
        >>> M = SparseMatrix(2, 2, {})
        >>> M[:, 1] = [1, 1]; M
        Matrix([
        [0, 1],
        [0, 1]])
        >>> M = SparseMatrix(2, 2, {})
        >>> M[1, :] = [[1, 1]]; M
        Matrix([
        [0, 0],
        [1, 1]])


        To replace row r you assign to position r*m where m
        is the number of columns:

        >>> M = SparseMatrix(4, 4, {})
        >>> m = M.cols
        >>> M[3*m] = ones(1, m)*2; M
        Matrix([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [2, 2, 2, 2]])

        And to replace column c you can assign to position c:

        >>> M[2] = ones(m, 1)*4; M
        Matrix([
        [0, 0, 4, 0],
        [0, 0, 4, 0],
        [0, 0, 4, 0],
        [2, 2, 4, 2]])
        """
        rv = self._setitem(key, value)
        if rv is not None:
            i, j, value = rv
            if value:
                self._smat[(i, j)] = value
            elif (i, j) in self._smat:
                del self._smat[(i, j)]

    __hash__ = None

    def row_del(self, k):
        """Delete the given row of the matrix.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> M = SparseMatrix([[0, 0], [0, 1]])
        >>> M
        Matrix([
        [0, 0],
        [0, 1]])
        >>> M.row_del(0)
        >>> M
        Matrix([[0, 1]])

        See Also
        ========

        col_del
        """
        newD = {}
        k = a2idx(k, self.rows)
        for (i, j) in self._smat:
            if i == k:
                pass
            elif i > k:
                newD[i - 1, j] = self._smat[i, j]
            else:
                newD[i, j] = self._smat[i, j]
        self._smat = newD
        self.rows -= 1

    def col_del(self, k):
        """Delete the given column of the matrix.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> M = SparseMatrix([[0, 0], [0, 1]])
        >>> M
        Matrix([
        [0, 0],
        [0, 1]])
        >>> M.col_del(0)
        >>> M
        Matrix([
        [0],
        [1]])

        See Also
        ========

        row_del
        """
        newD = {}
        k = a2idx(k, self.cols)
        for (i, j) in self._smat:
            if j == k:
                pass
            elif j > k:
                newD[i, j - 1] = self._smat[i, j]
            else:
                newD[i, j] = self._smat[i, j]
        self._smat = newD
        self.cols -= 1

    def row_swap(self, i, j):
        """Swap, in place, columns i and j.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> S = SparseMatrix.eye(3); S[2, 1] = 2
        >>> S.row_swap(1, 0); S
        Matrix([
        [0, 1, 0],
        [1, 0, 0],
        [0, 2, 1]])
        """
        if i > j:
            i, j = j, i
        rows = self.row_list()
        temp = []
        for ii, jj, v in rows:
            if ii == i:
                self._smat.pop((ii, jj))
                temp.append((jj, v))
            elif ii == j:
                self._smat.pop((ii, jj))
                self._smat[i, jj] = v
            elif ii > j:
                break
        for k, v in temp:
            self._smat[j, k] = v

    def col_swap(self, i, j):
        """Swap, in place, columns i and j.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> S = SparseMatrix.eye(3); S[2, 1] = 2
        >>> S.col_swap(1, 0); S
        Matrix([
        [0, 1, 0],
        [1, 0, 0],
        [2, 0, 1]])
        """
        if i > j:
            i, j = j, i
        rows = self.col_list()
        temp = []
        for ii, jj, v in rows:
            if jj == i:
                self._smat.pop((ii, jj))
                temp.append((ii, v))
            elif jj == j:
                self._smat.pop((ii, jj))
                self._smat[ii, i] = v
            elif jj > j:
                break
        for k, v in temp:
            self._smat[k, j] = v

    def row_join(self, other):
        """Returns B appended after A (column-wise augmenting)::

            [A B]

        Examples
        ========

        >>> from sympy import SparseMatrix, Matrix
        >>> A = SparseMatrix(((1, 0, 1), (0, 1, 0), (1, 1, 0)))
        >>> A
        Matrix([
        [1, 0, 1],
        [0, 1, 0],
        [1, 1, 0]])
        >>> B = SparseMatrix(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
        >>> B
        Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]])
        >>> C = A.row_join(B); C
        Matrix([
        [1, 0, 1, 1, 0, 0],
        [0, 1, 0, 0, 1, 0],
        [1, 1, 0, 0, 0, 1]])
        >>> C == A.row_join(Matrix(B))
        True

        Joining at row ends is the same as appending columns at the end
        of the matrix:

        >>> C == A.col_insert(A.cols, B)
        True
        """
        A, B = self, other
        if not A.rows == B.rows:
            raise ShapeError()
        A = A.copy()
        if not isinstance(B, SparseMatrix):
            k = 0
            b = B._mat
            for i in range(B.rows):
                for j in range(B.cols):
                    v = b[k]
                    if v:
                        A._smat[(i, j + A.cols)] = v
                    k += 1
        else:
            for (i, j), v in B._smat.items():
                A._smat[(i, j + A.cols)] = v
        A.cols += B.cols
        return A

    def col_join(self, other):
        """Returns B augmented beneath A (row-wise joining)::

            [A]
            [B]

        Examples
        ========

        >>> from sympy import SparseMatrix, Matrix, ones
        >>> A = SparseMatrix(ones(3))
        >>> A
        Matrix([
        [1, 1, 1],
        [1, 1, 1],
        [1, 1, 1]])
        >>> B = SparseMatrix.eye(3)
        >>> B
        Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]])
        >>> C = A.col_join(B); C
        Matrix([
        [1, 1, 1],
        [1, 1, 1],
        [1, 1, 1],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]])
        >>> C == A.col_join(Matrix(B))
        True

        Joining along columns is the same as appending rows at the end
        of the matrix:

        >>> C == A.row_insert(A.rows, Matrix(B))
        True
        """
        A, B = self, other
        if not A.cols == B.cols:
            raise ShapeError()
        A = A.copy()
        if not isinstance(B, SparseMatrix):
            k = 0
            b = B._mat
            for i in range(B.rows):
                for j in range(B.cols):
                    v = b[k]
                    if v:
                        A._smat[(i + A.rows, j)] = v
                    k += 1
        else:
            for (i, j), v in B._smat.items():
                A._smat[i + A.rows, j] = v
        A.rows += B.rows
        return A

    def copyin_list(self, key, value):
        if not is_sequence(value):
            raise TypeError("`value` must be of type list or tuple.")
        self.copyin_matrix(key, Matrix(value))

    def copyin_matrix(self, key, value):
        # include this here because it's not part of BaseMatrix
        rlo, rhi, clo, chi = self.key2bounds(key)
        shape = value.shape
        dr, dc = rhi - rlo, chi - clo
        if shape != (dr, dc):
            raise ShapeError(
                "The Matrix `value` doesn't have the same dimensions "
                "as the in sub-Matrix given by `key`.")
        if not isinstance(value, SparseMatrix):
            for i in range(value.rows):
                for j in range(value.cols):
                    self[i + rlo, j + clo] = value[i, j]
        else:
            if (rhi - rlo)*(chi - clo) < len(self):
                for i in range(rlo, rhi):
                    for j in range(clo, chi):
                        self._smat.pop((i, j), None)
            else:
                for i, j, v in self.row_list():
                    if rlo <= i < rhi and clo <= j < chi:
                        self._smat.pop((i, j), None)
            for k, v in value._smat.items():
                i, j = k
                self[i + rlo, j + clo] = value[i, j]

    def zip_row_op(self, i, k, f):
        """In-place operation on row ``i`` using two-arg functor whose args are
        interpreted as ``(self[i, j], self[k, j])``.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> M = SparseMatrix.eye(3)*2
        >>> M[0, 1] = -1
        >>> M.zip_row_op(1, 0, lambda v, u: v + 2*u); M
        Matrix([
        [2, -1, 0],
        [4,  0, 0],
        [0,  0, 2]])

        See Also
        ========
        row
        row_op
        col_op

        """
        self.row_op(i, lambda v, j: f(v, self[k, j]))

    def row_op(self, i, f):
        """In-place operation on row ``i`` using two-arg functor whose args are
        interpreted as ``(self[i, j], j)``.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> M = SparseMatrix.eye(3)*2
        >>> M[0, 1] = -1
        >>> M.row_op(1, lambda v, j: v + 2*M[0, j]); M
        Matrix([
        [2, -1, 0],
        [4,  0, 0],
        [0,  0, 2]])

        See Also
        ========
        row
        zip_row_op
        col_op

        """
        for j in range(self.cols):
            v = self._smat.get((i, j), S.Zero)
            fv = f(v, j)
            if fv:
                self._smat[(i, j)] = fv
            elif v:
                self._smat.pop((i, j))

    def col_op(self, j, f):
        """In-place operation on col j using two-arg functor whose args are
        interpreted as (self[i, j], i) for i in range(self.rows).

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> M = SparseMatrix.eye(3)*2
        >>> M[1, 0] = -1
        >>> M.col_op(1, lambda v, i: v + 2*M[i, 0]); M
        Matrix([
        [ 2, 4, 0],
        [-1, 0, 0],
        [ 0, 0, 2]])
        """
        for i in range(self.rows):
            v = self._smat.get((i, j), S.Zero)
            fv = f(v, i)
            if fv:
                self._smat[(i, j)] = fv
            elif v:
                self._smat.pop((i, j))

    def fill(self, value):
        """Fill self with the given value.

        Notes
        =====

        Unless many values are going to be deleted (i.e. set to zero)
        this will create a matrix that is slower than a dense matrix in
        operations.

        Examples
        ========

        >>> from sympy.matrices import SparseMatrix
        >>> M = SparseMatrix.zeros(3); M
        Matrix([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]])
        >>> M.fill(1); M
        Matrix([
        [1, 1, 1],
        [1, 1, 1],
        [1, 1, 1]])
        """
        if not value:
            self._smat = {}
        else:
            v = self._sympify(value)
            self._smat = dict([((i, j), v)
                for i in range(self.rows) for j in range(self.cols)])
