from operator import index as index_

from sympy.core.compatibility import is_sequence
from sympy.core.expr import Expr
from sympy.core.kind import NumberKind, UndefinedKind
from sympy.core.sympify import _sympify, SympifyError
from sympy.polys.domains import ZZ, QQ, EXRAW
from sympy.polys.matrices import DomainMatrix
from sympy.core.singleton import S

from .common import classof
from .matrices import MatrixBase, MatrixKind


class RepMatrix(MatrixBase):

    def __eq__(self, other):
        # Skip sympify for mutable matrices...
        if not isinstance(other, RepMatrix):
            try:
                other = _sympify(other)
            except SympifyError:
                return NotImplemented
            if not isinstance(other, RepMatrix):
                return NotImplemented

        return self._rep.unify_eq(other._rep)

    @property
    def _flat(self):
        return self._rep.to_sympy().to_list_flat()

    def _eval_tolist(self):
        return self._rep.to_sympy().to_list()

    def _eval_todok(self):
        return self._rep.to_sympy().to_dok()

    def _eval_values(self):
        return list(self.todok().values())

    def copy(self):
        return self._fromrep(self._rep.copy())

    @property
    def kind(self):
        domain = self._rep.domain
        if domain in (ZZ, QQ):
            element_kind = NumberKind
        elif domain == EXRAW:
            kinds = set(e.kind for e in self.values())
            if len(kinds) == 1:
                [element_kind] = kinds
            else:
                element_kind = UndefinedKind
        else: # pragma: no cover
            raise RuntimeError("Domain should only be ZZ, QQ or EXRAW")
        return MatrixKind(element_kind)

    def _eval_has(self, *patterns):
        # if the matrix has any zeros, see if S.Zero
        # has the pattern.  If _smat is full length,
        # the matrix has no zeros.
        zhas = False
        if len(self.todok()) != self.rows*self.cols:
            zhas = S.Zero.has(*patterns)
        return zhas or any(value.has(*patterns) for value in self.values())

    def _eval_is_Identity(self):
        if not all(self[i, i] == 1 for i in range(self.rows)):
            return False
        return len(self.todok()) == self.rows

    def _eval_is_symmetric(self, simpfunc):
        diff = (self - self.T).applyfunc(simpfunc)
        return len(diff.values()) == 0

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
        return self._fromrep(self._rep.transpose())

    def _eval_extract(self, rowsList, colsList):
        return self._fromrep(self._rep.extract(rowsList, colsList))

    def __getitem__(self, key):
        return _getitem_RepMatrix(self, key)

    @classmethod
    def _eval_zeros(cls, rows, cols):
        rep = DomainMatrix.zeros((rows, cols), ZZ)
        return cls._fromrep(rep)

    @classmethod
    def _eval_eye(cls, rows, cols):
        rep = DomainMatrix.eye((rows, cols), ZZ)
        return cls._fromrep(rep)

    def _eval_add(self, other):
        return classof(self, other)._fromrep(self._rep + other._rep)

    def _eval_matrix_mul(self, other):
        return classof(self, other)._fromrep(self._rep * other._rep)

    def _eval_matrix_mul_elementwise(self, other):
        rep = self._rep.mul_elementwise(other._rep)
        return classof(self, other)._fromrep(rep)

    def _eval_scalar_mul(self, other):
        rep, other = self._unify_element_sympy(self._rep, other)
        return self._fromrep(rep.scalarmul(other))

    def _eval_scalar_rmul(self, other):
        rep, other = self._unify_element_sympy(self._rep, other)
        return self._fromrep(rep.rscalarmul(other))

    def _eval_Abs(self):
        return self._fromrep(self._rep.applyfunc(abs))

    def _eval_conjugate(self):
        rep = self._rep
        domain = rep.domain
        if domain in (ZZ, QQ):
            return self.copy()
        else:
            return self._fromrep(rep.applyfunc(lambda e: e.conjugate()))

    def equals(self, other, failing_expression=False):
        """Applies ``equals`` to corresponding elements of the matrices,
        trying to prove that the elements are equivalent, returning True
        if they are, False if any pair is not, and None (or the first
        failing expression if failing_expression is True) if it cannot
        be decided if the expressions are equivalent or not. This is, in
        general, an expensive operation.

        Examples
        ========

        >>> from sympy.matrices import Matrix
        >>> from sympy.abc import x
        >>> A = Matrix([x*(x - 1), 0])
        >>> B = Matrix([x**2 - x, 0])
        >>> A == B
        False
        >>> A.simplify() == B.simplify()
        True
        >>> A.equals(B)
        True
        >>> A.equals(2)
        False

        See Also
        ========
        sympy.core.expr.Expr.equals
        """
        if self.shape != getattr(other, 'shape', None):
            return False

        rv = True
        for i in range(self.rows):
            for j in range(self.cols):
                ans = self[i, j].equals(other[i, j], failing_expression)
                if ans is False:
                    return False
                elif ans is not True and rv is True:
                    rv = ans
        return rv


class MutableRepMatrix(RepMatrix):

    def __setitem__(self, key, value):
        """

        Examples
        ========

        >>> from sympy import Matrix, I, zeros, ones
        >>> m = Matrix(((1, 2+I), (3, 4)))
        >>> m
        Matrix([
        [1, 2 + I],
        [3,     4]])
        >>> m[1, 0] = 9
        >>> m
        Matrix([
        [1, 2 + I],
        [9,     4]])
        >>> m[1, 0] = [[0, 1]]

        To replace row r you assign to position r*m where m
        is the number of columns:

        >>> M = zeros(4)
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
            self._rep, value = self._unify_element_sympy(self._rep, value)
            self._rep.rep.setitem(i, j, value)


def _getitem_RepMatrix(self, key):
    """Return portion of self defined by key. If the key involves a slice
    then a list will be returned (if key is a single slice) or a matrix
    (if key was a tuple involving a slice).

    Examples
    ========

    >>> from sympy import Matrix, I
    >>> m = Matrix([
    ... [1, 2 + I],
    ... [3, 4    ]])

    If the key is a tuple that doesn't involve a slice then that element
    is returned:

    >>> m[1, 0]
    3

    When a tuple key involves a slice, a matrix is returned. Here, the
    first column is selected (all rows, column 0):

    >>> m[:, 0]
    Matrix([
    [1],
    [3]])

    If the slice is not a tuple then it selects from the underlying
    list of elements that are arranged in row order and a list is
    returned if a slice is involved:

    >>> m[0]
    1
    >>> m[::2]
    [1, 3]
    """
    if isinstance(key, tuple):
        i, j = key
        try:
            i, j = self.key2ij(key)
            return self._rep.getitem_sympy(i, j)
        except (TypeError, IndexError):
            if (isinstance(i, Expr) and not i.is_number) or (isinstance(j, Expr) and not j.is_number):
                if ((j < 0) is True) or ((j >= self.shape[1]) is True) or\
                   ((i < 0) is True) or ((i >= self.shape[0]) is True):
                    raise ValueError("index out of boundary")
                from sympy.matrices.expressions.matexpr import MatrixElement
                return MatrixElement(self, i, j)

            if isinstance(i, slice):
                i = range(self.rows)[i]
            elif is_sequence(i):
                pass
            else:
                i = [i]
            if isinstance(j, slice):
                j = range(self.cols)[j]
            elif is_sequence(j):
                pass
            else:
                j = [j]
            return self.extract(i, j)

    else:
        # Index/slice like a flattened list
        rows, cols = self.shape

        # Raise the appropriate exception:
        if not rows * cols:
            return [][key]

        rep = self._rep.rep
        domain = rep.domain
        is_slice = isinstance(key, slice)

        if is_slice:
            values = [rep.getitem(*divmod(n, cols)) for n in range(rows * cols)[key]]
        else:
            values = [rep.getitem(*divmod(index_(key), cols))]

        if domain != EXRAW:
            to_sympy = domain.to_sympy
            values = [to_sympy(val) for val in values]

        if is_slice:
            return values
        else:
            return values[0]
