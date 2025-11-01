from __future__ import annotations

from typing import TYPE_CHECKING, overload

from collections import defaultdict
from collections.abc import Iterable, Sequence
from inspect import isfunction
from functools import reduce

from sympy.assumptions.refine import refine
from sympy.core import SympifyError, Add
from sympy.core.basic import Atom, Basic
from sympy.core.kind import UndefinedKind
from sympy.core.numbers import Integer, Float
from sympy.core.mod import Mod
from sympy.core.symbol import Symbol, Dummy
from sympy.core.sympify import sympify, _sympify
from sympy.core.function import diff
from sympy.polys import cancel
from sympy.functions.elementary.complexes import Abs, re, im
from sympy.printing import sstr
from sympy.functions.elementary.miscellaneous import Max, Min, sqrt
from sympy.functions.special.tensor_functions import KroneckerDelta, LeviCivita
from sympy.core.singleton import S
from sympy.printing.defaults import Printable
from sympy.printing.str import StrPrinter
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.combinatorial.factorials import binomial, factorial

import mpmath as mp
from collections.abc import Callable
from sympy.utilities.iterables import reshape
from sympy.core.expr import Expr
from sympy.core.power import Pow
from sympy.core.symbol import uniquely_named_symbol

from .utilities import _dotprodsimp, _simplify as _utilities_simplify
from sympy.polys.polytools import Poly
from sympy.utilities.iterables import flatten, is_sequence
from sympy.utilities.misc import as_int, filldedent
from sympy.core.decorators import call_highest_priority
from sympy.core.logic import fuzzy_and
from sympy.tensor.array import NDimArray
from sympy.utilities.iterables import NotIterable

from .utilities import _get_intermediate_simp_bool

from .kind import MatrixKind

from .exceptions import (
    MatrixError, ShapeError, NonSquareMatrixError, NonInvertibleMatrixError,
)

from .utilities import _iszero, _is_zero_after_expand_mul

from .determinant import (
    _find_reasonable_pivot, _find_reasonable_pivot_naive,
    _adjugate, _charpoly, _cofactor, _cofactor_matrix, _per,
    _det, _det_bareiss, _det_berkowitz, _det_bird, _det_laplace, _det_LU,
    _minor, _minor_submatrix)

from .reductions import _is_echelon, _echelon_form, _rank, _rref

from .solvers import (
    _diagonal_solve, _lower_triangular_solve, _upper_triangular_solve,
    _cholesky_solve, _LDLsolve, _LUsolve, _QRsolve, _gauss_jordan_solve,
    _pinv_solve, _cramer_solve, _solve, _solve_least_squares)

from .inverse import (
    _pinv, _inv_ADJ, _inv_GE, _inv_LU, _inv_CH, _inv_LDL, _inv_QR,
    _inv, _inv_block)

from .subspaces import _columnspace, _nullspace, _rowspace, _orthogonalize

from .eigen import (
    _eigenvals, _eigenvects,
    _bidiagonalize, _bidiagonal_decomposition,
    _is_diagonalizable, _diagonalize,
    _is_positive_definite, _is_positive_semidefinite,
    _is_negative_definite, _is_negative_semidefinite, _is_indefinite,
    _jordan_form, _left_eigenvects, _singular_values)

from .decompositions import (
    _rank_decomposition, _cholesky, _LDLdecomposition,
    _LUdecomposition, _LUdecomposition_Simple, _LUdecompositionFF,
    _singular_value_decomposition, _QRdecomposition, _upper_hessenberg_decomposition)

from .graph import (
    _connected_components, _connected_components_decomposition,
    _strongly_connected_components, _strongly_connected_components_decomposition)


if TYPE_CHECKING:
    from abc import ABCMeta, abstractmethod
else:
    from abc import abstractmethod
    ABCMeta = type


if TYPE_CHECKING:

    from typing import Literal, TypeVar, Iterator, Mapping, Any
    from typing_extensions import Self
    from sympy.combinatorics import Permutation
    from sympy.logic.boolalg import Boolean
    from sympy.integrals.integrals import SymbolLimits
    from sympy.matrices.expressions.matexpr import MatrixExpr
    from sympy.matrices.expressions.permutation import PermutationMatrix
    from sympy.matrices.expressions.blockmatrix import BlockDiagMatrix, BlockMatrix

    Tmat = TypeVar('Tmat', bound='MatrixBase')
    Tbasic = TypeVar('Tbasic', bound=Basic)
    SExpr = Expr | complex
    SBasic = Basic | complex
    Slice = slice | list[int]


__doctest_requires__ = {
    ('MatrixBase.is_indefinite',
     'MatrixBase.is_positive_definite',
     'MatrixBase.is_positive_semidefinite',
     'MatrixBase.is_negative_definite',
     'MatrixBase.is_negative_semidefinite'): ['matplotlib'],
}


class MatrixBase(Printable):
    """All common matrix operations including basic arithmetic, shaping,
    and special matrices like `zeros`, and `eye`."""

    _op_priority = 10.01

    # Added just for numpy compatibility
    __array_priority__ = 11

    is_Matrix = True
    _class_priority = 3
    zero = S.Zero
    one = S.One

    _diff_wrt: bool = True
    _simplify = None

    if TYPE_CHECKING:

        @property
        def rows(self) -> int:
            ...

        @property
        def cols(self) -> int:
            ...

    @overload
    @classmethod
    def _sympify(cls, expr: SExpr, /) -> Expr: ...
    @overload
    @classmethod
    def _sympify(cls, expr: Poly, /) -> Poly: ...

    @classmethod
    def _sympify(cls, expr: SExpr | Poly, /) -> Expr | Poly:
        if isinstance(expr, Poly):
            return expr
        return sympify(expr)

    @overload
    @classmethod
    def _new(cls, rows: int, cols: int, mat: Sequence[SExpr], /, copy: bool = False) -> Self: ...
    @overload
    @classmethod
    def _new(cls, rows: int, cols: int, func: Callable[[int, int], SExpr], /) -> Self: ...
    @overload
    @classmethod
    def _new(cls, mat: Sequence[Sequence[SExpr]] | Self, /) -> Self: ...
    @overload
    @classmethod
    def _new(cls, /) -> Self: ...
    @overload
    @classmethod
    def _new(cls, elements: Sequence[SExpr], /) -> Self: ...

    @classmethod
    @abstractmethod
    def _new(cls, *args, **kwargs) -> Self:
        """`_new` must, at minimum, be callable as
        `_new(rows, cols, mat) where mat is a flat list of the
        elements of the matrix."""
        raise NotImplementedError("Subclasses must implement this.")

    @classmethod
    def _as_type(cls, mat: MatrixBase) -> Self:
        if not isinstance(mat, cls):
            mat = cls._new(mat) # type: ignore
        return mat # type: ignore

    def as_explicit(self) -> Self:
        """Convert MatrixExpr to a matrix with explicit elements."""
        return self

    def __eq__(self, other):
        raise NotImplementedError("Subclasses must implement this.")

    @overload
    def __getitem__(self, key: tuple[int, int], /) -> Expr: ...
    @overload
    def __getitem__(self, key: tuple[int, Slice], /) -> Self: ...
    @overload
    def __getitem__(self, key: tuple[Slice, int], /) -> Self: ...
    @overload
    def __getitem__(self, key: tuple[Slice, Slice], /) -> Self: ...
    @overload
    def __getitem__(self, key: int, /) -> Expr: ...
    @overload
    def __getitem__(self, key: slice, /) -> list[Expr]: ...

    @abstractmethod
    def __getitem__(self, key: tuple[int | Slice, int | Slice] | int | slice, /
                    ) -> Expr | Self | list[Expr]:
        """Implementations of __getitem__ should accept ints, in which
        case the matrix is indexed as a flat list, tuples (i,j) in which
        case the (i,j) entry is returned, slices, or mixed tuples (a,b)
        where a and b are any combination of slices and integers."""
        raise NotImplementedError("Subclasses must implement this.")

    @property
    def shape(self) -> tuple[int, int]:
        """The shape (dimensions) of the matrix as the 2-tuple (rows, cols).

        Examples
        ========

        >>> from sympy import zeros
        >>> M = zeros(2, 3)
        >>> M.shape
        (2, 3)
        >>> M.rows
        2
        >>> M.cols
        3
        """
        return (self.rows, self.cols)

    def _eval_col_del(self, col: int, /) -> Self:
        def entry(i, j):
            return self[i, j] if j < col else self[i, j + 1]
        return self._new(self.rows, self.cols - 1, entry)

    def _eval_col_insert(self, pos: int, other: Self, /) -> Self:

        def entry(i, j):
            if j < pos:
                return self[i, j]
            elif pos <= j < pos + other.cols:
                return other[i, j - pos]
            return self[i, j - other.cols]

        return self._new(self.rows, self.cols + other.cols, entry)

    def _eval_col_join(self, other: Self, /) -> Self:
        rows = self.rows

        def entry(i, j):
            if i < rows:
                return self[i, j]
            return other[i - rows, j]

        return classof(self, other)._new(self.rows + other.rows, self.cols,
                                         entry)

    def _eval_extract(self, rowsList: Sequence[int], colsList: Sequence[int], /) -> Self:
        mat = self.flat()
        cols = self.cols
        indices = (i * cols + j for i in rowsList for j in colsList)
        return self._new(len(rowsList), len(colsList),
                         [mat[i] for i in indices])

    def _eval_get_diag_blocks(self) -> list[Self]:
        sub_blocks: list[Self] = []

        def recurse_sub_blocks(M: Self) -> None:
            for i in range(1, M.shape[0] + 1):
                if i == 1:
                    to_the_right = M[0, i:]
                    to_the_bottom = M[i:, 0]
                else:
                    to_the_right = M[:i, i:]
                    to_the_bottom = M[i:, :i]
                if any(to_the_right) or any(to_the_bottom): # type: ignore
                    continue
                sub_blocks.append(M[:i, :i])
                if M.shape != M[:i, :i].shape:
                    recurse_sub_blocks(M[i:, i:])
                return

        recurse_sub_blocks(self)
        return sub_blocks

    def _eval_row_del(self, row: int, /) -> Self:
        def entry(i, j):
            return self[i, j] if i < row else self[i + 1, j]
        return self._new(self.rows - 1, self.cols, entry)

    def _eval_row_insert(self, pos: int, other: Self, /) -> Self:
        entries = self.flat()
        insert_pos = pos * self.cols
        entries[insert_pos:insert_pos] = other.flat()
        return self._new(self.rows + other.rows, self.cols, entries)

    def _eval_row_join(self, other: Self, /) -> Self:
        cols = self.cols

        def entry(i, j):
            if j < cols:
                return self[i, j]
            return other[i, j - cols]

        return classof(self, other)._new(self.rows, self.cols + other.cols,
                                         entry)

    def _eval_tolist(self) -> list[list[Expr]]:
        return [self[i,:].flat() for i in range(self.rows)]

    def _eval_todok(self) -> dict[tuple[int, int], Expr]:
        dok = {}
        rows, cols = self.shape
        for i in range(rows):
            for j in range(cols):
                val = self[i, j]
                if val != self.zero:
                    dok[i, j] = val
        return dok

    @classmethod
    def _eval_from_dok(cls, rows: int, cols: int, dok: dict[tuple[int, int], Expr], /) -> Self:
        out_flat: list[Expr] = [cls.zero] * (rows * cols)
        for (i, j), val in dok.items():
            out_flat[i * cols + j] = val
        return cls._new(rows, cols, out_flat)

    def _eval_vec(self) -> Self:
        rows = self.rows

        def entry(n, _):
            # we want to read off the columns first
            j = n // rows
            i = n - j * rows
            return self[i, j]

        return self._new(len(self), 1, entry)

    def _eval_vech(self, diagonal: bool) -> Self:
        c = self.cols
        v = []
        if diagonal:
            for j in range(c):
                for i in range(j, c):
                    v.append(self[i, j])
        else:
            for j in range(c):
                for i in range(j + 1, c):
                    v.append(self[i, j])
        return self._new(len(v), 1, v)

    def col_del(self, col: int, /) -> Self:
        """Delete the specified column."""
        if col < 0:
            col += self.cols
        if not 0 <= col < self.cols:
            raise IndexError("Column {} is out of range.".format(col))
        return self._eval_col_del(col)

    def col_insert(self, pos: int, other: Self, /) -> Self:
        """Insert one or more columns at the given column position.

        Examples
        ========

        >>> from sympy import zeros, ones
        >>> M = zeros(3)
        >>> V = ones(3, 1)
        >>> M.col_insert(1, V)
        Matrix([
        [0, 1, 0, 0],
        [0, 1, 0, 0],
        [0, 1, 0, 0]])

        See Also
        ========

        col
        row_insert
        """
        # Allows you to build a matrix even if it is null matrix
        if not self:
            return self._new(other)

        pos = as_int(pos)

        if pos < 0:
            pos = self.cols + pos
        if pos < 0:
            pos = 0
        elif pos > self.cols:
            pos = self.cols

        if self.rows != other.rows:
            raise ShapeError(
                "The matrices have incompatible number of rows ({} and {})"
                .format(self.rows, other.rows))

        return self._eval_col_insert(pos, other)

    def col_join(self, other: Self, /) -> Self:
        """Concatenates two matrices along self's last and other's first row.

        Examples
        ========

        >>> from sympy import zeros, ones
        >>> M = zeros(3)
        >>> V = ones(1, 3)
        >>> M.col_join(V)
        Matrix([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [1, 1, 1]])

        See Also
        ========

        col
        row_join
        """
        # A null matrix can always be stacked (see  #10770)
        if self.rows == 0 and self.cols != other.cols:
            return self._new(0, other.cols, []).col_join(other)

        if self.cols != other.cols:
            raise ShapeError(
                "The matrices have incompatible number of columns ({} and {})"
                .format(self.cols, other.cols))
        return self._eval_col_join(other)

    def col(self, j: int, /) -> Self:
        """Elementary column selector.

        Examples
        ========

        >>> from sympy import eye
        >>> eye(2).col(0)
        Matrix([
        [1],
        [0]])

        See Also
        ========

        row
        col_del
        col_join
        col_insert
        """
        return self[:, j]

    def extract(self, rowsList: Sequence[int], colsList: Sequence[int], /) -> Self:
        r"""Return a submatrix by specifying a list of rows and columns.
        Negative indices can be given. All indices must be in the range
        $-n \le i < n$ where $n$ is the number of rows or columns.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(4, 3, range(12))
        >>> m
        Matrix([
        [0,  1,  2],
        [3,  4,  5],
        [6,  7,  8],
        [9, 10, 11]])
        >>> m.extract([0, 1, 3], [0, 1])
        Matrix([
        [0,  1],
        [3,  4],
        [9, 10]])

        Rows or columns can be repeated:

        >>> m.extract([0, 0, 1], [-1])
        Matrix([
        [2],
        [2],
        [5]])

        Every other row can be taken by using range to provide the indices:

        >>> m.extract(range(0, m.rows, 2), [-1])
        Matrix([
        [2],
        [8]])

        RowsList or colsList can also be a list of booleans, in which case
        the rows or columns corresponding to the True values will be selected:

        >>> m.extract([0, 1, 2, 3], [True, False, True])
        Matrix([
        [0,  2],
        [3,  5],
        [6,  8],
        [9, 11]])
        """

        if not is_sequence(rowsList) or not is_sequence(colsList):
            raise TypeError("rowsList and colsList must be iterable")
        # ensure rowsList and colsList are lists of integers
        if rowsList and all(isinstance(i, bool) for i in rowsList):
            rowsList = [index for index, item in enumerate(rowsList) if item]
        if colsList and all(isinstance(i, bool) for i in colsList):
            colsList = [index for index, item in enumerate(colsList) if item]

        # ensure everything is in range
        rowsList = [a2idx(k, self.rows) for k in rowsList]
        colsList = [a2idx(k, self.cols) for k in colsList]

        return self._eval_extract(rowsList, colsList)

    def get_diag_blocks(self) -> list[Self]:
        """Obtains the square sub-matrices on the main diagonal of a square matrix.

        Useful for inverting symbolic matrices or solving systems of
        linear equations which may be decoupled by having a block diagonal
        structure.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.abc import x, y, z
        >>> A = Matrix([[1, 3, 0, 0], [y, z*z, 0, 0], [0, 0, x, 0], [0, 0, 0, 0]])
        >>> a1, a2, a3 = A.get_diag_blocks()
        >>> a1
        Matrix([
        [1,    3],
        [y, z**2]])
        >>> a2
        Matrix([[x]])
        >>> a3
        Matrix([[0]])

        """
        return self._eval_get_diag_blocks()

    @classmethod
    def hstack(cls, *args: Self) -> Self:
        """Return a matrix formed by joining args horizontally (i.e.
        by repeated application of row_join).

        Examples
        ========

        >>> from sympy import Matrix, eye
        >>> Matrix.hstack(eye(2), 2*eye(2))
        Matrix([
        [1, 0, 2, 0],
        [0, 1, 0, 2]])
        """
        if len(args) == 0:
            return cls._new()

        kls = type(args[0])
        return reduce(kls.row_join, args)

    def reshape(self, rows: int, cols: int) -> Self:
        """Reshape the matrix. Total number of elements must remain the same.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(2, 3, lambda i, j: 1)
        >>> m
        Matrix([
        [1, 1, 1],
        [1, 1, 1]])
        >>> m.reshape(1, 6)
        Matrix([[1, 1, 1, 1, 1, 1]])
        >>> m.reshape(3, 2)
        Matrix([
        [1, 1],
        [1, 1],
        [1, 1]])

        """
        if self.rows * self.cols != rows * cols:
            raise ValueError("Invalid reshape parameters %d %d" % (rows, cols))
        dok = {divmod(i*self.cols + j, cols):
            v for (i, j), v in self.todok().items()}
        return self._eval_from_dok(rows, cols, dok)

    def row_del(self, row: int, /) -> Self:
        """Delete the specified row."""
        if row < 0:
            row += self.rows
        if not 0 <= row < self.rows:
            raise IndexError("Row {} is out of range.".format(row))

        return self._eval_row_del(row)

    def row_insert(self, pos: int, other: Self, /) -> Self:
        """Insert one or more rows at the given row position.

        Examples
        ========

        >>> from sympy import zeros, ones
        >>> M = zeros(3)
        >>> V = ones(1, 3)
        >>> M.row_insert(1, V)
        Matrix([
        [0, 0, 0],
        [1, 1, 1],
        [0, 0, 0],
        [0, 0, 0]])

        See Also
        ========

        row
        col_insert
        """
        # Allows you to build a matrix even if it is null matrix
        if not self:
            return self._new(other)

        pos = as_int(pos)

        if pos < 0:
            pos = self.rows + pos
        if pos < 0:
            pos = 0
        elif pos > self.rows:
            pos = self.rows

        if self.cols != other.cols:
            raise ShapeError(
                "The matrices have incompatible number of columns ({} and {})"
                .format(self.cols, other.cols))

        return self._eval_row_insert(pos, other)

    def row_join(self, other: Self) -> Self:
        """Concatenates two matrices along self's last and rhs's first column

        Examples
        ========

        >>> from sympy import zeros, ones
        >>> M = zeros(3)
        >>> V = ones(3, 1)
        >>> M.row_join(V)
        Matrix([
        [0, 0, 0, 1],
        [0, 0, 0, 1],
        [0, 0, 0, 1]])

        See Also
        ========

        row
        col_join
        """
        # A null matrix can always be stacked (see  #10770)
        if self.cols == 0 and self.rows != other.rows:
            return self._new(other.rows, 0, []).row_join(other)

        if self.rows != other.rows:
            raise ShapeError(
                "The matrices have incompatible number of rows ({} and {})"
                .format(self.rows, other.rows))
        return self._eval_row_join(other)

    def diagonal(self, k: int = 0) -> Self:
        """Returns the kth diagonal of self. The main diagonal
        corresponds to `k=0`; diagonals above and below correspond to
        `k > 0` and `k < 0`, respectively. The values of `self[i, j]`
        for which `j - i = k`, are returned in order of increasing
        `i + j`, starting with `i + j = |k|`.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(3, 3, lambda i, j: j - i); m
        Matrix([
        [ 0,  1, 2],
        [-1,  0, 1],
        [-2, -1, 0]])
        >>> _.diagonal()
        Matrix([[0, 0, 0]])
        >>> m.diagonal(1)
        Matrix([[1, 1]])
        >>> m.diagonal(-2)
        Matrix([[-2]])

        Even though the diagonal is returned as a Matrix, the element
        retrieval can be done with a single index:

        >>> Matrix.diag(1, 2, 3).diagonal()[1]  # instead of [0, 1]
        2

        See Also
        ========

        diag
        """
        k = as_int(k)
        rv = [self[r, r+k] for r in range(max(0, -k), min(self.rows, self.cols - k))]
        return self._new(1, len(rv), rv)

    def row(self, i: int, /) -> Self:
        """Elementary row selector.

        Examples
        ========

        >>> from sympy import eye
        >>> eye(2).row(0)
        Matrix([[1, 0]])

        See Also
        ========

        col
        row_del
        row_join
        row_insert
        """
        return self[i, :]

    def todok(self) -> dict[tuple[int, int], Expr]:
        """Return the matrix as dictionary of keys.

        Examples
        ========

        >>> from sympy import Matrix
        >>> M = Matrix.eye(3)
        >>> M.todok()
        {(0, 0): 1, (1, 1): 1, (2, 2): 1}
        """
        return self._eval_todok()

    @classmethod
    def from_dok(cls, rows: int, cols: int, dok: dict[tuple[int, int], Expr], /) -> Self:
        """Create a matrix from a dictionary of keys.

        Examples
        ========

        >>> from sympy import Matrix
        >>> d = {(0, 0): 1, (1, 2): 3, (2, 1): 4}
        >>> Matrix.from_dok(3, 3, d)
        Matrix([
        [1, 0, 0],
        [0, 0, 3],
        [0, 4, 0]])
        """
        dok = {ij: cls._sympify(val) for ij, val in dok.items()}
        return cls._eval_from_dok(rows, cols, dok)

    def tolist(self) -> list[list[Expr]]:
        """Return the Matrix as a nested Python list.

        Examples
        ========

        >>> from sympy import Matrix, ones
        >>> m = Matrix(3, 3, range(9))
        >>> m
        Matrix([
        [0, 1, 2],
        [3, 4, 5],
        [6, 7, 8]])
        >>> m.tolist()
        [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
        >>> ones(3, 0).tolist()
        [[], [], []]

        When there are no rows then it will not be possible to tell how
        many columns were in the original matrix:

        >>> ones(0, 3).tolist()
        []

        """
        if not self.rows:
            return []
        if not self.cols:
            return [[] for _ in range(self.rows)]
        return self._eval_tolist()

    def todod(M) -> dict[int, dict[int, Expr]]:
        """Returns matrix as dict of dicts containing non-zero elements of the Matrix

        Examples
        ========

        >>> from sympy import Matrix
        >>> A = Matrix([[0, 1],[0, 3]])
        >>> A
        Matrix([
        [0, 1],
        [0, 3]])
        >>> A.todod()
        {0: {1: 1}, 1: {1: 3}}


        """
        rowsdict = {}
        Mlol = M.tolist()
        for i, Mi in enumerate(Mlol):
            row = {j: Mij for j, Mij in enumerate(Mi) if Mij}
            if row:
                rowsdict[i] = row
        return rowsdict

    def vec(self) -> Self:
        """Return the Matrix converted into a one column matrix by stacking columns

        Examples
        ========

        >>> from sympy import Matrix
        >>> m=Matrix([[1, 3], [2, 4]])
        >>> m
        Matrix([
        [1, 3],
        [2, 4]])
        >>> m.vec()
        Matrix([
        [1],
        [2],
        [3],
        [4]])

        See Also
        ========

        vech
        """
        return self._eval_vec()

    def vech(self, diagonal: bool = True, check_symmetry: bool = True) -> Self:
        """Reshapes the matrix into a column vector by stacking the
        elements in the lower triangle.

        Parameters
        ==========

        diagonal : bool, optional
            If ``True``, it includes the diagonal elements.

        check_symmetry : bool, optional
            If ``True``, it checks whether the matrix is symmetric.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m=Matrix([[1, 2], [2, 3]])
        >>> m
        Matrix([
        [1, 2],
        [2, 3]])
        >>> m.vech()
        Matrix([
        [1],
        [2],
        [3]])
        >>> m.vech(diagonal=False)
        Matrix([[2]])

        Notes
        =====

        This should work for symmetric matrices and ``vech`` can
        represent symmetric matrices in vector form with less size than
        ``vec``.

        See Also
        ========

        vec
        """
        if not self.is_square:
            raise NonSquareMatrixError

        if check_symmetry and not self.is_symmetric():
            raise ValueError("The matrix is not symmetric.")

        return self._eval_vech(diagonal)

    @classmethod
    def vstack(cls, *args: Self) -> Self:
        """Return a matrix formed by joining args vertically (i.e.
        by repeated application of col_join).

        Examples
        ========

        >>> from sympy import Matrix, eye
        >>> Matrix.vstack(eye(2), 2*eye(2))
        Matrix([
        [1, 0],
        [0, 1],
        [2, 0],
        [0, 2]])
        """
        if len(args) == 0:
            return cls._new()

        kls = type(args[0])
        return reduce(kls.col_join, args)

    @classmethod
    def _eval_diag(cls, rows: int, cols: int, diag_dict: dict[tuple[int, int], SExpr], /) -> Self:
        """diag_dict is a defaultdict containing
        all the entries of the diagonal matrix."""
        def entry(i, j):
            return diag_dict[(i, j)]
        return cls._new(rows, cols, entry)

    @classmethod
    def _eval_eye(cls, rows: int, cols: int) -> Self:
        vals: list[Expr] = [cls.zero]*(rows*cols)
        vals[::cols+1] = [cls.one]*min(rows, cols)
        return cls._new(rows, cols, vals, copy=False)

    @classmethod
    def _eval_jordan_block(cls, size: int, eigenvalue: SExpr,
                                band: Literal['upper', 'lower'] = 'upper'
                           ) -> Self:
        if band == 'lower':
            def entry(i, j):
                if i == j:
                    return eigenvalue
                elif j + 1 == i:
                    return cls.one
                return cls.zero
        else:
            def entry(i, j):
                if i == j:
                    return eigenvalue
                elif i + 1 == j:
                    return cls.one
                return cls.zero
        return cls._new(size, size, entry)

    @classmethod
    def _eval_ones(cls, rows, cols) -> Self:
        def entry(i, j):
            return cls.one
        return cls._new(rows, cols, entry)

    @classmethod
    def _eval_zeros(cls, rows, cols) -> Self:
        return cls._new(rows, cols, [cls.zero]*(rows*cols), copy=False)

    @classmethod
    def _eval_wilkinson(cls, n) -> tuple[Self, Self]:
        def entry(i, j):
            return cls.one if i + 1 == j else cls.zero

        D = cls._new(2*n + 1, 2*n + 1, entry)

        wminus = cls.diag(list(range(-n, n + 1)), unpack=True) + D + D.T
        wplus = abs(cls.diag(list(range(-n, n + 1)), unpack=True)) + D + D.T

        return wminus, cls._as_type(wplus)

    @overload
    @classmethod
    def diag(kls, *args: SExpr | Sequence[SExpr] | MatrixBase,
                strict: bool = False,
                unpack: bool = True,
                rows: int | None = None,
                cols: int | None = None,
             cls: None = None) -> Self:
        ...

    @overload
    @classmethod
    def diag(kls, *args: SExpr | Sequence[SExpr] | MatrixBase,
                strict: bool = False,
                unpack: bool = True,
                rows: int | None = None,
                cols: int | None = None,
             cls: type[Tmat]) -> Tmat:
        ...

    @classmethod
    def diag(kls, *args: SExpr | Sequence[SExpr] | MatrixBase,
                strict: bool = False,
                unpack: bool = True,
                rows: int | None = None,
                cols: int | None = None,
             cls: type[Tmat] | None = None) -> Self | Tmat:
        """Returns a matrix with the specified diagonal.
        If matrices are passed, a block-diagonal matrix
        is created (i.e. the "direct sum" of the matrices).

        kwargs
        ======

        rows : rows of the resulting matrix; computed if
               not given.

        cols : columns of the resulting matrix; computed if
               not given.

        cls : class for the resulting matrix

        unpack : bool which, when True (default), unpacks a single
        sequence rather than interpreting it as a Matrix.

        strict : bool which, when False (default), allows Matrices to
        have variable-length rows.

        Examples
        ========

        >>> from sympy import Matrix
        >>> Matrix.diag(1, 2, 3)
        Matrix([
        [1, 0, 0],
        [0, 2, 0],
        [0, 0, 3]])

        The current default is to unpack a single sequence. If this is
        not desired, set `unpack=False` and it will be interpreted as
        a matrix.

        >>> Matrix.diag([1, 2, 3]) == Matrix.diag(1, 2, 3)
        True

        When more than one element is passed, each is interpreted as
        something to put on the diagonal. Lists are converted to
        matrices. Filling of the diagonal always continues from
        the bottom right hand corner of the previous item: this
        will create a block-diagonal matrix whether the matrices
        are square or not.

        >>> col = [1, 2, 3]
        >>> row = [[4, 5]]
        >>> Matrix.diag(col, row)
        Matrix([
        [1, 0, 0],
        [2, 0, 0],
        [3, 0, 0],
        [0, 4, 5]])

        When `unpack` is False, elements within a list need not all be
        of the same length. Setting `strict` to True would raise a
        ValueError for the following:

        >>> Matrix.diag([[1, 2, 3], [4, 5], [6]], unpack=False)
        Matrix([
        [1, 2, 3],
        [4, 5, 0],
        [6, 0, 0]])

        The type of the returned matrix can be set with the ``cls``
        keyword.

        >>> from sympy import ImmutableMatrix
        >>> from sympy.utilities.misc import func_name
        >>> func_name(Matrix.diag(1, cls=ImmutableMatrix))
        'ImmutableDenseMatrix'

        A zero dimension matrix can be used to position the start of
        the filling at the start of an arbitrary row or column:

        >>> from sympy import ones
        >>> r2 = ones(0, 2)
        >>> Matrix.diag(r2, 1, 2)
        Matrix([
        [0, 0, 1, 0],
        [0, 0, 0, 2]])

        See Also
        ========
        eye
        diagonal
        .dense.diag
        .expressions.blockmatrix.BlockMatrix
        .sparsetools.banded
       """
        from sympy.matrices.matrixbase import MatrixBase
        from sympy.matrices.dense import Matrix
        from sympy.matrices import SparseMatrix

        klass = cls if cls is not None else kls

        args2: tuple[SExpr | list[SExpr] | MatrixBase, ...]

        if unpack and len(args) == 1 and is_sequence(args[0]) and \
                not isinstance(args[0], MatrixBase):
            args2 = args[0] # type: ignore
        else:
            args2 = args # type: ignore

        # fill a default dict with the diagonal entries
        diag_entries: dict[tuple[int, int], SExpr] = defaultdict(int)
        rmax = cmax = 0  # keep track of the biggest index seen
        for m in args2:
            if isinstance(m, list):
                if strict:
                    # if malformed, Matrix will raise an error
                    _ = Matrix(m)
                    r, c = _.shape
                    m2 = _.tolist()
                else:
                    r, c, smat = SparseMatrix._handle_creation_inputs(m)
                    for (i, j), _ in smat.items():
                        diag_entries[(i + rmax, j + cmax)] = _
                    m2 = []  # to skip process below
            elif isinstance(m, MatrixBase):
                # convert to list of lists
                r, c = m.shape
                m2 = m.tolist()
            else:  # in this case, we're a single value
                diag_entries[(rmax, cmax)] = m
                rmax += 1
                cmax += 1
                continue
            # process list of lists
            for i, mi in enumerate(m2):
                for j, _ in enumerate(mi):
                    diag_entries[(i + rmax, j + cmax)] = _
            rmax += r
            cmax += c
        if rows is None:
            rows, cols = cols, rows
        if rows is None:
            rows, cols = rmax, cmax
        else:
            cols = rows if cols is None else cols
        if rows < rmax or cols < cmax:
            raise ValueError(filldedent('''
                The constructed matrix is {} x {} but a size of {} x {}
                was specified.'''.format(rmax, cmax, rows, cols)))
        return klass._eval_diag(rows, cols, diag_entries)

    @overload
    @classmethod
    def eye(kls, rows: int, cols: int | None = None, *,
                cls: None = None) -> Self:
        ...

    @overload
    @classmethod
    def eye(kls, rows: int, cols: int | None = None, *,
                cls: type[Tmat]) -> Tmat:
        ...

    @classmethod
    def eye(kls, rows: int, cols: int | None = None, *,
                 cls: type[Tmat] | None = None) -> Self | Tmat:
        """Returns an identity matrix.

        Parameters
        ==========

        rows : rows of the matrix
        cols : cols of the matrix (if None, cols=rows)

        kwargs
        ======
        cls : class of the returned matrix
        """
        if cols is None:
            cols2 = rows
        else:
            cols2 = cols
        if rows < 0 or cols2 < 0:
            raise ValueError("Cannot create a {} x {} matrix. "
                             "Both dimensions must be positive".format(rows, cols2))
        klass = cls if cls is not None else kls
        rows, cols2 = as_int(rows), as_int(cols2)

        return klass._eval_eye(rows, cols2)

    @overload
    @classmethod
    def jordan_block(kls, size: int | None = None,
                         eigenvalue: SExpr | None = None, *,
                         band: Literal['upper', 'lower'] = 'upper',
                         cls: None = None,
                         eigenval: SExpr | None = None) -> Self:
        ...

    @overload
    @classmethod
    def jordan_block(kls, size: int | None = None,
                         eigenvalue: SExpr | None = None, *,
                         band: Literal['upper', 'lower'] = 'upper',
                         cls: type[Tmat],
                         eigenval: SExpr | None = None) -> Tmat:
        ...

    @classmethod
    def jordan_block(kls, size: int | None = None,
                         eigenvalue: SExpr | None = None, *,
                         band: Literal['upper', 'lower'] = 'upper',
                         cls: type[Tmat] | None = None,
                         eigenval: SExpr | None = None) -> Self | Tmat:
        """Returns a Jordan block

        Parameters
        ==========

        size : Integer, optional
            Specifies the shape of the Jordan block matrix.

        eigenvalue : Number or Symbol
            Specifies the value for the main diagonal of the matrix.

            .. note::
                The keyword ``eigenval`` is also specified as an alias
                of this keyword, but it is not recommended to use.

                We may deprecate the alias in later release.

        band : 'upper' or 'lower', optional
            Specifies the position of the off-diagonal to put `1` s on.

        cls : Matrix, optional
            Specifies the matrix class of the output form.

            If it is not specified, the class type where the method is
            being executed on will be returned.

        Returns
        =======

        Matrix
            A Jordan block matrix.

        Raises
        ======

        ValueError
            If insufficient arguments are given for matrix size
            specification, or no eigenvalue is given.

        Examples
        ========

        Creating a default Jordan block:

        >>> from sympy import Matrix
        >>> from sympy.abc import x
        >>> Matrix.jordan_block(4, x)
        Matrix([
        [x, 1, 0, 0],
        [0, x, 1, 0],
        [0, 0, x, 1],
        [0, 0, 0, x]])

        Creating an alternative Jordan block matrix where `1` is on
        lower off-diagonal:

        >>> Matrix.jordan_block(4, x, band='lower')
        Matrix([
        [x, 0, 0, 0],
        [1, x, 0, 0],
        [0, 1, x, 0],
        [0, 0, 1, x]])

        Creating a Jordan block with keyword arguments

        >>> Matrix.jordan_block(size=4, eigenvalue=x)
        Matrix([
        [x, 1, 0, 0],
        [0, x, 1, 0],
        [0, 0, x, 1],
        [0, 0, 0, x]])

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Jordan_matrix
        """
        klass = cls if cls is not None else kls

        if eigenvalue is not None:
            if eigenval is not None and eigenvalue != eigenval:
                raise ValueError("Cannot supply both 'eigenvalue' and 'eigenval'")
            e = eigenvalue
        elif eigenval is not None:
            e = eigenval
        else:
            raise ValueError("Must supply an eigenvalue")

        if size is None:
            raise ValueError("Must supply a matrix size")
        else:
            size2 = as_int(size)

        return klass._eval_jordan_block(size2, e, band)

    @overload
    @classmethod
    def ones(kls, rows: int, cols: int | None = None, *,
                cls: None = None) -> Self:
        ...

    @overload
    @classmethod
    def ones(kls, rows: int, cols: int | None = None, *,
                cls: type[Tmat]) -> Tmat:
        ...

    @classmethod
    def ones(kls, rows: int, cols: int | None = None, *,
                cls: type[Tmat] | None = None) -> Self | Tmat:
        """Returns a matrix of ones.

        Parameters
        ==========

        rows : rows of the matrix
        cols : cols of the matrix (if None, cols=rows)

        kwargs
        ======
        cls : class of the returned matrix
        """
        if cols is None:
            cols = rows
        klass = cls if cls is not None else kls
        rows, cols = as_int(rows), as_int(cols)

        return klass._eval_ones(rows, cols)

    @overload
    @classmethod
    def zeros(kls, rows: int, cols: int | None = None, *,
                 cls: None = None) -> Self:
        ...

    @overload
    @classmethod
    def zeros(kls, rows: int, cols: int | None = None, *,
                 cls: type[Tmat]) -> Tmat:
        ...

    @classmethod
    def zeros(kls, rows: int, cols: int | None = None, *,
                 cls: type[Tmat] | None = None) -> Self | Tmat:
        """Returns a matrix of zeros.

        Parameters
        ==========

        rows : rows of the matrix
        cols : cols of the matrix (if None, cols=rows)

        kwargs
        ======
        cls : class of the returned matrix
        """
        if cols is None:
            cols = rows
        if rows < 0 or cols < 0:
            raise ValueError("Cannot create a {} x {} matrix. "
                             "Both dimensions must be positive".format(rows, cols))
        klass = cls if cls is not None else kls
        rows, cols = as_int(rows), as_int(cols)

        return klass._eval_zeros(rows, cols)

    @classmethod
    def companion(kls, poly: Poly) -> Self:
        """Returns a companion matrix of a polynomial.

        Examples
        ========

        >>> from sympy import Matrix, Poly, Symbol, symbols
        >>> x = Symbol('x')
        >>> c0, c1, c2, c3, c4 = symbols('c0:5')
        >>> p = Poly(c0 + c1*x + c2*x**2 + c3*x**3 + c4*x**4 + x**5, x)
        >>> Matrix.companion(p)
        Matrix([
        [0, 0, 0, 0, -c0],
        [1, 0, 0, 0, -c1],
        [0, 1, 0, 0, -c2],
        [0, 0, 1, 0, -c3],
        [0, 0, 0, 1, -c4]])
        """
        poly = kls._sympify(poly)
        if not isinstance(poly, Poly):
            raise ValueError("{} must be a Poly instance.".format(poly))
        if not poly.is_monic:
            raise ValueError("{} must be a monic polynomial.".format(poly))
        if not poly.is_univariate:
            raise ValueError(
                "{} must be a univariate polynomial.".format(poly))

        size = poly.degree()

        if size < 1:
            raise ValueError(f"{poly} must have degree not less than 1.")
        assert isinstance(size, int)

        coeffs = poly.all_coeffs()
        def entry(i, j):
            if j == size - 1:
                return -coeffs[-1 - i]
            elif i == j + 1:
                return kls.one
            return kls.zero
        return kls._new(size, size, entry)

    @overload
    @classmethod
    def wilkinson(kls, n: int, *, cls: None = None) -> tuple[Self, Self]:
        ...

    @overload
    @classmethod
    def wilkinson(kls, n: int, *, cls: type[Tmat]) -> tuple[Tmat, Tmat]:
        ...

    @classmethod
    def wilkinson(kls, n: int, *, cls: type[Tmat] | None = None
                  ) -> tuple[Self, Self] | tuple[Tmat, Tmat]:
        """Returns two square Wilkinson Matrix of size 2*n + 1
        $W_{2n + 1}^-, W_{2n + 1}^+ =$ Wilkinson(n)

        Examples
        ========

        >>> from sympy import Matrix
        >>> wminus, wplus = Matrix.wilkinson(3)
        >>> wminus
        Matrix([
        [-3,  1,  0, 0, 0, 0, 0],
        [ 1, -2,  1, 0, 0, 0, 0],
        [ 0,  1, -1, 1, 0, 0, 0],
        [ 0,  0,  1, 0, 1, 0, 0],
        [ 0,  0,  0, 1, 1, 1, 0],
        [ 0,  0,  0, 0, 1, 2, 1],
        [ 0,  0,  0, 0, 0, 1, 3]])
        >>> wplus
        Matrix([
        [3, 1, 0, 0, 0, 0, 0],
        [1, 2, 1, 0, 0, 0, 0],
        [0, 1, 1, 1, 0, 0, 0],
        [0, 0, 1, 0, 1, 0, 0],
        [0, 0, 0, 1, 1, 1, 0],
        [0, 0, 0, 0, 1, 2, 1],
        [0, 0, 0, 0, 0, 1, 3]])

        References
        ==========

        .. [1] https://blogs.mathworks.com/cleve/2013/04/15/wilkinsons-matrices-2/
        .. [2] J. H. Wilkinson, The Algebraic Eigenvalue Problem, Claredon Press, Oxford, 1965, 662 pp.

        """
        klass = kls if cls is None else cls
        n = as_int(n)
        return klass._eval_wilkinson(n)

    # The RepMatrix subclass uses more efficient sparse implementations of
    # _eval_iter_values and other things.

    def _eval_iter_values(self) -> Iterator[Expr]:
        return (i for i in self if i is not S.Zero) # type: ignore

    def _eval_values(self) -> list[Expr]:
        return list(self.iter_values())

    def _eval_iter_items(self) -> Iterator[tuple[tuple[int, int], Expr]]:
        for i in range(self.rows):
            for j in range(self.cols):
                if self[i, j]:
                    yield (i, j), self[i, j]

    @overload
    def _eval_atoms(self) -> set[Basic]: ...
    @overload
    def _eval_atoms(self, *types: Tbasic | type[Tbasic]) -> set[Tbasic]: ...

    def _eval_atoms(self, *types: Tbasic | type[Tbasic]) -> set[Basic] | set[Tbasic]:
        values = self.values()
        if len(values) < self.rows * self.cols:
            values.append(S.Zero)
        return set().union(*[v.atoms(*types) for v in values])

    def _eval_free_symbols(self) -> set[Basic]:
        return set().union(*(i.free_symbols for i in set(self.values())))

    def _eval_has(self, *patterns: SExpr | Basic | type[Basic]) -> bool:
        return any(a.has(*patterns) for a in self.iter_values())

    def _eval_is_symbolic(self) -> bool:
        return self.has(Symbol)

    # _eval_is_hermitian is called by some general SymPy
    # routines and has a different *args signature.  Make
    # sure the names don't clash by adding `_matrix_` in name.
    def _eval_is_matrix_hermitian(self, simpfunc: Callable[[Expr], Expr]) -> bool | None:
        herm = lambda i, j: simpfunc(self[i, j] - self[j, i].adjoint()).is_zero
        return fuzzy_and(herm(i, j) for (i, j), v in self.iter_items())

    def _eval_is_zero_matrix(self) -> bool | None:
        return fuzzy_and(v.is_zero for v in self.iter_values())

    def _eval_is_Identity(self) -> bool:
        one = self.one
        zero = self.zero
        ident = lambda i, j, v: v is one if i == j else v is zero
        return all(ident(i, j, v) for (i, j), v in self.iter_items())

    def _eval_is_diagonal(self) -> bool | None:
        return fuzzy_and(v.is_zero for (i, j), v in self.iter_items() if i != j)

    def _eval_is_lower(self) -> bool:
        return all(v.is_zero for (i, j), v in self.iter_items() if i < j)

    def _eval_is_upper(self) -> bool:
        return all(v.is_zero for (i, j), v in self.iter_items() if i > j)

    def _eval_is_lower_hessenberg(self) -> bool:
        return all(v.is_zero for (i, j), v in self.iter_items() if i + 1 < j)

    def _eval_is_upper_hessenberg(self) -> bool:
        return all(v.is_zero for (i, j), v in self.iter_items() if i > j + 1)

    def _eval_is_symmetric(self, simpfunc) -> bool | None:
        sym = lambda i, j: simpfunc(self[i, j] - self[j, i]).is_zero
        return fuzzy_and(sym(i, j) for (i, j), _ in self.iter_items())

    def _eval_is_anti_symmetric(self, simpfunc) -> bool | None:
        anti = lambda i, j: simpfunc(self[i, j] + self[j, i]).is_zero
        return fuzzy_and(anti(i, j) for (i, j), _ in self.iter_items())

    def _has_positive_diagonals(self) -> bool | None:
        diagonal_entries = (self[i, i] for i in range(self.rows))
        return fuzzy_and(x.is_positive for x in diagonal_entries)

    def _has_nonnegative_diagonals(self) -> bool | None:
        diagonal_entries = (self[i, i] for i in range(self.rows))
        return fuzzy_and(x.is_nonnegative for x in diagonal_entries)

    @overload
    def atoms(self) -> set[Basic]: ...
    @overload
    def atoms(self, *types: Tbasic | type[Tbasic]) -> set[Tbasic]: ...

    def atoms(self, *types: Tbasic | type[Tbasic]) -> set[Basic] | set[Tbasic]:
        """Returns the atoms that form the current object.

        Examples
        ========

        >>> from sympy.abc import x, y
        >>> from sympy import Matrix
        >>> Matrix([[x]])
        Matrix([[x]])
        >>> _.atoms()
        {x}
        >>> Matrix([[x, y], [y, x]])
        Matrix([
        [x, y],
        [y, x]])
        >>> _.atoms()
        {x, y}
        """
        types = tuple(t if isinstance(t, type) else type(t) for t in types)
        if not types:
            # XXX: .atoms(Atom) is not the same as .atoms()
            # This should be changed.
            return self._eval_atoms(Atom) # type: ignore
        else:
            return self._eval_atoms(*types)

    @property
    def free_symbols(self) -> set[Basic]:
        """Returns the free symbols within the matrix.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy import Matrix
        >>> Matrix([[x], [1]]).free_symbols
        {x}
        """
        return self._eval_free_symbols()

    def has(self, *patterns: SExpr | Basic | type[Basic]) -> bool:
        """Test whether any subexpression matches any of the patterns.

        Examples
        ========

        >>> from sympy import Matrix, SparseMatrix, Float
        >>> from sympy.abc import x, y
        >>> A = Matrix(((1, x), (0.2, 3)))
        >>> B = SparseMatrix(((1, x), (0.2, 3)))
        >>> A.has(x)
        True
        >>> A.has(y)
        False
        >>> A.has(Float)
        True
        >>> B.has(x)
        True
        >>> B.has(y)
        False
        >>> B.has(Float)
        True
        """
        return self._eval_has(*patterns)

    def is_anti_symmetric(self, simplify: bool | Callable[[Expr], Expr] = True) -> bool | None:
        """Check if matrix M is an antisymmetric matrix,
        that is, M is a square matrix with all M[i, j] == -M[j, i].

        When ``simplify=True`` (default), the sum M[i, j] + M[j, i] is
        simplified before testing to see if it is zero. By default,
        the SymPy simplify function is used. To use a custom function
        set simplify to a function that accepts a single argument which
        returns a simplified expression. To skip simplification, set
        simplify to False but note that although this will be faster,
        it may induce false negatives.

        Examples
        ========

        >>> from sympy import Matrix, symbols
        >>> m = Matrix(2, 2, [0, 1, -1, 0])
        >>> m
        Matrix([
        [ 0, 1],
        [-1, 0]])
        >>> m.is_anti_symmetric()
        True
        >>> x, y = symbols('x y')
        >>> m = Matrix(2, 3, [0, 0, x, -y, 0, 0])
        >>> m
        Matrix([
        [ 0, 0, x],
        [-y, 0, 0]])
        >>> m.is_anti_symmetric()
        False

        >>> from sympy.abc import x, y
        >>> m = Matrix(3, 3, [0, x**2 + 2*x + 1, y,
        ...                   -(x + 1)**2, 0, x*y,
        ...                   -y, -x*y, 0])

        Simplification of matrix elements is done by default so even
        though two elements which should be equal and opposite would not
        pass an equality test, the matrix is still reported as
        anti-symmetric:

        >>> m[0, 1] == -m[1, 0]
        False
        >>> m.is_anti_symmetric()
        True

        If ``simplify=False`` is used for the case when a Matrix is already
        simplified, this will speed things up. Here, we see that without
        simplification the matrix does not appear anti-symmetric:

        >>> print(m.is_anti_symmetric(simplify=False))
        None

        But if the matrix were already expanded, then it would appear
        anti-symmetric and simplification in the is_anti_symmetric routine
        is not needed:

        >>> m = m.expand()
        >>> m.is_anti_symmetric(simplify=False)
        True
        """
        # accept custom simplification
        simpfunc: Callable[[Any], Any]
        if not isfunction(simplify):
            simpfunc = _utilities_simplify if simplify else lambda x: x
        else:
            simpfunc = simplify

        if not self.is_square:
            return False
        return self._eval_is_anti_symmetric(simpfunc)

    def is_diagonal(self) -> bool | None:
        """Check if matrix is diagonal,
        that is matrix in which the entries outside the main diagonal are all zero.

        Examples
        ========

        >>> from sympy import Matrix, diag
        >>> m = Matrix(2, 2, [1, 0, 0, 2])
        >>> m
        Matrix([
        [1, 0],
        [0, 2]])
        >>> m.is_diagonal()
        True

        >>> m = Matrix(2, 2, [1, 1, 0, 2])
        >>> m
        Matrix([
        [1, 1],
        [0, 2]])
        >>> m.is_diagonal()
        False

        >>> m = diag(1, 2, 3)
        >>> m
        Matrix([
        [1, 0, 0],
        [0, 2, 0],
        [0, 0, 3]])
        >>> m.is_diagonal()
        True

        See Also
        ========

        is_lower
        is_upper
        sympy.matrices.matrixbase.MatrixBase.is_diagonalizable
        diagonalize
        """
        return self._eval_is_diagonal()

    @property
    def is_weakly_diagonally_dominant(self) -> bool | None:
        r"""Tests if the matrix is row weakly diagonally dominant.

        Explanation
        ===========

        A $n, n$ matrix $A$ is row weakly diagonally dominant if

        .. math::
            \left|A_{i, i}\right| \ge \sum_{j = 0, j \neq i}^{n-1}
            \left|A_{i, j}\right| \quad {\text{for all }}
            i \in \{ 0, ..., n-1 \}

        Examples
        ========

        >>> from sympy import Matrix
        >>> A = Matrix([[3, -2, 1], [1, -3, 2], [-1, 2, 4]])
        >>> A.is_weakly_diagonally_dominant
        True

        >>> A = Matrix([[-2, 2, 1], [1, 3, 2], [1, -2, 0]])
        >>> A.is_weakly_diagonally_dominant
        False

        >>> A = Matrix([[-4, 2, 1], [1, 6, 2], [1, -2, 5]])
        >>> A.is_weakly_diagonally_dominant
        True

        Notes
        =====

        If you want to test whether a matrix is column diagonally
        dominant, you can apply the test after transposing the matrix.
        """
        if not self.is_square:
            return False

        rows, cols = self.shape

        def test_row(i):
            summation = self.zero
            for j in range(cols):
                if i != j:
                    summation += Abs(self[i, j])
            return (Abs(self[i, i]) - summation).is_nonnegative

        return fuzzy_and(test_row(i) for i in range(rows))

    @property
    def is_strongly_diagonally_dominant(self) -> bool | None:
        r"""Tests if the matrix is row strongly diagonally dominant.

        Explanation
        ===========

        A $n, n$ matrix $A$ is row strongly diagonally dominant if

        .. math::
            \left|A_{i, i}\right| > \sum_{j = 0, j \neq i}^{n-1}
            \left|A_{i, j}\right| \quad {\text{for all }}
            i \in \{ 0, ..., n-1 \}

        Examples
        ========

        >>> from sympy import Matrix
        >>> A = Matrix([[3, -2, 1], [1, -3, 2], [-1, 2, 4]])
        >>> A.is_strongly_diagonally_dominant
        False

        >>> A = Matrix([[-2, 2, 1], [1, 3, 2], [1, -2, 0]])
        >>> A.is_strongly_diagonally_dominant
        False

        >>> A = Matrix([[-4, 2, 1], [1, 6, 2], [1, -2, 5]])
        >>> A.is_strongly_diagonally_dominant
        True

        Notes
        =====

        If you want to test whether a matrix is column diagonally
        dominant, you can apply the test after transposing the matrix.
        """
        if not self.is_square:
            return False

        rows, cols = self.shape

        def test_row(i):
            summation = self.zero
            for j in range(cols):
                if i != j:
                    summation += Abs(self[i, j])
            return (Abs(self[i, i]) - summation).is_positive

        return fuzzy_and(test_row(i) for i in range(rows))

    @property
    def is_hermitian(self) -> bool | None:
        """Checks if the matrix is Hermitian.

        In a Hermitian matrix element i,j is the complex conjugate of
        element j,i.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy import I
        >>> from sympy.abc import x
        >>> a = Matrix([[1, I], [-I, 1]])
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
        if not self.is_square:
            return False

        return self._eval_is_matrix_hermitian(_utilities_simplify)

    @property
    def is_Identity(self) -> bool:
        if not self.is_square:
            return False
        return self._eval_is_Identity()

    @property
    def is_lower_hessenberg(self) -> bool:
        r"""Checks if the matrix is in the lower-Hessenberg form.

        The lower hessenberg matrix has zero entries
        above the first superdiagonal.

        Examples
        ========

        >>> from sympy import Matrix
        >>> a = Matrix([[1, 2, 0, 0], [5, 2, 3, 0], [3, 4, 3, 7], [5, 6, 1, 1]])
        >>> a
        Matrix([
        [1, 2, 0, 0],
        [5, 2, 3, 0],
        [3, 4, 3, 7],
        [5, 6, 1, 1]])
        >>> a.is_lower_hessenberg
        True

        See Also
        ========

        is_upper_hessenberg
        is_lower
        """
        return self._eval_is_lower_hessenberg()

    @property
    def is_lower(self) -> bool:
        """Check if matrix is a lower triangular matrix. True can be returned
        even if the matrix is not square.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(2, 2, [1, 0, 0, 1])
        >>> m
        Matrix([
        [1, 0],
        [0, 1]])
        >>> m.is_lower
        True

        >>> m = Matrix(4, 3, [0, 0, 0, 2, 0, 0, 1, 4, 0, 6, 6, 5])
        >>> m
        Matrix([
        [0, 0, 0],
        [2, 0, 0],
        [1, 4, 0],
        [6, 6, 5]])
        >>> m.is_lower
        True

        >>> from sympy.abc import x, y
        >>> m = Matrix(2, 2, [x**2 + y, y**2 + x, 0, x + y])
        >>> m
        Matrix([
        [x**2 + y, x + y**2],
        [       0,    x + y]])
        >>> m.is_lower
        False

        See Also
        ========

        is_upper
        is_diagonal
        is_lower_hessenberg
        """
        return self._eval_is_lower()

    @property
    def is_square(self) -> bool:
        """Checks if a matrix is square.

        A matrix is square if the number of rows equals the number of columns.
        The empty matrix is square by definition, since the number of rows and
        the number of columns are both zero.

        Examples
        ========

        >>> from sympy import Matrix
        >>> a = Matrix([[1, 2, 3], [4, 5, 6]])
        >>> b = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> c = Matrix([])
        >>> a.is_square
        False
        >>> b.is_square
        True
        >>> c.is_square
        True
        """
        return self.rows == self.cols

    def is_symbolic(self) -> bool:
        """Checks if any elements contain Symbols.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.abc import x, y
        >>> M = Matrix([[x, y], [1, 0]])
        >>> M.is_symbolic()
        True

        """
        return self._eval_is_symbolic()

    def is_symmetric(self, simplify: Callable[[Expr], Expr] | bool = True) -> bool | None:
        """Check if matrix is symmetric matrix,
        that is square matrix and is equal to its transpose.

        By default, simplifications occur before testing symmetry.
        They can be skipped using 'simplify=False'; while speeding things a bit,
        this may however induce false negatives.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(2, 2, [0, 1, 1, 2])
        >>> m
        Matrix([
        [0, 1],
        [1, 2]])
        >>> m.is_symmetric()
        True

        >>> m = Matrix(2, 2, [0, 1, 2, 0])
        >>> m
        Matrix([
        [0, 1],
        [2, 0]])
        >>> m.is_symmetric()
        False

        >>> m = Matrix(2, 3, [0, 0, 0, 0, 0, 0])
        >>> m
        Matrix([
        [0, 0, 0],
        [0, 0, 0]])
        >>> m.is_symmetric()
        False

        >>> from sympy.abc import x, y
        >>> m = Matrix(3, 3, [1, x**2 + 2*x + 1, y, (x + 1)**2, 2, 0, y, 0, 3])
        >>> m
        Matrix([
        [         1, x**2 + 2*x + 1, y],
        [(x + 1)**2,              2, 0],
        [         y,              0, 3]])
        >>> m.is_symmetric()
        True

        If the matrix is already simplified, you may speed-up is_symmetric()
        test by using 'simplify=False'.

        >>> bool(m.is_symmetric(simplify=False))
        False
        >>> m1 = m.expand()
        >>> m1.is_symmetric(simplify=False)
        True
        """
        simpfunc: Callable[[Any], Any]
        if not isfunction(simplify):
            simpfunc = _utilities_simplify if simplify else lambda x: x
        else:
            simpfunc = simplify

        if not self.is_square:
            return False

        return self._eval_is_symmetric(simpfunc)

    @property
    def is_upper_hessenberg(self) -> bool:
        """Checks if the matrix is the upper-Hessenberg form.

        The upper hessenberg matrix has zero entries
        below the first subdiagonal.

        Examples
        ========

        >>> from sympy import Matrix
        >>> a = Matrix([[1, 4, 2, 3], [3, 4, 1, 7], [0, 2, 3, 4], [0, 0, 1, 3]])
        >>> a
        Matrix([
        [1, 4, 2, 3],
        [3, 4, 1, 7],
        [0, 2, 3, 4],
        [0, 0, 1, 3]])
        >>> a.is_upper_hessenberg
        True

        See Also
        ========

        is_lower_hessenberg
        is_upper
        """
        return self._eval_is_upper_hessenberg()

    @property
    def is_upper(self) -> bool:
        """Check if matrix is an upper triangular matrix. True can be returned
        even if the matrix is not square.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(2, 2, [1, 0, 0, 1])
        >>> m
        Matrix([
        [1, 0],
        [0, 1]])
        >>> m.is_upper
        True

        >>> m = Matrix(4, 3, [5, 1, 9, 0, 4, 6, 0, 0, 5, 0, 0, 0])
        >>> m
        Matrix([
        [5, 1, 9],
        [0, 4, 6],
        [0, 0, 5],
        [0, 0, 0]])
        >>> m.is_upper
        True

        >>> m = Matrix(2, 3, [4, 2, 5, 6, 1, 1])
        >>> m
        Matrix([
        [4, 2, 5],
        [6, 1, 1]])
        >>> m.is_upper
        False

        See Also
        ========

        is_lower
        is_diagonal
        is_upper_hessenberg
        """
        return self._eval_is_upper()

    @property
    def is_zero_matrix(self) -> bool | None:
        """Checks if a matrix is a zero matrix.

        A matrix is zero if every element is zero.  A matrix need not be square
        to be considered zero.  The empty matrix is zero by the principle of
        vacuous truth.  For a matrix that may or may not be zero (e.g.
        contains a symbol), this will be None

        Examples
        ========

        >>> from sympy import Matrix, zeros
        >>> from sympy.abc import x
        >>> a = Matrix([[0, 0], [0, 0]])
        >>> b = zeros(3, 4)
        >>> c = Matrix([[0, 1], [0, 0]])
        >>> d = Matrix([])
        >>> e = Matrix([[x, 0], [0, 0]])
        >>> a.is_zero_matrix
        True
        >>> b.is_zero_matrix
        True
        >>> c.is_zero_matrix
        False
        >>> d.is_zero_matrix
        True
        >>> e.is_zero_matrix
        """
        return self._eval_is_zero_matrix()

    def values(self) -> list[Expr]:
        """Return non-zero values of self.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix([[0, 1], [2, 3]])
        >>> m.values()
        [1, 2, 3]

        See Also
        ========

        iter_values
        tolist
        flat
        """
        return self._eval_values()

    def iter_values(self) -> Iterator[Expr]:
        """
        Iterate over non-zero values of self.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix([[0, 1], [2, 3]])
        >>> list(m.iter_values())
        [1, 2, 3]

        See Also
        ========

        values
        """
        return self._eval_iter_values()

    def iter_items(self) -> Iterator[tuple[tuple[int, int], Expr]]:
        """Iterate over indices and values of nonzero items.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix([[0, 1], [2, 3]])
        >>> list(m.iter_items())
        [((0, 1), 1), ((1, 0), 2), ((1, 1), 3)]

        See Also
        ========

        iter_values
        todok
        """
        return self._eval_iter_items()

    def _eval_adjoint(self) -> Self:
        return self.transpose().applyfunc(lambda x: x.adjoint())

    def _eval_applyfunc(self, f: Callable[[Expr], Expr]) -> Self:
        cols = self.cols
        size = self.rows*self.cols

        dok = self.todok()
        valmap = {v: f(v) for v in dok.values()}

        if len(dok) < size and ((fzero := f(S.Zero)) is not S.Zero):
            out_flat = [fzero]*size
            for (i, j), v in dok.items():
                out_flat[i*cols + j] = valmap[v]
            out = self._new(self.rows, self.cols, out_flat)
        else:
            fdok = {ij: valmap[v] for ij, v in dok.items()}
            out = self.from_dok(self.rows, self.cols, fdok)

        return out

    def _eval_as_real_imag(self) -> tuple[Self, Self]:
        return (self.applyfunc(re), self.applyfunc(im))

    def _eval_conjugate(self) -> Self:
        return self.applyfunc(lambda x: x.conjugate())

    def _eval_permute_cols(self, perm) -> Self:
        # apply the permutation to a list
        mapping = list(perm)

        def entry(i, j):
            return self[i, mapping[j]]

        return self._new(self.rows, self.cols, entry)

    def _eval_permute_rows(self, perm) -> Self:
        # apply the permutation to a list
        mapping = list(perm)

        def entry(i, j):
            return self[mapping[i], j]

        return self._new(self.rows, self.cols, entry)

    def _eval_trace(self) -> Expr:
        return Add(*(self[i, i] for i in range(self.rows)))

    def _eval_transpose(self) -> Self:
        return self._new(self.cols, self.rows, lambda i, j: self[j, i])

    def adjoint(self) -> Self:
        """Conjugate transpose or Hermitian conjugation."""
        return self._eval_adjoint()

    def applyfunc(self, f: Callable[[Expr], Expr]) -> Self:
        """Apply a function to each element of the matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix(2, 2, lambda i, j: i*2+j)
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

        return self._eval_applyfunc(f)

    def as_real_imag(self, deep: bool = True, **hints: Any) -> tuple[Self, Self]:
        """Returns a tuple containing the (real, imaginary) part of matrix."""
        # XXX: Ignoring deep and hints...
        return self._eval_as_real_imag()

    def conjugate(self) -> Self:
        """Return the by-element conjugation.

        Examples
        ========

        >>> from sympy import SparseMatrix, I
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
        sympy.matrices.matrixbase.MatrixBase.D: Dirac conjugation
        """
        return self._eval_conjugate()

    def doit(self, **hints) -> Self:
        return self.applyfunc(lambda x: x.doit(**hints))

    def evalf(self, n: int = 15,
              subs: Mapping[SExpr, SExpr] | None = None,
              maxn: int = 100,
              chop: bool | int = False,
              strict: bool = False,
              quad: Literal['osc'] | str | None = None,
              verbose: bool = False
        ) -> Self:
        """Apply evalf() to each element of self."""
        f = lambda i: i.evalf(n, subs=subs, maxn=maxn, chop=chop,
                              strict=strict, quad=quad, verbose=verbose)
        return self.applyfunc(f)

    def n(self, n: int = 15,
              subs: Mapping[SExpr, SExpr] | None = None,
              maxn: int = 100,
              chop: bool | int = False,
              strict: bool = False,
              quad: Literal['osc'] | str | None = None,
              verbose: bool = False
        ) -> Self:
        """Apply evalf() to each element of self."""
        return self.evalf(n, subs=subs, maxn=maxn, chop=chop,
                          strict=strict, quad=quad, verbose=verbose)

    def expand(self,
           deep: bool = True,
           modulus: int | None = None,
           power_base: bool = True,
           power_exp: bool = True,
           mul: bool = True,
           log: bool = True,
           multinomial: bool = True,
           basic: bool = True,
           **hints: Any,
       ) -> Self:
        """Apply core.function.expand to each entry of the matrix.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy import Matrix
        >>> Matrix(1, 1, [x*(x+1)])
        Matrix([[x*(x + 1)]])
        >>> _.expand()
        Matrix([[x**2 + x]])

        """
        return self.applyfunc(lambda x: x.expand(
            deep, modulus, power_base, power_exp, mul, log, multinomial, basic,
            **hints))

    @property
    def H(self) -> Self:
        """Return Hermite conjugate.

        Examples
        ========

        >>> from sympy import Matrix, I
        >>> m = Matrix((0, 1 + I, 2, 3))
        >>> m
        Matrix([
        [    0],
        [1 + I],
        [    2],
        [    3]])
        >>> m.H
        Matrix([[0, 1 - I, 2, 3]])

        See Also
        ========

        conjugate: By-element conjugation
        sympy.matrices.matrixbase.MatrixBase.D: Dirac conjugation
        """
        return self.adjoint()

    def permute(self,
            perm: list[int] | list[list[int]] | Permutation,
            orientation: Literal['rows', 'cols'] = 'rows',
            direction: Literal['forward', 'backward'] = 'forward',
        ) -> Self:
        r"""Permute the rows or columns of a matrix by the given list of
        swaps.

        Parameters
        ==========

        perm : Permutation, list, or list of lists
            A representation for the permutation.

            If it is ``Permutation``, it is used directly with some
            resizing with respect to the matrix size.

            If it is specified as list of lists,
            (e.g., ``[[0, 1], [0, 2]]``), then the permutation is formed
            from applying the product of cycles. The direction how the
            cyclic product is applied is described in below.

            If it is specified as a list, the list should represent
            an array form of a permutation. (e.g., ``[1, 2, 0]``) which
            would would form the swapping function
            `0 \mapsto 1, 1 \mapsto 2, 2\mapsto 0`.

        orientation : 'rows', 'cols'
            A flag to control whether to permute the rows or the columns

        direction : 'forward', 'backward'
            A flag to control whether to apply the permutations from
            the start of the list first, or from the back of the list
            first.

            For example, if the permutation specification is
            ``[[0, 1], [0, 2]]``,

            If the flag is set to ``'forward'``, the cycle would be
            formed as `0 \mapsto 2, 2 \mapsto 1, 1 \mapsto 0`.

            If the flag is set to ``'backward'``, the cycle would be
            formed as `0 \mapsto 1, 1 \mapsto 2, 2 \mapsto 0`.

            If the argument ``perm`` is not in a form of list of lists,
            this flag takes no effect.

        Examples
        ========

        >>> from sympy import eye
        >>> M = eye(3)
        >>> M.permute([[0, 1], [0, 2]], orientation='rows', direction='forward')
        Matrix([
        [0, 0, 1],
        [1, 0, 0],
        [0, 1, 0]])

        >>> from sympy import eye
        >>> M = eye(3)
        >>> M.permute([[0, 1], [0, 2]], orientation='rows', direction='backward')
        Matrix([
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 0]])

        Notes
        =====

        If a bijective function
        `\sigma : \mathbb{N}_0 \rightarrow \mathbb{N}_0` denotes the
        permutation.

        If the matrix `A` is the matrix to permute, represented as
        a horizontal or a vertical stack of vectors:

        .. math::
            A =
            \begin{bmatrix}
            a_0 \\ a_1 \\ \vdots \\ a_{n-1}
            \end{bmatrix} =
            \begin{bmatrix}
            \alpha_0 & \alpha_1 & \cdots & \alpha_{n-1}
            \end{bmatrix}

        If the matrix `B` is the result, the permutation of matrix rows
        is defined as:

        .. math::
            B := \begin{bmatrix}
            a_{\sigma(0)} \\ a_{\sigma(1)} \\ \vdots \\ a_{\sigma(n-1)}
            \end{bmatrix}

        And the permutation of matrix columns is defined as:

        .. math::
            B := \begin{bmatrix}
            \alpha_{\sigma(0)} & \alpha_{\sigma(1)} &
            \cdots & \alpha_{\sigma(n-1)}
            \end{bmatrix}
        """
        from sympy.combinatorics import Permutation

        # allow british variants and `columns`
        if direction == 'forwards':
            direction = 'forward'
        if direction == 'backwards':
            direction = 'backward'
        if orientation == 'columns':
            orientation = 'cols'

        if direction not in ('forward', 'backward'):
            raise TypeError("direction='{}' is an invalid kwarg. "
                            "Try 'forward' or 'backward'".format(direction))
        if orientation not in ('rows', 'cols'):
            raise TypeError("orientation='{}' is an invalid kwarg. "
                            "Try 'rows' or 'cols'".format(orientation))

        if not isinstance(perm, (Permutation, Iterable)):
            raise ValueError(
                "{} must be a list, a list of lists, "
                "or a SymPy permutation object.".format(perm))

        # ensure all swaps are in range
        max_index = self.rows if orientation == 'rows' else self.cols
        if not all(0 <= t <= max_index for t in flatten(list(perm))):
            raise IndexError("`swap` indices out of range.")

        if perm and not isinstance(perm, Permutation) and \
            isinstance(perm[0], Iterable):
            if direction == 'forward':
                perm = list(reversed(perm)) # type: ignore
            perm = Permutation(perm, size=max_index+1)
        else:
            perm = Permutation(perm, size=max_index+1)

        if orientation == 'rows':
            return self._eval_permute_rows(perm)
        if orientation == 'cols':
            return self._eval_permute_cols(perm)

    def permute_cols(self,
            swaps: list[int] | list[list[int]] | Permutation,
            direction: Literal['forward', 'backward'] = 'forward',
         ) -> Self:
        """Alias for
        ``self.permute(swaps, orientation='cols', direction=direction)``

        See Also
        ========

        permute
        """
        return self.permute(swaps, orientation='cols', direction=direction)

    def permute_rows(self,
            swaps: list[int] | list[list[int]] | Permutation,
            direction: Literal['forward', 'backward'] = 'forward',
         ) -> Self:
        """Alias for
        ``self.permute(swaps, orientation='rows', direction=direction)``

        See Also
        ========

        permute
        """
        return self.permute(swaps, orientation='rows', direction=direction)

    def refine(self, assumptions: Boolean | bool = True) -> Self:
        """Apply refine to each element of the matrix.

        Examples
        ========

        >>> from sympy import Symbol, Matrix, Abs, sqrt, Q
        >>> x = Symbol('x')
        >>> Matrix([[Abs(x)**2, sqrt(x**2)],[sqrt(x**2), Abs(x)**2]])
        Matrix([
        [ Abs(x)**2, sqrt(x**2)],
        [sqrt(x**2),  Abs(x)**2]])
        >>> _.refine(Q.real(x))
        Matrix([
        [  x**2, Abs(x)],
        [Abs(x),   x**2]])

        """
        return self.applyfunc(lambda x: refine(x, assumptions))

    def replace(self, F, G, map: bool=False, simultaneous: bool=True, exact=None):
        """Replaces Function F in Matrix entries with Function G.

        Examples
        ========

        >>> from sympy import symbols, Function, Matrix
        >>> F, G = symbols('F, G', cls=Function)
        >>> M = Matrix(2, 2, lambda i, j: F(i+j)) ; M
        Matrix([
        [F(0), F(1)],
        [F(1), F(2)]])
        >>> N = M.replace(F,G)
        >>> N
        Matrix([
        [G(0), G(1)],
        [G(1), G(2)]])
        """
        kwargs = {'map': map, 'simultaneous': simultaneous, 'exact': exact}

        if map:

            d = {}
            def func(eij):
                eij, dij = eij.replace(F, G, **kwargs)
                d.update(dij)
                return eij

            M = self.applyfunc(func)
            return M, d

        else:
            return self.applyfunc(lambda i: i.replace(F, G, **kwargs)) # type: ignore

    def rot90(self, k: int=1) -> Self:
        """Rotates Matrix by 90 degrees

        Parameters
        ==========

        k : int
            Specifies how many times the matrix is rotated by 90 degrees
            (clockwise when positive, counter-clockwise when negative).

        Examples
        ========

        >>> from sympy import Matrix, symbols
        >>> A = Matrix(2, 2, symbols('a:d'))
        >>> A
        Matrix([
        [a, b],
        [c, d]])

        Rotating the matrix clockwise one time:

        >>> A.rot90(1)
        Matrix([
        [c, a],
        [d, b]])

        Rotating the matrix anticlockwise two times:

        >>> A.rot90(-2)
        Matrix([
        [d, c],
        [b, a]])
        """
        mod = k%4
        if mod == 0:
            return self
        elif mod == 1:
            return self[::-1, ::].T
        elif mod == 2:
            return self[::-1, ::-1]
        elif mod == 3:
            return self[::, ::-1].T
        else:
            assert False

    def simplify(self, **kwargs: Any) -> Self:
        """Apply simplify to each element of the matrix.

        Examples
        ========

        >>> from sympy.abc import x, y
        >>> from sympy import SparseMatrix, sin, cos
        >>> SparseMatrix(1, 1, [x*sin(y)**2 + x*cos(y)**2])
        Matrix([[x*sin(y)**2 + x*cos(y)**2]])
        >>> _.simplify()
        Matrix([[x]])
        """
        return self.applyfunc(lambda x: x.simplify(**kwargs))

    @overload
    def subs(self, arg1: Mapping[SBasic, SBasic], arg2: None=None, **kwargs: Any) -> Self: ...
    @overload
    def subs(self, arg1: Iterable[tuple[SBasic, SBasic]], arg2: None=None, **kwargs: Any) -> Self: ...
    @overload
    def subs(self, arg1: SBasic, arg2: SBasic, **kwargs: Any) -> Self: ...

    def subs(self,
             arg1: Mapping[SBasic, SBasic] | Iterable[tuple[SBasic, SBasic]] | SBasic,
             arg2: SBasic | None = None,
             **kwargs: Any
        ) -> Self:
        """Return a new matrix with subs applied to each entry.

        Examples
        ========

        >>> from sympy.abc import x, y
        >>> from sympy import SparseMatrix, Matrix
        >>> SparseMatrix(1, 1, [x])
        Matrix([[x]])
        >>> _.subs(x, y)
        Matrix([[y]])
        >>> Matrix(_).subs(y, x)
        Matrix([[x]])
        """
        if arg2 is None and not isinstance(arg1, (dict, set)) and iter(arg1) and not is_sequence(arg1): # type: ignore
            arg1 = list(arg1) # type: ignore

        return self.applyfunc(lambda x: x.subs(arg1, arg2, **kwargs)) # type: ignore

    def trace(self) -> Expr:
        """
        Returns the trace of a square matrix i.e. the sum of the
        diagonal elements.

        Examples
        ========

        >>> from sympy import Matrix
        >>> A = Matrix(2, 2, [1, 2, 3, 4])
        >>> A.trace()
        5

        """
        if self.rows != self.cols:
            raise NonSquareMatrixError()
        return self._eval_trace()

    def transpose(self) -> Self:
        """
        Returns the transpose of the matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> A = Matrix(2, 2, [1, 2, 3, 4])
        >>> A.transpose()
        Matrix([
        [1, 3],
        [2, 4]])

        >>> from sympy import Matrix, I
        >>> m=Matrix(((1, 2+I), (3, 4)))
        >>> m
        Matrix([
        [1, 2 + I],
        [3,     4]])
        >>> m.transpose()
        Matrix([
        [    1, 3],
        [2 + I, 4]])
        >>> m.T == m.transpose()
        True

        See Also
        ========

        conjugate: By-element conjugation

        """
        return self._eval_transpose()

    @property
    def T(self) -> Self:
        '''Matrix transposition'''
        return self.transpose()

    @property
    def C(self) -> Self:
        '''By-element conjugation'''
        return self.conjugate()

    # should mirror core.basic.xreplace
    def xreplace(self, rule: Mapping[SBasic, SBasic]) -> Self:
        """Return a new matrix with xreplace applied to each entry.

        Examples
        ========

        >>> from sympy.abc import x, y
        >>> from sympy import SparseMatrix, Matrix
        >>> SparseMatrix(1, 1, [x])
        Matrix([[x]])
        >>> _.xreplace({x: y})
        Matrix([[y]])
        >>> Matrix(_).xreplace({y: x})
        Matrix([[x]])
        """
        return self.applyfunc(lambda x: x.xreplace(rule))

    def _eval_simplify(self, **kwargs: Any) -> Self:
        # XXX: We can't use self.simplify here as mutable subclasses will
        # override simplify and have it return None
        return self.applyfunc(lambda x: x.simplify(**kwargs))

    def _eval_trigsimp(self, **opts: Any) -> Self:
        from sympy.simplify.trigsimp import trigsimp
        return self.applyfunc(lambda x: trigsimp(x, **opts))

    def upper_triangular(self, k: int = 0) -> Self:
        """Return the elements on and above the kth diagonal of a matrix.
        If k is not specified then simply returns upper-triangular portion
        of a matrix

        Examples
        ========

        >>> from sympy import ones
        >>> A = ones(4)
        >>> A.upper_triangular()
        Matrix([
        [1, 1, 1, 1],
        [0, 1, 1, 1],
        [0, 0, 1, 1],
        [0, 0, 0, 1]])

        >>> A.upper_triangular(2)
        Matrix([
        [0, 0, 1, 1],
        [0, 0, 0, 1],
        [0, 0, 0, 0],
        [0, 0, 0, 0]])

        >>> A.upper_triangular(-1)
        Matrix([
        [1, 1, 1, 1],
        [1, 1, 1, 1],
        [0, 1, 1, 1],
        [0, 0, 1, 1]])

        """

        def entry(i, j):
            return self[i, j] if i + k <= j else self.zero

        return self._new(self.rows, self.cols, entry)

    def lower_triangular(self, k: int = 0) -> Self:
        """Return the elements on and below the kth diagonal of a matrix.
        If k is not specified then simply returns lower-triangular portion
        of a matrix

        Examples
        ========

        >>> from sympy import ones
        >>> A = ones(4)
        >>> A.lower_triangular()
        Matrix([
        [1, 0, 0, 0],
        [1, 1, 0, 0],
        [1, 1, 1, 0],
        [1, 1, 1, 1]])

        >>> A.lower_triangular(-2)
        Matrix([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [1, 0, 0, 0],
        [1, 1, 0, 0]])

        >>> A.lower_triangular(1)
        Matrix([
        [1, 1, 0, 0],
        [1, 1, 1, 0],
        [1, 1, 1, 1],
        [1, 1, 1, 1]])

        """

        def entry(i, j):
            return self[i, j] if i + k >= j else self.zero

        return self._new(self.rows, self.cols, entry)

    def _eval_Abs(self) -> Self:
        return self._new(self.rows, self.cols, lambda i, j: Abs(self[i, j]))

    def _eval_add(self, other: Self) -> Self:
        return self._new(self.rows, self.cols,
                         lambda i, j: self[i, j] + other[i, j])

    def _eval_matrix_mul(self, other: MatrixBase) -> Self:
        def entry(i, j):
            vec = [self[i,k]*other[k,j] for k in range(self.cols)]
            try:
                return Add(*vec)
            except (TypeError, SympifyError):
                # Some matrices don't work with `sum` or `Add`
                # They don't work with `sum` because `sum` tries to add `0`
                # Fall back to a safe way to multiply if the `Add` fails.
                return reduce(lambda a, b: a + b, vec)

        return self._new(self.rows, other.cols, entry)

    def _eval_matrix_mul_elementwise(self, other: MatrixBase) -> Self:
        return self._new(self.rows, self.cols, lambda i, j: self[i,j]*other[i,j])

    def _eval_matrix_rmul(self, other: MatrixBase) -> Self:
        def entry(i, j):
            return sum(other[i,k]*self[k,j] for k in range(other.cols))
        return self._new(other.rows, self.cols, entry)

    def _eval_pow_by_recursion(self, num: int) -> Self:
        if num == 1:
            return self

        if num % 2 == 1:
            a, b = self, self._eval_pow_by_recursion(num - 1)
        else:
            a = b = self._eval_pow_by_recursion(num // 2)

        return a.multiply(b)

    def _eval_pow_by_cayley(self, exp: int) -> Self:
        from sympy.discrete.recurrences import linrec_coeffs
        row = self.shape[0]
        p = self.charpoly()

        coeffs = (-p).all_coeffs()[1:]
        coeffs = linrec_coeffs(coeffs, exp)
        new_mat = self.eye(row)
        ans = self.zeros(row)

        for i in range(row):
            ans += coeffs[i]*new_mat
            new_mat *= self

        return ans

    def _eval_pow_by_recursion_dotprodsimp(self,
            num: int,
            prevsimp: list[bool] | None = None
        ) -> Self:

        if prevsimp is None:
            prevsimp = [True]*len(self)

        if num == 1:
            return self

        if num % 2 == 1:
            a, b = self, self._eval_pow_by_recursion_dotprodsimp(num - 1,
                    prevsimp=prevsimp)
        else:
            a = b = self._eval_pow_by_recursion_dotprodsimp(num // 2,
                    prevsimp=prevsimp)

        m     = a.multiply(b, dotprodsimp=False)
        lenm  = len(m)
        elems: list[Expr] = []

        for i in range(lenm):
            if prevsimp[i]:
                ei, prevsimp[i] = _dotprodsimp(m[i], withsimp=True)
            else:
                ei = m[i]
            elems.append(ei)

        return self._new(m.rows, m.cols, elems)

    def _eval_scalar_mul(self, other: SExpr) -> Self:
        return self._new(self.rows, self.cols, lambda i, j: self[i,j]*other)

    def _eval_scalar_rmul(self, other: SExpr) -> Self:
        return self._new(self.rows, self.cols, lambda i, j: other*self[i,j])

    def _eval_Mod(self, other: Expr) -> Self:
        return self._new(self.rows, self.cols, lambda i, j: Mod(self[i, j], other))

    # Python arithmetic functions
    def __abs__(self) -> Self:
        """Returns a new matrix with entry-wise absolute values."""
        return self._eval_Abs()

    @call_highest_priority('__radd__')
    def __add__(self: Tmat, other: Tmat) -> Tmat:
        """Return self + other, raising ShapeError if shapes do not match."""

        other_mat, _ = _coerce_operand(self, other)

        if not isinstance(other_mat, MatrixBase):
            return NotImplemented

        if self.shape != other_mat.shape:
            raise ShapeError(f"Matrix size mismatch: {self.shape} + {other.shape}.")

        # Unify matrix types
        a, b = self, other_mat
        if a.__class__ != classof(a, b):
            b, a = a, b

        return a._eval_add(b) # type: ignore

    @call_highest_priority('__rtruediv__')
    def __truediv__(self, other):
        return self * (self.one / other)

    @call_highest_priority('__rmatmul__')
    def __matmul__(self, other):
        self, other, T = _unify_with_other(self, other)

        if T != "is_matrix":
            return NotImplemented

        return self.__mul__(other)

    def __mod__(self, other) -> Self:
        return self.applyfunc(lambda x: x % other)

    @overload
    def __mul__(self, other: Self) -> Self: ...
    @overload
    def __mul__(self, other: MatrixBase) -> MatrixBase: ... # type: ignore
    @overload
    def __mul__(self, other: Expr) -> MatrixBase: ...

    @call_highest_priority('__rmul__')
    def __mul__(self, other: MatrixBase | Expr) -> MatrixBase | Expr:
        """Return self*other where other is either a scalar or a matrix
        of compatible dimensions.

        Examples
        ========

        >>> from sympy import Matrix
        >>> A = Matrix([[1, 2, 3], [4, 5, 6]])
        >>> 2*A == A*2 == Matrix([[2, 4, 6], [8, 10, 12]])
        True
        >>> B = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> A*B
        Matrix([
        [30, 36, 42],
        [66, 81, 96]])
        >>> B*A
        Traceback (most recent call last):
        ...
        ShapeError: Matrices size mismatch.
        >>>

        See Also
        ========

        matrix_multiply_elementwise
        """
        return self.multiply(other)

    @overload
    def multiply(self, other: Self, dotprodsimp: bool | None = None) -> Self: ...
    @overload
    def multiply(self, other: MatrixBase, dotprodsimp: bool | None = None) -> MatrixBase: ... # type: ignore
    @overload
    def multiply(self, other: Expr, dotprodsimp: bool | None = None) -> Expr: ...

    def multiply(self,
                 other: MatrixBase | Expr,
                 dotprodsimp: bool | None = None
             ) -> MatrixBase | Expr:
        """Same as __mul__() but with optional simplification.

        Parameters
        ==========

        dotprodsimp : bool, optional
            Specifies whether intermediate term algebraic simplification is used
            during matrix multiplications to control expression blowup and thus
            speed up calculation. Default is off.
        """

        isimpbool = _get_intermediate_simp_bool(False, dotprodsimp)

        self, other, T = _unify_with_other(self, other)

        if isinstance(other, MatrixBase):

            if self.shape[1] != other.shape[0]:
                raise ShapeError(f"Matrix size mismatch: {self.shape} * {other.shape}.")

            m = self._eval_matrix_mul(other)

            if isimpbool:
                m = m._new(m.rows, m.cols, [_dotprodsimp(e) for e in m.flat()])

            return m

        elif T == "possible_scalar":
            try:
                return self._eval_scalar_mul(other)
            except TypeError:
                return NotImplemented

        else:
            return NotImplemented

    def multiply_elementwise(self, other: MatrixBase) -> Self:
        """Return the Hadamard product (elementwise product) of A and B

        Examples
        ========

        >>> from sympy import Matrix
        >>> A = Matrix([[0, 1, 2], [3, 4, 5]])
        >>> B = Matrix([[1, 10, 100], [100, 10, 1]])
        >>> A.multiply_elementwise(B)
        Matrix([
        [  0, 10, 200],
        [300, 40,   5]])

        See Also
        ========

        sympy.matrices.matrixbase.MatrixBase.cross
        sympy.matrices.matrixbase.MatrixBase.dot
        multiply
        """
        if self.shape != other.shape:
            raise ShapeError("Matrix shapes must agree {} != {}".format(self.shape, other.shape))

        return self._eval_matrix_mul_elementwise(other)

    def __neg__(self) -> Self:
        return self._eval_scalar_mul(-1)

    @call_highest_priority('__rpow__')
    def __pow__(self, exp: SExpr) -> MatrixBase | MatrixExpr:
        """Return self**exp a scalar or symbol."""
        return self.pow(exp)

    @overload
    def pow(self, exp: int | Integer, method: str | None = None) -> Self: ...
    @overload
    def pow(self, exp: SExpr, method: str | None = None) -> MatrixBase | MatrixExpr: ...

    def pow(self, exp: SExpr, method: str | None = None) -> MatrixBase | MatrixExpr:
        r"""Return self**exp a scalar or symbol.

        Parameters
        ==========

        method : multiply, mulsimp, jordan, cayley
            If multiply then it returns exponentiation using recursion.
            If jordan then Jordan form exponentiation will be used.
            If cayley then the exponentiation is done using Cayley-Hamilton
            theorem.
            If mulsimp then the exponentiation is done using recursion
            with dotprodsimp. This specifies whether intermediate term
            algebraic simplification is used during naive matrix power to
            control expression blowup and thus speed up calculation.
            If None, then it heuristically decides which method to use.

        """

        if method is not None and method not in ['multiply', 'mulsimp', 'jordan', 'cayley']:
            raise TypeError('No such method')
        if self.rows != self.cols:
            raise NonSquareMatrixError()
        a = self
        jordan_pow = getattr(a, '_matrix_pow_by_jordan_blocks', None)
        exp = sympify(exp)

        if exp.is_zero:
            return a._new(a.rows, a.cols, lambda i, j: int(i == j))
        if exp == 1:
            return a

        diagonal = getattr(a, 'is_diagonal', None)
        if diagonal is not None and diagonal():
            return a._new(a.rows, a.cols, lambda i, j: a[i,j]**exp if i == j else 0)

        if exp.is_Number and exp % 1 == 0:
            if a.rows == 1:
                return a._new([[a[0]**exp]])
            if exp < 0:
                exp = -exp
                a = a.inv()
        # When certain conditions are met,
        # Jordan block algorithm is faster than
        # computation by recursion.
        if method == 'jordan' and jordan_pow is not None:
            return jordan_pow(exp)

        elif method == 'cayley':
            if not isinstance(exp, Integer):
                raise ValueError("cayley method is only valid for integer powers")
            return a._eval_pow_by_cayley(exp.p)

        elif method == "mulsimp":
            if not isinstance(exp, Integer):
                raise ValueError("mulsimp method is only valid for integer powers")
            return a._eval_pow_by_recursion_dotprodsimp(exp.p)

        elif method == "multiply":
            if not isinstance(exp, Integer):
                raise ValueError("multiply method is only valid for integer powers")
            return a._eval_pow_by_recursion(exp.p)

        elif method is None and isinstance(exp, (Integer, Float)) and exp % 1 == 0:
            if isinstance(exp, Float):
                expi = Integer(exp)
            else:
                expi = exp
            # Decide heuristically which method to apply
            if a.rows == 2 and expi > 100000 and jordan_pow is not None:
                return jordan_pow(expi)
            elif _get_intermediate_simp_bool(True, None):
                return a._eval_pow_by_recursion_dotprodsimp(expi.p)
            elif expi > 10000:
                return a._eval_pow_by_cayley(expi.p)
            else:
                return a._eval_pow_by_recursion(expi.p)

        if jordan_pow:
            try:
                return jordan_pow(exp)
            except NonInvertibleMatrixError:
                # Raised by jordan_pow on zero determinant matrix unless exp is
                # definitely known to be a non-negative integer.
                # Here we raise if n is definitely not a non-negative integer
                # but otherwise we can leave this as an unevaluated MatPow.
                if exp.is_integer is False or exp.is_nonnegative is False:
                    raise

        from sympy.matrices.expressions import MatPow
        return MatPow(a, exp)

    @call_highest_priority('__add__')
    def __radd__(self, other: MatrixBase) -> MatrixBase:
        return self.__add__(other)

    @call_highest_priority('__matmul__')
    def __rmatmul__(self, other: MatrixBase) -> MatrixBase:
        self, other, T = _unify_with_other(self, other)

        if T != "is_matrix":
            return NotImplemented

        return self.__rmul__(other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other: MatrixBase | SExpr) -> MatrixBase:
        return self.rmultiply(other)

    def rmultiply(self, other: MatrixBase | SExpr, dotprodsimp: bool | None = None) -> MatrixBase:
        """Same as __rmul__() but with optional simplification.

        Parameters
        ==========

        dotprodsimp : bool, optional
            Specifies whether intermediate term algebraic simplification is used
            during matrix multiplications to control expression blowup and thus
            speed up calculation. Default is off.
        """
        isimpbool = _get_intermediate_simp_bool(False, dotprodsimp)
        self, other_mat, T = _unify_with_other(self, other)

        if isinstance(other_mat, MatrixBase):

            if self.shape[0] != other_mat.shape[1]:
                raise ShapeError("Matrix size mismatch.")

            m = self._eval_matrix_rmul(other_mat)

            if isimpbool:
                return m._new(m.rows, m.cols, [_dotprodsimp(e) for e in m.flat()])

            return m

        elif T == "possible_scalar":
            try:
                return self._eval_scalar_rmul(other_mat)
            except TypeError:
                return NotImplemented

        else:
            return NotImplemented

    # XXX: No idea why mypy complains about this but type: ignore seems to be
    # needed here...
    @call_highest_priority('__sub__')
    def __rsub__(self, a: MatrixBase) -> MatrixBase: # type: ignore
        return (-self) + a

    @call_highest_priority('__rsub__')
    def __sub__(self: Tmat, a: Tmat) -> Tmat:
        return self + (-a)

    def _eval_det_bareiss(self, iszerofunc=_is_zero_after_expand_mul) -> Expr:
        return _det_bareiss(self, iszerofunc=iszerofunc)

    def _eval_det_berkowitz(self) -> Expr:
        return _det_berkowitz(self)

    def _eval_det_lu(self, iszerofunc=_iszero, simpfunc=None) -> Expr:
        return _det_LU(self, iszerofunc=iszerofunc, simpfunc=simpfunc)

    def _eval_det_bird(self) -> Expr:
        return _det_bird(self)

    def _eval_det_laplace(self) -> Expr:
        return _det_laplace(self)

    # for expressions.determinant.Determinant
    def _eval_determinant(self) -> Expr:
        return _det(self)

    def adjugate(self, method: str="berkowitz") -> Self:
        return _adjugate(self, method=method)

    def charpoly(self, x: str | Expr = 'lambda', simplify=_utilities_simplify) -> Poly:
        return _charpoly(self, x=x, simplify=simplify)

    def cofactor(self, i, j, method: str="berkowitz") -> Expr:
        return _cofactor(self, i, j, method=method)

    def cofactor_matrix(self, method: str="berkowitz") -> Self:
        return _cofactor_matrix(self, method=method)

    def det(self, method: str="bareiss", iszerofunc=None) -> Expr:
        return _det(self, method=method, iszerofunc=iszerofunc)

    def per(self) -> Expr:
        return _per(self)

    def minor(self, i, j, method: str="berkowitz") -> Expr:
        return _minor(self, i, j, method=method)

    def minor_submatrix(self, i, j) -> Self:
        return _minor_submatrix(self, i, j)

    _find_reasonable_pivot.__doc__       = _find_reasonable_pivot.__doc__
    _find_reasonable_pivot_naive.__doc__ = _find_reasonable_pivot_naive.__doc__
    _eval_det_bareiss.__doc__            = _det_bareiss.__doc__
    _eval_det_berkowitz.__doc__          = _det_berkowitz.__doc__
    _eval_det_bird.__doc__            = _det_bird.__doc__
    _eval_det_laplace.__doc__            = _det_laplace.__doc__
    _eval_det_lu.__doc__                 = _det_LU.__doc__
    _eval_determinant.__doc__            = _det.__doc__
    adjugate.__doc__                     = _adjugate.__doc__
    charpoly.__doc__                     = _charpoly.__doc__
    cofactor.__doc__                     = _cofactor.__doc__
    cofactor_matrix.__doc__              = _cofactor_matrix.__doc__
    det.__doc__                          = _det.__doc__
    per.__doc__                          = _per.__doc__
    minor.__doc__                        = _minor.__doc__
    minor_submatrix.__doc__              = _minor_submatrix.__doc__

    @overload
    def echelon_form(self,
                iszerofunc: Callable[[Expr], bool | None] = _iszero,
                simplify: bool | Callable[[Expr], Expr] = False,
                 *,
                with_pivots: Literal[False] = False,
            ) -> Self: ...
    @overload
    def echelon_form(self,
                iszerofunc: Callable[[Expr], bool | None] = _iszero,
                simplify: bool | Callable[[Expr], Expr] = False,
                *,
                with_pivots: Literal[True],
            ) -> tuple[Self, tuple[int]]: ...

    def echelon_form(self,
                iszerofunc: Callable[[Expr], bool | None] = _iszero,
                simplify: bool | Callable[[Expr], Expr] = False,
                *,
                with_pivots: bool = False,
            ) -> Self | tuple[Self, tuple[int]]:
        return _echelon_form(self, iszerofunc=iszerofunc, simplify=simplify,
                with_pivots=with_pivots)

    @property
    def is_echelon(self) -> bool:
        return _is_echelon(self)

    def rank(self,
             iszerofunc: Callable[[Expr], bool | None] = _iszero,
             simplify: bool | Callable[[Expr], Expr] = False) -> int:
        return _rank(self, iszerofunc=iszerofunc, simplify=simplify)

    def rref_rhs(self, rhs: Self) -> tuple[Self, Self]:
        """Return reduced row-echelon form of matrix, matrix showing
        rhs after reduction steps. ``rhs`` must have the same number
        of rows as ``self``.

        Examples
        ========

        >>> from sympy import Matrix, symbols
        >>> r1, r2 = symbols('r1 r2')
        >>> Matrix([[1, 1], [2, 1]]).rref_rhs(Matrix([r1, r2]))
        (Matrix([
        [1, 0],
        [0, 1]]), Matrix([
        [ -r1 + r2],
        [2*r1 - r2]]))
        """
        r, _ = _rref(self.hstack(self, self.eye(self.rows), rhs))
        return r[:, :self.cols], r[:, -rhs.cols:]

    @overload
    def rref(self,
             iszerofunc: Callable[[Expr], bool | None] = _iszero,
             simplify: bool | Callable[[Expr], Expr] = False,
             *,
             pivots: Literal[True] = True,
             normalize_last: bool = True
        ) -> tuple[Self, tuple[int]]: ...
    @overload
    def rref(self,
             iszerofunc: Callable[[Expr], bool | None] = _iszero,
             simplify: bool | Callable[[Expr], Expr] = False,
             *,
             pivots: Literal[False],
             normalize_last: bool = True
        ) -> Self: ...

    def rref(self,
             iszerofunc: Callable[[Expr], bool | None] = _iszero,
             simplify: bool | Callable[[Expr], Expr] = False,
             *,
             pivots: bool = True,
             normalize_last: bool = True
         ) -> Self | tuple[Self, tuple[int]]:
        return _rref(self, iszerofunc=iszerofunc, simplify=simplify,
            pivots=pivots, normalize_last=normalize_last)

    echelon_form.__doc__ = _echelon_form.__doc__
    is_echelon.__doc__   = _is_echelon.__doc__
    rank.__doc__         = _rank.__doc__
    rref.__doc__         = _rref.__doc__

    def _normalize_op_args(self, op: str,
                           col: int | None,
                           k: SExpr | None,
                           col1: int | None,
                           col2: int | None,
                           error_str: str = "col"
                        ) -> tuple[str, int, int, int, int]:
        """Validate the arguments for a row/column operation.  ``error_str``
        can be one of "row" or "col" depending on the arguments being parsed."""
        if op not in ["n->kn", "n<->m", "n->n+km"]:
            raise ValueError("Unknown {} operation '{}'. Valid col operations "
                             "are 'n->kn', 'n<->m', 'n->n+km'".format(error_str, op))

        # define self_col according to error_str
        self_cols = self.cols if error_str == 'col' else self.rows

        # normalize and validate the arguments
        if op == "n->kn":
            col = col if col is not None else col1
            if col is None or k is None:
                raise ValueError("For a {0} operation 'n->kn' you must provide the "
                                 "kwargs `{0}` and `k`".format(error_str))
            if not 0 <= col < self_cols:
                raise ValueError("This matrix does not have a {} '{}'".format(error_str, col))

        elif op == "n<->m":
            # we need two cols to swap. It does not matter
            # how they were specified, so gather them together and
            # remove `None`
            cols = {col, k, col1, col2}.difference([None])
            if len(cols) > 2:
                # maybe the user left `k` by mistake?
                cols = {col, col1, col2}.difference([None]) # type: ignore
            if len(cols) != 2:
                raise ValueError("For a {0} operation 'n<->m' you must provide the "
                                 "kwargs `{0}1` and `{0}2`".format(error_str))
            col1, col2 = cols # type: ignore
            if not 0 <= col1 < self_cols: # type: ignore
                raise ValueError("This matrix does not have a {} '{}'".format(error_str, col1))
            if not 0 <= col2 < self_cols: # type: ignore
                raise ValueError("This matrix does not have a {} '{}'".format(error_str, col2))

        elif op == "n->n+km":
            col = col1 if col is None else col
            col2 = col1 if col2 is None else col2
            if col is None or col2 is None or k is None:
                raise ValueError("For a {0} operation 'n->n+km' you must provide the "
                                 "kwargs `{0}`, `k`, and `{0}2`".format(error_str))
            if col == col2:
                raise ValueError("For a {0} operation 'n->n+km' `{0}` and `{0}2` must "
                                 "be different.".format(error_str))
            if not 0 <= col < self_cols:
                raise ValueError("This matrix does not have a {} '{}'".format(error_str, col))
            if not 0 <= col2 < self_cols:
                raise ValueError("This matrix does not have a {} '{}'".format(error_str, col2))

        else:
            raise ValueError('invalid operation %s' % repr(op))

        return op, col, k, col1, col2 # type: ignore

    def _eval_col_op_multiply_col_by_const(self, col: int, k: SExpr) -> Self:
        def entry(i, j):
            if j == col:
                return k * self[i, j]
            return self[i, j]
        return self._new(self.rows, self.cols, entry)

    def _eval_col_op_swap(self, col1: int, col2: int) -> Self:
        def entry(i, j):
            if j == col1:
                return self[i, col2]
            elif j == col2:
                return self[i, col1]
            return self[i, j]
        return self._new(self.rows, self.cols, entry)

    def _eval_col_op_add_multiple_to_other_col(self, col: int, k: SExpr, col2: int) -> Self:
        def entry(i, j):
            if j == col:
                return self[i, j] + k * self[i, col2]
            return self[i, j]
        return self._new(self.rows, self.cols, entry)

    def _eval_row_op_swap(self, row1: int, row2: int) -> Self:
        def entry(i, j):
            if i == row1:
                return self[row2, j]
            elif i == row2:
                return self[row1, j]
            return self[i, j]
        return self._new(self.rows, self.cols, entry)

    def _eval_row_op_multiply_row_by_const(self, row: int, k: SExpr) -> Self:
        def entry(i, j):
            if i == row:
                return k * self[i, j]
            return self[i, j]
        return self._new(self.rows, self.cols, entry)

    def _eval_row_op_add_multiple_to_other_row(self, row: int, k: SExpr, row2: int) -> Self:
        def entry(i, j):
            if i == row:
                return self[i, j] + k * self[row2, j]
            return self[i, j]
        return self._new(self.rows, self.cols, entry)

    def elementary_col_op(self, op: str = "n->kn",
                          col: int | None = None,
                          k: SExpr | None = None,
                          col1: int | None = None,
                          col2: int | None = None) -> Self:
        """Performs the elementary column operation `op`.

        `op` may be one of

            * ``"n->kn"`` (column n goes to k*n)
            * ``"n<->m"`` (swap column n and column m)
            * ``"n->n+km"`` (column n goes to column n + k*column m)

        Parameters
        ==========

        op : string; the elementary row operation
        col : the column to apply the column operation
        k : the multiple to apply in the column operation
        col1 : one column of a column swap
        col2 : second column of a column swap or column "m" in the column operation
               "n->n+km"
        """

        op, col, k, col1, col2 = self._normalize_op_args(op, col, k, col1, col2, "col")

        # now that we've validated, we're all good to dispatch
        if op == "n->kn":
            return self._eval_col_op_multiply_col_by_const(col, k)
        elif op == "n<->m":
            return self._eval_col_op_swap(col1, col2)
        elif op == "n->n+km":
            return self._eval_col_op_add_multiple_to_other_col(col, k, col2)
        else:
            raise ValueError(f'invalid operation {op!r}')

    def elementary_row_op(self, op: str = "n->kn",
                          row: int | None = None,
                          k: SExpr | None = None,
                          row1: int | None = None,
                          row2: int | None = None) -> Self:
        """Performs the elementary row operation `op`.

        `op` may be one of

            * ``"n->kn"`` (row n goes to k*n)
            * ``"n<->m"`` (swap row n and row m)
            * ``"n->n+km"`` (row n goes to row n + k*row m)

        Parameters
        ==========

        op : string; the elementary row operation
        row : the row to apply the row operation
        k : the multiple to apply in the row operation
        row1 : one row of a row swap
        row2 : second row of a row swap or row "m" in the row operation
               "n->n+km"
        """

        op, row, k, row1, row2 = self._normalize_op_args(op, row, k, row1, row2, "row")

        # now that we've validated, we're all good to dispatch
        if op == "n->kn":
            return self._eval_row_op_multiply_row_by_const(row, k)
        elif op == "n<->m":
            return self._eval_row_op_swap(row1, row2)
        elif op == "n->n+km":
            return self._eval_row_op_add_multiple_to_other_row(row, k, row2)
        else:
            raise ValueError(f'invalid operation {op!r}')

    def columnspace(self, simplify: bool | Callable[[Expr], Expr] = False) -> list[Self]:
        return _columnspace(self, simplify=simplify)

    def nullspace(self, simplify: bool | Callable[[Expr], Expr] = False, iszerofunc=_iszero) -> list[Self]:
        return _nullspace(self, simplify=simplify, iszerofunc=iszerofunc)

    def rowspace(self, simplify: bool | Callable[[Expr], Expr] = False) -> list[Self]:
        return _rowspace(self, simplify=simplify)

    # XXX: Somehow replacing this with an ordinary use of classmethod breaks
    # sphinx ...
    def orthogonalize(cls, *vecs: Self, **kwargs) -> list[Self]:
        return _orthogonalize(cls, *vecs, **kwargs)

    columnspace.__doc__   = _columnspace.__doc__
    nullspace.__doc__     = _nullspace.__doc__
    rowspace.__doc__      = _rowspace.__doc__
    orthogonalize.__doc__ = _orthogonalize.__doc__

    orthogonalize = classmethod(orthogonalize) # type: ignore

    def eigenvals(self,
                  error_when_incomplete: bool = True,
                  **flags) -> dict[Expr, int]:
        return _eigenvals(self, error_when_incomplete=error_when_incomplete, **flags)

    def eigenvects(self,
                   error_when_incomplete: bool = True,
                   iszerofunc: Callable[[Expr], bool | None] = _iszero,
                   **flags) -> list[tuple[Expr, int, list[Self]]]:
        return _eigenvects(self, error_when_incomplete=error_when_incomplete,
                iszerofunc=iszerofunc, **flags)

    def is_diagonalizable(self,
                          reals_only: bool = False,
                          **kwargs) -> bool:
        return _is_diagonalizable(self, reals_only=reals_only, **kwargs)

    def diagonalize(self,
                    reals_only: bool = False,
                    sort: bool = False,
                    normalize: bool = False
                    ) -> tuple[Self, Self]:
        return _diagonalize(self, reals_only=reals_only, sort=sort,
                normalize=normalize)

    def bidiagonalize(self, upper: bool=True) -> Self:
        return _bidiagonalize(self, upper=upper)

    def bidiagonal_decomposition(self, upper: bool=True) -> tuple[Self, Self, Self]:
        return _bidiagonal_decomposition(self, upper=upper)

    @property
    def is_positive_definite(self) -> bool | None:
        return _is_positive_definite(self)

    @property
    def is_positive_semidefinite(self) -> bool | None:
        return _is_positive_semidefinite(self)

    @property
    def is_negative_definite(self) -> bool | None:
        return _is_negative_definite(self)

    @property
    def is_negative_semidefinite(self) -> bool | None:
        return _is_negative_semidefinite(self)

    @property
    def is_indefinite(self) -> bool | None:
        return _is_indefinite(self)

    @overload
    def jordan_form(
            self: Tmat,
            calc_transform: Literal[True] = True,
            *,
            chop: bool = False
        ) -> tuple[Tmat, Tmat]: ...
    @overload
    def jordan_form(
            self: Tmat,
            calc_transform: Literal[False],
            *,
            chop: bool = False,
        ) -> Tmat: ...

    def jordan_form(self: Tmat,
             calc_transform: bool = True,
             *,
             chop: bool = False
         ) -> tuple[Tmat, Tmat] | Tmat:
        return _jordan_form(self, calc_transform=calc_transform, chop=chop)

    def left_eigenvects(self, **flags: Any) -> list[tuple[Expr, int, list[Self]]]:
        return _left_eigenvects(self, **flags)

    def singular_values(self) -> list[Expr]:
        return _singular_values(self)

    eigenvals.__doc__                  = _eigenvals.__doc__
    eigenvects.__doc__                 = _eigenvects.__doc__
    is_diagonalizable.__doc__          = _is_diagonalizable.__doc__
    diagonalize.__doc__                = _diagonalize.__doc__
    is_positive_definite.__doc__       = _is_positive_definite.__doc__
    is_positive_semidefinite.__doc__   = _is_positive_semidefinite.__doc__
    is_negative_definite.__doc__       = _is_negative_definite.__doc__
    is_negative_semidefinite.__doc__   = _is_negative_semidefinite.__doc__
    is_indefinite.__doc__              = _is_indefinite.__doc__
    jordan_form.__doc__                = _jordan_form.__doc__
    left_eigenvects.__doc__            = _left_eigenvects.__doc__
    singular_values.__doc__            = _singular_values.__doc__
    bidiagonalize.__doc__              = _bidiagonalize.__doc__
    bidiagonal_decomposition.__doc__   = _bidiagonal_decomposition.__doc__

    def diff(self, *args: Expr | int | tuple[Expr, int], evaluate: bool = True) -> Self:
        """Calculate the derivative of each element in the matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.abc import x, y
        >>> M = Matrix([[x, y], [1, 0]])
        >>> M.diff(x)
        Matrix([
        [1, 0],
        [0, 0]])

        See Also
        ========

        integrate
        limit
        """
        # XXX this should be handled here rather than in Derivative
        from sympy.tensor.array.array_derivatives import ArrayDerivative
        deriv = ArrayDerivative(self, *args, evaluate=evaluate)
        # XXX This can rather changed to always return immutable matrix
        if not isinstance(self, Basic) and evaluate:
            return deriv.as_mutable() # type: ignore
        return deriv # type: ignore

    def _eval_derivative(self, arg: Expr) -> Self:
        return self.applyfunc(lambda x: x.diff(arg))

    def integrate(self, *args: SymbolLimits, **kwargs):
        """Integrate each element of the matrix.  ``args`` will
        be passed to the ``integrate`` function.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.abc import x, y
        >>> M = Matrix([[x, y], [1, 0]])
        >>> M.integrate((x, ))
        Matrix([
        [x**2/2, x*y],
        [     x,   0]])
        >>> M.integrate((x, 0, 2))
        Matrix([
        [2, 2*y],
        [2,   0]])

        See Also
        ========

        limit
        diff
        """
        return self.applyfunc(lambda x: x.integrate(*args, **kwargs))

    def jacobian(self, X: MatrixBase | list[Expr]) -> Self:
        """Calculates the Jacobian matrix (derivative of a vector-valued function).

        Parameters
        ==========

        ``self`` : vector of expressions representing functions f_i(x_1, ..., x_n).
        X : set of x_i's in order, it can be a list or a Matrix

        Both ``self`` and X can be a row or a column matrix in any order
        (i.e., jacobian() should always work).

        Examples
        ========

        >>> from sympy import sin, cos, Matrix
        >>> from sympy.abc import rho, phi
        >>> X = Matrix([rho*cos(phi), rho*sin(phi), rho**2])
        >>> Y = Matrix([rho, phi])
        >>> X.jacobian(Y)
        Matrix([
        [cos(phi), -rho*sin(phi)],
        [sin(phi),  rho*cos(phi)],
        [   2*rho,             0]])
        >>> X = Matrix([rho*cos(phi), rho*sin(phi)])
        >>> X.jacobian(Y)
        Matrix([
        [cos(phi), -rho*sin(phi)],
        [sin(phi),  rho*cos(phi)]])

        See Also
        ========

        hessian
        wronskian
        """
        if not isinstance(X, MatrixBase):
            X = self._new(X)
        # Both X and ``self`` can be a row or a column matrix, so we need to make
        # sure all valid combinations work, but everything else fails:
        if self.shape[0] == 1:
            m = self.shape[1]
        elif self.shape[1] == 1:
            m = self.shape[0]
        else:
            raise TypeError("``self`` must be a row or a column matrix")
        if X.shape[0] == 1:
            n = X.shape[1]
        elif X.shape[1] == 1:
            n = X.shape[0]
        else:
            raise TypeError("X must be a row or a column matrix")

        # m is the number of functions and n is the number of variables
        # computing the Jacobian is now easy:
        return self._new(m, n, lambda j, i: self[j].diff(X[i]))

    def limit(self, x: Expr, xlim: Expr, dir: str = '+') -> Self:
        """Calculate the limit of each element in the matrix.
        ``args`` will be passed to the ``limit`` function.

        Examples
        ========

        >>> from sympy import Matrix
        >>> from sympy.abc import x, y
        >>> M = Matrix([[x, y], [1, 0]])
        >>> M.limit(x, 2)
        Matrix([
        [2, y],
        [1, 0]])

        See Also
        ========

        integrate
        diff
        """
        return self.applyfunc(lambda e: e.limit(x, xlim, dir))

    def berkowitz_charpoly(self, x: str | Expr = Dummy('lambda'),
               simplify: Callable[[Expr], Expr] = _utilities_simplify) -> Poly:
        return self.charpoly(x=x)

    def berkowitz_det(self) -> Expr:
        """Computes determinant using Berkowitz method.

        See Also
        ========

        det
        """
        return self.det(method='berkowitz')

    def berkowitz_eigenvals(self, **flags: Any) -> dict[Expr, int]:
        """Computes eigenvalues of a Matrix using Berkowitz method."""
        return self.eigenvals(**flags)

    def berkowitz_minors(self) -> tuple[Expr]:
        """Computes principal minors using Berkowitz method."""
        sign, minors = self.one, []

        for poly in self.berkowitz():
            minors.append(sign * poly[-1])
            sign = -sign

        return tuple(minors)

    def berkowitz(self) -> tuple[tuple[Expr], ...]:
        from sympy.matrices import zeros
        berk = ((S.One,),)
        if not self:
            return berk

        if not self.is_square:
            raise NonSquareMatrixError()

        A, N = self, self.rows
        transforms: list[int | MatrixBase] = [0] * (N - 1)

        for n in range(N, 1, -1):
            T, k = zeros(n + 1, n), n - 1

            R, C = -A[k, :k], A[:k, k]
            A, a = A[:k, :k], -A[k, k]

            items = [C]

            for i in range(0, n - 2):
                items.append(A * items[i])

            items2 = [self.one, a] + [(R * B)[0, 0] for B in items]

            for i in range(n):
                T[i:, i] = items2[:n - i + 1]

            transforms[k - 1] = T

        polys: list[MatrixBase] = [self._new([self.one, -A[0, 0]])]

        for i, T in enumerate(transforms):
            polys.append(T * polys[i])

        return berk + tuple((p[0,0],) for p in polys)

    def cofactorMatrix(self, method: str = "berkowitz") -> Self:
        return self.cofactor_matrix(method=method)

    def det_bareis(self) -> Expr:
        return _det_bareiss(self)

    def det_LU_decomposition(self) -> Expr:
        """Compute matrix determinant using LU decomposition.


        Note that this method fails if the LU decomposition itself
        fails. In particular, if the matrix has no inverse this method
        will fail.

        TODO: Implement algorithm for sparse matrices (SFF),
        http://www.eecis.udel.edu/~saunders/papers/sffge/it5.ps.

        See Also
        ========


        det
        berkowitz_det
        """
        return self.det(method='lu')

    def jordan_cell(self, eigenval: SExpr, n: int) -> Self:
        return self.jordan_block(n, eigenval)

    def jordan_cells(self, calc_transformation: bool = True) -> tuple[Self, list[Self]]:
        P, J = self.jordan_form()
        return P, J.get_diag_blocks()

    def minorEntry(self, i: int, j: int, method: str = "berkowitz") -> Expr:
        return self.minor(i, j, method=method)

    def minorMatrix(self, i: int, j: int) -> Self:
        return self.minor_submatrix(i, j)

    def permuteBkwd(self, perm: list[int] | list[list[int]] | Permutation) -> Self:
        """Permute the rows of the matrix with the given permutation in reverse."""
        return self.permute_rows(perm, direction='backward')

    def permuteFwd(self, perm: list[int] | list[list[int]] | Permutation) -> Self:
        """Permute the rows of the matrix with the given permutation."""
        return self.permute_rows(perm, direction='forward')

    @property
    def kind(self) -> MatrixKind:
        elem_kinds = {e.kind for e in self.flat()}
        if len(elem_kinds) == 1:
            elemkind, = elem_kinds
        else:
            elemkind = UndefinedKind
        return MatrixKind(elemkind)

    def flat(self) -> list[Expr]:
        """
        Returns a flat list of all elements in the matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> m = Matrix([[0, 2], [3, 4]])
        >>> m.flat()
        [0, 2, 3, 4]

        See Also
        ========

        tolist
        values
        """
        return [self[i, j] for i in range(self.rows) for j in range(self.cols)]

    def __array__(self, dtype=object, copy=None):
        if copy is not None and not copy:
            raise TypeError("Cannot implement copy=False when converting Matrix to ndarray")
        from .dense import matrix2numpy
        return matrix2numpy(self, dtype=dtype)

    def __len__(self) -> int:
        """Return the number of elements of ``self``.

        Implemented mainly so bool(Matrix()) == False.
        """
        return self.rows * self.cols

    def _matrix_pow_by_jordan_blocks(self, num) -> Self:
        from sympy.matrices import diag, MutableMatrix

        def jordan_cell_power(jc, n):
            N = jc.shape[0]
            l = jc[0,0]
            if l.is_zero:
                if N == 1 and n.is_nonnegative:
                    jc[0,0] = l**n
                elif not (n.is_integer and n.is_nonnegative):
                    raise NonInvertibleMatrixError("Non-invertible matrix can only be raised to a nonnegative integer")
                else:
                    for i in range(N):
                        jc[0,i] = KroneckerDelta(i, n)
            else:
                for i in range(N):
                    bn = binomial(n, i)
                    if isinstance(bn, binomial):
                        bn = bn._eval_expand_func()
                    jc[0,i] = l**(n-i)*bn
            for i in range(N):
                for j in range(1, N-i):
                    jc[j,i+j] = jc [j-1,i+j-1]

        P, J = self.jordan_form()
        jordan_cells = J.get_diag_blocks()
        # Make sure jordan_cells matrices are mutable:
        jordan_cells2 = [MutableMatrix(j) for j in jordan_cells]
        for j in jordan_cells2:
            jordan_cell_power(j, num)

        return self._as_type(P.multiply(diag(*jordan_cells2)).multiply(P.inv()))

    def __str__(self):
        if S.Zero in self.shape:
            return 'Matrix(%s, %s, [])' % (self.rows, self.cols)
        return "Matrix(%s)" % str(self.tolist())

    def _format_str(self, printer=None) -> str:
        if not printer:
            printer = StrPrinter()
        # Handle zero dimensions:
        if S.Zero in self.shape:
            return 'Matrix(%s, %s, [])' % (self.rows, self.cols)
        if self.rows == 1:
            return "Matrix([%s])" % self.table(printer, rowsep=',\n')
        return "Matrix([\n%s])" % self.table(printer, rowsep=',\n')

    @classmethod
    def irregular(cls, ntop: int, *matrices: Self) -> Self:
        """Return a matrix filled by the given matrices which
        are listed in order of appearance from left to right, top to
        bottom as they first appear in the matrix. They must fill the
        matrix completely.

        Examples
        ========

        >>> from sympy import ones, Matrix
        >>> Matrix.irregular(3, ones(2,1), ones(3,3)*2, ones(2,2)*3,
        ...   ones(1,1)*4, ones(2,2)*5, ones(1,2)*6, ones(1,2)*7)
        Matrix([
            [1, 2, 2, 2, 3, 3],
            [1, 2, 2, 2, 3, 3],
            [4, 2, 2, 2, 5, 5],
            [6, 6, 7, 7, 5, 5]])
        """
        ntop = as_int(ntop)
        # make sure we are working with explicit matrices
        b = [i.as_explicit() if hasattr(i, 'as_explicit') else i
            for i in matrices]
        q = list(range(len(b)))
        dat = [i.rows for i in b]
        active = [q.pop(0) for _ in range(ntop)]
        cols = sum(b[i].cols for i in active)
        rows = []
        while any(dat):
            r = []
            for a, j in enumerate(active):
                r.extend(b[j][-dat[j], :].flat())
                dat[j] -= 1
                if dat[j] == 0 and q:
                    active[a] = q.pop(0)
            if len(r) != cols:
                raise ValueError(filldedent('''
                    Matrices provided do not appear to fill
                    the space completely.'''))
            rows.append(r)
        return cls._new(rows)

    @classmethod
    def _handle_ndarray(cls, arg: Any) -> tuple[int, int, list[Expr]]:
        # NumPy array or matrix or some other object that implements
        # __array__. So let's first use this method to get a
        # numpy.array() and then make a Python list out of it.
        arr = arg.__array__()
        if len(arr.shape) == 2:
            rows, cols = arr.shape[0], arr.shape[1]
            flat_list = [cls._sympify(i) for i in arr.ravel()]
            return rows, cols, flat_list
        elif len(arr.shape) == 1:
            flat_list = [cls._sympify(i) for i in arr]
            return arr.shape[0], 1, flat_list
        else:
            raise NotImplementedError(
                "SymPy supports just 1D and 2D matrices")

    @classmethod
    def _handle_creation_inputs(cls, *args, **kwargs):
        """Return the number of rows, cols and flat matrix elements.

        Examples
        ========

        >>> from sympy import Matrix, I

        Matrix can be constructed as follows:

        * from a nested list of iterables

        >>> Matrix( ((1, 2+I), (3, 4)) )
        Matrix([
        [1, 2 + I],
        [3,     4]])

        * from un-nested iterable (interpreted as a column)

        >>> Matrix( [1, 2] )
        Matrix([
        [1],
        [2]])

        * from un-nested iterable with dimensions

        >>> Matrix(1, 2, [1, 2] )
        Matrix([[1, 2]])

        * from no arguments (a 0 x 0 matrix)

        >>> Matrix()
        Matrix(0, 0, [])

        * from a rule

        >>> Matrix(2, 2, lambda i, j: i/(j + 1) )
        Matrix([
        [0,   0],
        [1, 1/2]])

        See Also
        ========
        irregular - filling a matrix with irregular blocks
        """
        from sympy.matrices import SparseMatrix
        from sympy.matrices.expressions.matexpr import MatrixSymbol
        from sympy.matrices.expressions.blockmatrix import BlockMatrix

        flat_list = None

        if len(args) == 1:

            [arg1] = args

            # Matrix(SparseMatrix(...))
            if isinstance(arg1, SparseMatrix):
                return arg1.rows, arg1.cols, flatten(arg1.tolist())

            # Matrix(Matrix(...))
            elif isinstance(arg1, MatrixBase):
                return arg1.rows, arg1.cols, arg1.flat()

            # Matrix(MatrixSymbol('X', 2, 2))
            elif isinstance(arg1, Basic) and arg1.is_Matrix:
                return arg1.rows, arg1.cols, arg1.as_explicit().flat() # type: ignore

            elif isinstance(args[0], mp.matrix):
                M = args[0]
                flat_list = [cls._sympify(x) for x in M]
                return M.rows, M.cols, flat_list

            # Matrix(numpy.ones((2, 2)))
            elif hasattr(args[0], "__array__"):
                return cls._handle_ndarray(args[0])

            # Matrix([1, 2, 3]) or Matrix([[1, 2], [3, 4]])
            elif is_sequence(args[0]) and not isinstance(args[0], DeferredVector):
                dat = list(args[0])
                ismat = lambda i: isinstance(i, MatrixBase) and (
                    evaluate or isinstance(i, (BlockMatrix, MatrixSymbol)))
                raw = lambda i: is_sequence(i) and not ismat(i)
                evaluate = kwargs.get('evaluate', True)

                if evaluate:

                    def make_explicit(x):
                        """make Block and Symbol explicit"""
                        if isinstance(x, BlockMatrix):
                            return x.as_explicit()
                        elif isinstance(x, MatrixSymbol) and all(_.is_Integer for _ in x.shape):
                            return x.as_explicit()
                        else:
                            return x

                    def make_explicit_row(row):
                        # Could be list or could be list of lists
                        if isinstance(row, (list, tuple)):
                            return [make_explicit(x) for x in row]
                        else:
                            return make_explicit(row)

                    if isinstance(dat, (list, tuple)):
                        dat = [make_explicit_row(row) for row in dat]

                if len(dat) == 0:
                    rows = cols = 0
                    flat_list = []
                elif all(raw(i) for i in dat) and len(dat[0]) == 0: # type: ignore
                    if not all(len(i) == 0 for i in dat): # type: ignore
                        raise ValueError('mismatched dimensions')
                    rows = len(dat)
                    cols = 0
                    flat_list = []
                elif not any(raw(i) or ismat(i) for i in dat):
                    # a column as a list of values
                    flat_list = [cls._sympify(i) for i in dat] # type: ignore
                    rows = len(flat_list)
                    cols = 1 if rows else 0
                elif evaluate and all(ismat(i) for i in dat):
                    # a column as a list of matrices
                    ncol = {i.cols for i in dat if any(i.shape)} # type: ignore
                    if ncol:
                        if len(ncol) != 1:
                            raise ValueError('mismatched dimensions')
                        flat_list = [_ for i in dat for r in i.tolist() for _ in r] # type: ignore
                        cols = ncol.pop()
                        rows = len(flat_list)//cols
                    else:
                        rows = cols = 0
                        flat_list = []
                elif evaluate and any(ismat(i) for i in dat):
                    ncol = set()
                    flat_list = []
                    for i in dat:
                        if ismat(i):
                            flat_list.extend([k for j in i.tolist() for k in j]) # type: ignore
                            if any(i.shape): # type: ignore
                                ncol.add(i.cols) # type: ignore
                        elif raw(i):
                            if i:
                                ncol.add(len(i)) # type: ignore
                                flat_list.extend([cls._sympify(ij) for ij in i]) # type: ignore
                        else:
                            ncol.add(1)
                            flat_list.append(i)
                        if len(ncol) > 1:
                            raise ValueError('mismatched dimensions')
                    cols = ncol.pop()
                    rows = len(flat_list)//cols
                else:
                    # list of lists; each sublist is a logical row
                    # which might consist of many rows if the values in
                    # the row are matrices
                    flat_list = []
                    ncol = set()
                    rows = cols = 0
                    for row in dat:
                        if not is_sequence(row) and \
                                not getattr(row, 'is_Matrix', False):
                            raise ValueError('expecting list of lists')

                        if hasattr(row, '__array__'):
                            if 0 in row.shape: # type: ignore
                                continue

                        if evaluate and all(ismat(i) for i in row):
                            r, c, flatT = cls._handle_creation_inputs([i.T for i in row]) # type: ignore
                            T = reshape(flatT, [c])
                            flat = \
                                [T[i][j] for j in range(c) for i in range(r)]
                            r, c = c, r
                        else:
                            r = 1
                            if getattr(row, 'is_Matrix', False):
                                c = 1
                                flat = [row]
                            else:
                                c = len(row) # type: ignore
                                flat = [cls._sympify(i) for i in row] # type: ignore
                        ncol.add(c)
                        if len(ncol) > 1:
                            raise ValueError('mismatched dimensions')
                        flat_list.extend(flat)
                        rows += r
                    cols = ncol.pop() if ncol else 0

        elif len(args) == 3:
            rows = as_int(args[0])
            cols = as_int(args[1])

            if rows < 0 or cols < 0:
                raise ValueError("Cannot create a {} x {} matrix. "
                                 "Both dimensions must be positive".format(rows, cols))

            # Matrix(2, 2, lambda i, j: i+j)
            if len(args) == 3 and isinstance(args[2], Callable):
                op = args[2]
                flat_list = []
                for i in range(rows):
                    flat_list.extend(
                        [cls._sympify(op(cls._sympify(i), cls._sympify(j)))
                         for j in range(cols)])

            # Matrix(2, 2, [1, 2, 3, 4])
            elif len(args) == 3 and is_sequence(args[2]):
                flat_list = args[2]
                if len(flat_list) != rows * cols:
                    raise ValueError(
                        'List length should be equal to rows*columns')
                flat_list = [cls._sympify(i) for i in flat_list]


        # Matrix()
        elif len(args) == 0:
            # Empty Matrix
            rows = cols = 0
            flat_list = []

        else:
            raise TypeError("Need 0, 1 or 3 arguments")

        if flat_list is None:
            raise TypeError(filldedent('''
                Data type not understood; expecting list of lists
                or lists of values.'''))

        return rows, cols, flat_list # type: ignore

    def add(self, b: Self) -> Self:
        """Return self + b."""
        return self + b

    def condition_number(self) -> Expr:
        """Returns the condition number of a matrix.

        This is the maximum singular value divided by the minimum singular value

        Examples
        ========

        >>> from sympy import Matrix, S
        >>> A = Matrix([[1, 0, 0], [0, 10, 0], [0, 0, S.One/10]])
        >>> A.condition_number()
        100

        See Also
        ========

        singular_values
        """

        if not self:
            return self.zero
        singularvalues = self.singular_values()
        return Max(*singularvalues) / Min(*singularvalues)

    def copy(self) -> Self:
        """
        Returns the copy of a matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> A = Matrix(2, 2, [1, 2, 3, 4])
        >>> A.copy()
        Matrix([
        [1, 2],
        [3, 4]])

        """
        return self._new(self.rows, self.cols, self.flat())

    def cross(self, b: MatrixExpr) -> Self:
        r"""
        Return the cross product of ``self`` and ``b`` relaxing the condition
        of compatible dimensions: if each has 3 elements, a matrix of the
        same type and shape as ``self`` will be returned. If ``b`` has the same
        shape as ``self`` then common identities for the cross product (like
        `a \times b = - b \times a`) will hold.

        Parameters
        ==========
            b : 3x1 or 1x3 Matrix

        See Also
        ========

        dot
        hat
        vee
        multiply
        multiply_elementwise
        """
        from sympy.matrices.expressions.matexpr import MatrixExpr

        if not isinstance(b, (MatrixBase, MatrixExpr)):
            raise TypeError(
                "{} must be a Matrix, not {}.".format(b, type(b)))

        if not (self.rows * self.cols == b.rows * b.cols == 3):
            raise ShapeError("Dimensions incorrect for cross product: %s x %s" %
                             ((self.rows, self.cols), (b.rows, b.cols)))
        else:
            return self._new(self.rows, self.cols, (
                (self[1] * b[2] - self[2] * b[1]),
                (self[2] * b[0] - self[0] * b[2]),
                (self[0] * b[1] - self[1] * b[0])))

    def hat(self) -> Self:
        r"""
        Return the skew-symmetric matrix representing the cross product,
        so that ``self.hat() * b`` is equivalent to  ``self.cross(b)``.

        Examples
        ========

        Calling ``hat`` creates a skew-symmetric 3x3 Matrix from a 3x1 Matrix:

        >>> from sympy import Matrix
        >>> a = Matrix([1, 2, 3])
        >>> a.hat()
        Matrix([
        [ 0, -3,  2],
        [ 3,  0, -1],
        [-2,  1,  0]])

        Multiplying it with another 3x1 Matrix calculates the cross product:

        >>> b = Matrix([3, 2, 1])
        >>> a.hat() * b
        Matrix([
        [-4],
        [ 8],
        [-4]])

        Which is equivalent to calling the ``cross`` method:

        >>> a.cross(b)
        Matrix([
        [-4],
        [ 8],
        [-4]])

        See Also
        ========

        dot
        cross
        vee
        multiply
        multiply_elementwise
        """

        if self.shape != (3, 1):
            raise ShapeError("Dimensions incorrect, expected (3, 1), got " +
                             str(self.shape))
        else:
            x, y, z = self.flat()
            return self._new(3, 3, (
                 0, -z,  y,
                 z,  0, -x,
                -y,  x,  0))

    def vee(self) -> Self:
        r"""
        Return a 3x1 vector from a skew-symmetric matrix representing the cross product,
        so that ``self * b`` is equivalent to  ``self.vee().cross(b)``.

        Examples
        ========

        Calling ``vee`` creates a vector from a skew-symmetric Matrix:

        >>> from sympy import Matrix
        >>> A = Matrix([[0, -3, 2], [3, 0, -1], [-2, 1, 0]])
        >>> a = A.vee()
        >>> a
        Matrix([
        [1],
        [2],
        [3]])

        Calculating the matrix product of the original matrix with a vector
        is equivalent to a cross product:

        >>> b = Matrix([3, 2, 1])
        >>> A * b
        Matrix([
        [-4],
        [ 8],
        [-4]])

        >>> a.cross(b)
        Matrix([
        [-4],
        [ 8],
        [-4]])

        ``vee`` can also be used to retrieve angular velocity expressions.
        Defining a rotation matrix:

        >>> from sympy import rot_ccw_axis3, trigsimp
        >>> from sympy.physics.mechanics import dynamicsymbols
        >>> theta = dynamicsymbols('theta')
        >>> R = rot_ccw_axis3(theta)
        >>> R
        Matrix([
        [cos(theta(t)), -sin(theta(t)), 0],
        [sin(theta(t)),  cos(theta(t)), 0],
        [            0,              0, 1]])

        We can retrieve the angular velocity:

        >>> Omega = R.T * R.diff()
        >>> Omega = trigsimp(Omega)
        >>> Omega.vee()
        Matrix([
        [                      0],
        [                      0],
        [Derivative(theta(t), t)]])

        See Also
        ========

        dot
        cross
        hat
        multiply
        multiply_elementwise
        """

        if self.shape != (3, 3):
            raise ShapeError("Dimensions incorrect, expected (3, 3), got " +
                             str(self.shape))
        elif not self.is_anti_symmetric():
            raise ValueError("Matrix is not skew-symmetric")
        else:
            return self._new(3, 1, (
                 self[2, 1],
                 self[0, 2],
                 self[1, 0]))

    @property
    def D(self) -> Self:
        """Return Dirac conjugate (if ``self.rows == 4``).

        Examples
        ========

        >>> from sympy import Matrix, I, eye
        >>> m = Matrix((0, 1 + I, 2, 3))
        >>> m.D
        Matrix([[0, 1 - I, -2, -3]])
        >>> m = (eye(4) + I*eye(4))
        >>> m[0, 3] = 2
        >>> m.D
        Matrix([
        [1 - I,     0,      0,      0],
        [    0, 1 - I,      0,      0],
        [    0,     0, -1 + I,      0],
        [    2,     0,      0, -1 + I]])

        If the matrix does not have 4 rows an AttributeError will be raised
        because this property is only defined for matrices with 4 rows.

        >>> Matrix(eye(2)).D
        Traceback (most recent call last):
        ...
        AttributeError: Matrix has no attribute D.

        See Also
        ========

        sympy.matrices.matrixbase.MatrixBase.conjugate: By-element conjugation
        sympy.matrices.matrixbase.MatrixBase.H: Hermite conjugation
        """
        from sympy.physics.matrices import mgamma
        if self.rows != 4:
            # In Python 3.2, properties can only return an AttributeError
            # so we can't raise a ShapeError -- see commit which added the
            # first line of this inline comment. Also, there is no need
            # for a message since MatrixBase will raise the AttributeError
            raise AttributeError
        return self.H * self._as_type(mgamma(0))

    def dot(self, b: MatrixBase | Expr,
                  hermitian: bool | None = None,
                  conjugate_convention: str | None = None
            ) -> Expr:
        """Return the dot or inner product of two vectors of equal length.
        Here ``self`` must be a ``Matrix`` of size 1 x n or n x 1, and ``b``
        must be either a matrix of size 1 x n, n x 1, or a list/tuple of length n.
        A scalar is returned.

        By default, ``dot`` does not conjugate ``self`` or ``b``, even if there are
        complex entries. Set ``hermitian=True`` (and optionally a ``conjugate_convention``)
        to compute the hermitian inner product.

        Possible kwargs are ``hermitian`` and ``conjugate_convention``.

        If ``conjugate_convention`` is ``"left"``, ``"math"`` or ``"maths"``,
        the conjugate of the first vector (``self``) is used.  If ``"right"``
        or ``"physics"`` is specified, the conjugate of the second vector ``b`` is used.

        Examples
        ========

        >>> from sympy import Matrix
        >>> M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> v = Matrix([1, 1, 1])
        >>> M.row(0).dot(v)
        6
        >>> M.col(0).dot(v)
        12
        >>> v = [3, 2, 1]
        >>> M.row(0).dot(v)
        10

        >>> from sympy import I
        >>> q = Matrix([1*I, 1*I, 1*I])
        >>> q.dot(q, hermitian=False)
        -3

        >>> q.dot(q, hermitian=True)
        3

        >>> q1 = Matrix([1, 1, 1*I])
        >>> q.dot(q1, hermitian=True, conjugate_convention="maths")
        1 - 2*I
        >>> q.dot(q1, hermitian=True, conjugate_convention="physics")
        1 + 2*I


        See Also
        ========

        cross
        multiply
        multiply_elementwise
        """
        from .dense import Matrix

        if not isinstance(b, MatrixBase):

            # XXX: Figure what types b can actually be here...

            if is_sequence(b):
                if len(b) != self.cols and len(b) != self.rows: # type: ignore
                    raise ShapeError(
                        "Dimensions incorrect for dot product: %s, %s" % (
                            self.shape, len(b))) # type: ignore
                return self.dot(Matrix(b)) # type: ignore
            else:
                raise TypeError(
                    "`b` must be an ordered iterable or Matrix, not %s." %
                    type(b))

        if (1 not in self.shape) or (1 not in b.shape):
            raise ShapeError
        if len(self) != len(b):
            raise ShapeError(
                "Dimensions incorrect for dot product: %s, %s" % (self.shape, b.shape))

        mat = self
        n = len(mat)
        if mat.shape != (1, n):
            mat = mat.reshape(1, n)
        if b.shape != (n, 1):
            b = b.reshape(n, 1)

        # Now ``mat`` is a row vector and ``b`` is a column vector.

        # If it so happens that only conjugate_convention is passed
        # then automatically set hermitian to True. If only hermitian
        # is true but no conjugate_convention is not passed then
        # automatically set it to ``"maths"``

        if conjugate_convention is not None and hermitian is None:
            hermitian = True
        if hermitian and conjugate_convention is None:
            conjugate_convention = "maths"

        if hermitian == True:
            if conjugate_convention in ("maths", "left", "math"):
                mat = mat.conjugate()
            elif conjugate_convention in ("physics", "right"):
                b = b.conjugate()
            else:
                raise ValueError("Unknown conjugate_convention was entered."
                                 " conjugate_convention must be one of the"
                                 " following: math, maths, left, physics or right.")
        return (mat * b)[0]

    def dual(self) -> MatrixBase:
        """Returns the dual of a matrix.

        A dual of a matrix is:

        ``(1/2)*levicivita(i, j, k, l)*M(k, l)`` summed over indices `k` and `l`

        Since the levicivita method is anti_symmetric for any pairwise
        exchange of indices, the dual of a symmetric matrix is the zero
        matrix. Strictly speaking the dual defined here assumes that the
        'matrix' `M` is a contravariant anti_symmetric second rank tensor,
        so that the dual is a covariant second rank tensor.

        """
        from sympy.matrices import zeros

        M, n = self[:, :], self.rows
        work = zeros(n)
        if self.is_symmetric():
            return work

        for i in range(1, n):
            for j in range(1, n):
                acum = S.Zero
                for k in range(1, n):
                    acum += LeviCivita(i, j, 0, k) * M[0, k]
                work[i, j] = acum
                work[j, i] = -acum

        for l in range(1, n):
            acum = S.Zero
            for a in range(1, n):
                for b in range(1, n):
                    acum += LeviCivita(0, l, a, b) * M[a, b]
            acum /= 2
            work[0, l] = -acum
            work[l, 0] = acum

        return work

    def _eval_matrix_exp_jblock(self) -> Self:
        """A helper function to compute an exponential of a Jordan block
        matrix

        Examples
        ========

        >>> from sympy import Symbol, Matrix
        >>> l = Symbol('lamda')

        A trivial example of 1*1 Jordan block:

        >>> m = Matrix.jordan_block(1, l)
        >>> m._eval_matrix_exp_jblock()
        Matrix([[exp(lamda)]])

        An example of 3*3 Jordan block:

        >>> m = Matrix.jordan_block(3, l)
        >>> m._eval_matrix_exp_jblock()
        Matrix([
        [exp(lamda), exp(lamda), exp(lamda)/2],
        [         0, exp(lamda),   exp(lamda)],
        [         0,          0,   exp(lamda)]])

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Matrix_function#Jordan_decomposition
        """
        size = self.rows
        l = self[0, 0]
        exp_l = exp(l)

        bands = {i: exp_l / factorial(i) for i in range(size)}

        from .sparsetools import banded
        return self._as_type(banded(size, bands))


    def analytic_func(self, f, x) -> Self:
        """
        Computes f(A) where A is a Square Matrix
        and f is an analytic function.

        Examples
        ========

        >>> from sympy import Symbol, Matrix, S, log

        >>> x = Symbol('x')
        >>> m = Matrix([[S(5)/4, S(3)/4], [S(3)/4, S(5)/4]])
        >>> f = log(x)
        >>> m.analytic_func(f, x)
        Matrix([
        [     0, log(2)],
        [log(2),      0]])

        Parameters
        ==========

        f : Expr
            Analytic Function
        x : Symbol
            parameter of f

        """

        f, x = _sympify(f), _sympify(x)
        if not self.is_square:
            raise NonSquareMatrixError
        if not x.is_symbol:
            raise ValueError("{} must be a symbol.".format(x))
        if x not in f.free_symbols:
            raise ValueError(
                "{} must be a parameter of {}.".format(x, f))
        if x in self.free_symbols:
            raise ValueError(
                "{} must not be a parameter of {}.".format(x, self))

        eigen = self.eigenvals()
        max_mul = max(eigen.values())
        derivative = {}
        dd = f
        for ii in range(max_mul - 1):
            dd = diff(dd, x)
            derivative[ii + 1] = dd
        n = self.shape[0]
        r: list[list[Expr]] = [[S.Zero]*n for _ in range(n)]
        f_val = [S.Zero] * n
        row = 0

        for i in eigen:
            mul = eigen[i]
            f_val[row] = f.subs(x, i)
            if f_val[row].is_number and not f_val[row].is_complex:
                raise ValueError(
                    "Cannot evaluate the function because the "
                    "function {} is not analytic at the given "
                    "eigenvalue {}".format(f, f_val[row]))
            val = S.One
            for a in range(n):
                r[row][a] = val
                val *= i
            if mul > 1:
                coe = [1] * n
                deri = 1
                while mul > 1:
                    row = row + 1
                    mul -= 1
                    d_i = derivative[deri].subs(x, i)
                    if d_i.is_number and not d_i.is_complex:
                        raise ValueError(
                            "Cannot evaluate the function because the "
                            "derivative {} is not analytic at the given "
                            "eigenvalue {}".format(derivative[deri], d_i))
                    f_val[row] = d_i
                    for a in range(n):
                        if a - deri + 1 <= 0:
                            r[row][a] = S.Zero
                            coe[a] = 0
                            continue
                        coe[a] = coe[a]*(a - deri + 1)
                        r[row][a] = coe[a]*pow(i, a - deri)
                    deri += 1
            row += 1

        c = self._new(r).solve(self._new(f_val))
        ans = self.zeros(n)
        pre = self.eye(n)
        for i2 in range(n):
            ans = ans + c[i2]*pre
            pre *= self
        return ans


    def exp(self):
        """Return the exponential of a square matrix.

        Examples
        ========

        >>> from sympy import Symbol, Matrix

        >>> t = Symbol('t')
        >>> m = Matrix([[0, 1], [-1, 0]]) * t
        >>> m.exp()
        Matrix([
        [    exp(I*t)/2 + exp(-I*t)/2, -I*exp(I*t)/2 + I*exp(-I*t)/2],
        [I*exp(I*t)/2 - I*exp(-I*t)/2,      exp(I*t)/2 + exp(-I*t)/2]])
        """
        if not self.is_square:
            raise NonSquareMatrixError(
                "Exponentiation is valid only for square matrices")
        try:
            P, J = self.jordan_form()
            cells = J.get_diag_blocks()
        except MatrixError:
            raise NotImplementedError(
                "Exponentiation is implemented only for matrices for which the Jordan normal form can be computed")

        blocks = [cell._eval_matrix_exp_jblock() for cell in cells]
        eJ = self.diag(*blocks)
        # n = self.rows
        ret = P.multiply(eJ, dotprodsimp=None).multiply(P.inv(), dotprodsimp=None)

        if all(value.is_real for value in self.values()):
            return self._as_type(re(ret))
        else:
            return self._as_type(ret)

    def _eval_matrix_log_jblock(self) -> Self:
        """Helper function to compute logarithm of a jordan block.

        Examples
        ========

        >>> from sympy import Symbol, Matrix
        >>> l = Symbol('lamda')

        A trivial example of 1*1 Jordan block:

        >>> m = Matrix.jordan_block(1, l)
        >>> m._eval_matrix_log_jblock()
        Matrix([[log(lamda)]])

        An example of 3*3 Jordan block:

        >>> m = Matrix.jordan_block(3, l)
        >>> m._eval_matrix_log_jblock()
        Matrix([
        [log(lamda),    1/lamda, -1/(2*lamda**2)],
        [         0, log(lamda),         1/lamda],
        [         0,          0,      log(lamda)]])
        """
        size = self.rows
        l = self[0, 0]

        if l.is_zero:
            raise MatrixError(
                'Could not take logarithm or reciprocal for the given '
                'eigenvalue {}'.format(l))

        bands = {0: log(l)}
        for i in range(1, size):
            bands[i] = -((-l) ** -i) / i

        from .sparsetools import banded
        return self._as_type(banded(size, bands))

    def log(self, simplify: Callable[[Self], Self] | Literal[False] = cancel) -> Self:
        """Return the logarithm of a square matrix.

        Parameters
        ==========

        simplify : function, bool
            The function to simplify the result with.

            Default is ``cancel``, which is effective to reduce the
            expression growing for taking reciprocals and inverses for
            symbolic matrices.

        Examples
        ========

        >>> from sympy import S, Matrix

        Examples for positive-definite matrices:

        >>> m = Matrix([[1, 1], [0, 1]])
        >>> m.log()
        Matrix([
        [0, 1],
        [0, 0]])

        >>> m = Matrix([[S(5)/4, S(3)/4], [S(3)/4, S(5)/4]])
        >>> m.log()
        Matrix([
        [     0, log(2)],
        [log(2),      0]])

        Examples for non positive-definite matrices:

        >>> m = Matrix([[S(3)/4, S(5)/4], [S(5)/4, S(3)/4]])
        >>> m.log()
        Matrix([
        [         I*pi/2, log(2) - I*pi/2],
        [log(2) - I*pi/2,          I*pi/2]])

        >>> m = Matrix(
        ...     [[0, 0, 0, 1],
        ...      [0, 0, 1, 0],
        ...      [0, 1, 0, 0],
        ...      [1, 0, 0, 0]])
        >>> m.log()
        Matrix([
        [ I*pi/2,       0,       0, -I*pi/2],
        [      0,  I*pi/2, -I*pi/2,       0],
        [      0, -I*pi/2,  I*pi/2,       0],
        [-I*pi/2,       0,       0,  I*pi/2]])
        """
        if not self.is_square:
            raise NonSquareMatrixError(
                "Logarithm is valid only for square matrices")

        try:
            if simplify:
                P, J = simplify(self).jordan_form()
            else:
                P, J = self.jordan_form()

            cells = J.get_diag_blocks()
        except MatrixError:
            raise NotImplementedError(
                "Logarithm is implemented only for matrices for which "
                "the Jordan normal form can be computed")

        blocks = [cell._eval_matrix_log_jblock() for cell in cells]

        eJ = self.diag(*blocks)

        if simplify:
            ret = simplify(P * eJ * simplify(P.inv()))
            ret = self._new(ret)
        else:
            ret = P * eJ * P.inv()

        return ret

    def is_nilpotent(self) -> bool:
        """Checks if a matrix is nilpotent.

        A matrix B is nilpotent if for some integer k, B**k is
        a zero matrix.

        Examples
        ========

        >>> from sympy import Matrix
        >>> a = Matrix([[0, 0, 0], [1, 0, 0], [1, 1, 0]])
        >>> a.is_nilpotent()
        True

        >>> a = Matrix([[1, 0, 1], [1, 0, 0], [1, 1, 0]])
        >>> a.is_nilpotent()
        False
        """
        if not self:
            return True
        if not self.is_square:
            raise NonSquareMatrixError(
                "Nilpotency is valid only for square matrices")
        x = uniquely_named_symbol('x', self, modify=lambda s: '_' + s)
        p = self.charpoly(x)
        if p.args[0] == x ** self.rows:
            return True
        return False

    def key2bounds(self, keys: tuple[int | slice, int | slice]) -> tuple[int, int, int, int]:
        """Converts a key with potentially mixed types of keys (integer and slice)
        into a tuple of ranges and raises an error if any index is out of ``self``'s
        range.

        See Also
        ========

        key2ij
        """
        i, j = keys
        if isinstance(i, slice):
            if not self.rows:
                rlo = rhi = 0
            else:
                rlo, rhi = i.indices(self.rows)[:2]
        else:
            rlo = a2idx(i, self.rows)
            rhi = rlo + 1
        if isinstance(j, slice):
            if not self.cols:
                clo = chi = 0
            else:
                clo, chi = j.indices(self.cols)[:2]
        else:
            clo = a2idx(j, self.cols)
            chi = clo + 1
        return rlo, rhi, clo, chi

    def key2ij(self, key: int | slice | tuple[int | slice, int | slice]) -> tuple[int, int]:
        """Converts key into canonical form, converting integers or indexable
        items into valid integers for ``self``'s range or returning slices
        unchanged.

        See Also
        ========

        key2bounds
        """
        if is_sequence(key):
            if not len(key) == 2: # type: ignore
                raise TypeError('key must be a sequence of length 2')
            return [a2idx(i, n) if not isinstance(i, slice) else i for i, n in zip(key, self.shape)] # type: ignore
        elif isinstance(key, slice): # type: ignore
            return key.indices(len(self))[:2]
        else:
            return divmod(a2idx(key, len(self)), self.cols)

    def normalized(self, iszerofunc: Callable[[Expr], bool | None] = _iszero) -> Self:
        """Return the normalized version of ``self``.

        Parameters
        ==========

        iszerofunc : Function, optional
            A function to determine whether ``self`` is a zero vector.
            The default ``_iszero`` tests to see if each element is
            exactly zero.

        Returns
        =======

        Matrix
            Normalized vector form of ``self``.
            It has the same length as a unit vector. However, a zero vector
            will be returned for a vector with norm 0.

        Raises
        ======

        ShapeError
            If the matrix is not in a vector form.

        See Also
        ========

        norm
        """
        if self.rows != 1 and self.cols != 1:
            raise ShapeError("A Matrix must be a vector to normalize.")
        norm = self.norm()
        if iszerofunc(norm):
            out = self.zeros(self.rows, self.cols)
        else:
            out = self.applyfunc(lambda i: i / norm)
        return out

    def norm(self, ord: int | str | None = 2) -> Expr:
        """Return the Norm of a Matrix or Vector.

        In the simplest case this is the geometric size of the vector
        Other norms can be specified by the ord parameter


        =====  ============================  ==========================
        ord    norm for matrices             norm for vectors
        =====  ============================  ==========================
        None   Frobenius norm                2-norm
        'fro'  Frobenius norm                - does not exist
        inf    maximum row sum               max(abs(x))
        -inf   --                            min(abs(x))
        1      maximum column sum            as below
        -1     --                            as below
        2      2-norm (largest sing. value)  as below
        -2     smallest singular value       as below
        other  - does not exist              sum(abs(x)**ord)**(1./ord)
        =====  ============================  ==========================

        Examples
        ========

        >>> from sympy import Matrix, Symbol, trigsimp, cos, sin, oo
        >>> x = Symbol('x', real=True)
        >>> v = Matrix([cos(x), sin(x)])
        >>> trigsimp( v.norm() )
        1
        >>> v.norm(10)
        (sin(x)**10 + cos(x)**10)**(1/10)
        >>> A = Matrix([[1, 1], [1, 1]])
        >>> A.norm(1) # maximum sum of absolute values of A is 2
        2
        >>> A.norm(2) # Spectral norm (max of |Ax|/|x| under 2-vector-norm)
        2
        >>> A.norm(-2) # Inverse spectral norm (smallest singular value)
        0
        >>> A.norm() # Frobenius Norm
        2
        >>> A.norm(oo) # Infinity Norm
        2
        >>> Matrix([1, -2]).norm(oo)
        2
        >>> Matrix([-1, 2]).norm(-oo)
        1

        See Also
        ========

        normalized
        """
        # Row or Column Vector Norms
        vals = list(self.values()) or [S.Zero]
        if S.One in self.shape:
            if ord == 2:  # Common case sqrt(<x, x>)
                return sqrt(Add(*(abs(i) ** 2 for i in vals)))

            elif ord == 1:  # sum(abs(x))
                return Add(*(abs(i) for i in vals))

            elif ord is S.Infinity:  # max(abs(x))
                return Max(*[abs(i) for i in vals])

            elif ord is S.NegativeInfinity:  # min(abs(x))
                return Min(*[abs(i) for i in vals])

            # Otherwise generalize the 2-norm, Sum(x_i**ord)**(1/ord)
            # Note that while useful this is not mathematically a norm
            try:
                return Pow(Add(*(abs(i) ** ord for i in vals)), S.One / ord) # type: ignore
            except (NotImplementedError, TypeError):
                raise ValueError("Expected order to be Number, Symbol, oo")

        # Matrix Norms
        else:
            if ord == 1:  # Maximum column sum
                m = self.applyfunc(abs)
                return Max(*[sum(m.col(i).flat()) for i in range(m.cols)])

            elif ord == 2:  # Spectral Norm
                # Maximum singular value
                return Max(*self.singular_values())

            elif ord == -2:
                # Minimum singular value
                return Min(*self.singular_values())

            elif ord is S.Infinity:   # Infinity Norm - Maximum row sum
                m = self.applyfunc(abs)
                return Max(*[sum(m.row(i).flat()) for i in range(m.rows)])

            elif (ord is None or isinstance(ord,
                                            str) and ord.lower() in
                ['f', 'fro', 'frobenius', 'vector']):
                # Reshape as vector and send back to norm function
                return self.vec().norm(ord=2)

            else:
                raise NotImplementedError("Matrix Norms under development")

    def print_nonzero(self, symb: str = "X") -> None:
        """Shows location of non-zero entries for fast shape lookup.

        Examples
        ========

        >>> from sympy import Matrix, eye
        >>> m = Matrix(2, 3, lambda i, j: i*3+j)
        >>> m
        Matrix([
        [0, 1, 2],
        [3, 4, 5]])
        >>> m.print_nonzero()
        [ XX]
        [XXX]
        >>> m = eye(4)
        >>> m.print_nonzero("x")
        [x   ]
        [ x  ]
        [  x ]
        [   x]

        """
        s = []
        for i in range(self.rows):
            line = []
            for j in range(self.cols):
                if self[i, j] == 0:
                    line.append(" ")
                else:
                    line.append(str(symb))
            s.append("[%s]" % ''.join(line))
        print('\n'.join(s))

    def project(self, v: MatrixBase) -> MatrixBase:
        """Return the projection of ``self`` onto the line containing ``v``.

        Examples
        ========

        >>> from sympy import Matrix, S, sqrt
        >>> V = Matrix([sqrt(3)/2, S.Half])
        >>> x = Matrix([[1, 0]])
        >>> V.project(x)
        Matrix([[sqrt(3)/2, 0]])
        >>> V.project(-x)
        Matrix([[sqrt(3)/2, 0]])
        """
        return v * (self.dot(v) / v.dot(v))

    def table(self, printer, rowstart: str='[', rowend: str=']', rowsep: str='\n',
              colsep: str=', ', align: str='right'):
        r"""
        String form of Matrix as a table.

        ``printer`` is the printer to use for on the elements (generally
        something like StrPrinter())

        ``rowstart`` is the string used to start each row (by default '[').

        ``rowend`` is the string used to end each row (by default ']').

        ``rowsep`` is the string used to separate rows (by default a newline).

        ``colsep`` is the string used to separate columns (by default ', ').

        ``align`` defines how the elements are aligned. Must be one of 'left',
        'right', or 'center'.  You can also use '<', '>', and '^' to mean the
        same thing, respectively.

        This is used by the string printer for Matrix.

        Examples
        ========

        >>> from sympy import Matrix, StrPrinter
        >>> M = Matrix([[1, 2], [-33, 4]])
        >>> printer = StrPrinter()
        >>> M.table(printer)
        '[  1, 2]\n[-33, 4]'
        >>> print(M.table(printer))
        [  1, 2]
        [-33, 4]
        >>> print(M.table(printer, rowsep=',\n'))
        [  1, 2],
        [-33, 4]
        >>> print('[%s]' % M.table(printer, rowsep=',\n'))
        [[  1, 2],
        [-33, 4]]
        >>> print(M.table(printer, colsep=' '))
        [  1 2]
        [-33 4]
        >>> print(M.table(printer, align='center'))
        [ 1 , 2]
        [-33, 4]
        >>> print(M.table(printer, rowstart='{', rowend='}'))
        {  1, 2}
        {-33, 4}
        """
        # Handle zero dimensions:
        if S.Zero in self.shape:
            return '[]'
        # Build table of string representations of the elements
        table: list[list[str]] = []
        # Track per-column max lengths for pretty alignment
        maxlen = [0] * self.cols
        for i in range(self.rows):
            table.append([])
            for j in range(self.cols):
                s = printer._print(self[i, j])
                table[-1].append(s)
                maxlen[j] = max(len(s), maxlen[j])
        # Patch strings together
        align = {
            'left': 'ljust',
            'right': 'rjust',
            'center': 'center',
            '<': 'ljust',
            '>': 'rjust',
            '^': 'center',
        }[align]
        res = [""] * len(table)
        for i, row in enumerate(table):
            for j, elem in enumerate(row):
                row[j] = getattr(elem, align)(maxlen[j])
            res[i] = rowstart + colsep.join(row) + rowend
        return rowsep.join(res)

    def rank_decomposition(self,
                           iszerofunc: Callable[[Expr], bool | None] = _iszero,
                           simplify: bool | Callable[[Expr], Expr] = False
                        ) -> tuple[Self, Self]:
        return _rank_decomposition(self, iszerofunc=iszerofunc, simplify=simplify)

    @abstractmethod
    def cholesky(self, hermitian: bool = True) -> Self:
        raise NotImplementedError('This function is implemented in DenseMatrix or SparseMatrix')

    @abstractmethod
    def LDLdecomposition(self, hermitian: bool = True) -> tuple[Self, Self]:
        raise NotImplementedError('This function is implemented in DenseMatrix or SparseMatrix')

    def LUdecomposition(self,
                        iszerofunc: Callable[[Expr], bool | None] = _iszero,
                        simpfunc: Callable[[Expr], Expr] | None = None,
                        rankcheck: bool = False) -> tuple[Self, Self, list[list[int]]]:
        return _LUdecomposition(self, iszerofunc=iszerofunc, simpfunc=simpfunc,
                rankcheck=rankcheck)

    def LUdecomposition_Simple(self,
                               iszerofunc: Callable[[Expr], bool | None] = _iszero,
                               simpfunc: Callable[[Expr], Expr] | None = None,
                               rankcheck: bool = False) -> tuple[Self, list[list[int]]]:
        return _LUdecomposition_Simple(self, iszerofunc=iszerofunc,
                simpfunc=simpfunc, rankcheck=rankcheck)

    def LUdecompositionFF(self) -> tuple[Self, Self, Self, Self]:
        P, L, D, U = _LUdecompositionFF(self)
        return self._as_type(P), self._as_type(L), self._as_type(D), self._as_type(U)

    def singular_value_decomposition(self) -> tuple[Self, Self, Self]:
        return _singular_value_decomposition(self)

    def QRdecomposition(self) -> tuple[Self, Self]:
        return _QRdecomposition(self)

    def upper_hessenberg_decomposition(self) -> tuple[Self, Self]:
        return _upper_hessenberg_decomposition(self)

    def diagonal_solve(self, rhs) -> Self:
        return _diagonal_solve(self, rhs)

    @abstractmethod
    def lower_triangular_solve(self, rhs: MatrixBase) -> Self:
        raise NotImplementedError('This function is implemented in DenseMatrix or SparseMatrix')

    @abstractmethod
    def upper_triangular_solve(self, rhs: MatrixBase) -> Self:
        raise NotImplementedError('This function is implemented in DenseMatrix or SparseMatrix')

    def cholesky_solve(self, rhs) -> Self:
        return _cholesky_solve(self, rhs)

    def LDLsolve(self, rhs) -> Self:
        return _LDLsolve(self, rhs)

    def LUsolve(self, rhs, iszerofunc=_iszero) -> Self:
        return _LUsolve(self, rhs, iszerofunc=iszerofunc)

    def QRsolve(self, b) -> Self:
        return _QRsolve(self, b)

    @overload
    def gauss_jordan_solve(self, B: MatrixBase, freevar: Literal[False] = False,
                           ) -> tuple[Self, Self]: ...
    @overload
    def gauss_jordan_solve(self, B: MatrixBase, freevar: Literal[True],
                           ) -> tuple[Self, Self, list[int]]: ...

    def gauss_jordan_solve(self, B, freevar: bool=False
                           ) -> tuple[Self, Self] | tuple[Self, Self, list[int]]:
        return _gauss_jordan_solve(self, B, freevar=freevar)

    def pinv_solve(self, B, arbitrary_matrix=None) -> Self:
        return _pinv_solve(self, B, arbitrary_matrix=arbitrary_matrix)

    def cramer_solve(self, rhs: Self, det_method: str="laplace") -> Self:
        return _cramer_solve(self, rhs, det_method=det_method)

    def solve(self, rhs: Self, method: str='GJ') -> Self:
        return _solve(self, rhs, method=method)

    def solve_least_squares(self, rhs: Self, method: str='CH') -> Self:
        return _solve_least_squares(self, rhs, method=method)

    def pinv(self, method: str = 'RD') -> Self:
        return _pinv(self, method=method)

    def inverse_ADJ(self, iszerofunc: Callable[[Expr], bool | None] = _iszero) -> Self:
        return _inv_ADJ(self, iszerofunc=iszerofunc)

    def inverse_BLOCK(self, iszerofunc: Callable[[Expr], bool | None] = _iszero) -> Self:
        return self._as_type(_inv_block(self, iszerofunc=iszerofunc))

    def inverse_GE(self, iszerofunc=_iszero) -> Self:
        return _inv_GE(self, iszerofunc=iszerofunc)

    def inverse_LU(self, iszerofunc=_iszero) -> Self:
        return _inv_LU(self, iszerofunc=iszerofunc)

    def inverse_CH(self, iszerofunc=_iszero) -> Self:
        return _inv_CH(self, iszerofunc=iszerofunc)

    def inverse_LDL(self, iszerofunc=_iszero) -> Self:
        return _inv_LDL(self, iszerofunc=iszerofunc)

    def inverse_QR(self, iszerofunc=_iszero) -> Self:
        return _inv_QR(self, iszerofunc=iszerofunc)

    def inv(self, method: str | None = None,
                iszerofunc: Callable[[Expr], bool | None] = _iszero,
                try_block_diag: bool = False
            ) -> Self:
        return _inv(self, method=method, iszerofunc=iszerofunc,
                try_block_diag=try_block_diag)

    def connected_components(self) -> list[list[int]]:
        return _connected_components(self)

    def connected_components_decomposition(self
                            ) -> tuple[PermutationMatrix, BlockDiagMatrix]:
        return _connected_components_decomposition(self)

    def strongly_connected_components(self) -> list[list[int]]:
        return _strongly_connected_components(self)

    def strongly_connected_components_decomposition(self, lower: bool = True,
                            ) -> tuple[PermutationMatrix, BlockMatrix]:
        return _strongly_connected_components_decomposition(self, lower=lower)

    _sage_ = Basic._sage_

    rank_decomposition.__doc__     = _rank_decomposition.__doc__
    cholesky.__doc__               = _cholesky.__doc__
    LDLdecomposition.__doc__       = _LDLdecomposition.__doc__
    LUdecomposition.__doc__        = _LUdecomposition.__doc__
    LUdecomposition_Simple.__doc__ = _LUdecomposition_Simple.__doc__
    LUdecompositionFF.__doc__      = _LUdecompositionFF.__doc__
    singular_value_decomposition.__doc__ = _singular_value_decomposition.__doc__
    QRdecomposition.__doc__        = _QRdecomposition.__doc__
    upper_hessenberg_decomposition.__doc__ = _upper_hessenberg_decomposition.__doc__

    diagonal_solve.__doc__         = _diagonal_solve.__doc__
    lower_triangular_solve.__doc__ = _lower_triangular_solve.__doc__
    upper_triangular_solve.__doc__ = _upper_triangular_solve.__doc__
    cholesky_solve.__doc__         = _cholesky_solve.__doc__
    LDLsolve.__doc__               = _LDLsolve.__doc__
    LUsolve.__doc__                = _LUsolve.__doc__
    QRsolve.__doc__                = _QRsolve.__doc__
    gauss_jordan_solve.__doc__     = _gauss_jordan_solve.__doc__
    pinv_solve.__doc__             = _pinv_solve.__doc__
    cramer_solve.__doc__           = _cramer_solve.__doc__
    solve.__doc__                  = _solve.__doc__
    solve_least_squares.__doc__    = _solve_least_squares.__doc__

    pinv.__doc__                   = _pinv.__doc__
    inverse_ADJ.__doc__            = _inv_ADJ.__doc__
    inverse_GE.__doc__             = _inv_GE.__doc__
    inverse_LU.__doc__             = _inv_LU.__doc__
    inverse_CH.__doc__             = _inv_CH.__doc__
    inverse_LDL.__doc__            = _inv_LDL.__doc__
    inverse_QR.__doc__             = _inv_QR.__doc__
    inverse_BLOCK.__doc__          = _inv_block.__doc__
    inv.__doc__                    = _inv.__doc__

    connected_components.__doc__   = _connected_components.__doc__
    connected_components_decomposition.__doc__ = \
        _connected_components_decomposition.__doc__
    strongly_connected_components.__doc__   = \
        _strongly_connected_components.__doc__
    strongly_connected_components_decomposition.__doc__ = \
        _strongly_connected_components_decomposition.__doc__


def _convert_matrix(typ, mat):
    """Convert mat to a Matrix of type typ."""
    from sympy.matrices.matrixbase import MatrixBase
    if getattr(mat, "is_Matrix", False) and not isinstance(mat, MatrixBase):
        # This is needed for interop between Matrix and the redundant matrix
        # mixin types like _MinimalMatrix etc. If anyone should happen to be
        # using those then this keeps them working. Really _MinimalMatrix etc
        # should be deprecated and removed though.
        return typ(*mat.shape, list(mat))
    else:
        return typ(mat)


def _has_matrix_shape(other):
    shape = getattr(other, 'shape', None)
    if shape is None:
        return False
    return isinstance(shape, tuple) and len(shape) == 2


def _has_rows_cols(other):
    return hasattr(other, 'rows') and hasattr(other, 'cols')


def _coerce_operand(self, other: Any
    ) -> (tuple[None, Literal['invalid_type']]
          | tuple[MatrixBase, Literal['is_matrix']]
          | tuple[Expr, Literal['possible_scalar']]):
    """Convert other to a Matrix, or check for possible scalar."""

    # Disallow mixing Matrix and Array
    if isinstance(other, NDimArray):
        return None, 'invalid_type'

    is_Matrix = getattr(other, 'is_Matrix', None)

    # Return a Matrix as-is
    if is_Matrix:
        return other, 'is_matrix'

    # Try to convert numpy array, mpmath matrix etc.
    if is_Matrix is None:
        if _has_matrix_shape(other) or _has_rows_cols(other):
            return _convert_matrix(type(self), other), 'is_matrix'

    # Could be a scalar but only if not iterable...
    if not isinstance(other, Iterable):
        return other, 'possible_scalar'

    return None, 'invalid_type'


def classof(A: Tmat, B: Tmat) -> type[Tmat]:
    """
    Get the type of the result when combining matrices of different types.

    Currently the strategy is that immutability is contagious.

    Examples
    ========

    >>> from sympy import Matrix, ImmutableMatrix
    >>> from sympy.matrices.matrixbase import classof
    >>> M = Matrix([[1, 2], [3, 4]]) # a Mutable Matrix
    >>> IM = ImmutableMatrix([[1, 2], [3, 4]])
    >>> classof(M, IM)
    <class 'sympy.matrices.immutable.ImmutableDenseMatrix'>
    """
    priority_A = getattr(A, '_class_priority', None)
    priority_B = getattr(B, '_class_priority', None)
    if None not in (priority_A, priority_B):
        if A._class_priority > B._class_priority:
            return A.__class__
        else:
            return B.__class__

    try:
        import numpy
    except ImportError:
        pass
    else:
        if isinstance(A, numpy.ndarray):
            return B.__class__
        if isinstance(B, numpy.ndarray):
            return A.__class__

    raise TypeError("Incompatible classes %s, %s" % (A.__class__, B.__class__))

@overload
def _unify_with_other(self: Tmat, other: Tmat) -> tuple[Tmat, Tmat, Literal['is_matrix']]: ...
@overload
def _unify_with_other(self: Tmat, other: SExpr) -> tuple[Tmat, Expr, Literal['possible_scalar']]: ...

def _unify_with_other(self: MatrixBase, other: Any
            ) -> tuple[MatrixBase, MatrixBase | Expr, str]:
    """Unify self and other into a single matrix type, or check for scalar."""
    other, T = _coerce_operand(self, other)

    if isinstance(other, MatrixBase):
        typ = classof(self, other)
        if typ != self.__class__:
            self = _convert_matrix(typ, self)
        if typ != other.__class__:
            other = _convert_matrix(typ, other)

    return self, other, T


def a2idx(j, n=None):
    """Return integer after making positive and validating against n."""
    if not isinstance(j, int):
        jindex = getattr(j, '__index__', None)
        if jindex is not None:
            j = jindex()
        else:
            raise IndexError("Invalid index a[%r]" % (j,))
    if n is not None:
        if j < 0:
            j += n
        if not (j >= 0 and j < n):
            raise IndexError("Index out of range: a[%s]" % (j,))
    return int(j)


class DeferredVector(Symbol, NotIterable): # type: ignore
    """A vector whose components are deferred (e.g. for use with lambdify).

    Examples
    ========

    >>> from sympy import DeferredVector, lambdify
    >>> X = DeferredVector( 'X' )
    >>> X
    X
    >>> expr = (X[0] + 2, X[2] + 3)
    >>> func = lambdify( X, expr)
    >>> func( [1, 2, 3] )
    (3, 6)
    """

    def __getitem__(self, i):
        if i == -0:
            i = 0
        if i < 0:
            raise IndexError('DeferredVector index out of range')
        component_name = '%s[%d]' % (self.name, i)
        return Symbol(component_name)

    def __str__(self):
        return sstr(self)

    def __repr__(self):
        return "DeferredVector('%s')" % self.name
