from __future__ import annotations

from typing import TYPE_CHECKING, overload

if TYPE_CHECKING:
    from typing import Literal, Callable, TypeVar
    from sympy.core.expr import Expr
    from sympy.matrices.matrixbase import MatrixBase
    Tmat = TypeVar('Tmat', bound=MatrixBase)

from types import FunctionType

from sympy.polys.polyerrors import CoercionFailed
from sympy.polys.domains import ZZ, QQ

from .utilities import _get_intermediate_simp, _iszero, _dotprodsimp, _simplify
from .determinant import _find_reasonable_pivot


if TYPE_CHECKING:
    from typing import Any


def _row_reduce_list(mat, rows, cols, one, iszerofunc, simpfunc,
                normalize_last=True, normalize=True, zero_above=True):
    """Row reduce a flat list representation of a matrix and return a tuple
    (rref_matrix, pivot_cols, swaps) where ``rref_matrix`` is a flat list,
    ``pivot_cols`` are the pivot columns and ``swaps`` are any row swaps that
    were used in the process of row reduction.

    Parameters
    ==========

    mat : list
        list of matrix elements, must be ``rows`` * ``cols`` in length

    rows, cols : integer
        number of rows and columns in flat list representation

    one : SymPy object
        represents the value one, from ``Matrix.one``

    iszerofunc : determines if an entry can be used as a pivot

    simpfunc : used to simplify elements and test if they are
        zero if ``iszerofunc`` returns `None`

    normalize_last : indicates where all row reduction should
        happen in a fraction-free manner and then the rows are
        normalized (so that the pivots are 1), or whether
        rows should be normalized along the way (like the naive
        row reduction algorithm)

    normalize : whether pivot rows should be normalized so that
        the pivot value is 1

    zero_above : whether entries above the pivot should be zeroed.
        If ``zero_above=False``, an echelon matrix will be returned.
    """

    def get_col(i):
        return mat[i::cols]

    def row_swap(i, j):
        mat[i*cols:(i + 1)*cols], mat[j*cols:(j + 1)*cols] = \
            mat[j*cols:(j + 1)*cols], mat[i*cols:(i + 1)*cols]

    def cross_cancel(a, i, b, j):
        """Does the row op row[i] = a*row[i] - b*row[j]"""
        q = (j - i)*cols
        for p in range(i*cols, (i + 1)*cols):
            mat[p] = isimp(a*mat[p] - b*mat[p + q])

    isimp = _get_intermediate_simp(_dotprodsimp)
    piv_row, piv_col = 0, 0
    pivot_cols = []
    swaps = []

    # use a fraction free method to zero above and below each pivot
    while piv_col < cols and piv_row < rows:
        pivot_offset, pivot_val, \
        assumed_nonzero, newly_determined = _find_reasonable_pivot(
                get_col(piv_col)[piv_row:], iszerofunc, simpfunc)

        # _find_reasonable_pivot may have simplified some things
        # in the process.  Let's not let them go to waste
        for (offset, val) in newly_determined:
            offset += piv_row
            mat[offset*cols + piv_col] = val

        if pivot_offset is None:
            piv_col += 1
            continue

        pivot_cols.append(piv_col)
        if pivot_offset != 0:
            row_swap(piv_row, pivot_offset + piv_row)
            swaps.append((piv_row, pivot_offset + piv_row))

        # if we aren't normalizing last
        # or the pivot_val is non-commutative,
        # we normalize before we zero the other rows
        if normalize_last is False or not pivot_val.is_commutative:
            i, j = piv_row, piv_col
            mat[i*cols + j] = one
            for p in range(i*cols + j + 1, (i + 1)*cols):
                mat[p] = isimp(pivot_val**(-1) * mat[p])
            # after normalizing, the pivot value is 1
            pivot_val = one

        # zero above and below the pivot
        for row in range(rows):
            # don't zero our current row
            if row == piv_row:
                continue
            # don't zero above the pivot unless we're told.
            if zero_above is False and row < piv_row:
                continue
            # if we're already a zero, don't do anything
            val = mat[row*cols + piv_col]
            if iszerofunc(val):
                continue

            cross_cancel(pivot_val, row, val, piv_row)
        piv_row += 1

    # normalize each row
    if normalize_last is True and normalize is True:
        for piv_i, piv_j in enumerate(pivot_cols):
            pivot_val = mat[piv_i*cols + piv_j]
            mat[piv_i*cols + piv_j] = one
            for p in range(piv_i*cols + piv_j + 1, (piv_i + 1)*cols):
                mat[p] = isimp(pivot_val**(-1) * mat[p])

    return mat, tuple(pivot_cols), tuple(swaps)


# This functions is a candidate for caching if it gets implemented for matrices.
def _row_reduce(M, iszerofunc, simpfunc, normalize_last=True,
                normalize=True, zero_above=True):

    mat, pivot_cols, swaps = _row_reduce_list(list(M), M.rows, M.cols, M.one,
            iszerofunc, simpfunc, normalize_last=normalize_last,
            normalize=normalize, zero_above=zero_above)

    return M._new(M.rows, M.cols, mat), pivot_cols, swaps


def _is_echelon(M, iszerofunc=_iszero):
    """Returns `True` if the matrix is in echelon form. That is, all rows of
    zeros are at the bottom, and below each leading non-zero in a row are
    exclusively zeros."""

    if M.rows <= 0 or M.cols <= 0:
        return True

    zeros_below = all(iszerofunc(t) for t in M[1:, 0])

    if iszerofunc(M[0, 0]):
        return zeros_below and _is_echelon(M[:, 1:], iszerofunc)

    return zeros_below and _is_echelon(M[1:, 1:], iszerofunc)


@overload
def _echelon_form(M: Tmat,
        iszerofunc: Callable[[Expr], bool | None] = _iszero,
        simplify: bool | Callable[[Expr], Expr] = False,
        *,
        with_pivots: Literal[False] = False
    ) -> Tmat: ...
@overload
def _echelon_form(M: Tmat,
        iszerofunc: Callable[[Expr], bool | None] = _iszero,
        simplify: bool | Callable[[Expr], Expr] = False,
        *,
        with_pivots: Literal[True],
    ) -> tuple[Tmat, tuple[int]]: ...
@overload
def _echelon_form(M: Tmat,
        iszerofunc: Callable[[Expr], bool | None] = _iszero,
        simplify: bool | Callable[[Expr], Expr] = False,
        *,
        with_pivots: bool = False,
    ) -> Tmat | tuple[Tmat, tuple[int]]: ...

def _echelon_form(M: Tmat,
        iszerofunc: Callable[[Expr], bool | None] = _iszero,
        simplify: bool | Callable[[Expr], Expr] = False,
        *,
        with_pivots: bool = False,
    ) -> Tmat | tuple[Tmat, tuple[int]]:
    """Returns a matrix row-equivalent to ``M`` that is in echelon form. Note
    that echelon form of a matrix is *not* unique, however, properties like the
    row space and the null space are preserved.

    Examples
    ========

    >>> from sympy import Matrix
    >>> M = Matrix([[1, 2], [3, 4]])
    >>> M.echelon_form()
    Matrix([
    [1,  2],
    [0, -2]])
    """

    simpfunc = simplify if isinstance(simplify, FunctionType) else _simplify

    mat, pivots, _ = _row_reduce(M, iszerofunc, simpfunc,
            normalize_last=True, normalize=False, zero_above=False)

    if with_pivots:
        return mat, pivots

    return mat


# This functions is a candidate for caching if it gets implemented for matrices.
def _rank(M, iszerofunc=_iszero, simplify=False):
    """Returns the rank of a matrix.

    Examples
    ========

    >>> from sympy import Matrix
    >>> from sympy.abc import x
    >>> m = Matrix([[1, 2], [x, 1 - 1/x]])
    >>> m.rank()
    2
    >>> n = Matrix(3, 3, range(1, 10))
    >>> n.rank()
    2
    """

    def _permute_complexity_right(M, iszerofunc):
        """Permute columns with complicated elements as
        far right as they can go.  Since the ``sympy`` row reduction
        algorithms start on the left, having complexity right-shifted
        speeds things up.

        Returns a tuple (mat, perm) where perm is a permutation
        of the columns to perform to shift the complex columns right, and mat
        is the permuted matrix."""

        def complexity(i):
            # the complexity of a column will be judged by how many
            # element's zero-ness cannot be determined
            return sum(1 if iszerofunc(e) is None else 0 for e in M[:, i])

        complex = [(complexity(i), i) for i in range(M.cols)]
        perm    = [j for (i, j) in sorted(complex)]

        return (M.permute(perm, orientation='cols'), perm)

    simpfunc = simplify if isinstance(simplify, FunctionType) else _simplify

    # for small matrices, we compute the rank explicitly
    # if is_zero on elements doesn't answer the question
    # for small matrices, we fall back to the full routine.
    if M.rows <= 0 or M.cols <= 0:
        return 0

    if M.rows <= 1 or M.cols <= 1:
        zeros = [iszerofunc(x) for x in M]

        if False in zeros:
            return 1

    if M.rows == 2 and M.cols == 2:
        zeros = [iszerofunc(x) for x in M]

        if False not in zeros and None not in zeros:
            return 0

        d = M.det()

        if iszerofunc(d) and False in zeros:
            return 1
        if iszerofunc(d) is False:
            return 2

    mat, _       = _permute_complexity_right(M, iszerofunc=iszerofunc)
    _, pivots, _ = _row_reduce(mat, iszerofunc, simpfunc, normalize_last=True,
            normalize=False, zero_above=False)

    return len(pivots)


def _to_DM_ZZ_QQ(M):
    # We have to test for _rep here because there are tests that otherwise fail
    # with e.g. "AttributeError: 'SubspaceOnlyMatrix' object has no attribute
    # '_rep'." There is almost certainly no value in such tests. The
    # presumption seems to be that someone could create a new class by
    # inheriting some of the Matrix classes and not the full set that is used
    # by the standard Matrix class but if anyone tried that it would fail in
    # many ways.
    if not hasattr(M, '_rep'):
        return None

    rep = M._rep
    K = rep.domain

    if K.is_ZZ:
        return rep
    elif K.is_QQ:
        try:
            return rep.convert_to(ZZ)
        except CoercionFailed:
            return rep
    else:
        if not all(e.is_Rational for e in M):
            return None
        try:
            return rep.convert_to(ZZ)
        except CoercionFailed:
            return rep.convert_to(QQ)


def _rref_dm(dM):
    """Compute the reduced row echelon form of a DomainMatrix."""
    K = dM.domain

    if K.is_ZZ:
        dM_rref, den, pivots = dM.rref_den(keep_domain=False)
        dM_rref = dM_rref.to_field() / den
    elif K.is_QQ:
        dM_rref, pivots = dM.rref()
    else:
        assert False  # pragma: no cover

    M_rref = dM_rref.to_Matrix()

    return M_rref, pivots


@overload
def _rref(
        M: Tmat,
        iszerofunc: Callable[[Expr], bool | None] = _iszero,
        simplify: bool | Callable[[Expr], Expr] = False,
        *,
        pivots: Literal[False],
        normalize_last: bool = True,
        ) -> Tmat: ...
@overload
def _rref(
        M: Tmat,
        iszerofunc: Callable[[Expr], bool | None] = _iszero,
        simplify: bool | Callable[[Expr], Expr] = False,
        *,
        pivots: Literal[True] = True,
        normalize_last: bool = True,
    ) -> tuple[Tmat, tuple[int]]: ...
@overload
def _rref(
        M: Tmat,
        iszerofunc: Callable[[Expr], bool | None] = _iszero,
        simplify: bool | Callable[[Expr], Expr] = False,
        *,
        pivots: bool = True,
        normalize_last: bool = True,
    ) -> Tmat | tuple[Tmat, tuple[int]]: ...

def _rref(
        M: Tmat,
        iszerofunc: Callable[[Expr], bool | None] = _iszero,
        simplify: bool | Callable[[Expr], Expr] = False,
        pivots: bool = True,
        normalize_last: bool = True,
    ) -> Tmat | tuple[Tmat, tuple[int]]:
    """Return the reduced row-echelon form (RREF) of a matrix, together
    with the column indices of the pivot columns.

    **What is RREF?**

    A matrix is in *reduced row-echelon form* when all four of the
    following conditions hold:

    1. All rows consisting entirely of zeros are at the bottom.
    2. The first non-zero entry in every non-zero row (called the
       *pivot* or *leading entry*) is 1.
    3. Each pivot lies strictly to the right of the pivot in the row
       above it.
    4. Every entry directly *above and below* each pivot is 0.

    Unlike plain echelon form, the RREF of a matrix is **unique** -- there
    is exactly one RREF for any given matrix.

    Parameters
    ==========

    iszerofunc : callable, optional
        A function that accepts a single expression and returns:

        * ``True``  -- the entry is zero (may be used as a zero during
          row operations, but not chosen as a pivot).
        * ``False`` -- the entry is definitely non-zero (can be a pivot).
        * ``None``  -- cannot be determined; SymPy will try to simplify
          further before deciding.

        The default is ``lambda x: x.is_zero``, which works perfectly
        for exact symbolic and rational matrices.

        The most common reason to override ``iszerofunc`` is **floating-
        point rounding**: arithmetic on ``float`` entries can produce tiny
        residuals like ``2.3e-17`` that are logically zero but not detected
        as such by the default check.  Supplying a tolerance-based
        function like ``lambda x: abs(x) < 1e-9`` tells SymPy to treat
        those small values as zero.

        For purely integer, rational, or symbolic matrices you never need
        to change this argument.

    simplify : bool or callable, optional
        Controls whether algebraic simplification is applied to entries
        during the pivot search.

        * ``False`` (default) -- no extra simplification is performed
          while searching for a pivot.  SymPy's ``simplify`` is still
          called internally if ``iszerofunc`` returns ``None`` for a
          candidate pivot.
        * ``True`` -- use SymPy's built-in :func:`~sympy.simplify.simplify`
          during every pivot check.
        * A callable -- use that function (e.g. ``trigsimp``,
          ``radsimp``) instead of the default ``simplify``.

        You rarely need this for numeric matrices.  It is helpful when
        your matrix contains symbolic expressions whose zero-ness cannot
        be decided without transformation (e.g. ``sin(x)**2 + cos(x)**2
        - 1``).

    pivots : bool, optional
        Determines what the function returns.

        * ``True`` (default) -- return a 2-tuple
          ``(rref_matrix, pivot_columns)``.
        * ``False`` -- return only ``rref_matrix``, with no pivot
          information.

        Pass ``pivots=False`` when you only need the reduced matrix and
        the column indices are not important to you.

    normalize_last : bool, optional
        Controls the internal order in which row operations are applied.

        * ``True`` (default) -- use a **fraction-free** algorithm that
          postpones dividing each row by its pivot until *after* all
          entries above and below that pivot have been zeroed.  This
          avoids large intermediate fractions and is significantly faster
          for matrices with symbols or large integers.
        * ``False`` -- use the classic textbook algorithm where each
          pivot is normalized to 1 *before* the surrounding rows are
          cleared.  The final matrix is **identical**; only the
          intermediate steps differ.

        .. note::

           When the matrix contains only integers or rationals, SymPy
           routes the computation through a specialised ``DomainMatrix``
           backend and ``normalize_last`` has no effect.

    Returns
    =======

    (rref_matrix, pivot_columns) : tuple[Matrix, tuple[int, ...]]
        Returned when ``pivots=True`` (the default).

        - ``rref_matrix`` -- a new :class:`~sympy.matrices.Matrix` that
          is row-equivalent to the input and is in reduced row-echelon
          form.
        - ``pivot_columns`` -- a tuple of **0-based** column indices
          where the leading 1 of each non-zero row appears.  The length
          of this tuple equals the *rank* of the original matrix.

    rref_matrix : Matrix
        Returned when ``pivots=False``.  Only the reduced matrix.

    Examples
    ========

    **Basic usage**

    >>> from sympy import Matrix
    >>> M = Matrix([[1, 2, 3],
    ...             [4, 5, 6],
    ...             [7, 8, 9]])
    >>> rref_M, pivots = M.rref()
    >>> rref_M
    Matrix([
    [1, 0, -1],
    [0, 1,  2],
    [0, 0,  0]])
    >>> pivots
    (0, 1)

    Two pivot columns mean the rank of ``M`` is 2; column 2 is a free
    variable column.

    **Getting only the matrix with** ``pivots=False``

    >>> M.rref(pivots=False)
    Matrix([
    [1, 0, -1],
    [0, 1,  2],
    [0, 0,  0]])

    **Symbolic matrix**

    >>> from sympy.abc import x
    >>> S = Matrix([[1, 2],
    ...             [x, 1 - 1/x]])
    >>> S.rref()
    (Matrix([
    [1, 0],
    [0, 1]]), (0, 1))

    **Solving a linear system with an augmented matrix**

    Build ``[A | b]`` and reduce it to read off the solution directly::

        x = 2, y = 3, z = -1

    >>> A = Matrix([[2,  1, -1],
    ...             [-3, -1,  2],
    ...             [-2,  1,  2]])
    >>> b = Matrix([8, -11, -3])
    >>> A.row_join(b).rref()
    (Matrix([
    [1, 0, 0,  2],
    [0, 1, 0,  3],
    [0, 0, 1, -1]]), (0, 1, 2))

    **Using** ``iszerofunc`` **to handle floating-point rounding errors**

    Matrices built from ``float`` values can have tiny rounding residuals
    that cause a rank-deficient matrix to look full-rank.  The first call
    below incorrectly returns rank 3; the second, with a tolerance, gives
    the correct rank-2 result:

    >>> m = Matrix([[0.9, -0.1, -0.2, 0],
    ...             [-0.8, 0.9, -0.4, 0],
    ...             [-0.1, -0.8, 0.6, 0]])
    >>> m.rref()
    (Matrix([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0]]), (0, 1, 2))
    >>> m.rref(iszerofunc=lambda x: abs(x) < 1e-9)
    (Matrix([
    [1, 0, -0.301369863013699, 0],
    [0, 1, -0.712328767123288, 0],
    [0, 0,                  0, 0]]), (0, 1))

    **Using a custom** ``simplify`` **function**

    When entries involve trigonometric expressions, pass ``trigsimp`` so
    that identities like ``sin(x)**2 + cos(x)**2`` are recognised as 1:

    >>> from sympy import trigsimp, sin, cos, pi
    >>> T = Matrix([[sin(pi/4), cos(pi/4)],
    ...             [cos(pi/4), -sin(pi/4)]])
    >>> T.rref(simplify=trigsimp)
    (Matrix([
    [1, 0],
    [0, 1]]), (0, 1))

    **Effect of** ``normalize_last``

    Both settings produce the same final matrix; ``normalize_last=True``
    (the default) is typically faster for symbolic entries:

    >>> from sympy import symbols
    >>> a, b = symbols('a b')
    >>> N = Matrix([[a, 1],
    ...             [1, b]])
    >>> N.rref(normalize_last=True)   # default, fraction-free internally
    (Matrix([
    [1, 0],
    [0, 1]]), (0, 1))
    >>> N.rref(normalize_last=False)  # classic textbook order
    (Matrix([
    [1, 0],
    [0, 1]]), (0, 1))

    Notes
    =====

    The default ``normalize_last=True`` provides a significant speedup
    for matrices with symbolic entries because it keeps entries as
    integers (or polynomials) throughout the elimination phase and only
    divides at the very end.  If you need to observe exactly how
    intermediate rows look during reduction, set ``normalize_last=False``
    to follow the textbook algorithm step-by-step.

    See Also
    ========

    sympy.matrices.matrixbase.MatrixBase.echelon_form : row-echelon form
        (not fully reduced; pivots need not be 1 and entries above pivots
        need not be zero).
    sympy.matrices.matrixbase.MatrixBase.rank : number of pivot columns.
    sympy.matrices.matrixbase.MatrixBase.nullspace : the null space, which
        is determined by the free (non-pivot) columns of the RREF.
    """
    # Try to use DomainMatrix for ZZ or QQ
    dM = _to_DM_ZZ_QQ(M)

    if dM is not None:
        # Use DomainMatrix for ZZ or QQ
        mat, pivot_cols = _rref_dm(dM)
    else:
        # Use the generic Matrix routine.
        if isinstance(simplify, FunctionType):
            simpfunc: Callable[[Any], Any] = simplify
        else:
            simpfunc = _simplify

        mat, pivot_cols, _ = _row_reduce(M, iszerofunc, simpfunc,
                normalize_last, normalize=True, zero_above=True)

    if pivots:
        return mat, pivot_cols
    else:
        return mat
