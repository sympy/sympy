from __future__ import annotations
from typing import TYPE_CHECKING
from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.core.basic import Basic
from sympy.core.containers import Tuple
from sympy.core.singleton import S
from sympy.core.sympify import sympify
from sympy.functions.elementary.integers import ceiling
from sympy.functions.elementary.miscellaneous import Max, Min
from sympy.functions.elementary.complexes import sign

if TYPE_CHECKING:
    from sympy.core.expr import Expr


def delta(value: Expr) -> Expr:
    """
    Indicator delta function: delta(x) = (sign(x)^2 - sign(x)) / 2.
    Returns 1 if x < 0, and 0 if x >= 0.
    """
    return (sign(value)**2 - sign(value)) / 2


def decode_bound(bound: Expr, n: Expr) -> Expr:
    """
    Decodes stored canonical bound by adding n if bound is negative.
    """
    return bound + n*delta(bound)


def sympify_slice(n: Expr, value: slice | tuple | Tuple | Expr) -> Tuple:
    """
    Canonicalizes slice specification with respect to axis dimension n.
    Implements the complete symbolic formula for array slicing.
    """
    n = sympify(n)
    if isinstance(value, slice):
        value = (value.start, value.stop, value.step)
    if not isinstance(value, (tuple, list, Tuple)):
        v = sympify(value)
        if (v >= n) == S.true or (v < -n) == S.true:
            raise IndexError(
                f"Index {value} out of bounds for axis of size {n}"
            )
        v = v + n * delta(v)
        value = (v, v + 1, 1)
    if len(value) == 2:
        value = (*value, None)
    if len(value) != 3:
        raise TypeError("a slice tuple needs start, stop, and step")
    start = None if value[0] is None else sympify(value[0])
    stop = None if value[1] is None else sympify(value[1])
    step = S.One if value[2] is None else sympify(value[2])
    if step.is_zero is True:
        raise ValueError("slice step cannot be zero")

    p = delta(step)
    if start is None:
        start = p*(n - 1)
    else:
        start = start + n*delta(start)
        start = Max(-p, Min(n - p, start))
        start = start - n*delta(start)
    if stop is None:
        stop = (1 - 2*p)*n - p
    else:
        stop = stop + n*delta(stop)
        stop = Max(-p, Min(n - p, stop))
        stop = stop - n*delta(stop)
    return Tuple(start, stop, step)


normalize = sympify_slice


def slice_shape_1d(n: Expr, sl: tuple[Expr, Expr, Expr] | Tuple) -> Expr:
    start, stop, step = sl
    A = decode_bound(start, n)
    B = decode_bound(stop, n)
    return Max(0, ceiling((B - A) / step))


def compose_slices(
    s1: tuple[Expr, Expr, Expr] | Tuple,
    s2: tuple[Expr, Expr, Expr] | Tuple,
    n: Expr,
) -> Tuple:
    L1 = slice_shape_1d(n, s1)
    L2 = slice_shape_1d(L1, s2)
    A1 = decode_bound(s1[0], n)
    step1 = s1[2]
    A2 = decode_bound(s2[0], L1)
    step2 = s2[2]
    A_comp = A1 + A2 * step1
    step_comp = step1 * step2
    B_comp = A_comp + L2 * step_comp
    A_comp = A_comp - n * delta(A_comp)
    B_comp = B_comp - n * delta(B_comp)
    return Tuple(A_comp, B_comp, step_comp)


def mat_slice_of_slice(
    parent: MatrixSlice,
    rowslice: Tuple,
    colslice: Tuple,
) -> MatrixSlice:
    """ Collapse nested matrix slices """
    row = compose_slices(parent.rowslice, rowslice, parent.parent.rows)
    col = compose_slices(parent.colslice, colslice, parent.parent.cols)
    return MatrixSlice(parent.parent, row, col)


class MatrixSlice(MatrixExpr):
    """ A MatrixSlice of a Matrix Expression

    Examples
    ========

    >>> from sympy import MatrixSlice, ImmutableMatrix
    >>> M = ImmutableMatrix(4, 4, range(16))
    >>> M
    Matrix([
    [ 0,  1,  2,  3],
    [ 4,  5,  6,  7],
    [ 8,  9, 10, 11],
    [12, 13, 14, 15]])

    >>> B = MatrixSlice(M, (0, 2), (2, 4))
    >>> ImmutableMatrix(B)
    Matrix([
    [2, 3],
    [6, 7]])
    """
    parent = property(lambda self: self.args[0])
    rowslice = property(lambda self: self.args[1])
    colslice = property(lambda self: self.args[2])

    def __new__(cls, parent, rowslice, colslice):
        parent = sympify(parent)
        if not isinstance(parent, MatrixExpr):
            raise TypeError("parent must be a MatrixExpr")
        rowslice = sympify_slice(parent.rows, rowslice)
        colslice = sympify_slice(parent.cols, colslice)
        if isinstance(parent, MatrixSlice):
            return mat_slice_of_slice(parent, rowslice, colslice)
        return Basic.__new__(cls, parent, Tuple(*rowslice), Tuple(*colslice))

    @property
    def shape(self):
        rows = slice_shape_1d(self.parent.rows, self.rowslice)
        cols = slice_shape_1d(self.parent.cols, self.colslice)
        return rows, cols

    def _entry(self, i, j, **kwargs):
        row_start, _, row_step = self.rowslice
        col_start, _, col_step = self.colslice
        A_row = decode_bound(row_start, self.parent.rows)
        A_col = decode_bound(col_start, self.parent.cols)
        return self.parent._entry(
            A_row + i * row_step,
            A_col + j * col_step,
            **kwargs,
        )

    @property
    def on_diag(self):
        return self.rowslice == self.colslice
