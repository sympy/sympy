from __future__ import annotations
from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.core.basic import Basic
from sympy.core.containers import Tuple
from sympy.core.singleton import S
from sympy.core.sympify import sympify
from sympy.functions.elementary.integers import ceiling
from sympy.functions.elementary.miscellaneous import Max, Min
from sympy.functions.elementary.complexes import sign


def delta(value):
    """
    Indicator delta function: delta(x) = (sign(x)^2 - sign(x)) / 2.
    Returns 1 if x < 0, and 0 if x >= 0.
    """
    if hasattr(value, 'is_nonnegative') and value.is_nonnegative:
        return S.Zero
    if hasattr(value, 'is_negative') and value.is_negative:
        return S.One
    if hasattr(value, 'is_number') and value.is_number:
        if (value < 0) == S.true:
            return S.One
        return S.Zero
    return (sign(value)**2 - sign(value)) / 2


def decode_bound(bound, n):
    """
    Decodes stored canonical bound by adding n if bound is negative.
    """
    if (bound < 0) == S.true:
        return bound + n
    return bound


def sympify_slice(n, value):
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
        if (v < 0) == S.true:
            v += n
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

    is_negative_step = (
        (hasattr(step, 'is_negative') and step.is_negative is True) or
        (hasattr(step, 'is_number') and step.is_number and
         (step < 0) == S.true)
    )
    p = 1 if is_negative_step else 0

    is_concrete_n = hasattr(n, 'is_number') and n.is_number

    if start is None:
        start = p * (n - 1)
    else:
        if p == 1 and start == -n - 1:
            start = -1
        elif (start < 0) == S.true:
            if p == 1 and (start < -n) == S.true:
                start = -1
            else:
                start = Max(0, start + n)
        if is_concrete_n and hasattr(start, 'is_number') and start.is_number:
            start = Max(-p, Min(n - p, start))

    if stop is None:
        stop = -p if p == 1 else n
    else:
        if p == 1 and stop == -n - 1:
            stop = -1
        elif (stop < 0) == S.true:
            if (stop < -n) == S.true:
                stop = -1 if p == 1 else 0
            else:
                stop = stop + n
        if is_concrete_n and hasattr(stop, 'is_number') and stop.is_number:
            stop = Max(-p, Min(n - p, stop))

    if p == 1 and start == -1:
        start = -n - 1
    if p == 1 and stop == -1:
        stop = -n - 1

    return Tuple(start, stop, step)


normalize = sympify_slice


def slice_shape_1d(n, sl):
    start, stop, step = sl
    A = decode_bound(start, n)
    B = decode_bound(stop, n)
    diff = sympify(B - A)
    if step == 1:
        if diff == n:
            return n
        if hasattr(diff, 'is_number') and diff.is_number:
            return Max(0, diff)
        return diff
    if diff == n and (step > 0) == S.true:
        return ceiling(n / step)
    if hasattr(diff, 'is_number') and diff.is_number:
        return Max(0, ceiling(diff / step))
    return ceiling(diff / step)


def compose_slices(s1, s2, n):
    L1 = slice_shape_1d(n, s1)
    L2 = slice_shape_1d(L1, s2)
    A1 = decode_bound(s1[0], n)
    step1 = s1[2]
    A2 = decode_bound(s2[0], L1)
    step2 = s2[2]
    A_comp = A1 + A2 * step1
    step_comp = step1 * step2
    B_comp = A_comp + L2 * step_comp
    if (A_comp < 0) == S.true:
        A_comp = A_comp - n
    if (B_comp < 0) == S.true:
        B_comp = B_comp - n
    return Tuple(A_comp, B_comp, step_comp)


def mat_slice_of_slice(parent, rowslice, colslice):
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
