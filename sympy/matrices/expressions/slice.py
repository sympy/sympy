from __future__ import annotations
from sympy.matrices.expressions.matexpr import MatrixExpr
from sympy.core.basic import Basic
from sympy.core.containers import Tuple
from sympy.core.numbers import Integer
from sympy.functions.elementary.integers import ceiling


def _is_concrete_index(x):
    return isinstance(x, int) or (isinstance(x, Basic) and x.is_Integer)


def _slice_length(sl):
    start, stop, step = sl
    if all(_is_concrete_index(x) for x in (start, stop, step)):
        return Integer(len(range(int(start), int(stop), int(step))))
    length = stop - start
    return length if step == 1 else ceiling(length / step)


def normalize(i, parentsize):
    if isinstance(i, slice):
        i = (i.start, i.stop, i.step)
    if not isinstance(i, (tuple, list, Tuple)):
        # A single index is not clamped: an out-of-range index is an error,
        # just like indexing into a Python sequence.
        if (i < 0) == True:
            i += parentsize
        return (i, i + 1, 1)
    i = list(i)
    if len(i) == 2:
        i.append(1)
    start, stop, step = i
    # When the slice bounds and the parent dimension are concrete integers,
    # reuse Python's slice semantics (``slice.indices``) so that slicing a
    # symbolic matrix agrees with slicing an explicit one, including negative
    # steps and out-of-range bounds which are clamped rather than raising
    # (see issue #18411).
    if _is_concrete_index(parentsize) and all(
        x is None or _is_concrete_index(x) for x in (start, stop, step)
    ):
        return slice(
            None if start is None else int(start),
            None if stop is None else int(stop),
            None if step is None else int(step),
        ).indices(int(parentsize))
    start = start or 0
    if stop is None:
        stop = parentsize
    if (start < 0) == True:
        start += parentsize
    if (stop < 0) == True:
        stop += parentsize
    step = step or 1

    if ((stop - start) * step < 1) == True:
        raise IndexError()

    return (start, stop, step)

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
        rowslice = normalize(rowslice, parent.shape[0])
        colslice = normalize(colslice, parent.shape[1])
        if not (len(rowslice) == len(colslice) == 3):
            raise IndexError()
        if ((0 > rowslice[0]) == True or
            (parent.shape[0] < rowslice[1]) == True or
            (0 > colslice[0]) == True or
            (parent.shape[1] < colslice[1]) == True):
            raise IndexError()
        if isinstance(parent, MatrixSlice):
            return mat_slice_of_slice(parent, rowslice, colslice)
        return Basic.__new__(cls, parent, Tuple(*rowslice), Tuple(*colslice))

    @property
    def shape(self):
        return _slice_length(self.rowslice), _slice_length(self.colslice)

    def _entry(self, i, j, **kwargs):
        return self.parent._entry(i*self.rowslice[2] + self.rowslice[0],
                                  j*self.colslice[2] + self.colslice[0],
                                  **kwargs)

    @property
    def on_diag(self):
        return self.rowslice == self.colslice


def slice_of_slice(s, t):
    start1, stop1, step1 = s
    start2, stop2, step2 = t

    start = start1 + start2*step1
    step = step1 * step2
    stop = start1 + step1*stop2

    if stop > stop1:
        raise IndexError()

    return start, stop, step


def mat_slice_of_slice(parent, rowslice, colslice):
    """ Collapse nested matrix slices

    >>> from sympy import MatrixSymbol
    >>> X = MatrixSymbol('X', 10, 10)
    >>> X[:, 1:5][5:8, :]
    X[5:8, 1:5]
    >>> X[1:9:2, 2:6][1:3, 2]
    X[3:7:2, 4:5]
    """
    row = slice_of_slice(parent.rowslice, rowslice)
    col = slice_of_slice(parent.colslice, colslice)
    return MatrixSlice(parent.parent, row, col)
