from __future__ import print_function, division

from typing import Union

from sympy.core.containers import Tuple
from sympy.core.basic import Atom, Basic
from sympy.core.sympify import _sympify
from sympy.core.singleton import S, Singleton
from sympy.functions.elementary.integers import floor
from sympy.matrices.matrices import MatrixBase
from sympy.matrices.common import ShapeError

from .matexpr import MatrixExpr, MatrixSymbol, ZeroMatrix, OneMatrix, Identity


class SliceNone(Atom, metaclass=Singleton):
    """A corresponding sympy object for python singleton ``None`` used
    in ``slice``.

    Notes
    =====

    In some context of python array slicing, ``None`` cannot not be
    reduced to an integer.

    You cannot replace the second argument ``None`` with any integer
    to yield the same result, for following examples.

    >>> A = list(range(5))
    >>> A[slice(None, None, -1)]
    [4, 3, 2, 1, 0]
    >>> A[slice(None, None, -2)]
    [4, 2, 0]
    >>> A[slice(2, None, -1)]
    [2, 1, 0]
    >>> A[slice(2, None, -2)]
    [2, 0]

    This behavior is spcific to ``slice`` specification with
    ``step=None`` and ``step`` as a negative integer.
    """
    def __new__(cls):
        return super(SliceNone, cls).__new__(cls)


def _slice_sympify(i: Union[slice, tuple, list, Tuple]):
    if isinstance(i, slice):
        start, stop, step = i.start, i.stop, i.step
    else:
        if len(i) == 1:
            start, stop, step = i[0], S.SliceNone, S.SliceNone
        elif len(i) == 2:
            start, stop, step = i[0], i[1], S.SliceNone
        elif len(i) == 3:
            start, stop, step = i
        else:
            raise ValueError(
                "{} is not a valid slice specification because it has "
                "more than 3 elements.".format(i))

    def __sympify(item):
        if item is None:
            return S.SliceNone
        return _sympify(item)

    return __sympify(start), __sympify(stop), __sympify(step)


def _remove_slice_none(start, stop, step, parentsize):
    """Eliminate the arbitrary ``None``s from the slice specification
    if possible.
    """
    start_none = start is S.SliceNone
    stop_none = stop is S.SliceNone
    step_none = step is S.SliceNone

    if start_none and stop_none and step_none:
        return S.Zero, parentsize, S.One
    elif start_none and stop_none and not step_none:
        if step.is_positive:
            return S.Zero, parentsize, step
        elif step.is_negative:
            return parentsize-1, stop, step
    elif not start_none and stop_none and step_none:
        return start, parentsize, S.One
    elif start_none and not stop_none and step_none:
        return S.Zero, stop, S.One
    elif not start_none and not stop_none and step_none:
        return start, stop, S.One
    elif start_none and not stop_none and not step_none:
        return S.Zero, stop, step
    elif not start_none and stop_none and not step_none:
        if step.is_positive:
            return start, parentsize, step
    elif not start_none and not stop_none and step_none:
        return start, stop, S.One
    return start, stop, step


def normalize(i, parentsize: Basic):
    if not isinstance(i, (slice, tuple, list, Tuple)):
        if (i < -parentsize) | (i >= parentsize) == True:
            raise IndexError(
                "{} must be in the interval ({}, {}]."
                .format(i, -parentsize, parentsize))
        if (i < 0) == True:
            i += parentsize
        return i, i+1, 1

    start, stop, step = _slice_sympify(i)
    start, stop, step = _remove_slice_none(start, stop, step, parentsize)

    if (start <= -parentsize) == True:
        start = S.Zero
    elif (start > -parentsize) & (start < S.Zero) == True:
        start += parentsize
    elif start.is_zero:
        start = S.Zero
    elif (start > parentsize) == True:
        start = parentsize

    if (stop <= -parentsize) == True:
        stop = S.Zero
    elif (stop > -parentsize) & (stop < S.Zero) == True:
        stop += parentsize
    elif stop.is_zero:
        stop = S.Zero
    elif (stop > parentsize) == True:
        stop = parentsize

    return (start, stop, step)

class MatrixSlice(MatrixExpr):
    r"""A symbolic representation for slicing a matrix.

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

    How slicing works
    =================

    Slicing is one of the most popular feature of python arrays.
    However, to introduce some symbolic math or code generation
    features on top of it, we had to reinvent some wheels.

    Here, we take an example of slicing a one-dimensional vector,
    because extending this to any multidimensional array types can
    be done easily.

    Every specifications of python array slicing can be reduced either
    of these two cases.

    Case 1
    ======

    If *start*, *stop*, and *step* are defined explicitly and denoted
    by $i$, $j$, $k$ as:

    - $A$ is a vector with $n$ elements.
    - $i \in \mathbb{Z}$ denotes the *start*.
    - $j \in \mathbb{Z}$ denotes the *stop*.
    - $k \in \mathbb{Z} - \{0\}$ denotes the *step*.

    If $k > 0$, the indices $i$, $j$ can be normalized to
    $\bar{i}$, $\bar{j}$ as:

    .. math::
        \bar{i} = \begin{cases}
        0 & \text{ if } i \leq -n \\
        i+n & \text{ if } -n < i < 0 \\
        0 & \text{ if } i = 0 \\
        i & \text{ if } 0 < i < n
        \end{cases},
        \bar{j} = \begin{cases}
        0 & \text{ if } j \leq -n \\
        j+n & \text{ if } -n < j < 0 \\
        0 & \text{ if } j = 0 \\
        j & \text{ if } 0 < j < n \\
        n & \text{ if } j \geq n
        \end{cases}

    And every cases of $i \geq n$ can always return an empty vector.

    The slicing can be defined from the normalized indices as:

    - If $\bar{i} < \bar{j}$,
      $A[[\bar{i}:\bar{j}:k]] = [A[[\bar{i}]], A[[\bar{i}+k]], ...,
      A[[\bar{i} + (\lceil\frac{\bar{j}-\bar{i}}{k}\rceil-1) k]]]$
    - If $\bar{i} \geq \bar{j}$, $A[[\bar{i}:\bar{j}:k]] = []$

    If $k < 0$, the indices $i$, $j$ can be normalized to
    $\bar{i}$, $\bar{j}$ as:

    .. math::
        \bar{i} = \begin{cases}
        0 & \text{ if } i \leq -n \\
        i+n & \text{ if } -n < i < 0 \\
        0 & \text{ if } i = 0 \\
        i & \text{ if } 0 < i < n \\
        n-1 & \text{ if } i \geq n
        \end{cases},
        \bar{j} = \begin{cases}
        0 & \text{ if } j \leq -n \\
        j+n & \text{ if } -n < j < 0 \\
        0 & \text{ if } j = 0 \\
        j & \text{ if } 0 < j < n
        \end{cases}

    And every cases of $j \geq n$ can always return an empty vector.

    The slicing can be defined from the normalized indices as:

    - If $\bar{i} \leq \bar{j}$,
      $A[[\bar{i}:\bar{j}:k]] = []$
    - If $\bar{i} > \bar{j}$,
      $A[[\bar{i}:\bar{j}:k]] = [A[[\bar{i}]], A[[\bar{i}+k]], ...,
      A[[\bar{i} + (\lceil\frac{\bar{j}-\bar{i}}{k}\rceil-1) k]]]$

    Case 2
    ======

    This is a special case when *stop* is ``None`` and *step* is
    negative.

    - $A$ is a vector with $n$ elements.
    - $i \in \mathbb{Z}$ denotes the *start*.
    - $k \in \mathbb{Z}_{<0}$ denotes the *step*.

    The index $i$ can be normalized to $\bar{i}$ as:

    .. math::
        \bar{i} = \begin{cases}
        0 & \text{ if } i \leq -n \\
        i+n & \text{ if } -n < i < 0 \\
        0 & \text{ if } i = 0 \\
        i & \text{ if } 0 < i < n \\
        n-1 & \text{ if } i \geq n
        \end{cases}

    The slicing can be defined from the normalized indices as:

    $A[[\bar{i}::k]] = [A[[\bar{i}]], A[[\bar{i}+k]], ...,
    A[[\bar{i} + (\lceil\frac{-\bar{i}-1}{k}\rceil - 1) k]]]$

    How to comprehend None in slice
    ===============================

    Every python array slicing specifications with ``None`` can be
    reduced to one of the cases above because

    - $A[[::]] = A[[0:n:1]]$
    - $A[[::k]] = \begin{cases} A[[0:n:k]] & \text{ if } k > 0 \\
      A[[n-1::k]] & \text{ if } k < 0 \end{cases}$
    - $A[[:j:]] = A[[0:j:1]]$
    - $A[[i::]] = A[[i:n:1]]$
    - $A[[:j:k]] = A[[0:j:k]]$
    - $A[[i::k]] = \begin{cases} A[[i:n:k]] & \text{ if } k > 0 \\
      A[[i::k]] & \text{ if } k < 0 \end{cases}$
    - $A[[i:j:]] = A[[i:j:1]]$
    """
    parent = property(lambda self: self.args[0])
    rowslice = property(lambda self: self.args[1])
    colslice = property(lambda self: self.args[2])

    def __new__(cls, parent, rowslice, colslice):
        parent = _sympify(parent)
        if not parent.is_Matrix:
            raise ValueError("{} must be a matrix.".format(parent))

        rowslice = normalize(rowslice, parent.shape[0])
        colslice = normalize(colslice, parent.shape[1])

        if len(rowslice) != 3:
            raise IndexError(
                "{} is not a valid row slicing.".format(rowslice))
        if len(colslice) != 3:
            raise IndexError(
                "{} is not a valid column slicing.".format(colslice))

        if isinstance(parent, MatrixSlice):
            return mat_slice_of_slice(parent, rowslice, colslice)
        return Basic.__new__(cls, parent, Tuple(*rowslice), Tuple(*colslice))

    def _is_slice_numeric(self):
        def _pred(x):
            return x.is_Integer or x is S.SliceNone
        return all(_pred(x) for x in self.rowslice) and \
            all(_pred(x) for x in self.colslice)

    @property
    def shape(self):
        if not self._is_slice_numeric():
            raise ShapeError(
                "The matrix have an ambiguous slice specification to "
                "determine the shape.")

        rows = self.rowslice[1] - self.rowslice[0]
        rows = rows if self.rowslice[2] == 1 else floor(rows/self.rowslice[2])
        cols = self.colslice[1] - self.colslice[0]
        cols = cols if self.colslice[2] == 1 else floor(cols/self.colslice[2])
        return rows, cols

    def _entry(self, i, j, **kwargs):
        r_start, _, r_step = self.rowslice
        c_start, _, c_step = self.colslice

        return self.parent._entry(
            i*r_step + r_start, j*c_step + c_start, **kwargs)

    def doit(self, **kwargs):
        from sympy.combinatorics.permutations import Permutation
        from .blockmatrix import BlockMatrix, BlockDiagMatrix
        from .permutation import PermutationMatrix

        parent, rowslice, colslice = self.parent, self.rowslice, self.colslice

        are_slices_integers = all(x.is_Integer for x in rowslice.args) and \
            all(x.is_Integer for x in colslice.args)

        if not are_slices_integers:
            return self

        if rowslice == (0, parent.rows, 1) and \
            colslice == (0, parent.cols, 1):
            return parent

        if isinstance(parent, MatrixBase):
            rowslice = slice(*rowslice)
            colslice = slice(*colslice)
            return parent[rowslice, colslice]

        if isinstance(parent, MatrixSymbol):
            if rowslice[0] == 0 and rowslice[2] == 1 and \
                colslice[0] == 0 and colslice[2] == 1:
                return MatrixSymbol(parent.args[0], rowslice[1], colslice[1])

        if isinstance(parent, ZeroMatrix):
            return ZeroMatrix(*self.shape)

        if isinstance(parent, OneMatrix):
            return OneMatrix(*self.shape)

        if parent.is_Identity and self.on_diag:
            return Identity(self.rows)

        if isinstance(parent, BlockDiagMatrix):
            r_start, r_stop, r_step = rowslice
            c_start, c_stop, c_step = colslice
            if not r_step == c_step == 1:
                return self

            if not parent.rows.is_Integer or not parent.cols.is_Integer:
                return self

            cumulative_rows = parent._cumulative_rows()
            cumulative_cols = parent._cumulative_cols()
            blocks = parent.args

            if not r_start in cumulative_rows:
                return self
            if not r_stop in cumulative_rows:
                return self
            if not c_start in cumulative_cols:
                return self
            if not c_stop in cumulative_cols:
                return self

            r_start = cumulative_rows.index(r_start)
            r_stop = cumulative_rows.index(r_stop)
            c_start = cumulative_cols.index(c_start)
            c_stop = cumulative_cols.index(c_stop)
            if r_start == c_start and r_stop == c_stop:
                return BlockDiagMatrix(*blocks[r_start:r_stop])

        if isinstance(parent, BlockMatrix):
            r_start, r_stop, r_step = rowslice
            c_start, c_stop, c_step = colslice
            if not r_step == c_step == 1:
                return self

            if not parent.rows.is_Integer or not parent.cols.is_Integer:
                return self

            cumulative_rows = parent._cumulative_rows()
            cumulative_cols = parent._cumulative_cols()
            blocks = parent.blocks

            if not r_start in cumulative_rows:
                return self
            if not r_stop in cumulative_rows:
                return self
            if not c_start in cumulative_cols:
                return self
            if not c_stop in cumulative_cols:
                return self

            r_start = cumulative_rows.index(r_start)
            r_stop = cumulative_rows.index(r_stop)
            c_start = cumulative_cols.index(c_start)
            c_stop = cumulative_cols.index(c_stop)
            return BlockMatrix(blocks[r_start:r_stop, c_start:c_stop])

        if isinstance(parent, PermutationMatrix):
            if not self.on_diag:
                return self

            r_start, r_stop, r_step = rowslice
            c_start, c_stop, c_step = colslice

            if not r_step == c_step == 1:
                return self

            # Block decomposition of permutation matrix
            perm = parent.args[0]
            block_cyclic_form = perm._block_cyclic_form()

            cumulative_size = [0]
            for block in block_cyclic_form:
                cumulative_size += [
                    cumulative_size[-1] + sum(len(cycle) for cycle in block)]

            if not r_start in cumulative_size:
                return self
            if not r_stop in cumulative_size:
                return self

            block_r_start = cumulative_size.index(r_start)
            block_r_stop = cumulative_size.index(r_stop)

            norm = r_start
            new_cycles = []
            for block in block_cyclic_form[block_r_start:block_r_stop]:
                for cycle in block:
                    new_cycles.append([x-norm for x in cycle])

            return PermutationMatrix(Permutation(new_cycles))

        return self

    @property
    def on_diag(self):
        return self.rowslice == self.colslice

    def _eval_rewrite_as_MatrixPermute(self, *args, **kwargs):
        """Rewrite horizontal or vertical matrix flip as matrix
        permutation.
        """
        from sympy.combinatorics.permutations import Permutation
        from .permutation import MatrixPermute

        parent = self.parent
        r_start, r_stop, r_step = self.rowslice
        c_start, c_stop, c_step = self.colslice

        is_h_id = c_step == 1 and c_start == 0 and c_stop == parent.cols
        is_v_id = r_step == 1 and r_start == 0 and r_stop == parent.rows
        is_h_flip = \
            c_step == -1 and c_start == parent.cols and c_stop == 0
        is_v_flip = \
            r_step == -1 and r_start == parent.rows and r_stop == 0

        if is_h_flip and is_v_flip:
            af_rows = list(reversed(range(parent.rows)))
            af_cols = list(reversed(range(parent.cols)))
            return MatrixPermute(
                MatrixPermute(self.parent, Permutation(af_rows), 0),
                Permutation(af_cols), 1)
        elif is_h_flip and is_v_id:
            af_cols = list(reversed(range(parent.cols)))
            return MatrixPermute(self.parent, Permutation(af_cols), 1)
        elif is_v_flip and is_h_id:
            af_rows = list(reversed(range(parent.rows)))
            return MatrixPermute(self.parent, Permutation(af_rows), 0)


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
    X[3:7:2, 4]
    """
    row = slice_of_slice(parent.rowslice, rowslice)
    col = slice_of_slice(parent.colslice, colslice)
    return MatrixSlice(parent.parent, row, col)
