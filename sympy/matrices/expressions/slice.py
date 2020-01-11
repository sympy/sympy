from __future__ import print_function, division

from sympy.core.containers import Tuple
from sympy.core.basic import Basic
from sympy.core.sympify import _sympify
from sympy.functions.elementary.integers import floor
from sympy.matrices.matrices import MatrixBase

from .matexpr import MatrixExpr, MatrixSymbol, ZeroMatrix, OneMatrix, Identity



def normalize(i, parentsize):
    if isinstance(i, slice):
        i = (i.start, i.stop, i.step)
    if not isinstance(i, (tuple, list, Tuple)):
        if (i < 0) == True:
            i += parentsize
        i = (i, i+1, 1)
    i = list(i)
    if len(i) == 2:
        i.append(1)
    start, stop, step = i
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

        if (rowslice[0] < 0) == True:
            raise IndexError("{} should not be negative.".format(rowslice[0]))
        if (colslice[0] < 0) == True:
            raise IndexError("{} should not be negative.".format(colslice[0]))
        if (parent.shape[0] < rowslice[1]) == True:
            raise IndexError(
                "{} should not exceed {}."
                .format(rowslice[1], parent.shape[0]))
        if (parent.shape[1] < colslice[1]) == True:
            raise IndexError(
                "{} should not exceed {}."
                .format(colslice[1], parent.shape[1]))

        if isinstance(parent, MatrixSlice):
            return mat_slice_of_slice(parent, rowslice, colslice)
        return Basic.__new__(cls, parent, Tuple(*rowslice), Tuple(*colslice))

    @property
    def shape(self):
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
