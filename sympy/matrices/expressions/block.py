from sympy.matrices.expressions.matexpr  import MatrixExpr
from sympy import Tuple, Basic
from sympy.functions.elementary.integers import floor

def normalize(i, parentsize):
    if isinstance(i, slice):
        i = (i.start, i.stop, i.step)
    if not isinstance(i, (tuple, list, Tuple)):
        return (i, i+1, 1)
    i = list(i)
    if len(i) == 2:
        i.append(1)
    if i[0] == None:
        i[0] = 0
    if i[1] == None:
        i[1] = parentsize
    if i[2] == None:
        i[2] = 1
    return tuple(i)

class Block(MatrixExpr):
    """ A block of a Matrix Expression

    Examples

    >>> from sympy import Block, ImmutableMatrix
    >>> M = ImmutableMatrix(4, 4, range(16))
    >>> print M
    [ 0,  1,  2,  3]
    [ 4,  5,  6,  7]
    [ 8,  9, 10, 11]
    [12, 13, 14, 15]

    >>> B = Block(M, (0, 2), (2, 4))
    >>> print ImmutableMatrix(B)
    [2, 3]
    [6, 7]
    """
    parent = property(lambda self: self.args[0])
    rowslice = property(lambda self: self.args[1])
    colslice = property(lambda self: self.args[2])

    def __new__(cls, parent, rowslice, colslice):
        rowslice = normalize(rowslice, parent.shape[0])
        colslice = normalize(colslice, parent.shape[1])
        return Basic.__new__(cls, parent, Tuple(*rowslice), Tuple(*colslice))

    @property
    def shape(self):
        rows = self.rowslice[1] - self.rowslice[0]
        rows = rows if self.rowslice[2] == 1 else floor(rows/self.rowslice[2])
        cols = self.colslice[1] - self.colslice[0]
        cols = cols if self.colslice[2] == 1 else floor(cols/self.colslice[2])
        return rows, cols

    def _entry(self, i, j):
        return self.parent._entry(i*self.rowslice[2] + self.rowslice[0],
                                  j*self.colslice[2] + self.colslice[0])

    @property
    def on_diag(self):
        return self.rowslice == self.colslice
