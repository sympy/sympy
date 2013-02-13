from sympy.matrices.expressions.matexpr  import MatrixExpr
from sympy import Tuple, Basic

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
    rowbounds = property(lambda self: self.args[1])
    colbounds = property(lambda self: self.args[2])

    def __new__(cls, parent, rowbounds, colbounds):
        if not isinstance(rowbounds, (tuple, list, Tuple)):
            rowbounds = rowbounds, rowbounds + 1
        if not isinstance(colbounds, (tuple, list, Tuple)):
            colbounds = colbounds, colbounds + 1
        return Basic.__new__(cls, parent, Tuple(*rowbounds), Tuple(*colbounds))

    @property
    def shape(self):
        return (self.rowbounds[1] - self.rowbounds[0],
                self.colbounds[1] - self.colbounds[0])

    def _entry(self, i, j):
        return self.parent._entry(i + self.rowbounds[0], j + self.colbounds[0])

    @property
    def on_diag(self):
        return self.rowbounds == self.colbounds
