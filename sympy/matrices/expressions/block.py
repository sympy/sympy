from sympy.matrices.expressions.matexpr  import MatrixExpr

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

    @property
    def shape(self):
        return (self.rowbounds[1] - self.rowbounds[0],
                self.colbounds[1] - self.colbounds[0])

    def _entry(self, i, j):
        return self.parent._entry(i + self.rowbounds[0], j + self.colbounds[0])

    @property
    def on_diag(self):
        return self.rowbounds == self.colbounds
