from sympy.matrices.expressions.matexpr  import MatrixExpr

class Block(MatrixExpr):
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
