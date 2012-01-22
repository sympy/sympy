from matrices import MatrixBase
from expressions import MatrixExpr, Transpose
from sympy import Basic, Tuple

class ImmutableMatrix(MatrixExpr, MatrixBase):
    def __new__(cls, *args, **kwargs):
        rows, cols, mat = MatrixBase._handle_creation_inputs(*args, **kwargs)
        shape = Tuple(rows, cols)
        mat = Tuple(*mat)
        return Basic.__new__(cls, shape, mat)

    @property
    def shape(self):
        return self.args[0]

    @property
    def mat(self):
        return self.args[1]

    @property
    def rows(self):
        return self.shape[0]

    @property
    def cols(self):
        return self.shape[1]

    def _entry(self, i, j):
        return MatrixBase.__getitem__(self, (i,j))
