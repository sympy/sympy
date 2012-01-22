from matrices import Matrix
from expressions import MatrixExpr, Transpose
from sympy import Basic, Tuple

class ImmutableMatrix(MatrixExpr, Matrix):

    def __new__(cls, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], Matrix):
            mat = args[0]
        else:
            mat = Matrix(*args, **kwargs)
        shape = Tuple(*mat.shape)
        data = Tuple(*mat.mat)
        return Basic.__new__(cls, shape, data)

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

    def as_mutable(self):
        return Matrix(self.rows, self.cols, self.mat)

    def __setitem__(self, *args):
        raise TypeError('Can not set values in ImmutableMatrix')

    def __getitem__(self, *args):
        result = Matrix.__getitem__(self, *args)
        if result.is_Matrix:
            result = ImmutableMatrix(result)
        return result

    def transpose(self):
        return ImmutableMatrix(Matrix.transpose(self))

    def conjugate(self):
        return ImmutableMatrix(Matrix.conjugate(self))

    def expand(self, **hints):
        return ImmutableMatrix(Matrix.expand(self, **hints))

    def rref(self, *args, **kwargs):
        r, pivotlist = self.as_mutable().rref()
        return ImmutableMatrix(r), pivotlist

    def QRdecomposition(self):
        Q, R = self.as_mutable().QRdecomposition()
        return ImmutableMatrix(Q), ImmutableMatrix(R)

    def reshape(self, *args, **kwargs):
        return ImmutableMatrix(Matrix.reshape(self, *args, **kwargs))

