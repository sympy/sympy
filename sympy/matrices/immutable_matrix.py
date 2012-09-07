from matrices import MatrixBase, MutableMatrix
from expressions import MatrixExpr, Transpose
from sympy import Basic, Integer, Tuple

class ImmutableMatrix(MatrixExpr, MatrixBase):

    _class_priority = 8

    @classmethod
    def _new(cls, *args, **kwargs):
        if len(args)==1 and isinstance(args[0], ImmutableMatrix):
            return args[0]
        rows, cols, mat = MatrixBase._handle_creation_inputs(*args, **kwargs)
        rows = Integer(rows)
        cols = Integer(cols)
        mat = Tuple(*mat)
        return Basic.__new__(cls, rows, cols, mat)
    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    @property
    def shape(self):
        return self.args[0:2]

    @property
    def mat(self):
        return self.args[2]

    def _entry(self, i, j):
        return MatrixBase.__getitem__(self, (i,j))

    def __setitem__(self, *args):
        raise TypeError("Can not set values in Immutable Matrix. "
                        "Use Matrix instead.")

    __getitem__ = MatrixBase.__getitem__

    as_mutable = MatrixBase.as_mutable

    adjoint = MatrixBase.adjoint
    conjugate = MatrixBase.conjugate
    equals = MatrixBase.equals
    is_Identity = MatrixBase.is_Identity
    transpose = MatrixBase.transpose

    __add__ = MatrixBase.__add__
    __radd__ = MatrixBase.__radd__
    __mul__ = MatrixBase.__mul__
    __rmul__ = MatrixBase.__rmul__
    __pow__ = MatrixBase.__pow__
    __sub__ = MatrixBase.__sub__
    __rsub__ = MatrixBase.__rsub__
    __neg__ = MatrixBase.__neg__
    __div__ = MatrixBase.__div__
    __truediv__ = MatrixBase.__truediv__
