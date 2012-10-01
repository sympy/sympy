from matrices import MatrixBase
from mutable import MutableMatrix
from expressions import MatrixExpr
from sympy import Basic, Integer, Tuple

class ImmutableMatrix(MatrixExpr, MutableMatrix):

    _class_priority = 8

    @classmethod
    def _new(cls, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], ImmutableMatrix):
            return args[0]
        rows, cols, mat = MatrixBase._handle_creation_inputs(*args, **kwargs)
        rows = Integer(rows)
        cols = Integer(cols)
        mat = Tuple(*mat)
        return Basic.__new__(cls, rows, cols, mat)

    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    @property
    def _mat(self):
        return self.args[2]

    def _entry(self, i, j):
        return MutableMatrix.__getitem__(self, (i,j))

    __getitem__ = MutableMatrix.__getitem__

    def __setitem__(self, *args):
        raise TypeError("Cannot set values of ImmutableMatrix")

    as_mutable = MatrixBase.as_mutable

    adjoint = MatrixBase.adjoint
    conjugate = MatrixBase.conjugate
    equals = MatrixBase.equals
    is_Identity = MatrixBase.is_Identity
    _eval_trace = MutableMatrix._eval_trace
    _eval_transpose = MutableMatrix._eval_transpose
    _eval_conjugate = MutableMatrix._eval_conjugate
    _eval_inverse = MutableMatrix._eval_inverse

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
