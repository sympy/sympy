from sympy.core.cache import cacheit
from matrices import MatrixBase
from dense import Matrix
from sparse import MutableSparseMatrix as SparseMatrix
from expressions import MatrixExpr
from sympy import Basic, Integer, Tuple, Dict

class ImmutableMatrix(MatrixExpr, Matrix):
    """Create an immutable version of a matrix.

    Examples
    ========

    >>> from sympy import eye
    >>> from sympy.matrices import ImmutableMatrix
    >>> ImmutableMatrix(eye(3))
    [1, 0, 0]
    [0, 1, 0]
    [0, 0, 1]
    >>> _[0, 0] = 42
    Traceback (most recent call last):
    ...
    TypeError: Cannot set values of ImmutableDenseMatrix
    """

    _class_priority = 8

    @classmethod
    def _new(cls, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], ImmutableMatrix):
            return args[0]
        rows, cols, flat_list = MatrixBase._handle_creation_inputs(*args, **kwargs)
        rows = Integer(rows)
        cols = Integer(cols)
        mat = Tuple(*flat_list)
        return Basic.__new__(cls, rows, cols, mat)

    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    @property
    def shape(self):
        return self.args[:2]

    @property
    def _mat(self):
        return self.args[2]

    def _entry(self, i, j):
        return Matrix.__getitem__(self, (i,j))

    __getitem__ = Matrix.__getitem__

    def __setitem__(self, *args):
        raise TypeError("Cannot set values of ImmutableMatrix")

    as_mutable = MatrixBase.as_mutable

    adjoint = MatrixBase.adjoint
    conjugate = MatrixBase.conjugate
    _eval_trace = Matrix._eval_trace
    _eval_transpose = Matrix._eval_transpose
    _eval_conjugate = Matrix._eval_conjugate
    _eval_inverse = Matrix._eval_inverse

    equals = Matrix.equals
    is_Identity = Matrix.is_Identity

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

class ImmutableSparseMatrix(Basic, SparseMatrix):
    """Create an immutable version of a sparse matrix.

    Examples
    ========

    >>> from sympy import eye
    >>> from sympy.matrices.immutable import ImmutableSparseMatrix
    >>> ImmutableSparseMatrix(1, 1, {})
    [0]
    >>> ImmutableSparseMatrix(eye(3))
    [1, 0, 0]
    [0, 1, 0]
    [0, 0, 1]
    >>> _[0, 0] = 42
    Traceback (most recent call last):
    ...
    TypeError: Cannot set values of ImmutableSparseMatrix
    >>> _.shape
    (3, 3)
    """

    _class_priority = 9

    @classmethod
    def _new(cls, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], cls):
            return args[0]
        rows, cols, mat = MatrixBase._handle_creation_inputs(*args, **kwargs)
        mat = Dict(SparseMatrix(mat)._smat)
        rows = Integer(rows)
        cols = Integer(cols)
        return Basic.__new__(cls, rows, cols, mat)

    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    def __setitem__(self, *args):
        raise TypeError("Cannot set values of ImmutableSparseMatrix")

    def _entry(self, i, j):
        return SparseMatrix.__getitem__(self, (i, j))

    __getitem__ = SparseMatrix.__getitem__
