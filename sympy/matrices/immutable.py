from __future__ import print_function, division

from sympy.core import Basic, Integer, Tuple, Dict, S, sympify
from sympy.core.sympify import converter as sympify_converter

from sympy.matrices.matrices import MatrixBase
from sympy.matrices.dense import DenseMatrix
from sympy.matrices.sparse import SparseMatrix, MutableSparseMatrix
from sympy.matrices.expressions import MatrixExpr


def sympify_matrix(arg):
    return arg.as_immutable()
sympify_converter[MatrixBase] = sympify_matrix

class ImmutableMatrix(MatrixExpr, DenseMatrix):
    """Create an immutable version of a matrix.

    Examples
    ========

    >>> from sympy import eye
    >>> from sympy.matrices import ImmutableMatrix
    >>> ImmutableMatrix(eye(3))
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])
    >>> _[0, 0] = 42
    Traceback (most recent call last):
    ...
    TypeError: Cannot set values of ImmutableDenseMatrix
    """

    # MatrixExpr is set as NotIterable, but we want explicit matrices to be
    # iterable
    _iterable = True

    _class_priority = 8

    @classmethod
    def _new(cls, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], ImmutableMatrix):
            return args[0]
        rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)
        rows = Integer(rows)
        cols = Integer(cols)
        mat = Tuple(*flat_list)
        return Basic.__new__(cls, rows, cols, mat)

    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    @property
    def shape(self):
        return tuple([int(i) for i in self.args[:2]])

    @property
    def _mat(self):
        return list(self.args[2])

    def _entry(self, i, j):
        return DenseMatrix.__getitem__(self, (i, j))

    __getitem__ = DenseMatrix.__getitem__

    def __setitem__(self, *args):
        raise TypeError("Cannot set values of ImmutableMatrix")

    def _eval_Eq(self, other):
        """Helper method for Equality with matrices.

        Relational automatically converts matrices to ImmutableMatrix
        instances, so this method only applies here.  Returns True if the
        matrices are definitively the same, False if they are definitively
        different, and None if undetermined (e.g. if they contain Symbols).
        Returning None triggers default handling of Equalities.

        """
        if not hasattr(other, 'shape') or self.shape != other.shape:
            return S.false
        if isinstance(other, MatrixExpr) and not isinstance(
                other, ImmutableMatrix):
            return None
        diff = self - other
        return sympify(diff.is_zero)

    adjoint = MatrixBase.adjoint
    conjugate = MatrixBase.conjugate
    # C and T are defined in MatrixExpr...I don't know why C alone
    # needs to be defined here
    C = MatrixBase.C

    as_mutable = DenseMatrix.as_mutable
    _eval_trace = DenseMatrix._eval_trace
    _eval_transpose = DenseMatrix._eval_transpose
    _eval_conjugate = DenseMatrix._eval_conjugate
    _eval_adjoint = DenseMatrix._eval_adjoint
    _eval_inverse = DenseMatrix._eval_inverse
    _eval_simplify = DenseMatrix._eval_simplify

    equals = DenseMatrix.equals
    is_Identity = DenseMatrix.is_Identity

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
# This is included after the class definition as a workaround for issue 7213.
# See https://github.com/sympy/sympy/issues/7213
ImmutableMatrix.is_zero = DenseMatrix.is_zero


class ImmutableSparseMatrix(Basic, SparseMatrix):
    """Create an immutable version of a sparse matrix.

    Examples
    ========

    >>> from sympy import eye
    >>> from sympy.matrices.immutable import ImmutableSparseMatrix
    >>> ImmutableSparseMatrix(1, 1, {})
    Matrix([[0]])
    >>> ImmutableSparseMatrix(eye(3))
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])
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
        s = MutableSparseMatrix(*args)
        rows = Integer(s.rows)
        cols = Integer(s.cols)
        mat = Dict(s._smat)
        obj = Basic.__new__(cls, rows, cols, mat)
        obj.rows = s.rows
        obj.cols = s.cols
        obj._smat = s._smat
        return obj

    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    def __setitem__(self, *args):
        raise TypeError("Cannot set values of ImmutableSparseMatrix")

    subs = MatrixBase.subs

    def __hash__(self):
        return hash((type(self).__name__,) + (self.shape, tuple(self._smat)))

    _eval_Eq = ImmutableMatrix._eval_Eq
