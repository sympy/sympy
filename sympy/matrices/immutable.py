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

class ImmutableMatrix(DenseMatrix, MatrixExpr):
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

    _op_priority = 10.001
    _class_priority = 8

    def __new__(cls, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], ImmutableMatrix):
            return args[0]
        if not kwargs.get('copy', True):
            if len(args) != 3:
                raise TypeError("'copy' requires a matrix be initialized as rows,cols,[list]")
            rows, cols, flat_list = args
        else:
            rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)
        rows = Integer(rows)
        cols = Integer(cols)
        mat = Tuple(*flat_list)
        return Basic.__new__(cls, rows, cols, mat)

    __hash__ = MatrixExpr.__hash__

    @classmethod
    def _new(cls, *args, **kwargs):
        return cls(*args, **kwargs)

    def __setitem__(self, *args):
        raise TypeError("Cannot set values of ImmutableMatrix")

    def _entry(self, i, j):
        return DenseMatrix.__getitem__(self, (i, j))

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

    @property
    def _mat(self):
        return list(self.args[2])

    @property
    def shape(self):
        return tuple(int(i) for i in self.args[:2])
# This is included after the class definition as a workaround for issue 7213.
# if is_commutative == False, the assumptions system will declare
# the object is non-zero
# See https://github.com/sympy/sympy/issues/7213
ImmutableMatrix.is_zero = DenseMatrix.is_zero


class ImmutableSparseMatrix(SparseMatrix, MatrixExpr):
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

    def __new__(cls, *args, **kwargs):
        s = MutableSparseMatrix(*args, **kwargs)
        rows = Integer(s.rows)
        cols = Integer(s.cols)
        mat = Dict(s._smat)
        return Basic.__new__(cls, rows, cols, mat)

    __hash__ = MatrixExpr.__hash__

    @classmethod
    def _new(cls, *args, **kwargs):
        return cls(*args, **kwargs)

    @property
    def _smat(self):
        return {tuple(key): val for key,val in self.args[2].items()}

    def __setitem__(self, *args):
        raise TypeError("Cannot set values of ImmutableSparseMatrix")

    @property
    def shape(self):
        return tuple(int(i) for i in self.args[:2])
# This is included after the class definition as a workaround for issue 7213.
# if is_commutative == False, the assumptions system will declare
# the object is non-zero
# See https://github.com/sympy/sympy/issues/7213
ImmutableSparseMatrix.is_zero = SparseMatrix.is_zero

# these must be set after the classes are defined and one must be careful to
# avoid circular imports
DenseMatrix._immutable_variant = ImmutableMatrix
SparseMatrix._immutable_variant = ImmutableSparseMatrix
