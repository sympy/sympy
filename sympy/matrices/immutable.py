from mpmath.matrices.matrices import _matrix

from sympy.core import Basic, Dict, Expr, Integer, Tuple
from sympy.core.cache import cacheit
from sympy.core.sympify import converter as sympify_converter, _sympify
from sympy.matrices.dense import DenseMatrix, DenseDomainMatrix
from sympy.matrices.expressions import MatrixExpr
from sympy.matrices.matrices import MatrixBase
from sympy.matrices.sparse import SparseMatrix
from sympy.multipledispatch import dispatch
from sympy.polys.domainmatrix import DomainMatrix


def sympify_matrix(arg):
    return arg.as_immutable()


sympify_converter[MatrixBase] = sympify_matrix


def sympify_mpmath_matrix(arg):
    mat = [_sympify(x) for x in arg]
    return ImmutableDenseMatrix(arg.rows, arg.cols, mat)


sympify_converter[_matrix] = sympify_mpmath_matrix



class ImmutableDenseMatrix(DenseMatrix, MatrixExpr):  # type: ignore
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
    _op_priority = 10.001

    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    __hash__ = MatrixExpr.__hash__

    @classmethod
    def _new(cls, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], ImmutableDenseMatrix):
            return args[0]
        if kwargs.get('copy', True) is False:
            if len(args) != 3:
                raise TypeError("'copy=False' requires a matrix be initialized as rows,cols,[list]")
            rows, cols, flat_list = args
        else:
            rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)
            flat_list = list(flat_list) # create a shallow copy

        if all(isinstance(x, Expr) and cls._can_use_dm(x) for x in flat_list):
            return ImmutableDenseDomainMatrix(rows, cols, flat_list)

        obj = Basic.__new__(cls,
            Integer(rows),
            Integer(cols),
            Tuple(*flat_list))
        obj._rows = rows
        obj._cols = cols
        obj._mat = flat_list
        return obj

    def _entry(self, i, j, **kwargs):
        return DenseMatrix.__getitem__(self, (i, j))

    def __setitem__(self, *args):
        raise TypeError("Cannot set values of {}".format(self.__class__))

    def _eval_extract(self, rowsList, colsList):
        # self._mat is a Tuple.  It is slightly faster to index a
        # tuple over a Tuple, so grab the internal tuple directly
        mat = self._mat
        cols = self.cols
        indices = (i * cols + j for i in rowsList for j in colsList)
        return self._new(len(rowsList), len(colsList),
                         Tuple(*(mat[i] for i in indices), sympify=False), copy=False)

    @property
    def cols(self):
        return self._cols

    @property
    def rows(self):
        return self._rows

    @property
    def shape(self):
        return self._rows, self._cols

    def as_immutable(self):
        return self

    def is_diagonalizable(self, reals_only=False, **kwargs):
        return super().is_diagonalizable(
            reals_only=reals_only, **kwargs)

    is_diagonalizable.__doc__ = DenseMatrix.is_diagonalizable.__doc__
    is_diagonalizable = cacheit(is_diagonalizable)


# make sure ImmutableDenseMatrix is aliased as ImmutableMatrix
ImmutableMatrix = ImmutableDenseMatrix


class ImmutableDenseDomainMatrix(DenseDomainMatrix, ImmutableDenseMatrix):

    __hash__ = MatrixExpr.__hash__

    def __new__(cls, *args):
        if len(args) == 1:
            return ImmutableDenseMatrix(*args)
        rows, cols, flat_list = args
        flat_list = [_sympify(e) for e in flat_list]
        rows_list = [[flat_list[i*cols + j] for j in range(cols)] for i in range(rows)]
        rep = DomainMatrix.from_list_sympy(rows, cols, rows_list)

        obj = Basic.__new__(cls,
            Integer(rows),
            Integer(cols),
            Tuple(*flat_list))

        obj._rows = int(rows)
        obj._cols = int(cols)
        obj._rep = rep
        return obj

    @classmethod
    def _new(cls, *args, **kwargs):
        # This will come back to DenseDomainMatrix.__new__ if possible
        return ImmutableDenseMatrix._new(*args, **kwargs)

    @classmethod
    def from_DomainMatrix(cls, rep):
        rows, cols = rep.shape
        flat_list = [rep.domain.to_sympy(e) for row in rep.rep for e in row]
        return ImmutableDenseDomainMatrix(rows, cols, flat_list)

    #@property
    #def args(self):
    #    return (self.rows, self.cols, tuple(self._flat()))
#
#    def _hashable_content(self):
#        return self.args

    def _entry(self, i, j, **kwargs):
        return DenseDomainMatrix.__getitem__(self, (i, j))


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
    is_Matrix = True
    _class_priority = 9

    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    __hash__ = MatrixExpr.__hash__

    @classmethod
    def _new(cls, *args, **kwargs):
        rows, cols, smat = cls._handle_creation_inputs(*args, **kwargs)
        #if all(x.is_Rational for x in smat.values()):
        #    flat_list = [0] * (rows * cols)
        #    for (i, j), v in smat.items():
        #        flat_list[i*cols + j] = v
        #    return MutableDenseDomainMatrix(rows, cols, flat_list)
        obj = Basic.__new__(cls, Integer(rows), Integer(cols), Dict(smat))
        obj._rows = rows
        obj._cols = cols
        obj._smat = smat
        return obj

    def __setitem__(self, *args):
        raise TypeError("Cannot set values of ImmutableSparseMatrix")

    def _entry(self, i, j, **kwargs):
        return SparseMatrix.__getitem__(self, (i, j))

    @property
    def cols(self):
        return self._cols

    @property
    def rows(self):
        return self._rows

    @property
    def shape(self):
        return self._rows, self._cols

    def as_immutable(self):
        return self

    def is_diagonalizable(self, reals_only=False, **kwargs):
        return super().is_diagonalizable(
            reals_only=reals_only, **kwargs)

    is_diagonalizable.__doc__ = SparseMatrix.is_diagonalizable.__doc__
    is_diagonalizable = cacheit(is_diagonalizable)


@dispatch(ImmutableDenseMatrix, ImmutableDenseMatrix)
def _eval_is_eq(lhs, rhs): # noqa:F811
    """Helper method for Equality with matrices.sympy.

    Relational automatically converts matrices to ImmutableDenseMatrix
    instances, so this method only applies here.  Returns True if the
    matrices are definitively the same, False if they are definitively
    different, and None if undetermined (e.g. if they contain Symbols).
    Returning None triggers default handling of Equalities.

    """
    if lhs.shape != rhs.shape:
        return False
    return (lhs - rhs).is_zero_matrix
