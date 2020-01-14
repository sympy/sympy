from __future__ import division, print_function

from mpmath import mpf, mpc

from sympy import Integer, Float, Rational
from sympy.core.decorators import call_highest_priority
from sympy.core.sympify import sympify
from sympy.simplify import nsimplify

from .common import ShapeError, NonSquareMatrixError
from .dense import MutableDenseMatrix


def _is_integer(v):
    if isinstance(v, (int, Integer)):
        return True
    if isinstance(v, str):
        return False

    try:
        if int(v) == v:
            return True
    except (TypeError, ValueError):
        pass

    return False

def _is_scalar(v):
    if isinstance(v, (int, float, complex, mpf, mpc, Integer, Float, Rational)):
        return True
    if isinstance(v, str):
        return False

    try:
        complex(v)
    except (TypeError, ValueError):
        return False

    return True

def _detect_dtype(l): # Is there a better way to do this?
    import numpy as np

    try:
        [float(e) for e in l]
    except (TypeError, ValueError): # complex or symbolic source
        pass

    else:
        try:
            i = [int(e) for e in l]
        except (TypeError, ValueError): # have infinities and/or nans
            return np.float64

        else:
            if any(isinstance(e, Float) for e in l) or \
                    any(i != e for i, e in zip(i, l)):
                return np.float64 # this also catches integers which are too big
            else:
                return np.int64

    try:
        [complex(e) for e in l]
    except (TypeError, ValueError): # non-numeric source data
        return None

    return np.complex128


class NumPyMatrix(MutableDenseMatrix):
    """Matrix class which can wrap or create a new NumPy matrix (ndarray or
    matrix) for use in SymPy. Dispatches supported operations for faster
    execution to NumPy."""

    _op_priority          = 10.0005 # above Expr but below other matrix types
    _class_priority       = 2 # above _MatrixWrapper but below other matrix types
    _NDArrayFlatSympified = None

    # Wraps creation of a single NDArrayFlatSympified class to defer importing
    # numpy until it is needed.
    @property
    def NDArrayFlatSympified(self):
        if NumPyMatrix._NDArrayFlatSympified:
            return NumPyMatrix._NDArrayFlatSympified

        import numpy as np

        class NDArrayFlatSympified(np.ndarray):
            """Creates a view on an existing ``numpy.ndarray`` which sympifies
            all element reads so it can be used anywhere a sympified list of
            elements is required, such as ``M._mat``."""

            def __new__(cls, src):
                return np.ndarray.__new__(cls, src.size, dtype=src.dtype,
                        buffer=src)

            def __getitem__(self, key):
                return sympify(np.ndarray.__getitem__(self, key))

        NumPyMatrix._NDArrayFlatSympified = NDArrayFlatSympified

        return NumPyMatrix._NDArrayFlatSympified

    @classmethod
    def _new(cls, *args, dtype=None, copy=True, **kwargs):
        """Create a new NumPy store backed numeric matrix with either a new
        NumPy ``ndarray`` for the data or by creating a view on an existing
        ``ndarray``. If a view is being created on an existing ``ndarray`` then
        the dimensions can be reshaped by passing the new ``rows`` and ``cols``
        as long as the total number of elements remains the same.

        Parameters
        ==========

        dtype : numpy.dtype, optional
            This can be either a NumPy character dtype specification like
            ``'f4'`` or a class like ``numpy.float32``. Also, a value of
            ``None`` specifies that the type of matrix should be auto-detected
            which will result in either an ``int``, a ``float`` or a ``complex``
            matrix. Note that this could cause problems if you initially create
            a matrix of integers then carry out an in-place operation which
            results in floats or complexes. A value of ``True`` specifies the
            maximal information data type be used which is currently
            ``numpy.complex128``. Defaults to ``None``.

        copy : bool, optional
            Specifies whether the source data is to be copied or a view is
            created for the data. This only makes sense when the source data is
            a NumPy ``ndarray`` and the specified ``dtype`` is either ``None``
            or matches the source array ``dtype``. This allows wrapping of
            existing large NumPy matrices, otherwise a copy is always
            performed. If you specify ``copy=False`` for source data which is
            not an ``ndarray`` then the flag will be ignored and the data will
            be copied anyway. If the source data is not numeric then a
            ``MutableDenseMatrix`` will be created regardless of the ``copy``
            flag. Defaults to ``True``.
        """

        import numpy as np

        if dtype is True:
            dtype = np.complex128

        view_arr = None

        if copy is False:
            if len(args) == 1:
                view_arr = args[0]

                if not isinstance(view_arr, np.ndarray):
                    copy, view_arr = True, None
                elif len(view_arr.shape) != 2:
                    raise TypeError("'copy=False' for NumPyMatrix without rows and cols specified requires 2D source NumPy ndarray")
                else:
                    rows, cols = view_arr.shape

            elif len(args) == 3:
                rows, cols, view_arr = args

                if not isinstance(view_arr, np.ndarray):
                    copy, view_arr = True, None

            else:
                raise TypeError("'copy=False' requires a matrix be initialized as rows, cols, [list]")

            if copy is False:
                if not issubclass(view_arr.dtype.type, np.number):
                    raise TypeError("'copy=False' for NumPyMatrix with non-numeric source NumPy ndarray not allowed")
                if dtype is not None and dtype != view_arr.dtype:
                    raise TypeError("'copy=False' for NumPyMatrix with 'dtype' requires 'dtype' to match source NumPy ndarray")

                if isinstance(view_arr, np.matrix):
                    view_arr = np.ndarray(view_arr.shape, dtype=view_arr.dtype,
                            buffer=view_arr)

        if copy is True:
            if len(args) == 1 and isinstance(args[0], np.ndarray):
                if dtype is None and issubclass(args [0].dtype.type, np.number):
                    dtype = args[0].dtype

            elif len(args) == 3 and isinstance(args[2], np.ndarray):
                if dtype is None and issubclass(args [2].dtype.type, np.number):
                    dtype = args[2].dtype

                if len(args[2].shape) != 1:
                    args = args[:2] + (np.ndarray(args[2].size,
                            dtype=args[2].dtype, buffer=args[2]),)

            rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)

        if view_arr is None:
            if dtype is None:
                dtype = _detect_dtype(flat_list)

                if dtype is None:
                    if copy is True:
                        return MutableDenseMatrix(*args, copy=copy, **kwargs)
                    else:
                        raise TypeError("'copy=False' for NumPyMatrix requires numeric source data")

            view_arr = np.array([flat_list[cols * r : cols * (r + 1)]
                    for r in range (rows)], dtype=dtype)

        self       = object.__new__(cls)
        self.rows  = rows
        self.cols  = cols
        self.dtype = view_arr.dtype
        self._arr  = view_arr
        self._mat  = self.NDArrayFlatSympified(view_arr)

        return self


    # basic operations

    def __array__(self, *args):
        import numpy as np

        if len(args) > 1:
            raise TypeError('__array__() takes at most 1 argument (%d given)' % len(args))

        dtype = np.dtype(args[0]) if args else self.dtype

        if dtype != self.dtype:
            return np.array(self._arr, dtype)

        return self._arr

    def __abs__(self):
        return NumPyMatrix(self._arr.__abs__(), copy=False)

    def __neg__(self):
        return NumPyMatrix(self._arr.__neg__(), copy=False)

    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return (-self).__add__(other)

    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return self.__add__(-other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return self.__add__(other)

    @call_highest_priority('__radd__')
    def __add__(self, other):
        import numpy as np

        if isinstance(other, NumPyMatrix):
            other = other._arr
        elif not isinstance(other, np.ndarray):
            return super().__add__(other)

        if self.shape != other.shape:
            raise ShapeError("Matrix size mismatch: %s + %s" % (
                self.shape, other.shape))

        m = self._arr.__add__(other)

        return NumPyMatrix(m, copy=not issubclass(m.dtype.type, np.number))

    @call_highest_priority('__rmatmul__')
    def __matmul__(self, other):
        import numpy as np

        if isinstance(other, (NumPyMatrix, np.ndarray)):
            return self.__mul__(other)

        return super().__matmul__(other)

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        import numpy as np

        if isinstance(other, NumPyMatrix):
            other = other._arr

        elif not isinstance(other, np.ndarray):
            if _is_scalar(other):
                m = self._arr.__mul__(other)

                return NumPyMatrix(m, copy=not issubclass(m.dtype.type, np.number))

            return MutableDenseMatrix(self).__mul__(other)

        if self.cols != other.shape[0]:
            raise ShapeError("Matrix size mismatch: %s * %s." % (
                self.shape, other.shape))

        m = self._arr.__matmul__(other)

        return NumPyMatrix(m, copy=not issubclass(m.dtype.type, np.number))

    @call_highest_priority('__matmul__')
    def __rmatmul__(self, other):
        import numpy as np

        if isinstance(other, (NumPyMatrix, np.ndarray)):
            return self.__rmul__(other)

        return super().__rmatmul__(other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        import numpy as np

        if isinstance(other, NumPyMatrix):
            other = other._arr

        elif not isinstance(other, np.ndarray):
            if _is_scalar(other):
                m = self._arr.__rmul__(other)

                return NumPyMatrix(m, copy=not issubclass(m.dtype.type, np.number))

            return MutableDenseMatrix(self).__rmul__(other)

        if self.rows != other.shape[1]:
            raise ShapeError("Matrix size mismatch: %s * %s." % (
                other.shape, self.shape))

        m = self._arr.__rmatmul__(other)

        return NumPyMatrix(m, copy=not issubclass(m.dtype.type, np.number))

    @call_highest_priority('__rtruediv__')
    def __truediv__(self, other):
        return self.__div__(other)

    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        import numpy as np

        if _is_scalar(other):
            m = self._arr.__truediv__(other)

            return NumPyMatrix(m, copy=not issubclass(m.dtype.type, np.number))

        return MutableDenseMatrix(self).__div__(other)

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        import numpy as np
        import numpy.linalg as npl

        if not self.is_square:
            raise NonSquareMatrixError()

        if _is_integer(other):
            m = npl.linalg.matrix_power(self._arr, int(other))

            return NumPyMatrix(m, copy=not issubclass(m.dtype.type, np.number))

        if _is_scalar(other):
            return super().__pow__(other)

        return MutableDenseMatrix(self).__pow__(other)


    # housekeeping functions

    def as_mutable(self):
        """Returns a mutable version of this matrix, basically a copy since all
        NumPy matrices are mutable."""

        return NumPyMatrix(self)


    # higher level functions

    def det(self, *args, **kwargs):
        """Returns the determinant of this matrix using ``numpy.linalg.det``.
        This determinant may be imprecise for floating point or complex matrices
        but will always be precise for integer matrices (up to the limit of
        integer overflow where it will become an imprecise float)."""

        import numpy.linalg as npl

        if not self.is_square:
            raise NonSquareMatrixError()

        return nsimplify(npl.det(self._arr), full=True)

    det_bareiss   = det
    det_berkowitz = det
    det_LU        = det

    def rank(self, *args, **kwargs):
        """Returns the rank of the matrix using ``numpy.linalg.matrix_rank``."""

        import numpy.linalg as npl

        return sympify(npl.matrix_rank(self._arr))
