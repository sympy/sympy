from __future__ import division, print_function

from sympy import Float
from sympy.core.sympify import sympify

from .common import NonSquareMatrixError
from .dense import MutableDenseMatrix


def _detect_dtype(l): # Is there a better way to do this?
    import numpy as np

    try:
        f = [float(e) for e in l]
    except TypeError: # complex or symbolic source
        pass

    else:
        try:
            i = [int(e) for e in l]
        except TypeError: # have infinities and/or nans
            return np.float64

        else:
            if any(isinstance(e, Float) for e in l) or i != f:
                return np.float64
            else:
                return np.int64

    try:
        all(complex(e) for e in l)

        return np.complex128

    except TypeError: # non-numeric source data
        raise TypeError('source data for NumPyMatrix must be numeric') from None

    return np.complex128


class NumPyMatrix(MutableDenseMatrix):
    """Matrix class which can wrap or create a new NumPy matrix (ndarray) for
    use in SymPy. Dispatches supported operations for faster execution to NumPy.
    """

    _class_priority       = 2
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
                return np.ndarray.__new__(cls, src.size, dtype=src.dtype, buffer=src)

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
            results in floats or complexes. A value of ``'full'`` specifies the
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
            be copied anyway. Defaults to ``True``.
        """

        import numpy as np

        if dtype == 'full':
            dtype = np.complex128

        view_arr = None

        if copy is False:
            if len(args) == 1:
                view_arr = args[0]

                if not isinstance(view_arr, np.ndarray):
                    copy, view_arr = True, None
                elif len(view_arr.shape) != 2:
                    raise TypeError("'copy=False' for NumPyMatrix without rows and cols specified requires 2D NumPy matrix as source")
                else:
                    rows, cols = view_arr.shape

            elif len(args) == 3:
                rows, cols, view_arr = args

                if not isinstance(view_arr, np.ndarray):
                    copy, view_arr = True, None

            else:
                raise TypeError("'copy=False' requires a matrix be initialized as rows, cols, [list]")

            if copy is False:
                if dtype is not None and dtype != view_arr.dtype:
                    raise TypeError("'copy=False' for NumPyMatrix with 'dtype' requires 'dtype' to match source NumPy matrix")

                if isinstance(view_arr, np.matrix):
                    view_arr = np.ndarray(view_arr.shape, dtype=view_arr.dtype, buffer=view_arr)

        if copy is True:
            if len(args) == 1 and isinstance(args[0], np.ndarray):
                if dtype is None:
                    dtype = args[0].dtype

            elif len(args) == 3 and isinstance(args[2], np.ndarray):
                if dtype is None:
                    dtype = args[2].dtype

                if len(args[2].shape) != 1:
                    args = args[:2] + (np.ndarray(args[2].size, dtype=args[2].dtype, buffer=args[2]),)

            rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)

        self      = object.__new__(cls)
        self.rows = rows
        self.cols = cols

        if view_arr is None:
            if dtype is None:
                dtype = _detect_dtype(flat_list)

            view_arr = np.array([flat_list[cols * r : cols * (r + 1)]
                    for r in range (rows)], dtype=dtype)

        self.dtype = np.dtype(dtype)
        self._arr  = view_arr
        self._mat  = self.NDArrayFlatSympified(view_arr)

        return self

    def det(self, *args, **kwargs):
        import numpy.linalg as npla

        if not self.is_square:
            raise NonSquareMatrixError()

        if self.rows == 0:
            return self.one
            # sympy/matrices/tests/test_matrices.py contains a test that
            # suggests that the determinant of a 0 x 0 matrix is one, by
            # convention.

        return sympify(npla.det(self._arr))

    det_bareiss   = det
    det_berkowitz = det
    det_LU        = det

    def rank(self, *args, **kwargs):
        import numpy.linalg as npla

        return sympify(npla.matrix_rank(self._arr))
