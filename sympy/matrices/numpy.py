from __future__ import division, print_function

from sympy.core.sympify import sympify

from .common import NonSquareMatrixError
from .dense import MutableDenseMatrix

class NumPyMatrix(MutableDenseMatrix):
    """Matrix class which can wrap or create a new NumPy matrix (ndarray) for
    use in SymPy. Dispatches supported operations for faster execution to NumPy.
    """

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
    def _new(cls, *args, dtype=None, **kwargs):
        """
        Parameters
        ==========

        dtype : numpy.dtype, optional
            This can be either a NumPy character dtype specification like
            ``'f4'`` or a class like ``numpy.float32``. Defaults to
            ``numpy.complex128``.

        copy : bool, optional
            Specifies whether the source data is to be copied or a view is
            created for the data. This only makes sense when the source data is
            a NumPy ``ndarray`` and the specified ``dtype`` is either ``None``
            or matches the source array ``dtype``. This allows wrapping of
            existing large NumPy matrices, otherwise a copy is always
            performed. Defaults to ``True``.
        """

        import numpy as np

        view_arr = None

        if kwargs.get('copy', True) is False:
            if len(args) == 1:
                view_arr = args[0]

                if not isinstance(view_arr, (np.ndarray, np.matrix)):
                    raise TypeError("'copy=False' for NumPyMatrix requires NumPy matrix as source")
                if len(view_arr.shape) != 2:
                    raise TypeError("'copy=False' for NumPyMatrix without rows and cols specified requires 2D NumPy matrix as source")

                rows, cols = view_arr.shape

            elif len(args) == 3:
                rows, cols, view_arr = args

                if not isinstance(view_arr, (np.ndarray, np.matrix)):
                    raise TypeError("'copy=False' for NumPyMatrix requires NumPy matrix as source")

            else:
                raise TypeError("'copy=False' requires a matrix be initialized as rows, cols, [list]")

            if dtype is not None and dtype != view_arr.dtype:
                raise TypeError("'copy=False' for NumPyMatrix with 'dtype' requires 'dtype' to match source NumPy matrix")

            if isinstance(view_arr, np.matrix):
                view_arr = np.ndarray(view_arr.shape, dtype=view_arr.dtype, buffer=view_arr)

        else:
            if len(args) == 1 and isinstance(args[0], (np.ndarray, np.matrix)):
                if dtype is None:
                    dtype = args[0].dtype

            elif len(args) == 3 and isinstance(args[2], (np.ndarray, np.matrix)):
                if dtype is None:
                    dtype = args[2].dtype

                if len(args[2].shape) != 1:
                    args = args[:2] + (np.ndarray(args[2].size, dtype=args[2].dtype, buffer=args[2]),)

            rows, cols, flat_list = cls._handle_creation_inputs(*args, **kwargs)
            flat_list             = list(flat_list) # create a shallow copy

        self       = object.__new__(cls)
        self.rows  = rows
        self.cols  = cols
        self.dtype = np.complex128 if dtype is None else np.dtype(dtype)
        _arr       = np.ndarray((rows, cols), dtype=self.dtype) if view_arr is None else view_arr
        self._arr  = _arr
        self._mat  = self.NDArrayFlatSympified(_arr)

        if view_arr is None:
            for r in range(rows):
                for c in range(cols):
                    _arr[r, c] = flat_list[cols * r + c]

        return self

    def det(self, *args, **kwargs): # override det and calculate usign NumPy
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
