from sympy.tensor.array.dense_ndim_array import DenseNDimArray
from sympy.tensor.array.sparse_ndim_array import SparseNDimArray
from sympy.tensor.array.ndim_array import NDimArray
from sympy.core.compatibility import Iterable
from sympy.core.compatibility import SYMPY_INTS
from sympy.core.numbers import Integer


class nditer(object):

    def __init__(self, iterable):
        from sympy.matrices.matrices import MatrixBase

        if not isinstance(iterable, (Iterable, MatrixBase)):
            raise NotImplementedError("Data type not yet supported")
        self._iter = iterable
        self._idx = 0

    def __iter__(self):
        return self

    def __next__(self):
        from sympy.matrices.matrices import MatrixBase

        try:
            if isinstance(self._iter, list):
                temp = self._iter[self._idx]
                if hasattr(temp, '__len__') and len(temp) > 1:
                    result = next(nditer(temp))
                else:
                    result = temp

            elif isinstance(self._iter, DenseNDimArray):
                result = self._iter._array[self._idx]

            elif isinstance(self._iter, SparseNDimArray):
                if self._idx >= len(self._iter):
                    raise StopIteration

                if self._idx in self._iter._sparse_array:
                    result = self._iter._sparse_array[self._idx]
                else:
                    result = 0

            elif isinstance(self._iter, MatrixBase):
                result = self._iter[self._idx]

            elif hasattr(self._iter, '__next__'):
                result = next(self._iter)

            else:
                result = self._iter[self._idx]


        except IndexError:
            raise StopIteration

        self._idx += 1
        return result
        