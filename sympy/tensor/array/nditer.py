from sympy.tensor.array.dense_ndim_array import DenseNDimArray
from sympy.tensor.array.sparse_ndim_array import SparseNDimArray
from sympy.tensor.array.ndim_array import NDimArray
from sympy.core.compatibility import Iterable

class nditer(object):

    def __init__(self, iterable):
        from sympy.matrices.matrices import MatrixBase

        if not isinstance(iterable, (Iterable, MatrixBase)):
            raise NotImplementedError("Data type not yet supported")
        if isinstance(iterable, list) and isinstance(iterable[0], NDimArray):
            temp = []
            for i in range(len(iterable)):
                temp += iterable[0].tolist()
            iterable = temp
        self._iter = iterable
        self._idx = 0

    def __iter__(self):
        return self

    def __next__(self):
        from sympy.matrices.matrices import MatrixBase

        try:
            if isinstance(self._iter, list):
                result = self._iter[self._idx]

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
        