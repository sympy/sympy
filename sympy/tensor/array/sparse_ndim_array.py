from __future__ import print_function, division
import functools
from sympy import sympify, Matrix, S, Dict, flatten, SparseMatrix
from sympy.tensor.array.mutable_ndim_array import MutableNDimArray
from sympy.tensor.array.ndim_array import NDimArray


class SparseNDimArray(NDimArray):

    def __getitem__(self, index):
        """
        Get an element from a sparse N-dim array.

        Examples
        ========

        >>> from sympy.tensor.array import MutableSparseNDimArray
        >>> a = MutableSparseNDimArray(2, 2, range(4))
        >>> a
        [[0, 1], [2, 3]]
        >>> a[0, 0]
        0
        >>> a[1, 1]
        3
        >>> a[0]
        0
        >>> a[2]
        2

        """
        index = self._parse_index(index)

        if index in self._sparse_array:
            return self._sparse_array[index]
        else:
            return S.Zero

    @classmethod
    def zeros(cls, *shape):
        """
        Return a sparse N-dim array of zeros.
        """
        return cls(*(shape + ({},)))

    def tomatrix(self):
        """
        Converts MutableDenseNDimArray to Matrix. Can convert only 2-dim array, else will raise error.

        Examples
        ========

        >>> from sympy.tensor.array import MutableSparseNDimArray
        >>> a = MutableSparseNDimArray(3, 3, [1 for i in range(9)])
        >>> b = a.tomatrix()
        >>> b
        Matrix([
        [1, 1, 1],
        [1, 1, 1],
        [1, 1, 1]])
        """
        if self.rank() != 2:
            raise ValueError('Only rank 2 arrays can be converted to matrices')

        mat_sparse = {}
        for key, value in self._sparse_array.items():
            mat_sparse[self._get_tuple_index(key)] = value

        return SparseMatrix(self.shape[0], self.shape[1], mat_sparse)

    def __iter__(self):
        def iterator():
            for i in range(self._loop_size):
                yield self[i]
        return iterator()

    def reshape(self, *newshape):
        new_total_size = functools.reduce(lambda x,y: x*y, newshape)
        if new_total_size != self._loop_size:
            raise ValueError("Invalid reshape parameters " + newshape)

        return type(self)(*(newshape + (self._array,)))




class MutableSparseNDimArray(MutableNDimArray, SparseNDimArray):

    def __new__(cls, *args, **kwargs):

        shape, flat_list = cls._handle_ndarray_creation_inputs(*args, **kwargs)
        self = object.__new__(cls)
        self._shape = shape
        self._rank = len(shape)
        self._loop_size = functools.reduce(lambda x,y: x*y, shape) if shape else 0

        # Sparse array:
        if isinstance(flat_list, (dict, Dict)):
            self._sparse_array = args[-1]
            return self

        self._sparse_array = {}

        for i, el in enumerate(flatten(flat_list)):
            if el != 0:
                self._sparse_array[i] = sympify(el)

        return self

    def __setitem__(self, index, value):
        """Allows to set items to MutableDenseNDimArray.

        Examples
        ========

        >>> from sympy.tensor.array import MutableSparseNDimArray
        >>> a = MutableSparseNDimArray.zeros(2, 2)
        >>> a[0, 0] = 1
        >>> a[1, 1] = 1
        >>> a
        [[1, 0], [0, 1]]


        """
        index = self._parse_index(index)
        if not isinstance(value, MutableNDimArray):
            value = sympify(value)

        if isinstance(value, NDimArray):
            return NotImplementedError

        if value == 0 and index in self._sparse_array:
            self._sparse_array.pop(index)
        else:
            self._sparse_array[index] = value
