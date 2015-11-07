from __future__ import print_function, division
import functools
from sympy import sympify, Matrix, flatten
from sympy.tensor.array.mutable_ndim_array import MutableNDimArray
from sympy.tensor.array.ndim_array import NDimArray


class DenseNDimArray(NDimArray):

    def __getitem__(self, index):
        """
        Allows to get items from N-dim array.

        Examples
        ========

        >>> from sympy.tensor.array import MutableDenseNDimArray
        >>> a = MutableDenseNDimArray(2, 2, [0, 1, 2, 3])
        >>> a
        [[0, 1], [2, 3]]
        >>> a[0, 0]
        0
        >>> a[1, 1]
        3

        """
        index = self._parse_index(index)
        return self._array[index]

    @classmethod
    def zeros(cls, *shape):
        list_length = functools.reduce(lambda x, y: x*y, shape)
        return cls._new(*(shape + ([0]*list_length,)))

    def tomatrix(self):
        """
        Converts MutableDenseNDimArray to Matrix. Can convert only 2-dim array, else will raise error.

        Examples
        ========

        >>> from sympy.tensor.array import MutableDenseNDimArray
        >>> a = MutableDenseNDimArray(3, 3, [1 for i in range(9)])
        >>> b = a.tomatrix()
        >>> b
        Matrix([
        [1, 1, 1],
        [1, 1, 1],
        [1, 1, 1]])

        """
        if self.rank() != 2:
            raise ValueError('Dimensions must be of size of 2')

        return Matrix(self.shape[0], self.shape[1], self._array)

    def __iter__(self):
        return self._array.__iter__()

    def reshape(self, *newshape):
        """
        Returns MutableDenseNDimArray instance with new shape. Elements number
        must be        suitable to new shape. The only argument of method sets
        new shape.

        Examples
        ========

        >>> from sympy.tensor.array import MutableDenseNDimArray
        >>> a = MutableDenseNDimArray(2, 3, [1, 2, 3, 4, 5, 6])
        >>> a.shape
        (2, 3)
        >>> a
        [[1, 2, 3], [4, 5, 6]]
        >>> b = a.reshape(3, 2)
        >>> b.shape
        (3, 2)
        >>> b
        [[1, 2], [3, 4], [5, 6]]

        """
        new_total_size = functools.reduce(lambda x,y: x*y, newshape)
        if new_total_size != self._loop_size:
            raise ValueError("Invalid reshape parameters " + newshape)

        return type(self)(*(newshape + (self._array,)))


class MutableDenseNDimArray(DenseNDimArray, MutableNDimArray):

    def __new__(cls, *args, **kwargs):
        return cls._new(*args, **kwargs)

    @classmethod
    def _new(cls, *args, **kwargs):
        shape, flat_list = cls._handle_ndarray_creation_inputs(*args, **kwargs)
        flat_list = flatten(flat_list)
        self = object.__new__(cls)
        self._shape = shape
        self._array = list(flat_list)
        self._rank = len(shape)
        self._loop_size = functools.reduce(lambda x,y: x*y, shape) if shape else 0
        return self

    def __setitem__(self, index, value):
        """Allows to set items to MutableDenseNDimArray.

        Examples
        ========

        >>> from sympy.tensor.array import MutableDenseNDimArray
        >>> a = MutableDenseNDimArray.zeros(2,2)
        >>> a[0,0] = 1
        >>> a[1,1] = 1
        >>> a
        [[1, 0], [0, 1]]

        """
        index = self._parse_index(index)
        self._setter_iterable_check(value)
        value = sympify(value)

        self._array[index] = value
