from __future__ import print_function, division
import collections
from sympy import Matrix, Integer, sympify


class NDimArray(object):
    """

    Examples
    ========

    Create an N-dim array of zeros:

    >>> from sympy.tensor.array import MutableDenseNDimArray
    >>> a = MutableDenseNDimArray.zeros(2, 3, 4)
    >>> a
    [[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]]

    Create an N-dim array from a list;

    >>> a = MutableDenseNDimArray([[2, 3], [4, 5]])
    >>> a
    [[2, 3], [4, 5]]

    >>> b = MutableDenseNDimArray([[[1, 2], [3, 4], [5, 6]], [[7, 8], [9, 10], [11, 12]]])
    >>> b
    [[[1, 2], [3, 4], [5, 6]], [[7, 8], [9, 10], [11, 12]]]

    Create an N-dim array from a flat list with dimension shape:

    >>> a = MutableDenseNDimArray(2, 3, [1, 2, 3, 4, 5, 6])
    >>> a
    [[1, 2, 3], [4, 5, 6]]

    Create an N-dim array from a matrix:

    >>> from sympy import Matrix
    >>> a = Matrix([[1,2],[3,4]])
    >>> a
    Matrix([
    [1, 2],
    [3, 4]])
    >>> b = MutableDenseNDimArray(a)
    >>> b
    [[1, 2], [3, 4]]

    Arithmetic operations on N-dim arrays

    >>> a = MutableDenseNDimArray(2, 2, [1, 1, 1, 1])
    >>> b = MutableDenseNDimArray(2, 2, [4, 4, 4, 4])
    >>> c = a + b
    >>> c
    [[5, 5], [5, 5]]
    >>> a - b
    [[-3, -3], [-3, -3]]

    """
    def __new__(cls, *args, **kwargs):
        from sympy.tensor.array import MutableDenseNDimArray
        return MutableDenseNDimArray(*args, **kwargs)

    def _parse_index(self, index):

        if isinstance(index, (int, Integer)):
            if index >= self._loop_size:
                raise ValueError("index out of range")
            return index

        if len(index) != self._rank:
            raise ValueError('Wrong number of array axes')

        real_index = 0
        # check if input index can exist in current indexing
        for i in range(self._rank):
            if index[i] > self.shape[i]:
                raise ValueError('Index ' + str(index) + ' out of border')
            real_index = real_index*self.shape[i] + index[i]

        return real_index

    def _get_tuple_index(self, integer_index):
        index = []
        for i, sh in enumerate(reversed(self.shape)):
            index.append(integer_index % sh)
            integer_index //= sh
        index.reverse()
        return tuple(index)

    def _setter_iterable_check(self, value):
        if isinstance(value, (collections.Iterable, Matrix, NDimArray)):
            raise NotImplementedError

    def __len__(self):
        """Overload common function len(). Returns number of elements in array.

        Examples
        ========

        >>> from sympy.tensor.array.dense_ndim_array import MutableDenseNDimArray
        >>> a = MutableDenseNDimArray.zeros(3,3)
        >>> a
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        >>> len(a)
        9

        """
        return self._loop_size

    @property
    def shape(self):
        """
        Returns array shape (dimension).

        Examples
        ========

        >>> from sympy.tensor.array.dense_ndim_array import MutableDenseNDimArray
        >>> a = MutableDenseNDimArray.zeros(3,3)
        >>> a.shape
        (3, 3)

        """
        return self._shape

    def rank(self):
        """
        Returns rank of array.

        Examples
        ========

        >>> from sympy.tensor.array.dense_ndim_array import MutableDenseNDimArray
        >>> a = MutableDenseNDimArray.zeros(3,4,5,6,3)
        >>> a.rank()
        5

        """
        return self._rank

    def __str__(self):
        """Returns string, allows to use standard functions print() and str().

        Examples
        ========

        >>> from sympy.tensor.array import MutableDenseNDimArray
        >>> a = MutableDenseNDimArray.zeros(2, 2)
        >>> a
        [[0, 0], [0, 0]]

        """
        def f(sh, shape_left, i, j):
            if len(shape_left) == 1:
                return "["+", ".join([str(self[e]) for e in range(i, j)])+"]"

            sh //= shape_left[0]
            return "[" + ", ".join([f(sh, shape_left[1:], i+e*sh, i+(e+1)*sh) for e in range(shape_left[0])]) + "]" # + "\n"*len(shape_left)

        return f(self._loop_size, self.shape, 0, self._loop_size)

        out_str = ''

        # forming output string
        for i, el in enumerate(self):

            out_str += str(el) + '  '
            chidx = i+1
            for sh in reversed(self.shape):
                if chidx % sh == 0:
                    out_str += '\n'
                    chidx //= sh

        return out_str

    def __repr__(self):
        return self.__str__()

    def tolist(self):
        """
        Conveting MutableDenseNDimArray to one-dim list

        Examples
        ========

        >>> from sympy.tensor.array import MutableDenseNDimArray
        >>> a = MutableDenseNDimArray(2, 2, [1, 2, 3, 4])
        >>> a
        [[1, 2], [3, 4]]
        >>> b = a.tolist()
        >>> b
        [[1, 2], [3, 4]]
        """

        def f(sh, shape_left, i, j):
            if len(shape_left) == 1:
                return [self[e] for e in range(i, j)]
            result = []
            sh //= shape_left[0]
            for e in range(shape_left[0]):
                result.append(f(sh, shape_left[1:], i+e*sh, i+(e+1)*sh))
            return result

        return f(self._loop_size, self.shape, 0, self._loop_size)

    def __add__(self, other):
        if not isinstance(other, NDimArray):
            raise TypeError(str(other))

        if self.shape != other.shape:
            raise ValueError("array shape mismatch")
        result_list = [i+j for i,j in zip(self, other)]

        return type(self)(*(self.shape + (result_list,)))

    def __sub__(self, other):
        if not isinstance(other, NDimArray):
            raise TypeError(str(other))

        if self.shape != other.shape:
            raise ValueError("array shape mismatch")
        result_list = [i-j for i,j in zip(self, other)]

        return type(self)(*(self.shape + (result_list,)))

    def __mul__(self, other):
        if isinstance(other, (collections.Iterable,NDimArray, Matrix)):
            raise ValueError("scalar expected")
        other = sympify(other)
        result_list = [i*other for i in self]
        return type(self)(*(self.shape + (result_list,)))

    def __rmul__(self, other):
        if isinstance(other, (collections.Iterable,NDimArray, Matrix)):
            raise ValueError("scalar expected")
        other = sympify(other)
        result_list = [other*i for i in self]
        return type(self)(*(self.shape + (result_list,)))

    def __div__(self, other):
        if isinstance(other, (collections.Iterable,NDimArray, Matrix)):
            raise ValueError("scalar expected")
        other = sympify(other)
        result_list = [i/other for i in self]
        return type(self)(*(self.shape + (result_list,)))

    def __rdiv__(self, other):
        raise TypeError('unsupported operation on NDimArray')

    def __eq__(self, other):
        """
        NDimArray instances can be compared to each other.
        Instances equal if they have same shape and data.

        Examples
        ========

        >>> from sympy.tensor.array import MutableDenseNDimArray
        >>> a = MutableDenseNDimArray.zeros(2, 3)
        >>> b = MutableDenseNDimArray.zeros(2, 3)
        >>> a == b
        True
        >>> c = a.reshape(3, 2)
        >>> c == b
        False
        >>> a[0,0] = 1
        >>> b[0,0] = 2
        >>> a == b
        False
        """
        if not isinstance(other, NDimArray):
            return False
        return (self.shape == other.shape) and (list(self) == list(other))

    def __ne__(self, other):
        return not self.__eq__(other)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__
