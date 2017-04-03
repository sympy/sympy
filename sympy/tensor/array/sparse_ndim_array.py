from __future__ import print_function, division
import functools

import itertools

from sympy.core.sympify import _sympify

from sympy import S, Dict, Basic, Tuple
from sympy.tensor.array.mutable_ndim_array import MutableNDimArray
from sympy.tensor.array.ndim_array import NDimArray, ImmutableNDimArray


class SparseNDimArray(NDimArray):

    def __new__(self, *args, **kwargs):
        return ImmutableSparseNDimArray(*args, **kwargs)

    def __getitem__(self, index):
        """
        Get an element from a sparse N-dim array.

        Examples
        ========

        >>> from sympy import MutableSparseNDimArray
        >>> a = MutableSparseNDimArray(range(4), (2, 2))
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

        Symbolic indexing:

        >>> from sympy.abc import i, j
        >>> a[i, j]
        [[0, 1], [2, 3]][i, j]

        Replace `i` and `j` to get element `(0, 0)`:

        >>> a[i, j].subs({i: 0, j: 0})
        0

        """
        syindex = self._check_symbolic_index(index)
        if syindex is not None:
            return syindex

        # `index` is a tuple with one or more slices:
        if isinstance(index, tuple) and any([isinstance(i, slice) for i in index]):

            def slice_expand(s, dim):
                if not isinstance(s, slice):
                        return (s,)
                start, stop, step = s.indices(dim)
                return [start + i*step for i in range((stop-start)//step)]

            sl_factors = [slice_expand(i, dim) for (i, dim) in zip(index, self.shape)]
            eindices = itertools.product(*sl_factors)
            array = [self._sparse_array.get(self._parse_index(i), S.Zero) for i in eindices]
            nshape = [len(el) for i, el in enumerate(sl_factors) if isinstance(index[i], slice)]
            return type(self)(array, nshape)
        else:
            # `index` is a single slice:
            if isinstance(index, slice):
                start, stop, step = index.indices(self._loop_size)
                retvec = [self._sparse_array.get(ind, S.Zero) for ind in range(start, stop, step)]
                return retvec
            # `index` is a number or a tuple without any slice:
            else:
                index = self._parse_index(index)
                return self._sparse_array.get(index, S.Zero)

    @classmethod
    def zeros(cls, *shape):
        """
        Return a sparse N-dim array of zeros.
        """
        return cls({}, shape)

    def tomatrix(self):
        """
        Converts MutableDenseNDimArray to Matrix. Can convert only 2-dim array, else will raise error.

        Examples
        ========

        >>> from sympy import MutableSparseNDimArray
        >>> a = MutableSparseNDimArray([1 for i in range(9)], (3, 3))
        >>> b = a.tomatrix()
        >>> b
        Matrix([
        [1, 1, 1],
        [1, 1, 1],
        [1, 1, 1]])
        """
        from sympy.matrices import SparseMatrix
        if self.rank() != 2:
            raise ValueError('Dimensions must be of size of 2')

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


class ImmutableSparseNDimArray(SparseNDimArray, ImmutableNDimArray):

    def __new__(cls, iterable=None, shape=None, **kwargs):
        from sympy.utilities.iterables import flatten

        shape, flat_list = cls._handle_ndarray_creation_inputs(iterable, shape, **kwargs)
        shape = Tuple(*map(_sympify, shape))
        loop_size = functools.reduce(lambda x,y: x*y, shape) if shape else 0

        # Sparse array:
        if isinstance(flat_list, (dict, Dict)):
            sparse_array = Dict(flat_list)
        else:
            sparse_array = {}
            for i, el in enumerate(flatten(flat_list)):
                if el != 0:
                    sparse_array[i] = _sympify(el)

        sparse_array = Dict(sparse_array)

        self = Basic.__new__(cls, sparse_array, shape, **kwargs)
        self._shape = shape
        self._rank = len(shape)
        self._loop_size = loop_size
        self._sparse_array = sparse_array

        return self

    def __setitem__(self, index, value):
        raise TypeError("immutable N-dim array")


class MutableSparseNDimArray(MutableNDimArray, SparseNDimArray):

    def __new__(cls, iterable=None, shape=None, **kwargs):
        from sympy.utilities.iterables import flatten

        shape, flat_list = cls._handle_ndarray_creation_inputs(iterable, shape, **kwargs)
        self = object.__new__(cls)
        self._shape = shape
        self._rank = len(shape)
        self._loop_size = functools.reduce(lambda x,y: x*y, shape) if shape else 0

        # Sparse array:
        if isinstance(flat_list, (dict, Dict)):
            self._sparse_array = dict(flat_list)
            return self

        self._sparse_array = {}

        for i, el in enumerate(flatten(flat_list)):
            if el != 0:
                self._sparse_array[i] = _sympify(el)

        return self

    def __setitem__(self, index, value):
        """Allows to set items to MutableDenseNDimArray.

        Examples
        ========

        >>> from sympy import MutableSparseNDimArray
        >>> a = MutableSparseNDimArray.zeros(2, 2)
        >>> a[0, 0] = 1
        >>> a[1, 1] = 1
        >>> a
        [[1, 0], [0, 1]]


        """
        index = self._parse_index(index)
        if not isinstance(value, MutableNDimArray):
            value = _sympify(value)

        if isinstance(value, NDimArray):
            return NotImplementedError

        if value == 0 and index in self._sparse_array:
            self._sparse_array.pop(index)
        else:
            self._sparse_array[index] = value
