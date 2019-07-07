from sympy.core.sympify import _sympify
from sympy import Tuple, Basic
from sympy.tensor.array.ndim_array import NDimArray
from sympy.tensor.array.sparse_ndim_array import SparseNDimArray
from sympy.tensor.array.dense_ndim_array import DenseNDimArray
from sympy.core.compatibility import Iterable

import functools

# Row-based linked list sparse matrix(LIL)
class lil(Basic, NDimArray):
    """
    Create a sparse array with Row-based linked list sparse matrix(LIL) format.
    This data structure is efficient in incremental construction of sparse arrays.

    Note
    ====

    This is an experimental implementation of algorithm, the compatibility with
    existing format needs to be verified. Only some basic operations are implemented.

    Examples
    ========

    >>> from sympy.tensor.array.sp_utils import lil
    >>> a = lil([[0, 1, 0], [2, 0, 0]])
    >>> a._data
    [[1], [2]]
    >>> a._rows
    [[1], [0]]
    >>> a[0, 1]
    1
    >>> a[0, 1] = 3
    >>> a.tolist()
    [[0, 3, 0], [2, 0, 0]]
    """
    def __new__(cls, iterable=None, shape=None, **kwargs):
        shape, flat_list = cls._handle_ndarray_creation_inputs(iterable, shape, **kwargs)
        shape = Tuple(*map(_sympify, shape))
        cls._check_special_bounds(flat_list, shape)
        loop_size = functools.reduce(lambda x,y: x*y, shape) if shape else 0

        # Ideally we should use the Array module for this empty initialization, but
        # it is not yet supported in SymPy
        data = cls._empty(shape)
        rows = cls._empty(shape)

        self = Basic.__new__(cls, data, rows, shape, **kwargs)
        self._shape = shape
        self._rank = len(shape)
        self._loop_size = loop_size

        new_dict = {}
        if isinstance(iterable, Iterable):
            iterable = SparseNDimArray(iterable)

        if isinstance(iterable, SparseNDimArray):
            new_dict = iterable._sparse_array
        # if it is a dense array, it would be cast to a default sparse array
        # and then convert to a LIL format
        elif isinstance(iterable, DenseNDimArray):
            sp = SparseNDimArray(iterable)
            new_dict = sp._sparse_array
        else:
            raise NotImplementedError("Data type not yet supported")

        # Add non zero value to internal list.
        # This operation can be simplified once Array module supports tuple index
        for k, v in new_dict.items():
            idx = self._get_tuple_index(k)
            self._get_row(data, idx).append(v)
            self._get_row(rows, idx).append(idx[-1])

        self._data = data
        self._rows = rows
        return self

    @classmethod
    def _empty(cls, shape):
        def f(s):
            if len(s) == 0:
                return []
            if len(s) == 1:
                return [[] for i in range(s[0])]
            arr = f(s[1:])
            return [arr for i in range(s[0])]

        if not shape:
            raise ValueError("Shape must be defined")
        return f(shape[:-1])

    def _get_row(self, iterable, index):
        temp_iter = iterable
        for i in index[:-1]:
            temp_iter = temp_iter[i]
        return temp_iter


    def __getitem__(self, index):
        if not isinstance(index, (tuple, Tuple)):
            index = self._get_tuple_index(index)
        row_values = self._get_row(self._data, index)
        row_indices = self._get_row(self._rows, index)

        if index[-1] in row_indices:
            value_index = row_indices.index(index[-1])
            return row_values[value_index]
        else:
            return 0

    # For mutable arrays
    def __setitem__(self, index, value):
        if not isinstance(index, (tuple, Tuple)):
            index = self._get_tuple_index(index)
        row_values = self._get_row(self._data, index)
        row_indices = self._get_row(self._rows, index)

        if index[-1] in row_indices:
            value_index = row_indices.index(index[-1])
            row_values[value_index] = value
        else:
            row_values.append(value)
            row_indices.append(index[-1])
