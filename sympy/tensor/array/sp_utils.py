from sympy.core.sympify import _sympify
from sympy import Tuple, Basic
from sympy.tensor.array.ndim_array import NDimArray
from sympy.tensor.array.sparse_ndim_array import SparseNDimArray
from sympy.tensor.array.dense_ndim_array import DenseNDimArray
from sympy.core.compatibility import Iterable

import functools

# A set of functions that could be adde to the base class
# in order to avoid code repetition 
def check_bound(shape, index):
    if len(shape) != len(index):
        return False
    return all([shape[i] > index[i] for i in range(len(shape))])

# Row-based linked list sparse matrix(LIL)
class lil(Basic, SparseNDimArray):
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
        #shape = Tuple(*map(_sympify, shape))
        shape = Tuple.fromiter(_sympify(i) for i in shape)
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
        # if it is a dense array, it would be cast to a default sparse array
        # and then convert to a LIL format
        if isinstance(iterable, Iterable):
            iterable = SparseNDimArray(iterable)

        if isinstance(iterable, SparseNDimArray):
            new_dict = iterable._sparse_array

        # TODO: Enable the initialization with (data, rows) as a tuple
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

    def __iter__(self):
        def iterator():
            for i in range(self._loop_size):
                yield self[i]
        return iterator()

    def __getitem__(self, index):
        if not isinstance(index, (tuple, Tuple)):
            index = self._get_tuple_index(index)
        if not check_bound(self._shape, index):
            raise ValueError('Index ' + str(index) + ' out of border')

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
        if not check_bound(self._shape, index):
            raise ValueError('Index ' + str(index) + ' out of border')

        row_values = self._get_row(self._data, index)
        row_indices = self._get_row(self._rows, index)

        if index[-1] in row_indices:
            value_index = row_indices.index(index[-1])
            row_values[value_index] = value
        else:
            row_values.append(value)
            row_indices.append(index[-1])

    # TODO: convert to other formats


class coo(Basic, SparseNDimArray):
    '''
    Coordinate list format. This format is not very efficient for calculation. But it
    can be easily cast to csr/csc format. 
    '''
    def __new__(cls, iterable=None, shape=None, **kwargs):
        shape, flat_list = cls._handle_ndarray_creation_inputs(iterable, shape, **kwargs)
        shape = Tuple.fromiter(_sympify(i) for i in shape)
        cls._check_special_bounds(flat_list, shape)
        loop_size = functools.reduce(lambda x,y: x*y, shape) if shape else 0

        data = []
        coor = [[] for i in range(len(shape))]

        self = Basic.__new__(cls, data, coor, shape, **kwargs)
        self._data = data
        self._coor = coor
        self._shape = shape
        self._rank = len(shape)
        self._loop_size = loop_size

        if isinstance(iterable, Iterable):
            iterable = SparseNDimArray(iterable)

        if isinstance(iterable, SparseNDimArray):
            for k, v in sorted(iterable._sparse_array.items()):
                data.append(v)
                idx = self._get_tuple_index(k)
                for i in range(len(shape)):
                    coor[i].append(idx[i])
        # TODO: Enable the initialization with (data, coor) as a tuple
        else:
            raise NotImplementedError("Data type not yet supported")

        return self

    def __getitem__(self, index):
        if not isinstance(index, (tuple, Tuple)):
            index = self._get_tuple_index(index)
        if not check_bound(self._shape, index):
            raise ValueError('Index ' + str(index) + ' out of border')

        for i in range(len(self._data)):
            if all([index[j] == self._coor[j][i] for j in range(self._rank)]):
                return self._data[i]
        return 0

    def __iter__(self):
        def iterator():
            for i in range(self._loop_size):
                yield self[i]
        return iterator()

    # TODO: convert to other formats
    def tocsr(self):


class csr(Basic, SparseNDimArray):
    '''
    Compressed Sparse Row. This format is widely used. It has the following advantages:
    Efficient item access and slicing
    Fast arithmetics
    Fast vector products
    '''
    def __new__(cls, iterable=None, shape=None, **kwargs):
        shape, flat_list = cls._handle_ndarray_creation_inputs(iterable, shape, **kwargs)
        shape = Tuple.fromiter(_sympify(i) for i in shape)
        cls._check_special_bounds(flat_list, shape)
        loop_size = functools.reduce(lambda x,y: x*y, shape) if shape else 0

        data = []
        col_ind = []
        row_ptr = cls._empty(shape)

        self = Basic.__new__(cls, data, col_ind, row_ptr, shape, **kwargs)
        self._data = data
        self._col_ind = col_ind
        self._row_ptr = row_ptr
        self._shape = shape
        self._rank = len(shape)
        self._loop_size = loop_size

        counter = 0
        if isinstance(iterable, Iterable):
            iterable = SparseNDimArray(iterable)

        if isinstance(iterable, SparseNDimArray):
            for k, v in sorted(iterable._sparse_array.items()):
                data.append(v)
                idx = self._get_tuple_index(k)
                col_ind.append(idx[-1])
                self._get_row(row_ptr, idx) = counter
                counter += 1
            

        # TODO: Enable the initialization with (data, col_ind, row_ptr) as a tuple
        else:
            raise NotImplementedError("Data type not yet supported")

    # Repetitous code, should be regrouped in the base class
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
        if not check_bound(self._shape, index):
            raise ValueError('Index ' + str(index) + ' out of border')

        row_start = self._row_ptr[index[:-1]]
        if row_start == []:
            return 0
        