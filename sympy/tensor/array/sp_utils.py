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

    An example of storage:
    A = [[1, 0, 0]
         [0, 2, 0]
         [0, 0, 3]]
    shapeA = (3, 3)
    With lil, we will have a multidimensional list '._data' containing all non-zero
    value in a row majored order and another list '._rows' containing column index of
    each non-zero value in each row.

    It is like we remove all zero from the initial array and only store the non-zero
    value and its column index.

    Examples
    ========

    >>> from sympy.tensor.array.sp_utils import lil
    >>> a = lil([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
    >>> a._data
    [[1], [2], [3]]
    >>> a._rows
    [[1], [2], [3]]
    >>> a[0, 1]
    1
    >>> a[0, 1] = 3
    >>> a.tolist()
    [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
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

    Note
    ====

    An example of storage:
    A = [[1, 0, 0]
         [0, 2, 0]
         [0, 0, 3]]
    Here we have a flatten list '._data' for all non-zero values. Another list '._coor'
    contains as many as the number of rank lists, which represent the coordinate(the
    tuple index) of the value.
    In this case:
    >>> sp_A = coo(A)
    >>> sp_A._data
    [1, 2, 3]
    >>> sp_A._coor
    [[0, 1, 2]
      0, 1, 2]]
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


class csr(Basic, SparseNDimArray):
    '''
    Compressed Sparse Row. This format is widely used. It has the following advantages:
    - Efficient item access and slicing
    - Fast arithmetics
    - Fast vector products

    Note
    ====

    An example of storage:
    A = [[1, 0, 0]
         [0, 2, 0]
         [0, 0, 3]]
    This format has 3 list. '._data' is a flatten list for all non-zero value. '.col_ind'
    has the column index of each non-zero value in their row. '._row_ptr' is the pointer
    of the first non-zero value in each row.
    The notion of row is the last dimension, which means that for an array of shape (a, b, c, d),
    it will have (a, b, c) rows. And the (a, b, c) is also flatten for the purpose of simplicity.

    In this case, we have:
    >>> sp_A._data
    [1, 2, 3]
    >>> sp_A._col_ind
    [0, 1, 2]
    >>> sp_A._row_ptr
    [0, 1, 2, 4]

    For the last list, there is one more element in the end(we have 3 rows but 4 pointers).
    This is a convention for csr format, where we store number_of_non_zero_value+1 in the end
    of the list, so that the algorithm can handle the operation for the last row.
    In this case, last value 4 in '._row_ptr' represents 3+1, that is to say that there are 3
    non-zero values in this array.
    '''
    def __new__(cls, iterable=None, shape=None, **kwargs):
        shape, flat_list = cls._handle_ndarray_creation_inputs(iterable, shape, **kwargs)
        shape = Tuple.fromiter(_sympify(i) for i in shape)
        cls._check_special_bounds(flat_list, shape)
        loop_size = functools.reduce(lambda x,y: x*y, shape) if shape else 0

        data = []
        col_ind = []
        row_ptr = []

        self = Basic.__new__(cls, data, col_ind, row_ptr, shape, **kwargs)
        self._data = data
        self._col_ind = col_ind
        self._row_ptr = row_ptr
        self._shape = shape
        self._rank = len(shape)
        self._loop_size = loop_size

        if isinstance(iterable, Iterable):
            iterable = SparseNDimArray(iterable)

        if isinstance(iterable, SparseNDimArray):
            current_row = 0
            for i, (k, v) in enumerate(sorted(iterable._sparse_array.items())):
                data.append(v)
                idx = self._get_tuple_index(k)
                col_ind.append(idx[-1])
                row = self.calculate_integer_index(idx[:-1])
                if row == current_row:
                    self._row_ptr.append(i)
                    current_row += 1
                elif row > current_row:
                    last_value = self._row_ptr[-1]
                    for j in range(current_row, row):
                        self._row_ptr.append(last_value)
                    self._row_ptr.append(i)
                    current_row = row + 1
            self._row_ptr.append(len(iterable._sparse_array) + 1)

        # TODO: Enable the initialization with (data, col_ind, row_ptr) as a tuple
        else:
            raise NotImplementedError("Data type not yet supported")

        return self

    def calculate_integer_index(self, tuple_index):
        integer_idx = 0
        for i, idx in enumerate(tuple_index):
            integer_idx = integer_idx*self._shape[i] + tuple_index[i]
        return integer_idx

    def calculate_tuple_index(self, integer_idx):
        index = []
        for i, sh in enumerate(reversed(self._shape[:-1])):
            index.append(integer_index % sh)
            integer_index //= sh
        index.reverse()
        return tuple(index)

    def __getitem__(self, index):
        if not isinstance(index, (tuple, Tuple)):
            index = self._get_tuple_index(index)
        if not check_bound(self._shape, index):
            raise ValueError('Index ' + str(index) + ' out of border')

        row_idx = self.calculate_integer_index(index[:-1])
        row_start = self._row_ptr[row_idx]
        row_end = self._row_ptr[row_idx + 1]
        row_values = self._data[row_start:row_end]

        col_idx_start = self._row_ptr[row_idx]
        col_idx_end = self._row_ptr[row_idx + 1]
        col_idx = list(self._col_ind[col_idx_start:col_idx_end])

        if index[-1] in col_idx:
            value_index = col_idx.index(index[-1])
            return row_values[value_index]
        else:
            return 0

    def __iter__(self):
        def iterator():
            for i in range(self._loop_size):
                yield self[i]
        return iterator()
