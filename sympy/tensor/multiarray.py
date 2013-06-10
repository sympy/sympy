from sympy.core import sympify
from sympy.core.containers import Dict, Tuple
from sympy.core.expr import Expr
from sympy.matrices import Matrix
from sympy.core.singleton import S


class MultiArray(Expr):
    """
    Class to manage multidimensional arrays, it is a generalization of the
    Matrix class to any dimension.

    A MultiArray takes an optional string as first parameter upon construction,
    this is the name of the MultiArray object. It can be omitted, in which case
    the MultiArray will have no name.

    Examples
    ========

    There are many ways to create a MultiArray object, either specifying a
    Matrix, nested lists or tuples, or a function.

    >>> from sympy.tensor.multiarray import MultiArray
    >>> m = MultiArray.create([[2 ,3], [4, -1]])
    >>> m
    <MultiArray of rank 2, dimensions (2, 2)>

    ``m`` is two dimensional

    >>> print m.rank
    2

    Create a MultiArray from a Matrix element:

    >>> from sympy import Matrix
    >>> mat = Matrix([[1, 0], [2, 3]])
    >>> m2 = MultiArray.create(mat)

    get the matrix from MultiArray and check it is identical to the constructing one:

    >>> print m2.get_matrix() == mat
    True

    ``m2`` is two dimensional:

    >>> print m2.rank
    2

    Three dimensional MultiArray:

    >>> m3 = MultiArray.create([
    ...     [[173, 2], [3, 4], [5, 6]],
    ...     [[7, 8], [9, 10], [11, 12]]
    ... ])
    >>> print m3.rank
    3

    Access first element of three dimensional MultiArray:

    >>> print m3[0, 0, 0]
    173

    You can also store symbolic objects inside a MultiArray:

    >>> from sympy import var
    >>> x1, x2, x3 = var('x1 x2 x3')
    >>> m4 = MultiArray.create([x1, x2, x3])

    You can create sub-MultiArrays using slices as their indices,
    the same way as with a Matrix object:

    >>> m5 = m2[:, 1]
    >>> print m5.get_matrix() == mat[:, 1]
    True

    It is also possible to generate a MultiArray from a generating function,
    in this case MultiArray needs two parameters, the generating function
    and a list of dimensions (each corresponding to the range of the i-th index).
    The generating function takes a single parameter, which is a list
    of indices, to set the function result value as MultiArray value
    corresponding to those indices.

    >>> m6 = MultiArray.create(
    ...     lambda idx: -2 * (idx[0] + idx[1] + idx[2] + idx[3]),
    ...     (4, 5, 6, 7,)
    ... )
    >>> print m6[2, 3, 4, 5]
    -28

    Notes
    =====

    Elements can be accessed in various ways with the square bracket operator,
    i.e. [ ].

    Rank zero elements need [()] to access their unique element.

    Rank one elements may be given an integer (e.g. [3]), a tuple or a list
    (e.g. [(3,)]).

    Many indices may be passed both as tuples (e.g. [3, 2, 0]) or
    free indices (e.g. [(3, 2, 0)])
    """

    _op_priority = 12.0

    def __new__(cls, *args, **kwargs):
        if isinstance(args[0], (dict, Dict)):
            dict_data = args[0]
            dims = args[1]
            dict_data = Dict(dict_data)
            if not isinstance(dims, (tuple, list, Tuple,)):
                dims = (dims,)
            if not isinstance(dims, Tuple):
                dims = Tuple(*dims)
            obj = Expr.__new__(cls, dict_data, dims, **kwargs)
            from sympy.tensor.tensor import TensorSymmetry
            # TODO: detect symmetry from data
            # TODO: decide whether to use TensorSymmetry or something else
            obj._tensor_symmetry = TensorSymmetry.create([1]*obj.rank)
            return obj

        if isinstance(args[0], MultiArray):
            return Expr.__new__(cls, args[0].data_dict, args[0].dimensions, **kwargs)
        raise ValueError("Could not build MultiArray")

    @classmethod
    def create(cls, *args, **kwargs):
        # Dispatch parser according to various types of input:
        if hasattr(args[0], '__call__'):
            return MultiArray.from_function(args[0], args[1])
        elif isinstance(args[0], (list, tuple, Tuple,)):
            return MultiArray.from_nested(args[0])
        elif isinstance(args[0], Matrix):
            return MultiArray.from_matrix(args[0])
            # trailing_args = args[1:]
        elif isinstance(args[0], MultiArray):
            return MultiArray(args[0])

        raise ValueError("Unrecognized parameters.")

    @classmethod
    def from_matrix(cls, mat):
        rows, columns = mat.shape
        dimensions = Tuple(rows, columns,)

        temp_dict = dict()
        for i in xrange(rows):
            for j in xrange(columns):
                temp_dict[(i, j,)] = mat[i, j]

        return MultiArray(Dict(temp_dict), Tuple(*dimensions))

    @classmethod
    def from_nested(cls, nested):
        dimensions = []
        temp_dict = dict()
        try:
            mpointer = nested
            while True:
                dimensions.append(len(mpointer))
                mpointer = mpointer[0]
        except TypeError:
            pass
        current_position = []
        pointer = nested

        def tree_walker(pointer):
            if not isinstance(pointer, (tuple, list, Tuple,)):
                temp_dict[tuple(current_position)] = pointer
                return
            assert len(pointer) == dimensions[tree_walker.current_depth]
            for i, el in enumerate(pointer):
                current_position.append(i)
                tree_walker.current_depth += 1
                tree_walker(el)
                tree_walker.current_depth -= 1
                current_position.pop()
        tree_walker.current_depth = 0
        tree_walker(nested)

        return MultiArray(Dict(temp_dict), Tuple(*dimensions)) #, lcoeff, rcoeff)

    @classmethod
    def from_function(cls, func, dimensions):
        temp_dict = dict()

        def tree_walker_func(pointer):
            if tree_walker_func.depth == tree_walker_func.len_dimensions:
                temp_dict[tuple(pointer)] = func(pointer)
                return
            for i in xrange(dimensions[tree_walker_func.depth]):
                tree_walker_func.depth += 1
                tree_walker_func(pointer + [i])
                tree_walker_func.depth -= 1
        tree_walker_func.len_dimensions = len(dimensions)
        tree_walker_func.depth = 0
        tree_walker_func([])

        return MultiArray(Dict(temp_dict), Tuple(*dimensions))

    @property
    def data_dict(self):
        """
        Returns the basic constant dictionary-like object used to store data.
        """
        return self._args[0]

    @property
    def dimensions(self):
        """
        Returns the list of dimensions of the MultiArray's indices.
        """
        return self._args[1]

    def _populate_tensor_data(self, ranges, func, position, partly):
        if len(partly) == len(ranges):
            return func(partly)
        else:
            ret_val = []
            for i in range(0, 4):
                ret_val.append(
                    MultiArray.populate_tensor_data(
                        ranges,
                        func,
                        position,
                        partly + [i]
                    )
                )
            return ret_val

    @property
    def rank(self):
        """
            Returns the rank of the MultiArray, i.e. the number of indices
            needed to uniquely identify one of its internal members.
            The ``rank`` is sometimes called ``degree`` or ``order``.

            According to this definition of ``rank``, other objects can be
            given a ``rank`` value:
            A scalar has rank ``0``.
            A vector has rank ``1``.
            A matrix has rank ``2``.
        """
        lendim = len(self.dimensions)
        if lendim == 1 and self.dimensions[0] == 1:
            return 0
        return lendim

    @property
    def is_symmetric(self):
        return False
        # TODO: complete this step!!!

    def get_matrix(self):
        """
            If the rank of the MultiArray does not exceed the second order,
            it is possible to extract a MultiArray's data as a Matrix element.
        """
        # TODO: if rank == 0 or rank == 1 ???? Add in tests.
        if self.rank <= 2:
            if self.rank == 0:
                return Matrix([self[()]])
            rows = self.dimensions[0]
            columns = self.dimensions[1] if self.rank == 2 else 1
            if self.rank == 2:
                mat_list = [] * rows
                for i in xrange(rows):
                    mat_list.append([])
                    for j in xrange(columns):
                        mat_list[i].append(self.data_dict[(i, j)])
            else:
                mat_list = [None] * rows
                for i in xrange(rows):
                    mat_list[i] = self.data_dict[(i,)]
            return Matrix(mat_list)
        else:
            raise NotImplementedError(
                "missing multidimensional reduction to matrix.")

    def __str__(self):
        return "<MultiArray of rank %i, dimensions %s>" % (self.rank, str(self.dimensions),)

    def _pretty(self):
        if self.rank <= 2:
            return self.get_matrix()._pretty()
        return self.__str__()

    def __getitem__(self, item):
        """
        Get element uniquely identified by indices,
        otherwise get another MultiArray.

        If all indices are integers, this function returns get the corresponding element
        stored inside the MultiArray.

        If among the indices there are any slices, this function returns another
        MultiArray containing the data in the ranges defined by the slices.
        """
        other_array_types = (tuple, list, Tuple,)
        if not isinstance(item, other_array_types):
            item = [item, ]
        else:
            item = list(item)

        if isinstance(item, other_array_types):
            # check dimensions:
            new_dims = []
            basic_idx = []
            remap_idx = []
            for i, el in enumerate(item):
                idim = self.dimensions[i]
                basic_idx.append(None)
                if isinstance(el, slice):
                    sl = [el.start, el.stop, el.step]
                    # TODO: rewrite this code in order to have "step" working.
                    if sl[2]:
                        raise NotImplementedError("Unable to handle step in slice.")

                    if sl[0] is None:
                        sl[0] = 0
                    assert -idim <= sl[0] < idim
                    sl[0] %= idim

                    if sl[1] is None:
                        sl[1] = idim  # 0

                    assert -idim <= sl[1] <= idim
                    sl[1] = ((sl[1] - 1) % idim) + 1

                    new_dims.append(
#                         ((sl[1] % idim) - (sl[0] % idim)) % idim
                        (sl[1] - sl[0])
                    )
                    # ranges.append((i.start % idim, i.stop % idim, i.step))
                    remap_idx.append((i, sl[0],))
                    item[i] = slice(*sl)
                else:
                    assert -idim <= el < idim
                    basic_idx[-1] = el
                    # new_dims.append(1)
            for i in xrange(len(self.dimensions) - len(item)):
                new_dims.append(self.dimensions[i + len(item)])
                basic_idx.append(None)
                remap_idx.append((len(remap_idx) + 1, 0,))
            try:
                # TODO: create a flag instead of try-except.
                return self.data_dict[Tuple(*item)]
            except (KeyError, TypeError, AttributeError):
                def get_correct_element(x):
                    for i, el in enumerate(x):
                        basic_idx[remap_idx[i][0]] = remap_idx[i][1] + el
                    return self.data_dict[Tuple(*basic_idx)]
                return MultiArray.from_function(get_correct_element, new_dims)

        raise TypeError("Unrecognized type: {}".format(item))

    def get_indexwise_linear_transformation(self, *args):
        r"""
        Get a new MultiArray with same rank and same dimension as the original one,
        just with linear transformations applied on some of its indices.

        Given a transformation matrix `A_{ij}` acting on index `j`, this method performs the following operation:
        `\sum_j A_{ij} M_{\ldots, j, \ldots}`

        Examples
        ========

        >>> from sympy.tensor.multiarray import MultiArray
        >>> minkowski_metric = MultiArray.create([1, -1, -1, -1])
        >>> # TODO
        """
        trailing_ma = self
        for i, trma in enumerate(args):
            if trma in (None, False):
                continue
            if not isinstance(trma, MultiArray):
                trma = MultiArray.create(trma)
            if trma.rank == 1:
                if trma.dimensions[0] != self.dimensions[i]:
                    raise ValueError('Wrong dimension of index.')
            elif trma.rank == 2:
                if trma.dimensions[0] != trma.dimensions[1]:
                    raise ValueError('Transformation matrix is not square.')
                if trma.dimensions[0] != self.dimensions[i]:
                    raise ValueError('Wrong dimension of index.')
            else:
                raise ValueError('Transformation object has rank greater than two.')

            def fun(x):
                prevarg = [x[_] for _ in xrange(i)]
                postarg = [x[_] for _ in xrange(i + 1, self.rank)]
                summing = S.Zero
                if trma.rank == 2:
                    for j in xrange(self.dimensions[i]):
                        summing += trma[x[i], j] * trailing_ma[prevarg + [j] + postarg]

                elif trma.rank == 1:
                    summing = trma[x[i]] * trailing_ma[prevarg + [x[i]] + postarg]
                return summing

            trailing_ma = MultiArray.create(fun, self.dimensions)

        return trailing_ma

    def __setitem__(self, item, value):
        raise NotImplementedError()
        # TODO: should it create a new MultiArray or not be implemented at all?

    def _format_metric(self, metric, dim):
        if metric is None:
            metric = [1] * dim

        if not isinstance(metric, MultiArray):
            metric = MultiArray.create(metric)

        if metric.rank > 2:
            raise ValueError("Metric rank is greater than two.")

        return metric

    def self_extract(self, index_pos, index_value, metric = None):
        r"""
        This method extracts a rank `(n-1)` database from the rank `n` ``MultiArray``.

        Given a ``MultiArray`` `M_{i_1, \ldots, i_n}`, and a metric `L_{ij}`,
        an index position `P` and an index value `A`,
        this method performs an operation given by
        `\sum_k L_{A k} M_{i_1, \ldots, i_{P-1}, k, i_{P+1}, \ldots, i_n}`
        """
        if metric is None:
            metric = MultiArray.create([1] * self.dimensions[index_pos])
        params = [None] * index_pos + [metric] + [None] * (self.rank - index_pos - 1)
        pslice = [slice(None, None)] * index_pos + [index_value] + [slice(None, None)] * (self.rank - index_pos - 1)
        ma = self.get_indexwise_linear_transformation(*params)
        return ma[pslice]

    def self_contract(self, p1, metric=None):
        r"""
        This method extracts a rank `(n-1)` database from the rank `n` ``MultiArray``.

        Given a ``MultiArray`` `M_{i_1, \ldots, i_n}`, and a metric `L_{ij}`,
        an index position `P` and an index value `A`,
        this method performs an operation given by
        `\sum_{k, m} L_{k m} M_{i_1, \ldots, i_{P-1}, k, i_{P+1}, \ldots, i_n} M_{i_1, \ldots, i_{P-1}, m, i_{P+1}, \ldots, i_n}`
        """
        dim = self.dimensions[p1]
        metric = self._format_metric(metric, dim)

        def get_new_el(element):
            position = []
            for i, el in enumerate(element):
                position.append(el)
            position.insert(p1, None)
            summation = sympify(0)
            for i in xrange(dim):
                if metric.rank == 1:
                    position[p1] = i
                    summation += (self[position] ** 2) * metric[i]
                    continue
                for j in xrange(dim):
                    position[p1] = i
                    v1 = self[position]
                    position[p1] = j
                    v2 = self[position]
                    summation += v1 * metric[i, j] * v2
            return summation

        return MultiArray.create(get_new_el, [_ for dpos, _ in enumerate(self.dimensions) if dpos != p1])

    def contract_positions(self, p1, p2, metric=None):
        """
        Create a new MultiArray of rank n-2, with positions p1 and p2 summed over.

        Basically, given positions ``a`` and ``b`` representing some index,
        this function performs a summation
        `\sum_{i=0}^{d} \text{MultiArray}[\ldots, i, \ldots, i, \ldots]`,
        where the two ``i``s are placed in the positions of indices ``a`` and ``b``,
        and ``d`` stands for the dimension of ``a`` and ``b`` indices
        (they must have same dimension, otherwise a contraction is not possible).

        This function can be given an optional parameter ``metric``.
        This can be a list, Matrix, nested lists, or a MultiArray element.
        In any case its rank has to be one or two.
        In the case its rank is one, it will be considered as a diagonal matrix.

        The contraction with a metric is performed as
        `\sum_{i=0}^{d} \sum_{j=0}^d \text{MultiArrayInstance}[\ldots, i, \ldots, j, \ldots] g(i, j)`
        where $g(i, j)$ is the metric.

        Notice that when no metric is specified is the same as having $g(i, j)$
        to be a Kronecker delta, thus mirroring an Euclidean geometry.

        This function always returns a MultiArray. In case all indices are contracted,
        it will be a rank 0 instance of MultiArray, and its unique element
        is accessible as MultiArrayInstance[()]
        """
        dim1 = self.dimensions[p1]
        dim2 = self.dimensions[p2]

        if dim1 != dim2:
            raise ValueError("Not the same dimensions:"
                " Indices %i and %i." % (p1, p2,))

        metric = self._format_metric(metric, dim1)

        def get_new_el(element):
            position = []
            for i, el in enumerate(element):
                position.append(el)
            position.insert(p1, None)
            position.insert(p2, None)
            summation = sympify(0)
            for i in xrange(dim1):
                if metric.rank == 1:
                    position[p1] = i
                    position[p2] = i
                    summation += self[position] * metric[i]
                    continue
                for j in xrange(dim1):
                    position[p1] = j
                    position[p2] = i
                    summation += self[position] * metric[i, j]
            return summation

        if p1 == p2:
            raise ValueError("Cannot contract the same position.")
        if self.rank == 2:
            # TODO: decide whether this part is useful:
            summation = sympify(0)
            for i in range(dim1):
                if metric.rank == 1:
                    summation += self[(i, i)] * metric[i]
                    continue
                for j in range(dim1):
                    summation += self[(i, j)] * metric[i, j]
            return MultiArray.create((summation,))
        new_dims = [el for i, el in enumerate(self.dimensions) if i not in (p1, p2)]

        return MultiArray.from_function(get_new_el, new_dims)

    @property
    def tensor_symmetry(self):
        return self._tensor_symmetry

    def tensor_product(self, other_ma):
        """
        The tensor product create a new MultiArray by joining data from two factor MultiArray's.

        It is mathematically equivalent to:
        `C^{a_1 a_2 \cdots b_1 b_2 \cdots} = A^{a_1 a_2 \cdots} B^{b_1 b_2 \cdots}`
        The new MultiArray has all fields of the factors, returning their product.

        Usage
        =====

        >>> from sympy.tensor.multiarray import MultiArray
        >>> A = MultiArray.create([1, 2, 3])
        >>> B = MultiArray.create([2, -2, 0, 7])
        >>> tprod = A.tensor_product(B)
        >>> print tprod
        <MultiArray of rank 2, dimensions (3, 4)>

        >>> (A.tensor_product(B))[1, 2] == A[1]*B[2]
        True
        """
        if isinstance(other_ma, (tuple, list,)):
            other_ma = MultiArray.from_nested(other_ma)
        # TODO: join the multiarrays:
        # A[a,b,c], B[d, e]
        # ===> (A.tensor_product(B))[a, b, c, d, e] == A[a,b,c]*B[d,e]

        def get_tp_val(xval):
            return self[xval[:self.rank]] * other_ma[xval[self.rank:]]

        return MultiArray.from_function(get_tp_val, self.dimensions + other_ma.dimensions)

    def __mul__(self, other):
        # TODO: decide whether to implement as a tensor product
        if isinstance(other, MultiArray):
            raise NotImplementedError("Use tensor_product instead.")
        elif isinstance(other, Matrix):
            raise NotImplementedError("Use tensor_product instead.")
        else:
            return MultiArray.create(lambda x: self[x] * other, self.dimensions)

    def __rmul__(self, other):
        # TODO: implement
        if isinstance(other, MultiArray):
            raise NotImplementedError("Use tensor_product instead.")
        elif isinstance(other, Matrix):
            raise NotImplementedError("Use tensor_product instead.")
        else:
            return MultiArray.create(lambda x: other * self[x], self.dimensions)

    def __div__(self, other):
        if isinstance(other, (MultiArray, Matrix)):
            raise NotImplementedError()
        return MultiArray.create(lambda x: self[x] / other, self.dimensions)

    def __rdiv__(self, other):
        raise NotImplementedError()

    def __add__(self, other):
        """
        Addition of two MultiArray's: check that rank and dimensions are the same, then add elementwise.

        Addition of MultiArray and scalar:
        """
        if isinstance(other, MultiArray):
            if self.rank != other.rank:
                raise ValueError("rank is not the same")
            if self.dimensions != other.dimensions:
                raise ValueError("dimensions are not the same")
            return MultiArray.create(lambda x: self[x] + other[x], self.dimensions)
        return MultiArray.create(lambda x: self[x] + other, self.dimensions)

    def __radd__(self, other):
        if isinstance(other, MultiArray):
            if self.rank != other.rank:
                raise ValueError("rank is not the same")
            if self.dimensions != other.dimensions:
                raise ValueError("dimensions are not the same")
            return MultiArray.create(lambda x: other[x] + self[x], self.dimensions)
        return MultiArray.create(lambda x: other + self[x], self.dimensions)

    def __sub__(self, other):
        if isinstance(other, MultiArray):
            if self.rank != other.rank:
                raise ValueError("rank is not the same")
            if self.dimensions != other.dimensions:
                raise ValueError("dimensions are not the same")
            return MultiArray.create(lambda x: self[x] - other[x], self.dimensions)
        return MultiArray.create(lambda x: self[x] - other, self.dimensions)

    def __rsub__(self, other):
        if isinstance(other, MultiArray):
            if self.rank != other.rank:
                raise ValueError("rank is not the same")
            if self.dimensions != other.dimensions:
                raise ValueError("dimensions are not the same")
            return MultiArray.create(lambda x: other[x] - self[x], self.dimensions)
        return MultiArray.create(lambda x: other - self[x], self.dimensions)

    def applyfunc(self, func):
        """
        Applies a function to all elements of the MultiArray and returns a new
        MultiArray containing the resulting values.

        Examples
        ========

        >>> from sympy.tensor.multiarray import MultiArray
        >>> m = MultiArray.create([[2, 3], [4, -1]])
        >>> m.get_matrix()
        [2,  3]
        [4, -1]

        >>> m.applyfunc(lambda x: x**2).get_matrix()
        [ 4, 9]
        [16, 1]
        """
        return MultiArray.from_function(lambda x: func(self[x]), self.dimensions)

    def __iter__(self):
        key_list = sorted(self.data_dict.keys())

        for i in key_list:
            yield self.data_dict[i]

        raise StopIteration()


def multiempty(*args):
    marray = MultiArray.create(lambda x: 0, args)
    return marray


def multieye(rank, dim):
    marray = MultiArray.create(lambda x: 1 if len(set(x)) == 1 else 0, [dim] * rank)
    return marray


def multiones(rank, dim):
    marray = MultiArray.create(lambda x: 1, [dim] * rank)
    return marray
