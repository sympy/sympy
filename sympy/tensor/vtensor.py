from sympy.core import sympify
from sympy.core.basic import Basic
from sympy.core.numbers import Rational
from sympy.core.symbol import symbols
from sympy.tensor.multiarray import MultiArray
from sympy.tensor.tensor import TensorHead, TensorIndexType, \
    TensorType, TensMul, tensorsymmetry, TensAdd, TensExpr
from sympy.matrices import Matrix


def _contract_multiarray(free1, free2, ma):
    farray1 = [i for i in free1]
    farray2 = [i for i in free2]

    for freeidx1 in free1:
        for freeidx2 in free2:
            # we don't need to check wrong indices,
            # it has already been done by TensMul.__mul__
            if freeidx1[0] == -freeidx2[0]:
                pos1 = farray1.index(freeidx1)
                pos2 = farray2.index(freeidx2)
                ma = ma.contract_positions(
                    pos1,
                    pos2 + len(farray1),
                    freeidx1[0]._tensortype._vmetric
                )
                farray1.pop(pos1)
                farray2.pop(pos2)
    return ma


class VTensorHead(Basic):
    """
    This class is analogous to TensorHead, except that it adds indexed data.

    Indexed data are passed with parameter ``values`` upon construction,
    and are stored as a MultiArray inside of a VTensorHead instance.

    Examples
    ========

    >>> from sympy.tensor.vtensor import VTensorHead, VTensorIndexType, vtensorhead
    >>> from sympy.tensor.tensor import tensor_indices
    >>> from sympy import ones
    >>> Lorentz = VTensorIndexType('Lorentz', metric=[1, -1, -1, -1], dummy_fmt='L')
    >>> i0, i1 = tensor_indices('i0:2', Lorentz)
    >>> A = vtensorhead('A', [Lorentz] * 2, [[1], [1]], values=ones(4, 4))

    in order to retrieve data, it is also necessary to specify abstract indices
    enclosed by round brackets, then numerical indices inside square brackets.

    >>> A(i0, i1)[0, 0]
    1

    Notice that square brackets create a valued tensor expression instance:

    >>> A(i0, i1)
    A(i0, i1)
    """

    def __new__(cls, tensorhead, multiarray, **kwargs):
        obj = Basic.__new__(cls, tensorhead, multiarray)

        obj._multiarray = multiarray
        obj._tensorhead = tensorhead
        # tensorhead contains information about symmetry of elements in
        # multiarray, as of now multiarray data are not checked
        # to possess the symmetries specified in tensorhead,
        # so it is up to the user to check that data are consistent.
        return obj

    def __call__(self, *args, **kwargs):
        th = self._tensorhead(*args, **kwargs)
        if th.rank != self._multiarray.rank:
            # contract_positions = []
            free1 = []
            free2 = []
            if not isinstance(args, tuple):
                argsli = (args,)
            else:
                argsli = args
            for i, el1 in enumerate(argsli):
                for j, el2 in enumerate(argsli[i:]):
                    if el1 == -el2:
                        free1.append((el1, i))
                        free2.append((el2, j))
            ma = _contract_multiarray(free1, free2, self._multiarray)
            if ma.rank == 0:
                return ma[(0,)]
        else:
            ma = self._multiarray
        return VTensMul(th, ma)

    def __pow__(self, other):
        metrics = [_._vmetric for _ in self._tensorhead.args[1]._args[0]]
        marray = self._multiarray
        for i, metric in enumerate(metrics):
            marray = marray.self_contract(i, metric)
        pow2 = marray[()]

        return pow2 ** (Rational(1, 2) * other)

    def __iter__(self):
        return self._multiarray.__iter__()

    def __str__(self):
        return self._tensorhead.__str__()

    def __repr__(self):
        return self._tensorhead.__repr__()


def vtensorhead(name, typ, sym, values, comm=0):
    """
    Function generating vtensorhead(s).

    Parameters
    ==========

    name : name or sequence of names (as in ``symbol``)

    typ :  index types

    sym :  same as ``*args`` in ``tensorsymmetry``

    comm : commutation group number
    see ``_TensorManager.set_comm``

    values : the metric used to raise and lower indices,
    can be a tuple, list, matrix or MultiArray.

    Notes
    =====

    The passed symmetry information is not checked to be consistent
    with the symmetry structure of the passed data,
    so it is up to the user to check that data are consistent.

    There could be inconsistencies, such as

    >>> inconsistent = vtensorhead('A', [Lorentz]*2, [[1, 1]], [
    ...     [1, 2, 3, 4],
    ...     [1, 2, 3, 4],
    ...     [0, 0, 0, 0],
    ...     [1, 1, 1, 1],
    ... ])

    In this case the symmetry information shows that the indices
    are symmetric, i.e. ``[[1, 1]]``, while the matrix isn't.

    Examples
    ========

    >>> from sympy.tensor.vtensor import TensorIndexType, vtensorhead
    >>> from sympy.tensor.vtensor import VTensorIndexType
    >>> from sympy.tensor.tensor import tensor_indices
    >>> from sympy import eye
    >>> Lorentz = VTensorIndexType('Lorentz', [1, -1, -1, -1], dummy_fmt='L')
    >>> a, b = tensor_indices('a,b', Lorentz)
    >>> A = vtensorhead('A', [Lorentz]*2, [[1]*2], eye(4))
    >>> A(a, -b)
    A(a, -b)

    """
    sym = tensorsymmetry(*sym)
    S = TensorType(typ, sym)
    mvalues = MultiArray.create(values)
    if isinstance(name, (str, unicode)):
        names = [x.name for x in symbols(name, seq=True)]
    else:
        raise ValueError('expecting a string')
    if len(names) == 1:
        thead = TensorHead(name, S, comm)
        return VTensorHead(thead, mvalues)  # TensorHead(names[0], self, comm)
    else:
        return [VTensorHead(TensorHead(name_iter, S, comm), mvalues) for name_iter in names]


class VTensExpr(Basic):
    """
    Abstract base class for valued tensor expressions. This is an analogous to
    the unvalued ``TensExpr``.

    A tensor expression can be a ``VTensAdd`` or ``VTensMul`` object.
    """
    is_commutative = False
    _op_priority = 12.0

    def __new__(cls, tensexpr, multiarray, **kwargs):
        obj = Basic.__new__(cls, tensexpr, multiarray, **kwargs)
        obj._tensexpr = tensexpr
        obj._multiarray = multiarray
        return obj

    @property
    def abstract(self):
        """
        Returns the corresponding tensor expression stripped of values.
        """
        return self._tensexpr

    @property
    def rank(self):
        """
        Returns the tensor rank, i.e. the number of free indices.
        """
        return self._tensexpr.rank

    def __getitem__(self, item):
        if not isinstance(item, (tuple, list,)):
            item = (item,)

        if isinstance(self, VTensAdd):
            free_indices = [_[0] for _ in self._tensexpr._args[0].free]
        else:
            free_indices = [_[0] for _ in self._tensexpr.free]
        free_indices_types = [_._args[1] for _ in free_indices]
        free_indices_contravariant = [_._args[2] for _ in free_indices]
        # vtensorhead('', [] * 2, [[1], [1]], values=ma)
        def get_metric_covariant_values_right(pointer):
            current_pos = len(pointer)

            if current_pos >= self.rank:
                return self._multiarray[pointer]

            if free_indices_contravariant[current_pos]:
                return get_metric_covariant_values_right(pointer + [item[current_pos]])

            sumt = 0
            vmetric = free_indices_types[current_pos]._vmetric
            # TODO: verify if metric is diagonal, so to avoid this for-loop
            if vmetric.rank == 1:
                sumt = vmetric[item[current_pos]] * get_metric_covariant_values_right(pointer + [item[current_pos]])
            elif vmetric.rank == 2:
                for i in xrange(vmetric.dimensions[0]):
                    sumt += vmetric[item[current_pos], i] * get_metric_covariant_values_right(pointer + [i])
            return sumt

#         all_contravariant_value = self._multiarray[item]
#         for pos, metric in metric_list:
#             all_contravariant_value
        return get_metric_covariant_values_right([])

    def __setitem__(self, key, value):
        raise NotImplementedError()
        # self._multiarray[key] = value

    def __str__(self):
        return self._tensexpr.__str__()

    def _pretty(self):
        return self._tensexpr._pretty()

    def __repr__(self):
        return self._tensexpr.__repr__()

    def __iter__(self):
        return self._multiarray.__iter__()

    def __call__(self, *args, **kwargs):
        vte = VTensExpr(self._tensexpr(*args, **kwargs), self._multiarray)
        vte.__class__ = self.__class__
        return vte

    def __pow__(self, other):
        """
            Perform a self-contraction on all free indices.
        """
        free = self._tensexpr._free

        marray = self._multiarray
        for i, free_el in enumerate(free):
            marray = marray.self_contract(i, free_el[0].args[1]._vmetric)
        return marray[()]

    def get_matrix(self):
        if 0 < self.rank <= 2:
            rows = self._multiarray.dimensions[0]
            columns = self._multiarray.dimensions[1] if self.rank == 2 else 1
            if self.rank == 2:
                mat_list = [] * rows
                for i in xrange(rows):
                    mat_list.append([])
                    for j in xrange(columns):
                        mat_list[i].append(self[i, j])
            else:
                mat_list = [None] * rows
                for i in xrange(rows):
                    mat_list[i] = self[i]
            return Matrix(mat_list)
        else:
            raise NotImplementedError(
                "missing multidimensional reduction to matrix.")


class VTensAdd(VTensExpr):
    """
    Sum of valued tensor expressions.

    Attributes
    ==========

    ``args`` : tuple of addends
    ``rank`` : rank of the tensor

    Examples
    ========

    >>> from sympy.tensor.vtensor import VTensorIndexType, vtensorhead
    >>> from sympy.tensor.tensor import tensor_indices
    >>> from sympy import eye
    >>> Lorentz = VTensorIndexType('Lorentz', [1, -1, -1, -1], dummy_fmt='L')
    >>> a, b = tensor_indices('a, b', Lorentz)
    >>> p, q = vtensorhead('p, q', [Lorentz], [[1]], [2, 3, -2, 7])
    >>> t = p(a) + q(a); t
    p(a) + q(a)
    >>> t(b)
    p(b) + q(b)
    """

    def __new__(cls, *args, **kwargs):
        sympified_args = [sympify(x) for x in args if x]
        obj = VTensExpr.__new__(cls, *sympified_args, **kwargs)

        args_tensexpr = []
        for i in sympified_args:
            if isinstance(i, VTensExpr):
                args_tensexpr.append(i)
            elif isinstance(i, TensExpr):
                raise ValueError("Mixed VTensExpr and TensExpr not allowed.")

        tensexpr = TensAdd(*[i._tensexpr for i in args_tensexpr], **kwargs)

        if not len(args_tensexpr):
            return sum(sympified_args)

        marray = MultiArray.create(args_tensexpr[0]._multiarray)

        for argument in sympified_args[1:]:
            marray += argument._multiarray

        obj._multiarray = marray
        obj._tensexpr = tensexpr
        return obj

    def __add__(self, other):
        return VTensAdd(self, other)

        if isinstance(other, VTensExpr):
            mr = self._tensexpr + other._tensexpr
        elif isinstance(other, TensExpr):
            raise ValueError("Mixed VTensExpr and TensExpr not allowed.")
        else:
            mr = self._tensexpr + other
        mr

    def __radd__(self, other):
        return VTensAdd(other, self)

    def __sub__(self, other):
        return VTensAdd(self, -other)

    def __rsub__(self, other):
        return VTensAdd(other, -self)

    def __mul__(self, other):
        return VTensAdd(*[x * other for x in self.args])

    def __div__(self, other):
        if isinstance(other, TensExpr):
            raise ValueError('cannot divide by a tensor')
        return VTensAdd(*[x / other for x in self.args])

    def __rdiv__(self, other):
        raise ValueError('cannot divide by a tensor')

    def __str__(self):
        return self._tensexpr.__str__()

    def applyfunc(self, func):
        vta = VTensAdd(*self.args)
        vta._multiarray = vta._multiarray.applyfunc(func)
        return vta


class VTensMul(VTensExpr):
    """
    Product of valued tensors

    Parameters
    ==========

    ``tensmul`` : ``TensMul`` expression.
    ``multiarray`` : a ``MultiArray`` instance containing indexed data.

    Attributes
    ==========

    ``args`` : Basic args specifying this class.
    ``abstract`` : The abstract tensor stripped of all stored values.
    ``multiarray`` : the MultiArray object containing stored data.

    """

    def __new__(cls, tensmul, multiarray):
        autodrop = True

        obj = VTensExpr.__new__(cls, tensmul, multiarray)

        obj._tensexpr = tensmul

        # TODO: autodrop on rank == 0 ? Decide
        obj._autodrop = autodrop

        free = obj._tensexpr._free
        if len(free) != multiarray.rank:
            raise ValueError("MultiArray rank does not equal the number of free indices")

        obj._multiarray = multiarray
        return obj

    def __mul__(self, other):
        if isinstance(other, VTensMul):
            # After abstract index contraction in superclass,
            # perform the numerical index contraction:
            mr = self._tensexpr * other._tensexpr
            ma = self._multiarray.tensor_product(other._multiarray)

            free1 = self._tensexpr._free
            free2 = other._tensexpr._free

            ma = _contract_multiarray(free1, free2, ma)

            if ma.rank == 0:  # self._autodrop ??? (care for creation)
                return ma[0]
            mr._multiarray = ma
        elif isinstance(other, TensMul):
            raise TypeError("Need a VTensMul, not a TensMul.")
        else:
            mr = self._tensexpr * other
            ma = MultiArray(self._multiarray) * other

        return VTensMul(
            mr,
            ma
        )

    def __rmul__(self, other):
        if isinstance(other, TensExpr):
            raise ValueError('cannot mix VTensExpr with TensExpr.')

        if isinstance(other, VTensMul):
            mr = other._tensexpr * self._tensexpr
            ma = other._multiarray.tensor_product(self._multiarray)
        else:
            mr = other * self._tensexpr
            ma = other * self._multiarray
        return VTensMul(mr, ma)

    def __add__(self, other):
        return VTensAdd(self, other)

    def __radd__(self, other):
        return TensAdd(other, self)

    def __div__(self, other):
        if isinstance(other, VTensExpr):
            raise ValueError("Cannot divide by a tensor.")

        mr = self._tensexpr / other
        ma = self._multiarray / other
        return VTensMul(mr, ma)

    def __rdiv__(self, other):
        mr = other / self._tensexpr
        # TODO: should it fail here?
        ma = other / self._multiarray
        return VTensMul(mr, ma)

    def __sub__(self, other):
        return VTensAdd(self, -other)

    def __rsub__(self, other):
        return TensAdd(other, -self)

    def _pretty(self):
        return self._tensexpr._pretty()

    def __neg__(self):
        return -1 * self

    def applyfunc(self, func):
        return VTensMul(self._tensexpr, self._multiarray.applyfunc(func))


class VTensorIndexType(Basic):
    """
    An analogous to TensorIndexType, just with the addition of all values of
    the metric.

    The metric is then used to contract indices by standard Einstein summation.

    Parameters
    ==========

    name : name of the tensor type

    metric : metric object, can be a ``MultiArray``, a ``Matrix``, a list or
    a nested list. If the metric is diagonal it can just be a simple list
    containing the elements of the diagonal.

    eps_dim : dimension of the epsilon tensor

    dummy_fmt : name of the head of dummy indices

    Attributes
    ==========

    ``name``
    ``metric`` : the metric tensor
    ``delta`` : ``Kronecker delta``
    ``epsilon`` : the ``Levi-Civita epsilon`` tensor
    ``dim``
    ``dim_eps``
    ``dummy_fmt``

    Examples
    ========

    >>> from sympy.tensor.vtensor import VTensorIndexType
    >>> Lorentz = VTensorIndexType('Lorentz', [1, -1, -1, -1], dummy_fmt='L')
    >>> Lorentz
    Lorentz
    """

    def __new__(cls, name, metric, eps_dim=None,
                 dummy_fmt=None):
        if not isinstance(metric, MultiArray):
            metric = MultiArray.create(metric)
        if metric.rank == 2:
            if metric.dimensions[0] != metric.dimensions[1]:
                raise ValueError("Metric is not square.")
            dim = metric.dimensions[0]
        elif metric.rank == 1:
            dim = metric.dimensions[0]
        else:
            raise ValueError("Metric rank is greater than 2.")
        # keep metric None
        tensorindextype = TensorIndexType(name, None, dim, eps_dim, dummy_fmt)

        obj = Basic.__new__(cls, tensorindextype, metric)
        # if not metric.is_square:
        #    raise ValueError("Metric is not square.")
        obj._vmetric = metric
        obj._tensorindextype = tensorindextype
        return obj

    @property
    def metric(self):
        return self._vmetric

    @property
    def name(self):
        return self._tensorindextype.name

    @property
    def args(self):
        return (self._tensorindextype, self.metric)

    @property
    def dummy_fmt(self):
        return self._tensorindextype.dummy_fmt

    @property
    def _dummy_fmt(self):
        return self._tensorindextype._dummy_fmt

    @property
    def _name(self):
        return self._tensorindextype._name

    def __str__(self):
        return self._tensorindextype.__str__()

    def _pretty(self):
        return self._tensorindextype._pretty()

    __repr__ = __str__
