"""
This module defines tensors with abstract index notation.

The abstract index notation has been first formalized by Penrose.

Tensor indices are formal objects, with a tensor type; there is no
notion of index range, it is only possible to assign the dimension,
used to trace the Kronecker delta; the dimension can be a Symbol.

The Einstein summation convention is used.
The covariant indices are indicated with a minus sign in front of the index.

For instance the tensor ``t = p(a)*A(b,c)*q(-c)`` has the index ``c``
contracted.

A tensor expression ``t`` can be called; called with its
indices in sorted order it is equal to itself:
in the above example ``t(a, b) == t``;
one can call ``t`` with different indices; ``t(c, d) == p(c)*A(d,a)*q(-a)``.

The contracted indices are dummy indices, internally they have no name,
the indices being represented by a graph-like structure.

Tensors are put in canonical form using ``canon_bp``, which uses
the Butler-Portugal algorithm for canonicalization using the monoterm
symmetries of the tensors.

If there is a (anti)symmetric metric, the indices can be raised and
lowered when the tensor is put in canonical form.
"""

from __future__ import print_function, division

import functools
from collections import defaultdict
import itertools
from sympy import Matrix, Rational
from sympy.assumptions.ask import ask, Q
from sympy.combinatorics.tensor_can import get_symmetric_group_sgs, \
    bsgs_direct_product, canonicalize, riemann_bsgs
from sympy.core import Basic, sympify, S
from sympy.core.add import Add
from sympy.core.compatibility import string_types
from sympy.core.containers import Tuple
from sympy.core.decorators import deprecated
from sympy.core.mul import Mul
from sympy.core.symbol import Symbol, symbols
from sympy.core.sympify import CantSympify
from sympy.external import import_module
from sympy.utilities.decorator import doctest_depends_on
from sympy.matrices import eye
from sympy.tensor.tids import TIDS, VTIDS


class _TensorDataLazyEvaluator(CantSympify):
    """
    EXPERIMENTAL: do not rely on this class, it may change without deprecation
    warnings in future versions of SymPy.

    This object contains the logic to associate components data to a tensor
    expression. Components data are set via the ``.data`` property of tensor
    expressions, is stored inside this class as a mapping between the tensor
    expression and the ``ndarray``.

    Computations are executed lazily: whereas the tensor expressions can have
    contractions, tensor products, and additions, components data are not
    computed until they are accessed by reading the ``.data`` property
    associated to the tensor expression.
    """
    _substitutions_dict = dict()
    # TODO: remove this or not?
    _substitutions_dict_tensmul = dict()

    def __getitem__(self, key):
        dat = self._get(key)
        if dat is None:
            return None

        numpy = import_module("numpy")
        if not isinstance(dat, numpy.ndarray):
            return dat

        if dat.ndim == 0:
            return dat[()]
        elif dat.ndim == 1 and dat.size == 1:
            return dat[0]
        return dat

    def _get(self, key):
        """
        Retrieve ``data`` associated with ``key``.

        This algorithm looks into ``self._substitutions_dict`` for all
        ``TensorHead`` in the ``TensExpr`` (or just ``TensorHead`` if key is a
        TensorHead instance). It reconstructs the components data that the
        tensor expression should have by performing on components data the
        operations that correspond to the abstract tensor operations applied.

        Metric tensor is handled in a different manner: it is pre-computed in
        ``self._substitutions_dict_tensmul``.
        """
        if key in self._substitutions_dict:
            return self._substitutions_dict[key]

        if isinstance(key, TensorHead):
            return None

        if isinstance(key, Tensor):
            # special case to handle metrics. Metric tensors cannot be
            # constructed through contraction by the metric, their
            # components show if they are a matrix or its inverse.
            signature = tuple([i.is_up for i in key.get_indices()])
            srch = (key.component,) + signature
            if srch in self._substitutions_dict_tensmul:
                return self._substitutions_dict_tensmul[srch]
            return self.data_tensmul_from_tensorhead(key, key.component)

        if isinstance(key, TensMul):
            # TODO: should TensMul.split() include non-tensor expressions?
            tensmul_list = key.split()
            if len(tensmul_list) == 1 and len(tensmul_list[0].components) == 1:
                # special case to handle metrics. Metric tensors cannot be
                # constructed through contraction by the metric, their
                # components show if they are a matrix or its inverse.
                signature = tuple([i.is_up for i in tensmul_list[0].get_indices()])
                srch = (tensmul_list[0].components[0],) + signature
                if srch in self._substitutions_dict_tensmul:
                    return self._substitutions_dict_tensmul[srch]
            data_list = [self.data_tensmul_from_tensorhead(i, i.components[0]) for i in tensmul_list]
            if all([i is None for i in data_list]):
                return None
            if any([i is None for i in data_list]):
                raise ValueError("Mixing tensors with associated components "\
                                 "data with tensors without components data")
            data_result, tensmul_result = self.data_product_tensors(data_list, tensmul_list)
            return data_result

        if isinstance(key, TensAdd):
            sumvar = S.Zero
            data_list = [i.data for i in key.args]
            if all([i is None for i in data_list]):
                return None
            if any([i is None for i in data_list]):
                raise ValueError("Mixing tensors with associated components "\
                                 "data with tensors without components data")
            for i in data_list:
                sumvar += i
            return sumvar

        return None

    def data_tensorhead_from_tensmul(self, data, tensmul, tensorhead):
        """
        This method is used when assigning components data to a ``TensMul``
        object, it converts components data to a fully contravariant ndarray,
        which is then stored according to the ``TensorHead`` key.
        """
        if data is None:
            return None

        return self._correct_signature_from_indices(
            data,
            tensmul.get_indices(),
            # TODO: remove
            #tensmul.free,
            tensmul.dummy_indices_list,
            True)

    def data_tensmul_from_tensorhead(self, tensmul, tensorhead):
        """
        This method corrects the components data to the right signature
        (covariant/contravariant) using the metric associated with each
        ``TensorIndexType``.
        """
        if tensorhead.data is None:
            return None

        return self._correct_signature_from_indices(
            tensorhead.data,
            tensmul.get_indices(),
            #tensmul.free,
            tensmul.dummy_indices_list)

    def data_product_tensors(self, data_list, tensmul_list):
        """
        Given a ``data_list``, list of ``ndarray``'s and a ``tensmul_list``,
        list of ``TensMul`` instances, compute the resulting ``ndarray``,
        after tensor products and contractions.
        """
        # TODO: should tensor be assigned data only if their index order can be determined?
        def data_mul(f, g):
            """
            Multiplies two ``ndarray`` objects, it first calls ``TIDS.mul``,
            then checks which indices have been contracted, and finally
            contraction operation on data, according to the contracted indices.
            """
            data1, tensmul1 = f
            data2, tensmul2 = g
            components, free, dum = TIDS.mul(tensmul1, tensmul2)
            data = _TensorDataLazyEvaluator._contract_ndarray(tensmul1.free_indices_list, tensmul2.free_indices_list, data1, data2)
            # TODO: do this more efficiently... maybe by just passing an index list
            # to .data_product_tensor(...)
            return data, tensmul1*tensmul2

        return functools.reduce(data_mul, zip(data_list, tensmul_list))

    def _assign_data_to_tensor_expr(self, key, data):
        if isinstance(key, TensAdd):
            raise ValueError('cannot assign data to TensAdd')
        # here it is assumed that `key` is a `TensMul` instance.
        if len(key.components) != 1:
            raise ValueError('cannot assign data to TensMul with multiple components')
        tensorhead = key.components[0]
        newdata = self.data_tensorhead_from_tensmul(data, key, tensorhead)
        return tensorhead, newdata

    def __setitem__(self, key, value):
        """
        Set the components data of a tensor object/expression.

        Components data are transformed to the all-contravariant form and stored
        with the corresponding ``TensorHead`` object. If a ``TensorHead`` object
        cannot be uniquely identified, it will raise an error.
        """
        data = _TensorDataLazyEvaluator.parse_data(value)

        # TensorHead and TensorIndexType can be assigned data directly, while
        # TensMul must first convert data to a fully contravariant form, and
        # assign it to its corresponding TensorHead single component.
        if not isinstance(key, (TensorHead, TensorIndexType)):
            key, data = self._assign_data_to_tensor_expr(key, data)

        if isinstance(key, TensorHead):
            for dim, indextype in zip(data.shape, key.index_types):
                if indextype.data is None:
                    raise ValueError("index type {} has no components data"\
                    " associated (needed to raise/lower index)".format(indextype))
                if indextype.dim is None:
                    continue
                if dim != indextype.dim:
                    raise ValueError("wrong dimension of ndarray")
        self._substitutions_dict[key] = data

    def __delitem__(self, key):
        del self._substitutions_dict[key]

    def __contains__(self, key):
        return key in self._substitutions_dict

    @staticmethod
    def _contract_ndarray(free1, free2, ndarray1, ndarray2):
        numpy = import_module('numpy')

        axes1 = []
        axes2 = []
        for jpos, jindex in enumerate(free2):
            if -jindex not in free1:
                continue
            nidx = free1.index(-jindex)
            axes1.append(nidx)
            axes2.append(jpos)

        contracted_ndarray = numpy.tensordot(
            ndarray1,
            ndarray2,
            (axes1, axes2)
        )
        return contracted_ndarray

    @staticmethod
    def add_tensor_mul(prod, f, g):
        def mul_function():
            return _TensorDataLazyEvaluator._contract_ndarray(f.free, g.free, f.data, g.data)

        _TensorDataLazyEvaluator._substitutions_dict[prod] = mul_function()

    @staticmethod
    def add_tensor_add(addition, f, g):
        def add_function():
            return f.data + g.data

        _TensorDataLazyEvaluator._substitutions_dict[addition] = add_function()

    def add_metric_data(self, metric, data):
        """
        Assign data to the ``metric`` tensor. The metric tensor behaves in an
        anomalous way when raising and lowering indices.

        A fully covariant metric is the inverse transpose of the fully
        contravariant metric (it is meant matrix inverse). If the metric is
        symmetric, the transpose is not necessary and mixed
        covariant/contravariant metrics are Kronecker deltas.
        """
        # hard assignment, data should not be added to `TensorHead` for metric:
        # the problem with `TensorHead` is that the metric is anomalous, i.e.
        # raising and lowering the index means considering the metric or its
        # inverse, this is not the case for other tensors.
        self._substitutions_dict_tensmul[metric, True, True] = data
        inverse_transpose = self.inverse_transpose_matrix(data)
        # in symmetric spaces, the traspose is the same as the original matrix,
        # the full covariant metric tensor is the inverse transpose, so this
        # code will be able to handle non-symmetric metrics.
        self._substitutions_dict_tensmul[metric, False, False] = inverse_transpose
        # now mixed cases, these are identical to the unit matrix if the metric
        # is symmetric.
        m = Matrix(data)
        invt = Matrix(inverse_transpose)
        self._substitutions_dict_tensmul[metric, True, False] = m * invt
        self._substitutions_dict_tensmul[metric, False, True] = invt * m
        #TODO test that with antisymmetric metric the metric contraction relations hold.

    @staticmethod
    def _flip_index_by_metric(data, metric, pos):
        numpy = import_module('numpy')

        data = numpy.tensordot(
                metric,
                data,
                (1, pos))
        return numpy.rollaxis(data, 0, pos+1)

    @staticmethod
    def inverse_matrix(ndarray):
        m = Matrix(ndarray).inv()
        return _TensorDataLazyEvaluator.parse_data(m)

    @staticmethod
    def inverse_transpose_matrix(ndarray):
        m = Matrix(ndarray).inv().T
        return _TensorDataLazyEvaluator.parse_data(m)

    @staticmethod
    def _correct_signature_from_indices(data, indices, dum, inverse=False):
        """
        Utility function to correct the values inside the components data
        ndarray according to whether indices are covariant or contravariant.

        It uses the metric matrix to lower values of covariant indices.
        """
        numpy = import_module('numpy')
        # change the ndarray values according covariantness/contravariantness of the indices
        # use the metric
        for i, indx in enumerate(indices):
            if not indx.is_up and not inverse:
                data = _TensorDataLazyEvaluator._flip_index_by_metric(data, indx._tensortype.data, i)
            elif not indx.is_up and inverse:
                data = _TensorDataLazyEvaluator._flip_index_by_metric(
                    data,
                    _TensorDataLazyEvaluator.inverse_matrix(indx._tensortype.data),
                    i
                )

        if len(dum) > 0:
            ### perform contractions ###
            axes1 = []
            axes2 = []
            for i, indx1 in enumerate(indices):
                try:
                    nd = indices[:i].index(-indx1)
                except ValueError:
                    continue
                axes1.append(nd)
                axes2.append(i)

            for ax1, ax2 in zip(axes1, axes2):
                data = numpy.trace(data, axis1=ax1, axis2=ax2)
        return data

    @staticmethod
    def _sort_data_axes(old, new):
        numpy = import_module('numpy')

        new_data = old.data.copy()

        old_free = [i[0] for i in old.free]
        new_free = [i[0] for i in new.free]

        for i in range(len(new_free)):
            for j in range(i, len(old_free)):
                if old_free[j] == new_free[i]:
                    old_free[i], old_free[j] = old_free[j], old_free[i]
                    new_data = numpy.swapaxes(new_data, i, j)
                    break
        return new_data

    @staticmethod
    def add_rearrange_tensmul_parts(new_tensmul, old_tensmul):
        def sorted_compo():
            return _TensorDataLazyEvaluator._sort_data_axes(old_tensmul, new_tensmul)

        _TensorDataLazyEvaluator._substitutions_dict[new_tensmul] = sorted_compo()

    @staticmethod
    @doctest_depends_on(modules=('numpy',))
    def parse_data(data):
        """
        Transform ``data`` to a numpy ndarray. The parameter ``data`` may
        contain data in various formats, e.g. nested lists, sympy ``Matrix``,
        and so on.

        Examples
        ========

        >>> from sympy.tensor.tensor import _TensorDataLazyEvaluator
        >>> _TensorDataLazyEvaluator.parse_data([1, 3, -6, 12])
        [1 3 -6 12]

        >>> _TensorDataLazyEvaluator.parse_data([[1, 2], [4, 7]])
        [[1 2]
         [4 7]]
        """
        numpy = import_module('numpy')

        if (numpy is not None) and (not isinstance(data, numpy.ndarray)):
            if len(data) == 2 and hasattr(data[0], '__call__'):

                def fromfunction_sympify(*x):
                    return sympify(data[0](*x))

                data = numpy.fromfunction(fromfunction_sympify, data[1])
            else:
                vsympify = numpy.vectorize(sympify)
                data = vsympify(numpy.array(data))
        return data

_tensor_data_substitution_dict = _TensorDataLazyEvaluator()


class _TensorManager(object):
    """
    Class to manage tensor properties.

    Notes
    =====

    Tensors belong to tensor commutation groups; each group has a label
    ``comm``; there are predefined labels:

    ``0``   tensors commuting with any other tensor

    ``1``   tensors anticommuting among themselves

    ``2``   tensors not commuting, apart with those with ``comm=0``

    Other groups can be defined using ``set_comm``; tensors in those
    groups commute with those with ``comm=0``; by default they
    do not commute with any other group.
    """
    def __init__(self):
        self._comm_init()

    def _comm_init(self):
        self._comm = [{} for i in range(3)]
        for i in range(3):
            self._comm[0][i] = 0
            self._comm[i][0] = 0
        self._comm[1][1] = 1
        self._comm[2][1] = None
        self._comm[1][2] = None
        self._comm_symbols2i = {0:0, 1:1, 2:2}
        self._comm_i2symbol = {0:0, 1:1, 2:2}

    @property
    def comm(self):
        return self._comm

    def comm_symbols2i(self, i):
        """
        get the commutation group number corresponding to ``i``

        ``i`` can be a symbol or a number or a string

        If ``i`` is not already defined its commutation group number
        is set.
        """
        if i not in self._comm_symbols2i:
            n = len(self._comm)
            self._comm.append({})
            self._comm[n][0] = 0
            self._comm[0][n] = 0
            self._comm_symbols2i[i] = n
            self._comm_i2symbol[n] = i
            return n
        return self._comm_symbols2i[i]

    def comm_i2symbol(self, i):
        """
        Returns the symbol corresponding to the commutation group number.
        """
        return self._comm_i2symbol[i]

    def set_comm(self, i, j, c):
        """
        set the commutation parameter ``c`` for commutation groups ``i, j``

        Parameters
        ==========

        i, j : symbols representing commutation groups

        c  :  group commutation number

        Notes
        =====

        ``i, j`` can be symbols, strings or numbers,
        apart from ``0, 1`` and ``2`` which are reserved respectively
        for commuting, anticommuting tensors and tensors not commuting
        with any other group apart with the commuting tensors.
        For the remaining cases, use this method to set the commutation rules;
        by default ``c=None``.

        The group commutation number ``c`` is assigned in correspondence
        to the group commutation symbols; it can be

        0        commuting

        1        anticommuting

        None     no commutation property

        Examples
        ========

        ``G`` and ``GH`` do not commute with themselves and commute with
        each other; A is commuting.

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead, TensorManager
        >>> Lorentz = TensorIndexType('Lorentz')
        >>> i0,i1,i2,i3,i4 = tensor_indices('i0:5', Lorentz)
        >>> A = tensorhead('A', [Lorentz], [[1]])
        >>> G = tensorhead('G', [Lorentz], [[1]], 'Gcomm')
        >>> GH = tensorhead('GH', [Lorentz], [[1]], 'GHcomm')
        >>> TensorManager.set_comm('Gcomm', 'GHcomm', 0)
        >>> (GH(i1)*G(i0)).canon_bp()
        G(i0)*GH(i1)
        >>> (G(i1)*G(i0)).canon_bp()
        G(i1)*G(i0)
        >>> (G(i1)*A(i0)).canon_bp()
        A(i0)*G(i1)
        """
        if c not in (0, 1, None):
            raise ValueError('`c` can assume only the values 0, 1 or None')

        if i not in self._comm_symbols2i:
            n = len(self._comm)
            self._comm.append({})
            self._comm[n][0] = 0
            self._comm[0][n] = 0
            self._comm_symbols2i[i] = n
            self._comm_i2symbol[n] = i
        if j not in self._comm_symbols2i:
            n = len(self._comm)
            self._comm.append({})
            self._comm[0][n] = 0
            self._comm[n][0] = 0
            self._comm_symbols2i[j] = n
            self._comm_i2symbol[n] = j
        ni = self._comm_symbols2i[i]
        nj = self._comm_symbols2i[j]
        self._comm[ni][nj] = c
        self._comm[nj][ni] = c

    def set_comms(self, *args):
        """
        set the commutation group numbers ``c`` for symbols ``i, j``

        Parameters
        ==========

        args : sequence of ``(i, j, c)``
        """
        for i, j, c in args:
            self.set_comm(i, j, c)

    def get_comm(self, i, j):
        """
        Return the commutation parameter for commutation group numbers ``i, j``

        see ``_TensorManager.set_comm``
        """
        return self._comm[i].get(j, 0 if i == 0 or j == 0 else None)

    def clear(self):
        """
        Clear the TensorManager.
        """
        self._comm_init()


TensorManager = _TensorManager()


# TODO: should the first argument of TensMul always be a Mul(...) containing its coefficient?

@doctest_depends_on(modules=('numpy',))
class TensorIndexType(Basic):
    """
    A TensorIndexType is characterized by its name and its metric.

    Parameters
    ==========

    name : name of the tensor type

    metric : metric symmetry or metric object or ``None``


    dim : dimension, it can be a symbol or an integer or ``None``

    eps_dim : dimension of the epsilon tensor

    dummy_fmt : name of the head of dummy indices

    Attributes
    ==========

    ``name``
    ``metric_name`` : it is 'metric' or metric.name
    ``metric_antisym``
    ``metric`` : the metric tensor
    ``delta`` : ``Kronecker delta``
    ``epsilon`` : the ``Levi-Civita epsilon`` tensor
    ``dim``
    ``dim_eps``
    ``dummy_fmt``
    ``data`` : a property to add ``ndarray`` values, to work in a specified basis.

    Notes
    =====

    The ``metric`` parameter can be:
    ``metric = False`` symmetric metric (in Riemannian geometry)

    ``metric = True`` antisymmetric metric (for spinor calculus)

    ``metric = None``  there is no metric

    ``metric`` can be an object having ``name`` and ``antisym`` attributes.


    If there is a metric the metric is used to raise and lower indices.

    In the case of antisymmetric metric, the following raising and
    lowering conventions will be adopted:

    ``psi(a) = g(a, b)*psi(-b); chi(-a) = chi(b)*g(-b, -a)``

    ``g(-a, b) = delta(-a, b); g(b, -a) = -delta(a, -b)``

    where ``delta(-a, b) = delta(b, -a)`` is the ``Kronecker delta``
    (see ``TensorIndex`` for the conventions on indices).

    If there is no metric it is not possible to raise or lower indices;
    e.g. the index of the defining representation of ``SU(N)``
    is 'covariant' and the conjugate representation is
    'contravariant'; for ``N > 2`` they are linearly independent.

    ``eps_dim`` is by default equal to ``dim``, if the latter is an integer;
    else it can be assigned (for use in naive dimensional regularization);
    if ``eps_dim`` is not an integer ``epsilon`` is ``None``.

    Examples
    ========

    >>> from sympy.tensor.tensor import TensorIndexType
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> Lorentz.metric
    metric(Lorentz,Lorentz)

    Examples with metric components data added, this means it is working on a
    fixed basis:

    >>> Lorentz.data = [1, -1, -1, -1]
    >>> Lorentz
    TensorIndexType(Lorentz, 0)
    >>> Lorentz.data
    [[1 0 0 0]
    [0 -1 0 0]
    [0 0 -1 0]
    [0 0 0 -1]]
    """

    def __new__(cls, name, metric=False, dim=None, eps_dim=None,
                dummy_fmt=None):

        if isinstance(name, string_types):
            name = Symbol(name)
        obj = Basic.__new__(cls, name, S.One if metric else S.Zero)
        obj._name = str(name)
        if not dummy_fmt:
            obj._dummy_fmt = '%s_%%d' % obj.name
        else:
            obj._dummy_fmt = '%s_%%d' % dummy_fmt
        if metric is None:
            obj.metric_antisym = None
            obj.metric = None
        else:
            if metric in (True, False, 0, 1):
                metric_name = 'metric'
                obj.metric_antisym = metric
            else:
                metric_name = metric.name
                obj.metric_antisym = metric.antisym
            sym2 = TensorSymmetry(get_symmetric_group_sgs(2, obj.metric_antisym))
            S2 = TensorType([obj]*2, sym2)
            obj.metric = S2(metric_name)
            obj.metric._matrix_behavior = True

        obj._dim = dim
        obj._delta = obj.get_kronecker_delta()
        obj._eps_dim = eps_dim if eps_dim else dim
        obj._epsilon = obj.get_epsilon()
        obj._autogenerated = []
        return obj

    @property
    def auto_right(self):
        if not hasattr(self, '_auto_right'):
            self._auto_right = BindableTensorIndex("auto_right", self)
        return self._auto_right

    @property
    def auto_left(self):
        if not hasattr(self, '_auto_left'):
            self._auto_left = BindableTensorIndex("auto_left", self)
        return self._auto_left

    @property
    def auto_index(self):
        if not hasattr(self, '_auto_index'):
            self._auto_index = BindableTensorIndex("auto_index", self)
        return self._auto_index

    @property
    def data(self):
        return _tensor_data_substitution_dict[self]

    @data.setter
    def data(self, data):
        numpy = import_module('numpy')
        data = _TensorDataLazyEvaluator.parse_data(data)
        if data.ndim > 2:
            raise ValueError("data have to be of rank 1 (diagonal metric) or 2.")
        if data.ndim == 1:
            if self.dim is not None:
                nda_dim = data.shape[0]
                if nda_dim != self.dim:
                    raise ValueError("Dimension mismatch")

            dim = data.shape[0]
            newndarray = numpy.zeros((dim, dim), dtype=object)
            for i, val in enumerate(data):
                newndarray[i, i] = val
            data = newndarray
        dim1, dim2 = data.shape
        if dim1 != dim2:
            raise ValueError("Non-square matrix tensor.")
        if self.dim is not None:
            if self.dim != dim1:
                raise ValueError("Dimension mismatch")
        _tensor_data_substitution_dict[self] = data
        _tensor_data_substitution_dict.add_metric_data(self.metric, data)
        delta = self.get_kronecker_delta()
        i1 = TensorIndex('i1', self)
        i2 = TensorIndex('i2', self)
        delta(i1, -i2).data = _TensorDataLazyEvaluator.parse_data(eye(dim1))

    @data.deleter
    def data(self):
        if self in _tensor_data_substitution_dict:
            del _tensor_data_substitution_dict[self]
        if self.metric in _tensor_data_substitution_dict:
            del _tensor_data_substitution_dict[self.metric]

    @property
    def name(self):
        return self._name

    @property
    def dim(self):
        return self._dim

    @property
    def delta(self):
        return self._delta

    @property
    def eps_dim(self):
        return self._eps_dim

    @property
    def epsilon(self):
        return self._epsilon

    @property
    def dummy_fmt(self):
        return self._dummy_fmt

    def get_kronecker_delta(self):
        sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
        S2 = TensorType([self]*2, sym2)
        delta = S2('KD')
        delta._matrix_behavior = True
        return delta

    def get_epsilon(self):
        if not isinstance(self._eps_dim, int):
            return None
        sym = TensorSymmetry(get_symmetric_group_sgs(self._eps_dim, 1))
        Sdim = TensorType([self]*self._eps_dim, sym)
        epsilon = Sdim('Eps')
        return epsilon

    def __lt__(self, other):
        return self.name < other.name

    def __str__(self):
        return self.name

    __repr__ = __str__

    def _components_data_full_destroy(self):
        """
        EXPERIMENTAL: do not rely on this API method.

        This destroys components data associated to the ``TensorIndexType``, if
        any, specifically:

        * metric tensor data
        * Kronecker tensor data
        """
        if self in _tensor_data_substitution_dict:
            del _tensor_data_substitution_dict[self]

        def delete_tensmul_data(key):
            if key in _tensor_data_substitution_dict._substitutions_dict_tensmul:
                del _tensor_data_substitution_dict._substitutions_dict_tensmul[key]

        # delete metric data:
        delete_tensmul_data((self.metric, True, True))
        delete_tensmul_data((self.metric, True, False))
        delete_tensmul_data((self.metric, False, True))
        delete_tensmul_data((self.metric, False, False))

        # delete delta tensor data:
        delta = self.get_kronecker_delta()
        if delta in _tensor_data_substitution_dict:
            del _tensor_data_substitution_dict[delta]


@doctest_depends_on(modules=('numpy',))
class TensorIndex(Basic):
    """
    Represents an abstract tensor index.

    Parameters
    ==========

    name : name of the index, or ``True`` if you want it to be automatically assigned
    tensortype : ``TensorIndexType`` of the index
    is_up :  flag for contravariant index

    Attributes
    ==========

    ``name``
    ``tensortype``
    ``is_up``

    Notes
    =====

    Tensor indices are contracted with the Einstein summation convention.

    An index can be in contravariant or in covariant form; in the latter
    case it is represented prepending a ``-`` to the index name.

    Dummy indices have a name with head given by ``tensortype.dummy_fmt``


    Examples
    ========

    >>> from sympy.tensor.tensor import TensorIndexType, TensorIndex, TensorSymmetry, TensorType, get_symmetric_group_sgs
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> i = TensorIndex('i', Lorentz); i
    i
    >>> sym1 = TensorSymmetry(*get_symmetric_group_sgs(1))
    >>> S1 = TensorType([Lorentz], sym1)
    >>> A, B = S1('A,B')
    >>> A(i)*B(-i)
    A(L_0)*B(-L_0)

    If you want the index name to be automatically assigned, just put ``True``
    in the ``name`` field, it will be generated using the reserved character
    ``_`` in front of its name, in order to avoid conflicts with possible
    existing indices:

    >>> i0 = TensorIndex(True, Lorentz)
    >>> i0
    _i0
    >>> i1 = TensorIndex(True, Lorentz)
    >>> i1
    _i1
    >>> A(i0)*B(-i1)
    A(_i0)*B(-_i1)
    >>> A(i0)*B(-i0)
    A(L_0)*B(-L_0)
    """
    def __new__(cls, name, tensortype, is_up=True):
        if isinstance(name, string_types):
            name_symbol = Symbol(name)
        elif isinstance(name, Symbol):
            name_symbol = name
        elif name is True:
            name = "_i{0}".format(len(tensortype._autogenerated))
            name_symbol = Symbol(name)
            tensortype._autogenerated.append(name_symbol)
        else:
            raise ValueError("invalid name")

        obj = Basic.__new__(cls, name_symbol, tensortype, S.One if is_up else S.Zero)
        obj._name = str(name)
        obj._tensortype = tensortype
        obj._is_up = is_up
        return obj

    @property
    def is_free(self):
        return self.args[3] == 1

    @property
    def name(self):
        return self._name

    @property
    def tensortype(self):
        return self._tensortype

    @property
    def is_up(self):
        return self._is_up

    def _print(self):
        s = self._name
        if not self._is_up:
            s = '-%s' % s
        return s

    def __lt__(self, other):
        return (self._tensortype, self._name) < (other._tensortype, other._name)

    def __neg__(self):
        t1 = self.func(self.args[0], self.args[1],
                (not self.args[2]))
        return t1


FreeTensorIndex = TensorIndex
#class FreeTensorIndex(TensorIndex):
#    # TODO: finish
#    def __new__(cls, *args, **kw_args):
#        return TensorIndex.__new__(cls, *args, **kw_args)


class DummyTensorIndex(TensorIndex):
    # TODO: finish
    def __new__(cls, *args, **kw_args):
        return TensorIndex.__new__(cls, *args, **kw_args)


class BindableTensorIndex(TensorIndex):
    # TODO: FINISH
    pass


def tensor_indices(s, typ):
    """
    Returns list of tensor indices given their names and their types

    Parameters
    ==========

    s : string of comma separated names of indices

    typ : list of ``TensorIndexType`` of the indices

    Examples
    ========

    >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> a, b, c, d = tensor_indices('a,b,c,d', Lorentz)
    """
    if isinstance(s, str):
        a = [x.name for x in symbols(s, seq=True)]
    else:
        raise ValueError('expecting a string')

    # TODO remove
    # tilist = [FreeTensorIndex(i, typ) for i in a]
    tilist = [TensorIndex(i, typ) for i in a]
    if len(tilist) == 1:
        return tilist[0]
    return tilist


@doctest_depends_on(modules=('numpy',))
class TensorSymmetry(Basic):
    """
    Monoterm symmetry of a tensor

    Parameters
    ==========

    bsgs : tuple ``(base, sgs)`` BSGS of the symmetry of the tensor

    Attributes
    ==========

    ``base`` : base of the BSGS
    ``generators`` : generators of the BSGS
    ``rank`` : rank of the tensor

    Notes
    =====

    A tensor can have an arbitrary monoterm symmetry provided by its BSGS.
    Multiterm symmetries, like the cyclic symmetry of the Riemann tensor,
    are not covered.

    See Also
    ========

    sympy.combinatorics.tensor_can.get_symmetric_group_sgs

    Examples
    ========

    Define a symmetric tensor

    >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, TensorType, get_symmetric_group_sgs
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> V = S2('V')
    """
    def __new__(cls, *args, **kw_args):
        if len(args) == 1:
            base, generators = args[0]
        elif len(args) == 2:
            base, generators = args
        else:
            raise TypeError("bsgs required, either two separate parameters or one tuple")

        if not isinstance(base, Tuple):
            base = Tuple(*base)
        if not isinstance(generators, Tuple):
            generators = Tuple(*generators)
        obj = Basic.__new__(cls, base, generators, **kw_args)
        return obj

    @property
    def base(self):
        return self.args[0]

    @property
    def generators(self):
        return self.args[1]

    @property
    def rank(self):
        return self.args[1][0].size - 2


def tensorsymmetry(*args):
    """
    Return a ``TensorSymmetry`` object.

    One can represent a tensor with any monoterm slot symmetry group
    using a BSGS.

    ``args`` can be a BSGS
    ``args[0]``    base
    ``args[1]``    sgs

    Usually tensors are in (direct products of) representations
    of the symmetric group;
    ``args`` can be a list of lists representing the shapes of Young tableaux

    Notes
    =====

    For instance:
    ``[[1]]``       vector
    ``[[1]*n]``     symmetric tensor of rank ``n``
    ``[[n]]``       antisymmetric tensor of rank ``n``
    ``[[2, 2]]``    monoterm slot symmetry of the Riemann tensor
    ``[[1],[1]]``   vector*vector
    ``[[2],[1],[1]`` (antisymmetric tensor)*vector*vector

    Notice that with the shape ``[2, 2]`` we associate only the monoterm
    symmetries of the Riemann tensor; this is an abuse of notation,
    since the shape ``[2, 2]`` corresponds usually to the irreducible
    representation characterized by the monoterm symmetries and by the
    cyclic symmetry.

    Examples
    ========

    Symmetric tensor using a Young tableau

    >>> from sympy.tensor.tensor import TensorIndexType, TensorType, tensorsymmetry
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> sym2 = tensorsymmetry([1, 1])
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> V = S2('V')

    Symmetric tensor using a ``BSGS`` (base, strong generator set)

    >>> from sympy.tensor.tensor import TensorSymmetry, get_symmetric_group_sgs
    >>> sym2 = tensorsymmetry(*get_symmetric_group_sgs(2))
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> V = S2('V')
    """
    from sympy.combinatorics import Permutation

    def tableau2bsgs(a):
        if len(a) == 1:
            # antisymmetric vector
            n = a[0]
            bsgs = get_symmetric_group_sgs(n, 1)
        else:
            if all(x == 1 for x in a):
                # symmetric vector
                n = len(a)
                bsgs = get_symmetric_group_sgs(n)
            elif a == [2, 2]:
                bsgs = riemann_bsgs
            else:
                raise NotImplementedError
        return bsgs

    if not args:
        return TensorSymmetry(Tuple(), Tuple(Permutation(1)))

    if len(args) == 2 and isinstance(args[1][0], Permutation):
        return TensorSymmetry(args)
    base, sgs = tableau2bsgs(args[0])
    for a in args[1:]:
        basex, sgsx = tableau2bsgs(a)
        base, sgs = bsgs_direct_product(base, sgs, basex, sgsx)
    return TensorSymmetry(Tuple(base, sgs))


@doctest_depends_on(modules=('numpy',))
class TensorType(Basic):
    """
    Class of tensor types.

    Parameters
    ==========

    index_types : list of ``TensorIndexType`` of the tensor indices
    symmetry : ``TensorSymmetry`` of the tensor

    Attributes
    ==========

    ``index_types``
    ``symmetry``
    ``types`` : list of ``TensorIndexType`` without repetitions

    Examples
    ========

    Define a symmetric tensor

    >>> from sympy.tensor.tensor import TensorIndexType, tensorsymmetry, TensorType
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> sym2 = tensorsymmetry([1, 1])
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> V = S2('V')
    """
    is_commutative = False

    def __new__(cls, index_types, symmetry, **kw_args):
        assert symmetry.rank == len(index_types)
        obj = Basic.__new__(cls, Tuple(*index_types), symmetry, **kw_args)
        return obj

    @property
    def index_types(self):
        return self.args[0]

    @property
    def symmetry(self):
        return self.args[1]

    @property
    def types(self):
        return sorted(set(self.index_types), key=lambda x: x.name)

    def __str__(self):
        return 'TensorType(%s)' % ([str(x) for x in self.index_types])

    def __call__(self, s, comm=0, matrix_behavior=0):
        """
        Return a TensorHead object or a list of TensorHead objects.

        ``s``  name or string of names

        ``comm``: commutation group number
        see ``_TensorManager.set_comm``

        Examples
        ========

        Define symmetric tensors ``V``, ``W`` and ``G``, respectively
        commuting, anticommuting and with no commutation symmetry

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorsymmetry, TensorType, canon_bp
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> a, b = tensor_indices('a,b', Lorentz)
        >>> sym2 = tensorsymmetry([1]*2)
        >>> S2 = TensorType([Lorentz]*2, sym2)
        >>> V = S2('V')
        >>> W = S2('W', 1)
        >>> G = S2('G', 2)
        >>> canon_bp(V(a, b)*V(-b, -a))
        V(L_0, L_1)*V(-L_0, -L_1)
        >>> canon_bp(W(a, b)*W(-b, -a))
        0
        """
        if isinstance(s, str):
            names = [x.name for x in symbols(s, seq=True)]
        else:
            raise ValueError('expecting a string')
        if len(names) == 1:
            return TensorHead(names[0], self, comm, matrix_behavior=matrix_behavior)
        else:
            return [TensorHead(name, self, comm, matrix_behavior=matrix_behavior) for name in names]


def tensorhead(name, typ, sym, comm=0, matrix_behavior=0):
    """
    Function generating tensorhead(s).

    Parameters
    ==========

    name : name or sequence of names (as in ``symbol``)

    typ :  index types

    sym :  same as ``*args`` in ``tensorsymmetry``

    comm : commutation group number
    see ``_TensorManager.set_comm``


    Examples
    ========

    >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> a, b = tensor_indices('a,b', Lorentz)
    >>> A = tensorhead('A', [Lorentz]*2, [[1]*2])
    >>> A(a, -b)
    A(a, -b)

    """
    sym = tensorsymmetry(*sym)
    S = TensorType(typ, sym)
    th = S(name, comm, matrix_behavior=matrix_behavior)
    return th


@doctest_depends_on(modules=('numpy',))
class TensorHead(Basic):
    r"""
    Tensor head of the tensor

    Parameters
    ==========

    name : name of the tensor

    typ : list of TensorIndexType

    comm : commutation group number

    Attributes
    ==========

    ``name``
    ``index_types``
    ``rank``
    ``types``  :  equal to ``typ.types``
    ``symmetry`` : equal to ``typ.symmetry``
    ``comm`` : commutation group

    Notes
    =====

    A ``TensorHead`` belongs to a commutation group, defined by a
    symbol on number ``comm`` (see ``_TensorManager.set_comm``);
    tensors in a commutation group have the same commutation properties;
    by default ``comm`` is ``0``, the group of the commuting tensors.

    Examples
    ========

    >>> from sympy.tensor.tensor import TensorIndexType, tensorsymmetry, TensorType
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> sym2 = tensorsymmetry([1], [1])
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> A = S2('A')

    Examples with ndarray values, the components data assigned to the
    ``TensorHead`` object are assumed to be in a fully-contravariant
    representation. In case it is necessary to assign components data which
    represents the values of a non-fully covariant tensor, see the other
    examples.

    >>> from sympy.tensor.tensor import tensor_indices, tensorhead
    >>> Lorentz.data = [1, -1, -1, -1]
    >>> i0, i1 = tensor_indices('i0:2', Lorentz)
    >>> A.data = [[j+2*i for j in range(4)] for i in range(4)]

    in order to retrieve data, it is also necessary to specify abstract indices
    enclosed by round brackets, then numerical indices inside square brackets.

    >>> A(i0, i1)[0, 0]
    0
    >>> A(i0, i1)[2, 3] == 3+2*2
    True

    Notice that square brackets create a valued tensor expression instance:

    >>> A(i0, i1)
    A(i0, i1)

    To view the data, just type:

    >>> A.data
    [[0 1 2 3]
     [2 3 4 5]
     [4 5 6 7]
     [6 7 8 9]]

    Turning to a tensor expression, covariant indices get the corresponding
    components data corrected by the metric:

    >>> A(i0, -i1).data
    [[0 -1 -2 -3]
     [2 -3 -4 -5]
     [4 -5 -6 -7]
     [6 -7 -8 -9]]

    >>> A(-i0, -i1).data
    [[0 -1 -2 -3]
     [-2 3 4 5]
     [-4 5 6 7]
     [-6 7 8 9]]

    while if all indices are contravariant, the ``ndarray`` remains the same

    >>> A(i0, i1).data
     [[0 1 2 3]
     [2 3 4 5]
     [4 5 6 7]
     [6 7 8 9]]

    When all indices are contracted and components data are added to the tensor,
    accessing the data will return a scalar, no numpy object. In fact, numpy
    ndarrays are dropped to scalars if they contain only one element.

    >>> A(i0, -i0)
    A(L_0, -L_0)
    >>> A(i0, -i0).data
    -18

    It is also possible to assign components data to an indexed tensor, i.e. a
    tensor with specified covariant and contravariant components. In this
    example, the covariant components data of the Electromagnetic tensor are
    injected into `A`:

    >>> from sympy import symbols
    >>> Ex, Ey, Ez, Bx, By, Bz = symbols('E_x E_y E_z B_x B_y B_z')
    >>> c = symbols('c', positive=True)
    >>> A(-i0, -i1).data = [
    ... [0, Ex/c, Ey/c, Ez/c],
    ... [-Ex/c, 0, -Bz, By],
    ... [-Ey/c, Bz, 0, -Bx],
    ... [-Ez/c, -By, Bx, 0]]

    Now it is possible to retrieve the contravariant form of the Electromagnetic
    tensor:

    >>> A(i0, i1).data
    [[0 -E_x/c -E_y/c -E_z/c]
     [E_x/c 0 -B_z B_y]
     [E_y/c B_z 0 -B_x]
     [E_z/c -B_y B_x 0]]

    and the mixed contravariant-covariant form:

    >>> A(i0, -i1).data
    [[0 E_x/c E_y/c E_z/c]
     [E_x/c 0 B_z -B_y]
     [E_y/c -B_z 0 B_x]
     [E_z/c B_y -B_x 0]]

    To convert the numpy's ndarray to a sympy matrix, just cast:

    >>> from sympy import Matrix
    >>> Matrix(A.data)
    Matrix([
    [    0, -E_x/c, -E_y/c, -E_z/c],
    [E_x/c,      0,   -B_z,    B_y],
    [E_y/c,    B_z,      0,   -B_x],
    [E_z/c,   -B_y,    B_x,      0]])

    Still notice, in this last example, that accessing components data from a
    tensor without specifying the indices is equivalent to assume that all
    indices are contravariant.

    It is also possible to store symbolic components data inside a tensor, for
    example, define a four-momentum-like tensor:

    >>> from sympy import symbols
    >>> P = tensorhead('P', [Lorentz], [[1]])
    >>> E, px, py, pz = symbols('E p_x p_y p_z', positive=True)
    >>> P.data = [E, px, py, pz]

    The contravariant and covariant components are, respectively:

    >>> P(i0).data
    [E p_x p_y p_z]
    >>> P(-i0).data
    [E -p_x -p_y -p_z]

    The contraction of a 1-index tensor by itself is usually indicated by a
    power by two:

    >>> P(i0)**2
    E**2 - p_x**2 - p_y**2 - p_z**2

    As the power by two is clearly identical to `P_\mu P^\mu`, it is possible to
    simply contract the ``TensorHead`` object, without specifying the indices

    >>> P**2
    E**2 - p_x**2 - p_y**2 - p_z**2
    """
    is_commutative = False

    def __new__(cls, name, typ, comm=0, matrix_behavior=0, **kw_args):
        if isinstance(name, string_types):
            name_symbol = Symbol(name)
        elif isinstance(name, Symbol):
            name_symbol = name
        else:
            raise ValueError("invalid name")

        comm2i = TensorManager.comm_symbols2i(comm)

        obj = Basic.__new__(cls, name_symbol, typ, **kw_args)

        obj._matrix_behavior = matrix_behavior

        obj._name = obj.args[0].name
        #obj._rank = len(obj.index_types)
        obj._types = typ.types
        obj._symmetry = typ.symmetry
        obj._comm = comm2i
        return obj

    @property
    def name(self):
        return self._name

    @property
    def rank(self):
        return len(self.index_types)

    @property
    def types(self):
        return self._types[:]

    @property
    def symmetry(self):
        return self._symmetry

    @property
    def typ(self):
        return self.args[1]

    @property
    def comm(self):
        return self._comm

    @property
    def index_types(self):
        return self.args[1].index_types[:]

    def __lt__(self, other):
        return (self.name, self.index_types) < (other.name, other.index_types)

    def commutes_with(self, other):
        """
        Returns ``0`` if ``self`` and ``other`` commute, ``1`` if they anticommute.

        Returns ``None`` if ``self`` and ``other`` neither commute nor anticommute.
        """
        r = TensorManager.get_comm(self._comm, other._comm)
        return r

    def _print(self):
        return '%s(%s)' %(self.name, ','.join([str(x) for x in self.index_types]))

#     def _check_auto_matrix_indices_in_call(self, *indices):
#         # TODO: remove?
#         matrix_behavior_kinds = dict()
# 
#         if len(indices) != len(self.index_types):
#             if not self._matrix_behavior:
#                 raise ValueError('wrong number of indices')
# 
#             # _matrix_behavior is True, so take the last one or two missing
#             # indices as auto-matrix indices:
#             ldiff = len(self.index_types) - len(indices)
#             if ldiff > 2:
#                 raise ValueError('wrong number of indices')
#             if ldiff == 2:
#                 mat_ind = [len(indices), len(indices) + 1]
#             elif ldiff == 1:
#                 mat_ind = [len(indices)]
#             not_equal = True
#         else:
#             not_equal = False
#             mat_ind = [i for i, e in enumerate(indices) if e is True]
#             if mat_ind:
#                 not_equal = True
#             indices = tuple([_ for _ in indices if _ is not True])
# 
#             for i, el in enumerate(indices):
#                 if not isinstance(el, TensorIndex):
#                     not_equal = True
#                     break
#                 if el._tensortype != self.index_types[i]:
#                     not_equal = True
#                     break
# 
#         if not_equal:
#             for el in mat_ind:
#                 eltyp = self.index_types[el]
#                 if eltyp in matrix_behavior_kinds:
#                     elind = -self.index_types[el].auto_right
#                     matrix_behavior_kinds[eltyp].append(elind)
#                 else:
#                     elind = self.index_types[el].auto_left
#                     matrix_behavior_kinds[eltyp] = [elind]
#                 indices = indices[:el] + (elind,) + indices[el:]
# 
#         return indices, matrix_behavior_kinds

    def __call__(self, *indices):
        """
        Returns a tensor with indices.

        There is a special behavior in case of indices denoted by ``True``,
        they are considered auto-matrix indices, their slots are automatically
        filled, and confer to the tensor the behavior of a matrix or vector
        upon multiplication with another tensor containing auto-matrix indices
        of the same ``TensorIndexType``. This means indices get summed over the
        same way as in matrix multiplication. For matrix behavior, define two
        auto-matrix indices, for vector behavior define just one.

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> a, b = tensor_indices('a,b', Lorentz)
        >>> A = tensorhead('A', [Lorentz]*2, [[1]*2])
        >>> t = A(a, -b)
        >>> t
        A(a, -b)

        To use the auto-matrix index behavior, just put a ``True`` on the
        desired index position.

        >>> r = A(True, True)
        >>> r
        A(auto_left, -auto_right)

        Here ``auto_left`` and ``auto_right`` are automatically generated
        tensor indices, they are only two for every ``TensorIndexType`` and
        can be assigned to just one or two indices of a given type.

        Auto-matrix indices can be assigned many times in a tensor, if indices
        are of different ``TensorIndexType``

        >>> Spinor = TensorIndexType('Spinor', dummy_fmt='S')
        >>> B = tensorhead('B', [Lorentz, Lorentz, Spinor, Spinor], [[1]*4])
        >>> s = B(True, True, True, True)
        >>> s
        B(auto_left, -auto_right, auto_left, -auto_right)

        Here, ``auto_left`` and ``auto_right`` are repeated twice, but they are
        not the same indices, as they refer to different ``TensorIndexType``s.

        Auto-matrix indices are automatically contracted upon multiplication,

        >>> r*s
        A(auto_left, -L_0)*B(L_0, -auto_right, auto_left, -auto_right)

        The multiplication algorithm has found an ``auto_right`` index in ``A``
        and an ``auto_left`` index in ``B`` referring to the same
        ``TensorIndexType`` (``Lorentz``), so they have been contracted.

        Auto-matrix indices can be accessed from the ``TensorIndexType``:

        >>> Lorentz.auto_right
        auto_right
        >>> Lorentz.auto_left
        auto_left

        There is a special case, in which the ``True`` parameter is not needed
        to declare an auto-matrix index, i.e. when the matrix behavior has been
        declared upon ``TensorHead`` construction, in that case the last one or
        two tensor indices may be omitted, so that they automatically become
        auto-matrix indices:

        >>> C = tensorhead('C', [Lorentz, Lorentz], [[1]*2], matrix_behavior=True)
        >>> C()
        C(auto_left, -auto_right)

        """

        #indices, matrix_behavior_kinds = self._check_auto_matrix_indices_in_call(*indices)

        components = [self]
        tmul = Tensor(self, indices)
        #tmul._matrix_behavior_kinds = matrix_behavior_kinds
        return tmul

#        tids = TIDS.from_components_and_indices(components, indices)
#
#        tmul = TensMul.from_TIDS(S.One, tids)
#        tmul._matrix_behavior_kinds = matrix_behavior_kinds
#
#        return tmul

    def __pow__(self, other):
        if self.data is None:
            raise ValueError("No power on abstract tensors.")
        numpy = import_module('numpy')
        metrics = [_.data for _ in self.args[1].args[0]]

        marray = self.data
        for metric in metrics:
            marray = numpy.tensordot(marray, numpy.tensordot(metric, marray, (1, 0)), (0, 0))
        pow2 = marray[()]
        return pow2 ** (Rational(1, 2) * other)

    @property
    def data(self):
        return _tensor_data_substitution_dict[self]

    @data.setter
    def data(self, data):
        _tensor_data_substitution_dict[self] = data

    @data.deleter
    def data(self):
        if self in _tensor_data_substitution_dict:
            del _tensor_data_substitution_dict[self]

    def __iter__(self):
        return self.data.flatten().__iter__()

    def _components_data_full_destroy(self):
        """
        EXPERIMENTAL: do not rely on this API method.

        Destroy components data associated to the ``TensorHead`` object, this
        checks for attached components data, and destroys components data too.
        """
        # do not garbage collect Kronecker tensor (it should be done by
        # ``TensorIndexType`` garbage collection)
        if self.name == "KD":
            return

        # the data attached to a tensor must be deleted only by the TensorHead
        # destructor. If the TensorHead is deleted, it means that there are no
        # more instances of that tensor anywhere.
        if self in _tensor_data_substitution_dict:
            del _tensor_data_substitution_dict[self]


@doctest_depends_on(modules=('numpy',))
class TensExpr(Basic):
    """
    Abstract base class for tensor expressions

    Notes
    =====

    A tensor expression is an expression formed by tensors;
    currently the sums of tensors are distributed.

    A ``TensExpr`` can be a ``TensAdd`` or a ``TensMul``.

    ``TensAdd`` objects are put in canonical form using the Butler-Portugal
    algorithm for canonicalization under monoterm symmetries.

    ``TensMul`` objects are formed by products of component tensors,
    and include a coefficient, which is a SymPy expression.


    In the internal representation contracted indices are represented
    by ``(ipos1, ipos2, icomp1, icomp2)``, where ``icomp1`` is the position
    of the component tensor with contravariant index, ``ipos1`` is the
    slot which the index occupies in that component tensor.

    Contracted indices are therefore nameless in the internal representation.
    """

    _op_priority = 11.0
    is_commutative = False

    @property
    def has_index_order(self):
        return False

    def __neg__(self):
        return self*S.NegativeOne

    def __abs__(self):
        raise NotImplementedError

    def __add__(self, other):
        raise NotImplementedError

    def __radd__(self, other):
        raise NotImplementedError

    def __sub__(self, other):
        raise NotImplementedError

    def __rsub__(self, other):
        raise NotImplementedError

    def __mul__(self, other):
        raise NotImplementedError

    def __rmul__(self, other):
        raise NotImplementedError

    def __pow__(self, other):
        if self.data is None:
            raise ValueError("No power without ndarray data.")
        numpy = import_module('numpy')
        free = self.free

        marray = self.data
        for metric in free:
            marray = numpy.tensordot(
                marray,
                numpy.tensordot(
                    metric[0]._tensortype.data,
                    marray,
                    (1, 0)
                ),
                (0, 0)
            )
        pow2 = marray[()]
        return pow2 ** (Rational(1, 2) * other)

    def __rpow__(self, other):
        raise NotImplementedError

    def __div__(self, other):
        raise NotImplementedError

    def __rdiv__(self, other):
        raise NotImplementedError()

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    @doctest_depends_on(modules=('numpy',))
    def get_matrix(self):
        """
        Returns ndarray components data as a matrix, if components data are
        available and ndarray dimension does not exceed 2.

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensorsymmetry, TensorType
        >>> from sympy import ones
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> sym2 = tensorsymmetry([1]*2)
        >>> S2 = TensorType([Lorentz]*2, sym2)
        >>> A = S2('A')

        The tensor ``A`` is symmetric in its indices, as can be deduced by the
        ``[1, 1]`` Young tableau when constructing `sym2`. One has to be
        careful to assign symmetric component data to ``A``, as the symmetry
        properties of data are currently not checked to be compatible with the
        defined tensor symmetry.

        >>> from sympy.tensor.tensor import tensor_indices, tensorhead
        >>> Lorentz.data = [1, -1, -1, -1]
        >>> i0, i1 = tensor_indices('i0:2', Lorentz)
        >>> A.data = [[j+i for j in range(4)] for i in range(4)]
        >>> A(i0, i1).get_matrix()
        Matrix([
        [0, 1, 2, 3],
        [1, 2, 3, 4],
        [2, 3, 4, 5],
        [3, 4, 5, 6]])

        It is possible to perform usual operation on matrices, such as the
        matrix multiplication:

        >>> A(i0, i1).get_matrix()*ones(4, 1)
        Matrix([
        [ 6],
        [10],
        [14],
        [18]])
        """
        if 0 < self.rank <= 2:
            rows = self.data.shape[0]
            columns = self.data.shape[1] if self.rank == 2 else 1
            if self.rank == 2:
                mat_list = [] * rows
                for i in range(rows):
                    mat_list.append([])
                    for j in range(columns):
                        mat_list[i].append(self[i, j])
            else:
                mat_list = [None] * rows
                for i in range(rows):
                    mat_list[i] = self[i]
            return Matrix(mat_list)
        else:
            raise NotImplementedError(
                "missing multidimensional reduction to matrix.")

    def canon_bp(self):
        """
        Canonicalize using the Butler-Portugal algorithm for canonicalization
        under monoterm symmetries.

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2 = tensor_indices('m0,m1,m2', Lorentz)
        >>> A = tensorhead('A', [Lorentz]*2, [[2]])
        >>> t = A(m0,-m1)*A(m1,-m0)
        >>> t.canon_bp()
        -A(L_0, L_1)*A(-L_0, -L_1)
        >>> t = A(m0,-m1)*A(m1,-m2)*A(m2,-m0)
        >>> t.canon_bp()
        0
        """
        if not self.has_index_order:
            # TODO: Should this raise an error?
            return self
        if self.is_canon_bp:
            return self
        if not self.components:
            return self
        t = self.sorted_components()
        g, dummies, msym, v = t.canon_args()
        can = canonicalize(g, dummies, msym, *v)
        if can == 0:
            return S.Zero
        tmul = t.perm2tensor(can, True)
        tmul = tmul.renumber_dummies()
        return tmul

    def renumber_dummies(self):
        """
        TODO

        only if index order.

        TODO: test whether this correctly works with spinor dummies.
        """
        sign = S.One

        dummy_indices = self.virtual_dummy_indices_list
        repl_dummies = {}
        create_dummy = TensExpr._get_dummy_generator([])
        for index in dummy_indices:
            if -index in repl_dummies:
                repl_dummies[index] = -repl_dummies[-index]
                if index.tensortype.metric_antisym and index.is_up:
                    #import ipdb; ipdb.set_trace()
                    sign = -sign
                continue
            subs_dummy = create_dummy(index)
            if not index.is_up and index.tensortype.metric_antisym is not None:
                subs_dummy = -subs_dummy
            repl_dummies[index] = subs_dummy
        # TODO: correct bugs and add substitution rules.
        tmul = sign*self.xreplace(repl_dummies)
        tmul._is_canon_bp = self._is_canon_bp
        return tmul

    @property
    def types(self):
        # TODO: this should be independent from the indices.
        typs = set([i.tensortype for i in self.indices_set])
        typs = list(typs)
        typs.sort(key=lambda x: x.name)
        return typs

    def fun_eval(self, *index_tuples):
        """
        Return a tensor with free indices substituted according to ``index_tuples``

        Parameters
        ==========

        index_types : list of tuples ``(old_index, new_index)``

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> i, j, k, l = tensor_indices('i,j,k,l', Lorentz)
        >>> A, B = tensorhead('A,B', [Lorentz]*2, [[1]*2])
        >>> t = A(i, k)*B(-k, -j) + A(i, -j)
        >>> t.fun_eval((i, k),(-j, l))
        A(k, l) + A(k, L_0)*B(l, -L_0)
        """
        return self.xreplace(dict(index_tuples))

    def substitute_indices(self, *index_tuples):
        """
        TODO: add deprecation warning.
        TODO: add description of what advantages it gives compared to .xreplace/.subs: it substitutes also indices of opposite sign.

        Return a tensor with free indices substituted according to ``index_tuples``

        Parameters
        ==========

        ``index_types`` list of tuples ``(old_index, new_index)``

        Note: this method will neither raise or lower the indices, it will just replace their symbol.

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> i, j, k, l = tensor_indices('i,j,k,l', Lorentz)
        >>> A, B = tensorhead('A,B', [Lorentz]*2, [[1]*2])
        >>> t = A(i, k)*B(-k, -j); t
        A(i, L_0)*B(-L_0, -j)
        >>> t.substitute_indices((i,j), (j, k))
        A(j, L_0)*B(-L_0, -k)
        """
        repl_dict = dict(index_tuples)
        for old, new in index_tuples:
            repl_dict[-old] = -new
        return self.xreplace(repl_dict)

    # TODO: add notes, this is an experimental attempt, not linked to the main
    # code.
    def __perform_contract_metric(self, g):
        xrepl_dict = {}
        xrepl_inverse_dict = {}
        gfactor = [S.One]

        def add_parsed_couple(ind1, ind2):
            xrepl_dict[ind1] = ind2
            xrepl_inverse_dict[ind2] = ind1
            # TODO: add check if indices are both free, to save the metric.
            # TODO: restore FreeTensorIndex class?
            if not isinstance(ind1, DummyTensorIndex) and not isinstance(ind2, DummyTensorIndex):
                if ind1 == ind2:
                    pass
                import ipdb; ipdb.set_trace()
                #return None
                return g(ind2, -ind1)
            if -ind1 in xrepl_dict:
                # TODO: should this be popped?
                #import ipdb; ipdb.set_trace()
                indo1 = xrepl_dict[-ind1]
                return g(ind2, indo1)
                #return g(indo1, ind2)

        def determine_antisym_metric_sign(i1, i2, regu):
            a1 = i1
            a2 = -i2
            if regu:
                if a1.is_up == a2.is_up:
                    return 1
                else:
                    return -1
            else:
                if a1.is_up == a2.is_up:
                    return -1
                else:
                    return 1

        def add_substitution(ind1, ind2):
            if -ind1 in xrepl_dict:
                # g(ind1, indo)*g(-ind1, ind2)
                indo = xrepl_dict.pop(-ind1)
                xrepl_inverse_dict.pop(indo)
                if indo == -ind2:
                    # g(ind1, -ind2)*g(-ind1, ind2) ==> D
                    return determine_antisym_metric_sign(ind1, -ind2, True), None
                # indo could be free, check if ind2 is free and avoid this step?
                metric_subs = add_parsed_couple(-ind2, indo)
            elif ind2 in xrepl_dict:
                # g(-ind2, indo)*g(-ind1, ind2)
                indo = xrepl_dict.pop(ind2)
                xrepl_inverse_dict.pop(indo)
                if ind1 == indo:
                    import ipdb; ipdb.set_trace()
                    # g(-ind2, ind1)*g(-ind1, ind2) ==> D
                    return determine_antisym_metric_sign(ind1, -ind2, False), None
                metric_subs = add_parsed_couple(ind1, indo)
            elif -ind2 in xrepl_inverse_dict:
                # g(indo, -ind2)*g(-ind1, ind2)
                indo = xrepl_inverse_dict.pop(-ind2)
                xrepl_dict.pop(indo)
                if indo == ind1:
                    # g(ind1, -ind2)*g(-ind1, ind2)
                    return determine_antisym_metric_sign(ind1, -ind2, True), None
                metric_subs = add_parsed_couple(indo, -ind1)
            elif ind1 in xrepl_inverse_dict:
                # g(indo, ind1)*g(-ind1, ind2)
                indo = xrepl_inverse_dict.pop(ind1)
                xrepl_dict.pop(indo)
                if indo == -ind2:
                    # g(-ind2, ind1)*g(-ind1, ind2)
                    return determine_antisym_metric_sign(ind1, -ind2, False), None
                metric_subs = add_parsed_couple(indo, ind2)
            else:
                metric_subs = add_parsed_couple(ind1, ind2)

            return 0, metric_subs

        components = self.components
        antisym = g.index_types[0].metric_antisym
        #if not any(x == g for x in components):
        #    return self
        # list of positions of the metric ``g``
        if isinstance(self, Tensor):
            if self.component == g:
                gpos = [self]
            else:
                return self
        else:
            gpos = [i for i in self.args if isinstance(i, Tensor) and i.component == g] # TODO: rename gpos to `metric_list`
        if not gpos:
            return self
        mfactor = S.One
        elim = set()
        for tensl in gpos:
            if tensl in elim:
                continue
            typ = g.index_types[0]
            if typ.dim is None:
                raise ValueError('dimension not assigned')
            dim = typ.dim
            inde1, inde2 = tensl.indices_list
            free1 = tensl.free_indices_list
            dum1 = tensl.dummy_indices_list
            if not dum1:
                continue

            if len(dum1) == 1:
                dum10 = dum1[0]
                assert len(free1) == 1  # TODO: remove assertion
                free10 = free1[0]
                if not antisym:
                    factor_dim, repl = add_substitution(-dum10, free10)
                    if repl:
                        xrepl_dict[tensl] = repl
                    else:
                        elim.add(tensl)
                    if factor_dim:
                        mfactor *= dim
                else:
                    factor_dim, repl = add_substitution(-dum10, free10)
                    if repl:
                        xrepl_dict[tensl] = repl
                    else:
                        elim.add(tensl)
                    if factor_dim:
                        mfactor *= factor_dim*dim
                    if dum10.is_up:
                        # TODO: check position in tensl
                        import ipdb; ipdb.set_trace()
                        if dum10 == tensl.indices_list[0]:
                            # dummy is at first position
                            mfactor *= -1
                    else:
                        if dum10 == tensl.indices_list[1]:
                            mfactor *= -1

            elif len(dum1) == 2:
                assert len(dum1) == 2  # TODO: ...
                assert len(free1) == 0 # TODO: remove asserts
                dum11, dum12 = dum1

                if not antisym:
                    # dp0, dp1, dc0, dc1 = dum1[0]
                    if dum11 == -dum12:
                        # g(i, -i)
                        print("Case g(i, -i)")
                        mfactor *= dim
                        elim.add(tensl)
                        print("mfactor:\t", mfactor)
                    else:
                        # g(i0, i1)*p(-i1)
                        factor_dim, repl = add_substitution(-dum12, dum11)
                        if repl:
                            xrepl_dict[tensl] = repl
                        else:
                            elim.add(tensl)
                        if factor_dim:
                            mfactor *= dim
                else:
                    if dum11 == -dum12:
                        # g(i, -i) or g(-i, i)
                        mfactor *= dim
                        if dum11.is_up:
                            mfactor *= -1
                        elim.add(tensl)
                    else:
                        # A(i)*B(-i)*g(i, -i) or something like that
                        if dum11.is_up:
                            factor_dim, repl = add_substitution(-dum12, dum11)
                            if repl:
                                xrepl_dict[tensl] = repl
                            else:
                                elim.add(tensl)
                            if factor_dim:
                                mfactor *= factor_dim*dim
                        else:
                            factor_dim, repl = add_substitution(-dum11, dum12)
                            if repl:
                                xrepl_dict[tensl] = repl
                            else:
                                elim.add(tensl)
                            if factor_dim:
                                mfactor *= factor_dim*dim

        if isinstance(self, Tensor):
            res = self if self not in elim else S.One
        else:
            res = self.func(*[arg for arg in self.args if arg not in elim])
        return mfactor*res.xreplace(xrepl_dict)*gfactor[0]

    def contract_metric(self, metric):
        """
        Raise or lower indices with the metric ``g``

        Parameters
        ==========

        g : metric

        Notes
        =====

        see the ``TensorIndexType`` docstring for the contraction conventions

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2 = tensor_indices('m0,m1,m2', Lorentz)
        >>> g = Lorentz.metric
        >>> p, q = tensorhead('p,q', [Lorentz], [[1]])
        >>> t = p(m0)*q(m1)*g(-m0, -m1)
        >>> t.canon_bp()
        metric(L_0, L_1)*p(-L_0)*q(-L_1)
        >>> t.contract_metric(g).canon_bp()
        p(L_0)*q(-L_0)
        """
        # TODO: metric should possess info about dim without accessing indices.
        new_t = self.__perform_contract_metric_xreplace(metric)
        if isinstance(new_t, TensExpr):
            # TODO: renumber_dummies should become aware of spinor indices.
            # TODO: also add tests for this.
            pass
            #new_t = new_t.renumber_dummies()
        return new_t

    def __perform_contract_metric_xreplace(self, g):
        # TODO: documentation
        dim = g.index_types[0].dim

        # TODO: move to Tensor class.
        if isinstance(self, Tensor):
            args = [self]
        else:
            args = self.args

        components = self.components
        antisym = g.index_types[0].metric_antisym
        # find the first metric inside the expression arguments which has at
        # least one non-free index, then perform the contraction:
        target_metric = None
        for arg in args:
            if not isinstance(arg, TensExpr):
                continue
            if arg.component != g:
                continue
            if len(arg.free_indices_list) == 2:
                continue
            target_metric = arg
            break
        if target_metric is None:
            return self

        free = target_metric.free_indices_list
        dum = target_metric.dummy_indices_list
        if not antisym:
            # metric is not antisymmetric:
            if len(free) == 1:
                dum1 = dum[0]
                free1 = free[0]
                new_expr = self.xreplace({target_metric: S.One}).xreplace({-dum1: free1})
                return new_expr.contract_metric(g)
            else:
                dum1, dum2 = dum
                if dum1 == -dum2:
                    if dim is None:
                        raise ValueError("index does not specify its dimension")
                    new_expr = dim*self.xreplace({target_metric: S.One})
                    if not isinstance(new_expr, TensExpr):
                        return new_expr
                    return new_expr.contract_metric(g)
                ndum = self._get_dummy_generator(self.dummy_indices_set)(dum1)
                new_expr = self.xreplace({target_metric: S.One, -dum1: ndum, -dum2: -ndum})
                return new_expr.contract_metric(g)
        else:
            # antisymmetric metric:
            if len(free) == 1:
                dum1 = dum[0]
                free1 = free[0]
                # check if the free index is the first one:
                free_is_first = bool(target_metric.indices_list[0] == free1)

                # sign cases:
                # call `f` the free index, `d` the dummy index.
                # g( f,  d) ==> +
                # g( f, -d) ==> -
                # g(-f,  d) ==> +
                # g(-f, -d) ==> -
                # g( d,  f) ==> -
                # g( d, -f) ==> -
                # g(-d,  f) ==> +
                # g(-d, -f) ==> +
                sign = (1 if free_is_first else -1)*(1 if dum1.is_up else -1)
                new_expr = sign*self.xreplace({target_metric: S.One}).xreplace({-dum1: free1})
                return new_expr.contract_metric(g)
            else:
                dum1, dum2 = dum
                if dum1 == -dum2:
                    # TODO: find a better way to suppress an element from the args (using strategies?)

                    # sign cases:
                    # g( d, -d) ==> -
                    # g(-d,  d) ==> +
                    sign = -1 if dum1.is_up else 1
                    if dim is None:
                        raise ValueError("index does not specify its dimension")
                    new_expr = sign*dim*self.xreplace({target_metric: S.One})
                    if not isinstance(new_expr, TensExpr):
                        return new_expr
                    return new_expr.contract_metric(g)
                # sign cases:  # TODO: describe replacement process.
                # in the resulting expression:
                # resulting `d1` ==> flips valence
                # resulting `d2` ==> keeps valence # TODO: check
                # g( d1,  d2) ==> -
                # g( d1, -d2) ==> +
                # g(-d1,  d2) ==> -
                # g(-d1, -d2) ==> +
                sign = 1 if not dum2.is_up else -1
                ndum = self._get_dummy_generator(self.dummy_indices_set)(dum1)
                new_expr = sign * self.xreplace({target_metric: S.One, -dum1: ndum, -dum2: -ndum})
                return new_expr.contract_metric(g)

    @staticmethod
    def _get_dummy_generator(dummy_indices):
        # TODO documentation.
        dummy_type_counter = defaultdict(lambda: 0)

        def _create_dummy(index):
#            if isinstance(index, DummyTensorIndex):
#                # TODO check collisions!
#                return index
            dummy_fmt = index.tensortype.dummy_fmt[:-3]
            # TODO: this is repeated twice, find a way to write it just once:
            count = dummy_type_counter[index.tensortype]
            dummy_type_counter[index.tensortype] = count + 1
            dummy_name = dummy_fmt + "_" + str(count)
            dummy = DummyTensorIndex(*((dummy_name,)+index.args[1:]))
            # if dummy has already been used by a previous expression (indices
            # that were already contracted when passed to the TensMul
            # constructor), then choose another dummy:
            if dummy in dummy_indices:
                dummy = _create_dummy(index)
            return dummy

        #_create_dummy.dummy_type_counter = dummy_type_counter
        return _create_dummy

    @property
    def rank(self):
        return len(self.free_indices_set)

    @property
    def ext_rank(self):
        return len(self.indices_set)

    def _eval_simplify(self, ratio, measure):
        # this is a way to simplify a tensor expression.

        # This part walks for all `TensorHead`s appearing in the tensor expr
        # and looks for `simplify_this_type`, to specifically act on a subexpr
        # containing one type of `TensorHead` instance only:
        expr = self
        for i in list(set(self.components)):
            if hasattr(i, 'simplify_this_type'):
                expr = i.simplify_this_type(expr)
        # TODO: missing feature, perform metric contraction.
        return expr

    @staticmethod
    def _assert_automatrix_indices_structure(args):
        # TODO
        pass

    @staticmethod
    def _replace_automatrix_indices_in_args(args):
        # TODO

        # Replacement cases:

        # A(auto_left, -auto_right)*B(auto_left, -auto_right) 
        #   ===> A(auto_left, -s)*B(s, -auto_right)

        # A(auto_left)*B(auto_left, -auto_right)
        #   ===> A(-s)*B(s, -auto_left)

        # A(auto_left, -auto_right)*B(auto_left)
        #   ===> A(auto_left, -s)*B(s)

        # A(auto_left)*B(auto_left)
        #   ===> A(-s)*B(s)

        dummy_counter = defaultdict(lambda: 0)
        def get_dummy(indextype):
            c = dummy_counter[indextype]
            dummy_counter[indextype] += 1
            return FreeTensorIndex('dummy_bind_{0}'.format(c), indextype, True)

        d = {}
        last_bindable_by_indextype = {}
        for arg_i, arg in enumerate(args):
            if not isinstance(arg, TensExpr):
                continue
#           bind_indices = [i for i in arg.indices_set if isinstance(i, BindableTensorIndex)]
#           if len(bind_indices)>0:
#               import ipdb; ipdb.set_trace()
            #auto_left_indices = [i for i in bind_indices if i.name == 'auto_left']
            #auto_right_indices = [i for i in bind_indices if i.name == 'auto_right'] 
            # TODO: create a method to return this `bind_types` object
            bind_types = defaultdict(lambda: [None, None])
            for i in arg.indices_set:
                if not isinstance(i, BindableTensorIndex):
                    continue
                if i.name == 'auto_left':
                    bind_types[i.tensortype][0] = i
                else:
                    bind_types[i.tensortype][1] = i

            for key, indices in bind_types.items():
                al, ar = indices
                # TODO: remove these assertions:
                if ar is not None:
                    assert ar.name == 'auto_right'
                if al is not None:
                    assert al.name == 'auto_left'

                # TODO: bindable indices should count as free, add test.
                if key in last_bindable_by_indextype:
                    # TODO: remove
#                    if len(args) == 2 and str(args[0]).startswith('C(auto_left') and str(args[1]).startswith('C(auto_left'):
#                        import ipdb; ipdb.set_trace()
                    # TODO: replace indices and update
                    # last_bindable_by_indextype
                    dummy = get_dummy(key)
                    if ar is not None:
                        to_bind_left, pos = last_bindable_by_indextype[key]
                        last_bindable_by_indextype[key] = ar
                    else:
                        to_bind_left, pos = last_bindable_by_indextype.pop(key)
                    args[pos] = args[pos].xreplace({to_bind_left: -dummy})
                    args[arg_i] = args[arg_i].xreplace({al: dummy})
                    if to_bind_left == al and ar is not None:
                        #args[arg_i] = args[arg_i].xreplace({ar: al})
                        last_bindable_by_indextype[key] = ar
                        
                else:
                    # TODO: update last_bindable_by_indextype
                    if ar is not None:
                        last_bindable_by_indextype[key] = (ar, arg_i)
                    else:
                        last_bindable_by_indextype[key] = (al, arg_i)
#                d[al] = arg_i
#                if ar is not None:
#                    d[ar] = arg_i

        return args

    @staticmethod
    def _detect_contracted_indices_in_args(args):
        args = TensMul._replace_automatrix_indices_in_args(args)
        # TODO: relocate to TensExpr? or not?

        # this assumes that indices in args are contracted multiplicatively.

        # TODO: finish.
#        free_indices = set([])

        free_indices_dict = {}
        dummy_indices = set([])

        # TODO: the logic to detect wrong index signatures should be made
        # independent so it can be used with the unify module.

        xreplace_list_dict = [{} for i in args]
        all_indices_list = []
        #kept_back_dummy_indices = set([])
        dummy_repeated_counter = defaultdict(lambda: 0)

        # TODO: should this be a search through all of the args tree or just
        # through the first children?
        for arg_i, arg in enumerate(args):
            if not isinstance(arg, TensExpr):
                # if it is not a TensExpr, it cannot contain indices.
                continue
            # TODO: every TensExpr should possess free indices.

            # if it is possible to determine the index order, get the list,
            # otherwise an unsorted set:
            if arg.has_index_order:
                # TODO: rename has_index_order to are_indices_sorted/has_index_order ?
                index_walking_list = arg.indices_list
                # TODO: add tests for .free_indices_set, .free_indices_list, .dummy_indices_set, .dummy_indices_list
            else:
                index_walking_list = arg.indices_set

            # iterate over the just determined list or set:
            for index in index_walking_list:
                #paired_position = None
                if not isinstance(index, TensorIndex):
                    # TODO: how should it behave in cases like this one?
                    continue
#                 if index in arg.dummy_indices_set and -index in arg.dummy_indices_set:
#                     1/0
#                     # `index` and `-index` are a contracted pair inside `arg`,
#                     # the code will try not to renumber them, i.e. they are put
#                     # into a special variable (`kept_back_dummy_indices`). If
#                     # two args have two identical pairs, the code will be forced
#                     # to renumber one of the pairs.
#                     if index not in kept_back_dummy_indices:
#                         kept_back_dummy_indices.add(index)
#                         continue
#                     # TODO: SHOULD kept_back_dummy_indices be included into
#                     # dummy_indices of the object being constructed?
                if isinstance(index, DummyTensorIndex):
                    if index in dummy_indices:
                        # if another dummy index with the same name already
                        # exists, relabeling is necessary:

                        #drcount = dummy_repeated_counter[index]
                        #index = FreeTensorIndex(*(("DtFidxctr{0}X{1}".format(index.args[0], drcount),)+index.args[1:]))
                        #index = FreeTensorIndex(*index.args)
                        pass
                    else:
                        dummy_indices.add(index)
                        continue
                if index in free_indices_dict:
                    if not isinstance(index, BindableTensorIndex):
                    # TODO: remove if
                        raise ValueError('wrong index signature')

#                    # This is the case in which there are auto-matrix indices.
#                    other_arg_i = free_indices_dict[index]
#                # TODO: there's a better way of getting `auto_left` index:
#                    auto_left = BindableTensorIndex('auto_left', *index.args[1:])
#                    auto_right = BindableTensorIndex('auto_right', *index.args[1:])
#                    prev_arg = args[other_arg_i]
#                    #TODO: remove assertion
#                    assert auto_left.is_up
#                    #import ipdb; ipdb.set_trace()
##                    free_indices_dict.pop(auto_left, None)
##                    free_indices_dict.pop(-auto_right, None)
#                    fri = FreeTensorIndex(*(('dummy_contract',)+index.args[1:]))
#
#                    paired_position = free_indices_dict.pop(index)
#                    args[arg_i] = args[arg_i].xreplace({index: -index})
#                    opo = all_indices_list.index([index, paired_position, None])
#                    all_indices_list[opo][2] = arg_i

#                    # Case: A(auto_left)*B(auto_left)
#                    #  ===> A(d)*A(-d)
#                    if -auto_right not in prev_arg and -auto_right not in arg:
#                        args[other_arg_i] = args[other_arg_i].xreplace({auto_left: fri})
#                        args[arg_i] = args[arg_i].xreplace({auto_left: -fri})
#                    # Case: A(auto_left)*B(auto_left, -auto_right)
#                    #  ===> A(d)*B(-d, auto_left)
#                    elif auto_left in prev_arg and -auto_right not in arg:
#                        args[other_arg_i] = args[other_arg_i].xreplace({auto_left: fri})
#                        args[arg_i] = args[arg_i].xreplace({auto_left: -fri})
#                    # Case: A(auto_left, -auto_right)*B(auto_left)
#                    #  ===> A(auto_left, -d)*B(d)
#                    elif -auto_right not in arg:
#                        args[other_arg_i] = args[other_arg_i].xreplace({-auto_right: fri})
#                        args[arg_i] = args[arg_i].xreplace({auto_left: -fri})
#                    # Case: A(auto_left, -auto_right)*B(auto_left, -auto_right)
#                    #  ===> A(auto_left, -d)*B(d, -auto_right)
#                    else:
#                        args[other_arg_i] = args[other_arg_i].xreplace({-auto_right: fri})
#                        args[arg_i] = args[arg_i].xreplace({auto_left: -fri})

#                    free_indices_dict[fri] = arg_i
#                    free_indices_dict[-fri] = other_arg_i

#                        if index.name == "auto_left":
#                            fri = FreeTensorIndex(*(('dummy_contract',)+index.args[1:]))
#                            args[arg_i] = args[arg_i].xreplace({index: fri})
#                            index = fri
#                        elif index.name == "auto_right":
#                            free_indices_dict[fri] = other_pos
#                            args[other_pos] = args[other_pos].xreplace({index: -fri})
#                        else:
#                            raise ValueError('')

                if -index in free_indices_dict:
                    # TODO: wrong comment:::
                    # contract the index:
                    # this means substitute it by dummies.
                    paired_position = free_indices_dict.pop(-index)
#                    temp_dummy_identified.add(index)
#                    temp_dummy_identified.add(-index)
#                    free_indices.remove(-index)
                    opo = all_indices_list.index([-index, paired_position, None])
                    all_indices_list[opo][2] = arg_i
                else:
                    free_indices_dict[index] = arg_i
                #    continue
                all_indices_list.append([index, arg_i, None])
                #index_to_component_dict[index] = arg_i

        _create_dummy = TensExpr._get_dummy_generator(dummy_indices)

        for index, pos1, pos2 in all_indices_list:
            if pos2 is None:
                continue
#            if index not in temp_dummy_identified:
#                continue
#            temp_dummy_identified.remove(index)
#            temp_dummy_identified.remove(-index)

            dummy = _create_dummy(index)

#            pos1 = index_to_component_dict[index]
#            pos2 = index_to_component_dict[-index]

            xreplace_list_dict[pos1].update({index: dummy})
            xreplace_list_dict[pos2].update({-index: -dummy})
            # free_indices.remove(-index)
            dummy_indices.add(-dummy)
            dummy_indices.add(dummy)

        for i in range(len(args)):
            if not isinstance(args[i], TensExpr):
                continue

#             excluded_dummies = args[i].atoms(DummyTensorIndex) #- args[i].dummy_indices
#             pre_xreplace = {}
#             while excluded_dummies:
#                 j = excluded_dummies.pop()
#                 if -j in excluded_dummies:
#                     dummy = _create_dummy(j)
#                     pre_xreplace.update({j: dummy, -j: -dummy})
#                     excluded_dummies.remove(-j)
# #            for j in excluded_dummies:
# #                1/0
# #                xreplace_list_dict[i].update({j: _create_dummy(j)})
#             #new_dummies = xreplace_list_dict[i].values()
#             args[i] = args[i].xreplace(pre_xreplace)

            args[i] = args[i].xreplace(xreplace_list_dict[i])

#TODO: remove
        # Replace `-auto_right` with `auto_left` when there is no more `auto_left`:
#        for arg_i, arg in enumerate(args):
#            if not isinstance(arg, TensExpr):
#                continue
#            for index in arg.indices_set:
#                if not isinstance(index, BindableTensorIndex):
#                    continue
#                if index.name != 'auto_right':
#                    continue
#                # TODO: there's a better way of getting `auto_left` index:
#                auto_left = BindableTensorIndex('auto_left', index.args[1], 1)
#                # TODO: remove assert
#                assert auto_left.is_up
#                if auto_left not in arg.indices_set:
#                    print("test\n")
#                    args[arg_i] = args[arg_i].xreplace({index: auto_left})
            

        # TODO: add substitutions of dummy indices in subexpressions.
        # i.e., the numbering has to be consistent.

        # TODO: All TensExpr need to possess two sets, free and dummy indices.
        # Subclasses of TensExpr are divided into two classes: those with an
        # ordered free index list, and those with just a set.

        free_indices = set([])
        for i in free_indices_dict.keys():
            if isinstance(i, DummyTensorIndex):
                dummy_indices.add(i)
            else:
                free_indices.add(i)
        return args, free_indices, dummy_indices

    @property
    def indices_set(self):
        # TODO: add tests
        return self.atoms(TensorIndex)

    @property
    def indices_list(self):
        self._assert_has_index_order()
        # TODO: add tests
        ind = []
        for arg in self.args:
            if not isinstance(arg, TensExpr):
                continue
            ind.extend(arg.indices_list)
        return ind

#        indices_list = []
#        for i in self.args:
#            if not isinstance(i, TensExpr):
#                continue
#            if not isinstance(i, CanonicalizableTensor):
#                raise ValueError("No order to determine.")
#            indices_list.extend(i.indices)
#        return indices_list

    @property
    def free_indices_set(self):
        return self._free_indices

    @property
    def free_indices_list(self):
        self._assert_has_index_order()
        # TODO: wrong, they need to be ordered.
        indices = []
        for i in self.args:
            if not isinstance(i, TensExpr):
                continue
            # TODO: how should this behave in case of non-canonicalizable tensors?
            for j in i.free_indices_list:
                indices.append(j)
            #indices.extend(i.free_indices_list())
        return indices

    @property
    def dummy_indices_set(self):
        # TODO: define how TensAdd makes indices inheritable, write docstring.
        dis = set([])
        for arg in self.args:
            if not isinstance(arg, TensExpr):
                continue
            dis.update(arg.dummy_indices_set)
        return dis

    @property
    def dummy_indices_list(self):
        self._assert_has_index_order()
        d_indices = []
        for arg in self.args:
            if not isinstance(arg, TensExpr):
                continue
            d_indices.extend(arg.dummy_indices_list)
        return d_indices

    @property
    def virtual_free_indices_list(self):
        self._assert_has_index_order()
        indices = self.indices_list
        dummies = self.dummy_indices_set
        return [i for i in indices if -i not in dummies]

    @property
    def virtual_free_indices_set(self):
        return set(self.virtual_free_indices_list)

    @property
    def virtual_dummy_indices_list(self):
        self._assert_has_index_order()
        dummies = self.dummy_indices_list
        return [i for i in dummies if -i in dummies]

    @property
    def virtual_dummy_indices_set(self):
        return set(self.virtual_dummy_indices_list)

    # TODO: deprecate this:
    # TODO: use instead ==> free_indices_set
    #@deprecated()
    @property
    def free_args(self):
        return self.free_indices_set

    @property
    def is_canon_bp(self):
        self._assert_has_index_order()
        return self._is_canon_bp

    # TODO: add rank attribute to Tensor (and test it).

    # TODO: rename `get_tensor_index_permutation` ?
    def get_tensor_permutation(self):
        """
        Returns the list representing the permutation of the tensor indices.
        """
        self._assert_has_index_order()
        indices = self.indices_list
        sorted_indices = self.get_canon_sorted_indices()
        g = [sorted_indices.index(i) for i in indices]
        # TODO: when should these be flipped?
        g.extend([len(g), len(g)+1])
        return g

    def canon_args(self):
        """
        Returns ``(g, dummies, msym, v)``, the entries of ``canonicalize``

        see ``canonicalize`` in ``tensor_can.py``
        """
        self._assert_has_index_order()
        # to be called after sorted_components
        from sympy.combinatorics.permutations import _af_new

        virtual_dummies = self.virtual_dummy_indices_list

        # get `g` list:
        g = self.get_tensor_permutation()

        # TODO: this is wrong!!!
        # get `dummies` object:
        dummies_by_type = defaultdict(lambda: [])
        for i, ind in enumerate(self.indices_list):
            #if not isinstance(ind, DummyTensorIndex):
            if ind not in virtual_dummies:
                continue
            dummies_by_type[ind.tensortype].append(g[i])
        dummies = []
        ind_types = sorted(dummies_by_type.keys(), key=lambda a: a._name)
        for typ in ind_types:
            dummies.append(sorted(dummies_by_type[typ]))

            #return 1 if key in virtual_dummies else 0

#        if dummies != []:
#            dummies = [dummies]

        # get the expression for `msym`:
        msym = []
        prev = None
        for index in virtual_dummies:
            if not index.is_up:
                continue
            typ = index.tensortype
            if typ == prev:
                continue
            msym.append(typ.metric_antisym)
            prev = typ

        # now calculate `v`:
        numtyp = []
        prev = None
        for t in self.components:
            # TODO: add to documentation: canonicalizable tensors both have
            # ordered indices and have Tensor components.
            if t == prev:
                numtyp[-1][1] += 1
            else:
                prev = t
                numtyp.append([prev, 1])
        v = []
        for h, n in numtyp:
            if h._comm == 0 or h._comm == 1:
                comm = h._comm
            else:
                comm = TensorManager.get_comm(h._comm, h._comm)
            v.append((h._symmetry.base, h._symmetry.generators, n, comm))

        return _af_new(g), dummies, msym, v

    def get_canon_sorted_indices(self):
        self._assert_has_index_order()
        virtual_dummies = self.virtual_dummy_indices_list
        # TODO: rename `indices_list` to `extended_indices_list` ?

        def free_first(key):
            dummy_or_free = 1 if key in virtual_dummies else 0
            notup = not key.is_up
            return dummy_or_free, key.tensortype, key.name, notup

        indices = sorted(self.virtual_free_indices_list) + sorted(virtual_dummies, key=free_first)
        return indices

    def _assert_has_index_order(self):
        if not self.has_index_order:
            raise ValueError("tensor has no index order")

    def perm2tensor(self, g, canon_bp=False):
        """
        Returns the tensor corresponding to the permutation ``g``

        ``g``  permutation corresponding to the tensor in the representation
        used in canonicalization

        ``canon_bp``   if True, then ``g`` is the permutation
        corresponding to the canonical form of the tensor
        """
        self._assert_has_index_order()
        sign = 1
        if g[-1] != len(g) - 1:
            sign = -sign

        # TODO: rename `indices_list` to `extended_indices_list` ?

        indices = self.get_canon_sorted_indices()
        repl_indices = [None]*len(indices)
        for i in range(len(repl_indices)):
            repl_indices[i] = indices[g[i]]
        # TODO: this should just be an order replacement operation.
        # such order replacement operation should be acted upon by calling the TensExpr with the new indices.
        res = sign*self.xreplace(dict(zip(self.indices_list, repl_indices)))
        res._is_canon_bp = canon_bp
        return res

    @property
    def free(self):
        raise NotImplementedError("This feature is no longer available here,"\
            "use `sympy.tensor.deprecated_tensor` if you want it.")

    @property
    def dum(self):
        raise NotImplementedError("This feature is no longer available here,"\
            "use `sympy.tensor.deprecated_tensor` if you want it.")


@doctest_depends_on(modules=('numpy',))
class TensAdd(TensExpr):
    """
    Sum of tensors

    Parameters
    ==========

    free_args : list of the free indices

    Attributes
    ==========

    ``args`` : tuple of addends
    ``rank`` : rank of the tensor
    ``free_args`` : list of the free indices in sorted order

    Notes
    =====

    Sum of more than one tensor are put automatically in canonical form.

    Examples
    ========

    >>> from sympy.tensor.tensor import TensorIndexType, tensorhead, tensor_indices
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> a, b = tensor_indices('a,b', Lorentz)
    >>> p, q = tensorhead('p,q', [Lorentz], [[1]])
    >>> t = p(a) + q(a); t
    p(a) + q(a)
    >>> t(b)
    p(b) + q(b)

    Examples with components data added to the tensor expression:

    >>> from sympy import eye
    >>> Lorentz.data = [1, -1, -1, -1]
    >>> a, b = tensor_indices('a, b', Lorentz)
    >>> p.data = [2, 3, -2, 7]
    >>> q.data = [2, 3, -2, 7]
    >>> t = p(a) + q(a); t
    p(a) + q(a)
    >>> t(b)
    p(b) + q(b)

    The following are: 2**2 - 3**2 - 2**2 - 7**2 ==> -58

    >>> (p(a)*p(-a)).data
    -58
    >>> p(a)**2
    -58
    """

    def __new__(cls, *args, **kw_args):
        args = [sympify(x) for x in args if x]
        args = TensAdd._tensAdd_flatten(args)

        if not args:
            return S.Zero

        # TODO: remove
        # replace auto-matrix indices so that they are the same in all addends
        #args = TensAdd._tensAdd_check_automatrix(args)

        # now check that all addends have the same indices:
        free_indices = TensAdd._tensAdd_check(args)

        # if TensAdd has only 1 TensMul element in its `args`:
        if len(args) == 1 and isinstance(args[0], TensMul):
            # TODO: Should this just return a TensMul?
            #obj = Basic.__new__(cls, *args, **kw_args)
            #obj._free_indices = free_indices
            return args[0]

        # canonicalize all TensMul
        # TODO: remove canonicalization upon construction:
        args = [canon_bp(x) for x in args if x]
        # TODO: decide what to do
        is_canon_bp = kw_args.pop('is_canon_bp', True)
        args = [x for x in args if x]

        # if there are no more args (i.e. have cancelled out),
        # just return zero:
        if not args:
            return S.Zero

        # collect canonicalized terms
        # TODO: remove
        # args.sort(key=lambda x: (x.components, x.free_indices_set, x.dummy_indices_set))
        def sort_key(x):
            if not isinstance(x, TensExpr):
                return (-1, [], [], [])
            if x.has_index_order:
                return (x.coeff == 1, x.components, x.free_indices_list, x.dummy_indices_list)
            return (-1, [], [], [])
        args = TensAdd._tensAdd_collect_terms(args)
        args.sort(key=sort_key)
        if not args:
            return S.Zero
        # it there is only a component tensor return it
        if len(args) == 1:
            # TODO: REmove
            assert not isinstance(args[0], TensAdd)
            return args[0]

#        def _sort_key(x):
#            if not isinstance(x, TensExpr):
#                return [], []
#            return x.components, x.free_indices_list
        #args.sort()
        # TODO remove
        assert len(args) > 1

        #args = Tuple(*args)
        obj = Basic.__new__(cls, *args, **kw_args)
#        obj._args = tuple(args)  # TODO: is this required?
        obj._free_indices = free_indices
        # TODO: this will have to be changed after introducing operators:
        obj._dummy_indices = set([])  # args[0].atoms(TensorIndex).difference(free_indices)
        obj._is_canon_bp = is_canon_bp
        return obj

    @staticmethod
    def _tensAdd_flatten(args):
        # expand all TensAdd in args.
        # TODO: also expand Add
        new_arg = []
        for arg in args:
            if isinstance(arg, TensAdd):
                new_arg.extend(arg.args)
            else:
                new_arg.append(arg)
        return new_arg

        # TODO: OLD: flatten TensAdd, coerce terms which are not tensors to tensors

    @staticmethod
    def _tensAdd_check(args):
        # TODO: should this become args[i].free_indices?

        # check that all addends have the same free indices
        indices0 = _free_indices(args[0])  # TODO: set([x[0] for x in args[0].free])

        # TODO: if there are wildcards?
        # this is not optimal:
        #if not all([i.has_index_order for i in self.args if isinstance(i, TensExpr)]):
            #return indices0

        list_indices = [_free_indices(x) for x in args[1:]]
        if not all(x == indices0 for x in list_indices):
            raise ValueError('all tensors must have the same indices')
        return indices0

    @staticmethod
    def _tensAdd_collect_terms(args):
        # collect TensMul terms differing at most by their coefficient
        mapping = defaultdict(lambda: S.Zero)
        residual = S.Zero
        for arg in args:
            if isinstance(arg, Tensor):
                mapping[arg] += 1
            elif isinstance(arg, TensMul):
                if arg.coeff == 1:
                    mapping[arg] += 1
                elif arg.coeff == 0:
                    # ignore null arg
                    import ipdb; ipdb.set_trace()
                    continue
                else:
                    mapping[TensMul(*arg.args[1:])] += arg.coeff
            elif not isinstance(arg, TensExpr):
                residual += arg
            # TODO: handle other possible TensExpr cases, such as partial
            # derivative.
        coll = [coeff*t for coeff, t in mapping.items() if coeff]
        if residual != 0:
            coll = [residual] + coll
        return [t for t in coll if t != 0]

#         a = []
#         prev = args[0]
#         prev_coeff = prev._coeff
#         changed = False
#
#         for x in args[1:]:
#             # if x and prev have the same tensor, update the coeff of prev
#             if x.components == prev.components \
#                     and x.free == prev.free and x.dum == prev.dum:
#                 prev_coeff = prev_coeff + x._coeff
#                 changed = True
#                 op = 0
#             else:
#                 # x and prev are different; if not changed, prev has not
#                 # been updated; store it
#                 if not changed:
#                     a.append(prev)
#                 else:
#                     # get a tensor from prev with coeff=prev_coeff and store it
#                     if prev_coeff:
#                         t = TensMul.from_data(prev_coeff, prev.components,
#                             prev.free, prev.dum)
#                         a.append(t)
#                 # move x to prev
#                 op = 1
#                 pprev, prev = prev, x
#                 pprev_coeff, prev_coeff = prev_coeff, x._coeff
#                 changed = False
#         # if the case op=0 prev was not stored; store it now
#         # in the case op=1 x was not stored; store it now (as prev)
#         if op == 0 and prev_coeff:
#             prev = TensMul.from_data(prev_coeff, prev.components, prev.free, prev.dum)
#             a.append(prev)
#         elif op == 1:
#             a.append(prev)
#         return a

#    @property
#    def dummy_indices(self):
#        return self.atoms(DummyTensorIndex)

    @property
    def free_indices_list(self):
        raise ValueError('no order in TensAdd')

    @property
    def dummy_indices_list(self):
        raise ValueError('no order in TensAdd')

    @property
    def rank(self):
        # TODO: remove
        # return self.args[0].rank
        return len(_free_indices(self.args[0]))

#    @property
#    def free_args(self):
#        return self.args[0].free_args

    def __call__(self, *indices):
        """Returns tensor with ordered free indices replaced by ``indices``

        Parameters
        ==========

        indices

        Examples
        ========

        >>> from sympy import Symbol
        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
        >>> D = Symbol('D')
        >>> Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
        >>> i0,i1,i2,i3,i4 = tensor_indices('i0:5', Lorentz)
        >>> p, q = tensorhead('p,q', [Lorentz], [[1]])
        >>> g = Lorentz.metric
        >>> t = p(i0)*p(i1) + g(i0,i1)*q(i2)*q(-i2)
        >>> t(i0,i2)
        metric(i0, i2)*q(L_0)*q(-L_0) + p(i0)*p(i2)
        >>> t(i0,i1) - t(i1,i0)
        0
        """
        free_args = sorted(self.free_indices_set)
        indices = list(indices)
        if [x._tensortype for x in indices] != [x._tensortype for x in free_args]:
            raise ValueError('incompatible types')
        if indices == free_args:
            return self
        index_tuples = list(zip(free_args, indices))
        a = [x.fun_eval(*index_tuples) for x in self.args]
        res = TensAdd(*a)

        return res

    def canon_bp(self):
        """
        canonicalize using the Butler-Portugal algorithm for canonicalization
        under monoterm symmetries.
        """
        if not self.has_index_order:
            return self
        args = [x.canon_bp() for x in self.args]
        res = TensAdd(*args, is_canon_bp=True)
        return res

    def equals(self, other):
        other = sympify(other)
        if isinstance(other, TensMul) and other.coeff == 0:
            return all(x.coeff == 0 for x in self.args)
        if isinstance(other, TensExpr):
            if self.rank != other.rank:
                return False
        if isinstance(other, TensAdd):
            if set(self.args) != set(other.args):
                return False
            else:
                return True
        t = self - other
        if not isinstance(t, TensExpr):
            return t == 0
        else:
            if isinstance(t, TensMul):
                return t._coeff == 0
            else:
                return all(x._coeff == 0 for x in t.args)

    def expand(self):
        return TensAdd(*[i.expand() for i in self.args])

    def __add__(self, other):
        return TensAdd(self, other)

    def __radd__(self, other):
        return TensAdd(other, self)

    def __sub__(self, other):
        return TensAdd(self, -other)

    def __rsub__(self, other):
        return TensAdd(other, -self)

    def __mul__(self, other):
        return TensAdd(*(x*other for x in self.args))

    def __rmul__(self, other):
        return self*other

    def __div__(self, other):
        other = sympify(other)
        if isinstance(other, TensExpr):
            raise ValueError('cannot divide by a tensor')
        # TODO: remove
        # return TensAdd(*(x/other for x in self.args))
        return TensAdd(*[x/other for x in self.args])

    def __rdiv__(self, other):
        raise ValueError('cannot divide by a tensor')

    def __getitem__(self, item):
        return self.data[item]

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def contract_delta(self, delta):
        args = [x.contract_delta(delta) for x in self.args]
        t = TensAdd(*args)
        return canon_bp(t)

    def contract_metric(self, g):
        """
        Raise or lower indices with the metric ``g``

        Parameters
        ==========

        g :  metric

        contract_all : if True, eliminate all ``g`` which are contracted

        Notes
        =====

        see the ``TensorIndexType`` docstring for the contraction conventions
        """

        args = [x.contract_metric(g) if isinstance(x, TensExpr) else x for x in self.args]
        t = TensAdd(*args)
        # TODO: is canon_bp necessary?
        return canon_bp(t)

    def _print(self):
        a = []
        args = self.args
        for x in args:
            a.append(str(x))
        #a.sort()
        s = ' + '.join(a)
        s = s.replace('+ -', '- ')
        return s

    @staticmethod
    def from_TIDS_list(coeff, tids_list):
        """
        Given a list of coefficients and a list of ``TIDS`` objects, construct
        a ``TensAdd`` instance, equivalent to the one that would result from
        creating single instances of ``TensMul`` and then adding them.

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead, TensAdd
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> i, j = tensor_indices('i,j', Lorentz)
        >>> A, B = tensorhead('A,B', [Lorentz]*2, [[1]*2])
        >>> eA = 3*A(i, j)
        >>> eB = 2*B(j, i)
        >>> t1 = eA._tids
        >>> t2 = eB._tids
        >>> c1 = eA.coeff
        >>> c2 = eB.coeff
        >>> TensAdd.from_TIDS_list([c1, c2], [t1, t2])
        3*A(i, j) + 2*B(i, j)

        If the coefficient parameter is a scalar, then it will be applied
        as a coefficient on all ``TIDS`` objects.

        >>> TensAdd.from_TIDS_list(4, [t1, t2])
        4*A(i, j) + 4*B(i, j)

        """
        if not isinstance(coeff, (list, tuple, Tuple)):
            coeff = [coeff] * len(tids_list)
        tensmul_list = [TensMul.from_TIDS(c, t) for c, t in zip(coeff, tids_list)]
        return TensAdd(*tensmul_list)

    @property
    def data(self):
        return _tensor_data_substitution_dict[self]

    @data.setter
    def data(self, data):
        # TODO: check data compatibility with properties of tensor.
        _tensor_data_substitution_dict[self] = data

    @data.deleter
    def data(self):
        if self in _tensor_data_substitution_dict:
            del _tensor_data_substitution_dict[self]

    def __iter__(self):
        if not self.data:
            raise ValueError("No iteration on abstract tensors")
        return self.data.flatten().__iter__()


@doctest_depends_on(modules=('numpy',))
class Tensor(TensExpr):
    """
    """

    is_commutative = False

    def __new__(cls, tensor_head, indices, **kw_args):
        try:
            tids = TIDS.from_components_and_indices((tensor_head,), indices)
        except Exception as e:
            tids = None
        indices = Tensor._add_automatrix_indices(tensor_head, indices, tensor_head._matrix_behavior==True)
        Tensor._check_indices_compability(tensor_head, indices)
        indices, free_ind, dummy_ind = Tensor._detect_contracted_indices(indices)
        # TODO: replace contracted indices with dummy ones.
        obj = Basic.__new__(cls, tensor_head, Tuple(*indices), **kw_args)
        obj._tids = tids
        # TODO: do we need these two sets?
        obj._free_indices = free_ind
        obj._dummy_indices = dummy_ind
        obj._is_canon_bp = kw_args.get('is_canon_bp', False)
        return obj

    @staticmethod
    def _add_automatrix_indices(tensor_head, indices, add_trailing=False):
        ind_types = tensor_head.index_types
        if len(indices) > len(ind_types):
            raise ValueError('wrong number of indices')
        if (not add_trailing) and (len(indices) != len(ind_types)):
            raise ValueError('wrong number of indices')

        indices = [i for i in indices]
        typd = defaultdict(lambda: 0)

        for i, ind in enumerate(indices):
            typ_counter = typd[ind_types[i]]
            typd[ind_types[i]] += 1

            if ind != True:
                continue

            if typ_counter == 0:
                indices[i] = ind_types[i].auto_left
            elif typ_counter == 1:
                indices[i] = -ind_types[i].auto_right
            else:
                raise ValueError('wrong number of indices')

        if not add_trailing:
            return indices

        #indices = list(indices)

        for i in range(len(indices), len(ind_types)):
            typ_counter = typd[ind_types[i]]
            typd[ind_types[i]] += 1
            ind_typ = ind_types[i]

            if typ_counter == 0:
                indices.append(ind_typ.auto_left)
            elif typ_counter == 1:
                indices.append(-ind_typ.auto_right)
            else:
                raise ValueError('wrong number of indices')
        return indices

#        typd = collections.defaultdict(lambda: 0)
#        new_indices = [None]*len(indices)
#        index_pointer = 0
#        for i, ind_typ in enumerate(ind_types):
#            if indices[index_pointer].tensortype == ind_typ:
#                new_indices[i] = indices[index_pointer]
#                index_pointer += 1
#                continue
#            typ_counter = typd[ind_typ]
#            if typ_counter == 0:
#                new_indices[i] = ind_typ.auto_left
#            elif typ_counter == 1:
#                new_indices[i] = ind_typ.auto_right
#            else:
#                raise ValueError('wrong index signature')

    @staticmethod
    def _check_indices_compability(tensor_head, indices):
        ind_types = tensor_head.index_types
        # TODO: already asserted?
        if len(indices) != len(ind_types):
            raise ValueError('wrong number of indices')

        for index, ind_type in zip(indices, ind_types):
            if index.tensortype != ind_type:
                raise ValueError('incompatible index types')

# TODO: rename?
    @staticmethod
    def _detect_contracted_indices(indices):
        # TODO: finish.
        indices = list(indices)

        free_indices = set([])
        dummy_indices = set([])
        dum_type_counter = defaultdict(lambda: 0)
        for i, index in enumerate(indices):
            if not isinstance(index, TensorIndex):
                # TODO: remove to support wilds.
                raise ValueError('index is not a TensorIndex')
                continue
            if isinstance(index, DummyTensorIndex):
                # this means that a tensor with dummy indices is being
                # constructed, it has to be handled as a free index until the
                # contraction is detected:
                dummy_indices.add(index)
                continue
            if index in free_indices:
                raise ValueError('cannot repeat the same index')
            if index in dummy_indices:
                raise ValueError('wrong index signature')
            if -index in free_indices:
                free_indices.remove(-index)
                count = dum_type_counter[index.tensortype]
                dum_type_counter[index.tensortype] = count + 1
                dummy_fmt = index.tensortype.dummy_fmt[:-3]
                dummy_name = dummy_fmt + "_" + str(count)
                dummy_ind = DummyTensorIndex(*((dummy_name,)+index.args[1:]))
                dummy_indices.add(dummy_ind)
                dummy_indices.add(-dummy_ind)
                indices[i] = dummy_ind
                indices[indices.index(-index)] = -dummy_ind
            else:
                free_indices.add(index)
        return indices, set(free_indices), set(dummy_indices)

    @property
    def has_index_order(self):
        # TODO: add tests for `has_index_order`
        # TODO: check if this should always be true:
        return True

    @property
    def coeff(self):
        return S.One

    @property
    def component(self):
        return self.args[0]

    @property
    def components(self):
        return [self.args[0]]

    def split(self):
        return [self]

    def expand(self):
        return self

    def sorted_components(self):
        return self

    def get_indices(self):
        # TODO: Remove this.
        return list(self.args[1])

    @property
    def indices_list(self):
        return list(self.args[1])

    @property
    def free_indices_list(self):
        return [i for i in self.args[1] if i in self.free_indices_set]

#    @property
#    def free_indices_set(self):
#        return [i for i in self.args[1] if i in self.free_indices_set]

    @property
    def dummy_indices_list(self):
        return [i for i in self.args[1] if i not in self.free_indices_set]

    @property
    def dummy_indices_set(self):
        return set(self.dummy_indices_list)

    def as_base_exp(self):
        return self, S.One

    def __call__(self, *indices):
        # TODO: add docstring
        # TODO: is this correct?
        fi = self.free_indices_list[:len(indices)]
        fi.sort()
        if len(indices) > len(fi):
            raise ValueError("too many indices") # TODO create a unique function for this.
        return self.xreplace(dict(zip(fi, indices)))  # TODO: xreplace acts on all subexpressions, it could act on independent subexpressions in the args-tree.

    # TODO: put this into TensExpr?
    def __iter__(self):
        return self.data.flatten().__iter__()

    # TODO: put this into TensExpr?
    def __getitem__(self, item):
        return self.data[item]

    @property
    def data(self):
        return _tensor_data_substitution_dict[self]

    @data.setter
    def data(self, data):
        # TODO: check data compatibility with properties of tensor.
        _tensor_data_substitution_dict[self] = data

    @data.deleter
    def data(self):
        if self in _tensor_data_substitution_dict:
            del _tensor_data_substitution_dict[self]
        if self.metric in _tensor_data_substitution_dict:
            del _tensor_data_substitution_dict[self.metric]

    def __mul__(self, other):
        return TensMul(self, other)

    def __rmul__(self, other):
        return TensMul(other, self)

    def __div__(self, other):
        if isinstance(other, TensExpr):
            raise ValueError('cannot divide by a tensor')
        return TensMul(self, S.One/other, is_canon_bp=self.is_canon_bp)

    def __rdiv__(self, other):
        raise ValueError('cannot divide by a tensor')

    def __add__(self, other):
        return TensAdd(self, other)

    def __radd__(self, other):
        return TensAdd(other, self)

    def __sub__(self, other):
        return TensAdd(self, -other)

    def __rsub__(self, other):
        return TensAdd(other, self)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __neg__(self):
        return TensMul(S.NegativeOne, self)

    def _print(self):
        indices_str = str(tuple(self.get_indices())).replace(',)', ')')
        if indices_str == "()":
            return str(self.component.args[0])
        return str(self.component.args[0]) + indices_str

    def equals(self, other):
        if other == 0:
            return self.coeff == 0
        other = sympify(other)
        if not isinstance(other, TensExpr):
            assert not self.components
            return S.One == other

        def _get_compar_comp(self):
            t = self.canon_bp()
#             r = (t.coeff, tuple(t.components), \
#                     tuple(sorted(t.free)), tuple(sorted(t.dum)))
            r = (t.coeff, tuple(t.components), \
                    tuple(sorted(t.free_indices_set)), tuple(t.dummy_indices_set))
            return r

        return _get_compar_comp(self) == _get_compar_comp(other)

    def contract_delta(self, metric):
        return self.contract_metric(metric)

#     def canon_bp(self):
#         """
#         Canonicalize using the Butler-Portugal algorithm for canonicalization
#         under monoterm symmetries.
#
#         Examples
#         ========
#
#         >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
#         >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
#         >>> m0, m1, m2 = tensor_indices('m0,m1,m2', Lorentz)
#         >>> A = tensorhead('A', [Lorentz]*2, [[2]])
#         >>> t = A(m0,-m1)*A(m1,-m0)
#         >>> t.canon_bp()
#         -A(L_0, L_1)*A(-L_0, -L_1)
#         >>> t = A(m0,-m1)*A(m1,-m2)*A(m2,-m0)
#         >>> t.canon_bp()
#         0
#         """
#         if not self.has_index_order:
#             return self
#         if self.is_canon_bp:
#             return self
#         # t = self.sorted_components()
#         g, dummies, msym, v = self._tids.canon_args()
#         can = canonicalize(g, dummies, msym, *v)
#         if can == 0:
#             return S.Zero
#         # TODO remove
#         # tmul = self.func(self.component, self._tids.perm2tensor(can, True).to_indices())
#
#         tmul = self.perm2tensor(can, True)
#         return tmul

#        # TODO: remove
#        tmul = TensMul(*self._tids.perm2tensor(can, True).to_tensmul_args())
#        tmul._matrix_behavior_kinds = self._matrix_behavior_kinds
#        return tmul


@doctest_depends_on(modules=('numpy',))
class TensMul(TensExpr):
    """
    Product of tensors

    Parameters
    ==========

    coeff : SymPy coefficient of the tensor
    args

    Attributes
    ==========

    ``components`` : list of ``TensorHead`` of the component tensors
    ``types`` : list of nonrepeated ``TensorIndexType``
    ``free`` : list of ``(ind, ipos, icomp)``, see Notes
    ``dum`` : list of ``(ipos1, ipos2, icomp1, icomp2)``, see Notes
    ``ext_rank`` : rank of the tensor counting the dummy indices
    ``rank`` : rank of the tensor
    ``coeff`` : SymPy coefficient of the tensor
    ``free_args`` : list of the free indices in sorted order
    ``is_canon_bp`` : ``True`` if the tensor in in canonical form

    Notes
    =====

    ``args[0]``   list of ``TensorHead`` of the component tensors.

    ``args[1]``   list of ``(ind, ipos, icomp)``
    where ``ind`` is a free index, ``ipos`` is the slot position
    of ``ind`` in the ``icomp``-th component tensor.

    ``args[2]`` list of tuples representing dummy indices.
    ``(ipos1, ipos2, icomp1, icomp2)`` indicates that the contravariant
    dummy index is the ``ipos1``-th slot position in the ``icomp1``-th
    component tensor; the corresponding covariant index is
    in the ``ipos2`` slot position in the ``icomp2``-th component tensor.

    """

    def __new__(cls, *args, **kw_args):
        #components = []
        #indices = []
        # TODO: tensmul arguments are not flattened, i.e. there is still
        # a tensmul object inside the arguments when constructed. Add
        # test for this.

        def _parse_args(args):
            # TODO: this function should just collect the commutative scalars and put them into coeff
            # TODO: it was changed to include also non-commutative arguments.
            new_args = []
            coeff = S.One
            for arg in args:
                if isinstance(arg, TensMul):
                    if arg.args:
                        if not isinstance(arg.args[0], TensExpr):
                            coeff *= arg.args[0]
                            new_args.extend(arg.args[1:])
                        else:
                            new_args.extend(arg.args)
                    continue
                elif not isinstance(arg, TensExpr):
                    arg = sympify(arg)
                    coeff *= arg
#                    if arg.is_commutative:
#                        coeff *= arg
                    continue
                new_args.append(arg)
            return coeff, new_args

        coeff, new_args = _parse_args(args)
        if new_args == []:
            return coeff

        if coeff == S.Zero:
            return S.Zero

        if coeff != 1:
            # if coeff is not one, add it as first argument:
            new_args = [coeff] + new_args
        elif len(new_args) == 1:
            # if there is only one element and coefficient is one, don't create a TensMul:
            return new_args[0]

#        tids = None
#        if all([i.has_index_order for i in args if isinstance(i, TensExpr)]):
#            has_index_order = True
#        else:
#            has_index_order = False

            #try:
                # TODO: should tids be generated by multiplication?
                # TODO: make sure that tids only exists if TensMul is canonicalizable.
                #tids = TIDS.from_components_and_indices(components, indices)
#                seq = [i._tids for i in args if isinstance(i, TensExpr)]
#                if seq:
#                    tids = reduce(lambda x, y: x*y, seq)
                # TODO: this alternative does not work (makes canon_bp() get wrong results):
#                all_tids = [i._tids for i in new_args if isinstance(i, TensExpr)]
#                tids = all_tids[0]
#                for i in all_tids[1:]:
#                    tids *= i

#                if globals().get('flag_no', True):
#                    globals()['flag_no'] = False
#                    assert coeff*TensMul(*tids.to_tensmul_args()) == TensMul(*args)
#                    globals()['flag_no'] = True

#                has_index_order = True
#            except (ValueError, TypeError, AttributeError):
#                # TODO: how to handle index contraction with no TIDS?
#                # this applies to the class Tensor as well.
#                tids = None
#            except Exception:
#                tids = None

        # TODO: handle index contractions...
        new_args, free_indices, dummy_indices = TensExpr._detect_contracted_indices_in_args(new_args)

        obj = Basic.__new__(cls, *new_args)

        obj._types = []
#        for t in tids.components:
#            obj._types.extend(t._types)
#        obj._tids = tids
# TODO: this will fail:
#        obj._ext_rank = len(obj._tids.free) + 2*len(obj._tids.dum)
        obj._coeff = coeff
        obj._is_canon_bp = kw_args.get('is_canon_bp', False)
        obj._free_indices = free_indices
        obj._dummy_indices = dummy_indices
        return obj

    @property
    # TODO
    #@cacheit
    def has_index_order(self):
        return all([arg.has_index_order for arg in self.args if isinstance(arg, TensExpr)])

# TODO: rename get_indices() to @property indices and deprecate the old one?

## TODO: rename?
#    @staticmethod
#    def _detect_contracted_indices_in_args(args):
#        # TODO: finish.
#        #free = set([])
#        free_to_component_mapping = {}
#        dummy_to_arg_mapping = {}
#        dummy_type_counter = defaultdict(lambda: 0)
#        #dummies = set([])
#
#        for arg_i, arg in enumerate(args):
#            if not isinstance(arg, Tensor):
#                continue
#            for index in arg.get_indices():
#                if not isinstance(index, FreeTensorIndex):
#                    continue
#
#                if index in free_to_component_mapping:
#                    raise ValueError('wrong index signature')
#                if -index in free_to_component_mapping:
#                    # TODO: wrong comment:::
#                    # contract the index:
#                    # this means substitute it by dummies.
#                    dummy_to_arg_mapping[index] = arg_i
#                    dummy_to_arg_mapping[-index] = free_to_component_mapping.pop(-index)
#                free_to_component_mapping[index] = arg_i
#
#        for index, pos1 in dummy_to_arg_mapping.items():
#            dummy_fmt = index.tensortype.dummy_fmt[:-3]
#            # TODO: this is repeated twice, find a way to write it just once:
#            count = dummy_type_counter[index.tensortype]
#            dummy_type_counter[index.tensortype] = count + 1
#            dummy_name = dummy_fmt + "_" + str(count)
#            dummy = DummyTensorIndex(*((dummy_name,)+index.args[1:]))
#            args[pos1] = args[pos1].subs(index, dummy)
#            pos2 = dummy_to_arg_mapping.pop(-index)
#            args[pos2] = args[pos2].subs(-index, -dummy)
#
#        return args

# TODO: restore
#    @staticmethod
#    def from_data(coeff, components, free, dum, **kw_args):
#        tids = TIDS(components, free, dum)
#        return TensMul.from_TIDS(coeff, tids, **kw_args)

    # TODO: deprecate
    # TODO: or better remove?
    @staticmethod
    def from_TIDS(coeff, tids, **kw_args):
        # TODO: correct
        return coeff*TensMul(*tids.to_tensmul_args())
    # TODO: remove
#        components, free, dum = tids.to_tensmul_args()
#        free = map(lambda x: (FreeTensorIndex(*x[0].args), x[1], x[2]), free)
#        return TensMul(coeff, components, free, dum, **kw_args)

    @property
    def components(self):
        ret = []
        for i in self.args:
            if isinstance(i, Tensor):
                ret.append(i.component)
        return ret

    @property
    def coeff(self):
        return self._coeff

    @property
    def rank(self):
        return len(_free_indices(self))

    def equals(self, other):
        if other == 0:
            return self.coeff == 0
        other = sympify(other)
        if not isinstance(other, TensExpr):
            assert not self.components
            return self._coeff == other

#        def _get_compar_comp(self):
#            t = self.canon_bp()
#            r = (t.coeff, tuple(t.components), \
#                    tuple(sorted(t.free_indices_list)), \
#                    tuple(sorted(t.dummy_indices_list)))
#            return r

        s1 = self.expand()
        s2 = other.expand()

        if s1.has_index_order != s2.has_index_order:
            return False

        if s1.has_index_order:
            s1 = s1.canon_bp()
            s2 = s2.canon_bp()

        return s1 == s2

    def get_indices(self):
        """
        Returns the list of indices of the tensor

        The indices are listed in the order in which they appear in the
        component tensors.
        The dummy indices are given a name which does not collide with
        the names of the free indices.

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2 = tensor_indices('m0,m1,m2', Lorentz)
        >>> g = Lorentz.metric
        >>> p, q = tensorhead('p,q', [Lorentz], [[1]])
        >>> t = p(m1)*g(m0,m2)
        >>> t.get_indices()
        [m1, m0, m2]
        """
        # TODO: this does not return wilds
        indices = []
        # TODO: add tests for get_indices() for Tensor, TensMul.
        # get_indices() should work only when indices are sortable.
        for arg in self.args:
            if not isinstance(arg, TensExpr):
                continue
            indices.extend(arg.get_indices())

        return indices

#        return self.atoms(TensorIndex)

#        indices = [None]*self._ext_rank
#        start = 0
#        pos = 0
#        vpos = []
#        components = self.components
#        for t in components:
#            vpos.append(pos)
#            pos += t._rank
#        cdt = defaultdict(int)
#        # if the free indices have names with dummy_fmt, start with an
#        # index higher than those for the dummy indices
#        # to avoid name collisions
#        for indx, ipos, cpos in self.free:
#            if indx._name.split('_')[0] == indx._tensortype._dummy_fmt[:-3]:
#                cdt[indx._tensortype] = max(cdt[indx._tensortype], int(indx._name.split('_')[1]) + 1)
#            start = vpos[cpos]
#            indices[start + ipos] = indx
#        for ipos1, ipos2, cpos1, cpos2 in self.dum:
#            start1 = vpos[cpos1]
#            start2 = vpos[cpos2]
#            typ1 = components[cpos1].index_types[ipos1]
#            assert typ1 == components[cpos2].index_types[ipos2]
#            fmt = typ1._dummy_fmt
#            nd = cdt[typ1]
#            indices[start1 + ipos1] = TensorIndex(fmt % nd, typ1)
#            indices[start2 + ipos2] = TensorIndex(fmt % nd, typ1, False)
#            cdt[typ1] += 1
#        return indices

    def sorted_components(self):
        """
        TODO replace 'sorted_components'
        should this be a part of TensMul?

        TODO: test this with non-commutative scalars.

        Returns a tensor with sorted components.

        The sorting is done taking into account the commutation group
        of the component tensors.
        """
        #from sympy.combinatorics.permutations import _af_invert
        #cv = list(zip(self.args, range(len(self.components))))
        if self.coeff != 1:
            cv = list(zip(self.args[1:], range(len(self.components))))
        else:
            cv = list(zip(self.args, range(len(self.components))))
        sign = 1
        n = len(cv) - 1

        for i in range(n):
            for j in range(n, i, -1):
                comp0 = cv[j][0].component
                compm1 = cv[j-1][0].component
                c = compm1.commutes_with(comp0)
                if c not in [0, 1]:
                    continue
                if (compm1._types, compm1._name) > \
                        (comp0._types, comp0._name):
                    cv[j-1], cv[j] = cv[j], cv[j-1]
                    if c:
                        sign = -sign

        #cv.sort(key=lambda x: x[1])
        new_args = [x[0] for x in cv]
        return TensMul(sign*self.coeff, *new_args)

    def split(self):
        """
        Returns a list of tensors, whose product is ``self``

        Dummy indices contracted among different tensor components
        become free indices with the same name as the one used to
        represent the dummy indices.

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> a, b, c, d = tensor_indices('a,b,c,d', Lorentz)
        >>> A, B = tensorhead('A,B', [Lorentz]*2, [[1]*2])
        >>> t = A(a,b)*B(-b,c)
        >>> t
        A(a, L_0)*B(-L_0, c)
        >>> t.split()
        [A(a, L_0), B(-L_0, c)]
        """
        coeff = S.One
        args = []
        for i, arg in enumerate(self.args):
            if not isinstance(arg, TensExpr):
                coeff *= arg
                continue
            if not arg.has_index_order:
                # TODO: should this be skipped?
                continue
            dummies = arg.dummy_indices_list
            # TODO: Remove
            # frees = [FreeTensorIndex(*dummy.args) for dummy in dummies]
            frees = [FreeTensorIndex(*dummy.args) for dummy in dummies]
            args.append(coeff*arg.subs(dict(zip(dummies, frees))))
            coeff = S.One
        return args

        # TODO: remove
#        indices = self.get_indices()
#        pos = 0
#        components = self.components
#        if not components:
#            return [TensMul.from_data(self._coeff, [], [], [])]
#        res = []
#        for t in components:
#            t1 = t(*indices[pos:pos + t._rank])
#            pos += t._rank
#            res.append(t1)
#        res[0] = TensMul.from_data(self._coeff, res[0].components, res[0]._tids.free, res[0]._tids.dum, is_canon_bp=res[0].is_canon_bp)
#        return res

    # TODO
    def expand(self):
        """
        """
        # detect whether it has TensAdd subexpressions.
        #
        exp_pre = []
        for arg in self.args:
            if isinstance(arg, TensAdd):
                eargs = list(arg.args[:])
                exp_pre.append(eargs)
            else:
                exp_pre.append([arg])
        return TensAdd(*[TensMul(*item) for item in itertools.product(*exp_pre)])

    def __add__(self, other):
        return TensAdd(self, other)

    def __radd__(self, other):
        return TensAdd(other, self)

    def __sub__(self, other):
        return TensAdd(self, -other)

    def __rsub__(self, other):
        return TensAdd(other, -self)

    def __mul__(self, other):
        """
        Multiply two tensors using Einstein summation convention.

        If the two tensors have an index in common, one contravariant
        and the other covariant, in their product the indices are summed

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2 = tensor_indices('m0,m1,m2', Lorentz)
        >>> g = Lorentz.metric
        >>> p, q = tensorhead('p,q', [Lorentz], [[1]])
        >>> t1 = p(m0)
        >>> t2 = q(-m0)
        >>> t1*t2
        p(L_0)*q(-L_0)
        """
        other = sympify(other)
        if not isinstance(other, TensExpr):
            return TensMul(*(self.args + (other,)))
        if isinstance(other, TensAdd):
            return TensAdd(*[TensMul(self, x) for x in other.args])

        tmul = TensMul(self, other)
        return tmul.renumber_dummies()

    def __rmul__(self, other):
        other = sympify(other)
        return TensMul(other, *self.args)

    def __div__(self, other):
        other = sympify(other)
        if isinstance(other, TensExpr):
            raise ValueError('cannot divide by a tensor')
        #coeff = self._coeff/other
        #tmul = TensMul.from_TIDS(coeff, self._tids, is_canon_bp=self.is_canon_bp)
        tmul = TensMul(S.One/other, *self.args, is_canon_bp=self.is_canon_bp)
        return tmul

    def __rdiv__(self, other):
        raise ValueError('cannot divide by a tensor')

    def __getitem__(self, item):
        return self.data[item]

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def contract_delta(self, delta):
        t = self.contract_metric(delta)
        return t

    def __call__(self, *indices):
        """Returns tensor with ordered free indices replaced by ``indices``

        Examples
        ========

        >>> from sympy import Symbol
        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
        >>> D = Symbol('D')
        >>> Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
        >>> i0,i1,i2,i3,i4 = tensor_indices('i0:5', Lorentz)
        >>> g = Lorentz.metric
        >>> p, q = tensorhead('p,q', [Lorentz], [[1]])
        >>> t = p(i0)*q(i1)*q(-i1)
        >>> t(i1)
        p(i1)*q(L_0)*q(-L_0)
        """
        indices = list(indices)
        # get a set (not a list), so this applies even if the order is not
        # specified:
        sorted_free_indices = ([i for i in self.free_indices_set])
        sorted_free_indices.sort()
        if [x._tensortype for x in indices] != [x._tensortype for x in sorted_free_indices]:
            raise ValueError('incompatible types')
        if indices == sorted_free_indices:
            return self
        t = self.fun_eval(*list(zip(sorted_free_indices, indices)))

        # object is rebuilt in order to make sure that all contracted indices
        # get recognized as dummies, but only if there are contracted indices.
        if len(set(i if i.is_up else -i for i in indices)) != len(indices):
            return t.func(*t.args)
        return t

    def _print(self):
        def _format(o):
            if isinstance(o, TensAdd):
                return "({0})".format(o)
            if isinstance(o, Add):
                return "({0})".format(o)
            return str(o)
        outstr = '*'.join([_format(arg) for arg in self.args])
        if outstr.startswith('-1*'):
            outstr = '-'+outstr[3:]
        return outstr

# TODO: remove
#        if len(self.components) == 0:
#            return str(self._coeff)
#        indices = [str(ind) for ind in self.get_indices()]
#        pos = 0
#        a = []
#        for t in self.components:
#            if t._rank > 0:
#                a.append('%s(%s)' % (t.name, ', '.join(indices[pos:pos + t._rank])))
#            else:
#                a.append('%s' % t.name)
#            pos += t._rank
#        res = '*'. join(a)
#        if self._coeff == S.One:
#            return res
#        elif self._coeff == -S.One:
#            return '-%s' % res
#        if self._coeff.is_Atom:
#            return '%s*%s' % (self._coeff, res)
#        else:
#            return '(%s)*%s' %(self._coeff, res)

    @property
    def data(self):
        dat = _tensor_data_substitution_dict[self]
        if dat is None:
            return None
        return self.coeff * dat

    @data.setter
    def data(self, data):
        # TODO: check data compatibility with properties of tensor.
        _tensor_data_substitution_dict[self] = data

    @data.deleter
    def data(self):
        if self in _tensor_data_substitution_dict:
            del _tensor_data_substitution_dict[self]

    def __iter__(self):
        if self.data is None:
            raise ValueError("No iteration on abstract tensors")
        return (self.data.flatten()).__iter__()


def canon_bp(p):
    """
    Butler-Portugal canonicalization
    """
    if isinstance(p, TensExpr):
        return p.canon_bp()
    return p

def tensor_mul(*a):
    """
    product of tensors
    """
    if not a:
        return S.One
    t = a[0]
    for tx in a[1:]:
        t = t*tx
    return t


def riemann_cyclic_replace(t_r):
    """
    replace Riemann tensor with an equivalent expression

    ``R(m,n,p,q) -> 2/3*R(m,n,p,q) - 1/3*R(m,q,n,p) + 1/3*R(m,p,n,q)``

    """
    if not isinstance(t_r, TensExpr):
        # TODO: should it ever come to this place?
        return t_r
    # TODO: remove
    #free = sorted(t_r.free, key=lambda x: x[1])
    #m, n, p, q = [x[0] for x in free]
    m, n, p, q = t_r.free_indices_list
    t0 = S(2)/3*t_r
    t1 = - S(1)/3*t_r.substitute_indices((m,m),(n,q),(p,n),(q,p))
    t2 = S(1)/3*t_r.substitute_indices((m,m),(n,p),(p,n),(q,q))
    t3 = t0 + t1 + t2
    return t3

def riemann_cyclic(t2):
    """
    replace each Riemann tensor with an equivalent expression
    satisfying the cyclic identity.

    This trick is discussed in the reference guide to Cadabra.

    Examples
    ========

    >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead, riemann_cyclic
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> i, j, k, l = tensor_indices('i,j,k,l', Lorentz)
    >>> R = tensorhead('R', [Lorentz]*4, [[2, 2]])
    >>> t = R(i,j,k,l)*(R(-i,-j,-k,-l) - 2*R(-i,-k,-j,-l))
    >>> riemann_cyclic(t)
    0
    """
    t2 = t2.expand()
    if isinstance(t2, (TensMul, Tensor)):
        args = [t2]
    elif isinstance(t2, TensAdd):
        args = t2.args
    else:
        raise ValueError('wrong type')
    a1 = [x.split() for x in args]
    a2 = [[riemann_cyclic_replace(tx) for tx in y] for y in a1]
    a3 = [tensor_mul(*v) for v in a2]
    t3 = TensAdd(*a3)
    if not t3:
        return t3
    else:
        return canon_bp(t3)

# TODO: introduce get_component_list() and get_components_set() as for the indices?

def get_lines(ex, index_type):
    # TODO: add test for `get_lines`
    """
    returns ``(lines, traces, rest)`` for an index type,
    where ``lines`` is the list of list of positions of a matrix line,
    ``traces`` is the list of list of traced matrix lines,
    ``rest`` is the rest of the elements ot the tensor.
    """
    def _join_lines(a):
        i = 0
        while i < len(a):
            x = a[i]
            xend = x[-1]
            hit = True
            while hit:
                hit = False
                for j in range(i + 1, len(a)):
                    if j >= len(a):
                        break
                    if a[j][0] == xend:
                        hit = True
                        x.extend(a[j][1:])
                        xend = x[-1]
                        a.pop(j)
            i += 1
        return a

    tids = ex._tids
    components = tids.components
    dt = {}
    for c in components:
        if c in dt:
            continue
        index_types = c.index_types
        a = []
        for i in range(len(index_types)):
            if index_types[i] is index_type:
                a.append(i)
        if len(a) > 2:
            raise ValueError('at most two indices of type %s allowed' % index_type)
        if len(a) == 2:
            dt[c] = a
    dum = tids.dum
    lines = []
    traces = []
    traces1 = []
    for p0, p1, c0, c1 in dum:
        if components[c0] not in dt:
            continue
        if c0 == c1:
            traces.append([c0])
            continue
        ta0 = dt[components[c0]]
        ta1 = dt[components[c1]]
        if p0 not in ta0:
            continue
        if ta0.index(p0) == ta1.index(p1):
            # case gamma(i,s0,-s1)in c0, gamma(j,-s0,s2) in c1;
            # to deal with this case one could add to the position
            # a flag for transposition;
            # one could write [(c0, False), (c1, True)]
            raise NotImplementedError
        # if p0 == ta0[1] then G in pos c0 is mult on the right by G in c1
        # if p0 == ta0[0] then G in pos c1 is mult on the right by G in c0
        ta0 = dt[components[c0]]
        b0, b1 = (c0, c1) if p0 == ta0[1]  else (c1, c0)
        lines1 = lines[:]
        for line in lines:
            if line[-1] == b0:
                if line[0] == b1:
                    n = line.index(min(line))
                    traces1.append(line)
                    traces.append(line[n:] + line[:n])
                else:
                    line.append(b1)
                break
            elif line[0] == b1:
                line.insert(0, b0)
                break
        else:
            lines1.append([b0, b1])

        lines = [x for x in lines1 if x not in traces1]
        lines = _join_lines(lines)
    rest = []
    for line in lines:
        for y in line:
            rest.append(y)
    for line in traces:
        for y in line:
            rest.append(y)
    rest = [x for x in range(len(components)) if x not in rest]

    return lines, traces, rest


def _free_indices(expr):
    if isinstance(expr, TensExpr):
        return expr.free_indices_set
    return set([])

def _dummy_indices(expr):
    if isinstance(expr, TensExpr):
        return expr.dummy_indices_set
    return set([])
