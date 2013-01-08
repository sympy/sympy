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
the Butler-Portugal algorithm for canonicalization.

If there is a (anti)symmetric metric, the indices can be raised and
lowered when the tensor is put in canonical form.



TODO:

* better integration in SymPy

* term collection

* introduce scalars:
for example in ``p(a)*p(b)/(p**2 + m**2)`` the scalar ``p**2``
is equivalent to the tensor ``p(c)*p(-c)``, but can appear in the denominator.

* introduce tensor symmetrizers and algebraic operations on them

* Rule for contraction of Levi-Civita ``epsilon`` tensors:
one must intoduce generalized Kronecker deltas to do this properly

* develop a Young tableaux model:
one of its uses is in ``tensorsymmetry`` to provide the symmetry of the tensor
using the Young diagrams
"""


from collections import defaultdict
from sympy.core import Basic, sympify, Add, Mul, S
from sympy.core.symbol import Symbol, symbols
from sympy.combinatorics.tensor_can import get_symmetric_group_sgs, bsgs_direct_product, canonicalize, riemann_bsgs

class _TensorManager(object):
    def __init__(self):
        self._comm = defaultdict(dict)
        self._comm[0][0] = 0
        self._comm[0][1] = 0
        self._comm[1][0] = 1

    def set_comm(self, i, j, c):
        """
        set the commutation parameter ``c`` for commutation groups ``i, j``

        ``i, j`` they can be symbols or numbers, apart from ``0`` and ``1``
        which are reserved respectively for commuting and anticommuting
        tensors

        ``c``     commutation number
        0        commuting
        1        anticommuting
        None     no commutation property

        Tensors in the group ``i=0`` commute with any other tensor.
        Tensors in the group ``i=1`` anticommute within that group.
        For the remaining cases, use this method to set the commutation rules;
        by default ``c=None``.

        Examples
        ========

        ``G`` and ``GH`` do not commute with themselves and commute with
        each other; A is commuting.

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead, TensorManager
        >>> Lorentz = TensorIndexType('Lorentz')
        >>> i0,i1,i2,i3,i4 = tensor_indices('i0:5', Lorentz)
        >>> A = tensorhead('A', [Lorentz], [[1]])
        >>> G = tensorhead('G', [Lorentz], [[1]], 2)
        >>> GH = tensorhead('GH', [Lorentz], [[1]], 3)
        >>> TensorManager.set_comm(2, 3, 0)
        >>> (GH(i1)*G(i0)).canon_bp()
        G(i0)*GH(i1)
        >>> (G(i1)*G(i0)).canon_bp()
        G(i1)*G(i0)
        >>> (G(i1)*A(i0)).canon_bp()
        A(i0)*G(i1)
        """
        if i not in self._comm.keys():
            self._comm[i][0] = 0
            self._comm[0][i] = 0
        if j not in self._comm.keys():
            self._comm[0][j] = 0
        if j not in self._comm[i]:
            self._comm[i][j] = c
        if i not in self._comm[j]:
            self._comm[j][i] = c

    def set_comms(self, *args):
        for i, j, c in range(len(args)):
            set_comm(self, i, j, c)

    def get_comm(self, i, j):
        """
        Return the commutation parameter for commutation groups ``i, j``

        see ``_TensorManager.set_comm``
        """
        if i == 0 or j == 0:
            return 0
        return self._comm[i].get(j)


    def clear(self):
        for i in range(2, len(self._comm)):
            self._comm[i].clear()

TensorManager = _TensorManager()

class TensorIndexType(Basic):
    """
    A TensorIndexType is characterized by its name and its metric.

    ``metric = False`` symmetric metric (in Riemannian geometry)

    ``metric = True`` antisymmetric metric (for spinor calculus)

    In these two cases the metric is used to raise and lower indices.

    In the case of antisymmetric metric, the following raising and
    lowering conventions will be adopted:

    ``psi(a) = g(a, b)*psi(-b); chi(-a) = chi(b)*g(-b, -a)

    ``g(-a, b) = delta(-a, b); g(b, -a) = -delta(a, -b)``

    where ``delta(-a, b) = delta(b, -a)`` is the Kronecker delta

    ``metric = None``  there is no metric;
    it is not possible to raise or lower indices;
    e.g. the index of the defining representation of ``SU(N)``
    is 'covariant' and the conjugate representation is
    'contravariant'; for ``N > 2`` they are linearly independent.

    ``metric`` can be an object having ``name`` and ``antisym`` attributes.

    If a dimension ``dim`` is defined, it can be a symbol or an integer.
    """
    def __new__(cls, name, metric=False, dim=None, eps_dim = None,
                 dummy_fmt=None):
        """
        name   name of the tensor type

        ``metric``: it can be True, False, None or another object
        If it is True, False, None it gives its antisymmetry:
        False      symmetric metric
        True       antisymmetric
        None       no metric

        Otherwise, ``metric`` must have attributes ``name`` and ``antisym``

        ``dim``    dimension, it can be a symbol or a positive integer

        ``eps_dim``  dimension of the epsilon tensor; it is by default
                 equal to dim, if the latter is an integer; else
                 it can be assigned (for use in naive dimensional
                 regularization)

        ``dummy_fmt`` name of the head of dummy indices; by default it is
        the name of the tensor type

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        """
        obj = Basic.__new__(cls, name, metric)
        if not dummy_fmt:
            obj.dummy_fmt = '%s_%%d' % obj.name
        else:
            obj.dummy_fmt = '%s_%%d' % dummy_fmt
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

        obj.dim = dim
        obj.delta = obj.get_kronecker_delta()
        obj.eps_dim = eps_dim if eps_dim else dim
        obj.epsilon = obj.get_epsilon()
        return obj

    name = property(lambda self: self.args[0])


    def get_kronecker_delta(self):
        sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
        S2 = TensorType([self]*2, sym2)
        delta = S2('KD')
        return delta

    def get_epsilon(self):
        if not isinstance(self.eps_dim, int):
            return None
        sym = TensorSymmetry(get_symmetric_group_sgs(self.eps_dim, 1))
        Sdim = TensorType([self]*self.eps_dim, sym)
        epsilon = Sdim('Eps')
        return epsilon

    def __lt__(self, other):
        return self.name < other.name

    def __str__(self):
        return self.name

    __repr__ = __str__

class TensorIndex(Basic):
    """
    Tensor indices are contructed with the Einstein summation convention.

    A TensorIndex is chacterized by their ``name``, ``tensortype``
    and ``is_up``.

    An index can be in contravariant or in covariant form; in the latter
    case it is represented prepending a ``-`` to the index name.

    Dummy indices have a name with head given by ``tensortype.dummy_fmt``


    Examples
    ========

    >>> from sympy.tensor.tensor import TensorIndexType, TensorIndex, TensorSymmetry, TensorType, get_symmetric_group_sgs
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> i = TensorIndex('i', Lorentz); i
    i
    >>> sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
    >>> S1 = TensorType([Lorentz], sym1)
    >>> A, B = S1('A,B')
    >>> A(i)*B(-i)
    A(L_0)*B(-L_0)
    """
    def __new__(cls, name, tensortype, is_up=True):

        obj = Basic.__new__(cls, name, tensortype, is_up)
        obj.name = name
        obj.tensortype = tensortype
        obj.is_up = is_up
        return obj

    def _pretty(self):
        s = self.name
        if not self.is_up:
            s = '-%s' % s
        return s

    def __lt__(self, other):
        return (self.tensortype, self.name) < (other.tensortype, other.name)

    def __neg__(self):
        t1 = TensorIndex(self.name, self.tensortype,
                (not self.is_up))
        return t1

def tensor_indices(s, typ):
    """
    Returns list of tensor indices given their names and the type ``typ``

    ``s`` string of comma separated names of indices

    ``typ`` list of TensorIndexType of the indices

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

    return [TensorIndex(i, typ) for i in a]


class TensorSymmetry(Basic):
    """
    Symmetry of a tensor

    bsgs tuple (base, sgs) BSGS of the symmetry of the tensor

    Examples
    ========

    Define a symmetric tensor

    >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, TensorType, get_symmetric_group_sgs
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> V = S2('V')
    """
    def __new__(cls, bsgs, **kw_args):
        base, generators = bsgs
        obj = Basic.__new__(cls, base, generators, **kw_args)
        return obj

    base = property(lambda self: self.args[0])
    generators = property(lambda self: self.args[1])
    rank = property(lambda self: self.args[1][0].size)

def tensorsymmetry(*args):
    """
    return a ``TensorSymmetry`` object

    One can represent a tensor with any slot symmetry group using a BSGS
    ``args`` can be a BSGS
    ``args[0]``    base
    ``args[1]``    sgs

    Usually tensors are in (direct products of) irreducible representations
    of the symmetric group;
    ``args`` can be a list of lists representing Young tableaux
    ``[[1]]``       vector
    ``[[1]*n]``     symmetric tensor of rank ``n``
    ``[[n]``        antisymmetric tensor of rank ``n``
    ``[[2, 2]]``    slot symmetry of the Riemann tensor
    ``[[1],[1]]``   vector*vector

    TODO implement this with arbitrary Young tableaux

    Examples
    ========

    Symmetric tensor using a BSBS

    Symmetric tensor using a Young tableau

    >>> from sympy.tensor.tensor import TensorIndexType, TensorType, tensorsymmetry
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> sym2 = tensorsymmetry([1, 1])
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> V = S2('V')

    Symmetric tensor using a BSGS
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
        return TensorSymmetry([[], [Permutation(2)]])
    if len(args) == 2 and isinstance(args[1][0], Permutation):
        return TensorSymmetry(args)
    base, sgs = tableau2bsgs(args[0])
    for a in args[1:]:
        basex, sgsx = tableau2bsgs(a)
        base, sgs = bsgs_direct_product(base, sgs, basex, sgsx)
    return TensorSymmetry((base, sgs))


class TensorType(Basic):
    """
    A TensorType object is characterised by its index types and its symmetry

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
        assert symmetry.rank == len(index_types) + 2
        obj = Basic.__new__(cls, index_types, symmetry, **kw_args)
        return obj

    index_types = property(lambda self: self.args[0])

    symmetry = property(lambda self: self.args[1])

    types = property(lambda self: sorted(set(self.index_types), key=lambda x: x.name))

    def __str__(self):
        return 'TensorType(%s)' %([str(x) for x in self.index_types])

    def __call__(self, s, comm=0):
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
            return TensorHead(names[0], self, comm)
        else:
            return [TensorHead(name, self, comm) for name in names]

def tensorhead(name, typ, sym, comm=0):
    """
    Function generating tensorhead(s).

    ``name`` name or sequence of names (as in ``symbol``)

    ``typ``  index types

    ``sym``  same as ``*args`` in ``tensorsymmetry``

    ``comm``: commutation group number
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
    return S(name, comm)


class TensorHead(Basic):
    is_commutative = False

    def __new__(cls, name, typ, comm, **kw_args):
        """
        tensor with given name, index types, symmetry, commutation rule

        ``name`` name of the tensor

        ``typ`` list of TensorIndexType

        ``comm`` commutation group number
        see ``_TensorManager.set_comm``


        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensorsymmetry, TensorType
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> sym2 = tensorsymmetry([1]*2)
        >>> S2 = TensorType([Lorentz]*2, sym2)
        >>> A = S2('A')
        """
        assert isinstance(name, basestring)

        obj = Basic.__new__(cls, name, typ, **kw_args)
        obj.rank = len(obj.index_types)
        obj.types = typ.types
        obj.symmetry = typ.symmetry
        obj.comm = comm
        return obj

    name = property(lambda self: self.args[0])
    index_types = property(lambda self: self.args[1].index_types)

    def __lt__(self, other):
        return (self.name, self.index_types) < (other.name, other.index_types)

    def commutes_with(self, other):
        """
        Returns 0 (1) if self and other (anti)commute.

        Returns None if self and other do not (anti)commute.
        """
        r = TensorManager.get_comm(self.comm, other.comm)
        return r


    def __str__(self):
        return '%s(%s)' %(self.name, ','.join([str(x) for x in self.index_types]))
    __repr__ = __str__

    def __call__(self, *indices):
        """
        Returns a tensor with indices.

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, tensorhead
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> a, b = tensor_indices('a,b', Lorentz)
        >>> A = tensorhead('A', [Lorentz]*2, [[1]*2])
        >>> t = A(a, -b)
        """
        if not [indices[i].tensortype for i in range(len(indices))] == self.index_types:
            raise ValueError('wrong index type')
        components = [self]
        free, dum =  TensMul.from_indices(*indices)
        free.sort(key=lambda x: x[0].name)
        dum.sort()
        return TensMul(S.One, components, free, dum)


class TensExpr(Basic):
    """
    A tensor expression is an expression formed by tensors;
    currently the sums of tensors are distributed.

    A TensExpr can be a TensAdd or a TensMul.

    TensAdd objects are put in canonic form using the Butler-Portugal
    algorithm for canonicalization under monoterm symmetries.

    TensMul objects are formed by products of component tensors,
    and include a coefficient, which is a SymPy expression.


    In the internal representation contracted indices are represented
    by ``(ipos1, ipos2, icomp1, icomp2)``, where ``icomp1`` is the position
    of the component tensor with contravariant index, ``ipos1`` is the
    slot which the index occupies in that component tensor.

    Contracted indices are therefore nameless in the internal representation.
    """

    _op_priority = 11.0
    is_TensMul = False
    is_TensAdd = False
    is_commutative = False

    def __neg__(self):
        return (-1)*self

    def __abs__(self):
        raise NotImplementedError

    def __add__(self, other):
        other = sympify(other)
        if self.is_TensAdd:
            args = self.args + (other,)
            return TensAdd(*args)
        return TensAdd(self, other)

    def __radd__(self, other):
        return TensAdd(other, self)

    def __sub__(self, other):
        other = sympify(other)
        return TensAdd(self, -other)

    def __rsub__(self, other):
        return TensAdd(other, -self)

    def __mul__(self, other):
        if self.is_TensAdd:
            return TensAdd(*[x*other for x in self.args])
        if other.is_TensAdd:
            return TensAdd(*[self*x for x in other.args])
        return TensMul.__mul__(self, other)

    def __rmul__(self, other):
        if self.is_TensMul:
            coeff = other*self._coeff
            return TensMul(coeff, self._components, self._free, self._dum)
        return TensAdd(*[x*other for x in self.args])

    def __pow__(self, other):
        raise NotImplementedError

    def __rpow__(self, other):
        raise NotImplementedError

    def __div__(self, other):
        other = sympify(other)
        if isinstance(other, TensExpr):
            raise ValueError('cannot divide by a tensor')
        coeff = self._coeff/other
        return TensMul(coeff, self._components, self._free, self._dum, is_canon_bp=self._is_canon_bp)


    def __rdiv__(self, other):
        raise NotImplementedError()

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def substitute_indices(self, *index_tuples):
        """
        Return a tensor with free indices substituted according to ``index_tuples``

        ``index_types`` list of tuples ``(old_index, new_index)``

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
        if self.is_TensMul:
            free = self._free
            free1 = []
            for j, ipos, cpos in free:
                for i, v in index_tuples:
                    if i.name == j.name and i.tensortype == j.tensortype:
                        if i.is_up == j.is_up:
                            free1.append((v, ipos, cpos))
                        else:
                            free1.append((-v, ipos, cpos))
                        break
                else:
                    free1.append((j, ipos, cpos))

            return TensMul(self._coeff, self._components, free1, self._dum)
        if self.is_TensAdd:
            args = self.args
            args1 = []
            for x in args:
                y = x.substitute_indices(*index_tuples)
                args1.append(y)
            return TensAdd(*args1)


def _tensAdd_collect_terms(args):
    """
    collect TensMul terms differing at most by their coefficient
    """
    a = []
    pprev = None
    prev = args[0]
    prev_coeff = prev._coeff
    changed = False
    new = 0

    for x in args[1:]:
        # if x and prev have the same tensor, update the coeff of prev
        if x._components == prev._components \
                and x._free == prev._free and x._dum == prev._dum:
            prev_coeff = prev_coeff + x._coeff
            changed = True
            op = 0
        else:
            # x and prev are different; if not changed, prev has not
            # been updated; store it
            if not changed:
                a.append(prev)
            else:
                # get a tensor from prev with coeff=prev_coeff and store it
                if prev_coeff:
                    t = TensMul(prev_coeff, prev._components,
                        prev._free, prev._dum)
                    a.append(t)
            # move x to prev
            op = 1
            pprev, prev = prev, x
            pprev_coeff, prev_coeff = prev_coeff, x._coeff
            changed = False
    # if the case op=0 prev was not stored; store it now
    # in the case op=1 x was not stored; store it now (as prev)
    if op == 0 and prev_coeff:
        prev = TensMul(prev_coeff, prev._components, prev._free, prev._dum)
        a.append(prev)
    elif op == 1:
        a.append(prev)
    return a

def _tensAdd_flatten(args):
    """
    flatten TensAdd, coerce terms which are not tensors to tensors
    """
    if not all(isinstance(x, TensExpr) for x in args):
        args1 = []
        for x in args:
            if isinstance(x, TensExpr):
                if x.is_TensAdd:
                    args1.extend(list(x.args))
                else:
                    args1.append(x)
        args1 = [x for x in args1 if isinstance(x, TensExpr) and x._coeff]
        args2 = [x for x in args if not isinstance(x, TensExpr)]
        t0 = args1[0]
        if t0.is_TensAdd:
            t0 = t0.args[0]
        if t0._free:
           raise ValueError('all tensors must have the same indices')
        t1 = TensMul(Add(*args2), [], [], [])
        args = [t1] + args1
    a = []
    for x in args:
        if x.is_TensAdd:
            a.extend(list(x.args))
        else:
            a.append(x)
    args = [x for x in a if x._coeff]
    return args

def _tensAdd_check(args):
    """
    check that all addends have the same free indices
    """
    indices0 = sorted([x[0] for x in args[0]._free], key=lambda x: x.name)
    list_indices = [sorted([y[0] for y in x._free], key=lambda x: x.name) for x in args[1:]]
    if not all(x == indices0 for x in list_indices):
        raise ValueError('all tensors must have the same indices')


class TensAdd(TensExpr):
    """
    Sum of tensors

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
    """
    is_TensAdd = True

    def __new__(cls, *args, **kw_args):
        """
        args  tuple of addends
        """
        args = [sympify(x) for x in args if x]
        args = _tensAdd_flatten(args)
        if not args:
            return S.Zero

        _tensAdd_check(args)
        obj = Basic.__new__(cls, **kw_args)
        if len(args) == 1 and args[0].is_TensMul:
            obj._args = tuple(args)
            return obj
        args = [x.canon_bp() for x in args if x]
        args = [x for x in args if x]
        if not args:
            return S.Zero

        # collect canonicalized terms
        args.sort(key=lambda x: (x._components, x._free, x._dum))
        a = _tensAdd_collect_terms(args)
        if not a:
            return S.Zero
        # it there is only a component tensor return it
        if len(a) == 1:
            return a[0]
        obj._args = tuple(a)
        return obj

    @property
    def free_args(self):
        return self.args[0].free_args


    def __call__(self, *indices):
        """Returns tensor with ordered free indices replaced by ``indices``

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
        free_args = self.free_args
        indices = list(indices)
        if [x.tensortype for x in indices] != [x.tensortype for x in free_args]:
            raise ValueError('incompatible types')
        if indices == free_args:
            return self
        index_tuples = zip(free_args, indices)
        a = [x.fun_eval(*index_tuples) for x in self.args]
        res = TensAdd(*a)

        return res

    def canon_bp(self):
        """
        canonicalize using the Butler-Portugal algorithm for canonicalization
        under monoterm symmetries.
        """
        args = [x.canon_bp() for x in self.args]
        res = TensAdd(*args)
        return res

    def __eq__(self, other):
        other = sympify(other)
        if not isinstance(other, TensExpr):
            if len(self.args) == 1:
                return self.args[0]._coeff == other
        if isinstance(other, TensExpr) and other.is_TensMul and other._coeff == 0:
            return self == 0
        t = self - other
        if not isinstance(t, TensExpr):
            return t == 0
        else:
            if t.is_TensMul:
                return t._coeff == 0
            else:
                return all(x._coeff == 0 for x in t.args)

    def __ne__(self, other):
        return not (self == other)

    def contract_delta(self, delta):
        args = [x.contract_delta(delta) for x in self.args]
        if len(args) == 1:
            return args[0]
        t = TensAdd(*args)
        return canon_bp(t)

    def contract_metric(self, g, contract_all=False):
        """
        Raise or lower indices with the metric ``g``

        ``g``  metric

        ``contract_all`` if True, eliminate all ``g`` which are contracted
        """

        args = [x.contract_metric(g, contract_all) for x in self.args]
        if len(args) == 1:
            return args[0]
        t = TensAdd(*args)
        return canon_bp(t)


    def fun_eval(self, *index_tuples):
        """
        Return a tensor with free indices substituted according to ``index_tuples``

        ``index_types`` list of tuples ``(old_index, new_index)``

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
        args = self.args
        args1 = []
        for x in args:
            y = x.fun_eval(*index_tuples)
            args1.append(y)
        return TensAdd(*args1)

    def _pretty(self):
        a = []
        args = self.args
        for x in args:
            a.append(str(x))
        a.sort()
        s = ' + '.join(a)
        s = s.replace('+ -', '- ')
        return s

class TensMul(TensExpr):
    """
    Product of tensors
    """
    is_TensMul = True

    def __new__(cls, coeff, *args, **kw_args):
        """

        coeff SymPy expression coefficient of the tensor.

        ``args[0]``   list of TensorHead of the component tensors.

        ``args[1]``   list of ``(ind, ipos, icomp)``
        where ``ind`` is a free index, ``ipos`` is the slot position
        of ``ind`` in the ``icomp``-th component tensor.

        ``args[2]`` list of tuples representing dummy indices.
        (ipos1, ipos2, icomp1, icomp2) indicates that the contravariant
        dummy index is the ``ipos1``-th slot position in the ``icomp1``-th
        component tensor; the corresponding covariant index is
        in the ``ipos2`` slot position in the ``icomp2``-th component tensor.
        """
        obj = Basic.__new__(cls)
        obj._components = args[0]
        obj.types = []
        for t in obj._components:
            obj.types.extend(t.types)
        obj._free = args[1]
        obj._dum = args[2]
        obj.rank = len(obj._free) + 2*len(obj._dum)
        obj._coeff = coeff
        obj._is_canon_bp = kw_args.get('is_canon_bp', False)

        return obj

    @property
    def free_args(self):
        return sorted([x[0] for x in self._free])

    def __eq__(self, other):
        if other == 0:
            return self._coeff == 0
        other = sympify(other)
        if not isinstance(other, TensExpr):
            assert not self._components
            return self._coeff == other
        res = self - other
        return res == 0

    def __ne__(self, other):
        return not self == other

    @staticmethod
    def from_indices(*indices):
        """
        Convert ``indices`` into ``free``, ``dum`` for single component tensor

        ``free``     list of tuples ``(index, pos, 0)``,
                     where ``pos`` is the position of index in
                     the list of indices formed by the component tensors

        ``dum``      list of tuples ``(pos_contr, pos_cov, 0, 0)``

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensMul
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2, m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
        >>> TensMul.from_indices(m0, m1, -m1, m3)
        ([(m0, 0, 0), (m3, 3, 0)], [(1, 2, 0, 0)])
        """
        n = len(indices)
        if n == 1:
            return [(indices[0], 0, 0)], []

        # find the positions of the free indices and of the dummy indices
        free = [True]*len(indices)
        index_dict = {}
        dum = []
        for i, index in enumerate(indices):
            name = index.name
            typ = index.tensortype
            contr = index.is_up
            if (name, typ) in index_dict:
                # found a pair of dummy indices
                is_contr, pos = index_dict[(name, typ)]
                # check consistency and update free
                if is_contr:
                    if contr:
                        raise ValueError('two equal contravariant indices in slots %d and %d' %(pos, i))
                    else:
                        free[pos] = False
                        free[i] = False
                else:
                    if contr:
                        free[pos] = False
                        free[i] = False
                    else:
                        raise ValueError('two equal covariant indices in slots %d and %d' %(pos, i))
                if contr:
                    dum.append((i, pos, 0, 0))
                else:
                    dum.append((pos, i, 0, 0))
            else:
                index_dict[(name, typ)] = index.is_up, i

        free = [(index, i, 0) for i, index in enumerate(indices) if free[i]]
        free.sort()
        return free, dum

    def get_indices(self):
        """
        Returns the list of indices of the tensor

        The indices are listed in the order in which they appear in the
        component tensors.

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
        indices = [None]*self.rank
        start = 0
        pos = 0
        vpos = []
        components = self._components
        for t in components:
            vpos.append(pos)
            pos += t.rank
        cdt = defaultdict(int)
        # if the free indices have names with dummy_fmt, start with an
        # index higher than those for the dummy indices
        # to avoid name collisions
        for indx, ipos, cpos in self._free:
            if indx.name.split('_')[0] == indx.tensortype.dummy_fmt[:-3]:
                cdt[indx.tensortype] = max(cdt[indx.tensortype], int(indx.name.split('_')[1]) + 1)
            start = vpos[cpos]
            indices[start + ipos] = indx
        for ipos1, ipos2, cpos1, cpos2 in self._dum:
            start1 = vpos[cpos1]
            start2 = vpos[cpos2]
            typ1 = components[cpos1].index_types[ipos1]
            assert typ1 == components[cpos2].index_types[ipos2]
            fmt = typ1.dummy_fmt
            nd = cdt[typ1]
            indices[start1 + ipos1] = TensorIndex(fmt % nd, typ1)
            indices[start2 + ipos2] = TensorIndex(fmt % nd, typ1, False)
            cdt[typ1] += 1
        return indices

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
        indices = self.get_indices()
        pos = 0
        components = self._components
        if not components:
            return [TensMul(self._coeff, [], [], [])]
        res = []
        for t in components:
            t1 = t(*indices[pos:pos + t.rank])
            pos += t.rank
            res.append(t1)
        res[0] = TensMul(self._coeff, res[0]._components, res[0]._free, res[0]._dum, is_canon_bp=res[0]._is_canon_bp)
        return res

    def canon_args(self):
        """
        Returns ``(g, dummies, msym, v)``, the entries of ``canonicalize``

        see ``canonicalize`` in ``tensor_can.py``
        """
        # to be called after sorted_components
        from sympy.combinatorics.permutations import _af_new
        types = list(set(self.types))
        types.sort(key = lambda x: x.name)
        n = self.rank
        g = [None]*n + [n, n+1]
        pos = 0
        vpos = []
        components = self._components
        for t in components:
            vpos.append(pos)
            pos += t.rank
        # ordered indices: first the free indices, ordered by types
        # then the dummy indices, ordered by types and contravariant before
        # covariant
        # g[position in tensor] = position in ordered indices
        for i, (indx, ipos, cpos) in enumerate(self._free):
            pos = vpos[cpos] + ipos
            g[pos] = i
        pos = len(self._free)
        j = len(self._free)
        dummies = []
        prev = None
        a = []
        msym = []
        for ipos1, ipos2, cpos1, cpos2 in self._dum:
            pos1 = vpos[cpos1] + ipos1
            pos2 = vpos[cpos2] + ipos2
            g[pos1] = j
            g[pos2] = j + 1
            j += 2
            typ = components[cpos1].index_types[ipos1]
            if typ != prev:
                if a:
                    dummies.append(a)
                a = [pos, pos + 1]
                prev = typ
                msym.append(typ.metric_antisym)
            else:
                a.extend([pos, pos + 1])
            pos += 2
        if a:
            dummies.append(a)
        numtyp = []
        prev = None
        for t in components:
            if t == prev:
                numtyp[-1][1] += 1
            else:
                prev = t
                numtyp.append([prev, 1])
        v = []
        for h, n in numtyp:
            if h.comm == 0 or h.comm == 1:
                comm = h.comm
            else:
                comm = TensorManager.get_comm(h.comm, h.comm)

            v.append((h.symmetry.base, h.symmetry.generators, n, comm))
        return _af_new(g), dummies, msym, v

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
            coeff = self._coeff*other
            return TensMul(coeff, self._components, self._free, self._dum, is_canon_bp=self._is_canon_bp)
        if other.is_TensAdd:
            return TensAdd(*[self*x for x in other.args])

        components = self._components + other._components
        # find out which free indices of self and other are contracted
        free_dict1 = dict([(i.name, (pos, cpos, i)) for i, pos, cpos in self._free])
        free_dict2 = dict([(i.name, (pos, cpos, i)) for i, pos, cpos in other._free])

        free_names = set(free_dict1.keys()) & set(free_dict2.keys())
        # find the new `free` and `dum`
        nc1 = len(self._components)
        dum2 = [(i1, i2, c1 + nc1, c2 + nc1) for i1, i2, c1, c2 in other._dum]
        free1 = [(ind, i, c) for ind, i, c in self._free if ind.name not in free_names]
        free2 = [(ind, i, c + nc1) for ind, i, c in other._free if ind.name not in free_names]
        free = free1 + free2
        dum = self._dum + dum2
        for name in free_names:
            ipos1, cpos1, ind1 = free_dict1[name]
            ipos2, cpos2, ind2 = free_dict2[name]
            cpos2 += nc1
            if ind1.is_up == ind2.is_up:
                raise ValueError('wrong index contruction %s' % ind1)
            if ind1.is_up:
                new_dummy = (ipos1, ipos2, cpos1, cpos2)
            else:
                new_dummy = (ipos2, ipos1, cpos2, cpos1)
            dum.append(new_dummy)
        coeff = self._coeff*other._coeff
        return TensMul(coeff, components, free, dum)


    def sorted_components(self):
        """
        Returns a tensor with sorted components
        """
        from sympy.combinatorics.permutations import _af_invert
        cv = zip(self._components, range(len(self._components)))
        sign = 1
        n = len(cv) - 1
        for i in range(n):
            for j in range(n, i, -1):
                c = cv[j-1][0].commutes_with(cv[j][0])
                if c not in [0, 1]:
                    continue
                if (cv[j-1][0].types, cv[j-1][0].name) > \
                        (cv[j][0].types, cv[j][0].name):
                    cv[j-1], cv[j] = cv[j], cv[j-1]
                    if c:
                        sign = -sign

        # perm_inv[new_pos] = old_pos
        components = [x[0] for x in cv]
        perm_inv = [x[1] for x in cv]
        perm = _af_invert(perm_inv)
        free = [(ind, i, perm[c]) for ind, i, c in self._free]
        free.sort()
        dum = [(i1, i2, perm[c1], perm[c2]) for i1, i2, c1, c2 in self._dum]
        dum.sort(key = lambda x: components[x[2]].index_types[x[0]])
        coeff = -self._coeff if sign == -1 else self._coeff
        t = TensMul(coeff, components, free, dum)
        return t

    def perm2tensor(self, g, canon_bp=False):
        """
        Returns the tensor corresponding to the permutation ``g``

        ``g``  permutation corrisponding to the tensor in the representation
        used in canonicalization

        ``canon_bp``   if True, then ``g`` is the permutation
        corresponding to the canonical form of the tensor
        """
        from bisect import bisect_right
        vpos = []
        components = self._components
        pos = 0
        for t in components:
            vpos.append(pos)
            pos += t.rank
        sorted_free = [x[0] for x in self._free]
        sorted_free.sort()
        nfree = len(sorted_free)
        rank = self.rank
        indices = [None]*rank
        dum = [[None]*4 for i in range((rank - nfree)//2)]
        free = []
        icomp = -1
        for i in range(rank):
            if i in vpos:
                icomp += vpos.count(i)
                pos0 = i
            ipos = i - pos0
            gi = g[i]
            if gi < nfree:
                ind = sorted_free[gi]
                free.append((ind, ipos, icomp))
            else:
                j = gi - nfree
                idum, cov = divmod(j, 2)
                if cov:
                    dum[idum][1] = ipos
                    dum[idum][3] = icomp
                else:
                    dum[idum][0] = ipos
                    dum[idum][2] = icomp
        dum = [tuple(x) for x in dum]
        coeff = self._coeff
        if g[-1] != len(g) - 1:
            coeff = -coeff
        res = TensMul(coeff, components, free, dum, is_canon_bp=canon_bp)
        return res

    def canon_bp(self):
        """
        canonicalize using the Butler-Portugal algorithm for canonicalization
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
        from sympy.combinatorics.tensor_can import canonicalize
        if self._is_canon_bp:
            return self
        if not self._components:
            return self
        t = self.sorted_components()
        g, dummies, msym, v = t.canon_args()
        can = canonicalize(g, dummies, msym, *v)
        if can == 0:
            return S.Zero
        return t.perm2tensor(can, True)

    def _contract(self, g, antisym, contract_all=False):
        """
        helper method for ``contract_metric`` and ``contract_delta``

        ``g`` metric to be contracted

        ``antisym``:
        False  symmetric metric
        True   antisymmetric metric
        None   delta
        """
        if not self._components:
            return self
        free_indices = [x[0] for x in self._free]
        a = self.split()
        typ = g.index_types[0]
        for i, tg in enumerate(a):
            if tg._components[0] == g:
                tg_free = [x[0] for x in tg._free]
                if len(tg_free) == 0:
                    t = _contract_g_with_itself(a, i, tg, tg_free, g, antisym)
                    if contract_all == True and g in t._components:
                        return t._contract(g, antisym, True)
                    return t

                if all(indx in free_indices for indx in tg_free):
                    continue
                else:
                    break
        else:
            # all metric tensors have only free indices, there is no contraction
            return self

        # tg has one or two indices contracted with other tensors
        # i position of tg in a
        coeff = S.One
        tg_free = tg._free
        if antisym:
            # order by slot position
            tg_free = sorted(tg_free, key=lambda x: x[1])

        if tg_free[0][0] in free_indices or tg_free[1][0] in free_indices:
            # tg has one free index
            res = _contract_g_with_free_index(a, free_indices, i, tg, tg_free, g, antisym)
        else:
            # tg has two indices contracted with other tensors
            res = _contract_g_without_free_index(a, free_indices, i, tg, tg_free, g, typ, antisym)
        if contract_all == True and g in res._components:
            return res._contract(g, antisym, True)
        return res

    def contract_delta(self, delta):
        typ = delta.types[0]
        t = self._contract(delta, None, True)
        return t

    def contract_metric(self, g, contract_all=False):
        """
        Raise or lower indices with the metric ``g``

        ``g``  metric

        ``contract_all`` if True, eliminate all ``g`` which are contracted

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
        if g.index_types[0].metric_antisym is None:
            raise NotImplementedError
        return self._contract(g, g.index_types[0].metric_antisym, contract_all)


    def fun_eval(self, *index_tuples):
        """
        Return a tensor with free indices substituted according to ``index_tuples``

        ``index_types`` list of tuples ``(old_index, new_index)``

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
        free = self._free
        free1 = []
        for j, ipos, cpos in free:
            # search j in index_tuples
            for i, v in index_tuples:
                if i == j:
                    free1.append((v, ipos, cpos))
                    break
            else:
                free1.append((j, ipos, cpos))
        return TensMul(self._coeff, self._components, free1, self._dum)


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
        free_args = self.free_args
        indices = list(indices)
        if [x.tensortype for x in indices] != [x.tensortype for x in free_args]:
            raise ValueError('incompatible types')
        if indices == free_args:
            return self
        t = self.fun_eval(*zip(free_args, indices))
        return t


    def _pretty(self):
        if self._components == []:
            return str(self._coeff)
        indices = [str(ind) for ind in self.get_indices()]
        pos = 0
        a = []
        for t in self._components:
            if t.rank > 0:
                a.append('%s(%s)' % (t.name, ', '.join(indices[pos:pos + t.rank])))
            else:
                a.append('%s' % t.name)
            pos += t.rank
        res = '*'. join(a)
        if self._coeff == S.One:
            return res
        elif self._coeff == -S.One:
            return '-%s' % res
        if self._coeff.is_Atom:
            return '%s*%s' % (self._coeff, res)
        else:
            return '(%s)*%s' %(self._coeff, res)



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
        return TensMul(S.One, [], [], [])
    t = a[0]
    for tx in a[1:]:
        t = t*tx
    return t



def riemann_cyclic_replace(t_r):
    """
    replace Riemann tensor with an equivalent expression

    ``R(m,n,p,q) -> 2/3*R(m,n,p,q) - 1/3*R(m,q,n,p) + 1/3*R(m,p,n,q)``

    """
    free = sorted(t_r._free, key=lambda x: x[1])
    m, n, p, q = [x[0] for x in free]
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
    a = t2.args
    a1 = [x.split() for x in a]
    a2 = [[riemann_cyclic_replace(tx) for tx in y] for y in a1]
    a3 = [tensor_mul(*v) for v in a2]
    t3 = TensAdd(*a3)
    if not t3:
        return t3
    else:
        return canon_bp(t3)


def tensorlist_contract_metric(a, tg):
    """
    contract `tg` with a tensor in the list `a = t.split()`
    Only for symmetric metric.
    """
    ind1, ind2 = [x[0] for x in tg._free]
    mind1 = -ind1
    mind2 = -ind2
    for i in range(len(a)):
        t1 = a[i]
        for j in range(len(t1._free)):
            indx, ipos, _ = t1._free[j]
            if indx == mind1 or indx == mind2:
                ind3 = ind2 if indx == mind1 else ind1
                free1 = t1._free[:]
                free1[j] = (ind3, ipos, 0)
                t2 = TensMul(t1._coeff, t1._components, free1, t1._dum)
                a[i] = t2
                return a
    a.append(tg)
    return a

def _contract_g_with_itself(a, i, tg, tg_free, g, antisym):
    """
    helper function for _contract
    """
    typ = g.index_types[0]
    a1 = a[:i] + a[i + 1:]
    t11 = tensor_mul(*a1)
    if typ.dim is None:
        raise ValueError('dimension not assigned')
    coeff = typ.dim*a[i]._coeff
    if antisym and tg._dum[0][0] == 0:
        # g(i, -i) = -D
        coeff = -coeff
    t = tensor_mul(*a1)*coeff
    return t


def _contract_g_with_free_index(a, free_indices, i, tg, tg_free, g, antisym):
    """
    helper function for _contract
    """
    if tg_free[0][0] in free_indices:
        ind_free = tg_free[0][0]
        ind, ipos1, _ = tg_free[1]
    else:
        ind_free = tg_free[1][0]
        ind, ipos1, _ = tg_free[0]

    ind1 = -ind
    # search ind1 in the other component tensors
    for j, tx in enumerate(a):
        if ind1 in [x[0] for x in tx._free]:
            break
    # replace ind1 with ind_free
    free1 = []
    for indx, iposx, _ in tx._free:
        if indx == ind1:
            free1.append((ind_free, iposx, 0))
        else:
            free1.append((indx, iposx, 0))
    coeff = tx._coeff
    if antisym:
        if ind.is_up and ind == tg_free[0][0] or \
        (not ind.is_up) and ind == tg_free[1][0]:
            # g(i1, i0)*psi(-i1) = -psi(i0)
            # g(-i0, -i1)*psi(i1) = -psi(-i0)
            coeff = -coeff
    t1 = TensMul(coeff, tx._components, free1, tx._dum)
    a[j] = t1
    a = a[:i] + a[i + 1:]
    coeff = tg._coeff
    res = tensor_mul(*a)
    return coeff*res


def _contract_g_without_free_index(a, free_indices, i, tg, tg_free, g, typ, antisym):
    """
    helper function for _contract
    """
    coeff = S.One
    ind1 = tg_free[0][0]
    ind2 = tg_free[1][0]
    ind1m = -ind1
    ind2m = -ind2
    for k, ty in enumerate(a):
        if ind2m in [x[0] for x in ty._free]:
            break
    # ty has the index ind2m
    ty_free = ty._free[:]
    if ty._components == [g]:
        ty_indices = [x[0] for  x in ty._free]
        if all(x in [ind1m, ind2m] for x in ty_indices):
            # the two `g` are completely contracted
            # i < k always
            a = a[:i] + a[i+1:k] + a[k+1:]
            coeff = coeff*typ.dim*tg._coeff*ty._coeff
            if antisym:
                ty_free = sorted(ty_free, key=lambda x: x[1])
                if ind1.is_up == ind2.is_up:
                    # g(i,j)*g(-i,-j) = g(-i,-j)*g(i,j) = dim
                    # g(i,j)*g(-j,-i) = g(-i,-j)*g(j,i) = -dim
                    if ind1m == ty_free[1][0]:
                        coeff = -coeff
                else:
                    # g(-i,j)*g(i,-j) = g(i,-j)^g(-i,j) = -dim
                    # g(-i,j)*g(-j,i) = g(i,-j)*g(j,i) = dim
                    if ind1m == ty_free[0][0]:
                        coeff = -coeff

            if a:
                res = tensor_mul(*a)
                res = coeff*res
            else:
                res = TensMul(coeff, [],[],[], is_canon_bp=True)
            return res

    free2 = []
    ty_freeindices = [x[0] for x in ty_free]
    if ind1m in ty_freeindices:
        # tg has both indices contracted with ty
        free2 = [(indx, iposx, cposx) for indx, iposx, cposx in ty._free if indx != ind1m and indx != ind2m]
        dum2 = ty._dum[:]
        for indx, iposx, _ in ty._free:
            if indx == ind1m:
                iposx1 = iposx
            if indx == ind2m:
                iposx2 = iposx
        if antisym:
            if ind1.is_up == ind2.is_up:
                if iposx1 < iposx2:
                    coeff = -coeff
                    dum2.append((iposx1, iposx2, 0, 0))
                else:
                    dum2.append((iposx2, iposx1, 0, 0))
            else:
                if iposx1 > iposx2:
                    coeff = -coeff
                    dum2.append((iposx2, iposx1, 0, 0))
                else:
                    dum2.append((iposx1, iposx2, 0, 0))
        else:
            dum2.append((iposx1, iposx2, 0, 0))
    else:
        # replace ind2m with ind1 in the free indices of ty

        free2 = []
        if not antisym:
            for indx, iposx, _ in ty._free:
                if indx == ind2m:
                    free2.append((ind1, iposx, 0))
                else:
                    free2.append((indx, iposx, 0))
        else:
            for indx, iposx, _ in ty._free:
                if indx == ind2m:
                    free2.append((ind1, iposx, 0))
                    if indx.is_up:
                        coeff = -coeff
                else:
                    free2.append((indx, iposx, 0))
                    if not indx.is_up:
                        coeff = -coeff
        dum2 = ty._dum

    t2 = TensMul(ty._coeff, ty._components, free2, dum2)
    a[k] = t2
    a = a[:i] + a[i + 1:]
    coeff = coeff*tg._coeff
    res = tensor_mul(*a)
    return coeff*res
