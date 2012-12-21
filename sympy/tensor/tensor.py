from collections import defaultdict
from sympy.core import Basic, sympify, Add, Mul, S
from sympy.combinatorics.tensor_can import get_symmetric_group_sgs, bsgs_direct_product, canonicalize, riemann_bsgs

class TensorIndexType(object):
    """
    A TensorIndexType is characterized by its name, by metric_sym,
    giving the symmetry of its metric.

    If a dimension ``dim`` is defined, it can be a symbol or an integer.
    """
    def __init__(self, name, metric_sym=0, dim=None, eps_dim = None,
                 dummy_fmt=None):
        """
        name   name of the tensor type

        ``metric_sym``:
        0      symmetric
        1      antisymmetric
        None   no symmetry

        dim    dimension, it can be a symbol or a positive integer

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
        self.name = name
        self.metric_sym = metric_sym
        if not dummy_fmt:
            self.dummy_fmt = '%s_%%d' % self.name
        else:
            self.dummy_fmt = '%s_%%d' % dummy_fmt
        self.metric = self.get_metric()
        self.dim = dim
        if eps_dim:
            self.eps_dim = eps_dim
        else:
            self.eps_dim = dim
        if isinstance(self.eps_dim, int):
            self.epsilon = self.get_epsilon()
        else:
            self.epsilon = None

    def get_metric(self):
        if self.metric_sym is None:
            return None
        sym2 = TensorSymmetry(get_symmetric_group_sgs(2, self.metric_sym))
        S2 = TensorType([self]*2, sym2)
        metric = S2('metric')
        return metric

    def get_epsilon(self):
        sym = TensorSymmetry(get_symmetric_group_sgs(self.eps_dim, 1))
        Sdim = TensorType([self]*self.eps_dim, sym)
        epsilon = Sdim('Eps')
        return epsilon


    def __str__(self):
        return self.name

    __repr__ = __str__

class TensorIndex(object):
    """
    Tensor indices are contructed with the Einstein summation convention.

    A TensorIndex is chacterized by their ``name``, ``tensortype``
    and ``is_contravariant``.

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
    def __init__(self, name, tensortype, is_contravariant=True):
        self.name = name
        self.tensortype = tensortype
        self.is_contravariant = is_contravariant

    def __str__(self):
        s = self.name
        if not self.is_contravariant:
            s = '-%s' % s
        return s

    __repr__ = __str__

    def __eq__(self, other):
        return self.name == other.name and \
               self.tensortype == other.tensortype and \
               self.is_contravariant == other.is_contravariant

    def __ne__(self, other):
        return not (self == other)

    def __neg__(self):
        t1 = TensorIndex(self.name, self.tensortype,
                (not self.is_contravariant))
        return t1

def tensor_indices(s, typ):
    """
    Returns tensor indices given their names and the TensorIndexType ``typ``

    ``s`` string of comma separated names of indices

    ``typ`` list of TensorIndexType of the indices

    Examples
    ========

    >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> a, b, c, d = tensor_indices('a,b,c,d', Lorentz)
    """
    a = s.split(',')
    return [TensorIndex(i, typ) for i in a]


class TensorSymmetry(object):
    """
    Symmetry of a tensor

    bsgs tuple (base, sgs) BSGS of the symmetry of the tensor

    Examples
    ========

    Define a symmetric tensor

    >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, TensorType, get_symmetric_group_sgs
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> a, b, c, d = tensor_indices('a,b,c,d', Lorentz)
    >>> sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> V = S2('V')
    """
    def __init__(self, bsgs):
        self.base, self.generators = bsgs
        self.rank = self.generators[0].size


class TensorType(Basic):
    """
    A TensorType object is characterised by its index types and its symmetry

    Examples
    ========

    Define a symmetric tensor

    >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, TensorType, get_symmetric_group_sgs
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> a, b, c, d = tensor_indices('a,b,c,d', Lorentz)
    >>> sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
    >>> S2 = TensorType([Lorentz]*2, sym2)
    >>> V = S2('V')
    """
    is_commutative = False

    def __new__(cls, index_types, symmetry, **kw_args):
        obj = Basic.__new__(cls, **kw_args)
        obj.index_types = []
        obj.index_types = index_types
        obj.types = list(set(obj.index_types))
        obj.types.sort(key=lambda x: x.name)
        obj.symmetry = symmetry
        assert symmetry.rank == len(index_types) + 2
        return obj

    def __str__(self):
        return 'TensorType(%s)' %([str(x) for x in self.index_types])

    def __call__(self, s, commuting=0):
        """
        commuting:
        None no commutation rule
        0    commutes
        1    anticommutes

        Examples
        ========

        Define symmetric tensors ``V``, ``W`` and ``G``, respectively
        commuting, anticommuting and with no commutation symmetry

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, TensorType, get_symmetric_group_sgs, canon_bp
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> a, b = tensor_indices('a,b', Lorentz)
        >>> sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
        >>> S2 = TensorType([Lorentz]*2, sym2)
        >>> V = S2('V')
        >>> W = S2('W', 1)
        >>> G = S2('G', None)
        >>> canon_bp(V(a, b)*V(-b, -a))
        V(L_0, L_1)*V(-L_0, -L_1)
        >>> canon_bp(W(a, b)*W(-b, -a))
        0
        """
        names = s.split(',')
        if len(names) == 1:
            return TensorHead(names[0], self, commuting)
        else:
            return [TensorHead(name, self, commuting) for name in names]

class TensorHead(Basic):
    is_commutative = False

    def __new__(cls, name, typ, commuting, **kw_args):
        """
        tensor with given name, index types, symmetry, commutation rule

        ``name`` name of the tensor

        ``typ`` list of TensorIndexType

        ``commuting`` commutation property
        0     commuting tensor
        1     anticommuting tensor
        None  no commutation rule

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, TensorType, get_symmetric_group_sgs
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> a, b = tensor_indices('a,b', Lorentz)
        >>> sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
        >>> S2 = TensorType([Lorentz]*2, sym2)
        >>> A = S2('A')
        """
        assert isinstance(name, basestring)

        obj = Basic.__new__(cls, **kw_args)
        obj.name = name
        obj.index_types = typ.index_types
        obj.rank = len(obj.index_types)
        obj.types = typ.types
        obj.symmetry = typ.symmetry
        obj.commuting = commuting
        return obj

    def __eq__(self, other):
        if not isinstance(other, TensorHead):
            return False
        return self.name == other.name and self.index_types == other.index_types

    def __ne__(self, other):
        return not (self == other)

    def commutes_with(self, other):
        """
        Returns 0 (1) if self and other (anti)commute.

        Returns None if self and other do not (anti)commute.

        TODO: it should be possible to assign rules for commutations
        between tensors, to be used here.
        """
        if self.commuting == 0 or other.commuting == 0:
            return 0
        if self.commuting == 1 and other.commuting == 1:
            return 1
        return None


    def __str__(self):
        return '%s(%s)' %(self.name, ','.join([str(x) for x in self.index_types]))
    __repr__ = __str__

    def __call__(self, *indices):
        """
        Returns a tensor with indices.

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, TensorType, get_symmetric_group_sgs
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> a, b = tensor_indices('a,b', Lorentz)
        >>> sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
        >>> S2 = TensorType([Lorentz]*2, sym2)
        >>> A = S2('A')
        >>> t = A(a, -b)
        """
        assert self.rank == len(indices)
        components = [self]
        free, dum =  Tensor.from_indices(*indices)
        free.sort(key=lambda x: x[0].name)
        dum.sort()
        return Tensor(S.One, components, free, dum)


class TensExpr(Basic):
    """
    A tensor expression is an expression formed by tensors;
    currently the sums of tensors are distributed.

    A TensExpr can be a TensAdd or a Tensor.

    TensAdd objects are put in canonic form using the Butler-Portugal
    algorithm for canonicalization under monoterm symmetries.

    Tensor objects are formed by products of component tensors,
    and include a coefficient, which is a SymPy expression.


    In the internal representation contracted indices are represented
    by ``(ipos1, ipos2, icomp1, icomp2)``, where ``icomp1`` is the position
    of the component tensor with contravariant index, ``ipos1`` is the
    slot which the index occupies in that component tensor.

    Contracted indices are therefore nameless in the internal representation.
    """

    _op_priority = 11.0
    is_Tensor = True
    is_TensExpr = True
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
        return Tensor.__mul__(self, other)

    def __rmul__(self, other):
        if self.is_TensMul:
            coeff = other*self._coeff
            return Tensor(coeff, self._components, self._free, self._dum)
        return TensAdd(*[x*other for x in self.args])

    def __pow__(self, other):
        raise NotImplementedError

    def __rpow__(self, other):
        raise NotImplementedError

    def __div__(self, other):
        other = sympify(other)
        if other.is_Tensor:
            raise ValueError('cannot divide by a tensor')
        coeff = self._coeff/other
        return Tensor(coeff, self._components, self._free, self._dum, is_canon_bp=self._is_canon_bp)


    def __rdiv__(self, other):
        raise NotImplementedError()

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def substitute_indices(self, *index_tuples):
        """
        Return a tensor with free indices substituted according to ``index_tuples``

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, TensorType
        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, TensorType, get_symmetric_group_sgs
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> i, j, k, l = tensor_indices('i,j,k,l', Lorentz)
        >>> sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
        >>> S2 = TensorType([Lorentz]*2, sym2)
        >>> A, B = S2('A,B')
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
                        if i.is_contravariant == j.is_contravariant:
                            free1.append((v, ipos, cpos))
                        else:
                            free1.append((-v, ipos, cpos))
                        break
                else:
                    free1.append((j, ipos, cpos))

            return Tensor(self._coeff, self._components, free1, self._dum)
        if self.is_TensAdd:
            args = self.args
            args1 = []
            for x in args:
                y = x.substitute_indices(*index_tuples)
                args1.append(y)
            return TensAdd(*args1)

class TensAdd(TensExpr):
    """
    Sum of tensors

    Sum of more than one tensor are put automatically in canonical form.

    Examples
    ========

    >>> from sympy.tensor.tensor import TensorIndexType, TensorSymmetry, TensorType, get_symmetric_group_sgs, tensor_indices, TensAdd
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> a, b = tensor_indices('a,b', Lorentz)
    >>> sym = TensorSymmetry(get_symmetric_group_sgs(1))
    >>> S1 = TensorType([Lorentz], sym)
    >>> p, q = S1('p,q')
    >>> t1 = p(a)
    >>> t2 = q(a)
    >>> t1 + t2
    p(a) + q(a)
    """
    is_Tensor = True
    is_TensAdd = True

    def __new__(cls, *args, **kw_args):
        """
        args  tuple of addends
        """
        args = [sympify(x) for x in args if x]
        if not all(x.is_Tensor for x in args):
            args1 = []
            for x in args:
                if x.is_Tensor:
                    if x.is_TensAdd:
                        args1.extend(list(x.args))
                    else:
                        args1.append(x)
            args1 = [x for x in args1 if x.is_Tensor and x._coeff]
            args2 = [x for x in args if not x.is_Tensor]
            t0 = args1[0]
            if t0.is_TensAdd:
                t0 = t0.args[0]
            if t0._free:
               raise ValueError('all tensors must have the same indices')
            t1 = Tensor(Add(*args2), [], [], [])
            args = [t1] + args1
        a = []
        for x in args:
            if x.is_TensAdd:
                a.extend(list(x.args))
            else:
                a.append(x)
        args = a
        args = [x for x in args if x._coeff]
        if not args:
            return S.Zero

        indices0 = sorted([x[0] for x in args[0]._free], key=lambda x: x.name)
        list_indices = [sorted([y[0] for y in x._free], key=lambda x: x.name) for x in args[1:]]
        if not all(x == indices0 for x in list_indices):
            raise ValueError('all tensors must have the same indices')

        obj = Basic.__new__(cls, **kw_args)
        args.sort(key=lambda x: (x._components, x._free, x._dum))
        a = []
        pprev = None
        prev = args[0]
        prev_coeff = prev._coeff
        changed = False
        new = 0
        if len(args) == 1 and args[0].is_TensMul:
            obj._args = tuple(args)
            return obj

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
                        t = Tensor(prev_coeff, prev._components,
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
            prev = Tensor(prev_coeff, prev._components, prev._free, prev._dum)
            a.append(prev)
        elif op == 1:
            a.append(prev)
        if not a:
            return S.Zero

        # TODO introduce option not to use canon_bp automatically in TensAdd
        if all(x.is_TensMul for x in a):
            a = [canon_bp(x) for x in a]
        a = [x for x in a if x]
        if not a:
            return S.Zero
        obj._args = tuple(a)
        return obj

    def canon_bp(self):
        """
        canonicalize using the Butler-Portugal algorithm for canonicalization
        under monoterm symmetries.
        """
        args = [x.canon_bp() for x in self.args]
        args.sort(key=lambda x: (x._components, x._free, x._dum))
        res = TensAdd(*args)
        return res

    def __eq__(self, other):
        other = sympify(other)
        if not other.is_Tensor:
            if len(self.args) == 1:
                return self.args[0]._coeff == other
        if other.is_Tensor and other.is_TensMul and other._coeff == 0:
            return self == 0
        t = self - other
        if not t.is_Tensor:
            return t == 0
        else:
            if t.is_TensMul:
                return t._coeff == 0
            else:
                return all(x._coeff == 0 for x in t.args)

    def __ne__(self, other):
        return not (self == other)

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
        return t.canon_bp()

    def _pretty(self):
        a = []
        args = self.args
        for x in args:
            a.append(str(x))
        a.sort()
        s = ' + '.join(a)
        s = s.replace('+ -', '- ')
        return s

class Tensor(TensExpr):
    """
    Product of tensors
    """
    is_Tensor = True
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

    def __eq__(self, other):
        if other == 0:
            return self._coeff == 0
        other = sympify(other)
        if not other.is_Tensor :
            return False
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

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, TensorType, Tensor, get_symmetric_group_sgs
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2, m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
        >>> Tensor.from_indices(m0, m1, -m1, m3)
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
            contr = index.is_contravariant
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
                index_dict[(name, typ)] = index.is_contravariant, i

        free_indices = [(index, i, 0) for i, index in enumerate(indices) if free[i]]
        free = sorted(free_indices, key=lambda x: (x[0].tensortype, x[0].name))

        return free, dum

    def get_indices(self):
        """
        Returns the list of indices of the tensor

        The indices are listed in the order in which they appear in the
        component tensors.

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, get_symmetric_group_sgs, TensorType
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2 = tensor_indices('m0,m1,m2', Lorentz)
        >>> sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
        >>> S1 = TensorType([Lorentz], sym1)
        >>> g = Lorentz.metric
        >>> p, q = S1('p,q')
        >>> t = p(m1)*g(m0,m2)
        >>> t.get_indices()
        [m1, m0, m2]
        """
        indices = [None]*self.rank
        start = 0
        pos = 0
        vpos = []
        components = self._components
        for t in self._components:
            vpos.append(pos)
            pos += t.rank
        cdt = defaultdict(int)
        for indx, ipos, cpos in self._free:
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

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, TensorType, get_symmetric_group_sgs
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> a, b, c, d = tensor_indices('a,b,c,d', Lorentz)
        >>> sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
        >>> S2 = TensorType([Lorentz]*2, sym2)
        >>> A, B = S2('A,B')
        >>> t = A(a,b)*B(-b,c)
        >>> t
        A(a, L_0)*B(-L_0, c)
        >>> t.split()
        [A(a, L_0), B(-L_0, c)]
        """
        indices = self.get_indices()
        pos = 0
        components = self._components
        res = []
        for t in self._components:
            t1 = t(*indices[pos:pos + t.rank])
            pos += t.rank
            res.append(t1)
        res[0] = Tensor(self._coeff, res[0]._components, res[0]._free, res[0]._dum, is_canon_bp=res[0]._is_canon_bp)
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
        for t in self._components:
            vpos.append(pos)
            pos += t.rank
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
                msym.append(typ.metric_sym)
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
            v.append((h.symmetry.base, h.symmetry.generators, n, h.commuting))
        return _af_new(g), dummies, msym, v

    def __mul__(self, other):
        """
        Multiply two tensors using Einstein summation convention.

        If the two tensors have an index in common, one contravariant
        and the other covariant, in their product the indices are summed

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, get_symmetric_group_sgs, TensorType
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2 = tensor_indices('m0,m1,m2', Lorentz)
        >>> sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
        >>> S1 = TensorType([Lorentz], sym1)
        >>> g = Lorentz.metric
        >>> p, q = S1('p,q')
        >>> t1 = p(m0)
        >>> t2 = q(-m0)
        >>> t1*t2
        p(L_0)*q(-L_0)
        """
        other = sympify(other)
        if not other.is_Tensor:
            coeff = self._coeff*other
            return Tensor(coeff, self._components, self._free, self._dum, is_canon_bp=self._is_canon_bp)
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
            if ind1.is_contravariant == ind2.is_contravariant:
                raise ValueError('wrong index contruction %s' % ind1)
            if ind1.is_contravariant:
                new_dummy = (ipos1, ipos2, cpos1, cpos2)
            else:
                new_dummy = (ipos2, ipos1, cpos2, cpos1)
            dum.append(new_dummy)
        coeff = self._coeff*other._coeff
        return Tensor(coeff, components, free, dum)


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
                if c == None:
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
        dum = [(i1, i2, perm[c1], perm[c2]) for i1, i2, c1, c2 in self._dum]
        free.sort(key = lambda x: (x[0].tensortype, x[0].name))
        dum.sort(key = lambda x: components[x[2]].index_types[x[0]])
        if sign == -1:
            coeff = -self._coeff
        else:
            coeff = self._coeff
        t = Tensor(coeff, components, free, dum)
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
        for t in self._components:
            vpos.append(pos)
            pos += t.rank
        free_indices = [x[0] for x in self._free]
        sorted_free = sorted(free_indices, key=lambda x:(x.tensortype, x.name))
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
        res = Tensor(coeff, components, free, dum, is_canon_bp=canon_bp)
        return res

    def canon_bp(self):
        """
        canonicalize using the Butler-Portugal algorithm for canonicalization
        under monoterm symmetries.

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, get_symmetric_group_sgs, TensorType
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2 = tensor_indices('m0,m1,m2', Lorentz)
        >>> sym2a = TensorSymmetry(get_symmetric_group_sgs(2, 1))
        >>> S2 = TensorType([Lorentz]*2, sym2a)
        >>> A = S2('A')
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


    def contract_metric(self, g, contract_all=False):
        """
        Raise or lower indices with the metric ``g``

        ``g``  metric

        ``contract_all`` if True, eliminate all ``g`` which are contracted

        Examples
        ========

        >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, get_symmetric_group_sgs, TensorType
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2 = tensor_indices('m0,m1,m2', Lorentz)
        >>> sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
        >>> S1 = TensorType([Lorentz], sym1)
        >>> g = Lorentz.metric
        >>> p, q = S1('p,q')
        >>> t = p(m0)*q(m1)*g(-m0, -m1)
        >>> t.canon_bp()
        metric(L_0, L_1)*p(-L_0)*q(-L_1)
        >>> t.contract_metric(g).canon_bp()
        p(L_0)*q(-L_0)
        """
        if g.index_types[0].metric_sym != 0:
            # TODO case of antisymmetric metric
            raise NotImplementedError
        if not self._components:
            return self
        free_indices = [x[0] for x in self._free]
        a = self.split()
        typ = g.index_types[0]
        # if a component tensor of a has 2 dummy indices, it is g(d,-d) = dim
        for i, tx in enumerate(a):
            if tx._components[0] == g:
                free_indices_g = [x[0] for x in a[i]._free]
                if len(free_indices_g) == 0:
                    a1 = a[:i] + a[i + 1:]
                    t = tensor_mul(*a1)*(typ.dim*a[i]._coeff)
                    if contract_all == True and g in t._components:
                        return t.contract_metric(g, True)
                    return t

        # if all metric tensors have only free indices, there is no contraction
        for i, tg in enumerate(a):
            if tg._components[0] == g:
                tg_free_indices = [x[0] for x in tg._free]
                if all(indx in free_indices for indx in tg_free_indices):
                    continue
                break
        else:
            return self

        # tg has one or two indices contracted with other tensors
        # i position of tg in a
        coeff = S.One
        tg_free = tg._free
        if tg_free[0][0] in free_indices or tg_free[1][0] in free_indices:
            if tg_free[0][0] in free_indices:
                ind_free = tg_free[0][0]
                ind, ipos1, _ = tg_free[1]
            else:
                ind_free = tg_free[1][0]
                ind, ipos1, _ = tg_free[0]

            ind1 = -ind
            # search ind1 in self._dum
            for j, tx in enumerate(a):
                if ind1 in [x[0] for x in tx._free]:
                    break
            free1 = []
            for indx, iposx, _ in tx._free:
                if indx == ind1:
                    free1.append((ind_free, iposx, 0))
                else:
                    free1.append((indx, iposx, 0))
            t1 = Tensor(tx._coeff, tx._components, free1, tx._dum)
            a[j] = t1
            a = a[:i] + a[i + 1:]
            coeff = coeff*tg._coeff
            res = tensor_mul(*a)
        else:
            # tg has two indices contracted with other tensors
            ind1 = tg_free[0][0]
            ind2 = tg_free[1][0]
            ind1m = -ind1
            ind2m = -ind2
            for k, ty in enumerate(a):
                if ind2m in [x[0] for x in ty._free]:
                    break
            if ty._components == [g]:
                ty_indices = [x[0] for  x in ty._free]
                if all(x in [ind1m, ind2m] for x in ty_indices):
                    if i < k:
                        a = a[:i] + a[i+1:k] + a[k+1:]
                    else:
                        a = a[:k] + a[k+1:i] + a[k+1:]
                    if a:
                        res = tensor_mul(*a)
                        res = (coeff*typ.dim*tg._coeff*ty._coeff)*res
                    else:
                        res = coeff*typ.dim*tg._coeff*ty._coeff
                        res = Tensor(res, [],[],[], is_canon_bp=True)
                    if contract_all == True and g in res._components:
                        return res.contract_metric(g, True)
                    return res

            free2 = []
            for indx, iposx, _ in ty._free:
                if indx == ind2m:
                    free2.append((ind1, iposx, 0))
                else:
                    free2.append((indx, iposx, 0))
            t2 = Tensor(ty._coeff, ty._components, free2, ty._dum)
            a[k] = t2
            a = a[:i] + a[i + 1:]
            coeff = coeff*tg._coeff
            res = tensor_mul(*a)
        res = coeff*res
        if contract_all == True and g in res._components:
            return res.contract_metric(g, True)
        return res


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
    if p.is_Tensor:
        return p.canon_bp()
    return p

def tensor_mul(*a):
    """
    product of tensors
    """
    if not a:
        return Tensor(S.One, [], [], [])
    t = a[0]
    for tx in a[1:]:
        t = t*tx
    return t

def tensorlist_contract_metric(a, tg):
    """
    contract `tg` with a tensor in the list `a = t.split()`
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
                t2 = Tensor(t1._coeff, t1._components, free1, t1._dum)
                a[i] = t2
                return a
    a.append(tg)
    return a


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

    >>> from sympy.tensor.tensor import TensorIndexType, tensor_indices, TensorSymmetry, TensorType, riemann_cyclic, riemann_bsgs
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> i, j, k, l = tensor_indices('i,j,k,l', Lorentz)
    >>> symr = TensorSymmetry(riemann_bsgs)
    >>> R4 = TensorType([Lorentz]*4, symr)
    >>> R = R4('R')
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
