"""
TODO: legacy class to store deprecated classes.
"""
from sympy.core.sympify import CantSympify
from sympy.core.decorators import deprecated


class TIDS(CantSympify):
    """
    Tensor internal data structure. This contains internal data structures about
    components of a tensor expression, its free and dummy indices.

    To create a ``TIDS`` object via the standard constructor, the required
    arguments are

    WARNING: this class is meant as an internal representation of tensor data
    structures and should not be directly accessed by end users.

    Parameters
    ==========

    components : ``TensorHead`` objects representing the components of the tensor expression.

    free : Free indices in their internal representation.

    dum : Dummy indices in their internal representation.

    Examples
    ========

    >>> from sympy.tensor.poly_tensor import TensorIndexType, tensor_indices, TIDS, tensorhead
    >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    >>> m0, m1, m2, m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
    >>> T = tensorhead('T', [Lorentz]*4, [[1]*4])
    >>> TIDS([T], [(m0, 0, 0), (m3, 3, 0)], [(1, 2, 0, 0)])
    TIDS([T(Lorentz,Lorentz,Lorentz,Lorentz)], [(m0, 0, 0), (m3, 3, 0)], [(1, 2, 0, 0)])

    Notes
    =====

    In short, this has created the components, free and dummy indices for
    the internal representation of a tensor T(m0, m1, -m1, m3).

    Free indices are represented as a list of triplets. The elements of
    each triplet identify a single free index and are

    1. TensorIndex object
    2. position inside the component
    3. component number

    Dummy indices are represented as a list of 4-plets. Each 4-plet stands
    for couple for contracted indices, their original TensorIndex is not
    stored as it is no longer required. The four elements of the 4-plet
    are

    1. position inside the component of the first index.
    2. position inside the component of the second index.
    3. component number of the first index.
    4. component number of the second index.

    """

    def __init__(self, components, free, dum):
        self.components = components
        self.free = free
        self.dum = dum
        self._ext_rank = len(self.free) + 2*len(self.dum)

    def get_components_with_free_indices(self):
        """
        Get a list of components with their associated indices.

        Examples
        ========

        >>> from sympy.tensor.poly_tensor import TensorIndexType, tensor_indices, TIDS, tensorhead
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2, m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
        >>> T = tensorhead('T', [Lorentz]*4, [[1]*4])
        >>> A = tensorhead('A', [Lorentz], [[1]])
        >>> t = TIDS.from_components_and_indices([T], [m0, m1, -m1, m3])
        >>> t.get_components_with_free_indices()
        [(T(Lorentz,Lorentz,Lorentz,Lorentz), [(m0, 0, 0), (m3, 3, 0)])]
        >>> t2 = (A(m0)*A(-m0))._tids
        >>> t2.get_components_with_free_indices()
        [(A(Lorentz), []), (A(Lorentz), [])]
        >>> t3 = (A(m0)*A(-m1)*A(-m0)*A(m1))._tids
        >>> t3.get_components_with_free_indices()
        [(A(Lorentz), []), (A(Lorentz), []), (A(Lorentz), []), (A(Lorentz), [])]
        >>> t4 = (A(m0)*A(m1)*A(-m0))._tids
        >>> t4.get_components_with_free_indices()
        [(A(Lorentz), []), (A(Lorentz), [(m1, 0, 1)]), (A(Lorentz), [])]
        >>> t5 = (A(m0)*A(m1)*A(m2))._tids
        >>> t5.get_components_with_free_indices()
        [(A(Lorentz), [(m0, 0, 0)]), (A(Lorentz), [(m1, 0, 1)]), (A(Lorentz), [(m2, 0, 2)])]
        """
        components = self.components
        ret_comp = []

        free_counter = 0
        if len(self.free) == 0:
            return [(comp, []) for comp in components]

        for i, comp in enumerate(components):
            c_free = []
            while free_counter < len(self.free):
                if not self.free[free_counter][2] == i:
                    break

                c_free.append(self.free[free_counter])
                free_counter += 1

                if free_counter >= len(self.free):
                    break
            ret_comp.append((comp, c_free))

        return ret_comp

    @staticmethod
    def from_components_and_indices(components, indices):
        """
        Create a new ``TIDS`` object from ``components`` and ``indices``

        ``components``  ``TensorHead`` objects representing the components
                        of the tensor expression.

        ``indices``     ``TensorIndex`` objects, the indices. Contractions are
                        detected upon construction.

        Examples
        ========

        >>> from sympy.tensor.poly_tensor import TensorIndexType, tensor_indices, TIDS, tensorhead
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2, m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
        >>> T = tensorhead('T', [Lorentz]*4, [[1]*4])
        >>> TIDS.from_components_and_indices([T], [m0, m1, -m1, m3])
        TIDS([T(Lorentz,Lorentz,Lorentz,Lorentz)], [(m0, 0, 0), (m3, 3, 0)], [(1, 2, 0, 0)])

        In case of many components the same indices have slightly different
        indexes:

        >>> A = tensorhead('A', [Lorentz], [[1]])
        >>> TIDS.from_components_and_indices([A]*4, [m0, m1, -m1, m3])
        TIDS([A(Lorentz), A(Lorentz), A(Lorentz), A(Lorentz)], [(m0, 0, 0), (m3, 0, 3)], [(0, 0, 1, 2)])
        """
        tids = None
        cur_pos = 0
        for i in components:
            tids_sing = TIDS([i], *TIDS.free_dum_from_indices(*indices[cur_pos:cur_pos+i.rank]))
            if tids is None:
                tids = tids_sing
            else:
                tids *= tids_sing
            cur_pos += i.rank

        if tids is None:
            tids = TIDS([], [], [])

        tids.free.sort(key=lambda x: x[0].name)
        tids.dum.sort()

        return tids

    def to_indices(self):
        """
        Get a list of indices, creating new tensor indices to complete dummy indices.
        """
        component_indices = []
        for i in self.components:
            component_indices.append([None]*i.rank)

        for i in self.free:
            component_indices[i[2]][i[1]] = i[0]

        from sympy.tensor.poly_tensor import TensorIndex

        for i, dummy_pos in enumerate(self.dum):
            # TODO: add test for special case (Flavor index)
            tensor_index_type = self.components[dummy_pos[2]].args[1].args[0][dummy_pos[0]]
            dummy_index = TensorIndex('dummy_index_{0}'.format(i), tensor_index_type)
            component_indices[dummy_pos[2]][dummy_pos[0]] = dummy_index
            component_indices[dummy_pos[3]][dummy_pos[1]] = -dummy_index

        indices = []
        for i in component_indices:
            indices.extend(i)

        return indices

    @staticmethod
    def free_dum_from_indices(*indices):
        """
        Convert ``indices`` into ``free``, ``dum`` for single component tensor

        ``free``     list of tuples ``(index, pos, 0)``,
                     where ``pos`` is the position of index in
                     the list of indices formed by the component tensors

        ``dum``      list of tuples ``(pos_contr, pos_cov, 0, 0)``

        Examples
        ========

        >>> from sympy.tensor.poly_tensor import TensorIndexType, tensor_indices, TIDS
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2, m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
        >>> TIDS.free_dum_from_indices(m0, m1, -m1, m3)
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
            name = index._name
            typ = index._tensortype
            contr = index._is_up
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
                index_dict[(name, typ)] = index._is_up, i

        free = [(index, i, 0) for i, index in enumerate(indices) if free[i]]
        free.sort()
        return free, dum

    @staticmethod
    def mul(f, g):
        """
        The algorithms performing the multiplication of two ``TIDS`` instances.

        In short, it forms a new ``TIDS`` object, joining components and indices,
        checking that abstract indices are compatible, and possibly contracting
        them.

        Examples
        ========

        >>> from sympy.tensor.poly_tensor import TensorIndexType, tensor_indices, TIDS, tensorhead
        >>> Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
        >>> m0, m1, m2, m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
        >>> T = tensorhead('T', [Lorentz]*4, [[1]*4])
        >>> A = tensorhead('A', [Lorentz], [[1]])
        >>> tids_1 = TIDS.from_components_and_indices([T], [m0, m1, -m1, m3])
        >>> tids_2 = TIDS.from_components_and_indices([A], [m2])
        >>> tids_1 * tids_2
        TIDS([T(Lorentz,Lorentz,Lorentz,Lorentz), A(Lorentz)],\
            [(m0, 0, 0), (m3, 3, 0), (m2, 0, 1)], [(1, 2, 0, 0)])

        In this case no contraction has been performed.

        >>> tids_3 = TIDS.from_components_and_indices([A], [-m3])
        >>> tids_1 * tids_3
        TIDS([T(Lorentz,Lorentz,Lorentz,Lorentz), A(Lorentz)],\
            [(m0, 0, 0)], [(1, 2, 0, 0), (3, 0, 0, 1)])

        Free indices ``m3`` and ``-m3`` are identified as a contracted couple, and are
        therefore transformed into dummy indices.

        A wrong index construction (for example, trying to contract two
        contravariant indices or using indices multiple times) would result in
        an exception:

        >>> tids_4 = TIDS.from_components_and_indices([A], [m3])
        >>> # This raises an exception:
        >>> # tids_1 * tids_4
        """
        index_up = lambda u: u if u.is_up else -u

        # find out which free indices of f and g are contracted
        free_dict1 = dict([(i if i.is_up else -i, (pos, cpos, i)) for i, pos, cpos in f.free])
        free_dict2 = dict([(i if i.is_up else -i, (pos, cpos, i)) for i, pos, cpos in g.free])

        free_names = set(free_dict1.keys()) & set(free_dict2.keys())
        # find the new `free` and `dum`
        nc1 = len(f.components)
        dum2 = [(i1, i2, c1 + nc1, c2 + nc1) for i1, i2, c1, c2 in g.dum]
        free1 = [(ind, i, c) for ind, i, c in f.free if index_up(ind) not in free_names]
        free2 = [(ind, i, c + nc1) for ind, i, c in g.free if index_up(ind) not in free_names]
        free = free1 + free2
        dum = f.dum + dum2
        for name in free_names:
            ipos1, cpos1, ind1 = free_dict1[name]
            ipos2, cpos2, ind2 = free_dict2[name]
            cpos2 += nc1
            if ind1._is_up == ind2._is_up:
                raise ValueError('wrong index construction {0}'.format(ind1))
            if ind1._is_up:
                new_dummy = (ipos1, ipos2, cpos1, cpos2)
            else:
                new_dummy = (ipos2, ipos1, cpos2, cpos1)
            dum.append(new_dummy)
        return (f.components + g.components, free, dum)

    def __mul__(self, other):
        return TIDS(*self.mul(self, other))

    def __str__(self):
        return "TIDS({0}, {1}, {2})".format(self.components, self.free, self.dum)

    def __repr__(self):
        return self.__str__()

    def to_tensmul_args(self):
        """
        TODO

        EXPERIMENTAL

        Warning: anomalous method because ...
        """
        if len(self.components) == 0:
            return (S.One,)
        indices = self.to_indices()
        tensmul_args = []

        for component in self.components:
            comp_ind_size = len(component.index_types)
            ti = indices[:comp_ind_size]
            indices = indices[comp_ind_size:]
            tensmul_args.append(component(*ti))

        return tensmul_args

    def sorted_components(self):
        """
        Returns a ``TIDS`` with sorted components

        The sorting is done taking into account the commutation group
        of the component tensors.
        """
        from sympy.combinatorics.permutations import _af_invert
        cv = list(zip(self.components, range(len(self.components))))
        sign = 1
        n = len(cv) - 1
        for i in range(n):
            for j in range(n, i, -1):
                c = cv[j-1][0].commutes_with(cv[j][0])
                if c not in [0, 1]:
                    continue
                if (cv[j-1][0]._types, cv[j-1][0]._name) > \
                        (cv[j][0]._types, cv[j][0]._name):
                    cv[j-1], cv[j] = cv[j], cv[j-1]
                    if c:
                        sign = -sign

        # perm_inv[new_pos] = old_pos
        components = [x[0] for x in cv]
        perm_inv = [x[1] for x in cv]
        perm = _af_invert(perm_inv)
        free = [(ind, i, perm[c]) for ind, i, c in self.free]
        free.sort()
        dum = [(i1, i2, perm[c1], perm[c2]) for i1, i2, c1, c2 in self.dum]
        dum.sort(key=lambda x: components[x[2]].index_types[x[0]])

        return TIDS(components, free, dum), sign

    def canon_args(self):
        """
        Returns ``(g, dummies, msym, v)``, the entries of ``canonicalize``

        see ``canonicalize`` in ``tensor_can.py``
        """
        # to be called after sorted_components
        from sympy.combinatorics.permutations import _af_new
#         types = list(set(self._types))
#         types.sort(key = lambda x: x._name)
        n = self._ext_rank
        g = [None]*n + [n, n+1]
        pos = 0
        vpos = []
        components = self.components
        for t in components:
            vpos.append(pos)
            pos += t.rank
        # ordered indices: first the free indices, ordered by types
        # then the dummy indices, ordered by types and contravariant before
        # covariant
        # g[position in tensor] = position in ordered indices
        for i, (indx, ipos, cpos) in enumerate(self.free):
            pos = vpos[cpos] + ipos
            g[pos] = i
        pos = len(self.free)
        j = len(self.free)
        dummies = []
        prev = None
        a = []
        msym = []
        for ipos1, ipos2, cpos1, cpos2 in self.dum:
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
        from sympy.tensor.poly_tensor import TensorManager
        for h, n in numtyp:
            if h._comm == 0 or h._comm == 1:
                comm = h._comm
            else:
                comm = TensorManager.get_comm(h._comm, h._comm)
            v.append((h._symmetry.base, h._symmetry.generators, n, comm))
        return _af_new(g), dummies, msym, v

    def perm2tensor(self, g, canon_bp=False):
        """
        Returns a ``TIDS`` instance corresponding to the permutation ``g``

        ``g``  permutation corresponding to the tensor in the representation
        used in canonicalization

        ``canon_bp``   if True, then ``g`` is the permutation
        corresponding to the canonical form of the tensor
        """
        vpos = []
        components = self.components
        pos = 0
        for t in components:
            vpos.append(pos)
            pos += t.rank
        sorted_free = [x[0] for x in self.free]
        sorted_free.sort()
        nfree = len(sorted_free)
        rank = self._ext_rank
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

        return TIDS(components, free, dum)


class VTIDS(TIDS):
    """
    DEPRECATED: DO NOT USE.
    """

    @deprecated(useinstead="TIDS")
    def __init__(self, components, free, dum, data):
        super(VTIDS, self).__init__(components, free, dum)
        self.data = data

    @staticmethod
    @deprecated(useinstead="TIDS")
    def parse_data(data):
        """
        DEPRECATED: DO NOT USE.
        """
        return _TensorDataLazyEvaluator.parse_data(data)

    @deprecated(useinstead="TIDS")
    def correct_signature_from_indices(self, data, indices, free, dum):
        """
        DEPRECATED: DO NOT USE.
        """
        return _TensorDataLazyEvaluator._correct_signature_from_indices(data, indices, free, dum)

    @staticmethod
    @deprecated(useinstead="TIDS")
    def flip_index_by_metric(data, metric, pos):
        """
        DEPRECATED: DO NOT USE.
        """
        return _TensorDataLazyEvaluator._flip_index_by_metric(data, metric, pos)
