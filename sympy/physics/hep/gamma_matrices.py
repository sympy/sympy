from sympy import S
from sympy.tensor.tensor import TensorIndexType, TensorIndex,\
    TensMul, TensorHead, tensorsymmetry, TensorType,\
    TensAdd, tensor_mul, get_lines
from sympy.core.containers import Tuple


DiracSpinorIndex = TensorIndexType('DiracSpinorIndex', dim=4, dummy_fmt="S")


class _LorentzContainer(object):
    """
    Helper to collect LorentzIndex indices in various dimensions.

    It collects LorentzIndex TensorIndexType that have been implemented in the code,
    and stores them in a dict()
    """
    lorentz_types = dict()

    def __new__(cls, dim=4, eps_dim=None, dummy_fmt="L"):
        if (dim, eps_dim) in _LorentzContainer.lorentz_types:
            return _LorentzContainer.lorentz_types[(dim, eps_dim)]

        new_L = TensorIndexType("LorentzIndex", dim=dim, eps_dim=eps_dim, dummy_fmt=dummy_fmt)
        _LorentzContainer.lorentz_types[(dim, eps_dim)] = new_L
        return new_L


class GammaMatrixHead(TensorHead):
    r"""
    Class to wrap a ``TensorHead`` for gamma matrices.

    ``dim``       dimension of the gamma matrix.
    ``eps_dim``   correction for dimensional regularization, use None if not needed.

    Examples
    ========

    >>> from sympy.physics.hep.gamma_matrices import GammaMatrixHead
    >>> from sympy.tensor.tensor import tensor_indices
    >>> G = GammaMatrixHead()
    >>> i = tensor_indices('i', G.LorentzIndex)
    >>> G(i)
    gamma(i, auto_left, -auto_right)

    Note that there is already an instance of GammaMatrixHead in four dimensions:
    GammaMatrix, which is simply declare as

    ``GammaMatrix = GammaMatrixHead()``

    >>> from sympy.physics.hep.gamma_matrices import GammaMatrix
    >>> from sympy.tensor.tensor import tensor_indices
    >>> i = tensor_indices('i', GammaMatrix.LorentzIndex)
    >>> GammaMatrix(i)
    gamma(i, auto_left, -auto_right)

    To access the metric tensor

    >>> GammaMatrix.LorentzIndex.metric
    metric(LorentzIndex,LorentzIndex)

    """
    _gmhd = dict()

    def __new__(cls, dim=4, eps_dim=4):
        key = (dim, eps_dim)
        if key in GammaMatrixHead._gmhd:
            return GammaMatrixHead._gmhd[key]

        lorentz = _LorentzContainer(*key)

        gmh = TensorHead.__new__(cls, "gamma", TensorType(Tuple(lorentz, DiracSpinorIndex, DiracSpinorIndex), tensorsymmetry([1], [1], [1])), comm=2, matrix_behavior=True)
        GammaMatrixHead._gmhd[key] = gmh
        gmh.LorentzIndex = lorentz
        return gmh

    @staticmethod
    def extract_type_tens(expression):
        """
        Extract from a ``TensExpr`` all elements of this type.

        Returns two tensor expressions:

        * the first contains all ``TensorHead`` of this type.
        * the second contains all remaining.


        """
        sp = expression.split()

        # Collect all gamma matrices of the same dimension
        new_expr = S.One
        residual_expr = S.One
        for i in sp:
            if isinstance(i.args[1][0], GammaMatrixHead):
                new_expr *= i
            else:
                residual_expr *= i
        return new_expr, residual_expr

    @staticmethod
    def simplify_this_type(expression):
        extracted_expr, residual_expr = GammaMatrixHead.extract_type_tens(expression)
        res_expr = GammaMatrixHead._simplify_single_line(extracted_expr)
        return res_expr * residual_expr

    @staticmethod
    def simplify_gpgp(ex, sort=True):
        """
        simplify products ``G(i)*p(-i)*G(j)*p(-j) -> p(i)*p(-i)``

        Examples
        ========

        >>> from sympy.physics.hep.gamma_matrices import GammaMatrix as G
        >>> from sympy.tensor.tensor import tensor_indices, tensorhead
        >>> p, q = tensorhead('p, q', [G.LorentzIndex], [[1]])
        >>> i0,i1,i2,i3,i4,i5 = tensor_indices('i0:6', G.LorentzIndex)
        >>> ps = p(i0)*G(-i0)
        >>> qs = q(i0)*G(-i0)
        >>> G.simplify_gpgp(ps*qs*qs)
        gamma(-L_0, auto_left, -auto_right)*p(L_0)*q(L_1)*q(-L_1)
        """
        def _simplify_gpgp(ex):
            tids = ex._tids
            components = tids.components
            a = []
            for i in range(len(components)):
                if not isinstance(components[i], GammaMatrixHead):
                    continue
                dum = tids.dum
                for dx in dum:
                    if dx[2] == i:
                        p_pos1 = dx[3]
                    elif dx[3] == i:
                        p_pos1 = dx[2]
                    else:
                        continue
                    comp1 = components[p_pos1]
                    if comp1.comm == 0 and comp1.rank == 1:
                        a.append((i, p_pos1))
            if not a:
                return ex
            elim = set()
            tv = []
            hit = True
            coeff = S.One
            ta = None
            while hit:
                hit = False
                for i, ai in enumerate(a[:-1]):
                    if ai[0] in elim:
                        continue
                    if ai[0] != a[i + 1][0] - 1:
                        continue
                    if components[ai[1]] != components[a[i + 1][1]]:
                        continue
                    elim.add(ai[0])
                    elim.add(ai[1])
                    elim.add(a[i + 1][0])
                    elim.add(a[i + 1][1])
                    if not ta:
                        ta = ex.split()
                        mu = TensorIndex('mu', GammaMatrix.LorentzIndex)
                    ind1 = ta[ai[0]].args[-1][1]
                    ind2 = ta[ai[0] + 1].args[-1][2]
                    hit = True
                    if i == 0:
                        coeff = ex.coeff
                    tx = components[ai[1]](mu)*components[ai[1]](-mu)
                    tv.append(tx*DiracSpinorIndex.delta(ind1, ind2))
                    break

            if tv:
                a = [x for j, x in enumerate(ta) if j not in elim]
                a.extend(tv)
                t = tensor_mul(*a)*coeff
                t = t.contract_metric(DiracSpinorIndex.delta)
                return t
            else:
                return ex

        if sort:
            ex = ex.sorted_components()
        while 1:
            t = _simplify_gpgp(ex)
            if t != ex:
                ex = t
            else:
                return t

    @staticmethod
    def simplify_lines(ex):
        """
        simplify a product of gamma matrices

        Examples
        ========

        >>> from sympy.physics.hep.gamma_matrices import GammaMatrix, DiracSpinorIndex
        >>> from sympy.tensor.tensor import tensor_indices
        >>> i0,i1,i2,i3,i4,i5 = tensor_indices('i0:6', GammaMatrix.LorentzIndex)
        >>> s0,s1,s2,s3,s4,s5,s6,s7 = tensor_indices('s0:8', DiracSpinorIndex)
        >>> G = GammaMatrix
        >>> t = G(i1,s1,-s2)*G(i4,s7,-s6)*G(i2,s2,-s3)*G(i3,s4,-s5)*G(i5,s6,-s7)
        >>> G.simplify_lines(t)
        4*gamma(i3, s4, -s5)*gamma(i1, s1, -S_0)*gamma(i2, S_0, -s3)*metric(i4, i5)

        """
        lines, traces, rest = get_lines(ex, DiracSpinorIndex)
        a = ex.split()
        trest = tensor_mul(*[x for i, x in enumerate(a) if i in rest])
        tlines = []
        for line in lines:
            first = a[line[0]]
            last = a[line[-1]]
            first = [x[0] for x in first.free if x[1] == 1][0]
            last = [x[0] for x in last.free if x[1] == 2][0]
            tx = tensor_mul(*[x for i, x in enumerate(a) if i  in line])
            tx1 = GammaMatrixHead._simplify_single_line(tx)
            tlines.append(tx1)
        traces = [GammaMatrix._trace_single_line(tensor_mul(*[x for i, x in enumerate(a) if i  in line])) for line in traces]
        res = tensor_mul(*([trest] + tlines + traces))
        return res

    def gamma_trace(self, t):
        """
        trace of a single line of gamma matrices

        Examples
        ========

        >>> from sympy.physics.hep.gamma_matrices import GammaMatrix as G
        >>> from sympy.tensor.tensor import tensor_indices, tensorhead
        >>> p, q = tensorhead('p, q', [G.LorentzIndex], [[1]])
        >>> i0,i1,i2,i3,i4,i5 = tensor_indices('i0:6', G.LorentzIndex)
        >>> ps = p(i0)*G(-i0)
        >>> qs = q(i0)*G(-i0)
        >>> G.gamma_trace(G(i0)*G(i1))
        4*metric(i0, i1)
        >>> G.gamma_trace(ps*ps) - 4*p(i0)*p(-i0)
        0
        >>> G.gamma_trace(ps*qs + ps*ps) - 4*p(i0)*p(-i0) - 4*p(i0)*q(-i0)
        0

        """
        #assert any(x == DiracSpinorIndex.auto_right for x, p, c, in t._tids.free)
        if isinstance(t, TensAdd):
            res = TensAdd(*[self._trace_single_line(x) for x in t.args])
            return res
        t = self._simplify_single_line(t)
        res = self._trace_single_line(t)
        return res

    @staticmethod
    def _simplify_single_line(expression):
        """
        Simplify single-line product of gamma matrices.

        Examples
        ========

        >>> from sympy.physics.hep.gamma_matrices import GammaMatrix as G, DiracSpinorIndex as DS
        >>> from sympy.tensor.tensor import tensor_indices, tensorhead
        >>> p = tensorhead('p', [G.LorentzIndex], [[1]])
        >>> i0,i1 = tensor_indices('i0:2', G.LorentzIndex)
        >>> G._simplify_single_line(G(i0)*G(i1)*p(-i1)*G(-i0)) + 2*G(i0)*p(-i0)
        0

        """
        t1, t2 = GammaMatrixHead.extract_type_tens(expression)
        if t1 != 1:
            t1 = GammaMatrixHead._kahane_simplify(t1.coeff, t1._tids)
        res = t1*t2
        return res

    def _trace_single_line(self, t):
        """
        Evaluate the trace of a single gamma matrix line inside a ``TensExpr``.

        Notes
        =====

        If there are ``DiracSpinorIndex.auto_left`` and ``DiracSpinorIndex.auto_right``
        indices trace over them; otherwise traces are not implied (explain)


        Examples
        ========

        >>> from sympy.physics.hep.gamma_matrices import GammaMatrix as G
        >>> from sympy.tensor.tensor import tensor_indices, tensorhead
        >>> p = tensorhead('p', [G.LorentzIndex], [[1]])
        >>> i0,i1,i2,i3,i4,i5 = tensor_indices('i0:6', G.LorentzIndex)
        >>> G._trace_single_line(G(i0)*G(i1))
        4*metric(i0, i1)
        >>> G._trace_single_line(G(i0)*p(-i0)*G(i1)*p(-i1)) - 4*p(i0)*p(-i0)
        0

        """
        def _trace_single_line1(t):
            t = t.sorted_components()
            components = t.components
            ncomps = len(components)
            g = self.LorentzIndex.metric
            sg = DiracSpinorIndex.delta
            # gamma matirices are in a[i:j]
            hit = 0
            for i in range(ncomps):
                if isinstance(components[i], GammaMatrixHead):
                    hit = 1
                    break

            for j in range(i + hit, ncomps):
                if not isinstance(components[j], GammaMatrixHead):
                    break
            else:
                j = ncomps
            numG = j - i
            if numG == 0:
                spinor_free = [_[0] for _ in t._tids.free if _[0].tensortype is DiracSpinorIndex]
                tcoeff = t.coeff
                if spinor_free == [DiracSpinorIndex.auto_left, -DiracSpinorIndex.auto_right]:
                    t = t*DiracSpinorIndex.delta(-DiracSpinorIndex.auto_left, DiracSpinorIndex.auto_right)
                    t = t.contract_metric(sg)
                    return t/tcoeff if tcoeff else t
                else:
                    return t/tcoeff if tcoeff else t
            if numG % 2 == 1:
                return TensMul.from_data(S.Zero, [], [], [])
            elif numG > 4:
                t = t.substitute_indices((-DiracSpinorIndex.auto_right, -DiracSpinorIndex.auto_index), (DiracSpinorIndex.auto_left, DiracSpinorIndex.auto_index))
                a = t.split()
                ind1, lind1, rind1 = a[i].args[-1]
                ind2, lind2, rind2 = a[i + 1].args[-1]
                aa = a[:i] + a[i + 2:]
                t1 = tensor_mul(*aa)*g(ind1, ind2)*sg(lind1, rind1)*sg(lind2, rind2)
                t1 = t1.contract_metric(g)
                t1 = t1.contract_metric(sg)
                args = [t1]
                sign = 1
                for k in range(i + 2, j):
                    sign = -sign
                    ind2, lind2, rind2 = a[k].args[-1]
                    aa = a[:i] + a[i + 1:k] + a[k + 1:]
                    t2 = sign*tensor_mul(*aa)*g(ind1, ind2)*sg(lind1, rind1)*sg(lind2, rind2)
                    t2 = t2.contract_metric(g)
                    t2 = t2.contract_metric(sg)

                    t2 = GammaMatrixHead.simplify_gpgp(t2, False)
                    args.append(t2)
                t3 = TensAdd(*args)

                #aa = _tensorlist_contract_metric(aa, g(ind1, ind2))
                #t3 = t3.canon_bp()
                t3 = self._trace_single_line(t3)
                return t3
            else:
                a = t.split()
                if len(t.components) == 1:
                    if t.components[0] is DiracSpinorIndex.delta:
                        return 4  # FIXME only for D=4
                t1 = self._gamma_trace1(*a[i:j])
                a2 = a[:i] + a[j:]
                t2 = tensor_mul(*a2)
                t3 = t1*t2
                if not t3:
                    return t3
                t3 = t3.contract_metric(g)
                return t3

        if isinstance(t, TensAdd):
            a = [x.coeff*_trace_single_line1(x) for x in t.args]
            return TensAdd(*a)
        elif isinstance(t, TensMul):
            r = t.coeff*_trace_single_line1(t)
            return r
        else:
            return t

    def _gamma_trace1(self, *a):
        gctr = 4  # FIXME specific for d=4
        g = self.LorentzIndex.metric
        if not a:
            return gctr
        n = len(a)
        if n%2 == 1:
            #return TensMul.from_data(S.Zero, [], [], [])
            return S.Zero
        if n == 2:
            ind0 = a[0].args[-1][0]
            ind1 = a[1].args[-1][0]
            return gctr*g(ind0, ind1)
        if n == 4:
            ind0 = a[0].args[-1][0]
            ind1 = a[1].args[-1][0]
            ind2 = a[2].args[-1][0]
            ind3 = a[3].args[-1][0]

            return gctr*(g(ind0, ind1)*g(ind2, ind3) - \
               g(ind0, ind2)*g(ind1, ind3) + g(ind0, ind3)*g(ind1, ind2))

    @staticmethod
    def _kahane_simplify(coeff, tids):
        r"""
        This function cancels contracted elements in a product of four
        dimensional gamma matrices, resulting in an expression equal to the given
        one, without the contracted gamma matrices.

        Parameters
        ==========

        `coeff`     the coefficient of the tensor expression.
        `tids`      TIDS object representing the gamma matrix expression to simplify.

        Notes
        =====

        If spinor indices are given, the matrices must be given in
        the order given in the product.

        Algorithm
        =========

        The idea behind the algorithm is to use some well-known identities,
        i.e., for contractions enclosing an even number of `\gamma` matrices

        `\gamma^\mu \gamma_{a_1} \cdots \gamma_{a_{2N}} \gamma_\mu = 2 (\gamma_{a_{2N}} \gamma_{a_1} \cdots \gamma_{a_{2N-1}} + \gamma_{a_{2N-1}} \cdots \gamma_{a_1} \gamma_{a_{2N}} )`

        for an odd number of `\gamma` matrices

        `\gamma^\mu \gamma_{a_1} \cdots \gamma_{a_{2N+1}} \gamma_\mu = -2 \gamma_{a_{2N+1}} \gamma_{a_{2N}} \cdots \gamma_{a_{1}}`

        Instead of repeatedly applying these identities to cancel out all contracted indices,
        it is possible to recognize the links that would result from such an operation,
        the problem is thus reduced to a simple rearrangement of free gamma matrices.

        Examples
        ========

        When using, always remember that the original expression coefficient
        has to be handled separately

        >>> from sympy.physics.hep.gamma_matrices import GammaMatrix as G, DiracSpinorIndex as DS
        >>> from sympy.tensor.tensor import tensor_indices, tensorhead, TensMul, TensAdd
        >>> i0, i1, i2 = tensor_indices('i0:3', G.LorentzIndex)
        >>> s0,s1,s2,s3,s4,s5 = tensor_indices('s0:6', DS)
        >>> ta = G(i0)*G(-i0)
        >>> G._kahane_simplify(ta.coeff, ta._tids) - 4*DS.delta(DS.auto_left, -DS.auto_right)
        0
        >>> tb = G(i0)*G(i1)*G(-i0)
        >>> G._kahane_simplify(tb.coeff, tb._tids)
        -2*gamma(i1, auto_left, -auto_right)
        >>> t = G(i0, s0, -s1)*G(-i0,s1,-s2)
        >>> G._kahane_simplify(t.coeff, t._tids) - 4*DS.delta(s0, -s2)
        0
        >>> t = G(i0, s0, -s1)*G(-i0,s1,-s0)
        >>> G._kahane_simplify(t.coeff, t._tids)
        16

        If there are no contractions, the same expression is returned

        >>> tc = 3*G(i0)*G(i1)
        >>> G._kahane_simplify(tc.coeff, tc._tids)
        3*gamma(i0, auto_left, -S_0)*gamma(i1, S_0, -auto_right)

        References
        ==========

        [1] Algorithm for Reducing Contracted Products of gamma Matrices, Joseph Kahane, Journal of Mathematical Physics, Vol. 9, No. 10, October 1968.
        """

        for c in tids.components:
            if not(isinstance(tids.components[0], GammaMatrixHead)):
                raise ValueError('use only gamma matrices')
        n = len(tids.components)
        for p0, p1, c0, c1 in tids.dum:
            if p0 == 0:
                continue
            dc = abs(c0 - c1)
            if dc not in (1, n - 1):
                raise ValueError('wrong gamma matrix ordering')
        free = [_ for _ in tids.free if _[1] == 0]
        spinor_free = [_ for _ in tids.free if _[1] != 0]
        if len(spinor_free) == 2:
            spinor_free.sort(key=lambda x: x[2])
            assert spinor_free[0][1] == 1 and spinor_free[-1][1] == 2
            assert spinor_free[0][2] == 0
        elif spinor_free:
            raise ValueError('spinor indices do not match')

        dum = sorted([_ for _ in tids.dum if _[0] == 0 and _[1] == 0])

        if len(dum) == 0:  # or GammaMatrixHead:
            # no contractions in `expression`, just return it.
            return TensMul.from_TIDS(coeff, tids)

        # find the `first_dum_pos`, i.e. the position of the first contracted
        # gamma matrix, Kahane's algorithm as described in his paper requires the
        # gamma matrix expression to start with a contracted gamma matrix, this is
        # a workaround which ignores possible initial free indices, and re-adds
        # them later.
        dum_zip = list(zip(*dum))[2:]
        first_dum_pos = min(min(dum_zip[0]), min(dum_zip[1]))

        total_number = len(free) + len(dum)*2
        number_of_contractions = len(dum)

        free_pos = [None]*total_number
        for i in free:
            free_pos[i[2]] = i[0]

        # `index_is_free` is a list of booleans, to identify index position
        # and whether that index is free or dummy.
        index_is_free = [False]*total_number

        for i, indx in enumerate(free):
            if indx[1] != 0:
                raise ValueError("indx[1] should be equal to 0")
            index_is_free[indx[2]] = True

        # `links` is a dictionary containing the graph described in Kahane's paper,
        # to every key correspond one or two values, representing the linked indices.
        # All values in `links` are integers, negative numbers are used in the case
        # where it is necessary to insert gamma matrices between free indices, in
        # order to make Kahane's algorithm work (see paper).
        links = dict()
        for i in range(first_dum_pos, total_number):
            links[i] = []

        # `cum_sign` is a step variable to mark the sign of every index, see paper.
        cum_sign = -1
        # `cum_sign_list` keeps storage for all `cum_sign` (every index).
        cum_sign_list = [None]*total_number
        block_free_count = 0

        # multiply `resulting_coeff` by the coefficient parameter, the rest
        # of the algorithm ignores a scalar coefficient.
        resulting_coeff = S.One * coeff

        # initialize a lisf of lists of indices. The outer list will contain all
        # additive tensor expressions, while the inner list will contain the
        # free indices (rearranged according to the algorithm).
        resulting_indices = [[]]

        # start to count the `connected_components`, which together with the number
        # of contractions, determines a -1 or +1 factor to be multiplied.
        connected_components = 1

        # First loop: here we fill `cum_sign_list`, and draw the links
        # among consecutive indices (they are stored in `links`). Links among
        # non-consecutive indices will be drawn later.
        for i, is_free in enumerate(index_is_free):
            # if `expression` starts with free indices, they are ignored here;
            # they are later added as they are to the beginning of all
            # `resulting_indices` list of lists of indices.
            if i < first_dum_pos:
                continue

            if is_free:
                block_free_count += 1
                # if previous index was free as well, draw an arch in `links`.
                if block_free_count > 1:
                    links[i - 1].append(i)
                    links[i].append(i - 1)
            else:
                # Change the sign of the index (`cum_sign`) if the number of free
                # indices preceding it is even.
                cum_sign *= 1 if (block_free_count % 2) else -1
                if block_free_count == 0 and i != first_dum_pos:
                    # check if there are two consecutive dummy indices:
                    # in this case create virtual indices with negative position,
                    # these "virtual" indices represent the insertion of two
                    # gamma^0 matrices to separate consecutive dummy indices, as
                    # Kahane's algorithm requires dummy indices to be separated by
                    # free indices. The product of two gamma^0 matrices is unity,
                    # so the new expression being examined is the same as the
                    # original one.
                    if cum_sign == -1:
                        links[-1-i] = [-1-i+1]
                        links[-1-i+1] = [-1-i]
                if (i - cum_sign) in links:
                    if i != first_dum_pos:
                        links[i].append(i - cum_sign)
                    if block_free_count != 0:
                        if i - cum_sign < len(index_is_free):
                            if index_is_free[i - cum_sign]:
                                links[i - cum_sign].append(i)
                block_free_count = 0

            cum_sign_list[i] = cum_sign

        # The previous loop has only created links between consecutive free indices,
        # it is necessary to properly create links among dummy (contracted) indices,
        # according to the rules described in Kahane's paper. There is only one exception
        # to Kahane's rules: the negative indices, which handle the case of some
        # consecutive free indices (Kahane's paper just describes dummy indices
        # separated by free indices, hinting that free indices can be added without
        # altering the expression result).
        for i in dum:
            if i[0] != 0:
                raise ValueError("i[0] should be 0")
            if i[1] != 0:
                raise ValueError("i[1] should be 0")
            # get the positions of the two contracted indices:
            pos1 = i[2]
            pos2 = i[3]

            # create Kahane's upper links, i.e. the upper arcs between dummy
            # (i.e. contracted) indices:
            links[pos1].append(pos2)
            links[pos2].append(pos1)

            # create Kahane's lower links, this corresponds to the arcs below
            # the line described in the paper:

            # first we move `pos1` and `pos2` according to the sign of the indices:
            linkpos1 = pos1 + cum_sign_list[pos1]
            linkpos2 = pos2 + cum_sign_list[pos2]

            # otherwise, perform some checks before creating the lower arcs:

            # make sure we are not exceeding the total number of indices:
            if linkpos1 >= total_number:
                continue
            if linkpos2 >= total_number:
                continue

            # make sure we are not below the first dummy index in `expression`:
            if linkpos1 < first_dum_pos:
                continue
            if linkpos2 < first_dum_pos:
                continue

            # check if the previous loop created "virtual" indices between dummy
            # indices, in such a case relink `linkpos1` and `linkpos2`:
            if (-1-linkpos1) in links:
                linkpos1 = -1-linkpos1
            if (-1-linkpos2) in links:
                linkpos2 = -1-linkpos2

            # move only if not next to free index:
            if linkpos1 >= 0 and not index_is_free[linkpos1]:
                linkpos1 = pos1

            if linkpos2 >=0 and not index_is_free[linkpos2]:
                linkpos2 = pos2

            # create the lower arcs:
            if linkpos2 not in links[linkpos1]:
                links[linkpos1].append(linkpos2)
            if linkpos1 not in links[linkpos2]:
                links[linkpos2].append(linkpos1)

        # This loop starts from the `first_dum_pos` index (first dummy index)
        # walks through the graph deleting the visited indices from `links`,
        # it adds a gamma matrix for every free index in encounters, while it
        # completely ignores dummy indices and virtual indices.
        pointer = first_dum_pos
        previous_pointer = 0
        while True:
            if pointer in links:
                next_ones = links.pop(pointer)
            else:
                break

            if previous_pointer in next_ones:
                next_ones.remove(previous_pointer)

            previous_pointer = pointer

            if next_ones:
                pointer = next_ones[0]
            else:
                break

            if pointer == previous_pointer:
                break
            if pointer >=0 and free_pos[pointer] is not None:
                for ri in resulting_indices:
                    ri.append(free_pos[pointer])

        # The following loop removes the remaining connected components in `links`.
        # If there are free indices inside a connected component, it gives a
        # contribution to the resulting expression given by the factor
        # `gamma_a gamma_b ... gamma_z + gamma_z ... gamma_b gamma_a`, in Kahanes's
        # paper represented as  {gamma_a, gamma_b, ... , gamma_z},
        # virtual indices are ignored. The variable `connected_components` is
        # increased by one for every connected component this loop encounters.

        # If the connected component has virtual and dummy indices only
        # (no free indices), it contributes to `resulting_indices` by a factor of two.
        # The multiplication by two is a result of the
        # factor {gamma^0, gamma^0} = 2 I, as it appears in Kahane's paper.
        # Note: curly brackets are meant as in the paper, as a generalized
        # multi-element anticommutator!

        while links:
            connected_components += 1
            pointer = min(links.keys())
            previous_pointer = pointer
            # the inner loop erases the visited indices from `links`, and it adds
            # all free indices to `prepend_indices` list, virtual indices are
            # ignored.
            prepend_indices = []
            while True:
                if pointer in links:
                    next_ones = links.pop(pointer)
                else:
                    break

                if previous_pointer in next_ones:
                    if len(next_ones) > 1:
                        next_ones.remove(previous_pointer)

                previous_pointer = pointer

                if next_ones:
                    pointer = next_ones[0]

                if pointer >= first_dum_pos and free_pos[pointer] is not None:
                    prepend_indices.insert(0, free_pos[pointer])
            # if `prepend_indices` is void, it means there are no free indices
            # in the loop (and it can be shown that there must be a virtual index),
            # loops of virtual indices only contribute by a factor of two:
            if len(prepend_indices) == 0:
                resulting_coeff *= 2
            # otherwise, add the free indices in `prepend_indices` to
            # the `resulting_indices`:
            else:
                expr1 = prepend_indices
                expr2 = list(reversed(prepend_indices))
                resulting_indices = [expri + ri for ri in resulting_indices for expri in (expr1, expr2)]

        # sign correction, as described in Kahane's paper:
        resulting_coeff *= -1 if (number_of_contractions - connected_components + 1) % 2 else 1
        # power of two factor, as described in Kahane's paper:
        resulting_coeff *= 2**(number_of_contractions)

        # If `first_dum_pos` is not zero, it means that there are trailing free gamma
        # matrices in front of `expression`, so multiply by them:
        for i in range(0, first_dum_pos):
            [ri.insert(0, free_pos[i]) for ri in resulting_indices]

        resulting_expr = S.Zero
        for i in resulting_indices:
            temp_expr = S.One
            for j in i:
                temp_expr *= GammaMatrix(j)
            resulting_expr += temp_expr

        t = resulting_coeff * resulting_expr
        t1 = None
        if isinstance(t, TensAdd):
            t1 = t.args[0]
        elif isinstance(t, TensMul):
            t1 = t
        if t1:
            spinor_free1 = [_ for _ in t1._tids.free if _[1] != 0]
            if spinor_free1:
                if spinor_free:
                    t = t.substitute_indices((DiracSpinorIndex.auto_left, spinor_free[0][0]), (-DiracSpinorIndex.auto_right, spinor_free[-1][0]))
                else:
                    # FIXME trace
                    t = t*DiracSpinorIndex.delta(DiracSpinorIndex.auto_right, -DiracSpinorIndex.auto_left)
                    t = GammaMatrix.simplify_lines(t)
            else:
                if spinor_free:
                    t = t*DiracSpinorIndex.delta(spinor_free[0][0], spinor_free[-1][0])
                else:
                    t = t*4
        else:
            if spinor_free:
                t = t*DiracSpinorIndex.delta(spinor_free[0][0], spinor_free[-1][0])
            else:
                t = t*4
        return t

GammaMatrix = GammaMatrixHead()
