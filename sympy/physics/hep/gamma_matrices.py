from sympy import S
from sympy.tensor.tensor import TensorIndexType, tensorhead, TensorIndex,\
                                TensMul, TensorHead, tensorsymmetry, TensorType,\
                                TIDS, TensAdd, tensor_mul, TensExpr
from sympy.core.containers import Tuple
import collections


DiracSpinor = TensorIndexType('DiracSpinor', dim=4, dummy_fmt="S")


class _LorentzContainer(object):
    """
    Helper to collect Lorentz indices in various dimensions.

    It collects Lorentz TensorIndexType that have been implemented in the code,
    and stores them in a dict()
    """
    lorentz_types = dict()

    def __new__(cls, dim=4, eps_dim=None, dummy_fmt="L"):
        if (dim, eps_dim) in _LorentzContainer.lorentz_types:
            return _LorentzContainer.lorentz_types[(dim, eps_dim)]

        new_L = TensorIndexType("Lorentz", dim=dim, eps_dim=eps_dim, dummy_fmt=dummy_fmt)
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
    >>> i = tensor_indices('i', G.Lorentz)
    >>> G(i)
    gamma(i, auto_left, auto_right)

    Note that there is already an instance of GammaMatrixHead in four dimensions:
    GammaMatrix, which is simply declare as

    ``GammaMatrix = GammaMatrixHead()``

    >>> from sympy.physics.hep.gamma_matrices import GammaMatrix
    >>> from sympy.tensor.tensor import tensor_indices
    >>> i = tensor_indices('i', GammaMatrix.Lorentz)
    >>> GammaMatrix(i)
    gamma(i, auto_left, auto_right)

    To access the metric tensor

    >>> GammaMatrix.Lorentz.metric
    metric(Lorentz,Lorentz)

    """
    _gmhd = dict()

    def __new__(cls, dim=4, eps_dim=4):
        key = (dim, eps_dim)
        if key in GammaMatrixHead._gmhd:
            return GammaMatrixHead._gmhd[key]

        lorentz = _LorentzContainer(*key)

        gmh = TensorHead.__new__(cls, "gamma", TensorType(Tuple(lorentz, DiracSpinor, DiracSpinor), tensorsymmetry([1], [1], [1])), comm=0, matrix_behavior=True)
        GammaMatrixHead._gmhd[key] = gmh
        gmh.Lorentz = lorentz
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
    def _get_lines(ex):
        """
        returns ``(lines, traces)``
        where ``lines`` is the list of list of positions of a gamma line,
        ``traces`` is the list of list of traced gamma lines

        Examples
        ========
        from sympy.physics.hep.gamma_matrices import GammaMatrix
        from sympy.tensor.tensor import tensor_indices, TensorIndexType,  tensorhead
        >>> from sympy.physics.hep.gamma_matrices import GammaMatrix, DiracSpinor
        >>> from sympy.tensor.tensor import tensor_indices
        >>> i0,i1,i2,i3,i4,i5 = tensor_indices('i0:6', GammaMatrix.Lorentz)
        >>> s0,s1,s2,s3,s4,s5,s6,s7 = tensor_indices('s0:8', DiracSpinor)
        >>> G = GammaMatrix
        >>> t = G(i1,s1,-s2)*G(i4,s7,-s6)*G(i2,s2,-s3)*G(i3,s4,-s5)*G(i5,s6,-s7)
        >>> GammaMatrix._get_lines(t)
        ([[0, 2]], [[1, 4]])
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
        dum = tids.dum
        lines = []
        traces = []
        traces1 = []
        for p0, p1, c0, c1 in dum:
            # tr(G(i)) = 0
            if c0 == c1:
                return ([], [])
            if p0 == p1:
                if p0 > 0:
                    # case gamma(i,s0,-s1)in c0, gamma(j,-s0,s2) in c1;
                    # to deal with this case one could add to the position
                    # a flag for transposition;
                    # one could write [(c0, False), (c1, True)]
                    raise NotImplementedError
                else:
                    continue
            # if p0 == 2 then G in pos c0 is mult on the right by G in c1
            # if p0 == 1 then G in pos c1 is mult on the right by G in c0
            b0, b1 = (c0, c1) if p0 == 2 else (c1, c0)
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

        return lines, traces

    @staticmethod
    def simplify_lines(ex):
        """
        simplify a product of gamma matrices

        Examples
        ========
        >>> from sympy.physics.hep.gamma_matrices import GammaMatrix, DiracSpinor
        >>> from sympy.tensor.tensor import tensor_indices
        >>> i0,i1,i2,i3,i4,i5 = tensor_indices('i0:6', GammaMatrix.Lorentz)
        >>> s0,s1,s2,s3,s4,s5,s6,s7 = tensor_indices('s0:8', DiracSpinor)
        >>> G = GammaMatrix
        >>> t = G(i1,s1,-s2)*G(i4,s7,-s6)*G(i2,s2,-s3)*G(i3,s4,-s5)*G(i5,s6,-s7)
        >>> G.simplify_lines(t)
        4*gamma(i3, s4, -s5)*gamma(i1, s1, -S_0)*gamma(i2, S_0, -s3)*metric(i4, i5)

        """
        tids = ex._tids
        lines, traces = GammaMatrixHead._get_lines(ex)
        rest = []
        for line in lines:
            for y in line:
                rest.append(y)
        for line in traces:
            for y in line:
                rest.append(y)
        rest = [x for x in range(len(tids.components)) if x not in rest]
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
            if tx1.is_integer:
                tx1 = tx1*DiracSpinor.delta(first, last)
            tlines.append(tx1)
        traces = [GammaMatrix._trace_single_line(tensor_mul(*[x for i, x in enumerate(a) if i  in line])) for line in traces]
        res = tensor_mul(*([trest] + tlines + traces))
        return res


    @staticmethod
    def _simplify_single_line(expression):
        """
        Simplify single-line product of gamma matrices.
        """
        coeff = expression.coeff
        tids = expression._tids

        flag = False
        for i in expression.free:
            if (i[0]._tensortype == DiracSpinor):
                if not flag:
                    first_spinor_index = i[0]
                    flag = True
                last_spinor_index = i[0]

        tadd = GammaMatrixHead._kahane_simplify(coeff, tids)
        if not isinstance(tadd, TensExpr) and flag:
            tadd *= DiracSpinor.delta(first_spinor_index, last_spinor_index)

        return tadd

    def _trace_single_line(self, expression):
        """
        Evaluate the trace of a single gamma matrix line inside a ``TensExpr``.
        """
        tadd = GammaMatrix._simplify_single_line(expression)

        Lorentz = self.Lorentz
        D = Lorentz.dim

        if not isinstance(tadd, TensAdd):
            tadd = TensAdd(tadd)

        if tadd.args[0].get_indices() == [DiracSpinor.auto_left, DiracSpinor.auto_right] or tadd.rank == 0:
            # if rank of simplified expression is zero, return the dimension:
            assert len(tadd.args) == 1
            return D * tadd.args[0].coeff

        # Recurrence function to transform
        def even_trace_recursion(t):
            assert isinstance(t, TensMul)
            lorentz_indices = []
            for j in t.split():
                if not isinstance(j.components[0], GammaMatrixHead):
                    continue
                lorentz_indices += [i[0] for i in j.free if i[0]._tensortype == Lorentz]

            if len(lorentz_indices) == 0:
                return t
            elif len(lorentz_indices) == 2:
                return t.coeff * 4 * GammaMatrix.Lorentz.metric(lorentz_indices[0], lorentz_indices[1])

            a = t.split()
            for i in a:
                assert isinstance(i.components[0], GammaMatrixHead)
            ind1 = lorentz_indices[0]
            r_list = []
            metric_list = []
            sign = 1
            for k in range(1, len(a)):
                ind2 = lorentz_indices[k]
                aa1 = a[1:k] + a[k+1:]
                r_list.append(sign*tensor_mul(*aa1))
                metric_list.append(Lorentz.metric(ind1, ind2))
                sign *= -1
            p_list = [j*even_trace_recursion(i) for i, j in zip(r_list, metric_list)]
            return t.coeff * TensAdd(*p_list)

        res_add = []
        # in this for cycle the `TensMul` composing `tadd` get rearranged
        # according to their `DiracSpinor` free/dummy indices.
        for arg in tadd.args:
            # make the gamma matrix expression a circle with respect to its DiracSpinor indices:
            a = arg.substitute_indices((DiracSpinor.auto_right, -DiracSpinor.auto_index), (DiracSpinor.auto_left, DiracSpinor.auto_index))
            arg = S.One
            for i in a.split():
                arg *= i
            splits = arg.split()
            gamma_mats = [_ for _ in splits if isinstance(_.components[0], GammaMatrixHead)]
            gexpr = S.One
            while gamma_mats:
                gm1 = gamma_mats.pop(0)
                gexpr *= gm1
                if not gamma_mats:
                    break
                spn1, spn2 = [_[0] for _ in gm1.free if _[0]._tensortype == DiracSpinor]
                swappos = next(_ for _ in range(len(gamma_mats)) if gamma_mats[_].has(-spn1) or gamma_mats[_].has(-spn2))
                gamma_mats[0], gamma_mats[swappos] = gamma_mats[swappos], gamma_mats[0]

            res_add.append(even_trace_recursion(gexpr))

        t = TensAdd(*res_add)
        return t

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
        the order given in the product. If ``t`` is
        ``G(i0)*G(i1)`` or
        ``G(i0,s0,-s1)*G(i1,s1,-s2)*G(-i0,s2,-s3)`` or
        ``G(i0,True,-s1)*G(i1,s1,-s2)*G(-i0,s2,True)``,
        ``G._kahane_simplify(t.coeff, t._tids)`` gives
        ``-2*gamma(i1, auto_left, auto_right)``

        but if one reorders the terms as in
        ``t = G(-i0,s2,-s3)*G(i0,s0,-s1)*G(i1,s1,-s2)``
        one gets the wrong result ``4*gamma(i1, auto_left, auto_right)``

        If the result is a multiple of the identity, a number is given

        `` t = G(i0)*G(-i0)`` or
        ``t = G(i0,s0,-s1)*G(-i0,s1,-s2)``
        ``G._kahane_simplify(t.coeff, t._tids)`` gives 4

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

        >>> from sympy.physics.hep.gamma_matrices import GammaMatrix as G
        >>> from sympy.tensor.tensor import tensor_indices, tensorhead, TensMul, TensAdd
        >>> i0, i1, i2 = tensor_indices('i0:3', G.Lorentz)
        >>> ta = G(i0)*G(-i0)
        >>> sa = G._kahane_simplify(ta.coeff, ta._tids)
        >>> sa
        4
        >>> tb = G(i0)*G(i1)*G(-i0)
        >>> sb = G._kahane_simplify(tb.coeff, tb._tids)
        >>> sb
        -2*gamma(i1, auto_left, auto_right)

        If there are no contractions, the same expression is returned

        >>> tc = 3*G(i0)*G(i1)
        >>> sc = G._kahane_simplify(tc.coeff, tc._tids)
        >>> sc
        3*gamma(i0, auto_left, S_0)*gamma(i1, -S_0, auto_right)

        References
        ==========

        [1] Algorithm for Reducing Contracted Products of gamma Matrices, Joseph Kahane, Journal of Mathematical Physics, Vol. 9, No. 10, October 1968.
        """

        free = [_ for _ in tids.free if _[1] == 0]
        dum = sorted([_ for _ in tids.dum if _[0] == 0 and _[1] == 0])

        if len(dum) == 0:
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
            assert indx[1] == 0
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
        last_cum_sign = 1
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
                last_cum_sign = cum_sign

            cum_sign_list[i] = cum_sign

        # The previous loop has only created links between consecutive free indices,
        # it is necessary to properly create links among dummy (contracted) indices,
        # according to the rules described in Kahane's paper. There is only one exception
        # to Kahane's rules: the negative indices, which handle the case of some
        # consecutive free indices (Kahane's paper just describes dummy indices
        # separated by free indices, hinting that free indices can be added without
        # altering the expression result).
        for i in dum:
            assert i[0] == 0
            assert i[1] == 0
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

        return resulting_coeff * resulting_expr


GammaMatrix = GammaMatrixHead()
