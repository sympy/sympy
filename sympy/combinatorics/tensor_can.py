from sympy.combinatorics.permutations import Permutation, _af_rmul, _af_rmuln,\
  _af_invert, _af_new
from sympy.combinatorics.perm_groups import PermutationGroup, _orbit, \
  _orbit_transversal
from sympy.combinatorics.util import _distribute_gens_by_base, \
  _orbits_transversals_from_bsgs

"""
    References for tensor canonicalization:

    [1] R. Portugal "Algorithmic simplification of tensor expressions",
        J. Phys. A 32 (1999) 7779-7789

    [2] R. Portugal, B.F. Svaiter "Group-theoretic Approach for Symbolic
        Tensor Manipulation: I. Free Indices"
        arXiv:math-ph/0107031v1

    [3] L.R.U. Manssur, R. Portugal "Group-theoretic Approach for Symbolic
        Tensor Manipulation: II. Dummy Indices"
        arXiv:math-ph/0107032v1

    [4] xperm.c part of XPerm written by J. M. Martin-Garcia
        http://www.xact.es/index.html
"""


def dummy_sgs(dummies, sym, n):
    """
    Return the strong generators for dummy indices

    dummies   list of dummy indices,
              `dummies[2k], dummies[2k+1]` are paired indices
    sym       symmetry under interchange of contracted dummies
    sym =     None  no symmetry
              0 symmetric
              1 antisymmetric
    n         number of indices

    in base form the dummy indices are always in consecutive positions

    Examples
    ========
    >>> from sympy.combinatorics.tensor_can import dummy_sgs
    >>> dummy_sgs([2,3,4,5,6,7], 0, 8)
    [[0, 1, 3, 2, 4, 5, 6, 7, 8, 9], [0, 1, 2, 3, 5, 4, 6, 7, 8, 9], [0, 1, 2, 3, 4, 5, 7, 6, 8, 9], [0, 1, 4, 5, 2, 3, 6, 7, 8, 9], [0, 1, 2, 3, 6, 7, 4, 5, 8, 9]]
    """
    assert len(dummies) <= n
    res = []
    for j in dummies[::2]:
        a = range(n+2)
        if sym == 1:
            a[n] = n+1
            a[n+1] = n
        a[j], a[j+1] = a[j+1], a[j]
        res.append(a)
    for j in dummies[:-3:2]:
        a = range(n+2)
        a[j:j+4] = a[j+2],a[j+3],a[j],a[j+1]
        res.append(a)
    return res

def _min_dummies(dummies, sym, indices):
    """
    Return list of minima of the orbits of indices in group of dummies

    sym symmetry under interchange of contracted dummies
    sym = None  n symmetry
    sym = 0 symmetric
          1 antisymmetric

    Examples
    ========
    >>> from sympy.combinatorics.tensor_can import _min_dummies
    >>> _min_dummies([2,3,4,5,6,7], 0, range(10))
    [0, 1, 2, 2, 2, 2, 2, 2, 8, 9]
    """
    if sym is not None:
        if not dummies:
            return indices
        m = min(dummies)
        n = len(indices)
        res = indices[:]
        for c, i in enumerate(indices):
            if i in dummies:
                res[c] = m
        return res
    else:
        raise NotImplementedError


def _trace_S(s, j, b, S_cosets):
    """
    Return the representative h satisfying s[h[b]] == j

    If there is not such a representative return None
    """
    for h in S_cosets[b]:
        if s[h[b]] == j:
            return h
    return None

def _trace_D(gj, p_i, Dxtrav):
    """
    Return the representative h satisfying h[gj] == p_i

    If there is not such a representative return None
    """
    for h in Dxtrav:
        if h[gj] == p_i:
            return h
    return None

def _dumx_remove(dumx, p0):
    """
    remove p0 from dumx
    """
    if p0 not in dumx:
        return
    k = dumx.index(p0)
    if k % 2 == 0:
        p0_paired = dumx[k+1]
    else:
        p0_paired = dumx[k-1]
    dumx.remove(p0)
    dumx.remove(p0_paired)

def transversal2coset(size, base, transversal):
    a = []
    j = 0
    for i in range(size):
        if i in base:
            a.append(sorted(transversal[j].values()))
            j += 1
        else:
            a.append([range(size)])
    j = len(a)-1
    while a[j] == [range(size)]:
        j -= 1
    return a[:j+1]


def double_coset_can_rep(sym, sgens, g, sgs=None):
    """
    Butler-Portugal algorithm for tensor canonicalization with dummy indices

      sym     symmetry of the contructing metric
              0     symmetric
              1     antisymmetric
              None  no symmetry (not implemented yet)

      sgens   generators of the group of slot symmetries of a minimal BSGS
      g       permutation representing the tensor
      sgs     if it is not None, sgens is a strong generating set of S
              sgs is the tuple (S_transversals, S_base)

    A tensor with dummy indices can be represented in a number
    of equivalent ways which typically grows exponentially with
    the number of indices. To be able to establish if two tensors
    with many indices are equal becomes computationally very slow
    in absence of an efficient algorithm.

    The Butler-Portugal algorithm [3] is an efficient algorithm to
    put tensors in canonical form, solving the above problem.

    Let us explain the conventions by an example.

    Given a tensor T^{d3 d2 d1}_{d1 d2 d3}
    where T^{a0,a1,a2,a3,a4,a5} = -T^{a2,a1,a0,a3,a4,a5}
                                = -{Ta4,a1,a2,a3,a0,a5}
    and symmetric metric, find the tensor equivalent to it which
    is the lowest under the ordering of indices:
       lexicographic ordering d1, d2, d3
       contravariant index < covariant index
    that is the canonical form of the tensor.

    that is the canonical form of the tensor.

    The canonical form is -T^{d1 d2 d3}_{d1 d2 d3}
    obtained using T^{a0,a1,a2,a3,a4,a5} = -T^{a2,a1,a0,a3,a4,a5}

    To convert this problem in the input for this function,
    use the following labelling of the index names
    (- for covariant for short)
      d1, -d1, d2, -d2, d3, -d3
      0,  1.   2,  3,   4,  5

    T^{d3 d2 d1}_{d1 d2 d3} corresponds to
    g = [4,2,0,1,3,5,6,7]

    where the last two indices are for the sign

    sgens = [Permutation([2,1,0,3,4,5,7,6]), Permutation([4,1,2,3,0,5,7,6])]

    sgens[0] is the slot symmetry -(0,2)
    T^{a0,a1,a2,a3,a4,a5} = -T^{a2,a1,a0,a3,a4,a5};    -(0,2)

    sgens[1] is the slot symmetry -(0,4)
    T^{a0,a1,a2,a3,a4,a5} = -T^{a4,a1,a2,a3,a0,a5};    -(0,4)

    The dummy symmetry group D is generated by the strong base generators
    [(0,1),(2,3),(4,5),(0,1)(2,3),(2,3)(4,5)]

    The dummy symmetry acts from the left
    d = [1,0,2,3,4,5,6,7]  exchange d1 -> -d1
    T(d3,d2,d1,-d1,-d2,-d3) == T(d3,d2,-d1,d1,-d2,-d3)
    g=[4,2,0,1,3,5,6,7]     == [4,2,1,0,3,5,6,7] = _af_rmul(d, g)
    which differs from _af_rmul(g, d)

    The slot symmetry acts from the right
    s = [2,1,0,3,4,5,7,6]  exchanges slots 0 and 2 and changes sign
    T(d3,d2,d1,-d1,-d2,-d3) == -T(d1,d2,d3,-d1,-d2,-d3)
    g=[4,2,0,1,3,5,6,7]     == [0,2,4,1,3,5,7,6] = _af_rmul(g, s)

    The double coset D*g*S consists of permutations d*g*s corresponding
    to equal tensors; choose as representative the tensor with indices
    ordered lexicographically according to [d1, -d1, d2, -d2, d3, -d3]
    that is rep = min(D*g*S) = min([d*g*s for d in D for s in S])

    The indices are fixed one by one; first choose the lowest index
    for slot 0, then the lowest remaining index for slot 1, etc.
    Doing this one obtains a chain of stabilizers
    S -> S_b0 -> S_{b0,b1} -> ... and
    D -> D_p0 -> D_{p0,p1} -> ...
    where [b0, b1, ...] = range(b) is a base of the symmetric group;
    the strong base b_S of S is an ordered sublist of it;
    therefore it is sufficient to compute once the
    strong base generators of S using the Schreier-Sims algorithm;
    the stabilizers of the strong base generators are the
    strong base generators of the stabilizer subgroup.

    dbase = [p0,p1,...] is not in general in lexicographic order,
    so that one must recompute the strong base generators each time;
    however this is trivial, there is no need to use the Schreier-Sims
    algorithm for D.


    The algorithm keeps a TAB of elements (s_i, d_i, h_i)
    where h_i = d_i*g*s_i satisfying h_i*b_j = p_j for 0 <= j < i
    starting from s_0 = id, d_0 = id, h_0 = g
    The equations h_0*b_0 = p_0, h_1*b_1 = p_1 are solved in this order,
    choosing each time the lowest possible value of p_i
    For j < i
    d_i*g*s_i*S_{b_0,...,b_{i-1}}*b_j = D_{p_0,...,p_{i-1}}*p_j
    so that for dx in D_{p_0,...,p_{i-1}} and sx in S_{base[0],...,base[i-1]}
    one has dx*d_i*g*s_i*sx*b_j = p_j
    Search for dx, sx such that this equation holds for j = i;
    it can be written as s_i*sx*b_j = J, dx*d_i*g*J = p_j
    sx*b_j = s_i**-1*J; sx = trace(s_i**-1, S_{b_0,...,b_{i-1}})
    dx**-1*p_j = d_i*g*J; dx = trace(d_i*g*J, D_{p_0,...,p_{i-1}})

    s_{i+1} = s_i*trace(s_i**-1*J, S_{b_0,...,b_{i-1}})
    d_{i+1} = trace(d_i*g*J, D_{p_0,...,p_{i-1}})**-1*d_i
    h_{i+1}*b_i = d_{i+1}*g*s_{i+1}*b_i = p_i

    h_n*b_j = p_j for all j, so that h_n is the solution

    Examples
    ========

    >>> from sympy.combinatorics.permutations import Permutation
    >>> from sympy.combinatorics.perm_groups import PermutationGroup
    >>> from sympy.combinatorics.tensor_can import double_coset_can_rep, get_transversals
    >>> gens = [Permutation(x) for x in [[2,1,0,3,4,5,7,6], [4,1,2,3,0,5,7,6]]]
    >>> base = [0, 2]
    >>> g = Permutation([4,2,0,1,3,5,6,7])
    >>> transversals = get_transversals(base, gens)
    >>> double_coset_can_rep(0, gens, g, (transversals, base))
    [0, 1, 2, 3, 4, 5, 7, 6]

    T^{d3}_{d1 d2}^{d1}_{d3}^{d2}
    = -T^{d3}_{d1 d3}^{d1}_{d2}^{d2}   under slot symmetry -(2,4)
    = T_{d3 d1}^{d3}^{d1}_{d2}^{d2}    under slot symmetry -(0,2)
    = T^{d3}_{d1 d3}^{d1}_{d2}^{d2}    symmetric metric
    = 0
    >>> g = [4,1,3,0,5,2,6,7]
    >>> double_coset_can_rep(0, gens, g, (transversals, base))
    0
    """
    if sym is None:
        raise NotImplementedError
    if isinstance(g, Permutation):
        size = g.size
        g = g.array_form
    else:
        size = len(g)
        sortedg = sorted(g)
        if sortedg != range(size):
            raise ValueError('g must be a permutation')
        size1 = size - 1
    dummies = range(size - 2)
    num_dummies = size - 2
    dumx = range(num_dummies)
    base = range(num_dummies)
    # strong base generators for Sx; start with Sx=S
    if not sgs:
        S = PermutationGroup(sgens)
        b_S, strong_gens = S.schreier_sims_incremental()
        if sorted(b_S) != b_S:
            # TODO use baseswap to get a BSGS with ordered strong baswe
            raise NotImplementedError
        strong_gens_distr = _distribute_gens_by_base(b_S, strong_gens)
        b_S, S_transversals = _orbits_transversals_from_bsgs(b_S,
                               strong_gens_distr)
        sgensx = [h._array_form for h in strong_gens]
    else:
        sgensx = [h.array_form for h in sgens]
        S_transversals, b_S = sgs
    b_S = b_S[:]
    if b_S:
        S_transversals = transversal2coset(size, b_S, S_transversals)
    # strong generating set for D
    dsgsx = dummy_sgs(dumx, sym, num_dummies)
    ginv = _af_invert(g)
    idn = range(size)
    # TAB = list of entries (s, d, h) where h = _af_rmuln(d,g,s)
    # for short, in the following d*g*s means _af_rmuln(d,g,s)
    TAB = [(idn, idn, g)]
    for i in range(size - 2):
        b = base[i]
        testb = b in b_S and sgensx
        sgensx1 = [_af_new(_) for _ in sgensx]
        if testb:
            deltab = _orbit(size, sgensx1, b)
        else:
            deltab = set([b])
        # p1 = min(IMAGES) = min(Union D_p*h*deltab for h in TAB)
        md = _min_dummies(dumx, 0, dummies)
        p_i = min([min([md[h[x]] for x in deltab]) for s,d,h in TAB])
        dsgsx1 = [_af_new(_) for _ in dsgsx]
        Dxtrav = _orbit_transversal(size, dsgsx1, p_i, False, af=True) \
                if dumx else None
        if Dxtrav:
            Dxtrav = [_af_invert(x) for x in Dxtrav]
        deltap = dumx if p_i in dumx else [p_i]
        TAB1 = []
        nTAB = len(TAB)
        while TAB:
            s, d, h = TAB.pop()
            if min([md[h[x]] for x in deltab]) != p_i:
                continue
            deltab1 = [x for x in deltab if md[h[x]] == p_i]
            # NEXT = s*deltab1 intersection (d*g)**-1*deltap
            dg = _af_rmul(d, g)
            dginv = _af_invert(dg)
            sdeltab = [s[x] for x in deltab1]
            gdeltap = [dginv[x] for x in deltap]
            NEXT = [x for x in sdeltab if x in gdeltap]
            # d, s satisfy
            # d*g*s*base[i-1] = p_{i-1}; using the stabilizers
            # d*g*s*S_{base[0],...,base[i-1]}*base[i-1] =
            # D_{p_0,...,p_{i-1}}*p_{i-1}
            # so that to find d1, s1 satisfying d1*g*s1*b = p_i
            # one can look for dx in D_{p_0,...,p_{i-1}} and
            # sx in S_{base[0],...,base[i-1]}
            # d1 = dx*d; s1 = s*sx
            # d1*g*s1*b = dx*d*g*s*sx*b = p_i
            for j in NEXT:
                if testb:
                    # solve s1*b = j with s1 = s*sx for some element sx
                    # of the stabilizer of ..., base[i-1]
                    # sx*b = s**-1*j; sx = _trace_S(s, j,...)
                    # s1 = s*trace_S(s**-1*j,...)
                    s1 = _trace_S(s, j, b, S_transversals)
                    if not s1:
                        continue
                    else:
                        s1 = [s[ix] for ix in s1]
                else:
                    s1 = s

                #assert s1[b] == j  # invariant
                # solve d1*g*j = p_i with d1 = dx*d for some element dg
                # of the stabilizer of ..., p_{i-1}
                # dx**-1*p_i = d*g*j; dx**-1 = trace_D(d*g*j,...)
                # d1 = trace_D(d*g*j,...)**-1*d
                # to save an inversion in the inner loop; notice we did
                # Dxtrav = [perm_af_invert(x) for x in Dxtrav] out of the loop
                if dumx:
                    d1 = _trace_D(dg[j], p_i, Dxtrav)
                    if not d1:
                        continue
                else:
                    if p_i != dg[j]:
                        continue
                    d1 = idn
                #assert d1[p_i] == dg[j]  # invariant
                d1 = [d1[ix] for ix in d]
                h1 = [d1[g[ix]] for ix in s1]
                #assert h1[b] == p_i  # invariant
                TAB1.append((s1, d1, h1))

        # if TAB contains equal permutations, keep only one of them;
        # if TAB contains equal permutations up to the sign, return 0
        TAB1.sort(key=lambda x: x[-1])
        nTAB1 = len(TAB1)
        prev = [0]*size
        while TAB1:
            s, d, h = TAB1.pop()
            if h[:-2] == prev[:-2]:
                if h[-1] != prev[-1]:
                    return 0
            else:
                TAB.append((s, d, h))
            prev = h

        # stabilize the SGS
        sgensx = [h for h in sgensx if h[b] == b]
        if b in b_S:
            b_S.remove(b)
        _dumx_remove(dumx, p_i)
        dsgsx = dummy_sgs(dumx, sym, num_dummies)
    return TAB[0][-1]


def canonical_free(base, gens, g, num_free):
    """
    canonicalization of a tensor with respect to free indices
    choosing the minimum with respect to lexicographical ordering

    base, gens  BSGS for slot permutation group
    g           permutation representing the tensor
    num_free    number of free indices
    The indices must be ordered with first the free indices

    see explanation in double_coset_can_rep
    The algorithm is given in [2]

    Examples
    ========
    >>> from sympy.combinatorics import Permutation
    >>> from sympy.combinatorics.tensor_can import canonical_free
    >>> gens = [[1,0,2,3,5,4], [2,3,0,1,4,5],[0,1,3,2,5,4]]
    >>> gens = [Permutation(h) for h in gens]
    >>> base = [0, 2]
    >>> g = Permutation([2, 1, 0, 3, 4, 5])
    >>> canonical_free(base, gens, g, 4)
    ([0, 3, 1, 2, 5, 4], [])

    Consider the product of Riemann tensors
    `T = R^{a}_{d0}^{d1,d2}*R_{d2,d1}^{d0,b}`
    The order of the indices is [a,b,d0,-d0,d1,-d1,d2,-d2]
    The permutation corresponding to the tensor is
    g = [0,3,4,6,7,5,2,1,8,9]
    Use the slot symmetries to get `T` is the form which is the minimum
    in lexicographic order
    `R^{a}_{d0}^{d1,d2}*R^{b,d0}_{d1,d2}` corresponding to
    `[0, 3, 4, 6, 1, 2, 5, 7, 8, 9]`
    >>> from sympy.combinatorics.tensor_can import riemann_bsgs, tensor_gens
    >>> base, gens = riemann_bsgs
    >>> size, sbase, sgens = tensor_gens(base, gens, [[],[]], 0)
    >>> g = Permutation([0,3,4,6,7,5,2,1,8,9])
    >>> canonical_free(sbase, [Permutation(h) for h in sgens], g, 2)
    ([0, 3, 4, 6, 1, 2, 5, 7, 8, 9], [])
    """
    g = g.array_form
    size = len(g)
    if not base:
        return g[:], [range(size)]
    K = [x.array_form for x in gens]

    transversals = get_transversals(base, K)
    cosets = transversal2coset(size, base, transversals)
    m = len(base)
    for x in sorted(g[:-2]):
        if x not in base:
            base.append(x)
    lambd = g
    K1 = []
    for i in range(m):
        b = base[i]
        delta = [x[b] for x in cosets[b]]
        delta1 = [lambd[x] for x in delta]
        delta1min = min(delta1)

        k = delta1.index(delta1min)
        p = delta[k]
        for omega in cosets[b]:
            if omega[b] == p:
                break
        lambd = _af_rmul(lambd, omega)
        K = [px for px in K if px[b] == b]
    return lambd, K


def _get_map_slots(size, fixed_slots):
    res = range(size)
    pos = 0
    for i in range(size):
      if i in fixed_slots:
          continue
      res[i] = pos
      pos += 1
    return res

def _lift_sgens(size, fixed_slots, free, s):
    a = []
    j = k = 0
    fd = zip(fixed_slots, free)
    fd = [y for x,y in sorted(fd)]
    num_free = len(free)
    for i in range(size):
        if i in fixed_slots:
            a.append(fd[k])
            k += 1
        else:
            a.append(s[j] + num_free)
            j += 1
    return a

def canonicalize(g, dummies, msym, *v):
    """
    canonicalize tensor formed by tensors of the different types

    g  permutation representing the tensor
    dummies  list of dummy indices; the dummy indices must come
             after the free indices, and put in order
             contravariant, covariant
             [d0, -d0, d1,-d1,...]

    msym symmetry of the dummy index metric

    v is a list of (base_i, gens_i, n_i, sym_i) for tensors of type `i`
         base_i, gens_i BSGS for tensors of this type
         The BSGS should have minimal base under lexicographic ordering;
         if not, canonicalize_naive is used, whichis much slower.

    n_i  number ot tensors of type `i`

    sym_i symmetry under exchange of two component tensors of type `i`
          None  no symmetry
          0     commuting
          1     anticommuting

    Return 0 if the tensor is zero, else return the array form of
    the permutation representing the canonical form of the tensor.

    Algorithm
    =========

    First one uses canonical_free to get the minimum tensor under
    lexicographic order, using only the slot symmetries.
    If the component tensors have not minimal BSGS, canonicalize_naive
    is used instead.
    Compute the residual slot symmetry keeping fixed the free indices
    using tensor_gens(base, gens, list_free_indices, sym)
    Reduce the problem eliminating the free indices.
    Then use double_coset_can_rep and lift back the result reintroducing
    the free indices.

    Examples
    ========
    >>> from sympy.combinatorics.tensor_can import get_symmetric_group_sgs, canonicalize
    >>> from sympy.combinatorics import Permutation, PermutationGroup
    >>> g = Permutation([1,3,2,0,4,5])
    >>> base2, gens2 = get_symmetric_group_sgs(2)
    >>> canonicalize(g, [2, 3], 0, (base2, gens2, 2, 0))
    [0, 2, 1, 3, 4, 5]
    """
    from sympy.combinatorics.testutil import canonicalize_naive
    if not msym in [0, 1, None]:
        raise ValueError('msym must be 0, 1 or None')
    size = g.size
    num_tensors = 0
    v1 = []
    for i in range(len(v)):
        base_i, gens_i, n_i, sym_i = v[i]
        # check that the BSGS is minimal;
        # this property is used in double_coset_can_rep;
        # if it is not minimal use canonicalize_naive
        # TODO use baseswap to find a minimal BSGS if this occurs.
        if not _is_minimal_bsgs(base_i, gens_i):
            mbsgs = get_minimal_bsgs(base_i, gens_i)
            if not mbsgs:
                can = canonicalize_naive(g, dummies, msym, *v)
                return can
            base_i, gens_i = mbsgs
        v1.append((base_i, gens_i, [[]]*n_i, sym_i))
        num_tensors += n_i
    # slot symmetry of the tensor
    size1, sbase, sgens = gens_products(*v1)
    assert size == size1
    free = [i for i in range(size-2) if i not in dummies]
    num_free = len(free)

    sgens = [_af_new(h) for h in sgens]
    # g1 minimal tensor under slot symmetry
    g1, gens1 = canonical_free(sbase, sgens, g, num_free)
    if not dummies:
        return g1
    # save the sign of g1
    sign = 0 if g1[-1] == size-1 else 1

    # the free indices are kept fixed.
    # Determine free_i, the list of slots of tensors which are fixed
    # since they are occupied by free indices, which are fixed.
    start = 0
    for i in range(len(v)):
        free_i = []
        base_i, gens_i, n_i, sym_i = v[i]
        # TODO deal with case gens_i == []
        len_tens = len(gens_i[0]) - 2
        # for each component tensor get a list od fixed islots
        for j in range(n_i):
            # get the elements corresponding to the component tensor
            h = g1[start:(start + len_tens)]
            fr = []
            # get the positions of the fixed elements in h
            for k in free:
                if k in h:
                    fr.append(h.index(k))
            free_i.append(fr)
            start += len_tens
        v1[i] = (base_i, gens_i, free_i, sym_i)
    # BSGS of the tensor with fixed free indices
    # if tensor_gens fails in gens_product, use canonicalize_naive
    size, sbase, sgens = gens_products(*v1)

    # reduce the permutations getting rid of the free indices
    pos_dummies = [g1.index(x) for x in dummies]
    pos_free = [g1.index(x) for x in range(num_free)]
    size_red = size - num_free
    g1_red = [x - num_free for x in g1 if x in dummies]
    if sign:
        g1_red.extend([size_red - 1, size_red - 2])
    else:
        g1_red.extend([size_red - 2, size_red - 1])
    map_slots = _get_map_slots(size, pos_free)
    sbase_red = [map_slots[i] for i in sbase if i not in pos_free]
    sgens_red = [[map_slots[i] for i in y if i not in pos_free] for y in sgens]
    transv_red = get_transversals(sbase_red, sgens_red)
    g1_red = _af_new(g1_red)
    sgens_red = [_af_new(h) for h in sgens_red]
    g2 = double_coset_can_rep(msym, sgens_red, g1_red, (transv_red, sbase_red))
    if g2 == 0:
        return 0
    # lift to the case with the free indices
    g3 = _lift_sgens(size, pos_free, free, g2)
    return g3


def perm_af_direct_product(gens1, gens2, signed=False):
    """
    direct products of the generators gens1 and gens2

    Examples
    ========
    >>> from sympy.combinatorics.tensor_can import perm_af_direct_product
    >>> gens1 = [[1,0,2,3], [0,1,3,2]]
    >>> gens2 = [[1,0]]
    >>> perm_af_direct_product(gens1, gens2, False)
    [[1, 0, 2, 3, 4, 5], [0, 1, 3, 2, 4, 5], [0, 1, 2, 3, 5, 4]]
    >>> gens1 = [[1,0,2,3,5,4], [0,1,3,2,4,5]]
    >>> gens2 = [[1,0,2,3]]
    >>> perm_af_direct_product(gens1, gens2, True)
    [[1, 0, 2, 3, 4, 5, 7, 6], [0, 1, 3, 2, 4, 5, 6, 7], [0, 1, 2, 3, 5, 4, 6, 7]]
    """
    gens1 = [list(x) for x in gens1]
    gens2 = [list(x) for x in gens2]
    s = 2 if signed else 0
    n1 = len(gens1[0]) - s
    n2 = len(gens2[0]) - s
    start = list(range(n1))
    end = list(range(n1, n1 + n2))
    if signed:
        gens1 = [gen[:-2] + end + [gen[-2]+n2, gen[-1]+n2] for gen in gens1]
        gens2 = [start + [x + n1 for x in gen] for gen in gens2]
    else:
        gens1 = [gen + end for gen in gens1]
        gens2 = [start + [x + n1 for x in gen] for gen in gens2]

    res = gens1 + gens2

    return res

def bsgs_direct_product(base1, gens1, base2, gens2, signed=False):
    """
    direct product of two BSGS

    base1    base of the first BSGS
    gens1    strong generating sequence of the first BSGS
    base2, gens2   similarly for the second BSGS
    signed   flag for signed permutations

    Examples
    ========
    >>> from sympy.combinatorics.tensor_can import (get_symmetric_group_sgs, bsgs_direct_product)
    >>> base1, gens1 = get_symmetric_group_sgs(1)
    >>> base2, gens2 = get_symmetric_group_sgs(2)
    >>> bsgs_direct_product(base1, gens1, base2, gens2, 1)
    ([1], [[0, 2, 1, 3, 4]])
    """
    s = 2 if signed else 0
    n1 = len(gens1[0]) - s
    base = base1[:]
    base += [x + n1 for x in base2]
    gens = perm_af_direct_product(gens1, gens2, signed)
    size = len(gens[0])
    id_af = range(size)
    gens = [h for h in gens if h != id_af]
    if not gens:
        gens = [id_af]
    return base, gens

def get_symmetric_group_sgs(n, sym=0):
    """
    Return base, gens of the minimal BSGS for (anti)symmetric group
    with n elements

    Examples
    ========
    >>> from sympy.combinatorics.tensor_can import get_symmetric_group_sgs
    >>> get_symmetric_group_sgs(3)
    ([0, 1], [[1, 0, 2, 3, 4], [0, 2, 1, 3, 4]])
    """
    if n == 1:
        return [], [range(3)]
    gens = [Permutation(n-1)(i, i+1).array_form for i in range(n-1)]
    if sym == 0:
        gens = [x + [n, n+1] for x in gens]
    else:
        gens = [x + [n+1, n] for x in gens]
    base = range(n-1)
    return base, gens

riemann_bsgs = [0, 2], [[1,0,2,3,5,4], [0,1,3,2,5,4], [2,3,0,1,4,5]]

def get_transversals(base, gens):
    """
    Return transversals for the group with BSGS base, gens
    """
    if not base:
        return []
    if not isinstance(gens[0], Permutation):
        gens = [_af_new(h) for h in gens]
    stabs =  _distribute_gens_by_base(base, gens)
    orbits, transversals = _orbits_transversals_from_bsgs(base, stabs)
    transversals = [dict((x, h._array_form) for x, h in y.items()) for y in \
            transversals]
    return transversals


def _is_minimal_bsgs(base, gens):
    """
    Check if the BSGS has minimal base under lexigographic order.

    base, gens BSGS

    Examples
    ========
    >>> from sympy.combinatorics.tensor_can import riemann_bsgs, _is_minimal_bsgs
    >>> _is_minimal_bsgs(*riemann_bsgs)
    True
    >>> riemann_bsgs1 = ([2, 0], [[1,0,2,3,5,4], [2,3,0,1,4,5]])
    >>> _is_minimal_bsgs(*riemann_bsgs1)
    False
    """
    base1 = []
    sgs1 = gens[:]
    size = len(gens[0])
    for i in range(size):
        if not all(h[i] == i for h in sgs1):
            base1.append(i)
            sgs1 = [h for h in sgs1 if h[i] == i]
    return base1 == base

def get_minimal_bsgs(base, gens):
    """
    Compute a minimal GSGS

    base, gens BSGS

    If base, gens is a minimal BSGS return it; else return a minimal BSGS
    if it fails in finding one, it returns None

    TODO: use baseswap in the case in which if it fails in finding a
    minimal BSGS

    Examples
    ========

    >>> from sympy.combinatorics.tensor_can import get_minimal_bsgs
    >>> riemann_bsgs1 = ([2, 0], [[1,0,2,3,5,4], [2,3,0,1,4,5]])
    >>> get_minimal_bsgs(*riemann_bsgs1)
    ([0, 2], [[1, 0, 2, 3, 5, 4], [2, 3, 0, 1, 4, 5], [0, 1, 3, 2, 5, 4]])
    """
    if _is_minimal_bsgs(base, gens):
        return base, gens
    if isinstance(gens[0], list):
        gens = [Permutation(h) for h in gens]
    G = PermutationGroup(gens)
    base, gens = G.schreier_sims_incremental()
    gens = [h._array_form for h in gens]
    if not  _is_minimal_bsgs(base, gens):
        raise None
    return base, gens

def tensor_gens(base, gens, list_free_indices, sym=0):
    """
    Returns size, res_base, res_gens BSGS for n tensors of the same type

    base, gens BSGS for tensors of this type
    list_free_indices  list of the slots occupied by fixed indices
                       for each of the tensors

    sym symmetry under commutation of two tensors
    sym   None  no symmetry
    sym   0     commuting
    sym   1     anticommuting

    Examples
    ========
    >>> from sympy.combinatorics.tensor_can import tensor_gens, get_symmetric_group_sgs

    two symmetric tensors with 3 indices without free indices
    >>> base, gens = get_symmetric_group_sgs(3)
    >>> tensor_gens(base, gens, [[], []])
    (8, [0, 1, 3, 4], [[1, 0, 2, 3, 4, 5, 6, 7], [0, 2, 1, 3, 4, 5, 6, 7], [0, 1, 2, 4, 3, 5, 6, 7], [0, 1, 2, 3, 5, 4, 6, 7], [3, 4, 5, 0, 1, 2, 6, 7]])

    two symmetric tensors with 3 indices with free indices in slot 1 and 0
    >>> tensor_gens(base, gens, [[1],[0]])
    (8, [0, 4], [[2, 1, 0, 3, 4, 5, 6, 7], [0, 1, 2, 3, 5, 4, 6, 7]])

    four symmetric tensors with 3 indices, two of which with free indices

    """
    def _get_bsgs(G, base, gens, free_indices):
        """
        return the BSGS for G.pointwise_stabilizer(free_indices)
        """
        if not free_indices:
            return base[:], gens[:]
        else:
            H = G.pointwise_stabilizer(free_indices)
            base, sgs = H.schreier_sims_incremental()
            return base, [h._array_form for h in sgs]

    # if not base there is no slot symmetry for the component tensors
    # if list_free_indices.count([]) < 2 there is no commutation symmetry
    # so there is no resulting slot symmetry
    if not base and list_free_indices.count([]) < 2:
        n = len(list_free_indices)
        size = len(gens[0])
        size = n*(len(gens[0]) - 2) + 2
        return size, [], [range(size)]

    # if any(list_free_indices) one needs to compute the pointwise
    # stabilizer, so G is needed
    if any(list_free_indices):
        G = PermutationGroup([_af_new(h) for h in gens])
    else:
        G = None

    # no_free list of lists of indices for component tensors without fixed
    # indices
    no_free = []
    size = len(gens[0])
    id_af = range(size)
    num_indices = size - 2
    if not list_free_indices[0]:
        no_free.append(range(num_indices))
    res_base, res_gens = _get_bsgs(G, base, gens, list_free_indices[0])
    for i in range(1, len(list_free_indices)):
        base1, gens1 = _get_bsgs(G, base, gens, list_free_indices[i])
        res_base, res_gens = bsgs_direct_product(res_base, res_gens,
            base1, gens1, 1)
        if not list_free_indices[i]:
            no_free.append(range(size-2, size - 2 + num_indices))
        size += num_indices
    nr = size - 2
    res_gens = [h for h in res_gens if h != id_af]
    # if sym == None there are no generators for commuting tensors
    if sym == None:
        if not res_gens:
            res_gens = id_af
        return size, res_base, res_gens

    # if the component tensors have moinimal BSGS, so is their direct
    # product P; the slot symmetry group is S = P*C, where C is the group
    # to (anti)commute the component tensors with no free indices
    # a stabilizer has the property S_i = P_i*C_i;
    # the BSGS of P*C has SGS_P + SGS_C and the base is
    # the ordered union of the bases of P and C.
    # If P has minimal BSGS, so has S with this base.
    base_comm = []
    for i in range(len(no_free) - 1):
        ind1 = no_free[i]
        ind2 = no_free[i + 1]
        a = range(ind1[0])
        a.extend(ind2)
        a.extend(ind1)
        base_comm.append(ind1[0])
        a.extend(range(ind2[-1]+1, nr))
        if sym == 0:
            a.extend([nr, nr + 1])
        else:
            a.extend([nr + 1, nr])
        res_gens.append(a)
    # each base is ordered; order the union of the two bases
    for i in base_comm:
        if i not in res_base:
            res_base.append(i)
    res_base.sort()
    # if commutator generators have been added, res_base might not fix
    # them, so add elements to it.
    fixed_gens = res_gens
    for i in res_base:
        fixed_gens = [h for h in fixed_gens if h[i] == i]
    if fixed_gens:
        for i in range(size):
            if i in res_base:
                continue
            if any(h[i] != i for h in fixed_gens):
                res_base.append(i)
                fixed_gens = [h for h in fixed_gens if h[i] == i]
    if not res_gens:
        res_gens = [id_af]

    return size, res_base, res_gens


def gens_products(*v):
    """
    Returns size, res_base, res_gens BSGS for n tensors of different types

    v is a sequence of (base_i, gens_i, free_i, sym_i)
    where
    base_i, gens_i  BSGS of tensor of type `i`
    free_i          list of the fixed slots for each of the tensors
                    of type `i`; if there are `n_i` tensors of type `i`
                    and none of them have fixed slots, `free = [[]]*n_i`
    sym   0 (1) if the tensors of type `i` (anti)commute among themselves

    Examples
    ========

    >>> from sympy.combinatorics.tensor_can import get_symmetric_group_sgs, gens_products
    >>> base, gens = get_symmetric_group_sgs(2)
    >>> gens_products((base,gens,[[],[]],0))
    (6, [0, 2], [[1, 0, 2, 3, 4, 5], [0, 1, 3, 2, 4, 5], [2, 3, 0, 1, 4, 5]])
    >>> gens_products((base,gens,[[1],[]],0))
    (6, [2], [[0, 1, 3, 2, 4, 5]])
    """
    res_size, res_base, res_gens = tensor_gens(*v[0])
    for i in range(1, len(v)):
        size, base, gens = tensor_gens(*v[i])
        res_base, res_gens = bsgs_direct_product(res_base, res_gens, base,
                                               gens, 1)
    res_size = len(res_gens[0])
    id_af = range(res_size)
    res_gens = [h for h in res_gens if h != id_af]
    if not res_gens:
        res_gens = [id_af]
    return res_size, res_base, res_gens
