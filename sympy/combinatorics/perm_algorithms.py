from sympy.combinatorics.permutations import Permutation, perm_af_mul, perm_af_muln, perm_af_invert, _new_from_array_form, cyclic
from sympy.combinatorics.perm_groups import PermutationGroup

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

def orbit(genv, alpha):
    r"""
    Compute the orbit of alpha `\{g(\alpha) | g \in G\}` as a set.

    genv   generators of the Group G in array form
    alpha  starting point of the orbit

    The time complexity of the algorithm used here is `O(|Orb|*r)` where
    `|Orb|` is the size of the orbit and `r` is the number of generators
    of the group. For a proof of correctness, see [1], p.78.

    Examples
    ========

    >>> from sympy.combinatorics.perm_algorithms import orbit
    >>> from sympy.combinatorics.permutations import Permutation
    >>> genv = [[1,2,0,4,5,6,3]]
    >>> orbit(genv, 0)
    set([0, 1, 2])
    >>> orbit(genv, 4)
    set([3, 4, 5, 6])

    See Also
    ========

    orbit_transversal

    References
    ==========

    [1] Holt, D., Eick, B., O'Brien, E.
    "Handbook of computational group theory"

    """
    n = len(genv[0])
    orb = [alpha]
    used = [False]*n
    used[alpha] = True
    for b in orb:
        for gen in genv:
            temp = gen[b]
            if used[temp] == False:
                orb.append(temp)
                used[temp] = True
    return set(orb)

def orbit_transversal(genv, alpha, pairs=False, af=False):
    r"""
    Computes a transversal for the orbit of ``alpha`` as a set.

    genv   generators of the Group G
    alpha  starting point of the orbit
    pairs  = True return the list of pairs
           `(\beta, g_\beta)`. For a proof of correctness, see [1], p.79
    af     = True genv and the transversal elements are in array form

    For a permutation group `G`, a transversal for the orbit
    `Orb = \{g(\alpha) | g \in G\}` is a set
    `\{g_\beta | g_\beta(\alpha) = \beta\}` for `\beta \in Orb`.
    Note that there may be more than one possible transversal.

    Examples
    ========
    >>> from sympy.combinatorics.permutations import _new_from_array_form
    >>> from sympy.combinatorics.named_groups import DihedralGroup
    >>> from sympy.combinatorics.perm_algorithms import orbit_transversal
    >>> G = DihedralGroup(6)
    >>> orbit_transversal(G.generators, 0)
    [Permutation([0, 1, 2, 3, 4, 5]), Permutation([1, 2, 3, 4, 5, 0]),
    Permutation([5, 4, 3, 2, 1, 0]), Permutation([2, 3, 4, 5, 0, 1]),
    Permutation([4, 3, 2, 1, 0, 5]), Permutation([3, 4, 5, 0, 1, 2])]

    See Also
    ========

    orbit

    References
    ==========

    [1] Holt, D., Eick, B., O'Brien, E.
    "Handbook of computational group theory"

    """
    if af:
        n = len(genv[0])
    else:
        n = genv[0].size
    tr = [(alpha, range(n))]
    used = [False]*n
    used[alpha] = True
    if not af:
        genv = [h.array_form for h in genv]
    for pair in tr:
        for gen in genv:
            temp = gen[pair[0]]
            if used[temp] == False:
                tr.append((temp, perm_af_mul(gen, pair[1])))
                used[temp] = True
    if pairs:
        if af:
            return tr
        else:
            tr = [(v, _new_from_array_form(h)) for v, h in tr]
            return tr
    if af:
        return [pair[1] for pair in tr]
    else:
        return [_new_from_array_form(pair[1]) for pair in tr]

def dummy_sgs(d, sym, n):
    """
    Return the strong generators for dummy indices

    d   list of dummy indices
    sym symmetry under interchange of contracted dummies
    sym = None  n symmetry
    sym = 0 symmetric
          1 antisymmetric
    n   number of indices

    in base form the dummy indices are always in consecutive positions

    Examples
    ========
    >>> from sympy.combinatorics.perm_algorithms import dummy_sgs
    >>> dummy_sgs([2,3,4,5,6,7], 0, 8)
    [[0, 1, 3, 2, 4, 5, 6, 7, 8, 9], [0, 1, 2, 3, 5, 4, 6, 7, 8, 9], [0, 1, 2, 3, 4, 5, 7, 6, 8, 9], [0, 1, 4, 5, 2, 3, 6, 7, 8, 9], [0, 1, 2, 3, 6, 7, 4, 5, 8, 9]]
    """
    assert len(d) <= n
    res = []
    for j in d[::2]:
        a = range(n+2)
        if sym == 1:
            a[n] = n+1
            a[n+1] = n
        a[j], a[j+1] = a[j+1], a[j]
        res.append(a)
    for j in d[:-3:2]:
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

def double_coset_can_rep(sym, sgens, g, sgs=None):
    """
    Butler-Portugal algorithm for tensor canonicalization with dummy indices

      sym     symmetry of the contructing metric
              0     symmetric
              1     antisymmetric
              None  no symmetry (not implemented yet)

      sgens   generators of the group of slot symmetries
      g       permutation representing the tensor
      sgs     if it is not None, sgens is a strong generating set of S
              sgs is the tuple (S.coset_repr(), S.base)

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
    g=[4,2,0,1,3,5,6,7]     == [4,2,1,0,3,5,6,7] = perm_af_mul(d, g)
    which differs from perm_af_mul(g, d)

    The slot symmetry acts from the right
    s = [2,1,0,3,4,5,7,6]  exchanges slots 0 and 2 and changes sign
    T(d3,d2,d1,-d1,-d2,-d3) == -T(d1,d2,d3,-d1,-d2,-d3)
    g=[4,2,0,1,3,5,6,7]     == [0,2,4,1,3,5,7,6] = perm_af_mul(g, s)

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
    >>> from sympy.combinatorics.perm_algorithms import double_coset_can_rep
    >>> sgens = [Permutation(x) for x in [[2,1,0,3,4,5,7,6], [4,1,2,3,0,5,7,6]]]
    >>> g = [4,2,0,1,3,5,6,7]
    >>> double_coset_can_rep(0, sgens, g)
    [0, 1, 2, 3, 4, 5, 7, 6]

    T^{d3}_{d1 d2}^{d1}_{d3}^{d2}
    = -T^{d3}_{d1 d3}^{d1}_{d2}^{d2}   under slot symmetry -(2,4)
    = T_{d3 d1}^{d3}^{d1}_{d2}^{d2}    under slot symmetry -(0,2)
    = T^{d3}_{d1 d3}^{d1}_{d2}^{d2}    symmetric metric
    = 0
    >>> g = [4,1,3,0,5,2,6,7]
    >>> double_coset_can_rep(0, sgens, g)
    0
    """
    if sym is None:
        raise NotImplementedError
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
        # strong generating set for S
        sgens = [h.array_form for h in sgens]
        sgensx = sgens + S.stabilizers_gens()
        b_S = S.base
        S_cosets = S.coset_repr()
    else:
        # strong generating set for S
        sgensx = [h.array_form for h in sgens]
        S_cosets, b_S = sgs
    # strong generating set for D
    dsgsx = dummy_sgs(dumx, sym, num_dummies)
    ginv = perm_af_invert(g)
    # md[j] = min(D_p.orbit(j)
    idn = range(size)
    # TAB = list of entries (s, d, h) where h = d*g*s
    TAB = [(idn, idn, g)]
    for i in range(size - 2):
        b = base[i]
        testb = b in b_S and sgensx
        if testb:
            deltab = orbit(sgensx, b)
        else:
            deltab = set([b])
        # p1 = min(IMAGES) = min(Union D_p*h*deltab for h in TAB)
        md = _min_dummies(dumx, 0, dummies)
        p_i = min([min([md[h[x]] for x in deltab]) for s,d,h in TAB])
        Dxtrav = orbit_transversal(dsgsx, p_i, af=True) if dumx else None
        if Dxtrav:
            Dxtrav = [perm_af_invert(x) for x in Dxtrav]
        deltap = dumx if p_i in dumx else [p_i]
        TAB1 = []
        nTAB = len(TAB)
        while TAB:
            s, d, h = TAB.pop()
            if min([md[h[x]] for x in deltab]) != p_i:
                continue
            deltab1 = [x for x in deltab if md[h[x]] == p_i]
            # NEXT = s*deltab1 intersection (d*g)**-1*deltap
            dg = perm_af_mul(d, g)
            dginv = perm_af_invert(dg)
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
                    s1 = _trace_S(s, j, b, S_cosets)
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
        _dumx_remove(dumx, p_i)
        dsgsx = dummy_sgs(dumx, sym, num_dummies)
    return TAB[0][-1]

def canonical_free(G, p, num_free):
    """
    canonization of a tensor with respect to  free indices

    G   slot permutation group
    p   permutation representing the tensor
    num_free  number of free indices
    The indices must be ordered with first the free indices

    see explanation in double_coset_can_rep
    The algorithm is given in [2]

    Examples
    ========

    >>> from sympy.combinatorics.permutations import Permutation
    >>> from sympy.combinatorics.perm_groups import PermutationGroup
    >>> from sympy.combinatorics.perm_algorithms import canonical_free
    >>> sgens = [Permutation([1,0,2,3,5,4]), Permutation([2,3,0,1,4,5])]
    >>> p = Permutation([2, 1, 0, 3, 4, 5])
    >>> S = PermutationGroup(sgens)
    >>> canonical_free(S, p, 4)
    ([0, 3, 2, 1, 4, 5], [])
    """
    p = p.array_form
    K = [px.array_form for px in G.generators] + G.stabilizers_gens()
    K1 = []
    for x in K:
        if x not in K1:
            K1.append(x)
    K = K1
    cosets = G.coset_repr()
    base = list(G.base)
    m = len(base)
    for x in sorted(p[:-2]):
        if x not in base:
            base.append(x)
    lambd = p
    K1 = []
    for i in range(m):
        b = base[i]
        if b == num_free:
            return lambd, K

        delta = [x[b] for x in cosets[b]]
        delta1 = [lambd[x] for x in delta]
        delta1min = min(delta1, key=lambda x: base.index(x))

        k = delta1.index(delta1min)
        p = delta[k]
        for omega in cosets[b]:
            if omega[b] == p:
                break
        lambd = perm_af_mul(lambd, omega)
        K = [px for px in K if px[b] == b]
    return lambd, K

def canonicalization(dummies, sym, sgens, g):
    """
    canonicalization of a tensor

      dummies dummy indices

      sym     symmetry of the contructing metric
              0     symmetric
              1     antisymmetric
              None  no symmetry (not implemented yet)

      sgens   generators of the group of slot symmetries
      g       permutation representing the tensor


    Using the slot and dummy symmetries find the tensor with
    minimal indices according to free index < dummy index and
    then lexicographically

    Example:
    slot symmetry
    T(a0,a1,a2,a3,a4,a5) = -T(a2,a1,a0,a3,a4,a5);   -(0,2)
    T(a0,a1,a2,a3,a4,a5) = -T(a4,a1,a2,a3,a0,a5);   -(0,4)
    order of the indices:
    [a,b,d1,-d1,d2,-d2]
     0 1 2   3  4   5

    dummy indices: 2,3,4,5
    T(d2,-d2,d1,-d1,a,b)  corresponds to the permutation
    g = [4,5,2,3,0,1,6,7]

    Find the dummy group element `d` and the slot group element `s`
    such that g_can = d*g*s
    One gets g_can = [0, 2, 3, 4, 5, 1, 6, 7] corresponding to
    T(a,d1,-d1,d2,-d2,b)

    Let us verify it by hand:
    T(d2,-d2,d1,-d1,a,b) = T^{d2}_{d2}^{d1}_{d1}^{a b}
    symmetry -(0,4):
    T^{d2}_{d2}^{d1}_{d1}^{a b} = -T^{a}_{d2}^{d1}_{d1}^{d2 b}
    = T^{a}_{d2}^{d2}_{d1}^{d1 b}  by symmetry -(2,4)
    = T^{a}_{d1}^{d1}_{d2}^{d2 b}  by symmetry d1 <-> d2
    = T^{a d1}_{d1}^{d2}_{d2}^b    with (anti)symmetric metric

    For a more detailed expanation see double_coset_can_rep

    Examples
    ========

    >>> from sympy.combinatorics.permutations import Permutation
    >>> from sympy.combinatorics.perm_groups import PermutationGroup
    >>> from sympy.combinatorics.perm_algorithms import canonicalization
    >>> sgens = [Permutation([2,1,0,3,4,5,7,6]), Permutation([4,1,2,3,0,5,7,6])]
    >>> g = [4,5,2,3,0,1,6,7]
    >>> dummies = [2,3,4,5]
    >>> canonicalization(dummies, 0, sgens, g)
    [0, 2, 3, 4, 5, 1, 6, 7]

    """
    if sym is None:
        raise NotImplementedError
    size = len(g)
    num_free = size - 2 - len(dummies)
    if not num_free:
        return double_coset_can_rep(sym, sgens, g)
    S = PermutationGroup(sgens)
    if not dummies:
        return canonical_free(S, Permutation(g), num_free)[0]
    sgens = [x.array_form for x in sgens] + S.stabilizers_gens()
    g1, sgens1 = canonical_free(S, Permutation(g), num_free)
    sign = 1 if g1[-1] == size-1 else -1
    pos_dummies = [g1.index(x) for x in dummies]
    pos_free = [g1.index(x) for x in range(num_free)]
    # permutation representing the tensor corresponding to g1
    # without the free indices
    gd = [x - num_free for x in g1 if x in dummies]
    sgens1 = [[x-num_free for x in y][num_free:] for y in sgens1]
    size1 = size - num_free
    dummies = range(size1-2)
    sgens1 = [Permutation(x) for x in sgens1]
    g2 =  double_coset_can_rep(sym, sgens1, gd + [size1-2,size1-1])
    if g2 == 0:
        return 0
    if g2[-1] != len(g2) - 1:
        sign = -sign
    g3 = range(size)
    if sign == -1:
        g3[-2] = size - 1
        g3[-1] = size - 2
    for i in pos_free:
        g3[i] = g1[i]
    pos_dummies.sort()
    for c, i in enumerate(pos_dummies):
        g3[i] = g2[c] + num_free
    return g3

# SGS for the slot symmetry group of the Riemann tensor
riemann_gens = [[1,0,2,3,5,4], [0,1,3,2,5,4], [2,3,0,1,4,5]]

def get_symmetric_group_sgs(n, sym=True):
    """
    SGS for (anti)symmetric group with n elements

    Examples
    ========
    >>> from sympy.combinatorics.perm_algorithms import get_symmetric_group_sgs
    >>> get_symmetric_group_sgs(3)
    [[1, 0, 2, 3, 4], [0, 2, 1, 3, 4]]
    """
    gens = [Permutation(cyclic([(i,i+1)], n)).array_form for i in range(1, n)]
    if sym:
        gens = [x + [n, n+1] for x in gens]
    else:
        gens = [x + [n+1, n] for x in gens]
    return gens

def perm_af_direct_product(gens1, gens2, signed=False):
    """
    direct products of the generators gens1 and gens2

    Examples
    ========
    >>> from sympy.combinatorics.perm_algorithms import perm_af_direct_product
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

def get_coset_repr(sgens):
    """
    return (coset_repr, strong_base)

    sgens is a strong generating set for a permutation group of
    signed permutations
    """
    sgensx = sgens[:]
    n = len(sgens[0]) - 2
    b_S = []
    S_cosets = []
    # generate S_cosets, b_S
    for ii in range(n-1):
        if not sgensx:
            S_cosets.append([range(n+2)])
            continue
        Sxtrav = orbit_transversal(sgensx, ii, af=True)
        if len(Sxtrav) > 1:
            b_S.append(ii)
        S_cosets.append(Sxtrav)
        sgensx = [h for h in sgensx if h[ii] == ii]
    return S_cosets, b_S


def tensor_gens(basic_gens, num_indices, n):
    """
    Returns the strong generating set for n tensors of the same type

    n             number of tensors
    num_indices   number of indices for each tensor
    basic_gens    strong generating set for the slot symmetries
                  of a tensor
    The tensors are assumed to be commuting.

    The strong generating set for the invariant tensor made by the
    product of n tensors of the same type is composed of
      n*len(basic_gens) generators coming from basic_gens
      n-1   generators for the commutation of the tensors
    The size of each generator in n*sum_indices + 2
    each generator having the last two indices dedicated to the sign

    Examples
    ========
    >>> from sympy.combinatorics.perm_algorithms import tensor_gens, get_symmetric_group_sgs
    >>> sym3_gens = get_symmetric_group_sgs(3)
    >>> tensor_gens(sym3_gens, 3, 2)
    [[1, 0, 2, 3, 4, 5, 6, 7], [0, 2, 1, 3, 4, 5, 6, 7], [0, 1, 2, 4, 3, 5, 6, 7], [0, 1, 2, 3, 5, 4, 6, 7], [3, 4, 5, 0, 1, 2, 6, 7]]

    """
    res = basic_gens[:]
    ngens = len(basic_gens)
    size = num_indices + 2
    size2 = num_indices + 2
    for i in range(2, n+1):
        size1 = num_indices*(i-1) + 2
        res = perm_af_direct_product(res, basic_gens, 1)
        ngens += len(basic_gens)
        size += num_indices
        nr = size - 2
        a = list(range(nr-2*num_indices))
        for j in range(num_indices):
            a.append(nr-num_indices+j)
        for j in range(num_indices):
            a.append(nr-2*num_indices+j)
        a.append(nr)
        a.append(nr+1)
        res.append(a)
    return res

def gens_products(vgens, nvgens):
    """
    SGS slot generators for tensors of different types

    vgens  list of SGS of component tensors
    nvgens list of the number of component tensors of a given type
    """
    sgens = tensor_gens(vgens[0], len(vgens[0][0])-2, nvgens[0])
    for i in range(1, len(nvgens)):
        tmp = tensor_gens(vgens[i], len(vgens[i][0])-2, nvgens[i])
        sgens = perm_af_direct_product(sgens, tmp, 1)
    return sgens
