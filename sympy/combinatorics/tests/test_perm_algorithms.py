from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.permutations import Permutation, perm_af_muln, cyclic,perm_af_invert
from sympy.combinatorics.perm_algorithms import double_coset_can_rep, canonicalization, tensor_gens, get_coset_repr, riemann_gens, gens_products, get_symmetric_group_sgs

"""
References:
[1] test_xperm.cc in cadabra by Kasper Peeters, http://cadabra.phi-sci.com/
"""

def test_double_coset_can_rep():
    sgens = [Permutation([2,1,0,3,4,5,7,6]), Permutation([4,1,2,3,0,5,7,6])]
    t = [([4,2,0,1,3,5,6,7], [0,1,2,3,4,5,7,6]),
         ([4,2,0,5,1,3,6,7], 0),
         ([4,2,0,1,5,3,6,7], 0),
         ([1,2,3,4,5,0,6,7], [0,1,2,3,4,5,6,7]),
         ([5,2,3,4,0,1,6,7], [0,1,2,3,4,5,7,6])]
    for g, r in t:
        assert double_coset_can_rep(0, sgens, g) == r

    sgens = [Permutation([2,1,0,3,4,5,6,7]), Permutation([4,1,2,3,0,5,6,7])]
    t = [([4,2,0,1,3,5,6,7], [0,1,2,3,4,5,6,7]),
         ([4,2,0,5,1,3,6,7], [0, 2, 1, 4, 5, 3, 6, 7]),
         ([4,2,0,1,5,3,6,7], [0, 2, 1, 4, 5, 3, 6, 7]),
         ([1,2,3,4,5,0,6,7], [0,1,2,3,4,5,6,7]),
         ([5,2,3,4,0,1,6,7], [0,1,2,3,4,5,6,7])]
    for g, r in t:
        assert double_coset_can_rep(0, sgens, g) == r

    # comparison with test_xperm.cc test2 in cadabra [1]
    # write the generators in cyclic form, to save space
    sgens = [[[25,26],[1,2]], [[25,26],[1,3]], [[25,26],[1,4]], [[25,26],[2,3]],
      [[25,26],[2,4]], [[25,26],[3,4]], [[25,26],[5,6]], [[25,26],[5,7]],
      [[25,26],[5,8]], [[25,26],[6,7]], [[25,26],[6,8]], [[25,26],[7,8]],
      [[25,26],[9,10]], [[25,26],[9,11]], [[25,26],[9,12]],
      [[25, 26], [10, 11]], [[25, 26], [10, 12]], [[25, 26], [11, 12]],
      [[25,26],[13,14]], [[25,26],[13,15]], [[25,26],[13,16]],
      [[25,26],[14,15]], [[25,26],[14,16]], [[25,26],[15,16]],
      [[25,26],[17,18]], [[25,26],[19,20]], [[18,20],[17,19]],
      [[25,26],[21,22]], [[25,26],[23,24]], [[22,24],[21,23]],
      [[4,8],[3,7],[2,6],[1,5]], [[4,12],[3,11],[2,10],[1,9]],
      [[4,16],[3,15],[2,14],[1,13]], [[8,12],[7,11],[6,10],[5,9]],
      [[8,16],[7,15],[6,14],[5,13]], [[12,16],[11,15],[10,14],[9,13]],
      [[20,24],[19,23],[18,22],[17,21]]]

    sgens = [Permutation(cyclic(x, 26)) for x in sgens]

    # perm element representing a tensor in the notation of test_xperm.c test2
    g0 = [1,23,5,13,8,16,11,17,9,18,12,20,10,24,2,21,6,14,3,22,7,15,4,19]
    g1 = [1,17,6,13,14,24,7,15,3,12,10,16,18,22,5,23,19,21,2,11,4,9,8,20]
    g2 = [24,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,1,2]
    g3 = [2,24,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,1,3]
    gv = [g0,g1,g2,g3]
    r0 = [1,5,2,6,3,7,4,8,9,17,10,18,11,19,12,21,13,20,14,22,15,23,16,24,25,26]
    r1 = [1,5,2,6,3,7,4,9,8,17,10,13,11,14,12,21,15,19,16,22,18,23,20,24,26,25]
    rv = [r0, r1, 0, 0]

    # conversion from xperm notation to our notation
    # convert to 0-base
    for i in range(4):
        g = [x-1 for x in gv[i]] + [24,25]
        # invert to have Portugal notation, giving slots from indices
        g = perm_af_invert(g)
        r = double_coset_can_rep(0, sgens, g)
        if r:
            # convert back
            r = perm_af_invert(r)
            r = [x+1 for x in r]
        assert r == rv[i]

def test_get_symmetric_group_sgs():
    assert get_symmetric_group_sgs(2) == [[1,0,2,3]]
    assert get_symmetric_group_sgs(2, 0) == [[1,0,3,2]]
    assert get_symmetric_group_sgs(3) == [[1,0,2,3,4],[0,2,1,3,4]]
    assert get_symmetric_group_sgs(3, 0) == [[1,0,2,4,3],[0,2,1,4,3]]
    assert get_symmetric_group_sgs(4) == [[1,0,2,3,4,5], [0,2,1,3,4,5], [0,1,3,2,4,5]]
    assert get_symmetric_group_sgs(4, 0) == [[1,0,2,3,5,4], [0,2,1,3,5,4], [0,1,3,2,5,4]]

def test_riemann_invariants():
    """
    Riemann invariant with 6 Riemann tensors
    R_{\mu_{11}}{}^{\mu_1}{}_{\mu_0}{}^{\mu_5}R^{\mu_6 \mu_4 \mu_0}{}_{\mu_5}R_{\mu_7 \mu_2 \mu_8 \mu_9}R_{\mu_{10} \mu_3 \mu_6 \mu_4}R^{\mu_2 \mu_7 \mu_{11}}{}_{\mu_1}R^{\mu_8 \mu_9 \mu_3 \mu_{10}}
    Replacing contravariant \mu_k with 2*k and
    covariant \mu_k with 2*k+1 one gets the corresponding permutation g
    """
    g = [23,2,1,10,12,8,0,11,15,5,17,19,21,7,13,9,4,14,22,3,16,18,6,20,24,25]
    sgens = tensor_gens(riemann_gens, 4, 6)
    sgs = get_coset_repr(sgens)
    sgens = [Permutation(x) for x in sgens]
    gcan = double_coset_can_rep(0, sgens, g, sgs)
    assert gcan == [0,2,4,6,1,3,8,10,5,7,12,14,9,11,16,18,13,15,20,22,17,19,21,23,24,25]
    """
    from gcan one gets the canonicalized tensor
    R^{\mu_0 \mu_1 \mu_2 \mu_3}R_{\mu_0 \mu_1}{}^{\mu_4 \mu_5}R_{\mu_2 \mu_3}{}^{\mu_6 \mu_7}R_{\mu_4 \mu_5}{}^{\mu_8 \mu_9}R_{\mu_6 \mu_7}{}^{\mu_{10} \mu_{11}}R_{\mu_8 \mu_9 \mu_{10} \mu_{11}}
    """

    sgens = tensor_gens(riemann_gens, 4, 10)
    sgs = get_coset_repr(sgens)
    sgens = [Permutation(x) for x in sgens]
    g = [0,2,5,7,4,6,9,11,8,10,13,15,12,14,17,19,16,18,21,23,20,22,25,27,24,26,29,31,28,30,33,35,32,34,37,39,36,38,1,3,40,41]
    gcan = double_coset_can_rep(0, sgens, g, sgs)
    assert gcan == [0,2,4,6,1,3,8,10,5,7,12,14,9,11,16,18,13,15,20,22,17,19,24,26,21,23,28,30,25,27,32,34,29,31,36,38,33,35,37,39,40,41]

    sgens = tensor_gens(riemann_gens, 4, 12)
    sgs = get_coset_repr(sgens)
    sgens = [Permutation(x) for x in sgens]
    g = [17, 44, 11, 3, 0, 19, 23, 15, 38, 4, 25, 27, 43, 36, 22, 14, 8, 30, 41, 20, 2, 10, 12, 28, 18, 1, 29, 13, 37, 42, 33, 7, 9, 31, 24, 26, 39, 5, 34, 47, 32, 6, 21, 40, 35, 46, 45, 16, 48, 49]
    gcan = double_coset_can_rep(0, sgens, g, sgs)
    assert gcan == [0, 2, 4, 6, 1, 3, 8, 10, 5, 7, 12, 14, 9, 11, 16, 18, 13, 15, 20, 22, 17, 19, 24, 26, 21, 23, 28, 30, 25, 27, 32, 34, 29, 31, 36, 38, 33, 35, 40, 42, 37, 39, 44, 46, 41, 43, 45, 47, 48, 49]

    sgens = tensor_gens(riemann_gens, 4, 20)
    sgs = get_coset_repr(sgens)
    sgens = [Permutation(x) for x in sgens]
    g = [0,2,4,6, 7,8,10,12, 14,16,18,20, 19,22,24,26, 5,21,28,30, 32,34,36,38, 40,42,44,46, 13,48,50,52, 15,49,54,56, 17,33,41,58, 9,23,60,62, 29,35,63,64, 3,45,66,68, 25,37,47,57, 11,31,69,70, 27,39,53,72, 1,59,73,74, 55,61,67,76, 43,65,75,78, 51,71,77,79, 80,81]
    gcan = double_coset_can_rep(0, sgens, g, sgs)
    assert gcan == [0,2,4,6, 1,8,10,12, 3,14,16,18, 5,20,22,24, 7,26,28,30, 9,15,32,34, 11,36,23,38, 13,40,42,44, 17,39,29,46, 19,48,43,50, 21,45,52,54, 25,56,33,58, 27,60,53,62, 31,51,64,66, 35,65,47,68, 37,70,49,72, 41,74,57,76, 55,67,59,78, 61,69,71,75, 63,79,73,77, 80,81]

def test_mixed_tensor_invariant():
    sym2_gens = get_symmetric_group_sgs(2)
    sgens = gens_products([riemann_gens, sym2_gens], [2, 3])
    sgs = get_coset_repr(sgens)
    sgens = [Permutation(x) for x in sgens]
    g = [12,10,5,2,8,0,4,6,13,1,7,3,9,11,14,15]
    gcan = double_coset_can_rep(0, sgens, g, sgs)
    assert gcan == [0, 2, 4, 6, 1, 8, 10, 12, 3, 9, 5, 11, 7, 13, 15, 14]


def test_canonicalization():
    sgens = [Permutation([1,0,2,3,5,4]), Permutation([2,3,0,1,4,5])]
    a = (([1,2,3,0],[3,0,1,2],[2,1,0,3],[0,3,2,1],[0,3,2,1,4,5]),
         ([1,3,0,2],[2,0,3,1],[0,2,1,3],[3,1,2,0],[0,2,1,3,4,5]),
         ([1,0,3,2],[2,3,0,1],[3,2,1,0],[0,1,2,3],[0,1,2,3,4,5]),
         ([3,2,0,1],[2,3,1,0],[1,0,2,3],[0,1,3,2],[0,1,2,3,5,4]),
         ([1,2,0,3],[3,0,2,1],[2,1,3,0],[0,3,1,2],[0,3,2,1,5,4]),
         ([2,0,1,3],[1,3,2,0],[3,1,0,2],[0,2,3,1],[0,2,1,3,5,4])
         )
    for t in a:
        for p in t[:-1]:
            assert canonicalization([], 0, sgens, p+[4,5]) == t[-1]

    sgens = [Permutation(x) for x in [[2,1,0,3,4,5,7,6], [4,1,2,3,0,5,7,6]] ]
    g = [5,4,2,3,0,1,6,7]
    dummies = [2,3,4,5]
    assert canonicalization(dummies, 0, sgens, g) == [0, 2, 3, 4, 5, 1, 6, 7]
    assert canonicalization(dummies, 1, sgens, g) == [0, 2, 3, 4, 5, 1, 7, 6]
