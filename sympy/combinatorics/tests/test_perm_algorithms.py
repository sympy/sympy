from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.permutations import Permutation, perm_af_muln, cyclic,perm_af_invert
from sympy.combinatorics.perm_algorithms import double_coset_can_rep, canonicalization

"""
References:
[1] test_xperm.c in cadabra by Kasper Peeters, http://cadabra.phi-sci.com/
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
    # g0 leads to sagfault or hanging in test_xperm.cc
    # g1, g2, g3 give r1, 0, 0
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
