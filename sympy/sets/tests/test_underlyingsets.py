from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.permutations import Permutation
from sympy.sets.sets import FiniteSet
from sympy.sets.underlyingsets import UnderlyingSetOf


def test_PermutationGroup():
    a = Permutation([2, 0, 1])
    b = Permutation([2, 1, 0])
    G = PermutationGroup([a, b])

    g = G.underlying_set()
    assert g == UnderlyingSetOf(G)

    assert a in g
    assert b in g
    assert a*b in g

    g = g.rewrite(FiniteSet)
    expected_members = [
        Permutation([0, 1, 2]),
        Permutation([0, 2, 1]),
        Permutation([1, 0, 2]),
        Permutation([1, 2, 0]),
        Permutation([2, 0, 1]),
        Permutation([2, 1, 0]),
    ]
    assert g == FiniteSet(*expected_members)
