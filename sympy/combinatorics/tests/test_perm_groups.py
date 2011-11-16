from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.permutations import Permutation

def test_permgroups():
    a = Permutation([0,1])
    b = Permutation([1,0])
    c = PermutationGroup([a,b], 2)
    assert c.perm_size == 2
    assert c.generators == [Permutation([0, 1]), \
                            Permutation([1, 0])]
    assert c.is_abelian
    assert c.coset_repr(0) == set([Permutation([0, 1]), Permutation([1, 0])])
    assert c.coset_repr(1) == set([Permutation([0, 1])])
    assert list(c.generate()) == [Permutation([1, 0])]
    assert c.pointwise_stabilizers([0]) == []
    assert c.pointwise_stabilizers([1]) == []
    assert c.setwise_stabilizers([0], 1) == []
    assert c.setwise_stabilizers([1], 1) == [Permutation([1, 0])]
    assert c.center == [Permutation([1, 0])]
    assert c.commutators == [Permutation([0, 1])]

    a = Permutation([4, 1, 2, 3, 5, 0])
    b = Permutation([4, 1, 2, 3, 0, 5])
    c = Permutation([1, 4, 2, 5, 3, 0])
    d = Permutation([3, 1, 2, 4, 5, 0])
    e = Permutation([4, 2, 3, 1, 0, 5])
    f = Permutation([1, 2, 3, 4, 5, 0])
    g = Permutation([0, 1, 3, 2, 4, 5])
    h = Permutation([0, 2, 1, 3, 5, 4])
    i = PermutationGroup([a, b, c, d, e, f, g, h], 6)
    assert i.generators == [Permutation([4, 1, 2, 3, 5, 0]), Permutation([4, 1, 2, 3, 0, 5]), \
                            Permutation([1, 4, 2, 5, 3, 0]), Permutation([3, 1, 2, 4, 5, 0]), \
                            Permutation([4, 2, 3, 1, 0, 5]), Permutation([1, 2, 3, 4, 5, 0]), \
                            Permutation([0, 1, 3, 2, 4, 5]), Permutation([0, 2, 1, 3, 5, 4])]
    assert i.is_abelian == False
    assert i.coset_repr(3) == set([Permutation([1, 0, 2, 4, 3, 5]), Permutation([1, 4, 2, 5, 3, 0]), \
                                   Permutation([0, 1, 2, 3, 4, 5]), Permutation([1, 5, 2, 0, 3, 4]), \
                                   Permutation([0, 4, 3, 2, 5, 1]), Permutation([4, 0, 2, 1, 5, 3])])
    assert i.coset_repr(5) == set([Permutation([1, 0, 2, 4, 3, 5]), Permutation([0, 1, 2, 3, 4, 5]), \
                                   Permutation([0, 4, 3, 2, 5, 1]), Permutation([4, 0, 2, 1, 5, 3]), \
                                   Permutation([1, 5, 2, 0, 3, 4])])
    assert len(set(i.generate())) == 719


