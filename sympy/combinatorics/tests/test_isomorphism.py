from sympy.combinatorics import Permutation
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.homomorphisms import homomorphism
from sympy.combinatorics.free_groups import free_group
from sympy.combinatorics.fp_groups import FpGroup
from sympy.combinatorics.named_groups import AlternatingGroup, DihedralGroup
from sympy.combinatorics.named_groups import CyclicGroup
from sympy.combinatorics.isomorphisms import group_isomorphism, is_isomorphism

def test_isomorphisms():
    # Trivial Case
    # FpGroup -> FpGroup
    F, a, b = free_group("a, b")
    H = FpGroup(F, [a**3, b**3, (a*b)**2])
    F, c, d = free_group("c, d")
    G = FpGroup(F, [c**3, d**3, (c*d)**2])
    assert is_isomorphism(G, H)

    # FpGroup -> PermutationGroup
    # FpGroup is converted to the equivalent isomorphic group.
    F, a, b = free_group("a, b")
    G = FpGroup(F, [a**3, b**3, (a*b)**2])
    H = AlternatingGroup(4)
    T = group_isomorphism(G, H)
    assert T.kernel().order() == G.order()
    assert T.is_trivial()
    assert is_isomorphism(G, H)

    # PermutationGroup -> PermutationGroup
    D = DihedralGroup(8)
    p = Permutation(0, 1, 2, 3, 4, 5, 6, 7)
    P = PermutationGroup(p)
    assert not is_isomorphism(D, P)

    # Cyclic Groups of prime order are isomorphic to each other.
    A = CyclicGroup(5)
    B = CyclicGroup(7)
    assert is_isomorphism(A, B)
