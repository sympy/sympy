from sympy.combinatorics import Permutation
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.grouphomomorphism import homomorphism
from sympy.combinatorics.free_groups import free_group
from sympy.combinatorics.fp_groups import FpGroup
from sympy.combinatorics.named_groups import AlternatingGroup

def test_homomorphism():
    F, a, b = free_group("a, b")
    G = FpGroup(F, [a**3, b**3, (a*b)**2])

    c = Permutation(3)(0, 1, 2)
    d = Permutation(3)(1, 2, 3)
    A = AlternatingGroup(4)
    T = homomorphism(G, A, [a, b], [c, d])
    assert T(a*b**2*a**-1) == c*d**2*c**-1
    assert T.is_isomorphism()
    assert T(T.invert(Permutation(3)(0, 2, 3))) == Permutation(3)(0, 2, 3)

    T = homomorphism(G, AlternatingGroup(4), G.generators)
    assert T.is_trivial()
    assert T.kernel().order() == G.order()

    F, a = free_group("a")
    G = FpGroup(F, [a**8])
    P = PermutationGroup([Permutation(0, 1, 2, 3), Permutation(0, 2)])
    T = homomorphism(G, P, [a], [Permutation(0, 1, 2, 3)])
    assert T.image().order() == 4
    assert T(T.invert(Permutation(0, 2)(1, 3))) == Permutation(0, 2)(1, 3)
