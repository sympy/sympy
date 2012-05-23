from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.permutations import Permutation, perm_af_muln, cyclic
from sympy.utilities.pytest import raises, skip, XFAIL
from sympy.combinatorics.generators import rubik_cube_generators
import random

def test_new():
    a = Permutation([1, 0])
    G =  PermutationGroup([a])
    assert G.is_abelian
    a = Permutation([2, 0, 1])
    b = Permutation([2, 1, 0])
    G =  PermutationGroup([a, b])
    assert not G.is_abelian

    def test1():
        a = Permutation([2, 0, 1, 3, 4, 5])
        b = Permutation([0, 2, 1, 3, 4])
        g = PermutationGroup([a, b])
    raises(ValueError, 'test1()')

def test_generate():
    a = Permutation([1, 0])
    g = PermutationGroup([a]).generate()
    assert list(g) == [Permutation([0, 1]), Permutation([1, 0])]
    g = PermutationGroup([a]).generate(method='dimino')
    assert list(g) == [Permutation([0, 1]), Permutation([1, 0])]

    a = Permutation([2, 0, 1])
    b = Permutation([2, 1, 0])
    G = PermutationGroup([a, b])
    g = G.generate()
    v1 = [p.array_form for p in list(g)]
    v1.sort()
    assert v1 == [[0,1,2], [0,2,1], [1,0,2], [1,2,0], [2,0,1], [2,1,0]]

    v2 = list(G.generate(method='dimino', af=True))
    assert v1 == sorted(v2)

    a = Permutation([2, 0, 1, 3, 4, 5])
    b = Permutation([2, 1, 3, 4, 5, 0])
    g = PermutationGroup([a, b]).generate(af=True)
    assert len(list(g)) == 360

def test_order():
    a = Permutation([2,0,1,3,4,5,6,7,8,9])
    b = Permutation([2,1,3,4,5,6,7,8,9,0])
    g = PermutationGroup([a, b])
    assert g.order() == 1814400

def test_stabilizer():
    a = Permutation([2,0,1,3,4,5])
    b = Permutation([2,1,3,4,5,0])
    G = PermutationGroup([a,b])
    G0 = G.stabilizer(0)
    assert G0.order() == 60

    gens_cube = [[1, 3, 5, 7, 0, 2, 4, 6], [1, 3, 0, 2, 5, 7, 4, 6]]
    gens = [Permutation(p) for p in gens_cube]
    G = PermutationGroup(gens)
    G2 = G.stabilizer(2)
    assert G2.order() == 6
    G2_1 = G2.stabilizer(1)
    v = list(G2_1.generate(af=True))
    assert v == [[0, 1, 2, 3, 4, 5, 6, 7], [3, 1, 2, 0, 7, 5, 6, 4]]

    gens = ((1,2,0,4,5,3,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
            (0,1,2,3,4,5,19,6,8,9,10,11,12,13,14,15,16,7,17,18),
            (0,1,2,3,4,5,6,7,9,18,16,11,12,13,14,15,8,17,10,19))
    gens = [Permutation(p) for p in gens]
    G = PermutationGroup(gens)
    G2 = G.stabilizer(2)
    assert G2.order() == 181440

def test_coset_repr():
    a = Permutation([0, 2, 1])
    b = Permutation([1, 0, 2])
    G = PermutationGroup([a, b])
    assert G.coset_repr() == [[[0,1,2], [1,0,2], [2,0,1]], [[0,1,2], [0,2,1]]]
    assert G.stabilizers_gens() == [[0, 2, 1]]

def test_coset_rank():
    gens_cube = [[1, 3, 5, 7, 0, 2, 4, 6], [1, 3, 0, 2, 5, 7, 4, 6]]
    gens = [Permutation(p) for p in gens_cube]
    G = PermutationGroup(gens)
    i = 0
    for h in G.generate(af=True):
        rk = G.coset_rank(h)
        assert rk == i
        h1 = G.coset_unrank(rk, af=True)
        assert h == h1
        i += 1
    assert G.coset_unrank(48) == None
    assert G.coset_rank(gens[0]) == 6
    assert G.coset_unrank(6) == gens[0]

def test_coset_decomposition():
    a = Permutation([2,0,1,3,4,5])
    b = Permutation([2,1,3,4,5,0])
    g = PermutationGroup([a, b])
    assert g.order() == 360
    rep = g.coset_repr()
    d = Permutation([1,0,2,3,4,5])
    assert not g.coset_decomposition(d.array_form)
    assert not g.has_element(d)
    c = Permutation([1,0,2,3,5,4])
    v = g.coset_decomposition(c)
    assert perm_af_muln(*v) == [1,0,2,3,5,4]
    assert g.has_element(c)

    a = Permutation([0,2,1])
    g = PermutationGroup([a])
    c = Permutation([2,1,0])
    assert not g.coset_decomposition(c)
    assert g.coset_rank(c) == None

def test_orbits():
    a = Permutation([2, 0, 1])
    b = Permutation([2, 1, 0])
    g = PermutationGroup([a, b])
    assert g.orbit(0) == set([0, 1, 2])
    assert g.orbits() == [set([0, 1, 2])]
    assert g.is_transitive()
    assert g.orbits(rep=True) == [0]
    assert g.orbit_traversal(0) == \
        [Permutation([0,1,2]), Permutation([1,2,0]), Permutation([2,0,1])]
    orbt = g.orbit_traversal(1, af=True)
    assert g.orbit_traversal(1, af=True) == [[2, 0, 1], [0, 1, 2], [1, 2, 0]]

    a = Permutation(range(1, 100) + [0])
    G = PermutationGroup([a])
    assert G.orbits(rep=True) == [0]
    gens = rubik_cube_generators()
    g = PermutationGroup(gens, 48)
    assert g.orbits(rep=True) == [0, 1]
    assert not g.is_transitive()

def test_is_normal():
    gens_s5 = [Permutation(p) for p in [[1,2,3,4,0], [2,1,4,0,3]]]
    G1 = PermutationGroup(gens_s5)
    assert G1.order() == 120
    gens_a5 = [Permutation(p) for p in [[1,0,3,2,4], [2,1,4,3,0]]]
    G2 = PermutationGroup(gens_a5)
    assert G2.order() == 60
    assert G2.is_normal(G1)
    gens3 = [Permutation(p) for p in [[2,1,3,0,4], [1,2,0,3,4]]]
    G3 = PermutationGroup(gens3)
    assert not G3.is_normal(G1)
    assert G3.order() == 12
    G4 = G1.normal_closure(G3.generators)
    assert G4.order() == 60
    gens5 = [Permutation(p) for p in [[1,2,3,0,4], [1,2,0,3,4]]]
    G5 = PermutationGroup(gens5)
    assert G5.order() == 24
    G6 = G1.normal_closure(G5.generators)
    assert G6.order() == 120
    assert G1 == G6
    assert G1 != G4
    assert G2 == G4

def test_eq():
    a = [[1,2,0,3,4,5], [1,0,2,3,4,5], [2,1,0,3,4,5], [1,2,0,3,4,5]]
    a = [Permutation(p) for p in a + [[1,2,3,4,5,0]]]
    g = Permutation([1,2,3,4,5,0])
    G1, G2, G3 = [PermutationGroup(x) for x in [a[:2],a[2:4],[g, g**2]]]
    assert G1.order() == G2.order() == G3.order() == 6
    assert G1 == G2
    assert G1 != G3
    G4 = PermutationGroup([Permutation([0,1])])
    assert G1 != G4
    assert not G4.is_subgroup(G1)

def test_commutators():
    a = Permutation([1, 0, 2, 4, 3])
    b = Permutation([0, 1, 3, 2, 4])
    G = PermutationGroup([a,b])
    C = G.commutator()
    assert C.order() == 3
    assert C.is_normal(G)
    assert C.is_subgroup(G)
    assert not G.is_subgroup(C)
    gens_cube = [[1, 3, 5, 7, 0, 2, 4, 6], [1, 3, 0, 2, 5, 7, 4, 6]]
    gens = [Permutation(p) for p in gens_cube]
    G = PermutationGroup(gens)
    C = G.commutator()
    assert C.order() == 12

def test_is_solvable():
    a = Permutation([1,2,0])
    b = Permutation([1,0,2])
    G = PermutationGroup([a, b])
    assert G.is_solvable()
    a = Permutation([1,2,3,4,0])
    b = Permutation([1,0,2,3,4])
    G = PermutationGroup([a, b])
    assert not G.is_solvable()

def test_rubik1():
    gens = rubik_cube_generators()
    gens1 = [gens[0]] + [p**2 for p in gens[1:]]
    G1 = PermutationGroup(gens1)
    assert G1.order() == 19508428800
    gens2 = [p**2 for p in gens]
    G2 = PermutationGroup(gens2)
    assert G2.order() == 663552
    assert G2.is_subgroup(G1)
    C1 = G1.commutator()
    assert C1.order() == 4877107200
    assert C1.is_subgroup(G1)
    assert not G2.is_subgroup(C1)

@XFAIL
def test_rubik():
    skip('takes too much time')
    gens = rubik_cube_generators()
    G = PermutationGroup(gens)
    assert G.order() == 43252003274489856000
    G1 = PermutationGroup(gens[:3])
    assert G1.order() == 170659735142400
    assert not G1.is_normal(G)
    G2 = G.normal_closure(G1.generators)
    assert G2 == G
