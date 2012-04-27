from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.permutations import Permutation, perm_af_muln, cyclic
from sympy.utilities.pytest import raises, skip, XFAIL
from sympy.combinatorics.generators import rubik_cube_generators
import random

def test_new():
    g = PermutationGroup([])
    assert g.order() == 0
    def test1():
        a = Permutation([2, 0, 1, 3, 4, 5])
        b = Permutation([0, 2, 1, 3, 4])
        g = PermutationGroup([a, b])
    raises(ValueError, 'test1()')
    

def test_generate():
    a = Permutation([1, 0])
    g = PermutationGroup([a]).generate()
    assert list(g) == [Permutation([0, 1]), Permutation([1, 0])]

    a = Permutation([2, 0, 1])
    b = Permutation([2, 1, 0])
    g = PermutationGroup([a, b]).generate()
    v1 = [p.array_form for p in list(g)]
    v1.sort()
    assert v1 == [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]

    a = Permutation([2, 0, 1, 3, 4, 5])
    b = Permutation([2, 1, 3, 4, 5, 0])
    g = PermutationGroup([a, b]).generate(af=True)
    assert len(list(g)) == 360

def test_order():
    a = Permutation([2,0,1,3,4,5,6,7,8,9])
    b = Permutation([2,1,3,4,5,6,7,8,9,0])
    g = PermutationGroup([a, b])
    assert g.order() == 1814400

def test_center():
    # TODO to be changed; center should return a group
    a = Permutation([1, 2, 0, 3, 5, 4])
    b = Permutation([2, 0, 1, 4, 3, 5])
    g = PermutationGroup([a, b])
    assert g.order() == 18
    v = [x.array_form for x in g.center]
    v.sort()
    assert v == [[0, 1, 2, 3, 4, 5], [1, 2, 0, 3, 4, 5], [2, 0, 1, 3, 4, 5]]

    N = 12
    a = range(N)
    v = []
    for i in range(0, N, 2):
        b = a[:]
        b[i+1], b[i] = b[i], b[i+1]
        v.append(Permutation(b))
    g = PermutationGroup(v)
    c = g.center
    assert len(c) == 2**(N//2)


def test_stabilizers():
    a = Permutation([2,0,1,3,4,5])
    b = Permutation([2,1,3,4,5,0])
    G = PermutationGroup([a,b])
    stab = G.pointwise_stabilizers([0])
    assert len(stab) == 60
    G0 = G.stabilizer_group(0)
    assert G0.order() == 60

    a = Permutation([[0], [1], [2,3], [4]])
    b = Permutation([[0, 2], [1], [3], [4]])
    g = PermutationGroup([a, b])
    s = g.pointwise_stabilizers([0])
    assert s == [Permutation(range(5)), Permutation([0, 1, 3, 2, 4])]
    s = g.pointwise_stabilizers([2])
    assert s == [Permutation(range(5)), Permutation([3, 1, 2, 0, 4])]

    gens_cube = [[1, 3, 5, 7, 0, 2, 4, 6], [1, 3, 0, 2, 5, 7, 4, 6]]
    gens = [Permutation(p) for p in gens_cube]
    G = PermutationGroup(gens)
    G2 = G.stabilizer_group(2)
    assert G2.order() == 6
    G2_1 = G2.stabilizer_group(1)
    v = list(G2_1.generate(af=True))
    s = G.pointwise_stabilizers([2, 1], af=True)
    assert s == v

    gens = ((1,2,0,4,5,3,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
            (0,1,2,3,4,5,19,6,8,9,10,11,12,13,14,15,16,7,17,18),
            (0,1,2,3,4,5,6,7,9,18,16,11,12,13,14,15,8,17,10,19))
    gens = [Permutation(p) for p in gens]
    G = PermutationGroup(gens)
    G2 = G.stabilizer_group(2)
    assert G2.order() == 181440

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

def test_random():
    random.seed(10)
    gens_cube = [[1, 3, 5, 7, 0, 2, 4, 6], [1, 3, 0, 2, 5, 7, 4, 6]]
    gens = [Permutation(p) for p in gens_cube]
    G = PermutationGroup(gens)
    h = G.random(af=True)
    assert h == [4, 5, 6, 7, 0, 1, 2, 3]

def test_has():
    a = Permutation([2,0,1,3,4,5])
    b = Permutation([2,1,3,4,5,0])
    g = PermutationGroup([a, b])
    assert g.order() == 360
    rep = g.coset_repr()
    d = Permutation([1,0,2,3,4,5])
    v1 = g.coset_decomposition(d.array_form)
    assert not v1
    v = g.coset_decomposition([1,0,2,3,5,4])
    assert perm_af_muln(*v) == [1,0,2,3,5,4]

def test_normal_closure():
    a = Permutation([2,0,1,3,4,5,6])
    b = Permutation([2,1,3,4,5,6,0])
    g = PermutationGroup([a, b])
    c = a*b*~a*~b
    g2 = PermutationGroup([c], a.size)
    assert g.order() == 5040
    assert g2.order() == 2
    g3 = g.normal_closure(g2)
    assert g3.order() == 2520

def test_commutators():
    a = Permutation([1, 0, 2, 4, 3])
    b = Permutation([0, 1, 3, 2, 4])
    g = PermutationGroup([a,b])
    assert g.order() == 12
    cm = g.commutators()
    assert cm.order() == 3
    cm1 = g.commutators(method="comm")
    assert len(cm1) == 3

    # http://math.stackexchange.com/questions/7811/derived-subgroups-and-commutators
    # lowest order permutation group for which the commutator group
    # contain elements which are not commutators
    a = [[(3,8,6),(4,7,5),(9,27,17),(10,28,18),(11,30,22),(12,29,21),(13,26,23),
     (14,25,24),(15,31,20),(16,32,19)],
    [(1,17,7,23),(2,18,8,24),(3,19,5,21),(4,20,6,22),(9,26,15,32),
     (10,25,16,31),(11,28,13,30),(12,27,14,29)],
    [(1,9,5,13),(2,10,6,14),(3,11,7,15),(4,12,8,16),(17,25,21,29),(18,26,22,30),
     (19,27,23,31),(20,28,24,32)],
    [(1,5),(2,6),(3,7),(4,8),(9,13),(10,14),(11,15),(12,16),(17,21),(18,22),
     (19,23),(20,24),(25,29),(26,30),(27,31),(28,32)],
    [(1,3),(2,4),(5,7),(6,8),(9,11),(10,12),(13,15),(14,16),(17,19),
     (18,20),(21,23),(22,24),(25,27),(26,28),(29,31),(30,32)],
    [(1,2),(3,4),(5,6),(7,8),(9,10),(11,12),(13,14),(15,16),(17,18),(19,20),
     (21,22),(23,24),(25,26),(27,28),(29,30),(31,32)]]
    va = [Permutation(cyclic(x, 32)) for x in a]
    g = PermutationGroup(va)
    assert g.order() == 96
    gcv = g.commutators(af=True)
    assert gcv.order() == 32
    gcv1 = g.commutators(af=True, method="comm")
    assert len(gcv1) == 29

def test_orbits():

    a = Permutation([2, 0, 1])
    b = Permutation([2, 1, 0])
    g = PermutationGroup([a,b])
    assert g.orbit(0) == set([0, 1, 2])
    assert g.orbits() == [set([0, 1, 2])]
    assert g.orbits(rep=True) == [0]

    gens = rubik_cube_generators()
    g = PermutationGroup(gens, 48)
    assert g.orbits(rep=True) == [0, 1]

@XFAIL
def test_rubik():
    skip('takes too much time')
    gens = rubik_cube_generators()
    G = PermutationGroup(gens)
    assert G.order() == 43252003274489856000


