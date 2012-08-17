from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics.group_constructs import DirectProduct
from sympy.combinatorics.named_groups import SymmetricGroup, CyclicGroup,\
DihedralGroup, AlternatingGroup, AbelianGroup
from sympy.combinatorics.permutations import Permutation, perm_af_muln, cyclic
from sympy.utilities.pytest import raises, skip, XFAIL
from sympy.combinatorics.generators import rubik_cube_generators
import random
from sympy.combinatorics.testutil import _verify_bsgs, _verify_centralizer,\
_cmp_perm_lists, _verify_normal_closure
from sympy.combinatorics.util import _distribute_gens_by_base

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
    raises(ValueError, lambda: test1())

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

def test_center():
    # the center of the dihedral group D_n is of order 2 for even n
    for i in (4, 6, 10):
        D = DihedralGroup(i)
        assert (D.center()).order() == 2
    # the center of the dihedral group D_n is of order 1 for odd n>2
    for i in (3, 5, 7):
        D = DihedralGroup(i)
        assert (D.center()).order() == 1
    # the center of an abelian group is the group itself
    for i in (2, 3, 5):
        for j in (1, 5, 7):
            for k in (1, 1, 11):
                G = AbelianGroup(i, j, k)
                assert G.center() == G
    # the center of a nonabelian simple group is trivial
    for i in(1, 5, 9):
        A = AlternatingGroup(i)
        assert (A.center()).order() == 1
    # brute-force verifications
    D = DihedralGroup(5)
    A = AlternatingGroup(3)
    C = CyclicGroup(4)
    G = D*A*C
    assert _verify_centralizer(G, G)

def test_centralizer():
    # the centralizer of the trivial group is the entire group
    S = SymmetricGroup(2)
    assert S.centralizer(Permutation(range(2))) == S
    A = AlternatingGroup(5)
    assert A.centralizer(Permutation(range(5))) == A
    # a centralizer in the trivial group is the trivial group itself
    triv = PermutationGroup([Permutation([0,1,2,3])])
    D = DihedralGroup(4)
    assert triv.centralizer(D) == triv
    # brute-force verifications for centralizers of groups
    for i in (4, 5, 6):
        S = SymmetricGroup(i)
        A = AlternatingGroup(i)
        C = CyclicGroup(i)
        D = DihedralGroup(i)
        for gp in (S, A, C, D):
            for gp2 in (S, A, C, D):
                if gp2 != gp:
                    assert _verify_centralizer(gp, gp2)
    # verify the centralizer for all elements of several groups
    S = SymmetricGroup(5)
    elements = list(S.generate_dimino())
    for element in elements:
        assert _verify_centralizer(S, element)
    A = AlternatingGroup(5)
    elements = list(A.generate_dimino())
    for element in elements:
        assert _verify_centralizer(A, element)
    D = DihedralGroup(7)
    elements = list(D.generate_dimino())
    for element in elements:
        assert _verify_centralizer(D, element)
    # verify centralizers of small groups within small groups
    small = []
    for i in (1, 2, 3):
        small.append(SymmetricGroup(i))
        small.append(AlternatingGroup(i))
        small.append(DihedralGroup(i))
        small.append(CyclicGroup(i))
    for gp in small:
        for gp2 in small:
            if gp.degree == gp2.degree:
                assert _verify_centralizer(gp, gp2)

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
    assert g.is_transitive
    assert g.orbits(rep=True) == [0]
    assert g.orbit_transversal(0) == \
        [Permutation([0, 1, 2]), Permutation([2, 0, 1]), Permutation([1, 2, 0])]
    assert g.orbit_transversal(0, True) == \
        [(0, Permutation([0, 1, 2])), (2, Permutation([2, 0, 1])), \
        (1, Permutation([1, 2, 0]))]

    a = Permutation(range(1, 100) + [0])
    G = PermutationGroup([a])
    assert G.orbits(rep=True) == [0]
    gens = rubik_cube_generators()
    g = PermutationGroup(gens, 48)
    assert g.orbits(rep=True) == [0, 1]
    assert not g.is_transitive

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

def test_derived_subgroup():
    a = Permutation([1, 0, 2, 4, 3])
    b = Permutation([0, 1, 3, 2, 4])
    G = PermutationGroup([a,b])
    C = G.derived_subgroup()
    assert C.order() == 3
    assert C.is_normal(G)
    assert C.is_subgroup(G)
    assert not G.is_subgroup(C)
    gens_cube = [[1, 3, 5, 7, 0, 2, 4, 6], [1, 3, 0, 2, 5, 7, 4, 6]]
    gens = [Permutation(p) for p in gens_cube]
    G = PermutationGroup(gens)
    C = G.derived_subgroup()
    assert C.order() == 12

def test_is_solvable():
    a = Permutation([1,2,0])
    b = Permutation([1,0,2])
    G = PermutationGroup([a, b])
    assert G.is_solvable
    a = Permutation([1,2,3,4,0])
    b = Permutation([1,0,2,3,4])
    G = PermutationGroup([a, b])
    assert not G.is_solvable

def test_rubik1():
    gens = rubik_cube_generators()
    gens1 = [gens[0]] + [p**2 for p in gens[1:]]
    G1 = PermutationGroup(gens1)
    assert G1.order() == 19508428800
    gens2 = [p**2 for p in gens]
    G2 = PermutationGroup(gens2)
    assert G2.order() == 663552
    assert G2.is_subgroup(G1)
    C1 = G1.derived_subgroup()
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

def test_direct_product():
    C = CyclicGroup(4)
    D = DihedralGroup(4)
    G = C*C*C
    assert G.order() == 64
    assert G.degree == 12
    assert len(G.orbits()) == 3
    assert G.is_abelian == True
    H = D*C
    assert H.order() == 32
    assert H.is_abelian == False

def test_orbit_rep():
    G = DihedralGroup(6)
    assert G.orbit_rep(1,3) in [Permutation([2, 3, 4, 5, 0, 1]),\
    Permutation([4, 3, 2, 1, 0, 5])]
    H = CyclicGroup(4)*G
    assert H.orbit_rep(1, 5) == False

def test_schreier_vector():
    G = CyclicGroup(50)
    v = [0]*50
    v[23] = -1
    assert G.schreier_vector(23) == v
    H = DihedralGroup(8)
    assert H.schreier_vector(2) == [0, 1, -1, 0, 0, 1, 0, 0]
    L = SymmetricGroup(4)
    assert L.schreier_vector(1) == [1, -1, 0, 0]

def test_random_pr():
    D = DihedralGroup(6)
    r = 11
    n = 3
    _random_prec_n = {}
    _random_prec_n[0] = {'s': 7, 't': 3, 'x': 2, 'e': -1}
    _random_prec_n[1] = {'s': 5, 't': 5, 'x': 1, 'e': -1}
    _random_prec_n[2] = {'s': 3, 't': 4, 'x': 2, 'e': 1}
    D._random_pr_init(r, n, _random_prec_n = _random_prec_n)
    assert D._random_gens[11] == Permutation([0, 1, 2, 3, 4, 5])
    _random_prec = {'s': 2, 't': 9, 'x': 1, 'e': -1}
    assert D.random_pr(_random_prec = _random_prec) == \
    Permutation([0, 5, 4, 3, 2, 1])

def test_is_alt_sym():
    G = DihedralGroup(10)
    assert G.is_alt_sym() == False
    S = SymmetricGroup(10)
    N_eps = 10
    _random_prec = {'N_eps': N_eps,
    0: Permutation([[2], [1, 4], [0, 6, 7, 8, 9, 3, 5]]),
    1: Permutation([[1, 8, 7, 6, 3, 5, 2, 9], [0, 4]]),
    2: Permutation([[5, 8], [4, 7], [0, 1, 2, 3, 6, 9]]),
    3: Permutation([[3], [0, 8, 2, 7, 4, 1, 6, 9, 5]]),
    4: Permutation([[8], [4, 7, 9], [3, 6], [0, 5, 1, 2]]),
    5: Permutation([[6], [0, 2, 4, 5, 1, 8, 3, 9, 7]]),
    6: Permutation([[6, 9, 8], [4, 5], [1, 3, 7], [0, 2]]),
    7: Permutation([[4], [0, 2, 9, 1, 3, 8, 6, 5, 7]]),
    8: Permutation([[1, 5, 6, 3], [0, 2, 7, 8, 4, 9]]),
    9: Permutation([[8], [6, 7], [2, 3, 4, 5], [0, 1, 9]])}
    assert S.is_alt_sym(_random_prec = _random_prec) == True
    A = AlternatingGroup(10)
    _random_prec = {'N_eps': N_eps,
    0: Permutation([[1, 6, 4, 2, 7, 8, 5, 9, 3], [0]]),
    1: Permutation([[1], [0, 5, 8, 4, 9, 2, 3, 6, 7]]),
    2: Permutation([[1, 9, 8, 3, 2, 5], [0, 6, 7, 4]]),
    3: Permutation([[6, 8, 9], [4, 5], [1, 3, 7, 2], [0]]),
    4: Permutation([[8], [5], [4], [2, 6, 9, 3], [1], [0, 7]]),
    5: Permutation([[3, 6], [0, 8, 1, 7, 5, 9, 4, 2]]),
    6: Permutation([[5], [2, 9], [1, 8, 3], [0, 4, 7, 6]]),
    7: Permutation([[1, 8, 4, 7, 2, 3], [0, 6, 9, 5]]),
    8: Permutation([[5, 8, 7], [3], [1, 4, 2, 6], [0, 9]]),
    9: Permutation([[4, 9, 6], [3, 8], [1, 2], [0, 5, 7]])}
    assert A.is_alt_sym(_random_prec = _random_prec) == False

def test_minimal_block():
    D = DihedralGroup(6)
    block_system = D.minimal_block([0,3])
    for i in range(3):
        assert block_system[i] == block_system[i+3]
    S = SymmetricGroup(6)
    assert S.minimal_block([0, 1]) == [0, 0, 0, 0, 0, 0]

def test_max_div():
    S = SymmetricGroup(10)
    assert S.max_div == 5

def test_is_primitive():
    S = SymmetricGroup(5)
    assert S.is_primitive() == True
    C = CyclicGroup(7)
    assert C.is_primitive() == True

def test_random_stab():
    S = SymmetricGroup(5)
    _random_el = Permutation([1, 3, 2, 0, 4])
    _random_prec = {'rand': _random_el}
    g = S.random_stab(2, _random_prec = _random_prec)
    assert g == Permutation([1, 3, 2, 0, 4])
    h = S.random_stab(1)
    assert h(1) == 1

def test_transitivity_degree():
    perm = Permutation([1, 2, 0])
    C = PermutationGroup([perm])
    assert C.transitivity_degree == 1
    gen1 = Permutation([1, 2, 0, 3, 4])
    gen2 = Permutation([1, 2, 3, 4, 0])
    # alternating group of degree 5
    Alt = PermutationGroup([gen1, gen2])
    assert Alt.transitivity_degree == 3

def test_schreier_sims_random():
    S = SymmetricGroup(3)
    base = [0, 1]
    strong_gens = [Permutation([1, 2, 0]), Permutation([1, 0, 2]),\
                  Permutation([0, 2, 1])]
    assert S.schreier_sims_random(base, strong_gens, 5) == (base, strong_gens)
    D = DihedralGroup(3)
    _random_prec = {'g': [Permutation([2, 0, 1]), Permutation([1, 2, 0]),\
                         Permutation([1, 0, 2])]}
    base = [0, 1]
    strong_gens = [Permutation([1, 2, 0]), Permutation([2, 1, 0]),\
                  Permutation([0, 2, 1])]
    assert D.schreier_sims_random([], D.generators, 2,\
           _random_prec=_random_prec) == (base, strong_gens)

def test_baseswap():
    S = SymmetricGroup(4)
    S.schreier_sims()
    base = S.base
    strong_gens = S.strong_gens
    assert base == [0, 1, 2]
    deterministic = S.baseswap(base, strong_gens, 1, randomized=False)
    randomized = S.baseswap(base, strong_gens, 1)
    assert deterministic[0] == [0, 2, 1]
    assert _verify_bsgs(S, deterministic[0], deterministic[1]) == True
    assert randomized[0] == [0, 2, 1]
    assert _verify_bsgs(S, randomized[0], randomized[1]) == True

def test_schreier_sims_incremental():
    identity = Permutation([0, 1, 2, 3, 4])
    TrivialGroup = PermutationGroup([identity])
    base, strong_gens = TrivialGroup.schreier_sims_incremental(base=[0, 1, 2])
    assert _verify_bsgs(TrivialGroup, base, strong_gens) == True
    S = SymmetricGroup(5)
    base, strong_gens = S.schreier_sims_incremental(base=[0,1,2])
    assert _verify_bsgs(S, base, strong_gens) == True
    D = DihedralGroup(2)
    base, strong_gens = D.schreier_sims_incremental(base=[1])
    assert _verify_bsgs(D, base, strong_gens) == True
    A = AlternatingGroup(7)
    gens = A.generators[:]
    gen0 = gens[0]
    gen1 = gens[1]
    gen1 = gen1*(~gen0)
    gen0 = gen0*gen1
    gen1 = gen0*gen1
    base, strong_gens = A.schreier_sims_incremental(base=[0,1], gens=gens)
    assert _verify_bsgs(A, base, strong_gens) == True
    C = CyclicGroup(11)
    gen = C.generators[0]
    base, strong_gens = C.schreier_sims_incremental(gens=[gen**3])
    assert _verify_bsgs(C, base, strong_gens) == True

def test_subgroup_search():
    prop_true = lambda x: True
    prop_fix_points = lambda x: [x(point) for point in points] == points
    prop_comm_g = lambda x: x*g == g*x
    prop_even = lambda x: x.is_even
    for i in range(10, 17, 2):
        S = SymmetricGroup(i)
        A = AlternatingGroup(i)
        C = CyclicGroup(i)
        Sym = S.subgroup_search(prop_true)
        assert Sym == S
        Alt = S.subgroup_search(prop_even)
        assert Alt == A
        Sym = S.subgroup_search(prop_true, init_subgroup=C)
        assert Sym == S
        points = [7]
        assert S.stabilizer(7) == S.subgroup_search(prop_fix_points)
        points = [3, 4]
        assert S.stabilizer(3).stabilizer(4) == S.subgroup_search(prop_fix_points)
        points = [3, 5]
        fix35 = A.subgroup_search(prop_fix_points)
        points = [5]
        fix5 = A.subgroup_search(prop_fix_points)
        assert A.subgroup_search(prop_fix_points, init_subgroup=fix35) == fix5
        base, strong_gens = A.schreier_sims_incremental()
        g = A.generators[0]
        comm_g = A.subgroup_search(prop_comm_g, base=base, strong_gens=strong_gens)
        assert _verify_bsgs(comm_g, base, comm_g.generators) == True
        assert [prop_comm_g(gen) == True for gen in comm_g.generators]

def test_normal_closure():
    # the normal closure of the trivial group is trivial
    S = SymmetricGroup(3)
    identity = Permutation([0, 1, 2])
    closure = S.normal_closure(identity)
    assert closure.is_trivial
    # the normal closure of the entire group is the entire group
    A = AlternatingGroup(4)
    assert A.normal_closure(A) == A
    # brute-force verifications for subgroups
    for i in (3, 4, 5):
        S = SymmetricGroup(i)
        A = AlternatingGroup(i)
        D = DihedralGroup(i)
        C = CyclicGroup(i)
        for gp in (A, D, C):
            assert _verify_normal_closure(S, gp)
    # brute-force verifications for all elements of a group
    S = SymmetricGroup(5)
    elements = list(S.generate_dimino())
    for element in elements:
        assert _verify_normal_closure(S, element)
    # small groups
    small = []
    for i in (1, 2, 3):
        small.append(SymmetricGroup(i))
        small.append(AlternatingGroup(i))
        small.append(DihedralGroup(i))
        small.append(CyclicGroup(i))
    for gp in small:
        for gp2 in small:
            if gp2.is_subgroup(gp):
                assert _verify_normal_closure(gp, gp2)

def test_derived_series():
    # the derived series of the trivial group consists only of the trivial group
    triv = PermutationGroup([Permutation([0, 1, 2])])
    assert triv.derived_series() == [triv]
    # the derived series for a simple group consists only of the group itself
    for i in (5, 6, 7):
        A = AlternatingGroup(i)
        assert A.derived_series() == [A]
    # the derived series for S_4 is S_4 > A_4 > K_4 > triv
    S = SymmetricGroup(4)
    series = S.derived_series()
    assert series[1] == AlternatingGroup(4)
    assert series[2] == DihedralGroup(2)
    assert series[3].is_trivial

def test_lower_central_series():
    # the lower central series of the trivial group consists only of the trivial group
    triv = PermutationGroup([Permutation([0, 1, 2])])
    assert triv.lower_central_series() == [triv]
    # the lower central series of a simple group consists only of the group itself
    for i in (5, 6, 7):
        A = AlternatingGroup(i)
        assert A.lower_central_series() == [A]
    # GAP-verified example
    S = SymmetricGroup(6)
    series = S.lower_central_series()
    assert len(series) == 2
    assert series[1] == AlternatingGroup(6)

def test_commutator():
    # the commutator of the trivial group and the trivial group is trivial
    S = SymmetricGroup(3)
    triv = PermutationGroup([Permutation([0, 1, 2])])
    assert S.commutator(triv, triv) == triv
    # the commutator of the trivial group and any other group is again trivial
    A = AlternatingGroup(3)
    assert S.commutator(triv, A) == triv
    # the commutator is commutative
    for i in (3, 4, 5):
        S = SymmetricGroup(i)
        A = AlternatingGroup(i)
        D = DihedralGroup(i)
        assert S.commutator(A, D) == S.commutator(D, A)
    # the commutator of an abelian group is trivial
    S = SymmetricGroup(7)
    A1 = AbelianGroup(2, 5)
    A2 = AbelianGroup(3, 4)
    triv = PermutationGroup([Permutation([0, 1, 2, 3, 4, 5, 6])])
    assert S.commutator(A1, A1) == triv
    assert S.commutator(A2, A2) == triv
    # examples calculated by hand
    S = SymmetricGroup(3)
    A = AlternatingGroup(3)
    assert S.commutator(A, S) == A

def test_is_nilpotent():
    # every abelian group is nilpotent
    for i in (1, 2, 3):
        C = CyclicGroup(i)
        Ab = AbelianGroup(i, i + 2)
        assert C.is_nilpotent
        assert Ab.is_nilpotent
    Ab = AbelianGroup(5, 7, 10)
    assert Ab.is_nilpotent
    # A_5 is not solvable and thus not nilpotent
    assert AlternatingGroup(5).is_nilpotent == False

def test_is_trivial():
    for i in range(5):
        triv = PermutationGroup([Permutation(range(i))])
        assert triv.is_trivial

def test_pointwise_stabilizer():
    S = SymmetricGroup(5)
    points = []
    stab = S
    for point in (2, 0, 3, 4, 1):
        stab = stab.stabilizer(point)
        points.append(point)
        assert S.pointwise_stabilizer(points) == stab
