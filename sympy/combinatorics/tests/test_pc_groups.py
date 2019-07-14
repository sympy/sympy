from sympy.combinatorics.pc_groups import PolycyclicGroup, Collector
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.free_groups import free_group
from sympy.combinatorics.named_groups import SymmetricGroup

def test_collected_word():
    F, x0, x1, x2, x3 = free_group("x0, x1, x2, x3")

    # Polycyclic relators for SymmetricGroup(4)
    pc_relators = { x0**2: (), x1**3: (), x2**2: (), x3**2: (),
                    x0**-1*x1*x0: x1**2, x0**-1*x2*x0: x2*x3,
                    x0**-1*x3*x0: x3, x1**-1*x2*x1: x3,
                    x1**-1*x3*x1: x2*x3, x2**-1*x3*x2: x3
                  }

    word = x3*x2*x1*x0
    relative_order = [2, 3, 2, 2]
    group = word.group
    collector = Collector(pc_relators, relative_order, group)
    collected_word_ = collector.collected_word(word)

    assert collected_word_ == x0*x1**2*x2*x3

    # Polycyclic Generators of SymmetricGroup(4)
    x0 = Permutation(0, 1)
    x1 = Permutation(0, 1, 2)
    x2 = Permutation(0, 2)(1, 3)
    x3 = Permutation(0, 1)(2, 3)

    word = x3*x2*x1*x0
    collected_word_ = x0*x1**2*x2*x3
    assert word == collected_word_



    F, x0, x1 = free_group("x0, x1")
    # polycyclic relators for Symmetricgroup(3)
    pc_relators = {x0**2: (), x1**3: (), x0**-1*x1*x0: x1**2}
    relative_order = [2, 3]
    group = F
    collector = Collector(pc_relators, relative_order, group)

    a = Permutation(0, 1) # x0
    b = Permutation(0, 1, 2) # x1

    word = x1*x0
    assert collector.collected_word(word) == x0*x1**2
    assert b*a == a*b**2

    word = x1*x0**2
    assert collector.collected_word(word) == x1
    assert b*a**2 == b

    word = x1**2*x0
    assert collector.collected_word(word) == x0*x1
    assert b**2*a == a*b

    word = x1**4*x0**6
    assert collector.collected_word(word) == x1
    assert b**4*a**6 == b

    word = x0*x1
    # The word is already collected
    assert collector.collected_word(word) == x0*x1
    assert a*b == a*b

    word = x0**2*x1
    assert collector.collected_word(word) == x1
    assert a**2*b == b

    word = x0**2*x1**3
    # Handle Identity case
    assert collector.collected_word(word) == F.identity
    assert a**2*b**3 == Permutation(2)

    word = x1**-2*x0
    assert collector.collected_word(word) == x0*x1**2
    assert b**-2*a == a*b**2


def test_pc_presentation():
    F1, x0, x1 = free_group("x0, x1")
    F2, x0, x1, x2, x3 = free_group("x0, x1, x2, x3")
    F3, x0, x1, x2, x3, x4, x5, x6 = free_group("x0, x1, x2, x3, x4, x5, x6")

    l = [(SymmetricGroup(3), F1), (SymmetricGroup(4), F2),
         (SymmetricGroup(9).sylow_subgroup(3), F2), (SymmetricGroup(9).sylow_subgroup(2), F3),
         (SymmetricGroup(8).sylow_subgroup(2), F3)]

    for t in l:
        pc_group = t[0].polycyclic_group()
        pc_presentation = pc_group.pc_presentation(t[1])

        pcgs = pc_group.pcgs
        free_to_perm = {}
        for s, g in zip(t[1].symbols, pcgs):
            free_to_perm[s] = g

        for k, v in pc_presentation.items():
            k_array = k.array_form
            if v != ():
                v_array = v.array_form

            lhs = Permutation()
            for gen in k_array:
                s = gen[0]
                e = gen[1]
                lhs = lhs*free_to_perm[s]**e

            if v == ():
                assert lhs.is_identity
                continue

            rhs = Permutation()
            for gen in v_array:
                s = gen[0]
                e = gen[1]
                rhs = rhs*free_to_perm[s]**e

            assert lhs == rhs


def test_exponent_vector():
    F1, x0, x1 = free_group("x0, x1")
    F2, x0, x1, x2, x3 = free_group("x0, x1, x2, x3")
    F3, x0, x1, x2, x3, x4, x5, x6 = free_group("x0, x1, x2, x3, x4, x5, x6")

    l = [(SymmetricGroup(3), F1), (SymmetricGroup(4), F2),
         (SymmetricGroup(9).sylow_subgroup(3), F2), (SymmetricGroup(9).sylow_subgroup(2), F3),
         (SymmetricGroup(8).sylow_subgroup(2), F3)]

    for t in l:
        PcGroup = t[0].polycyclic_group()
        pcgs = PcGroup.pcgs

        for gen in t[0].generators:
            exp = PcGroup.exponent_vector(gen, t[1])
            g = Permutation()
            for i in range(len(exp)):
                g = g*pcgs[i] if exp[i] else g

            assert g == gen
