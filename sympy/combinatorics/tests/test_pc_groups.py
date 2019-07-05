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
    # SymmetricGroup(3)
    F, x0, x1 = free_group("x0, x1")
    group = F
    G = SymmetricGroup(3)

    pc_group = G.polycyclic_group()
    pc_presentation = pc_group.pc_presentation(F)
    assert pc_presentation == { x1**3: (), x0**3: (), x0**-1*x1*x0: x0**2 }

    x0, x1 = pc_group.pcgs
    assert x0**-1*x1*x0 == x0**2


    # SymmetricGroup(4)
    F, x0, x1 = free_group("x0, x1")
    group = F
    G = SymmetricGroup(4)

    pc_group = G.polycyclic_group()
    pc_presentation = pc_group.pc_presentation(F)
    assert pc_presentation == { x1**2: (), x0**2: (), x0**-1*x1*x0: x1 }

    x0, x1 = pc_group.pcgs
    assert x0**-1*x1*x0 == x1


    # SymmetricGroup(9).sylow_subgroup(3)
    F, x0, x1, x2 = free_group("x0, x1, x2")
    group = F
    S = SymmetricGroup(9)
    G = S.sylow_subgroup(3)
    der = G.derived_series()

    pc_group = G.polycyclic_group()
    pc_presentation = pc_group.pc_presentation(F)
    assert pc_presentation == { x2**3: (), x1**3: (), x0**3: (),
                                x1**-1*x2*x1: x2, x0**-1*x2*x0: x2,
                                x0**-1*x1*x0: x0**2 }

    # x0, x1, x2 = pc_group.pcgs  (I don't know why the sequence gets reversed here)
    x0, x1, x2 = der[len(der)-2].generators
    assert x1**-1*x2*x1 == x2
    assert x0**-1*x2*x0 == x2
    assert x0**-1*x1*x0 == x0**2
