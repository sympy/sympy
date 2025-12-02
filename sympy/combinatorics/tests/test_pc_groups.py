from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.named_groups import (SymmetricGroup,
    AlternatingGroup, DihedralGroup)
from sympy.matrices import Matrix
from sympy.combinatorics import free_group

def test_pc_presentation():
    Groups = [SymmetricGroup(3), SymmetricGroup(4),
         SymmetricGroup(9).sylow_subgroup(3),
         SymmetricGroup(9).sylow_subgroup(2),
         SymmetricGroup(8).sylow_subgroup(2), DihedralGroup(10)]

    S = SymmetricGroup(125).sylow_subgroup(5)
    G = S.derived_series()[2]
    Groups.append(G)

    G = SymmetricGroup(25).sylow_subgroup(5)
    Groups.append(G)

    S = SymmetricGroup(11**2).sylow_subgroup(11)
    G = S.derived_series()[2]
    Groups.append(G)

    for G in Groups:
        PcGroup = G.polycyclic_group()
        collector = PcGroup.collector
        pc_presentation = collector.pc_presentation

        pcgs = PcGroup.pcgs
        free_group = collector.free_group
        free_to_perm = {}
        for s, g in zip(free_group.symbols, pcgs):
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

    Groups = [SymmetricGroup(3), SymmetricGroup(4),
         SymmetricGroup(9).sylow_subgroup(3),
         SymmetricGroup(9).sylow_subgroup(2),
         SymmetricGroup(8).sylow_subgroup(2)]

    for G in Groups:
        PcGroup = G.polycyclic_group()
        collector = PcGroup.collector

        pcgs = PcGroup.pcgs

        for gen in G.generators:
            exp = collector.exponent_vector(gen)
            g = Permutation()
            for i in range(len(exp)):
                g = g*pcgs[i]**exp[i] if exp[i] else g
            assert g == gen


def test_induced_pcgs():
    G = [SymmetricGroup(9).sylow_subgroup(3),
    SymmetricGroup(20).sylow_subgroup(2), AlternatingGroup(4),
    DihedralGroup(4), DihedralGroup(10), DihedralGroup(9),
    SymmetricGroup(3), SymmetricGroup(4)]

    for g in G:
        PcGroup = g.polycyclic_group()
        collector = PcGroup.collector
        gens = list(g.generators)
        ipcgs = collector.induced_pcgs(gens)
        m = []
        for i in ipcgs:
            m.append(collector.exponent_vector(i))
        assert Matrix(m).is_upper


def test_collected_word_issue_28637():
    """
    Test that collected_word preserves permutation information
    for negative exponents.
    
    This addresses issue #28637 where collector.collected_word
    was dropping permutation data for certain words with negative
    exponents.
    
    References
    ==========
    
    .. [1] https://github.com/sympy/sympy/issues/28637
    """
    
    def perm_from_word(pcgs, word):
        """Convert a free group word to a permutation."""
        if pcgs:
            perm = Permutation(pcgs[0].size - 1)
        else:
            perm = Permutation()
        
        for sym, exp in word.array_form:
            idx = int(str(sym)[1:])
            perm *= pcgs[idx] ** exp
        return perm
    
    G = SymmetricGroup(4)
    pc_group = G.polycyclic_group()
    collector = pc_group.collector
    pcgs = collector.pcgs
    F, *gens = free_group(','.join([f'x{i}' for i in range(len(pcgs))]))
    
    word = gens[2] ** -2
    collected = collector.collected_word(word)
    word_perm = perm_from_word(pcgs, word)
    collected_perm = perm_from_word(pcgs, collected)
    assert word_perm == collected_perm, (
        f"collected_word dropped permutation info: "
        f"{word_perm} != {collected_perm}")
    
    for i in range(len(pcgs)):
        for exp in [-4, -3, -2, -1]:
            word = gens[i] ** exp
            collected = collector.collected_word(word)
            word_perm = perm_from_word(pcgs, word)
            collected_perm = perm_from_word(pcgs, collected)
            assert word_perm == collected_perm, (
                f"x{i}^{exp}: collected_word changed permutation: "
                f"{word_perm} != {collected_perm}")
    
    test_words = [
        gens[0] ** -1 * gens[2] ** -2,
        gens[1] ** -3 * gens[3] ** -1,
        gens[2] * gens[2] ** -3,
    ]
    
    for word in test_words:
        collected = collector.collected_word(word)
        word_perm = perm_from_word(pcgs, word)
        collected_perm = perm_from_word(pcgs, collected)
        assert word_perm == collected_perm, (
            f"collected_word changed permutation for {word}: "
            f"{word_perm} != {collected_perm}")
    
    for G in [SymmetricGroup(3), SymmetricGroup(4), DihedralGroup(8)]:
        pc_group = G.polycyclic_group()
        collector = pc_group.collector
        pcgs = collector.pcgs
        F, *gens = free_group(','.join([f'x{i}' for i in range(len(pcgs))]))
        
        for i in range(len(pcgs)):
            word = gens[i] ** -2
            collected = collector.collected_word(word)
            word_perm = perm_from_word(pcgs, word)
            collected_perm = perm_from_word(pcgs, collected)
            assert word_perm == collected_perm
