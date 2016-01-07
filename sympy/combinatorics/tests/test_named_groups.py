from sympy.combinatorics.named_groups import (SymmetricGroup, CyclicGroup,
DihedralGroup, AlternatingGroup, AbelianGroup)


def test_SymmetricGroup():
    s4 = SymmetricGroup(4)
    s1 = SymmetricGroup(1)
    s2 = SymmetricGroup(2)
    s5 = SymmetricGroup(5)
    elements = list(s5.generate())
    assert (s5.generators[0]).size == 5
    assert len(elements) == 120
    assert s5.is_solvable is False
    assert s5.is_abelian is False
    assert s5.is_nilpotent is False
    assert s5.is_transitive() is True
    assert s1.order() == 1
    assert s2.order() == 2

    assert s4.stabilizer(2) == PermGroup(Perm(0, 1, 3), Perm(3)(0, 1))
    assert s4.stabilizer(3) == PermGroup(Perm(3)(0, 1, 2), Perm(3)(0, 1))

    assert s4.stabilizer(1) == PermGroup(Perm(0, 2, 3), Perm(3)(0, 2))
    assert s4.stabilizer(0) == PermGroup(Perm(1, 2, 3), Perm(3)(1, 2))
    assert s1.stabilizer(0) == SymmetricGroup(1)
    assert s2.stabilizer(0) == PermGroup(Perm(1))


def test_CyclicGroup():
    G = CyclicGroup(10)
    elements = list(G.generate())
    assert len(elements) == 10
    assert (G.derived_subgroup()).order() == 1
    assert G.is_abelian is True
    assert G.is_solvable is True
    assert G.is_nilpotent is True
    H = CyclicGroup(1)
    assert H.order() == 1
    L = CyclicGroup(2)
    assert L.order() == 2


def test_DihedralGroup():
    G = DihedralGroup(6)
    elements = list(G.generate())
    assert len(elements) == 12
    assert G.is_transitive() is True
    assert G.is_abelian is False
    assert G.is_solvable is True
    assert G.is_nilpotent is False
    H = DihedralGroup(1)
    assert H.order() == 2
    L = DihedralGroup(2)
    assert L.order() == 4
    assert L.is_abelian is True
    assert L.is_nilpotent is True


def test_AlternatingGroup():
    G = AlternatingGroup(5)
    elements = list(G.generate())
    assert len(elements) == 60
    assert [perm.is_even for perm in elements] == [True]*60
    H = AlternatingGroup(1)
    assert H.order() == 1
    L = AlternatingGroup(2)
    assert L.order() == 1


def test_AbelianGroup():
    A = AbelianGroup(3, 3, 3)
    assert A.order() == 27
    assert A.is_abelian is True
