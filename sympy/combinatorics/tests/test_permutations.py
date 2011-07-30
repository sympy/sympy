from sympy.combinatorics.permutations import Permutation

def test_Permutation():
    p = Permutation([2,5,1,6,3,0,4])
    q = Permutation([[1,4,5],[2,0,6],[3]])

    assert q.cycles == 3
    assert p*q == Permutation([4, 6, 1, 2, 5, 3, 0])
    assert q*p == Permutation([6, 5, 3, 0, 2, 4, 1])

    assert q.array_form == [3, 1, 4, 5, 0, 6, 2]
    assert p.cyclic_form == [[3, 6, 4], [0, 2, 1, 5]]

    assert p**13 == p
    assert q**2 == Permutation([5, 1, 0, 6, 3, 2, 4])

    assert p+q == Permutation([5, 6, 3, 1, 2, 4, 0])
    assert q+p == p+q

    assert p-q == Permutation([6, 3, 5, 1, 2, 4, 0])
    assert q-p == Permutation([1, 4, 2, 6, 5, 3, 0])

    a = p-q
    b = q-p
    assert (a+b).is_Identity

    assert len(p.atoms()) == 7
    assert q.atoms() == set([0, 1, 2, 3, 4, 5, 6])

    assert p.inversion_vector == [2, 4, 1, 3, 1, 0]
    assert q.inversion_vector == [3, 1, 2, 2, 0, 1]

    assert Permutation.from_inversion_vector(p.inversion_vector) == p
    assert Permutation.from_inversion_vector(q.inversion_vector).array_form\
           == q.array_form

    s = Permutation([0])

    assert s.is_Singleton

    r = Permutation([3,2,1,0])
    assert (r**2).is_Identity

    assert (p*(~p)).is_Identity
    assert (~p)**13 == Permutation([5, 2, 0, 4, 6, 1, 3])
    assert ~(r**2).is_Identity
    assert p.max == 6
    assert p.min == 0

    q = Permutation([[4,1,2,3],[0,5,6]])

    assert q.max == 4
    assert q.min == 0

    assert p.rank_nonlex() == 14830
    assert q.rank_nonlex() == 8441
    assert Permutation.unrank_nonlex(7, 41) == Permutation([4, 2, 3, 5, 1, 0, 6])

    assert q.rank == 870
    assert p.rank == 1964

    p = Permutation([1,5,2,0,3,6,4])
    q = Permutation([[2,3,5],[1,0,6],[4]])

    assert p.ascents == [0, 3, 4]
    assert q.ascents == [1, 2, 4]
    assert r.ascents == []

    assert p.descents == [1, 2, 5]
    assert q.descents == [0, 3, 5]
    assert Permutation(r.descents).is_Identity

    assert p.inversions == 7
    assert p.signature == -1
    assert q.inversions == 11
    assert q.signature == -1
    assert (p*(~p)).inversions == 0
    assert (p*(~p)).signature == 1

    assert p.order == 6
    assert q.order == 3
    assert (p**(p.order)).is_Identity

    assert p.length == 6
    assert q.length == 7
    assert r.length == 4

    assert not p.is_Positive
    assert p.is_Negative
    assert not q.is_Positive
    assert q.is_Negative
    assert r.is_Positive
    assert not r.is_Negative

    assert p.runs() == [[1, 5], [2], [0, 3, 6], [4]]
    assert q.runs() == [[4], [2, 3, 5], [0, 6], [1]]
    assert r.runs() == [[3], [2], [1], [0]]

    assert p.index == 8
    assert q.index == 8
    assert r.index == 3

    assert q.rank_trotterjohnson == 259
    assert p.rank_trotterjohnson == 1087

    assert p.get_precedence_distance(q) == q.get_precedence_distance(p)
    assert p.get_adjacency_distance(q) == p.get_adjacency_distance(q)
    assert p.get_positional_distance(q) == p.get_positional_distance(q)
    p = Permutation([0, 1, 2, 3])
    q = Permutation([3, 2, 1, 0])
    assert p.get_precedence_distance(q) == 6
    assert p.get_adjacency_distance(q) == 3
    assert p.get_positional_distance(q) == 8

def test_josephus():
    assert Permutation.josephus(4, 6, 1) == Permutation([3, 1, 0, 2, 5, 4])
    assert Permutation.josephus(1, 5, 1).is_Identity

def test_unrank_lex():
    assert Permutation.unrank_lex(5, 10).rank == 10
    assert Permutation.unrank_lex(15, 225).rank == 225
    assert Permutation.unrank_lex(10, 0).is_Identity
    assert Permutation.unrank_lex(4, 23).array_form == [0, 3, 2, 1]
