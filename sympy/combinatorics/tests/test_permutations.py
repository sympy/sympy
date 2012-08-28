from sympy.core import FiniteSet
from sympy.combinatorics.permutations import (Permutation, _af_parity,
    _af_mul, _af_mul, DisjointCycle as DC)
from sympy.core.compatibility import permutations

from sympy.utilities.pytest import raises

def test_Permutation():
    assert Permutation([1,2,3]) == Permutation([0, 1, 2, 3])
    assert Permutation([[1,2]], size=4) == Permutation([[1, 2], [0], [3]])
    p = Permutation([2, 5, 1, 6, 3, 0, 4])
    q = Permutation([[1], [0, 3, 5, 6, 2, 4]])
    r = Permutation([1,3,2,0,4,6,5])
    ans = Permutation(_af_mul(*[w.array_form for w in (p, q, r)])).array_form
    assert (p*q*r).array_form == ans
    # make sure no other permutation of p, q, r could have given
    # that answer
    for a, b, c in permutations((p, q, r)):
        if (a, b, c) == (p, q, r): continue
        assert (a*b*c).array_form != ans

    assert p.support() == range(7)
    assert q.support() == [0, 2, 3, 4, 5, 6]
    assert Permutation(p.cyclic_form).array_form == p.array_form
    assert p.cardinality == 5040
    assert q.cardinality == 5040
    assert q.cycles == 2
    assert q*p == Permutation([4, 6, 1, 2, 5, 3, 0])
    assert p*q == Permutation([6, 5, 3, 0, 2, 4, 1])
    assert _af_mul(p.array_form, q.array_form) == \
        [6, 5, 3, 0, 2, 4, 1]

    assert (Permutation([[1,2,3],[0,4]])*Permutation([[1,2,4],[0],[3]])).cyclic_form == \
        [[0, 4, 2], [1, 3]]
    assert q.array_form == [3, 1, 4, 5, 0, 6, 2]
    assert q.cyclic_form == [[0, 3, 5, 6, 2, 4]]
    assert q.full_cyclic_form == [[0, 3, 5, 6, 2, 4], [1]]
    assert p.cyclic_form == [[0, 2, 1, 5], [3, 6, 4]]
    assert p.transpositions() == FiniteSet([(0, 1), (0, 2), (0, 5), (3, 4), (3, 6)])

    assert p**13 == p
    assert q**2 == Permutation([5, 1, 0, 6, 3, 2, 4])

    assert p+q == Permutation([5, 6, 3, 1, 2, 4, 0])
    assert q+p == p+q

    assert p-q == Permutation([6, 3, 5, 1, 2, 4, 0])
    assert q-p == Permutation([1, 4, 2, 6, 5, 3, 0])

    a = p-q
    b = q-p
    assert (a+b).is_Identity

    pq = p.conjugate(q)
    assert pq == Permutation([5, 3, 0, 4, 6, 2, 1])
    assert pq == ~q*p*q
    assert pq == p**q
    qp = q.conjugate(p)
    assert qp == Permutation([6, 3, 2, 0, 1, 4, 5])
    assert qp == ~p*q*p
    assert qp == q**p

    assert p.commutator(q) == Permutation([1, 4, 5, 6, 3, 0, 2])
    assert q.commutator(p) == Permutation([5, 0, 6, 4, 1, 2, 3])
    assert p.commutator(q) == ~ q.commutator(p)

    assert len(p.atoms()) == 7
    assert q.atoms() == set([0, 1, 2, 3, 4, 5, 6])

    assert p.inversion_vector() == [2, 4, 1, 3, 1, 0]
    assert q.inversion_vector() == [3, 1, 2, 2, 0, 1]

    assert Permutation.from_inversion_vector(p.inversion_vector()) == p
    assert Permutation.from_inversion_vector(q.inversion_vector()).array_form\
           == q.array_form
    assert Permutation([i for i in range(500,-1,-1)]).inversions() == 125250

    assert Permutation([0, 4, 1, 3, 2]).parity() == 0
    assert Permutation([0, 1, 4, 3, 2]).parity() == 1
    assert _af_parity([0, 4, 1, 3, 2]) == 0
    assert _af_parity([0, 1, 4, 3, 2]) == 1

    s = Permutation([0])

    assert s.is_Singleton

    r = Permutation([3, 2, 1, 0])
    assert (r**2).is_Identity

    assert (p*(~p)).is_Identity
    assert (~p)**13 == Permutation([5, 2, 0, 4, 6, 1, 3])
    assert ~(r**2).is_Identity
    assert p.max() == 6
    assert p.min() == 0

    q = Permutation([[6], [5], [0, 1, 2, 3, 4]])

    assert q.max() == 4
    assert q.min() == 0

    p = Permutation([1, 5, 2, 0, 3, 6, 4])
    q = Permutation([[1, 2, 3, 5, 6], [0, 4]])

    assert p.ascents() == [0, 3, 4]
    assert q.ascents() == [1, 2, 4]
    assert r.ascents() == []

    assert p.descents() == [1, 2, 5]
    assert q.descents() == [0, 3, 5]
    assert Permutation(r.descents()).is_Identity

    assert p.inversions() == 7
    assert p.signature() == -1
    assert q.inversions() == 11
    assert q.signature() == -1
    assert (p*(~p)).inversions() == 0
    assert (p*(~p)).signature() == 1

    assert p.order() == 6
    assert q.order() == 10
    assert (p**(p.order())).is_Identity

    assert p.length() == 6
    assert q.length() == 7
    assert r.length() == 4

    assert p.runs() == [[1, 5], [2], [0, 3, 6], [4]]
    assert q.runs() == [[4], [2, 3, 5], [0, 6], [1]]
    assert r.runs() == [[3], [2], [1], [0]]

    assert p.index() == 8
    assert q.index() == 8
    assert r.index() == 3

    assert p.get_precedence_distance(q) == q.get_precedence_distance(p)
    assert p.get_adjacency_distance(q) == p.get_adjacency_distance(q)
    assert p.get_positional_distance(q) == p.get_positional_distance(q)
    p = Permutation([0, 1, 2, 3])
    q = Permutation([3, 2, 1, 0])
    assert p.get_precedence_distance(q) == 6
    assert p.get_adjacency_distance(q) == 3
    assert p.get_positional_distance(q) == 8

    a = [Permutation.unrank_nonlex(4, i) for i in range(5)]
    iden = Permutation([0, 1, 2, 3])
    for i in range(5):
        for j in range(i+1, 5):
            assert a[i].commutes_with(a[j]) == (a[i]*a[j] == a[j]*a[i])
            if a[i].commutes_with(a[j]):
                assert a[i].commutator(a[j]) == iden
                assert a[j].commutator(a[i]) == iden

def test_josephus():
    assert Permutation.josephus(4, 6, 1) == Permutation([3, 1, 0, 2, 5, 4])
    assert Permutation.josephus(1, 5, 1).is_Identity

def test_ranking():
    assert Permutation.unrank_lex(5, 10).rank() == 10
    p = Permutation.unrank_lex(15, 225)
    assert p.rank() == 225
    p1 = p.next_lex()
    assert p1.rank() == 226
    assert Permutation.unrank_lex(15, 225).rank() == 225
    assert Permutation.unrank_lex(10, 0).is_Identity
    p = Permutation.unrank_lex(4, 23)
    assert p.rank() == 23
    assert p.array_form == [3, 2, 1, 0]
    assert p.next_lex() == None

    p = Permutation([1, 5, 2, 0, 3, 6, 4])
    q = Permutation([[1, 2, 3, 5, 6], [0, 4]])
    a = [Permutation.unrank_trotterjohnson(4, i).array_form for i in range(5)]
    assert a == [[0,1,2,3], [0,1,3,2], [0,3,1,2], [3,0,1,2], [3,0,2,1] ]
    assert [Permutation(pa).rank_trotterjohnson() for pa in a] == range(5)
    assert Permutation([0,1,2,3]).next_trotterjohnson() == \
        Permutation([0,1,3,2])

    assert q.rank_trotterjohnson() == 2283
    assert p.rank_trotterjohnson() == 3389

    p = Permutation([2, 5, 1, 6, 3, 0, 4])
    q = Permutation([[6], [5], [0, 1, 2, 3, 4]])
    assert p.rank() == 1964
    assert q.rank() == 870
    assert Permutation([]).rank_nonlex() == 0
    prank = p.rank_nonlex()
    assert prank == 1600
    assert Permutation.unrank_nonlex(7, 1600) == p
    qrank = q.rank_nonlex()
    assert qrank == 41
    assert Permutation.unrank_nonlex(7, 41) == Permutation(q.array_form)

    a = [Permutation.unrank_nonlex(4, i).array_form for i in range(24)]
    assert a == \
    [[1, 2, 3, 0], [3, 2, 0, 1], [1, 3, 0, 2], [1, 2, 0, 3], [2, 3, 1, 0], \
     [2, 0, 3, 1], [3, 0, 1, 2], [2, 0, 1, 3], [1, 3, 2, 0], [3, 0, 2, 1], \
     [1, 0, 3, 2], [1, 0, 2, 3], [2, 1, 3, 0], [2, 3, 0, 1], [3, 1, 0, 2], \
     [2, 1, 0, 3], [3, 2, 1, 0], [0, 2, 3, 1], [0, 3, 1, 2], [0, 2, 1, 3], \
     [3, 1, 2, 0], [0, 3, 2, 1], [0, 1, 3, 2], [0, 1, 2, 3]]

    assert Permutation([3, 2, 0, 1]).next_nonlex() == Permutation([1, 3, 0, 2])
    assert [Permutation(pa).rank_nonlex() for pa in a] == range(24)

def test_mul():
    a, b = [0,2,1,3], [0,1,3,2]
    assert _af_mul(a, b) == [0, 2, 3, 1]
    assert _af_mul(a, b, range(4)) == [0, 2, 3, 1]
    assert (Permutation(a)*Permutation(b)).array_form == [0, 2, 3, 1]
    assert _af_mul([0,2,1,3], [0,1,3,2], reverse=True) == [0, 3, 1, 2]

    a = Permutation([0, 2, 1, 3])
    b = (0, 1, 3, 2)
    c = (3, 1, 2, 0)
    assert Permutation.mul(a, b, c) == Permutation([1, 2, 3, 0])
    assert Permutation.mul(a, b, c, reverse=True) == Permutation([3, 0, 1, 2])
    assert Permutation.mul(b, c) == [2, 1, 3, 0]

    n = 6
    m = 8
    a = [Permutation.unrank_nonlex(n, i).array_form for i in range(m)]
    h = range(n)
    for i in range(m):
        h = _af_mul(h, a[i])
        h2 = _af_mul(*a[:i+1])
        assert h == h2

def test_args():
    p = Permutation([(0, 3, 1, 2), (4, 5)])
    assert p.cyclic_form == [[0, 3, 1, 2], [4, 5]]
    assert p._array_form == None
    p = Permutation((0, 3, 1, 2))
    assert p._cyclic_form == None
    assert p._array_form == [0, 3, 1, 2]
    assert Permutation([0]) == Permutation((0, ))
    assert Permutation([[0], [1]]) == Permutation(((0, ), (1, ))) == Permutation(((0, ), [1]))
    raises(TypeError, lambda: Permutation(0, 1, 2)) # enclosing brackets needed
    raises(TypeError, lambda: Permutation([1, 2], [0])) # enclosing brackets needed
    raises(ValueError, lambda: Permutation([[1, 2], 0])) # enclosing brackets needed on 0

def test_DisjointCycle():
    a, b = DC(1, 2), DC(2, 3)
    ans = DC(1, 3, 2)
    assert a*b == ans
    assert a*(2, 3) == ans
    assert a*[(2, 3), (4, 5)] == DC(1, 3, 2)*DC(4, 5)
    assert a*Permutation([2, 1, 0, 3]) == DC(0, 2, 1)
    raises(ValueError, lambda: DC().as_list())
    assert DC(1,2).as_list() == [0, 2, 1]
    assert DC(1,2).as_list(4) == [0, 2, 1, 3]
    assert DC(1,2).as_permutation(4).array_form == \
        Permutation([0, 2, 1, 3]).array_form
    assert str(DC(1,2)*DC(4,5)) == '[(1, 2), (4, 5)]'
    assert str(DC(1,2)) == '[(1, 2)]'
    assert DC(Permutation(range(3))) == DC()
    assert DC(1,2).as_list() == [0, 2, 1]
    assert DC(1,2).as_list(4) == [0, 2, 1, 3]
    raises(TypeError, lambda: DC((1,2)))
    raises(ValueError, lambda: DC(1,2,1))
