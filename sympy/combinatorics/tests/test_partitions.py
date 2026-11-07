from __future__ import annotations
from sympy.core.sorting import ordered, default_sort_key
from sympy.combinatorics.partitions import (Partition, IntegerPartition,
                                            RGS_enum, RGS_unrank, RGS_rank,
                                            random_integer_partition)
from sympy.testing.pytest import raises
from sympy.utilities.iterables import partitions
from sympy.sets.sets import Set, FiniteSet


def test_partition_constructor():
    raises(ValueError, lambda: Partition([1, 1, 2]))
    raises(ValueError, lambda: Partition([1, 2, 3], [2, 3, 4]))
    raises(ValueError, lambda: Partition(1, 2, 3))
    raises(ValueError, lambda: Partition(*list(range(3))))

    assert Partition([1, 2, 3], [4, 5]) == Partition([4, 5], [1, 2, 3])
    assert Partition({1, 2, 3}, {4, 5}) == Partition([1, 2, 3], [4, 5])

    a = FiniteSet(1, 2, 3)
    b = FiniteSet(4, 5)
    assert Partition(a, b) == Partition([1, 2, 3], [4, 5])
    assert Partition({a, b}) == Partition(FiniteSet(a, b))
    assert Partition({a, b}) != Partition(a, b)

def test_partition():
    from sympy.abc import x

    a = Partition([1, 2, 3], [4])
    b = Partition([1, 2], [3, 4])
    c = Partition([x])
    l = [a, b, c]
    l.sort(key=default_sort_key)
    assert l == [c, a, b]
    l.sort(key=lambda w: default_sort_key(w, order='rev-lex'))
    assert l == [c, a, b]

    assert (a == b) is False
    assert a <= b
    assert (a > b) is False
    assert a != b
    assert a < b

    assert (a + 2).partition == [[1, 2], [3, 4]]
    assert (b - 1).partition == [[1, 2, 4], [3]]

    assert (a - 1).partition == [[1, 2, 3, 4]]
    assert (a + 1).partition == [[1, 2, 4], [3]]
    assert (b + 1).partition == [[1, 2], [3], [4]]

    assert a.rank == 1
    assert b.rank == 3

    assert a.RGS == (0, 0, 0, 1)
    assert b.RGS == (0, 0, 1, 1)


def test_partition_aliasing():
    p = Partition([1, 2], [3, 4])
    blocks = p.partition
    blocks[0].append(99)
    assert p.partition == [[1, 2], [3, 4]]


def test_integer_partition():
    # no zeros in partition
    raises(ValueError, lambda: IntegerPartition(list(range(3))))
    # check fails since 1 + 2 != 100
    raises(ValueError, lambda: IntegerPartition(100, list(range(1, 3))))
    a = IntegerPartition(8, [1, 3, 4])
    b = a.next_lex()
    c = IntegerPartition([1, 3, 4])
    d = IntegerPartition(8, {1: 3, 3: 1, 2: 1})
    assert a == c
    assert a.integer == d.integer
    assert a.conjugate == [3, 2, 2, 1]
    assert (a == b) is False
    assert a <= b
    assert (a > b) is False
    assert a != b

    for i in range(1, 11):
        next = set()
        prev = set()
        a = IntegerPartition([i])
        ans = {IntegerPartition(p) for p in partitions(i)}
        n = len(ans)
        for j in range(n):
            next.add(a)
            a = a.next_lex()
            IntegerPartition(i, a.partition)  # check it by giving i
        for j in range(n):
            prev.add(a)
            a = a.prev_lex()
            IntegerPartition(i, a.partition)  # check it by giving i
        assert next == ans
        assert prev == ans

    assert IntegerPartition([1, 2, 3]).as_ferrers() == '###\n##\n#'
    assert IntegerPartition([1, 1, 3]).as_ferrers('o') == 'ooo\no\no'
    assert str(IntegerPartition([1, 1, 3])) == '[3, 1, 1]'
    assert IntegerPartition([1, 1, 3]).partition == [3, 1, 1]

    raises(ValueError, lambda: random_integer_partition(-1))
    assert random_integer_partition(1) == [1]
    assert random_integer_partition(10, seed=[1, 3, 2, 1, 5, 1]
            ) == [5, 2, 1, 1, 1]


def test_integer_partition_as_dict_aliasing():
    p = IntegerPartition([1, 1, 2, 3])
    d = p.as_dict()
    d[99] = 1
    assert p.as_dict() == {1: 2, 2: 1, 3: 1}


def test_rgs():
    raises(ValueError, lambda: RGS_unrank(-1, 3))
    raises(ValueError, lambda: RGS_unrank(3, 0))
    raises(ValueError, lambda: RGS_unrank(10, 1))

    raises(ValueError, lambda: Partition.from_rgs(list(range(3)), list(range(2))))
    raises(ValueError, lambda: Partition.from_rgs(list(range(1, 3)), list(range(2))))
    assert RGS_enum(-1) == 0
    assert RGS_enum(1) == 1
    assert RGS_unrank(7, 5) == [0, 0, 1, 0, 2]
    assert RGS_unrank(23, 14) == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 2]
    assert RGS_rank(RGS_unrank(40, 100)) == 40

def test_ordered_partition_9608():
    a = Partition([1, 2, 3], [4])
    b = Partition([1, 2], [3, 4])
    assert list(ordered([a, b], Set._infimum_key))


def test_integer_partition_dominance():
    p1 = IntegerPartition([4, 1, 1])
    p2 = IntegerPartition([3, 2, 1])
    p3 = IntegerPartition([3, 3])
    p4 = IntegerPartition([2, 2, 2])

    # Basic dominance relations for n=6
    assert p1.dominates(p2) is True
    assert p2.is_dominated_by(p1) is True
    assert p2.dominates(p1) is False
    assert p1.is_dominated_by(p2) is False

    assert p3.dominates(p4) is True
    assert p4.dominates(p3) is False

    # Incomparable partitions for n=6
    assert p1.dominates(p3) is False
    assert p3.dominates(p1) is False
    assert p1.is_dominated_by(p3) is False
    assert p3.is_dominated_by(p1) is False

    # Top [n] dominates all, bottom [1^n] is dominated by all
    top = IntegerPartition([6])
    bottom = IntegerPartition([1, 1, 1, 1, 1, 1])
    for p in [p1, p2, p3, p4, top, bottom]:
        assert top.dominates(p) is True
        assert p.is_dominated_by(top) is True
        assert p.dominates(bottom) is True
        assert bottom.is_dominated_by(p) is True

    # Reflexivity
    assert p1.dominates(p1) is True
    assert p1.is_dominated_by(p1) is True

    # Coercion with list / tuple / dict
    assert p1.dominates([3, 2, 1]) is True
    assert p1.dominates((3, 2, 1)) is True
    assert p1.dominates({3: 1, 2: 1, 1: 1}) is True
    assert p2.is_dominated_by([4, 1, 1]) is True

    # Different integer partitions are incomparable
    p_diff = IntegerPartition([5])
    assert p1.dominates(p_diff) is False
    assert p_diff.dominates(p1) is False
    assert p1.is_dominated_by(p_diff) is False

    # Conjugation Duality Theorem: lambda >= mu <=> mu' >= lambda'
    for n in range(1, 8):
        parts = [IntegerPartition(p) for p in partitions(n)]
        for a in parts:
            for b in parts:
                ca = IntegerPartition(a.conjugate)
                cb = IntegerPartition(b.conjugate)
                assert a.dominates(b) == cb.dominates(ca)
                assert a.is_dominated_by(b) == b.dominates(a)
                # Antisymmetry
                if a.dominates(b) and b.dominates(a):
                    assert a.partition == b.partition
