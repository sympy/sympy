from sympy.combinatorics.partitions import (Partition, IntegerPartition,
                                            RGS_generalized,
                                            RGS_enum, RGS_unrank, RGS_rank)
from sympy.utilities.pytest import raises
from sympy.utilities.iterables import default_sort_key, partitions

def test_partition():
    from sympy.abc import x

    raises(ValueError, lambda: Partition(range(3)))
    raises(ValueError, lambda: Partition([[1, 1, 2]]))

    a = Partition([[1,2,3], [4]])
    b = Partition([[1,2], [3,4]])
    c = Partition([[x]])
    l = [a, b, c]
    l.sort(key=default_sort_key)
    assert l == [c, a, b]
    l.sort(key=lambda w: default_sort_key(w, order='rev-lex'))
    assert l == [c, a, b]

    assert (a == b) == False
    assert a < b
    assert (a > b) == False
    assert a != b

    assert (a + b).as_list() == [[1, 2], [3], [4]]
    assert (a - b).as_list() == [[1], [2], [3, 4]]
    assert (a + 2).as_list() == [[1, 2], [3, 4]]
    assert (b - 1).as_list() == [[1, 2, 4], [3]]

    assert (a - 1).as_list() == [[1, 2, 3, 4]]
    assert (a + 1).as_list() == [[1, 2, 4], [3]]
    assert (b + 1).as_list() == [[1, 2], [3], [4]]

    assert a.rank == 1
    assert b.rank == 3

    assert a.RGS == [0, 0, 0, 1]
    assert b.RGS == [0, 0, 1, 1]


def test_integer_partition():
    raises(ValueError, lambda: IntegerPartition(range(3)))
    raises(ValueError, lambda: IntegerPartition(100, range(1, 3)))
    a = IntegerPartition(8, [1,3,4])
    b = IntegerPartition(8, [1,1,2,4])
    c = IntegerPartition([1,3,4])
    d = IntegerPartition(8, {1:3, 3:1, 2:1})
    assert a == c
    assert a.integer == d.integer
    assert a.conjugate == [3, 2, 2, 1]
    assert (a == b) == False
    assert a > b
    assert (a < b) == False
    assert a != b

    from sympy.utilities.iterables import partitions
    for i in range(1, 11):
        next = set()
        prev = set()
        a = IntegerPartition([i])
        ans = set([IntegerPartition(p).partition for p in partitions(i)])
        n = len(ans)
        for j in range(n):
            next.add(a.partition)
            a = a.next_lex()
            IntegerPartition(i, a.partition) # check it
        for j in range(n):
            prev.add(a.partition)
            a = a.prev_lex()
            IntegerPartition(i, a.partition) # check it
        assert next == ans
        assert prev == ans

def test_rgs():
    assert RGS_unrank(7, 5) == [0, 0, 1, 0, 2]
    assert RGS_unrank(23, 14) == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 2]
    assert RGS_rank(RGS_unrank(40, 100)) == 40
