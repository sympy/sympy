from sympy.combinatorics.partitions import (Partition, IntegerPartition,
                                            RGS_generalized,
                                            RGS_enum, RGS_unrank, RGS_rank)

def test_partition():
    a = Partition([[1,2,3], [4]])
    b = Partition([[1,2], [3,4]])

    assert (a == b) == False
    assert a < b
    assert (a > b) == False
    assert a != b

    assert (a + b).as_list() == [[1, 2], [3], [4]]
    assert (a - b).as_list() == [[1], [2], [3, 4]]
    assert (a + 2).as_list() == [[1, 2], [3, 4]]
    assert (b - 1).as_list() == [[1, 2, 4], [3]]

    assert a.previous().as_list() == [[1, 2, 3, 4]]
    assert a.next().as_list() == [[1, 2, 4], [3]]
    assert b.previous().as_list() == (b - 1).as_list()
    assert b.next().as_list() == [[1, 2], [3], [4]]

    assert a.rank == 1
    assert b.rank == 3

    assert a.RGS == [0, 0, 0, 1]
    assert b.RGS == [0, 0, 1, 1]


def test_integer_partition():
    a = IntegerPartition(8, [1,3,4])
    b = IntegerPartition(8, [1,1,2,4])
    assert a.conjugate == [3, 2, 2, 1]

    assert (a == b) == False
    assert a > b
    assert (a < b) == False
    assert a != b

def test_rgs():
    assert RGS_unrank(7, 5) == [0, 0, 1, 0, 2]
    assert RGS_unrank(23, 14) == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 2]
    assert RGS_rank(RGS_unrank(40, 100)) == 40
