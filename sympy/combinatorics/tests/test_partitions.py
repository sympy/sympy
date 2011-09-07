from sympy.combinatorics.partitions import (Partition, IntegerPartition,
                                            from_RGS, RGS_generalized,
                                            RGS_enum, RGS_unrank, RGS_rank)

def test_partition():
    a = Partition([[1,2,3],[4]],[1,2,3,4])
    b = Partition([[1,2], [3,4]], [1,2,3,4])

    assert (a == b) == False
    assert a > b
    assert (a < b) == False
    assert a != b

    assert (a + b).partition == [[1, 2], [3], [4]]
    assert (a - b).partition == [[1], [2], [3, 4]]
    assert (a + 2).partition == [[1, 2], [3, 4]]
    assert (b - 1).partition == [[1, 2, 4], [3]]

    assert a.previous().partition == [[1, 2, 3, 4]]
    assert a.next().partition ==  [[1, 2, 4], [3]]
    assert b.previous().partition == (b - 1).partition
    assert b.next().partition == [[1, 2], [3], [4]]

    assert a.rank == 1
    assert b.rank == 3

    assert a.RGS == [0, 0, 0, 1]
    assert b.RGS == [0, 0, 1, 1]

def test_integer_partition():
    a = IntegerPartition([1,3,4], 8)
    b = IntegerPartition([1,1,2,4], 8)
    assert a.conjugate == [3, 2, 2, 1]
    assert a.conjugate_partition == a.conjugate

    assert (a == b) == False
    assert a > b
    assert (a < b) == False
    assert a != b

def test_rgs():
    assert RGS_unrank(7, 5) == [0, 0, 1, 0, 2]
    assert RGS_unrank(23, 14) == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 2]
    assert RGS_rank(RGS_unrank(40, 100)) == 40
