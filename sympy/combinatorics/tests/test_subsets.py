from sympy.combinatorics import (Subset, unrank_graycode_subset,
                                 unrank_binary)

def test_subset():
    a = Subset(['c','d'],['a','b','c','d'])
    assert a.next_binary() == ['b']
    assert a.prev_binary() == ['c']
    assert a.next_graycode() == ['c']
    assert a.prev_graycode() == ['c']
    assert a.rank_binary == 3
    assert a.rank_lexicographic == 14
    assert a.rank_graycode == 8
    assert a.cardinality == 16

    a = Subset([2,5,7],[1,2,3,4,5,6,7])
    assert a.next_binary() == [2, 5, 6]
    assert a.prev_binary() == [2, 5]
    assert a.next_graycode() == [2, 5, 6, 7]
    assert a.prev_graycode() == [1, 2, 5, 7]
    assert a.rank_binary == 37
    assert a.rank_lexicographic == 93
    assert a.rank_graycode == 99
    assert a.cardinality == 128

    superset = ['a','b','c','d']
    assert Subset(unrank_binary(4,superset),
                  superset).rank_binary == 4
    assert Subset(unrank_graycode_subset(10,superset),
                  superset).rank_graycode == 10

    superset = [1,2,3,4,5,6,7,8,9]
    assert Subset(unrank_binary(33,superset),
                  superset).rank_binary == 33
    assert Subset(unrank_graycode_subset(25,superset),
                  superset).rank_graycode == 25
