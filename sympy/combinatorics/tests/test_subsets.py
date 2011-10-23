from sympy.combinatorics import Subset

def test_subset():
    a = Subset(['c','d'],['a','b','c','d'])
    assert a.next_binary() == Subset(['b'], ['a','b','c','d'])
    assert a.prev_binary() == Subset(['c'], ['a','b','c','d'])
    assert a.next_graycode() == Subset(['a','c','d'], ['a','b','c','d'])
    assert a.prev_graycode() == Subset(['c'], ['a','b','c','d'])
    assert a.rank_binary == 3
    assert a.rank_lexicographic == 14
    assert a.rank_graycode == 8
    assert a.cardinality == 16

    a = Subset([2,5,7],[1,2,3,4,5,6,7])
    assert a.next_binary() == Subset([2, 5, 6],[1,2,3,4,5,6,7])
    assert a.prev_binary() == Subset([2, 5],[1,2,3,4,5,6,7])
    assert a.next_graycode() == Subset([2, 3, 5, 7],[1,2,3,4,5,6,7])
    assert a.prev_graycode() == Subset([1, 2, 5, 7],[1,2,3,4,5,6,7])
    assert a.rank_binary == 37
    assert a.rank_lexicographic == 93
    assert a.rank_graycode == 99
    assert a.cardinality == 128

    superset = ['a','b','c','d']
    assert Subset.unrank_binary(4,superset).rank_binary == 4
    assert Subset.unrank_graycode_subset(10,superset).rank_graycode == 10

    superset = [1,2,3,4,5,6,7,8,9]
    assert Subset.unrank_binary(33,superset).rank_binary == 33
    assert Subset.unrank_graycode_subset(25,superset).rank_graycode == 25
