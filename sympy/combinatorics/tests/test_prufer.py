from sympy.combinatorics.prufer import Prufer

def test_prufer():
    a = Prufer([[0, 1], [0, 2], [0, 3], [0, 4]], 5)
    assert a.rank == 0
    assert a.nodes == 5
    assert a.prufer_repr == [0, 0, 0]

    a = Prufer([[2, 4], [1, 4], [1, 3], [0, 5], [0, 4]], 6)
    assert a.rank == 1008
    assert a.nodes == 6
    assert a.tree_repr == [[2, 4], [1, 4], [1, 3], [0, 5], [0, 4]]
    assert a.prufer_repr == [4, 4, 0, 0]
