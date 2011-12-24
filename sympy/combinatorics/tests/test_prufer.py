from sympy.combinatorics.prufer import Prufer
from sympy.utilities.pytest import raises

def test_prufer():
    # number of nodes is optional
    assert Prufer([[0, 1], [0, 2], [0, 3], [0, 4]], 5).nodes == 5
    assert Prufer([[0, 1], [0, 2], [0, 3], [0, 4]]).nodes == 5

    a = Prufer([[0, 1], [0, 2], [0, 3], [0, 4]])
    assert a.rank == 0
    assert a.nodes == 5
    assert a.prufer_repr == [0, 0, 0]

    a = Prufer([[2, 4], [1, 4], [1, 3], [0, 5], [0, 4]])
    assert a.rank == 64
    assert a.nodes == 6
    assert a.tree_repr == [[2, 4], [1, 4], [1, 3], [0, 5], [0, 4]]
    assert a.prufer_repr == [0, 1, 4, 4]

    raises(ValueError, 'Prufer([[1, 2], [3, 4]])') # 0 is missing
