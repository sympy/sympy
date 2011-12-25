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

    assert Prufer.edges([0, 1, 2, 3], [1, 4, 5], [1, 4, 6]) == \
        [[0, 1], [1, 2], [4, 6], [4, 5], [1, 4], [2, 3]]

    # accept iterables but convert to list of lists
    tree = [(0, 1), (1, 5), (0, 3), (0, 2), (2, 6), (4, 7), (2, 4)]
    Prufer(tree).tree_repr == tree
    Prufer(set(tree)).tree_repr == tree

    raises(ValueError, 'Prufer([[1, 2], [3, 4]])') # 0 is missing
