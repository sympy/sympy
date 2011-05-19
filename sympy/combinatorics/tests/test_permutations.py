from sympy.combinatorics.permutations import Permutation

def test_Permutation():
    p = Permutation([2,5,1,6,3,0,4])
    q = Permutation([[1,4,5],[2,0,6],[3]])

    assert p.is_ArrayForm
    assert q.is_CyclicForm

    assert p*q == Permutation([4, 6, 1, 2, 5, 3, 0])
    assert q*p == Permutation([6, 5, 3, 0, 2, 4, 1])

    assert q.to_array() == Permutation([3, 1, 4, 5, 0, 6, 2])
    assert p.to_cycles() == Permutation([[3, 6, 4], [0, 2, 1, 5]])

    assert p**13 == p
    assert q**2 == Permutation([[3, 6, 4], [1], [0, 5, 2]])

    assert len(p.atoms()) == 7
    assert q.atoms() == set([0, 1, 2, 3, 4, 5, 6])

