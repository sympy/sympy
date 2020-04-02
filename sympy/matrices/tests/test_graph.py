from sympy.combinatorics import Permutation
from sympy.core.symbol import symbols
from sympy.matrices import Matrix
from sympy.matrices.expressions import PermutationMatrix, BlockDiagMatrix


def test_connected_components():
    a, b, c, d, e, f, g, h, i, j, k, l, m = symbols('a:m')

    M = Matrix([
        [a, 0, 0, 0, b, 0, 0, 0, 0, 0, c, 0, 0],
        [0, d, 0, 0, 0, e, 0, 0, 0, 0, 0, f, 0],
        [0, 0, g, 0, 0, 0, h, 0, 0, 0, 0, 0, i],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [m, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, m, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, m, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [j, 0, 0, 0, k, 0, 0, 1, 0, 0, l, 0, 0],
        [0, j, 0, 0, 0, k, 0, 0, 1, 0, 0, l, 0],
        [0, 0, j, 0, 0, 0, k, 0, 0, 1, 0, 0, l],
        [0, 0, 0, 0, d, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, d, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, d, 0, 0, 0, 0, 0, 1]])
    cc = M.connected_components()
    assert cc == [[0, 4, 7, 10], [1, 5, 8, 11], [2, 6, 9, 12], [3]]

    P, B = M.connected_components_decomposition()
    p = Permutation([0, 4, 7, 10, 1, 5, 8, 11, 2, 6, 9, 12, 3])
    assert P == PermutationMatrix(p)

    B0 = Matrix([
        [a, b, 0, c],
        [m, 1, 0, 0],
        [j, k, 1, l],
        [0, d, 0, 1]])
    B1 = Matrix([
        [d, e, 0, f],
        [m, 1, 0, 0],
        [j, k, 1, l],
        [0, d, 0, 1]])
    B2 = Matrix([
        [g, h, 0, i],
        [m, 1, 0, 0],
        [j, k, 1, l],
        [0, d, 0, 1]])
    B3 = Matrix([[1]])
    assert B == BlockDiagMatrix(B0, B1, B2, B3)
