from sympy.combinatorics import Permutation
from sympy.core.symbol import Symbol, symbols
from sympy.matrices import Matrix
from sympy.matrices.expressions import PermutationMatrix, BlockDiagMatrix


def test_connected_components():
    kp, ki, kd = symbols('kp ki kd')
    Ix, Iy, Iz = symbols('Ix Iy Iz')
    d = Symbol('d')

    M = Matrix([
        [1 + d*kd/Ix,           0,           0, 0, d*kp/Ix,       0,       0, 0, 0, 0, d*ki/Ix,       0,       0],
        [          0, 1 + d*kd/Iy,           0, 0,       0, d*kp/Iy,       0, 0, 0, 0,       0, d*ki/Iy,       0],
        [          0,           0, 1 + d*kd/Iz, 0,       0,       0, d*kp/Iz, 0, 0, 0,       0,       0, d*ki/Iz],
        [          0,           0,           0, 1,       0,       0,       0, 0, 0, 0,       0,       0,       0],
        [        d/2,           0,           0, 0,       1,       0,       0, 0, 0, 0,       0,       0,       0],
        [          0,         d/2,           0, 0,       0,       1,       0, 0, 0, 0,       0,       0,       0],
        [          0,           0,         d/2, 0,       0,       0,       1, 0, 0, 0,       0,       0,       0],
        [      -d*kd,           0,           0, 0,   -d*kp,       0,       0, 1, 0, 0,   -d*ki,       0,       0],
        [          0,       -d*kd,           0, 0,       0,   -d*kp,       0, 0, 1, 0,       0,   -d*ki,       0],
        [          0,           0,       -d*kd, 0,       0,       0,   -d*kp, 0, 0, 1,       0,       0,   -d*ki],
        [          0,           0,           0, 0,       d,       0,       0, 0, 0, 0,       1,       0,       0],
        [          0,           0,           0, 0,       0,       d,       0, 0, 0, 0,       0,       1,       0],
        [          0,           0,           0, 0,       0,       0,       d, 0, 0, 0,       0,       0,       1]])
    cc = M.connected_components()
    assert cc == [[0, 4, 7, 10], [1, 5, 8, 11], [2, 6, 9, 12], [3]]

    P, B = M.connected_components_decomposition()
    p = Permutation([0, 4, 7, 10, 1, 5, 8, 11, 2, 6, 9, 12, 3])
    assert P == PermutationMatrix(p)

    B0 = Matrix([
        [1 + d*kd/Ix, d*kp/Ix, 0, d*ki/Ix],
        [        d/2,       1, 0,       0],
        [      -d*kd,   -d*kp, 1,   -d*ki],
        [          0,       d, 0,       1]])
    B1 = Matrix([
        [1 + d*kd/Iy, d*kp/Iy, 0, d*ki/Iy],
        [        d/2,       1, 0,       0],
        [      -d*kd,   -d*kp, 1,   -d*ki],
        [          0,       d, 0,       1]])
    B2 = Matrix([
        [1 + d*kd/Iz, d*kp/Iz, 0, d*ki/Iz],
        [        d/2,       1, 0,       0],
        [      -d*kd,   -d*kp, 1,   -d*ki],
        [          0,       d, 0,       1]])
    B3 = Matrix([[1]])
    assert B == BlockDiagMatrix(B0, B1, B2, B3)
