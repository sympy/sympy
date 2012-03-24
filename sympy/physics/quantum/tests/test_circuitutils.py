from sympy import Wild
from sympy.physics.quantum.circuitutils import *
from sympy.physics.quantum.gate import (X, Y, Z, H, S, T, CNOT,
        CGate)

def test_kmp_table():
    word = ('a', 'b', 'c', 'd', 'a', 'b', 'd')
    expected_table = [-1, 0, 0, 0, 0, 1, 2]
    assert expected_table == kmp_table(word)

    word = ('P', 'A', 'R', 'T', 'I', 'C', 'I', 'P', 'A', 'T', 'E', ' ',
            'I', 'N', ' ', 'P', 'A', 'R', 'A', 'C', 'H', 'U', 'T', 'E')
    expected_table = [-1, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0,
                       0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0]
    assert expected_table == kmp_table(word)

    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    word = (x, y, y, x, z)
    expected_table = [-1, 0, 0, 0, 1]
    assert expected_table == kmp_table(word)

    word = (x, x, y, h, z)
    expected_table = [-1, 0, 1, 0, 0]
    assert expected_table == kmp_table(word)

def test_find_subcircuit_with_seq():
    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    x1 = X(1)
    y1 = Y(1)

    i0 = Wild('i0')
    x_i0 = X(i0)
    y_i0 = Y(i0)
    z_i0 = Z(i0)

    circuit = (x, y, z)

    assert find_subcircuit_with_seq(circuit, (x,)) == 0
    assert find_subcircuit_with_seq(circuit, (x1,)) == -1
    assert find_subcircuit_with_seq(circuit, (y,)) == 1
    assert find_subcircuit_with_seq(circuit, (h,)) == -1
    assert find_subcircuit_with_seq(circuit, (x, h)) == -1
    assert find_subcircuit_with_seq(circuit, (x, y, z)) == 0
    assert find_subcircuit_with_seq(circuit, (y, z)) == 1
    assert find_subcircuit_with_seq(circuit, (x, y, z, h)) == -1
    assert find_subcircuit_with_seq(circuit, (z, y, x)) == -1
    assert find_subcircuit_with_seq(circuit, (x,), start=2, end=1) == -1

    circuit = (x, y, x, y, z)
    assert find_subcircuit_with_seq(circuit, (x, y, z)) == 2
    assert find_subcircuit_with_seq(circuit, (x,), start=1) == 2
    assert find_subcircuit_with_seq(circuit, (x, y), start=1, end=2) == -1
    assert find_subcircuit_with_seq(circuit, (x, y), start=1, end=3) == -1
    assert find_subcircuit_with_seq(circuit, (x, y), start=1, end=4) == 2
    assert find_subcircuit_with_seq(circuit, (x, y), start=2, end=4) == 2

    circuit = (x, y, z, x1, x, y, z, h, x, y, x1,
               x, y, z, h, y1, h)
    assert find_subcircuit_with_seq(circuit, (x, y, z, h, y1)) == 11

    circuit = (x, y, x_i0, y_i0, z_i0, z)
    assert find_subcircuit_with_seq(circuit, (x_i0, y_i0, z_i0)) == 2

def test_remove_subcircuit_with_seq():
    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    cnot = CNOT(1,0)
    cgate_z = CGate((0,), Z(1))

    # Standard cases
    circuit = (z, y, x, x)
    remove = (z, y, x)
    assert remove_subcircuit_with_seq(circuit, remove) == (x,)
    assert remove_subcircuit_with_seq(circuit, remove + (x,)) == ()
    assert remove_subcircuit_with_seq(circuit, remove, pos=1) == circuit
    assert remove_subcircuit_with_seq(circuit, remove, pos=0) == (x,)
    assert remove_subcircuit_with_seq(circuit, (x, x), pos=2) == (z, y)
    assert remove_subcircuit_with_seq(circuit, (h,)) == circuit

    circuit = (x, y, x, y, z)
    remove = (x, y, z)
    assert remove_subcircuit_with_seq(circuit, remove) == (x, y)
    remove = (x, y, x, y)
    assert remove_subcircuit_with_seq(circuit, remove) == (z,)

    circuit = (x, h, cgate_z, h, cnot)
    remove = (x, h, cgate_z)
    assert remove_subcircuit_with_seq(circuit, remove, pos=-1) == (h, cnot)
    assert remove_subcircuit_with_seq(circuit, remove, pos=1) == circuit
    remove = (h, h)
    assert remove_subcircuit_with_seq(circuit, remove) == circuit
    remove = (h, cgate_z, h, cnot)
    assert remove_subcircuit_with_seq(circuit, remove) == (x,)
