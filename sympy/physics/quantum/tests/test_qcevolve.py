from sympy.physics.quantum.qcevolve import *
from sympy.physics.quantum.gate import (X, Y, Z, H, S, T, CNOT,
        IdentityGate, CGate, gate_simp)
from sympy.physics.quantum.identitysearch import (
        GateIdentity, bfs_identity_search)

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

def test_find_subcircuit():
    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    x1 = X(1)
    y1 = Y(1)

    circuit = (x, y, z)

    assert find_subcircuit(circuit, (x,)) == 0
    assert find_subcircuit(circuit, (x1,)) == -1
    assert find_subcircuit(circuit, (y,)) == 1
    assert find_subcircuit(circuit, (h,)) == -1
    assert find_subcircuit(circuit, (x, h)) == -1
    assert find_subcircuit(circuit, (x, y, z)) == 0
    assert find_subcircuit(circuit, (y, z)) == 1
    assert find_subcircuit(circuit, (x, y, z, h)) == -1
    assert find_subcircuit(circuit, (z, y, x)) == -1
    assert find_subcircuit(circuit, (x,), start=2, end=1) == -1

    circuit = (x, y, x, y, z)
    assert find_subcircuit(circuit, (x, y, z)) == 2
    assert find_subcircuit(circuit, (x,), start=1) == 2
    assert find_subcircuit(circuit, (x, y), start=1, end=2) == -1
    assert find_subcircuit(circuit, (x, y), start=1, end=3) == -1
    assert find_subcircuit(circuit, (x, y), start=1, end=4) == 2
    assert find_subcircuit(circuit, (x, y), start=2, end=4) == 2

    circuit = (x, y, z, x1, x, y, z, h, x, y, x1,
               x, y, z, h, y1, h)
    assert find_subcircuit(circuit, (x, y, z, h, y1)) == 11

def test_qc_reduce():
    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    cnot = CNOT(1,0)

    gate_list = [x, y, z]
    id_set = bfs_identity_search(gate_list, 1, max_depth=4)

    circuit = (z, y, x, x)
    # assert qc_reduce(circuit, list(id_set)) == (x,)

    circuit = (x, y, x, y, z)
    #assert qc_reduce(circuit, list(id_set)) == (x, y)

    circuit = (x, y, z, y, x)
    #assert qc_reduce(circuit, list(id_set)) == (y, x)

    circuit = (x, y, z, x, y, x, y)
    #assert qc_reduce(circuit, list(id_set), homogeneous=False) == ()
