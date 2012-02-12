from sympy.physics.quantum.qcevolve import *
from sympy.physics.quantum.gate import (X, Y, Z, H, S, T, CNOT,
        CGate)
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

def test_qc_remove_subcircuit():
    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    cnot = CNOT(1,0)
    cgate_z = CGate((0,), Z(1))

    # Standard cases
    circuit = (z, y, x, x)
    remove = (z, y, x)
    assert qc_remove_subcircuit(circuit, remove) == (x,)
    assert qc_remove_subcircuit(circuit, remove + (x,)) == ()
    assert qc_remove_subcircuit(circuit, remove, pos=1) == circuit
    assert qc_remove_subcircuit(circuit, remove, pos=0) == (x,)
    assert qc_remove_subcircuit(circuit, (x, x), pos=2) == (z, y)
    assert qc_remove_subcircuit(circuit, (h,)) == circuit

    circuit = (x, y, x, y, z)
    remove = (x, y, z)
    assert qc_remove_subcircuit(circuit, remove) == (x, y)
    remove = (x, y, x, y)
    assert qc_remove_subcircuit(circuit, remove) == (z,)

    circuit = (x, h, cgate_z, h, cnot)
    remove = (x, h, cgate_z)
    assert qc_remove_subcircuit(circuit, remove, pos=-1) == (h, cnot)
    assert qc_remove_subcircuit(circuit, remove, pos=1) == circuit
    remove = (h, h)
    assert qc_remove_subcircuit(circuit, remove) == circuit
    remove = (h, cgate_z, h, cnot)
    assert qc_remove_subcircuit(circuit, remove) == (x,)

def test_qc_random_reduce():
    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    cnot = CNOT(1,0)
    cgate_z = CGate((0,), Z(1))

    seed = 1
    gate_list = [x, y, z]
    ids = list(bfs_identity_search(gate_list, 1, max_depth=4))
    ids.sort()

    circuit = (x, y, h, z, cnot)
    assert qc_random_reduce(circuit, []) == circuit
    assert qc_random_reduce(circuit, ids) == circuit

    circuit = (x, y, z, x, y, h)
    # With seed = 1, randint(0, 0) = 0
    assert qc_random_reduce(circuit, ids, seed=seed) == (x, y, h)

    circuit = (x, x, y, y, z, z)
    # randint(0, 2) = 0
    assert qc_random_reduce(circuit, ids, seed=seed) == (y, y, z, z)

    seed = 2
    # randint(0, 2) = 2
    assert qc_random_reduce(circuit, ids, seed=seed) == (x, x, y, y)

    gate_list = [x, y, z, h, cnot, cgate_z]
    ids = list(bfs_identity_search(gate_list, 2, max_depth=4))
    ids.sort()

    circuit = (x, y, z, y, h, y, h, cgate_z, h, cnot)
    expected = (x, y, z, y, h, y)
    # randint(0, 2) == 2
    assert qc_random_reduce(circuit, ids, seed=seed) == expected

    seed = 5
    expected = (x, y, z, cgate_z, h, cnot)
    # randint(0, 1) == 1
    assert qc_random_reduce(circuit, ids, seed=seed) == expected

def test_qc_random_insert():
    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    cnot = CNOT(1,0)
    cgate_z = CGate((0,), Z(1))

    seed = 1
    choices = [(x, x)]
    circuit = (y, y)
    # insert location: 0; 
    assert qc_random_insert(circuit, choices, seed=seed) == (x, x, y, y)

    seed = 8
    circuit = (x, y, z, h)
    choices = [(h, h), (x, y, z)]
    expected = (x, x, y, z, y, z, h)
    # insert location: 1; circuit choice: 1
    assert qc_random_insert(circuit, choices, seed=seed) == expected

    gate_list = [x, y, z, h, cnot, cgate_z]
    ids = list(bfs_identity_search(gate_list, 2, max_depth=4))
    ids.sort()

    collapse_func = lambda acc, an_id: acc + list(an_id.eq_identities)
    ids = reduce(collapse_func, ids, [])

    circuit = (x, y, h, cnot, cgate_z)
    expected = (x, cnot, h, cgate_z, h, y, h, cnot, cgate_z)
    # insert location: 1; circuit choice: 31
    actual = qc_random_insert(circuit, ids, seed=seed)
    assert actual == expected

