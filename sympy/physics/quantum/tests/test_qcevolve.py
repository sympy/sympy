from sympy.physics.quantum.qcevolve import *
from sympy.physics.quantum.gate import (X, Y, Z, H, S, T, CNOT,
        CGate)
from sympy.physics.quantum.identitysearch import (
        GateIdentity, bfs_identity_search)

def test_random_reduce():
    x = X(0)
    y = Y(0)
    z = Z(0)
    h = H(0)
    cnot = CNOT(1,0)
    cgate_z = CGate((0,), Z(1))

    seed = 1
    gate_list = [x, y, z]
    ids = list(bfs_identity_search(gate_list, 1, max_depth=4))

    circuit = (x, y, h, z, cnot)
    assert random_reduce(circuit, []) == circuit
    assert random_reduce(circuit, ids) == circuit

    circuit = (x, y, z, x, y, h)
    # With seed = 1, randint(0, 0) = 0
    actual = random_reduce(circuit, ids, seed=seed)
    assert actual == (x, y, h)

    circuit = (x, x, y, y, z, z)
    # randint(0, 2) = 0
    actual = random_reduce(circuit, ids, seed=seed)
    assert actual == (x, x, y, y)

    seed = 2
    # randint(0, 2) = 2
    assert random_reduce(circuit, ids, seed=seed) == (x, x, z, z)

    gate_list = [x, y, z, h, cnot, cgate_z]
    ids = list(bfs_identity_search(gate_list, 2, max_depth=4))

    circuit = (x, y, z, y, h, y, h, cgate_z, h, cnot)
    expected = (x, y, z, y, h, y)
    expected = (y, h, y, h, cgate_z, h, cnot)
    # randint(0, 2) == 2
    assert random_reduce(circuit, ids, seed=seed) == expected

    seed = 5
    expected = (x, y, z, cgate_z, h, cnot)
    expected = (x, y, z, y, h, y)
    # randint(0, 2) == 1
    assert random_reduce(circuit, ids, seed=seed) == expected

def test_random_insert():
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
    assert random_insert(circuit, choices, seed=seed) == (x, x, y, y)

    seed = 8
    circuit = (x, y, z, h)
    choices = [(h, h), (x, y, z)]
    expected = (x, x, y, z, y, z, h)
    # insert location: 1; circuit choice: 1
    assert random_insert(circuit, choices, seed=seed) == expected

    gate_list = [x, y, z, h, cnot, cgate_z]
    ids = list(bfs_identity_search(gate_list, 2, max_depth=4))

    collapse_func = lambda acc, an_id: acc + list(an_id.equivalent_ids)
    ids = reduce(collapse_func, ids, [])

    circuit = (x, y, h, cnot, cgate_z)
    expected = (x, y, z, y, z, y, h, cnot, cgate_z)
    # insert location: 1; circuit choice: 30
    actual = random_insert(circuit, ids, seed=seed)
    assert actual == expected
