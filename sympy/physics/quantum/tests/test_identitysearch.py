from sympy.physics.quantum.gate import (X, Y, Z, H, S, T, CNOT,
        IdentityGate, CGate, gate_simp)
from sympy.physics.quantum.identitysearch import *
from sympy.physics.quantum.dagger import Dagger

def test_generate_equivalent_ids():
    x = X(0)
    y = Y(0)
    z = Z(0)

    assert generate_equivalent_ids(x) == set([(x,)])
    assert generate_equivalent_ids(x, y) == set([(x, y), (y, x)])

    gate_seq = (x, y, z)
    gate_rules = set([(x, y, z), (y, z, x), (z, x, y), (z, y, x),
                      (y, x, z), (x, z, y)])
    assert generate_equivalent_ids(*gate_seq) == gate_rules

    h = H(0)
    gate_seq = (x, y, z, h)
    gate_rules = set([(x, y, z, h), (y, z, h, x),
                      (h, x, y, z), (h, z, y, x),
                      (z, y, x, h), (y, x, h, z),
                      (z, h, x, y) ,(x, h, z, y)])
    assert generate_equivalent_ids(*gate_seq) == gate_rules

    gate_seq = (x, y, x, y)
    gate_rules = set([(x, y, x, y), (y, x, y, x)])
    assert generate_equivalent_ids(*gate_seq) == gate_rules

    cgate_y = CGate((1,), y)
    gate_seq = (y, cgate_y, y, cgate_y)
    gate_rules = set([(y, cgate_y, y, cgate_y), (cgate_y, y, cgate_y, y)])
    assert generate_equivalent_ids(*gate_seq) == gate_rules

    cnot = CNOT(1,0)
    cgate_z = CGate((0,), Z(1))
    gate_seq = (cnot, h, cgate_z, h)
    gate_rules = set([(cnot, h, cgate_z, h), (h, cgate_z, h, cnot),
                      (h, cnot, h, cgate_z), (cgate_z, h, cnot, h)])
    assert generate_equivalent_ids(*gate_seq) == gate_rules

def test_is_scalar_matrix():
    numqubits = 2
    id_only = False

    id_gate = (IdentityGate(1),)
    assert is_scalar_matrix(id_gate, numqubits, id_only) == True
    assert is_scalar_sparse_matrix(id_gate, numqubits, id_only) == True

    x0 = X(0)
    xx_circuit = (x0, x0)
    assert is_scalar_matrix(xx_circuit, numqubits, id_only) == True
    assert is_scalar_sparse_matrix(xx_circuit, numqubits, id_only) == True

    x1 = X(1)
    y1 = Y(1)
    xy_circuit = (x1, y1)
    assert is_scalar_matrix(xy_circuit, numqubits, id_only) == False
    assert is_scalar_sparse_matrix(xy_circuit, numqubits, id_only) == False

    z1 = Z(1)
    xyz_circuit = (x1, y1, z1)
    assert is_scalar_matrix(xyz_circuit, numqubits, id_only) == True
    assert is_scalar_sparse_matrix(xyz_circuit, numqubits, id_only) == True

    cnot = CNOT(1,0)
    cnot_circuit = (cnot, cnot)
    assert is_scalar_matrix(cnot_circuit, numqubits, id_only) == True
    assert is_scalar_sparse_matrix(cnot_circuit, numqubits, id_only) == True

    h = H(0)
    hh_circuit = (h, h)
    assert is_scalar_matrix(hh_circuit, numqubits, id_only) == True
    assert is_scalar_sparse_matrix(hh_circuit, numqubits, id_only) == True

    id_only = True
    assert is_scalar_matrix(xyz_circuit, numqubits, id_only) == False
    assert is_scalar_matrix(cnot_circuit, numqubits, id_only) == True
    assert is_scalar_matrix(hh_circuit, numqubits, id_only) == True

    assert is_scalar_sparse_matrix(xyz_circuit, numqubits, id_only) == False
    assert is_scalar_sparse_matrix(cnot_circuit, numqubits, id_only) == True
    assert is_scalar_sparse_matrix(hh_circuit, numqubits, id_only) == True

def test_is_degenerate():
    x = X(0)
    y = Y(0)
    z = Z(0)

    gate_id = GateIdentity(x, y, z)
    ids = set([gate_id])

    another_id = (z, y, x)
    assert is_degenerate(ids, another_id) == True

def test_is_reducible():
    nqubits = 2

    x = X(0)
    y = Y(0)
    z = Z(0)

    circuit = (x, y, y)
    assert is_reducible(circuit, nqubits, 1, 3) == True

    circuit = (x, y, x)
    assert is_reducible(circuit, nqubits, 1, 3) == False

    circuit = (x, y, y, x)
    assert is_reducible(circuit, nqubits, 0, 4) == True

    circuit = (x, y, y, x)
    assert is_reducible(circuit, nqubits, 1, 3) == True

    circuit = (x, y, z, y, y)
    assert is_reducible(circuit, nqubits, 1, 5) == True

def test_bfs_identity_search():
    assert bfs_identity_search([], 1) == set()

    x = X(0)
    y = Y(0)
    z = Z(0)

    gate_list = [x]
    id_set = set([GateIdentity(x, x)])
    assert bfs_identity_search(gate_list, 1, max_depth=2) == id_set

    # Set should not contain degenerate quantum circuits
    gate_list = [x, y, z]
    id_set = set([GateIdentity(x, x), 
                  GateIdentity(y, y),
                  GateIdentity(z, z),
                  GateIdentity(x, y, z)])
    assert bfs_identity_search(gate_list, 1) == id_set

    id_set = set([GateIdentity(x, x),
                  GateIdentity(y, y),
                  GateIdentity(z, z),
                  GateIdentity(x, y, z),
                  GateIdentity(x, y, x, y),
                  GateIdentity(x, z, x, z),
                  GateIdentity(y, z, y, z)])
    assert bfs_identity_search(gate_list, 1, max_depth=4) == id_set
    assert bfs_identity_search(gate_list, 1, max_depth=5) == id_set

    h = H(0)
    gate_list = [x, y, z, h]
    id_set = set([GateIdentity(x, x),
                  GateIdentity(y, y),
                  GateIdentity(z, z),
                  GateIdentity(h, h),
                  GateIdentity(x, y, z),
                  GateIdentity(x, y, x, y),
                  GateIdentity(x, z, x, z),
                  GateIdentity(x, h, z, h),
                  GateIdentity(y, z, y, z),
                  GateIdentity(y, h, y, h)])
    assert bfs_identity_search(gate_list, 1) == id_set

    id_set = set([GateIdentity(x, x),
                  GateIdentity(y, y),
                  GateIdentity(z, z),
                  GateIdentity(h, h)])
    assert id_set == bfs_identity_search(gate_list, 1, max_depth=3,
                                         identity_only=True)

    id_set = set([GateIdentity(x, x),
                  GateIdentity(y, y),
                  GateIdentity(z, z),
                  GateIdentity(h, h),
                  GateIdentity(x, y, z),
                  GateIdentity(x, y, x, y),
                  GateIdentity(x, z, x, z),
                  GateIdentity(x, h, z, h),
                  GateIdentity(y, z, y, z),
                  GateIdentity(y, h, y, h),
                  GateIdentity(x, y, h, x, h),
                  GateIdentity(x, z, h, y, h),
                  GateIdentity(y, z, h, z, h)])
    assert bfs_identity_search(gate_list, 1, max_depth=5) == id_set

    cnot = CNOT(1,0)
    gate_list = [x, cnot]
    id_set = set([GateIdentity(x, x),
                  GateIdentity(cnot, cnot),
                  GateIdentity(x, cnot, x, cnot)])
    assert bfs_identity_search(gate_list, 2, max_depth=4) == id_set

    cgate_x = CGate((1,), x)
    gate_list = [x, cgate_x]
    id_set = set([GateIdentity(x, x),
                  GateIdentity(cgate_x, cgate_x),
                  GateIdentity(x, cgate_x, x, cgate_x)])
    assert bfs_identity_search(gate_list, 2, max_depth=4) == id_set

    cgate_z = CGate((0,), Z(1))
    gate_list = [cnot, cgate_z, h]
    id_set = set([GateIdentity(h, h),
                  GateIdentity(cgate_z, cgate_z),
                  GateIdentity(cnot, cnot),
                  GateIdentity(cnot, h, cgate_z, h)])
    assert bfs_identity_search(gate_list, 2, max_depth=4) == id_set
