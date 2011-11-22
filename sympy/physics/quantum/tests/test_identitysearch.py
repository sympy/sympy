from sympy.physics.quantum.gate import (X, Y, Z, H, S, T, CNOT,
        IdentityGate, gate_simp)
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.identitysearch import *
from sympy.physics.quantum.dagger import Dagger

def test_generate_gate_rules():
    x = X(0)
    y = Y(0)
    z = Z(0)

    assert generate_gate_rules(x) == [(x,)]
    assert generate_gate_rules(x, y) == [(x, y), (y, x)]

    gate_seq = (x, y, z)
    gate_rules = [(x, y, z), (y, z, x), (z, x, y), (z, y, x),
                  (y, x, z), (x, z, y)]
    assert generate_gate_rules(*gate_seq) == gate_rules

    h = H(0)
    h_dag = Dagger(h)
    gate_seq = (x, y, z, h)
    gate_rules = [(x, y, z, h), (y, z, h, x),
                  (h, x, y, z), (h_dag, z, y, x),
                  (z, y, x, h_dag), (y, x, h_dag, z),
                  (z, h, x, y) ,(x, h_dag, z, y)]
    assert generate_gate_rules(*gate_seq) == gate_rules

def test_is_scalar_matrix():
    numqubits = 2

    id_gate = (IdentityGate(1),)
    assert is_scalar_matrix(id_gate, numqubits) == True

    x0 = X(0)
    xx_circuit = (x0, x0)
    assert is_scalar_matrix(xx_circuit, numqubits) == True

    x1 = X(1)
    y1 = Y(1)
    xy_circuit = (x1, y1)
    assert is_scalar_matrix(xy_circuit, numqubits) == False

    z1 = Z(1)
    xyz_circuit = (x1, y1, z1)
    assert is_scalar_matrix(xyz_circuit, numqubits) == True

    cnot = CNOT(1,0)
    cnot_circuit = (cnot, cnot)
    assert is_scalar_matrix(cnot_circuit, numqubits) == True

    h = H(0)
    hh_circuit = (h, h)
    assert is_scalar_matrix(hh_circuit, numqubits) == True

def test_is_degenerate():
    gate_id = GateIdentity(X(0), Y(0), Z(0))
    ids = set([gate_id])

    another_id = (Z(0), Y(0), X(0))
    assert is_degenerate(ids, another_id) == True

def test_is_reducible():
    nqubits = 2

    circuit = (X(0), Y(0), Y(0))
    assert is_reducible(circuit, nqubits, 1, 3) == True

    circuit = (X(0), Y(0), X(0))
    assert is_reducible(circuit, nqubits, 1, 3) == False

    circuit = (X(0), Y(0), Y(0), X(0))
    assert is_reducible(circuit, nqubits, 0, 4) == True

    circuit = (X(0), Y(0), Y(0), X(0))
    assert is_reducible(circuit, nqubits, 1, 3) == True

def test_bfs_identity_search():
    assert bfs_identity_search([], 1) == set()

    gate_list = [X(0)]
    id_set = set([GateIdentity(X(0), X(0))])
    assert bfs_identity_search(gate_list, 1, 2) == id_set

    # Set should not contain degenerate quantum circuits
    gate_list = [X(0), Y(0), Z(0)]
    id_set = set([GateIdentity(X(0), X(0)), 
                  GateIdentity(Y(0), Y(0)),
                  GateIdentity(Z(0), Z(0)),
                  GateIdentity(X(0), Y(0), Z(0))])
    assert bfs_identity_search(gate_list, 1) == id_set

    id_set = set([GateIdentity(X(0), X(0)),
                  GateIdentity(Y(0), Y(0)),
                  GateIdentity(Z(0), Z(0)),
                  GateIdentity(X(0), Y(0), Z(0)),
                  GateIdentity(X(0), Y(0), X(0), Y(0)),
                  GateIdentity(X(0), Z(0), X(0), Z(0)),
                  GateIdentity(Y(0), Z(0), Y(0), Z(0))])
    assert bfs_identity_search(gate_list, 1, 4) == id_set
    assert bfs_identity_search(gate_list, 1, 5) == id_set

    gate_list = [X(0), Y(0), Z(0), H(0)]
    id_set = set([GateIdentity(X(0), X(0)),
                  GateIdentity(Y(0), Y(0)),
                  GateIdentity(Z(0), Z(0)),
                  GateIdentity(H(0), H(0)),
                  GateIdentity(X(0), Y(0), Z(0)),
                  GateIdentity(X(0), Y(0), X(0), Y(0)),
                  GateIdentity(X(0), Z(0), X(0), Z(0)),
                  GateIdentity(X(0), H(0), Z(0), H(0)),
                  GateIdentity(Y(0), Z(0), Y(0), Z(0)),
                  GateIdentity(Y(0), H(0), Y(0), H(0))])
    assert bfs_identity_search(gate_list, 1) == id_set

    id_set = set([GateIdentity(X(0), X(0)),
                  GateIdentity(Y(0), Y(0)),
                  GateIdentity(Z(0), Z(0)),
                  GateIdentity(H(0), H(0)),
                  GateIdentity(X(0), Y(0), Z(0)),
                  GateIdentity(X(0), Y(0), X(0), Y(0)),
                  GateIdentity(X(0), Z(0), X(0), Z(0)),
                  GateIdentity(X(0), H(0), Z(0), H(0)),
                  GateIdentity(Y(0), Z(0), Y(0), Z(0)),
                  GateIdentity(Y(0), H(0), Y(0), H(0)),
                  GateIdentity(Z(0), H(0), X(0), H(0))])
    #assert bfs_identity_search(gate_list, 1, 5) == id_set

