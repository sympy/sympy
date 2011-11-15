from sympy.physics.quantum.gate import (X, Y, Z, H, S, T, CNOT,
        IdentityGate, gate_simp)
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.identitysearch import *

def permutations_recursive():
    assert permutations_recursive([], 0, 0, 1) == []
    assert permutations_recursive([1], 0, 0, -1) == []
    assert permutations_recursive([1], 0, -1, 1) == []
    assert permutations_recursive([1], -1, 0, 1) == []
    assert permutations_recursive([1], 0, 0, 2) == []
    assert permutations_recursive([1], 0, 2, 1) == []

def test_generate_gate_rules():
    assert generate_gate_rules((X(0),)) == [(X(0),)]
    assert generate_gate_rules((X(0), Y(0))) == [(X(0), Y(0)), (Y(0), X(0))]

    gate_seq = (X(0), Y(0), Z(0))
    gate_rules = [(X(0), Y(0), Z(0)), (Y(0), Z(0), X(0)),
                  (Z(0), X(0), Y(0)), (Z(0), Y(0), X(0)),
                  (Y(0), X(0), Z(0)), (X(0), Z(0), Y(0))]
    assert generate_gate_rules(gate_seq) == gate_rules


def test_is_scalar_matrix():
    numqubits = 2

    id_gate = IdentityGate(1)
    id_matrix = represent(id_gate, nqubits=numqubits, format='scipy.sparse')
    assert is_scalar_matrix(id_matrix) == True

    xy_matrix = represent(X(1)*Y(1), nqubits=numqubits, format='scipy.sparse')
    assert is_scalar_matrix(xy_matrix) == False

    cnot_matrix = represent(CNOT(1,0)*CNOT(1,0), nqubits=numqubits,
                          format='scipy.sparse')
    assert is_scalar_matrix(cnot_matrix) == True

    # Fails - wondering if it might be a floating point issue
    h_matrix = represent(H(0)*H(0), nqubits=numqubits, format='scipy.sparse')
    assert is_scalar_matrix(h_matrix) == True

def test_is_degenerate():
    gate_id = GateIdentity((X(0), Y(0), Z(0)))
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
    id_set = set([GateIdentity((X(0), X(0)))])
    assert bfs_identity_search(gate_list, 1, 2) == id_set

    # Set should not contain degenerate quantum circuits
    gate_list = [X(0), Y(0), Z(0)]
    id_set = set([GateIdentity((X(0), X(0))), 
                  GateIdentity((X(0), Y(0), Z(0))),
                  GateIdentity((Y(0), Y(0))),
                  GateIdentity((Z(0), Z(0)))])
    assert bfs_identity_search(gate_list, 1) == id_set

    id_set = set([GateIdentity((X(0), X(0))),
                  GateIdentity((X(0), Y(0), X(0), Y(0))),
                  GateIdentity((X(0), Y(0), Z(0))),
                  GateIdentity((X(0), Z(0), X(0), Z(0))),
                  GateIdentity((Y(0), Y(0))),
                  GateIdentity((Y(0), Z(0), Y(0), Z(0))),
                  GateIdentity((Z(0), Z(0)))])
    assert bfs_identity_search(gate_list, 1, 4) == id_set

