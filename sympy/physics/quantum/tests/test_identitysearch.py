from sympy.physics.quantum.gate import (X, Y, Z, H, S, T, CNOT,
        IdentityGate, gate_simp)
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.identitysearch import (GateIdentity,
        is_scalar_matrix, bfs_identity_search)

def test_is_scalar_matrix():
    numqubits = 2
    id_gate = IdentityGate(1)
    id_matrix = represent(id_gate, nqubits=numqubits, format='scipy.sparse')
    assert is_scalar_matrix(id_matrix) == True

    id_matrix = represent(X(1)*Y(1), nqubits=numqubits, format='scipy.sparse')
    assert is_scalar_matrix(id_matrix) == False

    id_matrix = represent(CNOT(1,0)*CNOT(1,0), nqubits=numqubits,
                          format='scipy.sparse')
    assert is_scalar_matrix(id_matrix) == True

def test_bfs_identity_search():
    assert bfs_identity_search([], 1) == set()

    gate_list = [X(0)]
    id_set = set([GateIdentity((X(0), X(0)))])
    assert bfs_identity_search(gate_list, 1, 2) == id_set

    gate_list = gate_list + [Y(0), Z(0)]
    additional_set = set([GateIdentity((X(0), Y(0), Z(0))),
                          GateIdentity((X(0), Z(0), Y(0))),
                          GateIdentity((Y(0), X(0), Z(0))),
                          GateIdentity((Y(0), Y(0))),
                          GateIdentity((Y(0), Z(0), X(0))),
                          GateIdentity((Z(0), X(0), Y(0))),
                          GateIdentity((Z(0), Y(0), X(0))),
                          GateIdentity((Z(0), Z(0)))])
    id_set = id_set | additional_set
    assert bfs_identity_search(gate_list, 1) == id_set
