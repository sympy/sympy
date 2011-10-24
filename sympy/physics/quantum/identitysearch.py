from collections import deque
from random import randint

from sympy import Mul
from sympy.matrices import Matrix, eye
from sympy.physics.quantum.gate import X, Y, Z, H, S, T, CNOT, IdentityGate, gate_simp
from sympy.physics.quantum.represent import represent

class GateIdentity(Basic):
    """Wrapper class for circuits that reduce to a scalar value."""

    def __new__(cls, circuit):
        obj = Basic.__new__(cls, circuit)
        obj._circuit = circuit
        return obj

    @property
    def circuit(self):
        return self._circuit

    def generate_gate_rules():
        raise NotImplementedError(
            "Generation of gate rules from identity not implemented."
        )

    def __str__(self):
        return str(self.circuit)

# Dynamically generate and concatenate a list of all possible sympy gate objects in given space.
def construct_gate_list(numqubits):
    Xs = [X(i) for i in xrange(numqubits)]
    Ys = [Y(i) for i in xrange(numqubits)]
    Zs = [Z(i) for i in xrange(numqubits)]
    Hs = [H(i) for i in xrange(numqubits)]
    Ss = [S(i) for i in xrange(numqubits)]
    Ts = [T(i) for i in xrange(numqubits)]
    CNOTs = [CNOT(i,j) for i in xrange(numqubits) for j in xrange(numqubits) if i != j]
    
    gate_list = Xs+Ys+Zs+Hs+Ss+Ts+CNOTs

    return gate_list

# Dynamically generate and concatenate a list of all possible scipy.sparse gate matrices in given space.
def construct matrix_list(numqubits):
    xs = [represent(X(i), nqubits=numqubits, format='scipy.sparse') for i in xrange(numqubits)]
    ys = [represent(Y(i), nqubits=numqubits, format='scipy.sparse') for i in xrange(numqubits)]
    zs = [represent(Z(i), nqubits=numqubits, format='scipy.sparse') for i in xrange(numqubits)]
    hs = [represent(H(i), nqubits=numqubits, format='scipy.sparse') for i in xrange(numqubits)]
    ss = [represent(S(i), nqubits=numqubits, format='scipy.sparse') for i in xrange(numqubits)]
    ts = [represent(T(i), nqubits=numqubits, format='scipy.sparse') for i in xrange(numqubits)]
    cnots = [represent(CNOT(i,j), nqubits=numqubits, format='scipy.sparse') for i in xrange(numqubits) for j in xrange(numqubits) if i != j]
    
    matrix_list = xs+ys+zs+hs+ss+ts+cnots

    return matrix_list

def is_scalar_matrix(matrix):
    """Checks if given scipy.sparse matrix is a scalar matrix."""

    if list(matrix.nonzero()[0]) == list(matrix.nonzero()[1]):
        diag = list(matrix.diagonal())
        if diag.count(diag[0]) == len(diag):
            return True
    return False

def write_to_file(identityset, filename, append=False):
    """Write gate identities to a specified file."""

    # File is truncated unless given a flag to not truncate.
    identities = open(filename, 'a') if append else open(filename, 'w+')

    for identity in identityset:
        identities.write(str(identity) + "\n")

    identities.close()

def bfs_identity_search(gate_list, numqubits):
    # Breadth first search might be more efficient because it eliminates
    # a search down paths like ZZZZZZ or XXXXYY.
    # Returns the set of gate identities from the list.

    # For now, limit size of identity search based on size of list
    max_length = len(gate_list)

    # Root of BFS tree is an IdentityGate(0)
    id_gate = IdentityGate(0)
    queue = deque([id_gate])

    # Create an empty set of gate identities
    ids = set()

    # Begin searching for gate identities in given space.
    while (queue.__len__() > 0):
        current_circuit = queue.popleft()

        for next_gate in gate_list:
            new_circuit = current_circuit*next_gate
            matrix_version = represent(new_circuit, nqubits=numqubits,
                                       format='scipy.sparse')

            # If a matrix equivalent to a scalar value is found
            if (is_scalar_matrix(matrix_version)):
                ids.add(GateIdentity(new_circuit))

            # Number of operators in the circuit gives the
            # number of gates in the circuit
            elif (new_circuit.count_op < max_length):
                queue.append(circuit)

    return ids

def random_identity_search(gate_list, numgates, numqubits):
    """Randomly selects numgates from gate_list and checks if it is
       a gate identity.

       If the circuit is a gate identity, the circuit is returned;
       Otherwise, None is returned
    """

    gate_size = len(gate_list)
    circuit = IdentityGate(0)

    for i in range(numgates):
        next_gate = gate_list[randint(0, gate_size - 1)]
        circuit = initial_circuit*next_gate

    matrix_version = represent(circuit, nqubits=numqubits,
                               format='scipy.sparse')

    return circuit if is_scalar_matrix(matrix_version) else None
