from collections import deque
from random import randint

from sympy import Mul, Basic
from sympy.matrices import Matrix, eye
from sympy.physics.quantum.gate import (X, Y, Z, H, S, T, CNOT,
        IdentityGate, gate_simp)
from sympy.physics.quantum.represent import represent

__all__ = [
    'generate_gate_rules'
    'GateIdentity',
    'is_scalar_matrix',
    'bfs_identity_search',
    'random_identity_search'
]

def generate_gate_rules_recursive(gate_seq, recurse_pt):

    seq = list(gate_seq)

    if (recurse_pt < 0):
        return []

    if (recurse_pt + 1 == len(gate_seq)):
        return [seq]

    permutations = []

    for i in range(recurse_pt, len(gate_seq)):
        print seq, recurse_pt
        current_permutes = generate_gate_rules_recursive(seq, recurse_pt + 1)
        permutations = permutations + current_permutes

        if (i + 1 < len(gate_seq)):            
            temp_gate = seq[i + 1]
            for j in reversed(range(recurse_pt, i + 1)):
                seq[j + 1] = seq[j]
            seq[recurse_pt] = temp_gate

    return permutations

def generate_gate_rules(gate_seq):
    # Not recursive version - unsure if interpreter
    # optimizes recursion 
    return None 

class GateIdentity(Basic):
    """Wrapper class for circuits that reduce to a scalar value."""

    def __new__(cls, circuit):
        # circuit should be a tuple
        obj = Basic.__new__(cls, circuit)
        obj._circuit = circuit
        obj._gate_rules = generate_gate_rules(circuit)

        return obj

    @property
    def circuit(self):
        return self._circuit

    @property
    def gate_rules(self):
        return self._gate_rules

    def __str__(self):
        """Returns the string of gates in a tuple."""
        return str(self.circuit)

# Dynamically generate and concatenate a list of all possible
# sympy gate objects in given space.
def construct_gate_list(numqubits):
    Xs = [X(i) for i in xrange(numqubits)]
    Ys = [Y(i) for i in xrange(numqubits)]
    Zs = [Z(i) for i in xrange(numqubits)]
    Hs = [H(i) for i in xrange(numqubits)]
    Ss = [S(i) for i in xrange(numqubits)]
    Ts = [T(i) for i in xrange(numqubits)]
    CNOTs = [CNOT(i,j) for i in xrange(numqubits)
                           for j in xrange(numqubits) if i != j]
   
    gate_list = Xs+Ys+Zs+Hs+Ss+Ts+CNOTs

    return gate_list

# Dynamically generate and concatenate a list of all possible
# scipy.sparse gate matrices in given space.
def construct_matrix_list(numqubits):
    xs = [represent(X(i), nqubits=numqubits, format='scipy.sparse')
          for i in xrange(numqubits)]
    ys = [represent(Y(i), nqubits=numqubits, format='scipy.sparse')
          for i in xrange(numqubits)]
    zs = [represent(Z(i), nqubits=numqubits, format='scipy.sparse')
          for i in xrange(numqubits)]
    hs = [represent(H(i), nqubits=numqubits, format='scipy.sparse')
          for i in xrange(numqubits)]
    ss = [represent(S(i), nqubits=numqubits, format='scipy.sparse')
          for i in xrange(numqubits)]
    ts = [represent(T(i), nqubits=numqubits, format='scipy.sparse')
          for i in xrange(numqubits)]
    cnots = [represent(CNOT(i,j), nqubits=numqubits, format='scipy.sparse')
             for i in xrange(numqubits) for j in xrange(numqubits) if i != j]

    matrix_list = xs+ys+zs+hs+ss+ts+cnots

    return matrix_list

def is_scalar_matrix(matrix):
    """Checks if given scipy.sparse matrix is a scalar matrix."""

    if (list(matrix.nonzero()[0]) == list(matrix.nonzero()[1])):
        diag = list(matrix.diagonal())
        if (diag.count(diag[0]) == len(diag)):
            return True
    return False

def bfs_identity_search(gate_list, numqubits, max_depth=0):
    # Breadth first search might be more efficient because it eliminates
    # a search down paths like ZZZZZZ or XXXXYY.
    # Returns the set of gate identities from the list.

    # If max depth of a path isn't given, use the length of the gate_list
    if (max_depth == 0):
        max_depth = len(gate_list)

    # Root of BFS tree is an IdentityGate(0)
    id_gate = IdentityGate(0)
    queue = deque([(id_gate,)])

    # Create an empty set of gate identities
    ids = set()

    # Begin searching for gate identities in given space.
    while (len(queue) > 0):
        current_circuit = queue.popleft()

        for next_gate in gate_list:
            new_circuit = current_circuit + (next_gate,)
            matrix_version = represent(Mul(*new_circuit), nqubits=numqubits,
                                       format='scipy.sparse')

            # In many cases when the matrix is a scalar value,
            # the evaluated matrix will actually be an integer          
            if (isinstance(matrix_version, int)):
                # When adding a gate identity, remove the
                # identity gate at the beginning of the tuple
                ids.add(GateIdentity(new_circuit[1:]))

            # If a matrix is equivalent to a scalar value is found
            elif (is_scalar_matrix(matrix_version)):
                # When adding a gate identity, remove the
                # identity gate at the beginning of the tuple
                ids.add(GateIdentity(new_circuit[1:]))

            # Number of operators in the circuit gives the
            # number of gates in the circuit
            elif (len(new_circuit) < max_depth + 1):
                queue.append(new_circuit)

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
