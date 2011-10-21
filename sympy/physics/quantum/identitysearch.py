from sets import Set

from sympy import Mul
from sympy.matrices import Matrix, eye
from sympy.physics.quantum.gate import X, Y, Z, H, S, T, CNOT, IdentityGate, gate_simp
from sympy.physics.quantum.represent import represent

# Possible inherit Python List
class GateIdentity():
    
# Save the identities to a set
class GateIdentitySet(Set):
    def write_to_file(filename, append=False):
        '''Write gate identities to a specified file.'''

        # File is truncated unless given a flag to not truncate.
        identities = open(filename, 'a') if append else open(filename, 'w+')
        identities.write(identity+'\n')
        identities.close()

# Number of qubits in the spaces (to be operated on).
def default_numqubits():
    return 2

# Number base the counter will use or number of gates possible in given space.
def base(gate_list):
    return len(gate_list)

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

# Initial number list - start identity search with 2 gates
def initial_num():
    return [0, 1]

def count(base, number):
    # Uses a list of numbers (as digits) to count in a specified base.
    for i in xrange(len(number)):
        if number[i] == (base-1):
            number[i] = 0
        else:
            number[i] += 1
            break
    if number[i] == 0:
        number.append(1)

def number_to_gates(number, gate_list):
    # Converts a list of numbers into a list of corresponding gate objects.
    gates = []
    for digit in number:
        gates.append(gate_list[digit])
    return gates

def number_to_matrices(number, gate_list):
    # Converts a list of numbers into a list of corresponding matrices.
    matrices = []
    if isinstance(number, int):
        return matrix_list[number]
    for digit in number:
        matrices.append(matrix_list[digit])
    return matrices

def matrix_mul(matrices):
    # Multiplies given scipy.sparse matrices.
    mul = represent(IdentityGate(0), nqubits=numqubits, format='scipy.sparse')
    for matrix in matrices:
        mul = mul*matrix
    return mul

def is_scalar_matrix(matrix):
    # Checks if given scipy.sparse matrix is a scalar matrix.
    if list(matrix.nonzero()[0]) == list(matrix.nonzero()[1]):
        diag = list(matrix.diagonal())
        if diag.count(diag[0]) == len(diag):
            return True
    return False

def iterative_identity_search(numqubits):
    # Begin searching for gate identities in given space.

    # Create the list of gates and matrices
    gate_list = construct_gate_list(numqubits)
    matrix_list = construct_matrix_list(numqubits)

    num = initial_num()

    # Iteratively search for identities by increasing the
    # number of gates in the list
    while True:
        if num[0] is 0:
            matrices = number_to_matrices(num[1:], matrix_list)
            cached = matrix_mul(matrices)

        gates = number_to_gates(num, gate_list)

        if len(gate_simp(Mul(*gates)).args) == len(gates):
            circuit = number_to_matrices(num[0], matrix_list)*cached
            if is_scalar_matrix(circuit) and check_identity(filename, str(gate_simp(Mul(*gates)).args)):
                write_identity(filename, str(gate_simp(Mul(*gates)).args))
                print num
                print gate_simp(Mul(*gates)).args
                print circuit
                print ''
        count(base, num)

def random_identity_search(numqubits):

