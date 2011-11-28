from collections import deque
from random import randint

from sympy import Mul, Basic, Number
from sympy.matrices import Matrix, eye
from sympy.physics.quantum.gate import (X, Y, Z, H, S, T, CNOT,
        IdentityGate, gate_simp)
from sympy.physics.quantum.represent import represent
from sympy.physics.quantum.operator import (UnitaryOperator,
        HermitianOperator)
from sympy.physics.quantum.dagger import Dagger

__all__ = [
    'generate_gate_rules',
    'GateIdentity',
    'is_scalar_matrix',
    'is_degenerate',
    'is_reducible',
    'bfs_identity_search',
    'random_identity_search'
]

def generate_gate_rules(*gate_seq):
    """Returns a list of equivalent gate identities.

    The function uses the four operations (LL, LR, RL, RR)
    to find equivalent gate identities.

    Note that all equivalent identities are reachable in n operations
    from the starting gate identity, where n is the number of gates
    in the sequence.

    Parameters
    ----------

    Examples
    --------

    """

    # Each item in queue is a 3-tuple:
    #     i)   first item is the left side of an equality
    #    ii)   second item is the right side of an equality
    #   iii)   third item is the number of operations performed
    # The argument, gate_seq, will start on the left side, and
    # the right side will be empty, implying the presence of an
    # identity.
    queue = deque()
    # visited is a list of equalities that's been visited
    vis = []
    # A list of equivalent gate identities
    gate_rules = []
    # Maximum number of operations to perform
    max_ops = len(gate_seq)

    queue.append((gate_seq, (), 0))
    vis.append((gate_seq, ()))
    gate_rules.append(gate_seq)

    while (len(queue) > 0):
        rule = queue.popleft()

        left = rule[0]
        rite = rule[1]
        ops = rule[2]

        if (len(left) > 0):
            ll_op = left[0]
            lr_op = left[len(left)-1]

        if (len(rite) > 0):
            rl_op = rite[0]
            rr_op = rite[len(rite)-1]

        # Do a LL, if possible
        if (len(left) > 0 and Dagger(ll_op)*ll_op == 1):
            # Get the new left side w/o the leftmost gate
            new_left = left[1:len(left)]
            # Add the leftmost gate to the left position on the right side
            new_rite = (Dagger(ll_op),) + rite

            new_rule = (new_left, new_rite)
            # If the left side is empty (left side is scalar)
            if (len(new_left) == 0 and new_rite not in gate_rules):
                gate_rules.append(new_rite)                
            # If the equality has not been seen and has not reached the
            # max limit on operations
            elif (new_rule not in vis and ops + 1 < max_ops):
                queue.append(new_rule + (ops + 1,))

            vis.append(new_rule)

        # Do a LR, if possible
        if (len(left) > 0 and Dagger(lr_op)*lr_op == 1):
            # Get the new left side w/o the rightmost gate
            new_left = left[0:len(left)-1]
            # Add the rightmost gate to the right position on the right side
            new_rite = rite + (Dagger(lr_op),)

            new_rule = (new_left, new_rite)
            if (len(new_left) == 0 and new_rite not in gate_rules):
                gate_rules.append(new_rite)
            elif (new_rule not in vis and ops + 1 < max_ops):
                queue.append(new_rule + (ops + 1,))

            vis.append(new_rule)

        # Do a RL, if possible
        if (len(rite) > 0 and Dagger(rl_op)*rl_op == 1):
            # Get the new right side w/o the leftmost gate
            new_rite = rite[1:len(rite)]
            # Add the leftmost gate to the left position on the left side
            new_left = (Dagger(rl_op),) + left

            new_rule = (new_left, new_rite)
            if (len(new_rite) == 0 and new_left not in gate_rules):
                gate_rules.append(new_left)
            elif (new_rule not in vis and ops + 1 < max_ops):
                queue.append(new_rule + (ops + 1,))

            vis.append(new_rule)

        # Do a RR, if possible
        if (len(rite) > 0 and Dagger(rr_op)*rr_op == 1):
            # Get the new right side w/o the rightmost gate
            new_rite = rite[0:len(rite)-1]
            # Add the rightmost gate to the right position on the right side
            new_left = left + (Dagger(rr_op),)

            new_rule = (new_left, new_rite)
            if (len(new_rite) == 0 and new_left not in gate_rules):
                gate_rules.append(new_left)
            elif (new_rule not in vis and ops + 1 < max_ops):
                queue.append(new_rule + (ops + 1,))

            vis.append(new_rule)

    return gate_rules

class GateIdentity(Basic):
    """Wrapper class for circuits that reduce to a scalar value.

    Parameters
    ----------

    Examples
    --------
    """
    def __new__(cls, *args):
        # args should be a tuple - a variable length argument list
        obj = Basic.__new__(cls, *args)
        obj._gate_rules = generate_gate_rules(*args)

        return obj

    @property
    def circuit(self):
        return self.args

    @property
    def gate_rules(self):
        return self._gate_rules

    def __str__(self):
        """Returns the string of gates in a tuple."""
        return str(self.circuit)

def construct_gate_list(numqubits):
    """Dynamically generate and concatenate a list of all possible
       sympy gate objects in given space."""

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

def construct_matrix_list(numqubits):
    """ Dynamically generate and concatenate a list of all possible
        scipy.sparse gate matrices in given space."""

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

def is_scalar_sparse_matrix(matrix):
    """Checks if a given scipy.sparse matrix is a scalar matrix."""

    if (list(matrix.nonzero()[0]) == list(matrix.nonzero()[1])):
        diag = list(matrix.diagonal())
        if (diag.count(diag[0]) == len(diag)):
            return True
    return False

def is_scalar_matrix(circuit, numqubits, identity_only):
    """Checks if a given circuit, in matrix form, is equivalent to
       a scalar value."""

    # A sparse matrix is faster but there's a few problems with it,
    # such as not being able to determine H(0)*H(0) is the identity matrix.

    matrix_version = represent(Mul(*circuit), nqubits=numqubits)

    # In some cases, represent returns a 1D scalar value in place
    # of a multi-dimensional scalar matrix
    if (isinstance(matrix_version, Number)):
        return matrix_version == 1 if identity_only else True

    # If represent returns a matrix, check if the matrix is diagonal
    # and if every item along the diagonal is the same
    else:
        # Added up the diagonal elements
        matrix_trace = matrix_version.trace()
        # Divide the trace by the first element in the matrix
        # if matrix is not required to be the identity matrix
        adjusted_matrix_trace = (matrix_trace/matrix_version[0] 
                                 if not identity_only
                                 else matrix_trace)
        # The matrix is scalar if it's diagonal and the adjusted trace
        # value is equal to 2^numqubits
        return (matrix_version.is_diagonal() and
                adjusted_matrix_trace == pow(2, numqubits))

def is_degenerate(identity_set, gate_identity):
    """Checks if a gate identity is a permutation of another identity.

    Parameters
    ----------

    Examples
    --------
    """

    # For now, just iteratively go through the set and check if the current
    # gate_identity is a permutation of an identity in the set
    for an_id in identity_set:
        if (gate_identity in an_id.gate_rules):
            return True
    return False

def is_reducible(circuit, numqubits, begin, end):
    """Determines if a subcircuit in circuit is reducible to a scalar value.

    Parameters
    ----------

    Examples
    --------
    """

    current_circuit = ()
    # Start from the gate at "end" and go down to almost the gate at "begin"
    for ndx in reversed(range(begin, end)):
        next_gate = circuit[ndx]
        current_circuit = (next_gate,) + current_circuit

        # If a circuit as a matrix is equivalent to a scalar value
        if (is_scalar_matrix(current_circuit, numqubits, False)):
            return True

    return False

def bfs_identity_search(gate_list, numqubits, **kwargs):
    """Constructs a set of gate identities from the list of possible gates.

    Performs a breadth first search over the space of gate identities.
    This allows the finding of the shortest gate identities first.

    Parameters
    ----------

    Examples
    --------
    """

    # If max depth of a path isn't given, use the length of the gate_list
    if ("max_depth" in kwargs and kwargs["max_depth"] > 0):
        max_depth = kwargs["max_depth"]
    else:        
        max_depth = len(gate_list)

    identity_only = (kwargs["identity_only"] if "identity_only" in kwargs else
                     False)

    # Start with an empty sequence (implicitly contains an IdentityGate)
    queue = deque([()])

    # Create an empty set of gate identities
    ids = set()

    # Begin searching for gate identities in given space.
    while (len(queue) > 0):
        current_circuit = queue.popleft()

        for next_gate in gate_list:
            new_circuit = current_circuit + (next_gate,)
            #matrix_version = represent(Mul(*new_circuit), nqubits=numqubits,
            #                           format='scipy.sparse')

            # Determines if a (strict) subcircuit is a scalar matrix
            circuit_reducible = is_reducible(new_circuit, numqubits, 
                                             1, len(new_circuit))

            # In many cases when the matrix is a scalar value,
            # the evaluated matrix will actually be an integer          
            if (is_scalar_matrix(new_circuit, numqubits, identity_only) and
                not is_degenerate(ids, new_circuit) and
                not circuit_reducible):
                ids.add(GateIdentity(*new_circuit))

            elif (len(new_circuit) < max_depth and
                not circuit_reducible):
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

    return (circuit if is_scalar_matrix(matrix_version, numqubits, False)
            else None)
