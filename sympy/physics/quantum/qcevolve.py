"""Evolution of quantum circuits.

This module provides utilities to evolve quantum circuits in SymPy.
It uses Pyevolve, a Python library for developing genetic algorithms.
More information about Pyevolve available at:

http://pyevolve.sourceforge.net/index.html

"""

from random import Random
from sympy import Basic
from pyevolve.GenomeBase import GenomeBase

__all__ = [
    'GQCBase',
    'kmp_table',
    'find_subcircuit',
    'qc_remove_subcircuit',
    'qc_random_reduce'
]

class GQCBase(Basic, GenomeBase):
    """A base representation of quantum circuits for genetic algorithms.

    Specifically, GQCBase is used for genetic programming,
    a specialization of genetic algorithms that evolve computer
    programs, or generally, algorithms.  Each quantum circuit
    represents a quantum algorithm.

    The base class provides no default behavior other than those
    provided by Basic and GenomeBase; all subclasses are expected
    to provide class specific implementations.
    """
    pass

class GQCLinear(GQCBase):
    """A linear program representation of quantum circuits.

    This class also does not provide behavior for the following:

        * Initializer operator
        * Mutator operator
        * Crossover operator
        * Evaluation function

    These functions must be set before an instance of the
    class can be used.

    GQCLinear was created to provide a meaningful name
    for a representation of quantum circuits.
    """
    pass

def kmp_table(word):
    """Build the 'partial match' table of the
       Knuth-Morris-Pratt algorithm.

    Note: This is applicable to strings or quantum circuits.
    """

    # Current position in subcircuit
    pos = 2
    # Beginning position of candidate substring that
    # may reappear later in word
    cnd = 0
    # The 'partial match' table that helps one determine
    # the next location to start substring search
    table = list()
    table.append(-1)
    table.append(0)

    while pos < len(word):
        if word[pos-1] == word[cnd]:
            cnd = cnd + 1
            table.append(cnd)
            pos = pos + 1
        elif cnd > 0:
            cnd = table[cnd]
        else:
            table.append(0)
            pos = pos + 1

    return table

def find_subcircuit(circuit, subcircuit, start=0, end=0):
    """Finds the subcircuit in circuit, if it exists.

    If the subcircuit exists, the index of the start of
    the subcircuit in circuit is returned; otherwise,
    -1 is returned.  The algorithm that is implemented
    is the Knuth-Morris-Pratt algorithm.

    Parameters
    ==========
    circuit : tuple, Gate
        A tuple of Gates representing a quantum circuit
    subcircuit : tuple, Gate
        A tuple of Gates to check if circuit contains
    start : int
        The location to start looking for subcircuit.
        If start is the same or past end, -1 is returned.
    end : int
        The last place to look for a subcircuit.  If end
        is less than 1 (one), then the length of circuit
        is taken to be end.
    """

    if len(subcircuit) == 0 or len(subcircuit) > len(circuit):
        return -1

    if end < 1:
        end = len(circuit)

    # Location in circuit
    pos = start
    # Location in the subcircuit
    index = 0
    # 'Partial match' table
    table = kmp_table(subcircuit)

    while (pos + index) < end:
        if subcircuit[index] == circuit[pos + index]:
            index = index + 1
        else:
            pos = pos + index - table[index]
            index = table[index] if table[index] > -1 else 0

        if index == len(subcircuit):
            return pos

    return -1

def qc_remove_subcircuit(circuit, subcircuit, pos=0):
    """Removes subcircuit from circuit, if it exists.

    If multiple instances of subcircuit exists, the
    first instance is removed.  A location to check may
    be optionally given.  If subcircuit can't be found,
    circuit is returned.

    Parameters
    ==========
    circuit : tuple, Gate
        A quantum circuit represented by a tuple of Gates
    subcircuit : tuple, Gate
        A quantum circuit to remove from circuit
    pos : int
        The location to remove subcircuit, if it exists.
        This may be used if it is known beforehand that
        multiple instances exist, and it is desirable
        to remove a specific instance.  If a negative number
        is given, pos will be defaulted to 0.
    """

    if pos < 0:
        pos = 0

    # Look for the subcircuit starting at pos
    loc = find_subcircuit(circuit, subcircuit, start=pos)

    # If subcircuit was found
    if loc > -1:
        # Get the gates to the left of subcircuit
        left = circuit[0:loc]
        # Get the gates to the right of subcircuit
        right = circuit[loc + len(subcircuit):len(circuit)]
        # Recombine the left and right side gates into a circuit
        circuit = left + right

    return circuit

def qc_random_reduce(circuit, gate_ids, seed=None):
    """Shorten the length of a quantum circuit.

    qc_random_reduce looks for circuit identities
    in circuit, randomly chooses one to remove,
    and returns a shorter yet equivalent circuit.
    If no identities are found, the same circuit
    is returned.

    Parameters
    ==========
    circuit : tuple, Gate
        A tuple of Gates representing a quantum circuit
    gate_ids : list, GateIdentity
        List of gate identities to find in circuit
    seed : int
        Seed value for the random number generator
    """

    if len(gate_ids) < 1:
        return circuit

    # Create the random integer generator with the seed
    int_gen = Random()
    int_gen.seed(seed)

    # Flatten the GateIdentity objects (with gate rules)
    # into one single list
    collapse_func = lambda acc, an_id: acc + an_id.gate_rules
    ids = reduce(collapse_func, gate_ids, [])

    # List of identities found in circuit
    ids_found = []

    # Look for identities in circuit
    for an_id in ids:
        if find_subcircuit(circuit, an_id) > -1:
            ids_found.append(an_id)

    # Randomly choose an identity to remove
    remove_id = int_gen.randint(0, len(ids_found)-1)

    # Remove the identity
    new_circuit = qc_remove_subcircuit(circuit, ids_found[remove_id])

    return new_circuit
    
# For now leaving this code in as a reminder of past progress
# Since using genetic programming, mutation by reduction
# should be a simple procedure.
# ===========================================================
def conflict_indices(pair, pairs):
    """Returns a list of positions in pairs where
    the endpoints in pair intersect with another set
    of endpoints in pairs.
    """
    pass

def qc_reduce(circuit, ids, quant=0, homogeneous=True):
    """Shorten the length of a quantum circuit.

    qc_reduce looks for subcircuits - circuit identities -
    in circuit, removes them, and returns the new shorter circuit.

    Parameters
    ==========
    circuit : tuple, Gate
        A tuple of Gates representing a quantum circuit
    ids : list, GateIdentity
        List of gate identities to find in circuit
    quant : int
        Number of identities to remove.  The default
        is to remove the maximum number of identities
        that can be found.
    homogeneous : bool
        True to remove up to quant number of the
        same identities in circuit;
        False to remove up to the first quant number of
        identities found

    """

    # Greedy approach - look for all locations of a
    # subcircuit in circuit and remove the subcircuits
    # that produce the shortest resulting circuit.

    # The circuits to remove from circuit, and the list
    # may contain more than one copy of the same identity.
    # Each item in the list is a 2-tuple:
    #     i) the identity to remove
    #    ii) the endpoints the identity occupies
    # The endpoints are a 2-tuple giving the location of the
    # first gate and the location after the last gate
    remove_circuits = []
    # Index into the list of gate identities
    i = 0
    # Location to start finding identities in circuit
    pos = 0

    # Overview of following code
    # collapse GateIdentity objects (with gate rules) into one list
    # for each identity
    #     while pos less than end of circuit
    #         get location of identity in circuit
    #         if location doesn't conflict with another identity
    #             save endpoints of identity in circuit
    #             save the circuit
    #         else
    #             if 'new' identity is longer than 'old' identities
    #                 remove 'old' identities
    #                 save 'new' identity
                
    #while i < len(ids):
    #    identity = ids[i]
    #    if find_subcircuit
    pass
# ===========================================================

def qc_random_insert(circuit, choices, identity=False):
    """Insert a circuit into another quantum circuit.

    qc_insert looks for subcircuits - circuit identities -
    in circuit, makes the reduction, and returns the new
    smaller circuit.

    Parameters
    ==========
    circuit : tuple, Gate
        A tuple of Gates representing a quantum circuit
    choices : list
        Set of circuit choices
    identity : boolean
        Indicates whether choices are of type GateIdentity
        or tuples
    """
    pass

"""
For the current problem of optimizing quantum circuits,
a semi-genetic programming approach is used where only
the mutation operation is applied to create future populations.

The first approach keeps all members of the population
within the same equivalence class of the original circuit.
"""
def qcopt_linear_init(orig_circuit, ids, popsize=100):
    """Initialization function for optimizing quantum circuits
       problem using an linear circuit representataion."""
    pass

def qcopt_linear_mutator():
    """Mutator function for optimizing quantum circuits
       problem using an linear circuit representataion."""

    # The mutator may either look for identities in
    # the circuit and make a reduction or it may
    # insert new identities into the circuit.
    pass

def qcopt_linear_eval():
    """Evaluation function for optimizing quantum circuits
       problem using an linear circuit representataion."""
    pass

