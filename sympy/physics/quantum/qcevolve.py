"""Evolution of quantum circuits.

This module provides utilities to evolve quantum circuits in SymPy.
It uses Pyevolve, a Python library for developing genetic algorithms.
More information about Pyevolve available at:

http://pyevolve.sourceforge.net/index.html

"""

from sympy import Basic
from pyevolve.GenomeBase import GenomeBase

__all__ = [
    'GQCBase',
    'kmp_table',
    'find_subcircuit'
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

def find_subcircuit(circuit, subcircuit):
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

    """

    if len(subcircuit) == 0 or len(subcircuit) > len(circuit):
        return -1

    # Location in circuit
    pos = 0
    # Location in the subcircuit
    index = 0
    # 'Partial match' table
    table = kmp_table(subcircuit)

    while (pos + index) < len(circuit):
        if subcircuit[index] == circuit[pos + index]:
            index = index + 1
        else:
            pos = pos + index - table[index]
            index = table[index] if table[index] > -1 else 0

        if index == len(subcircuit):
            return pos

    return -1

def qc_reduce(circuit, ids):
    """Shorten the length of a quantum circuit.

    qc_reduce looks for subcircuits - circuit identities -
    in circuit, removes them, and returns the new shorter circuit.

    Parameters
    ==========
    circuit : tuple, Gate
        A tuple of Gates representing a quantum circuit
    ids : list, GateIdentity
        List of gate identities to find in circuit
    """

    # Try greedy approach - find longest subcircuits first
    # before optimizing.  Theoretically, it's possible that
    # circuit may contain two shorter identities and one
    # longer identity and removing the two shorter identities
    # may be more optimal than removing the one longer identity.

    # Sort places the shortest length circuits at the top of the list
    ids.sort().reverse()

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

