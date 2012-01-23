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

def find_subcircuit(circuit, subcircuit):
    """Finds the subcircuit in circuit, if it exists.

    If the subcircuit exists, the index of the start of
    the subcircuit in circuit is returned; otherwise,
    -1 is returned.

    Parameters
    ==========
    circuit : tuple, Gate
        A tuple of Gates representing a quantum circuit
    subcircuit : tuple, Gate
        A tuple of Gates to check if circuit contains

    """    
    if len(subcircuit) == 0 or len(subcircuit) > len(circuit):
        return -1

    # Location of subcircuit that this identity contains
    index = 0
    loc = -1
    for gate in circuit:
        if index < len(subcircuit):
            index = index + 1 if gate == subcircuit[index] else 0
            # if index == 1, then the start of the subcircuit is found
            loc = circuit.index(gate) if index == 1 else loc

        # If the value of index reaches the length of subcircuit
        # then this identity contains subcircuit
        if index == len(subcircuit):
            return loc

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

