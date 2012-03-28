"""Evolution of quantum circuits.

This module provides utilities to evolve quantum circuits in SymPy.
It uses Pyevolve, a Python library for developing genetic algorithms.
More information about Pyevolve available at:

http://pyevolve.sourceforge.net/index.html

TODO or THINGS TO CONSIDER:
* Implement a mutator operator
* Implement an evaluator
* Implement a random replace (similar to reduce)
* Sorting the circuits to deal with 'trivial identities'
  i.e. X(0)X(1)X(0)X(1)
  Gate sorting in gate.py limited because can't swap gates that act on
  same qubit even if they do commute.
* When inserting circuits, use the dictionary of symbolic to real
  indices (from the circuit being manipulated) to create a circuit
  with real indices.
"""

from random import Random
from sympy import Basic
from pyevolve.GenomeBase import GenomeBase
from sympy.physics.quantum.circuitutils import (
         find_subcircuit_with_seq, replace_subcircuit_with_seq,
         conv2_symbolic_qubits_with_seq, conv2_real_qubits_with_seq)

__all__ = [
    'GQCBase',
    'GQCLinear',
    'random_reduce',
    'random_insert'
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

    def __new__(cls, *args):
        # args is the quantum circuit representing a genome
        obj = Basic.__new__(cls, *args)
        return obj

    @property
    def circuit(self):
        return self.args

    def __repr__(self):
        """Return a string representation of GQCBase"""
        rep = GenomeBase.__repr__(self)
        rep += "- GQCBase\n"
        rep += "\tCircuit:\t\t %s\n\n" % (self.circuit,)
        return rep

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

    def __new__(cls, *args, **kargs):
        # args is the quantum circuit representing a genome
        # kargs should a variable length dictionary
        obj = GQCBase.__new__(cls, *args)

        obj._genome_circuit = args
        obj._gate_identities = False
        obj._insert_choices = []

        # Doing this is a questionable design
        if 'GateIdentity' in kargs:
            obj._gate_identities = kargs['GateIdentity']

        if 'choices' in kargs:
            obj._insert_choices = kargs['choices']

        if obj._gate_identities:
            collapse_func = lambda acc, an_id: acc + list(an_id.symbolic_ids)
            ids_flat = reduce(collapse_func, obj._insert_choices, [])
            obj._insert_choices = ids_flat

        return obj

    @property
    def genome_circuit(self):
        return self._genome_circuit

    @genome_circuit.setter
    def genome_circuit(self, new_circuit):
        self._genome_circuit = new_circuit

    @property
    def insert_choices(self):
        # List of circuits that could be inserted into another circuit
        return self._insert_choices

    def __repr__(self):
        """Return a string representation of GQCBase"""
        rep = GQCBase.__repr__(self)
        rep += "- GQCLinear\n"
        rep += "\tIdentities:\t\t %s\n\n" % (self.identities,)
        return rep

def random_reduce(circuit, gate_ids, seed=None):
    """Shorten the length of a quantum circuit.

    random_reduce looks for circuit identities
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
    collapse_func = lambda acc, an_id: acc + list(an_id.symbolic_ids)
    ids = reduce(collapse_func, gate_ids, [])

    # List of identities found in circuit
    ids_found = []

    sym_circuit, mapping, ndx = conv2_symbolic_qubits_with_seq(*circuit)
    # Look for identities in circuit
    for an_id in ids:
        if find_subcircuit_with_seq(sym_circuit, an_id) > -1:
            id_real = conv2_real_qubits_with_seq(*an_id, qubit_map=mapping)
            ids_found.append(id_real)

    if len(ids_found) < 1:
        return circuit

    # Randomly choose an identity to remove
    remove_id = int_gen.randint(0, len(ids_found)-1)

    # Remove the identity
    new_circuit = replace_subcircuit_with_seq(circuit, ids_found[remove_id])

    return new_circuit

def random_insert(circuit, choices, symbolic=True, seed=None):
    """Insert a circuit into another quantum circuit.

    random_insert randomly selects a circuit from
    choices and randomly chooses a location to insert
    into circuit.

    Parameters
    ==========
    circuit : tuple, Gate
        A tuple of Gates representing a quantum circuit
    choices : list
        Set of circuit choices
    symbolic : bool
        Indicates whether the circuit choices use symbolic indices or
        not.  By default, assumes the choices are symbolic.
    seed : int
        Seed value for the random number generator
    """

    if len(choices) < 1:
        return circuit

    # Create the random integer generator with the seed
    int_gen = Random()
    int_gen.seed(seed)

    insert_loc = int_gen.randint(0, len(circuit))
    insert_circuit_loc = int_gen.randint(0, len(choices)-1)
    insert_circuit = choices[insert_circuit_loc]

    if symbolic:
        sym_circuit, mapping, ndx = conv2_symbolic_qubits_with_seq(*circuit)
        insert_circuit = conv2_real_qubits_with_seq(*insert_circuit,
                                 qubit_map=mapping)

    left = circuit[0:insert_loc]
    right = circuit[insert_loc:len(circuit)]

    return left + insert_circuit + right

"""
For the current problem of optimizing quantum circuits,
a semi-genetic programming approach is used where only
the mutation operation is applied to create future populations.

The first approach keeps all members of the population
within the same equivalence class of the original circuit.
"""
def linear_init(genome, **args):
    """Initialization function for optimizing quantum circuits
       problem using an linear circuit representataion.

    Parameters
    ==========
    genome : GQCLinear
        The genome in the population
    args : any (variable length)
        In many cases, will include an instance of GSimpleGA
    """

    new_circuit = random_insert(
                      genome.genome_circuit,
                      genome.insert_choices
                  )

    genome.genome_circuit = new_circuit

def linear_mutator(genome, **args):
    """Mutator function for optimizing quantum circuits
       problem using an linear circuit representataion.

    Parameters
    ==========
    genome : GQCLinear
        The genome in the population
    args : any (variable length)
        In many cases, will include an instance of GSimpleGA
    """

    # The mutator may either look for identities in
    # the circuit and make a reduction or it may
    # insert new identities into the circuit.

    # Return the number of mutations that occurred
    # with this genome (following convention)
    pass

def linear_eval(chromosome):
    """Evaluation function for optimizing quantum circuits
       problem using an linear circuit representataion."""
    pass

