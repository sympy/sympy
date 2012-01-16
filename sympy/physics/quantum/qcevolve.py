"""Evolution of quantum circuits.

This module provides utilities to evolve quantum circuits in SymPy.
It uses Pyevolve, a Python library for developing genetic algorithms.
More information about Pyevolve available at:

http://pyevolve.sourceforge.net/index.html

"""

from pyevolve.GenomeBase import GenomeBase

__all__ = [
    'GQuantumCircuit'
]

class GQuantumCircuit(GenomeBase):
    """A representation of quantum circuits for genetic algorithms.

    Specifically, GQuantumCircuit is used for genetic programming,
    a specialization of genetic algorithms that evolve computer
    programs, or generally, algorithms.  Each quantum circuit
    represents a quantum algorithm.
    """
    pass
