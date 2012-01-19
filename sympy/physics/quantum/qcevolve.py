"""Evolution of quantum circuits.

This module provides utilities to evolve quantum circuits in SymPy.
It uses Pyevolve, a Python library for developing genetic algorithms.
More information about Pyevolve available at:

http://pyevolve.sourceforge.net/index.html

"""

from sympy import Basic
from pyevolve.GenomeBase import GenomeBase

__all__ = [
    'GQCBase'
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

# At the moment, this function would construct circuits of popsize
# by randomly choosing gates.  The number of gates in a circuit range
# from minsize to maxsize.
def qclinear_initializator(gates, nqubits, popsize=100,
        minsize=1, maxsize=1):
    pass

class GQCLinear(GQCBase):
    """A linear program representation of quantum circuits.

    TODO:
    * Initializer operator
    * Mutator operator
    * Crossover operator
    * Evaluator function
    """
    pass
