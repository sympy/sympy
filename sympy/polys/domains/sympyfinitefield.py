"""Implementation of :class:`SymPyFiniteField` class. """

from sympy.polys.domains.finitefield import FiniteField
from sympy.polys.domains.sympyintegerring import SymPyIntegerRing

class SymPyFiniteField(FiniteField):
    """Finite field based on SymPy's integers. """

    dom = SymPyIntegerRing()
    alias = 'FF_sympy'

