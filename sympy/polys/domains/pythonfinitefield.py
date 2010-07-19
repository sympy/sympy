"""Implementation of :class:`PythonFiniteField` class. """

from sympy.polys.domains.finitefield import FiniteField
from sympy.polys.domains.pythonintegerring import PythonIntegerRing

class PythonFiniteField(FiniteField):
    """Finite field based on Python's integers. """

    dom = PythonIntegerRing()
    alias = 'FF_python'

