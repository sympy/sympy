"""Implementation of :class:`PythonFiniteField` class. """

__all__ = ["PythonFiniteField"]

from sympy.polys.domains.finitefield import FiniteField
from sympy.polys.domains.pythonintegerring import PythonIntegerRing

class PythonFiniteField(FiniteField):
    """Finite field based on Python's integers. """

    alias = 'FF_python'

    def __init__(self, mod, symmetric=True):
        return super(PythonFiniteField, self).__init__(mod, PythonIntegerRing(), symmetric)
