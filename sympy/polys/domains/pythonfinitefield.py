"""Implementation of :class:`PythonFiniteField` class. """
from __future__ import annotations


from sympy.polys.domains.finitefield import FiniteField
from sympy.polys.domains.pythonintegerring import PythonIntegerRing

from sympy.utilities import public

@public
class PythonFiniteField(FiniteField):
    """Finite field based on Python's integers. """

    alias = 'FF_python'

    def __init__(self, mod, symmetric=True):
        super().__init__(mod, PythonIntegerRing(), symmetric)
