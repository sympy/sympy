"""Implementation of :class:`GMPYFiniteField` class. """

from sympy.polys.domains.finitefield import FiniteField
from sympy.polys.domains.gmpyintegerring import GMPYIntegerRing

class GMPYFiniteField(FiniteField):
    """Finite field based on Python's integers. """

    alias = 'FF_gmpy'

    def __init__(self, mod, symmetric=True):
        return super(GMPYFiniteField, self).__init__(mod, GMPYIntegerRing(), symmetric)
