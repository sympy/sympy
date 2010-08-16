"""Implementation of :class:`GMPYFiniteField` class. """

from sympy.polys.domains.finitefield import FiniteField
from sympy.polys.domains.gmpyintegerring import GMPYIntegerRing

class GMPYFiniteField(FiniteField):
    """Finite field based on Python's integers. """

    dom = GMPYIntegerRing()
    alias = 'FF_gmpy'

