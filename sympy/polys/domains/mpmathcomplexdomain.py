"""Implementation of :class:`MPmathComplexDomain` class. """

from sympy.polys.domains.realdomain import RealDomain

class MPmathComplexDomain(RealDomain): # XXX: tmp solution
    """Complex domain. """

    alias = 'CC_mpmath'

    def __init__(self):
        pass

