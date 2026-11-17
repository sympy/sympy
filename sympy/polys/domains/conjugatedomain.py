"""Implementation of :class:`ConjugateDomain` class. """
from __future__ import annotations


from sympy.polys.domains.domain import Domain, Er
from sympy.utilities import public

@public
class ConjugateDomain(Domain[Er]):
    """Base class for domains supporting conjugation. """

    is_ConjugateDomain = True

    def conjugate(self, a: Er) -> Er:
        """Returns the complex conjugate of ``a``."""
        raise NotImplementedError
