"""Implementation of :class:`SimpleDomain` class. """


from sympy.polys.domains.domain import Domain, Er
from sympy.utilities import public

@public
class SimpleDomain(Domain[Er]):
    """Base class for simple domains, e.g. ZZ, QQ. """

    is_Simple = True

    def inject(self, *gens):
        """Inject generators into this domain. """
        return self.poly_ring(*gens)
