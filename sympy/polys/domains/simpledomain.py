"""Implementation of :class:`SimpleDomain` class. """


from sympy.polys.domains.domain import Domain, Er
from sympy.utilities import public
from typing_extensions import Self

@public
class SimpleDomain(Domain[Er]):
    """Base class for simple domains, e.g. ZZ, QQ. """

    is_Simple: bool = True

    def inject(self, *gens) -> Self:
        """Inject generators into this domain. """
        return self.poly_ring(*gens)
