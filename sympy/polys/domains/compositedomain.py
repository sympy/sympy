"""Implementation of :class:`CompositeDomain` class. """

from sympy.polys.domains.domain import Domain
from sympy.polys.polyerrors import GeneratorsError

class CompositeDomain(Domain):
    """Base class for composite domains, e.g. ZZ[x], ZZ(X). """

    is_Composite = True

    def inject(self, *gens):
        """Inject generators into this domain. """
        if not (set(self.gens) & set(gens)):
            return self.__class__(self.dom, *(self.gens + gens))
        else:
            raise GeneratorsError("common generators in %s and %s" % (self.gens, gens))
