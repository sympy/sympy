"""Implementation of :class:`CompositeDomain` class. """

from sympy.polys.domains.domain import Domain
from sympy.polys.polyerrors import GeneratorsError

class CompositeDomain(Domain):
    """Base class for composite domains, e.g. ZZ[x], ZZ(X). """

    is_Composite = True

    gens, ngens, symbols, domain = [None]*4

    def inject(self, *symbols):
        """Inject generators into this domain. """
        if not (set(self.symbols) & set(symbols)):
            return self.__class__.init(self.dom, *(self.symbols + symbols))
        else:
            raise GeneratorsError("common generators in %s and %s" % (self.symbols, symbols))
