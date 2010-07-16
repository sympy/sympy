"""Implementation of :class:`CompositeDomain` class. """

from sympy.polys.domains.domain import Domain

class CompositeDomain(Domain):
    """Base class for composite domains, e.g. ZZ[x]. """

    is_Composite = True

