"""Algebra between functions such as f+g or 2*f"""

from sympy.core import S
from sympy.sets import Set

__all__ = [
    "FunctionSet",
]

class FunctionSet(Set):
    """
    Set of functions who have same domains or codomains.

    Examples
    ========

    >>> from sympy import FunctionSet, Map, S

    >>> fs = FunctionSet(domain=S.Reals)
    >>> f = Map('f', domain=S.Integers)

    >>> f in fs
    True

    """

    def __new__(cls, domain=None, codomain=None, **kwargs):
        if domain is None:
            domain = S.UniversalSet
        if codomain is None:
            codomain = S.UniversalSet
        return super().__new__(cls, domain, codomain)

    @property
    def domain(self):
        return self.args[0]

    @property
    def codomain(self):
        return self.args[1]

    def _contains(self, other):
        return (
            self.domain.is_superset(other.domain)
            and self.codomain.is_superset(other.codomain)
        )

    def _element(self, name):
        from sympy.map import Map
        return Map(name, domain=self.domain, codomain=self.codomain)
