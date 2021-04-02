"""

Module for the DomainScalar class.

A DomainScalar represents an element which is in a particular
Domain. The idea is that the DomainScalar class provides the
convenience routines for unifying elements with different domains.

It assists in Scalar Multiplication and getitem for DomainMatrix.

"""
from ..constructor import construct_domain


class DomainScalar:
    r"""
    docstring
    """

    def __init__(self, element, domain):
        self.element = element
        self.domain = domain

    @classmethod
    def new(cls, element, domain):
        return cls(element, domain)

    def __repr__(self):
        return repr(self.element)

    @classmethod
    def from_sympy(cls, expr):
        [domain, [element]] = construct_domain([expr])
        return cls.new(element, domain)

    def to_domain(self, domain):
        element = domain.convert_from(self.element, self.domain)
        return self.new(element, domain)

    def unify(self, other):
        domain = self.domain.unify(other.domain)
        return self.to_domain(domain), other.to_domain(domain)

    def __add__(self, other):
        if not isinstance(other, DomainScalar):
            return NotImplemented
        self, other = self.unify(other)
        return self.new(self.element + other.element, self.domain)

    def __sub__(self, other):
        if not isinstance(other, DomainScalar):
            return NotImplemented
        self, other = self.unify(other)
        return self.new(self.element - other.element, self.domain)

    def __mul__(self, other):
        if not isinstance(other, DomainScalar):
            from .domainmatrix import DomainMatrix
            if isinstance(other, DomainMatrix):
                self, other = self.unify(other)
                return other.scalarmul(self)
            return NotImplemented

        self, other = self.unify(other)
        return self.new(self.element * other.element, self.domain)

    def __floordiv__(self, other):
        if not isinstance(other, DomainScalar):
            return NotImplemented
        self, other = self.unify(other)
        return self.new(self.domain.quo(self.element, other.element), self.domain)

    def __pow__(self, n):
        if not isinstance(n, int):
            return NotImplemented
        return self.new(self.element**n, self.domain)

    def __pos__(self):
        return self.new(+self.element, self.domain)

    def __eq__(self, other):
        if not isinstance(other, DomainScalar):
            return NotImplemented
        return self.element == other.element and self.domain == other.domain
