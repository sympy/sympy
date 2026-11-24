"""

Module for the DomainScalar class.

A DomainScalar represents an element which is in a particular
Domain. The idea is that the DomainScalar class provides the
convenience routines for unifying elements with different domains.

It assists in Scalar Multiplication and getitem for DomainMatrix.

"""

from __future__ import annotations
from typing import TYPE_CHECKING, Generic, overload
from ..constructor import construct_domain

from sympy.polys.domains import Domain, ZZ
from sympy.polys.domains.domain import Er, Es

if TYPE_CHECKING:
    from typing import Any
    from sympy.core.expr import Expr


class DomainScalar(Generic[Er]):
    r"""
    docstring
    """

    element: Er
    domain: Domain[Er]

    def __new__(cls, element: Er, domain: Domain[Er]) -> DomainScalar[Er]:
        if not isinstance(domain, Domain):
            raise TypeError("domain should be of type Domain")
        if not domain.of_type(element):
            raise TypeError("element %s should be in domain %s" % (element, domain))
        return cls.new(element, domain)

    @classmethod
    def new(cls, element: Es, domain: Domain[Es]) -> DomainScalar[Es]:
        obj: DomainScalar[Es] = object.__new__(DomainScalar)
        obj.element = element
        obj.domain = domain
        return obj

    def __repr__(self) -> str:
        return repr(self.element)

    @classmethod
    def from_sympy(cls, expr: Expr | complex) -> DomainScalar[Any]:
        [domain, [element]] = construct_domain([expr])
        return cls.new(element, domain)

    def to_sympy(self) -> Expr:
        return self.domain.to_sympy(self.element)

    def to_domain(self, domain: Domain[Es]) -> DomainScalar[Es]:
        element = domain.convert_from(self.element, self.domain)
        return self.new(element, domain)

    def convert_to(self, domain: Domain[Es]) -> DomainScalar[Es]:
        return self.to_domain(domain)

    @overload
    def unify(
        self, other: DomainScalar[Er]
    ) -> tuple[DomainScalar[Er], DomainScalar[Er]]: ...
    @overload
    def unify(
        self, other: DomainScalar[Any]
    ) -> tuple[DomainScalar[Any], DomainScalar[Any]]: ...

    def unify(
        self, other: DomainScalar[Any]
    ) -> tuple[DomainScalar[Any], DomainScalar[Any]]:
        domain = self.domain.unify(other.domain)
        return self.to_domain(domain), other.to_domain(domain)

    def __bool__(self) -> bool:
        return bool(self.element)

    @overload
    def __add__(self, other: DomainScalar[Er]) -> DomainScalar[Er]: ...
    @overload
    def __add__(self, other: DomainScalar[Any]) -> DomainScalar[Any]: ...

    def __add__(self, other: object) -> DomainScalar[Any]:
        if not isinstance(other, DomainScalar):
            return NotImplemented
        self, other = self.unify(other)
        return self.new(self.element + other.element, self.domain)

    @overload
    def __sub__(self, other: DomainScalar[Er]) -> DomainScalar[Er]: ...
    @overload
    def __sub__(self, other: DomainScalar[Any]) -> DomainScalar[Any]: ...

    def __sub__(self, other: object) -> DomainScalar[Any]:
        if not isinstance(other, DomainScalar):
            return NotImplemented
        self, other = self.unify(other)
        return self.new(self.element - other.element, self.domain)

    @overload
    def __mul__(self, other: DomainScalar[Er] | int) -> DomainScalar[Er]: ...
    @overload
    def __mul__(self, other: DomainScalar[Any]) -> DomainScalar[Any]: ...

    def __mul__(self, other: object) -> DomainScalar[Any]:
        if not isinstance(other, DomainScalar):
            if isinstance(other, int):
                other = DomainScalar(ZZ(other), ZZ)
            else:
                return NotImplemented

        self, other = self.unify(other)
        return self.new(self.element * other.element, self.domain)

    @overload
    def __floordiv__(self, other: DomainScalar[Er]) -> DomainScalar[Er]: ...
    @overload
    def __floordiv__(self, other: DomainScalar[Any]) -> DomainScalar[Any]: ...

    def __floordiv__(self, other: object) -> DomainScalar[Any]:
        if not isinstance(other, DomainScalar):
            return NotImplemented
        self, other = self.unify(other)
        return self.new(self.domain.quo(self.element, other.element), self.domain)

    @overload
    def __mod__(self, other: DomainScalar[Er]) -> DomainScalar[Er]: ...
    @overload
    def __mod__(self, other: DomainScalar[Any]) -> DomainScalar[Any]: ...

    def __mod__(self, other: object) -> DomainScalar[Any]:
        if not isinstance(other, DomainScalar):
            return NotImplemented
        self, other = self.unify(other)
        return self.new(self.domain.rem(self.element, other.element), self.domain)

    @overload
    def __divmod__(
        self, other: DomainScalar[Er]
    ) -> tuple[DomainScalar[Er], DomainScalar[Er]]: ...
    @overload
    def __divmod__(
        self, other: DomainScalar[Any]
    ) -> tuple[DomainScalar[Any], DomainScalar[Any]]: ...

    def __divmod__(self, other: object) -> tuple[DomainScalar[Any], DomainScalar[Any]]:
        if not isinstance(other, DomainScalar):
            return NotImplemented
        self, other = self.unify(other)
        q, r = self.domain.div(self.element, other.element)
        return (self.new(q, self.domain), self.new(r, self.domain))

    def __pow__(self, n: int) -> DomainScalar[Er]:
        if not isinstance(n, int):
            return NotImplemented
        return self.new(self.element**n, self.domain)

    def __pos__(self) -> DomainScalar[Er]:
        return self.new(+self.element, self.domain)

    def __neg__(self) -> DomainScalar[Er]:
        return self.new(-self.element, self.domain)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, DomainScalar):
            return NotImplemented
        return self.element == other.element and self.domain == other.domain

    def is_zero(self) -> bool:
        return self.element == self.domain.zero

    def is_one(self) -> bool:
        return self.element == self.domain.one
