from __future__ import annotations

from functools import cached_property

from sympy.core.expr import Expr
from sympy.polys.domains.ring import Ring
from sympy.polys.domains.compositedomain import CompositeDomain
from sympy.polys.domains.domain import Er, Domain
from sympy.polys.series.ring import (
    power_series_ring,
    PowerSeriesElement,
)

from typing import TYPE_CHECKING, Protocol

if TYPE_CHECKING:
    from typing import TypeIs
    from sympy.polys.densebasic import dup
    from sympy.polys.series.tring import TSeriesElement
    from sympy.polys.series.base import PowerSeriesRingProto


class SeriesRingProto(Protocol[Er]):
    """Generic protocol for power series rings."""

    ring: PowerSeriesRingProto[TSeriesElement[Er], Er]
    dtype: type[PowerSeriesElement[Er]]
    domain: Domain[Er]
    symbol: Expr
    prec: int

    one: PowerSeriesElement[Er]
    zero: PowerSeriesElement[Er]
    gen: PowerSeriesElement[Er]

    def __init__(self, domain: Domain[Er], symbol: str | Expr = "x", prec: int = 6): ...

    def __repr__(self) -> str: ...

    def __eq__(self, other) -> bool: ...

    def __hash__(self) -> int: ...

    def is_element(self, element: object) -> TypeIs[PowerSeriesElement[Er]]: ...

    def order_term(self) -> PowerSeriesElement[Er]: ...

    def from_expr(self, expr: Expr) -> PowerSeriesElement[Er]: ...

    def from_list(
        self, lst: list[Er], prec: int | None = None
    ) -> PowerSeriesElement[Er]: ...

    def from_element(self, element: TSeriesElement[Er]) -> PowerSeriesElement[Er]: ...

    def from_int(self, arg: int) -> PowerSeriesElement[Er]: ...

    def from_ground(self, arg: Er) -> PowerSeriesElement[Er]: ...

    def to_expr(self, element: PowerSeriesElement[Er]) -> Expr: ...

    def to_list(self, element: PowerSeriesElement[Er]) -> list[Er]: ...

    def to_dense(self, element: PowerSeriesElement[Er]) -> dup[Er]: ...

    def domain_new(self, arg: Er | int) -> Er: ...

    def ring_new(self, arg: Expr | Er | int) -> PowerSeriesElement[Er]: ...


class PowerSeriesRing(Ring[PowerSeriesElement[Er]], CompositeDomain):
    """A Domain class for representing univariate power series rings.

    Notes
    =====

    This class is at experimental stage. Proper domain methods should be added to
    integrate with SymPy's existing domain framework.
    """

    is_PowerSeriesRing = is_Series = True

    has_assoc_Ring = True
    has_assoc_Field = False

    ring: SeriesRingProto[Er]
    dtype: type[PowerSeriesElement[Er]]
    gen: PowerSeriesElement[Er]
    symbol: Expr
    domain: Domain[Er]

    def __init__(self, domain: Domain[Er], symbol: Expr | str = "x", prec: int = 6):
        ring, gen = power_series_ring(str(symbol), domain, prec)

        self.ring = ring
        self.gen = gen

        self.dtype = ring.dtype
        self.domain = ring.domain
        self.symbol = ring.symbol

    def __eq__(self, other: object) -> bool:
        if isinstance(other, PowerSeriesRing):
            return self.ring == other.ring
        else:
            return NotImplemented

    def __hash__(self) -> int:
        return hash((self.__class__.__name__, self.ring, self.domain, self.symbol))

    def __repr__(self) -> str:
        return f"{self.domain}[[{self.symbol}], {self.prec}]"

    def new(  # type: ignore
        self, element: Expr | Er | int
    ) -> PowerSeriesElement[Er]:
        return self.ring.ring_new(element)

    def of_type(self, element) -> TypeIs[PowerSeriesElement[Er]]:
        """Check if ``a`` is of type ``dtype``."""
        return self.ring.is_element(element)

    @cached_property
    def zero(self) -> PowerSeriesElement[Er]:  # type: ignore
        return self.ring.zero

    @cached_property
    def one(self) -> PowerSeriesElement[Er]:  # type: ignore
        return self.ring.one

    @cached_property
    def prec(self) -> int:
        return self.ring.prec

    def is_unit(self, a: PowerSeriesElement[Er]) -> bool:
        """Returns ``True`` if ``constant coefficient`` of series is a unit of ``self``"""
        if not a.is_ground:
            return False
        K = self.domain
        return K.is_unit(a.constant_coefficient())

    def canonical_unit(self, a: PowerSeriesElement[Er]) -> PowerSeriesElement[Er]:
        u = self.domain.canonical_unit(a.constant_coefficient())
        return self.ring.from_ground(u)

    def to_sympy(self, a: PowerSeriesElement[Er]) -> Expr:
        """Convert `a` to a SymPy object."""
        return self.ring.to_expr(a)

    def from_sympy(self, a: Expr) -> PowerSeriesElement[Er]:
        """Convert SymPy's expression to `dtype`."""
        return self.ring.from_expr(a)

    def from_ZZ(K1, a, K0):
        """Convert a Python `int` object to `dtype`."""
        return K1.ring.from_ground(K1.domain.convert(a, K0))

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`."""
        return K1.ring.from_ground(K1.domain.convert(a, K0))

    def from_QQ(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`."""
        return K1.ring.from_ground(K1.domain.convert(a, K0))

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`."""
        return K1.ring.from_ground(K1.domain.convert(a, K0))

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`."""
        return K1.ring.from_ground(K1.domain.convert(a, K0))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`."""
        return K1.ring.from_ground(K1.domain.convert(a, K0))

    def is_positive(self, a: PowerSeriesElement[Er]) -> bool:
        """Returns True if `constant coefficient(a)` is positive."""
        c = a.constant_coefficient()
        return self.domain.is_positive(c)

    def is_negative(self, a: PowerSeriesElement[Er]) -> bool:
        """Returns True if `constant coefficient(a)` is negative."""
        c = a.constant_coefficient()
        return self.domain.is_negative(c)

    def is_nonpositive(self, a: PowerSeriesElement[Er]) -> bool:
        """Returns True if `constant coefficient(a)` is non-positive."""
        c = a.constant_coefficient()
        return self.domain.is_nonpositive(c)

    def is_nonnegative(self, a: PowerSeriesElement[Er]) -> bool:
        """Returns True if `constant coefficient(a)` is non-negative."""
        c = a.constant_coefficient()
        return self.domain.is_nonnegative(c)
