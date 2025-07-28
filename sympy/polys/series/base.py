from __future__ import annotations

from typing import Protocol, Sequence, TypeVar
from sympy.polys.domains import Domain
from sympy.polys.domains.domain import Er
from sympy.polys.densebasic import dup, dup_pretty


TSeries = TypeVar("TSeries")


def series_pprint(
    series: dup[Er],
    prec: int | None,
    *,
    sym: str = "x",
    ascending: bool = True,
) -> str:
    """Convert a list of coefficients into a string representation of a power series."""
    if prec == 0:
        return f"O({sym}**{prec})"

    poly = dup_pretty(series, sym, ascending=ascending)

    if not poly or poly == "0":
        if prec is not None:
            return f"0 + O({sym}**{prec})"
        else:
            return "0"

    if prec is not None:
        return f"{poly} + O({sym}**{prec})"
    return poly


class PowerSeriesRingProto(Protocol[TSeries, Er]):
    """A protocol for a power series ring."""

    def __init__(self, prec: int = 6, /) -> None:
        """Initialize a power series ring over ZZ."""
        ...

    def __repr__(self, /) -> str:
        """Return string representation of the ring."""
        ...

    def __call__(self, coeffs: Sequence, prec: int | None = None) -> TSeries:
        """Return a power series element."""
        ...

    @property
    def domain(self, /) -> Domain[Er]:
        """Return the ground domain of the power series ring."""
        ...

    @property
    def prec(self, /) -> int:
        """Return the default precision for power series operations."""
        ...

    @property
    def one(self, /) -> TSeries:
        """Return the multiplicative identity (1) as a power series."""
        ...

    @property
    def zero(self, /) -> TSeries:
        """Return the additive identity (0) as a power series."""
        ...

    @property
    def gen(self, /) -> TSeries:
        """Return the generator (x) as a power series."""
        ...

    def pretty(
        self, s: TSeries, /, *, symbol: str = "x", ascending: bool = True
    ) -> str:
        """Return a pretty string representation of a power series."""
        ...

    def print(
        self, s: TSeries, /, *, symbol: str = "x", ascending: bool = True
    ) -> None:
        """Return a printable string representation of a power series."""
        ...

    def from_list(self, coeffs: list[Er], prec: int | None = None, /) -> TSeries:
        """Create a power series from a list of coefficients."""
        ...

    def to_list(self, s: TSeries, /) -> list[Er]:
        """Return the coefficients of a power series as a list."""
        ...

    def equal(self, s1: TSeries, s2: TSeries, /) -> bool | None:
        """Check if two power series are equal."""
        ...

    def equal_repr(self, s1: TSeries, s2: TSeries, /) -> bool:
        """Check if two power series have the same representation."""
        ...

    def positive(self, s: TSeries, /) -> TSeries:
        """Return the positive of a power series."""
        ...

    def negative(self, s: TSeries, /) -> TSeries:
        """Return the additive inverse of a power series."""
        ...

    def add(self, s1: TSeries, s2: TSeries, /) -> TSeries:
        """Add two power series."""
        ...

    def subtract(self, s1: TSeries, s2: TSeries, /) -> TSeries:
        """Subtract two power series."""
        ...

    def multiply(self, s1: TSeries, s2: TSeries, /) -> TSeries:
        """Multiply two power series."""
        ...

    def multiply_ground(self, s: TSeries, n: Er, /) -> TSeries:
        """Multiply a power series by a ground element (integer or rational)."""
        ...

    def pow_int(self, s: TSeries, n: int, /) -> TSeries:
        """Raise a power series to an integer power."""
        ...

    def square(self, s: TSeries, /) -> TSeries:
        """Return the square of a power series."""
        ...

    def truncate(self, s: TSeries, n: int, /) -> TSeries:
        """Truncate a power series to the first n terms."""
        ...
