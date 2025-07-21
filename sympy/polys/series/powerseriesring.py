from __future__ import annotations

from typing import Protocol, TypeVar
from sympy.core.symbol import Symbol
from sympy.polys.domains import Domain
from sympy.polys.domains.domain import Er


_x: Symbol = Symbol("x")
TSeries = TypeVar("TSeries")


def series_pprint(series: list[Er], prec: int | None) -> str:
    """Convert a list of coefficients into a string representation of a power series."""
    terms = []

    if prec == 0:
        return f"O({_x}**{prec})"

    for i, coeff in enumerate(series):
        if coeff == 0:
            continue
        if coeff == 1 and i != 0:
            term = f"{_x}" if i == 1 else f"{_x}**{i}"
        elif coeff == -1 and i != 0:
            term = f"-{_x}" if i == 1 else f"-{_x}**{i}"
        else:
            if i == 0:
                term = f"{coeff}"
            elif i == 1:
                term = f"{coeff}*{_x}"
            else:
                term = f"{coeff}*{_x}**{i}"
        terms.append((i, term))

    if not terms:
        if prec is not None:
            return f"0 + O({_x}**{prec})"
        else:
            return "0"

    poly = " + ".join(term for _, term in terms).replace("+ -", "- ")
    if prec is not None:
        return f"{poly} + O({_x}**{prec})"
    return poly


class PowerSeriesRing(Protocol[TSeries, Er]):
    """A protocol for a power series ring."""

    def __init__(self, prec: int = 6, /) -> None:
        """Initialize a power series ring over ZZ."""
        ...

    def __repr__(self, /) -> str:
        """Return string representation of the ring."""
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

    def pretty(self, s: TSeries, /) -> str:
        """Return a pretty string representation of a power series."""
        ...

    def print(self, s: TSeries, /) -> None:
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
