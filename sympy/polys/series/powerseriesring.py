from __future__ import annotations

from typing import Any, Protocol, TypeVar
from sympy.core.symbol import Symbol
from sympy.polys.domains import Domain


_x: Symbol = Symbol("x")
Sr = TypeVar("Sr")


def _series_from_list(series: list[Any], prec: int | None) -> str:
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


class PowerSeriesRing(Protocol[Sr]):
    """Protocol for power series ring over the integer ring (ZZ).

    Attributes
    ----------
    domain : Domain
        The ground domain (ZZ).
    prec : int
        The default precision for power series operations.
    """

    def __init__(self, prec: int = 6, /) -> None:
        """Initialize a power series ring over ZZ."""
        ...

    def __repr__(self, /) -> str:
        """Return string representation of the ring."""
        ...

    @property
    def domain(self, /) -> Domain:
        """Return the ground domain of the power series ring."""
        ...

    @property
    def prec(self, /) -> int:
        """Return the default precision for power series operations."""
        ...

    @property
    def one(self, /) -> Sr:
        """Return the multiplicative identity (1) as a power series."""
        ...

    @property
    def zero(self, /) -> Sr:
        """Return the additive identity (0) as a power series."""
        ...

    @property
    def gen(self, /) -> Sr:
        """Return the generator (x) as a power series."""
        ...

    def pretty(self, s: Sr, /) -> str:
        """Return a pretty string representation of a power series."""
        ...

    def print(self, s: Sr, /) -> None:
        """Return a printable string representation of a power series."""
        ...

    def from_list(self, coeffs: list[Any], prec: int | None = None, /) -> Sr:
        """Create a power series from a list of coefficients."""
        ...

    def to_list(self, s: Sr, /) -> list[Any]:
        """Return the coefficients of a power series as a list."""
        ...

    def equal(self, s1: Sr, s2: Sr, /) -> bool | None:
        """Check if two power series are equal."""
        ...

    def equal_repr(self, s1: Sr, s2: Sr, /) -> bool:
        """Check if two power series have the same representation."""
        ...

    def positive(self, s: Sr, /) -> Sr:
        """Return the positive of a power series."""
        ...

    def negative(self, s: Sr, /) -> Sr:
        """Return the additive inverse of a power series."""
        ...

    def add(self, s1: Sr, s2: Sr, /) -> Sr:
        """Add two power series."""
        ...

    def subtract(self, s1: Sr, s2: Sr, /) -> Sr:
        """Subtract two power series."""
        ...

    def multiply(self, s1: Sr, s2: Sr, /) -> Sr:
        """Multiply two power series."""
        ...

    def multiply_ground(self, s: Sr, n: Any, /) -> Sr:
        """Multiply a power series by a ground element (integer or rational)."""
        ...

    def pow_int(self, s: Sr, n: int, /) -> Sr:
        """Raise a power series to an integer power."""
        ...

    def square(self, s: Sr, /) -> Sr:
        """Return the square of a power series."""
        ...

    def truncate(self, s: Sr, n: int, /) -> Sr:
        """Truncate a power series to the first n terms."""
        ...
