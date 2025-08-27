from __future__ import annotations

from typing import Protocol, Sequence, TypeVar
from sympy.polys.domains import Domain
from sympy.polys.domains.domain import Ef, Er
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
            return f"O({sym}**{prec})"
        else:
            return "0"

    if prec is not None:
        return f"{poly} + O({sym}**{prec})"
    return poly


class PowerSeriesRingProto(Protocol[TSeries, Er]):
    """A protocol for a power series ring."""

    def __init__(self, prec: int = 6, /) -> None:
        """Initialize a power series ring over a ring."""
        ...

    def __repr__(self, /) -> str:
        """Return string representation of the ring."""
        ...

    def __call__(self, coeffs: Sequence[Er | int], prec: int | None = None) -> TSeries:
        """Return a power series element over a ring."""
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

    def to_dense(self, s: TSeries, /) -> dup[Er]:
        """Return the coefficients of a power series as a dense list."""
        ...

    def series_prec(self, s: TSeries, /) -> int | None:
        """Return the precision of a power series."""
        ...

    def equal(self, s1: TSeries, s2: TSeries, /) -> bool | None:
        """Check if two power series are equal."""
        ...

    def equal_repr(self, s1: TSeries, s2: TSeries, /) -> bool:
        """Check if two power series have the same representation."""
        ...

    def is_ground(self, arg: TSeries) -> bool | None:
        """Check if a arg is a ground element of the power series ring."""
        ...

    def constant_coefficient(self, s: TSeries) -> Er:
        """Return the constant coefficient of a power series."""
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

    def add_ground(self, s: TSeries, n: Er, /) -> TSeries:
        """Add a ground element to a power series."""
        ...

    def subtract(self, s1: TSeries, s2: TSeries, /) -> TSeries:
        """Subtract two power series."""
        ...

    def subtract_ground(self, s: TSeries, n: Er, /) -> TSeries:
        """Subtract a ground element from a power series."""
        ...

    def rsubtract_ground(self, s: TSeries, n: Er, /) -> TSeries:
        """Subtract a power series from a ground element."""
        ...

    def multiply(self, s1: TSeries, s2: TSeries, /) -> TSeries:
        """Multiply two power series."""
        ...

    def multiply_ground(self, s: TSeries, n: Er, /) -> TSeries:
        """Multiply a power series by a ground element (integer or rational)."""
        ...

    def divide(self, s1: TSeries, s2: TSeries, /) -> TSeries:
        """Divide a power series by another power series."""
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

    def compose(self, s1: TSeries, s2: TSeries, /) -> TSeries:
        """Compose two power series."""
        ...

    def inverse(self, s: TSeries, /) -> TSeries:
        """Return the inverse of a power series."""
        ...

    def reversion(self, s: TSeries, /) -> TSeries:
        """Return the reversion of a power series."""
        ...

    def differentiate(self, s: TSeries, /) -> TSeries:
        """Return the derivative of a power series."""
        ...


class PowerSeriesRingFieldProto(
    PowerSeriesRingProto[TSeries, Ef], Protocol[TSeries, Ef]
):
    """A protocol for a power series ring over a field."""

    def __call__(
        self, coeffs: Sequence[Ef | int | tuple[int, int]], prec: int | None = None
    ) -> TSeries:
        """Return a power series element over a field."""
        ...

    def integrate(self, s: TSeries, /) -> TSeries:
        """Return the integral of a power series."""
        ...

    def sqrt(self, s: TSeries, /) -> TSeries:
        """Return the square root of a power seris."""
        ...

    def log(self, s: TSeries, /) -> TSeries:
        """Return the logarithm of a power series."""
        ...

    def log1p(self, s: TSeries, /) -> TSeries:
        """Return the logarithm of (1 + s) for a power series with constant term."""
        ...

    def exp(self, s: TSeries, /) -> TSeries:
        """Return the exponential of a power series."""
        ...

    def expm1(self, s: TSeries, /) -> TSeries:
        """Return the exponential of a power series minus 1."""
        ...

    def atan(self, s: TSeries, /) -> TSeries:
        """Return the arctangent of a power series."""
        ...

    def atanh(self, s: TSeries, /) -> TSeries:
        """Return the hyperbolic arctangent of a power series."""
        ...

    def asin(self, s: TSeries, /) -> TSeries:
        """Return the arcsine of a power series."""
        ...

    def asinh(self, s: TSeries, /) -> TSeries:
        """Return the hyperbolic arcsine of a power series."""
        ...

    def tan(self, s: TSeries, /) -> TSeries:
        """Return the tangent of a power series."""
        ...

    def tanh(self, s: TSeries, /) -> TSeries:
        """Return the hyperbolic tangent of a power series."""
        ...

    def sin(self, s: TSeries, /) -> TSeries:
        """Return the sine of a power series."""
        ...

    def sinh(self, s: TSeries, /) -> TSeries:
        """Return the hyperbolic sine of a power series."""
        ...

    def cos(self, s: TSeries, /) -> TSeries:
        """Return the cosine of a power series."""
        ...

    def cosh(self, s: TSeries, /) -> TSeries:
        """Return the hyperbolic cosine of a power series."""
        ...

    def hypot(self, s1: TSeries, s2: TSeries, /) -> TSeries:
        """Return the hypotenuse of two power series."""
        ...
