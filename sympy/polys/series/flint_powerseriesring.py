from __future__ import annotations

from contextlib import contextmanager
from typing import Any, TypeVar

from sympy.polys.domains import Domain, QQ, ZZ
from sympy.polys.series.powerseriesring import _series_from_list, PowerSeriesRing

try:
    from flint import ctx, fmpq_poly, fmpq_series, fmpz_poly, fmpz_series
except ImportError:
    pass

ZZSeries = TypeVar("ZZSeries", bound='fmpz_series | fmpz_poly')
QQSeries = TypeVar("QQSeries", bound='fmpq_series | fmpq_poly')


def _get_series_precision(s: ZZSeries | QQSeries) -> int:
    """Helper function to get the precision of a series. By using the
    representation of the ring"""

    # XXX: This approach is inefficient, but as of python-flint 0.7.1, there is
    # no alternative method to extract the precision from a series element.
    rep = s.repr()
    prec = int(rep.split('prec=')[1].split(')')[0].strip())
    return prec

@contextmanager
def _global_cap(cap: int):
    """Temporarily set the global series cap within a context."""
    if ctx is None:
        raise RuntimeError("Flint is not available")
    old_cap, ctx.cap = ctx.cap, cap
    try:
        yield
    finally:
        ctx.cap = old_cap

class FlintPowerSeriesRingZZ(PowerSeriesRing[ZZSeries]):
    """Flint implementation of power series ring over integer ring."""

    _domain = ZZ

    def __init__(self, prec: int = 6) -> None:
        self._prec = prec

    def __repr__(self) -> str:
        return f"Flint Power Series Ring over {self.domain} with precision {self.prec}"

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, FlintPowerSeriesRingZZ):
            return NotImplemented
        return self.prec == other.prec

    def __hash__(self) -> int:
        return hash((self.domain, self.prec))

    @property
    def domain(self) -> Domain:
        """Return the ground domain of the power series ring."""
        return self._domain

    @property
    def prec(self) -> int:
        """Return the default precision for power series operations."""
        return self._prec

    @property
    def one(self) -> ZZSeries:
        if self.prec == 0:
            return fmpz_series([1], prec=0)
        return fmpz_poly([1])

    @property
    def zero(self) -> ZZSeries:
        if self.prec == 0:
            return fmpz_series([0], prec=0)
        return fmpz_poly([0])

    @property
    def gen(self) -> ZZSeries:
        if self.prec < 2:
            return fmpz_series([0, 1], prec=self.prec)
        return fmpz_poly([0, 1])

    def pretty(self, s: ZZSeries) -> str:
        """Pretty print a power series with improved formatting."""
        if isinstance(s, fmpz_poly):
            return _series_from_list(s.coeffs(), None)

        prec = _get_series_precision(s)
        return _series_from_list(s.coeffs(), prec)

    def print(self, s: ZZSeries) -> str:
        return self.pretty(s)

    def from_list(self, coeffs: list[Any], prec: int | None = None) -> ZZSeries:
        """Create a power series from a list of coefficients."""
        if prec is None:
            if len(coeffs) < self.prec:
                return fmpz_poly(coeffs)
            prec = self.prec

        return fmpz_series(coeffs, prec=prec)

    def to_list(self, s: ZZSeries) -> list[Any]:
        """Returns coeffs list."""
        return s.coeffs()

    def equal(self, s1: ZZSeries, s2: ZZSeries) -> bool | None:
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            return s1 == s2
        elif isinstance(s1, fmpz_series) and isinstance(s2, fmpz_series):
            if s1._equal_repr(s2):
                return None
        return False

    def negative(self, s: ZZSeries) -> ZZSeries:
        """Return the negative of a power series."""
        return -s

    def add(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Add two power series."""
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            poly = s1 + s2
            if poly.degree() < self.prec:
                return poly
            return fmpz_series(poly, prec=self.prec)

        with _global_cap(self.prec):
            return s1 + s2

    def subtract(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Subtract two power series."""
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            poly = s1 - s2
            if poly.degree() < self.prec:
                return poly
            return fmpz_series(poly, prec=self.prec)

        with _global_cap(self.prec):
            return s1 - s2

    def multiply(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Multiply two power series."""
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            poly = s1 * s2
            if poly.degree() < self.prec:
                return poly
            return fmpz_series(poly, prec=self.prec)

        with _global_cap(self.prec):
            return s1 * s2

    def multiply_ground(self, s1: ZZSeries, n: Any) -> ZZSeries:
        """Multiply a power series by an integer or a polynomial."""
        if isinstance(s1, fmpz_poly):
            poly = s1 * n
            if poly.degree() < self.prec:
                return poly
            return fmpz_series(poly, prec=self.prec)

        with _global_cap(self.prec):
            return s1 * n

    def pow_int(self, s: ZZSeries, n: int) -> ZZSeries:
        """Raise a power series to an integer power."""
        if n < 0:
            raise ValueError("Power must be non-negative")

        if isinstance(s, fmpz_poly):
            poly = s ** n
            if poly.degree() < self.prec:
                return poly
            return fmpz_series(poly, prec=self.prec)

        with _global_cap(self.prec):
            return s ** n

    def square(self, s: ZZSeries) -> ZZSeries:
        """Return the square of a power series."""
        return self.pow_int(s, 2)

    def truncate(self, s: ZZSeries, n: int) -> ZZSeries:
        """Truncate a power series to the first n terms."""
        if n < 0:
            raise ValueError("Truncation precision must be non-negative")

        if len(s) < n:
            return s

        coeffs = s.coeffs()[:n]
        with _global_cap(self.prec):
            return fmpz_series(coeffs, prec=n)


class FlintPowerSeriesRingQQ(PowerSeriesRing[QQSeries]):
    """Flint implementation of power series ring over rational fields."""

    _domain = QQ

    def __init__(self, prec: int = 6) -> None:
        self._prec = prec

    def __repr__(self) -> str:
        return f"Flint Power Series Ring over {self.domain} with precision {self.prec}"

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, FlintPowerSeriesRingQQ):
            return NotImplemented
        return self.prec == other.prec

    def __hash__(self) -> int:
        return hash((self.domain, self.prec))

    @property
    def domain(self) -> Domain:
        """Return the ground domain of the power series ring."""
        return self._domain

    @property
    def prec(self) -> int:
        """Return the default precision for power series operations."""
        return self._prec

    @property
    def one(self) -> QQSeries:
        if self.prec == 0:
            return fmpq_series([1], prec=0)
        return fmpq_poly([1])

    @property
    def zero(self) -> QQSeries:
        if self.prec == 0:
            return fmpq_series([0], prec=0)
        return fmpq_poly([0])

    @property
    def gen(self) -> QQSeries:
        if self.prec < 2:
            return fmpq_series([0, 1], prec=self.prec)
        return fmpq_poly([0, 1])

    def pretty(self, s: QQSeries) -> str:
        """Pretty print a power series with improved formatting."""
        if isinstance(s, fmpq_poly):
            return _series_from_list(s.coeffs(), None)

        prec = _get_series_precision(s)
        return _series_from_list(s.coeffs(), prec)

    def print(self, s: QQSeries) -> str:
        return self.pretty(s)

    def from_list(self, coeffs: list[Any], prec: int | None = None) -> QQSeries:
        """Create a power series from a list of coefficients."""
        if prec is None:
            if len(coeffs) < self.prec:
                return fmpq_poly(coeffs)
            prec = self.prec

        return fmpq_series(coeffs, prec=prec)

    def to_list(self, s: QQSeries) -> list[Any]:
        """Returns coeffs list."""
        return s.coeffs()

    def equal(self, s1: QQSeries, s2: QQSeries) -> bool | None:
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpq_poly):
            return s1 == s2
        elif isinstance(s1, fmpz_series) and isinstance(s2, fmpz_series):
            if s1._equal_repr(s2):
                return None
        return False

    def negative(self, s: QQSeries) -> QQSeries:
        """Return the negative of a power series."""
        return -s

    def add(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Add two power series."""
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            poly = s1 + s2
            if poly.degree() < self.prec:
                return poly
            return fmpq_series(poly, prec=self.prec)

        with _global_cap(self.prec):
            return s1 + s2

    def subtract(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Subtract two power series."""
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            poly = s1 - s2
            if poly.degree() < self.prec:
                return poly
            return fmpq_series(poly, prec=self.prec)

        with _global_cap(self.prec):
            return s1 - s2

    def multiply(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Multiply two power series."""
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            poly = s1 * s2
            if poly.degree() < self.prec:
                return poly
            return fmpq_series(poly, prec=self.prec)

        with _global_cap(self.prec):
            return s1 * s2

    def multiply_ground(self, s1: QQSeries, n: Any) -> QQSeries:
        """Multiply a power series by a rational number."""
        if isinstance(s1, fmpq_poly):
            poly = s1 * n
            if poly.degree() < self.prec:
                return poly
            return fmpq_series(poly, prec=self.prec)

        with _global_cap(self.prec):
            return s1 * n

    def pow_int(self, s: QQSeries, n: int) -> QQSeries:
        """Raise a power series to an integer power."""
        if n < 0:
            raise ValueError("Power must be non-negative")

        if isinstance(s, fmpq_poly):
            poly = s ** n
            if poly.degree() < self.prec:
                return poly
            return fmpq_series(poly, prec=self.prec)

        with _global_cap(self.prec):
            return s ** n

    def square(self, s: QQSeries) -> QQSeries:
        """Return the square of a power series."""
        return self.pow_int(s, 2)

    def truncate(self, s: QQSeries, n: int) -> QQSeries:
        """Truncate a power series to the first n terms."""
        if n < 0:
            raise ValueError("Truncation precision must be non-negative")

        if len(s) < n:
            return s

        coeffs = s.coeffs()[:n]
        with _global_cap(self.prec):
            return fmpq_series(coeffs, prec=n)
