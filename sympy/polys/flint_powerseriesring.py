from __future__ import annotations

from contextlib import contextmanager
from typing import Any, TypeAlias

from sympy.external.gmpy import GROUND_TYPES
from sympy.polys.domains import QQ, ZZ
from sympy.polys.powerseriesring import (PowerSeriesRingQQ, PowerSeriesRingZZ,
                                        series_from_list)

if GROUND_TYPES == 'flint':
    from flint import ctx, fmpq_poly, fmpq_series, fmpz_poly, fmpz_series
else:
    ctx = None
    fmpq_poly = None
    fmpq_series = None
    fmpz_poly = None
    fmpz_series = None


ZZSeries: TypeAlias = 'fmpz_series | fmpz_poly'
QQSeries: TypeAlias = 'fmpq_series | fmpq_poly'


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

class FlintPowerSeriesRingZZ(PowerSeriesRingZZ):
    """Flint implementation of power series ring over integer ring."""

    def __init__(self, prec: int = 6) -> None:
        self.domain = ZZ
        self.prec = prec

    def __repr__(self) -> str:
        return f"Flint Power Series Ring over {self.domain} with precision {self.prec}"

    def __eq__(self, other: Any) -> bool:
        return (isinstance(other, FlintPowerSeriesRingZZ) and
                self.prec == other.prec)

    def __ne__(self, other: Any) -> bool:
        return not self.__eq__(other)

    def __hash__(self) -> int:
        return hash((self.domain, self.prec))

    @property
    def one(self) -> ZZSeries:
        return fmpz_poly([1])

    @property
    def zero(self) -> ZZSeries:
        return fmpz_poly([0])

    @property
    def gen(self) -> ZZSeries:
        return fmpz_poly([0, 1])

    def pretty(self, s: ZZSeries) -> str:
        """Pretty print a power series with improved formatting."""
        if isinstance(s, fmpz_poly):
            return series_from_list(s.coeffs(), None)

        prec = int(s.repr().split('prec=')[1].split(')')[0].strip())
        return series_from_list(s.coeffs(), prec)

    def print(self, s: ZZSeries) -> str:
        return self.pretty(s)

    def equal(self, s1: ZZSeries, s2: ZZSeries) -> bool:
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            return s1==s2
        elif isinstance(s1, fmpz_series) and isinstance(s2, fmpz_series):
            return s1._equal_repr(s2)

        return False

    def negative(self, s: ZZSeries) -> ZZSeries:
        """Return the negative of a power series."""
        return -s

    def add(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Add two power series."""
        with _global_cap(self.prec):
            if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
                poly = s1 + s2
                if poly.degree() < self.prec:
                    return poly
                return fmpz_series(poly, prec=self.prec)

            return s1 + s2

    def subtract(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Subtract two power series."""
        with _global_cap(self.prec):
            if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
                poly = s1 - s2
                if poly.degree() < self.prec:
                    return poly
                return fmpz_series(poly, prec=self.prec)

            return s1 - s2

    def multiply(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Multiply two power series."""
        with _global_cap(self.prec):
            if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
                poly = s1 * s2
                if poly.degree() < self.prec:
                    return poly
                return fmpz_series(poly, prec=self.prec)

            return s1 * s2

    def multiply_ground(self, s1: ZZSeries, n: Any) -> ZZSeries:
        """Multiply a power series by an integer or a polynomial."""
        with _global_cap(self.prec):
            if isinstance(s1, fmpz_poly):
                poly = s1 * n
                if poly.degree() < self.prec:
                    return poly
                return fmpz_series(poly, prec=self.prec)

            return s1 * n

    def trunc(self, s: ZZSeries, n: int) -> ZZSeries:
        """Truncate a power series to the first n terms."""
        if n < 0:
            raise ValueError("Truncation precision must be non-negative")
        with _global_cap(self.prec):
            coeffs = s.coeffs()[:n]

            if isinstance(s, fmpz_poly):
                return fmpz_poly(coeffs)
            return fmpz_series(coeffs, prec=n)


class FlintPowerSeriesRingQQ(PowerSeriesRingQQ):
    """Flint implementation of power series ring over rational fields."""

    def __init__(self, prec: int = 6) -> None:
        self.domain = QQ
        self.prec = prec

    def __repr__(self) -> str:
        return f"Flint Power Series Ring over {self.domain} with precision {self.prec}"

    def __eq__(self, other: Any) -> bool:
        return (isinstance(other, FlintPowerSeriesRingQQ) and
                self.prec == other.prec)

    def __ne__(self, other: Any) -> bool:
        return not self.__eq__(other)

    def __hash__(self) -> int:
        return hash((self.domain, self.prec))

    @property
    def one(self) -> QQSeries:
        return fmpq_poly([1])

    @property
    def zero(self) -> QQSeries:
        return fmpq_poly([0])

    @property
    def gen(self) -> QQSeries:
        return fmpq_poly([0, 1])

    def pretty(self, s: QQSeries) -> str:
        """Pretty print a power series with improved formatting."""
        if isinstance(s, fmpq_poly):
            return series_from_list(s.coeffs(), None)

        prec = int(s.repr().split('prec=')[1].split(')')[0].strip())
        return series_from_list(s.coeffs(), prec)

    def print(self, s: QQSeries) -> str:
        return self.pretty(s)

    def equal(self, s1: QQSeries, s2: QQSeries) -> bool:
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            return s1==s2
        elif isinstance(s1, fmpq_series) and isinstance(s2, fmpq_series):
            return s1._equal_repr(s2)

        return False

    def negative(self, s: QQSeries) -> QQSeries:
        """Return the negative of a power series."""
        return -s

    def add(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Add two power series."""
        with _global_cap(self.prec):
            if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
                poly = s1 + s2
                if poly.degree() < self.prec:
                    return poly
                return fmpq_series(poly, prec=self.prec)

            return s1 + s2

    def subtract(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Subtract two power series."""
        with _global_cap(self.prec):
            if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
                poly = s1 - s2
                if poly.degree() < self.prec:
                    return poly
                return fmpq_series(poly, prec=self.prec)

            return s1 - s2

    def multiply(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Multiply two power series."""
        with _global_cap(self.prec):
            if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
                poly = s1 * s2
                if poly.degree() < self.prec:
                    return poly
                return fmpq_series(poly, prec=self.prec)

            return s1 * s2

    def multiply_ground(self, s1: QQSeries, n: Any) -> QQSeries:
        """Multiply a power series by a rational number."""
        with _global_cap(self.prec):
            if isinstance(s1, fmpq_poly):
                poly = s1 * n
                if poly.degree() < self.prec:
                    return poly
                return fmpq_series(poly, prec=self.prec)

            return s1 * n

    def trunc(self, s: QQSeries, n: int) -> QQSeries:
        """Truncate a power series to the first n terms."""
        if n < 0:
            raise ValueError("Truncation precision must be non-negative")
        with _global_cap(self.prec):
            coeffs = s.coeffs()[:n]

            if isinstance(s, fmpq_poly):
                return fmpq_poly(coeffs)
            return fmpq_series(coeffs, prec=n)
