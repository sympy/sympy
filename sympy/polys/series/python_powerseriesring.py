from __future__ import annotations

from typing import Any, TypeVar
from sympy.polys.densearith import (dup_add, dup_mul, dup_mul_ground, dup_neg, dup_pow,
                                dup_sub)
from sympy.polys.densebasic import dup_degree, dup_reverse, dup_slice
from sympy.polys.domains import Domain, QQ, ZZ
from sympy.polys.series.powerseriesring import _series_from_list, PowerSeriesRing


T = TypeVar("T")
DUP = list[T]
USeries = tuple[DUP[T], int | None]

def _unify_prec(
    s1: USeries[T], s2: USeries[T], dom: Domain, ring_prec: int
) -> tuple[DUP[T], DUP[T], int]:
    """Unify the precision of two series."""
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if prec1 is None and prec2 is None:
        unified_prec = ring_prec
    elif prec1 is None:
        unified_prec = prec2
    elif prec2 is None:
        unified_prec = prec1
    else:
        unified_prec = min(prec1, prec2)

    coeffs1 = dup_slice(coeffs1, 0, unified_prec, dom)
    coeffs2 = dup_slice(coeffs2, 0, unified_prec, dom)

    return coeffs1, coeffs2, unified_prec

def _useries_equality(
    s1: USeries[T], s2: USeries[T]
) -> bool | None:
    """Check if two power series are equal."""
    if s1 == s2:
        return None
    return False

def _useries_add(
    s1: USeries[T], s2: USeries[T], dom: Domain, ring_prec: int
) -> USeries[T]:
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2
    max_degree = max(dup_degree(coeffs1), dup_degree(coeffs2))

    if prec1 is None and prec2 is None and (max_degree < ring_prec):
        series = dup_add(coeffs1, coeffs2, dom)

        return series, None

    coeffs1, coeffs2, prec = _unify_prec(s1, s2, dom, ring_prec)
    return dup_add(coeffs1, coeffs2, dom), prec

def _useries_sub(
    s1: USeries[T], s2: USeries[T], dom: Domain, ring_prec: int
) -> USeries[T]:
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2
    max_degree = max(dup_degree(coeffs1), dup_degree(coeffs2))

    if prec1 is None and prec2 is None and (max_degree < ring_prec):
        series = dup_sub(coeffs1, coeffs2, dom)
        return series, None

    coeffs1, coeffs2, prec = _unify_prec(s1, s2, dom, ring_prec)
    return dup_sub(coeffs1, coeffs2, dom), prec

def _useries_mul(
    s1: USeries[T], s2: USeries[T], dom: Domain, ring_prec: int
) -> USeries[T]:
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if prec1 is None and prec2 is None:
        d = dup_degree(coeffs1) + dup_degree(coeffs2)
        if d < ring_prec:
            return dup_mul(coeffs1, coeffs2, dom), None

    coeffs1, coeffs2, min_prec = _unify_prec(s1, s2, dom, ring_prec)
    series = dup_slice(dup_mul(coeffs1, coeffs2, dom), 0, min_prec, dom)
    return (series, min_prec)

def _useries_mul_ground(
    s: USeries[T], n: Any, dom: Domain, ring_prec: int
) -> USeries[T]:
    coeffs, prec = s

    if n == 0:
        if prec is not None:
            return [dom(0)], prec
        elif prec == 0:
            return [dom(0)], None
        else:
            return [], None
    if n == 1:
        return s
    if n == -1:
        neg_coeffs = dup_neg(coeffs, dom)
        return neg_coeffs, prec

    series = dup_mul_ground(coeffs, n, dom)
    deg = dup_degree(series)

    if prec is None:
        if deg < ring_prec:
            return series, None
        prec = ring_prec
    return dup_slice(series, 0, prec, dom), prec

def _useries_pow_int(s: USeries[T], n: int, dom: Domain, ring_prec: int) -> USeries[T]:
    """Raise a power series to an integer power."""
    if n < 0:
        raise ValueError("Power must be non-negative")

    coeffs, prec = s
    series = dup_pow(coeffs, n, dom)
    deg = dup_degree(series)
    if prec is None:
        if deg < ring_prec:
            return series, None
        prec = ring_prec
    series = dup_slice(series, 0, prec, dom)
    return series, prec

def _useries_truncate(s: USeries[T], n: int, dom: Domain) -> USeries[T]:
    """Truncate a power series to the first n terms."""
    coeffs, prec = s

    if n < 0:
        raise ValueError("Truncation precision must be non-negative")

    if prec is None:
        return dup_slice(coeffs, 0, n, dom), None

    prec = min(n, prec)
    return dup_slice(coeffs, 0, prec, dom), prec


class PythonPowerSeriesRingZZ(PowerSeriesRing[USeries[T]]):
    """Python implementation of power series ring over integers ring."""

    _domain = ZZ

    def __init__(self, prec: int = 6) -> None:
        self._prec = prec

    def __repr__(self) -> str:
        return f"Python Power Series Ring over {self.domain} with precision {self.prec}"

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, PythonPowerSeriesRingZZ):
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
    def one(self) -> USeries:
        if self.prec == 0:
            return ([], 0)
        return ([self.domain.one], None)

    @property
    def zero(self) -> USeries:
        if self.prec == 0:
            return ([], 0)
        return ([self.domain.zero], None)

    @property
    def gen(self) -> USeries:
        if self.prec == 0:
            return ([], 0)
        if self.prec == 1:
            return ([self.domain.zero], 1)
        return ([self.domain.one, self.domain.zero], None)

    def pretty(self, series: USeries[T]) -> str:
        coeffs, prec = series
        return _series_from_list(coeffs[::-1], prec)

    def print(self, series: USeries[T]) -> str:
        return self.pretty(series)

    def from_list(self, coeffs: DUP[T], prec: int | None = None) -> USeries[T]:
        """Create a power series from a list of coefficients."""
        if prec is None and len(coeffs) > self.prec:
            prec = self.prec
            coeffs = dup_slice(coeffs, 0, prec, self.domain)

        return dup_reverse(coeffs), prec

    def to_list(self, series: USeries[T]) -> DUP[T]:
        """Returns coeffs list."""
        coeffs, _ = series
        return coeffs[::-1]

    def equal(self, s1: USeries[T], s2: USeries[T]) -> bool | None:
        """Check if two power series are equal."""
        return _useries_equality(s1, s2)

    def negative(self, s: USeries[T]) -> USeries[T]:
        """Return the negative of a power series."""
        coeffs, prec = s
        neg_coeffs = dup_neg(coeffs, self.domain)
        return neg_coeffs, prec

    def add(self, s1: USeries[T], s2: USeries[T]) -> USeries[T]:
        """Add two power series."""
        return _useries_add(s1, s2, self.domain, self.prec)

    def subtract(self, s1: USeries[T], s2: USeries[T]) -> USeries[T]:
        """Subtract two power series."""
        return _useries_sub(s1, s2, self.domain, self.prec)

    def multiply(self, s1: USeries[T], s2: USeries[T]) -> USeries[T]:
        """Multiply two power series."""
        return _useries_mul(s1, s2, self.domain, self.prec)

    def multiply_ground(self, s: USeries[T], n: Any) -> USeries[T]:
        """Multiply a power series by a ground element."""
        return _useries_mul_ground(s, n, self.domain, self.prec)

    def pow_int(self, s: USeries[T], n: int) -> USeries[T]:
        """Raise a power series to an integer power."""
        return _useries_pow_int(s, n, self.domain, self.prec)

    def square(self, s: USeries[T]) -> USeries[T]:
        """Return the square of a power series."""
        return _useries_mul(s, s, self.domain, self.prec)

    def truncate(self, s: USeries[T], n: int) -> USeries[T]:
        """Truncate a power series to the first n terms."""
        return _useries_truncate(s, n, self.domain)


class PythonPowerSeriesRingQQ(PowerSeriesRing[USeries[T]]):
    """Python implementation of power series ring over rational field."""

    _domain = QQ

    def __init__(self, prec: int = 6) -> None:
        self._prec = prec

    def __repr__(self) -> str:
        return f"Python Power Series Ring over {self.domain} with precision {self.prec}"

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, PythonPowerSeriesRingQQ):
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
    def one(self) -> USeries[T]:
        if self.prec == 0:
            return ([], 0)
        return ([QQ(1)], None)

    @property
    def zero(self) -> USeries[T]:
        if self.prec == 0:
            return ([], 0)
        return ([QQ(0)], None)

    @property
    def gen(self) -> USeries[T]:
        if self.prec == 0:
            return ([], 0)
        if self.prec == 1:
            return ([QQ(0)], 1)
        return ([QQ(1), QQ(0)], None)

    def pretty(self, series: USeries[T]) -> str:
        coeffs, prec = series
        return _series_from_list(coeffs[::-1], prec)

    def print(self, series: USeries[T]) -> str:
        return self.pretty(series)

    def from_list(self, coeffs: DUP[T], prec: int | None = None) -> USeries[T]:
        """Create a power series from a list of coefficients."""
        if prec is None and len(coeffs) > self.prec:
            prec = self.prec
            coeffs = dup_slice(coeffs, 0, prec, self.domain)

        return dup_reverse(coeffs), prec

    def to_list(self, series: USeries[T]) -> DUP[T]:
        """Returns coeffs list."""
        coeffs, _ = series
        return coeffs[::-1]

    def equal(self, s1: USeries[T], s2: USeries[T]) -> bool | None:
        """Check if two power series are equal."""
        return _useries_equality(s1, s2)

    def negative(self, s: USeries[T]) -> USeries[T]:
        """Return the negative of a power series."""
        coeffs, prec = s
        neg_coeffs = dup_neg(coeffs, self.domain)
        return neg_coeffs, prec

    def add(self, s1: USeries[T], s2: USeries[T]) -> USeries[T]:
        """Add two power series."""
        return _useries_add(s1, s2, self.domain, self.prec)

    def subtract(self, s1: USeries[T], s2: USeries[T]) -> USeries[T]:
        """Subtract two power series."""
        return _useries_sub(s1, s2, self.domain, self.prec)

    def multiply(self, s1: USeries[T], s2: USeries[T]) -> USeries[T]:
        """Multiply two power series."""
        return _useries_mul(s1, s2, self.domain, self.prec)

    def multiply_ground(self, s: USeries[T], n: Any) -> USeries[T]:
        """Multiply a power series by a ground element."""
        return _useries_mul_ground(s, n, self.domain, self.prec)

    def pow_int(self, s: USeries[T], n: int) -> USeries[T]:
        """Raise a power series to an integer power."""
        return _useries_pow_int(s, n, self.domain, self.prec)

    def square(self, s: USeries[T]) -> USeries[T]:
        """Return the square of a power series."""
        return _useries_mul(s, s, self.domain, self.prec)

    def truncate(self, s: USeries[T], n: int) -> USeries[T]:
        """Truncate a power series to the first n terms."""
        return _useries_truncate(s, n, self.domain)
