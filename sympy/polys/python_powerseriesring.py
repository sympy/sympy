from __future__ import annotations

from typing import Any, TypeAlias

from sympy.polys.densearith import dup_add, dup_mul, dup_mul_ground, dup_sub
from sympy.polys.domains import QQ, ZZ
from sympy.polys.polyerrors import CoercionFailed
from sympy.polys.powerseriesring import (PowerSeriesRingQQ, PowerSeriesRingZZ,
                                        series_from_list)


DUPZ: TypeAlias = list[Any]
DUPQ: TypeAlias = list[Any]

ZZSeries: TypeAlias = tuple[DUPZ, int | None]
QQSeries: TypeAlias = tuple[DUPQ, int | None]


def _first_nonzero_index(coeffs):
    """Find the index of the first nonzero coefficient."""
    for i, coeff in enumerate(coeffs):
        if coeff != 0:
            return i
    return len(coeffs)


def _unify_prec(s1: ZZSeries | QQSeries,
                s2: ZZSeries | QQSeries
            ) -> tuple[DUPZ | DUPQ, DUPZ | DUPQ, int | None]:

    """Unify the precision of two series."""
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if prec1 is None:
        unified_prec = prec2
    elif prec2 is None:
        unified_prec = prec1
    else:
        unified_prec = min(prec1, prec2)

    if unified_prec is not None:
        coeffs1 = coeffs1[:unified_prec]
        coeffs2 = coeffs2[:unified_prec]

    return coeffs1, coeffs2, unified_prec


class PythonPowerSeriesRingZZ(PowerSeriesRingZZ):
    """Python implementation of power series ring over integers ring."""

    def __init__(self, prec: int = 6) -> None:
        self.domain = ZZ
        self.prec = prec

    def __repr__(self) -> str:
        return f"Python Power Series Ring over {self.domain} with precision {self.prec}"

    def __eq__(self, other: Any) -> bool:
        return (isinstance(other, PythonPowerSeriesRingZZ) and
               self.prec == other.prec)

    def __ne__(self, other: Any) -> bool:
        return not self.__eq__(other)

    def __hash__(self) -> int:
        return hash((self.domain, self.prec))

    @property
    def one(self) -> ZZSeries:
        return ([1], None)

    @property
    def zero(self) -> ZZSeries:
        return ([0], None)

    @property
    def gen(self) -> ZZSeries:
        return ([0, 1], None)

    def pretty(self, series: ZZSeries) -> str:
        coeffs, prec = series
        return series_from_list(coeffs, prec)

    def print(self, series: ZZSeries) -> str:
        return self.pretty(series)

    def equal(self, s1: ZZSeries, s2: ZZSeries) -> bool:
        """Check if two power series are equal."""

        coeffs1, prec1 = s1
        coeffs2, prec2 = s2

        return (prec1 == prec2 and
                coeffs1 == coeffs2)

    def negative(self, s: ZZSeries) -> ZZSeries:
        """Return the negative of a power series."""

        coeffs, prec = s

        neg_coeffs = [-c for c in coeffs]
        return neg_coeffs, prec

    def add(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Add two power series."""
        coeffs1, prec1 = s1
        coeffs2, prec2 = s2

        # case where both are DUPs
        if not prec1 and not prec2:
            coeffs1 = coeffs1[:self.prec]
            coeffs2 = coeffs2[:self.prec]

            series = dup_add(coeffs1[::-1], coeffs2[::-1], self.domain)
            series = series[::-1][:self.prec]
            if len(series) < self.prec:
                return series, None
            return series, self.prec

        coeffs1, coeffs2, prec = _unify_prec(s1, s2)
        series = dup_add(coeffs1[::-1], coeffs2[::-1], self.domain)
        series = series[::-1]

        return series, prec


    def subtract(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Subtract two power series."""
        coeffs1, prec1 = s1
        coeffs2, prec2 = s2

        # case where both are DUPs
        if not prec1 and not prec2:
            coeffs1 = coeffs1[:self.prec]
            coeffs2 = coeffs2[:self.prec]

            series = dup_sub(coeffs1[::-1], coeffs2[::-1], self.domain)
            series = series[::-1][:self.prec]
            if len(series) < self.prec:
                return series, None
            return series, self.prec

        coeffs1, coeffs2, prec = _unify_prec(s1, s2)
        series = dup_sub(coeffs1[::-1], coeffs2[::-1], self.domain)
        series = series[::-1]

        return series, prec

    def multiply(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Multiply two power series."""
        coeffs1, prec1 = s1
        coeffs2, prec2 = s2

        if not prec1 and not prec2:
            coeffs1 = coeffs1[:self.prec]
            coeffs2 = coeffs2[:self.prec]
            series = dup_mul(coeffs1[::-1], coeffs2[::-1], self.domain)
            series = series[::-1][:self.prec]
            if len(series) < self.prec:
                return series, None
            return series, self.prec

        if not prec1:
            min_prec = prec2
        elif not prec2:
            min_prec = prec1
        else:
            low_expv1 = _first_nonzero_index(coeffs1)
            low_expv2 = _first_nonzero_index(coeffs2)
            min_prec = min(prec1+low_expv2, prec2+low_expv1)

        coeffs1, coeffs2 = coeffs1[:min_prec], coeffs2[:min_prec]
        series = dup_mul(coeffs1[::-1], coeffs2[::-1], self.domain)
        series = series[::-1][:min_prec]

        return (series, min_prec)

    def multiply_ground(self, s: ZZSeries, n: Any) -> ZZSeries:
        """Multiply a power series by a ground element."""

        try:
            domain_ground = ZZ.convert(n)
            if n == domain_ground:
                n = domain_ground
            else:
                raise CoercionFailed(f"Cannot coerce {n} to {self.domain}")
        except Exception:
            raise CoercionFailed(f"Cannot coerce {n} to {self.domain}")

        coeffs, prec = s
        if n == 0:
            return ([0], None)
        if n == 1:
            return s
        if n == -1:
            return self.negative(s)

        series = dup_mul_ground(coeffs[::-1], n, self.domain)
        series = series[::-1]

        if prec is None:
            if len(series) < self.prec:
                return series, None
            prec = self.prec
        return series, prec

    def trunc(self, s: ZZSeries, n: int) -> ZZSeries:
        if n < 0:
            raise ValueError("Truncation precision must be non-negative")

        coeffs, prec = s
        if prec is None:
                return coeffs[:n], None

        return coeffs[:n], min(n, prec)


class PythonPowerSeriesRingQQ(PowerSeriesRingQQ):
    """Python implementation of power series ring over rational field."""

    def __init__(self, prec: int = 6) -> None:
        self.domain = QQ
        self.prec = prec

    def __repr__(self) -> str:
        return f"Python Power Series Ring over {self.domain} with precision {self.prec}"

    def __eq__(self, other: Any) -> bool:
        return (isinstance(other, PythonPowerSeriesRingQQ) and
                self.prec == other.prec)

    def __ne__(self, other: Any) -> bool:
        return not self.__eq__(other)

    def __hash__(self) -> int:
        return hash((self.domain, self.prec))

    @property
    def one(self) -> QQSeries:
        return ([QQ(1)], None)

    @property
    def zero(self) -> QQSeries:
        return ([QQ(0)], None)

    @property
    def gen(self) -> QQSeries:
        return ([QQ(0), QQ(1)], None)

    def pretty(self, series: QQSeries) -> str:
        coeffs, prec = series
        return series_from_list(coeffs, prec)

    def print(self, series: QQSeries) -> str:
        return self.pretty(series)

    def equal(self, s1: QQSeries, s2: QQSeries) -> bool:
        """Check if two power series are equal."""
        coeffs1, prec1 = s1
        coeffs2, prec2 = s2

        return (prec1 == prec2 and
                coeffs1 == coeffs2)

    def negative(self, s: QQSeries) -> QQSeries:
        """Return the negative of a power series."""

        coeffs, prec = s

        neg_coeffs = [-c for c in coeffs]
        return neg_coeffs, prec

    def add(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Add two power series."""
        coeffs1, prec1 = s1
        coeffs2, prec2 = s2

        # case where both are DUPs
        if not prec1 and not prec2:
            coeffs1 = coeffs1[:self.prec]
            coeffs2 = coeffs2[:self.prec]

            series = dup_add(coeffs1[::-1], coeffs2[::-1], self.domain)
            series = series[::-1][:self.prec]
            if len(series) < self.prec:
                return series, None
            return series, self.prec

        coeffs1, coeffs2, prec = _unify_prec(s1, s2)
        series = dup_add(coeffs1[::-1], coeffs2[::-1], self.domain)
        series = series[::-1]

        return series, prec

    def subtract(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Subtract two power series."""
        coeffs1, prec1 = s1
        coeffs2, prec2 = s2

        # case where both are DUPs
        if not prec1 and not prec2:
            coeffs1 = coeffs1[:self.prec]
            coeffs2 = coeffs2[:self.prec]

            series = dup_sub(coeffs1[::-1], coeffs2[::-1], self.domain)
            series = series[::-1][:self.prec]
            if len(series) < self.prec:
                return series, None
            return series, self.prec

        coeffs1, coeffs2, prec = _unify_prec(s1, s2)
        series = dup_sub(coeffs1[::-1], coeffs2[::-1], self.domain)
        series = series[::-1]

        return series, prec

    def multiply(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Multiply two power series."""
        coeffs1, prec1 = s1
        coeffs2, prec2 = s2

        if not prec1 and not prec2:
            coeffs1 = coeffs1[:self.prec]
            coeffs2 = coeffs2[:self.prec]
            series = dup_mul(coeffs1[::-1], coeffs2[::-1], self.domain)
            series = series[::-1][:self.prec]
            if len(series) < self.prec:
                return series, None
            return series, self.prec

        if not prec1:
            min_prec = prec2
        elif not prec2:
            min_prec = prec1
        else:
            low_expv1 = _first_nonzero_index(coeffs1)
            low_expv2 = _first_nonzero_index(coeffs2)
            min_prec = min(prec1 + low_expv2, prec2 + low_expv1)

        coeffs1, coeffs2 = coeffs1[:min_prec], coeffs2[:min_prec]
        series = dup_mul(coeffs1[::-1], coeffs2[::-1], self.domain)
        series = series[::-1][:min_prec]

        return (series, min_prec)

    def multiply_ground(self, s: QQSeries, n: Any) -> QQSeries:
        """Multiply a power series by a ground element."""

        try:
            domain_ground = QQ.convert(n)
            if n == domain_ground:
                n = domain_ground
            else:
                raise CoercionFailed(f"Cannot coerce {n} to {self.domain}")
        except Exception:
            raise CoercionFailed(f"Cannot coerce {n} to {self.domain}")

        coeffs, prec = s
        if n == 0:
            return ([QQ(0)], None)
        if n == 1:
            return s
        if n == -1:
            return self.negative(s)

        series = dup_mul_ground(coeffs[::-1], n, self.domain)
        series = series[::-1]

        if prec is None:
            if len(series) < self.prec:
                return series, None
            prec = self.prec
        return series, prec

    def trunc(self, s: QQSeries, n: int) -> QQSeries:
        if n < 0:
            raise ValueError("Truncation precision must be non-negative")

        coeffs, prec = s
        if prec is None:
            return coeffs[:n], None

        return coeffs[:n], min(n, prec)
