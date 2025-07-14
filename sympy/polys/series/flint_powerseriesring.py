from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Union

from sympy.polys.densebasic import dup_reverse
from sympy.polys.domains import Domain, QQ, ZZ
from sympy.polys.series.powerseriesring import (
    _series_from_list,
    PowerSeriesRing,
)
from sympy.external.gmpy import MPZ, MPQ
from sympy.polys.polyerrors import NotReversible
from sympy.polys.densetools import dup_revert, dup_series_compose
from sympy.polys.densetools import dup_truncate
from sympy.polys.series.python_powerseriesring import _unify_prec

from flint import fmpq_poly, fmpq_series, fmpz_poly, fmpz_series, ctx  # type: ignore


ZZSeries = Union[fmpz_series, fmpz_poly]
QQSeries = Union[fmpq_series, fmpq_poly]


def _get_series_precision(s: ZZSeries | QQSeries) -> int:
    """Helper function to get the precision of a series. By using the
    representation of the ring"""

    # XXX: This approach is inefficient, but as of python-flint 0.7.1, there is
    # no alternative method to extract the precision from a series element.
    rep = s.repr()
    prec = int(rep.split("prec=")[1].split(")")[0].strip())
    return prec


@contextmanager
def _global_cap(cap: int):
    """Temporarily set the global series cap within a context."""
    old_cap, ctx.cap = ctx.cap, cap
    try:
        yield
    finally:
        ctx.cap = old_cap


class FlintPowerSeriesRingZZ(PowerSeriesRing):
    """
    Flint implementation of power series ring over integers (ZZ).

    This class provides high-performance power series operations over the integer ring,
    leveraging the FLINT library for optimized arithmetic and series manipulations
    precision handling and truncation.

    Parameters
    ==========

    prec : int, optional
        The default precision for power series operations. Default is 6.

    Examples
    ========

    >>> from sympy.polys.series.flint_powerseriesring import FlintPowerSeriesRingZZ
    >>> R = FlintPowerSeriesRingZZ(5)
    >>> R
    Flint Power Series Ring over ZZ with precision 5

    >>> s = R([1, 2, 3])
    >>> R.print(s)
    1 + 2*x + 3*x**2

    >>> x = R.gen
    >>> squared = R.square(R.add(R.one, x))  # (1+x)^2
    >>> R.print(squared)
    1 + 2*x + x**2

    >>> # Geometric series: 1/(1-x)
    >>> geom = R.inverse(R.subtract(R.one, x))
    >>> R.print(geom)
    1 + x + x**2 + x**3 + x**4 + O(x**5)
    """

    _domain = ZZ

    def __init__(self, prec: int = 6) -> None:
        if prec < 0:
            raise ValueError("Power series precision must be non-negative")
        self._prec = prec

    def __repr__(self) -> str:
        return f"Flint Power Series Ring over {self.domain} with precision {self._prec}"

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, FlintPowerSeriesRingZZ):
            return NotImplemented
        return self._prec == other.prec

    def __hash__(self) -> int:
        return hash((self.domain, self._prec))

    def __call__(self, coeffs: list[MPZ], prec: int | None = None) -> ZZSeries:
        """Create a power series from a list of coefficients. If `prec` is not specified,
        it defaults to the ring's precision."""
        return self.from_list(coeffs, prec)

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
        if self._prec == 0:
            return fmpz_series([1], prec=0)
        return fmpz_poly([1])

    @property
    def zero(self) -> ZZSeries:
        if self._prec == 0:
            return fmpz_series([0], prec=0)
        return fmpz_poly([0])

    @property
    def gen(self) -> ZZSeries:
        if self._prec < 2:
            return fmpz_series([0, 1], prec=self._prec)
        return fmpz_poly([0, 1])

    def pretty(self, s: ZZSeries) -> str:
        """Pretty print a power series with improved formatting."""
        if isinstance(s, fmpz_poly):
            return _series_from_list(s.coeffs(), None)

        prec = _get_series_precision(s)
        return _series_from_list(s.coeffs(), prec)

    def print(self, s: ZZSeries) -> None:
        print(self.pretty(s))

    def from_list(self, coeffs: list[MPZ], prec: int | None = None) -> ZZSeries:
        """
        Create a power series from a list of coefficients in ascending order of
        exponents. If `prec` is not specified, it defaults to the ring's precision.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> s = R.from_list([1, 2, 3, 4, 5])
        >>> R.print(s)
        1 + 2*x + 3*x**2 + 4*x**3 + 5*x**4
        """
        if prec is None:
            if len(coeffs) <= self._prec:
                return fmpz_poly(coeffs)
            prec = self._prec

        return fmpz_series(coeffs, prec=prec)

    def to_list(self, s: ZZSeries) -> list[MPZ]:
        """
        Returns the list of coefficients.

        Examples
        ========
        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> x = R.gen
        >>> R.to_list(x)
        [0, 1]
        """
        return s.coeffs()

    def equal(self, s1: ZZSeries, s2: ZZSeries) -> bool | None:
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            return s1 == s2
        elif isinstance(s1, fmpz_poly):
            min_prec = _get_series_precision(s2)
        elif isinstance(s2, fmpz_poly):
            min_prec = _get_series_precision(s1)
        else:
            min_prec = min(_get_series_precision(s1), _get_series_precision(s2))

        coeffs1 = s1.coeffs()[:min_prec]
        coeffs2 = s2.coeffs()[:min_prec]

        if coeffs1 != coeffs2:
            return False
        return None

    def equal_repr(self, s1: ZZSeries, s2: ZZSeries) -> bool:
        """
        Check if two power series have the same representation.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> s1 = R.from_list([1, 2, 1])
        >>> s2 = R.square(R.add(R.one, R.gen))
        >>> R.equal_repr(s1, s2)
        True
        """
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            return s1 == s2
        elif isinstance(s1, fmpz_series) and isinstance(s2, fmpz_series):
            return s1._equal_repr(s2)
        else:
            return False

    def positive(self, s: ZZSeries) -> ZZSeries:
        """
        Return the positive of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> x = R.gen
        >>> R.print(R.positive(x))
        x
        """
        if isinstance(s, fmpz_poly):
            poly = s
            if poly.degree() < self._prec:
                return poly
            return fmpz_series(poly, prec=self._prec)
        return s

    def negative(self, s: ZZSeries) -> ZZSeries:
        """
        Return the negative of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> x = R.gen
        >>> R.print(R.negative(x))
        -x
        """
        return self.positive(-s)

    def add(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """
        Add two power series.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 3)
        >>> s1 = R.from_list([1, 2, 3, 4])
        >>> s2 = R.from_list([3, 4, 5, 6])
        >>> R.print(R.add(s1, s2))
        4 + 6*x + 8*x**2 + O(x**3)
        """
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            poly = s1 + s2
            if poly.degree() < self._prec:
                return poly
            return fmpz_series(poly, prec=self._prec)

        with _global_cap(self._prec):
            return s1 + s2

    def subtract(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """
        Subtract two power series.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 3)
        >>> s1 = R.from_list([1, 2, 3, 4])
        >>> s2 = R.from_list([3, 4, 5, 6])
        >>> R.print(R.subtract(s1, s2))
        -2 - 2*x - 2*x**2 + O(x**3)
        """
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            poly = s1 - s2
            if poly.degree() < self._prec:
                return poly
            return fmpz_series(poly, prec=self._prec)

        with _global_cap(self._prec):
            return s1 - s2

    def multiply(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """
        Multiply two power series.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> s1 = R.from_list([1, 3, 3, 1])
        >>> s2 = R.from_list([1, 4, 6, 4, 1])
        >>> R.print(R.multiply(s1, s2))
        1 + 7*x + 21*x**2 + 35*x**3 + 35*x**4 + O(x**5)
        """
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            poly = s1 * s2
            if poly.degree() < self._prec:
                return poly
            return fmpz_series(poly, prec=self._prec)

        with _global_cap(self._prec):
            return s1 * s2

    def multiply_ground(self, s: ZZSeries, n: MPZ) -> ZZSeries:
        """
        Multiply a power series by an integer or a polynomial.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> x = R.gen
        >>> R.print(R.multiply_ground(x, 3))
        3*x
        """
        if isinstance(s, fmpz_poly):
            poly = s * n
            if poly.degree() < self._prec:
                return poly
            return fmpz_series(poly, prec=self._prec)

        with _global_cap(self._prec):
            return s * n

    def pow_int(self, s: ZZSeries, n: int) -> ZZSeries:
        """
        Raise a power series to an integer power.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> s = R.from_list([1, 2, 1])
        >>> R.print(R.pow_int(s, 5))
        1 + 10*x + 45*x**2 + 120*x**3 + 210*x**4 + O(x**5)
        """
        if n < 0:
            raise ValueError("Power must be non-negative")

        if isinstance(s, fmpz_poly):
            poly = s**n
            if poly.degree() < self._prec:
                return poly
            return fmpz_series(poly, prec=self._prec)

        with _global_cap(self._prec):
            return s**n

    def square(self, s: ZZSeries) -> ZZSeries:
        """
        Return the square of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> s = R.from_list([1, 2, 1])
        >>> R.print(R.square(s))
        1 + 4*x + 6*x**2 + 4*x**3 + x**4
        """
        return self.pow_int(s, 2)

    def compose(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """
        Compose two power series.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> s1 = R.from_list([1, 2])
        >>> s2 = R.from_list([0, 3])
        >>> R.print(R.compose(s1, s2))
        1 + 6*x
        """
        dom = self._domain
        ring_prec = self._prec
        dup1 = dup_reverse(s1.coeffs())
        dup2 = dup_reverse(s2.coeffs())

        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            comp = dup_series_compose(dup1, dup2, ring_prec, dom)
            comp = comp[::-1]

            deg1 = s1.degree()
            deg2 = s2.degree()
            if deg1 + deg2 < ring_prec:
                return fmpz_poly(comp)
            else:
                return fmpz_series(comp, prec=ring_prec)

        if dup2 and not dom.is_zero(dup2[-1]):
            raise ValueError(
                "Series composition requires the second series to have a zero constant term."
            )

        prec1 = _get_series_precision(s1) if isinstance(s1, fmpz_series) else ring_prec
        prec2 = _get_series_precision(s2) if isinstance(s2, fmpz_series) else ring_prec

        DUP1 = (dup1, prec1)
        DUP2 = (dup2, prec2)

        coeffs1, coeffs2, min_prec = _unify_prec(DUP1, DUP2, dom, ring_prec)
        comp = dup_series_compose(coeffs1, coeffs2, min_prec, dom)
        comp = comp[::-1]
        return fmpz_series(comp, prec=min_prec)

    def inverse(self, s: ZZSeries) -> ZZSeries:
        """
        Compute the inverse of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> s = R.from_list([1, 2, 3])
        >>> R.print(R.inverse(s))
        1 - 2*x + x**2 + 4*x**3 - 11*x**4 + O(x**5)
        """
        dom = self._domain
        coeffs = s.coeffs()

        if not coeffs or not dom.is_unit(coeffs[0]):
            raise NotReversible(
                "Series inverse requires the constant term to be a unit"
            )

        prec = _get_series_precision(s) if isinstance(s, fmpz_series) else self._prec

        dup = dup_reverse(coeffs)
        inv = dup_revert(dup, prec, dom)
        inv = dup_truncate(inv, prec, dom)[::-1]
        return fmpz_series(inv, prec=prec)

    def reversion(self, s: ZZSeries) -> ZZSeries:
        """
        Compute the reversion of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> s = R.from_list([0, 1, 1])
        >>> R.print(R.reversion(s))
        x - x**2 + 2*x**3 - 5*x**4 + O(x**5)
        """
        dom = self._domain
        coeffs = s.coeffs()

        if not coeffs or not dom.is_zero(coeffs[0]):
            raise NotReversible(
                "Series reversion requires the constant term to be zero."
            )

        if len(coeffs) >= 2 and not dom.is_unit(coeffs[1]):
            raise NotReversible("Series reversion requires the linear term to be unit.")

        if isinstance(s, fmpz_poly):
            s = fmpz_series(coeffs, prec=self._prec)

        with _global_cap(self._prec):
            return s.reversion()

    def truncate(self, s: ZZSeries, n: int) -> ZZSeries:
        """
        Truncate a power series to the first n terms.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> s = R.from_list([1, 2, 3, 4, 5, 6])
        >>> R.print(s)
        1 + 2*x + 3*x**2 + 4*x**3 + 5*x**4 + O(x**5)
        >>> t = R.truncate(s, 3)
        >>> R.print(t)
        1 + 2*x + 3*x**2 + O(x**3)
        """
        if n < 0:
            raise ValueError("Truncation precision must be non-negative")

        if len(s) <= n:
            return s

        coeffs = s.coeffs()[:n]
        return fmpz_series(coeffs, prec=n)

    def differentiate(self, s: ZZSeries) -> ZZSeries:
        """
        Compute the derivative of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import ZZ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(ZZ, 5)
        >>> s = R.from_list([1, 2, 1])
        >>> R.print(R.differentiate(s))
        2 + 2*x
        """
        if isinstance(s, fmpz_poly):
            poly = s.derivative()
            if poly.degree() < self._prec:
                return poly
            return fmpz_series(poly, prec=self._prec)

        poly = fmpz_poly(s.coeffs())
        derivative = poly.derivative()
        prec = _get_series_precision(s)
        return fmpz_series(derivative, prec=prec - 1)


class FlintPowerSeriesRingQQ(PowerSeriesRing):
    """
    Flint implementation of power series ring over rational field (QQ).

    This class provides high-performance power series operations over the rational field,
    leveraging the FLINT library for optimized arithmetic and series manipulations
    with  precision handling and truncation. It extends the integer ring functionality
    with support for rational coefficients and integration.

    Parameters
    ==========

    prec : int, optional
        The default precision for power series operations. Default is 6.

    Examples
    ========

    >>> from sympy import QQ
    >>> from sympy.polys.series.flint_powerseriesring import FlintPowerSeriesRingQQ
    >>> R = FlintPowerSeriesRingQQ(5)
    >>> R
    Flint Power Series Ring over QQ with precision 5

    >>> s = R([QQ(1,2), QQ(2,3)])
    >>> R.print(s)
    1/2 + 2/3*x

    >>> x = R.gen
    >>> integrated = R.integrate(R.add(R.one, x))
    >>> R.print(integrated)
    x + 1/2*x**2

    >>> inv = R.inverse(s)
    >>> R.print(inv)
    2 - 8/3*x + 32/9*x**2 - 128/27*x**3 + 512/81*x**4 + O(x**5)
    """

    _domain = QQ

    def __init__(self, prec: int = 6) -> None:
        if prec < 0:
            raise ValueError("Power series precision must be non-negative")
        self._prec = prec

    def __repr__(self) -> str:
        return f"Flint Power Series Ring over {self.domain} with precision {self._prec}"

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, FlintPowerSeriesRingQQ):
            return NotImplemented
        return self._prec == other.prec

    def __hash__(self) -> int:
        return hash((self.domain, self._prec))

    def __call__(self, coeffs: list[MPQ], prec: int | None = None) -> QQSeries:
        """Create a power series from a list of coefficients. If `prec` is not specified,
        it defaults to the ring's precision."""
        return self.from_list(coeffs, prec)

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
        if self._prec == 0:
            return fmpq_series([1], prec=0)
        return fmpq_poly([1])

    @property
    def zero(self) -> QQSeries:
        if self._prec == 0:
            return fmpq_series([0], prec=0)
        return fmpq_poly([0])

    @property
    def gen(self) -> QQSeries:
        if self._prec < 2:
            return fmpq_series([0, 1], prec=self._prec)
        return fmpq_poly([0, 1])

    def pretty(self, s: QQSeries) -> str:
        """Pretty print a power series with improved formatting."""
        if isinstance(s, fmpq_poly):
            return _series_from_list(s.coeffs(), None)

        prec = _get_series_precision(s)
        return _series_from_list(s.coeffs(), prec)

    def print(self, s: QQSeries) -> None:
        print(self.pretty(s))

    def from_list(self, coeffs: list[MPQ], prec: int | None = None) -> QQSeries:
        """
        Create a power series from a list of coefficients in ascending order of
        exponents. If `prec` is not specified, it defaults to the ring's precision.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> s = R.from_list([QQ(1,2), QQ(3,4)])
        >>> R.print(s)
        1/2 + 3/4*x
        """
        if prec is None:
            if len(coeffs) <= self._prec:
                return fmpq_poly(coeffs)
            prec = self._prec

        return fmpq_series(coeffs, prec=prec)

    def to_list(self, s: QQSeries) -> list[MPQ]:
        """
        Returns the list of coefficients.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> x = R.gen
        >>> R.to_list(x)
        [0, 1]
        """
        return s.coeffs()

    def equal(self, s1: QQSeries, s2: QQSeries) -> bool | None:
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            return s1 == s2
        elif isinstance(s1, fmpq_poly):
            min_prec = _get_series_precision(s2)
        elif isinstance(s2, fmpq_poly):
            min_prec = _get_series_precision(s1)
        else:
            min_prec = min(_get_series_precision(s1), _get_series_precision(s2))

        coeffs1 = s1.coeffs()[:min_prec]
        coeffs2 = s2.coeffs()[:min_prec]

        if coeffs1 != coeffs2:
            return False
        return None

    def equal_repr(self, s1: QQSeries, s2: QQSeries) -> bool:
        """
        Check if two power series have the same representation.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> s1 = R.from_list([QQ(1), QQ(2), QQ(1)])
        >>> s2 = R.square(R.add(R.one, R.gen))
        >>> R.equal_repr(s1, s2)
        True
        """
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            return s1 == s2
        elif isinstance(s1, fmpq_series) and isinstance(s2, fmpq_series):
            return s1._equal_repr(s2)
        else:
            return False

    def positive(self, s: QQSeries) -> QQSeries:
        """
        Return the positive of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> x = R.gen
        >>> R.print(R.positive(x))
        x
        """
        if isinstance(s, fmpq_poly):
            poly = s
            if poly.degree() < self._prec:
                return poly
            return fmpq_series(poly, prec=self._prec)
        return s

    def negative(self, s: QQSeries) -> QQSeries:
        """
        Return the negative of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> x = R.gen
        >>> R.print(R.negative(x))
        -x
        """
        return self.positive(-s)

    def add(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """
        Add two power series.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 3)
        >>> s1 = R.from_list([QQ(1,2), QQ(2,3)])
        >>> s2 = R.from_list([QQ(3,4), QQ(4,5)])
        >>> R.print(R.add(s1, s2))
        5/4 + 22/15*x
        """
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            poly = s1 + s2
            if poly.degree() < self._prec:
                return poly
            return fmpq_series(poly, prec=self._prec)

        with _global_cap(self._prec):
            return s1 + s2

    def subtract(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """
        Subtract two power series.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 3)
        >>> s1 = R.from_list([QQ(1,2), QQ(2,3)])
        >>> s2 = R.from_list([QQ(3,4), QQ(4,5)])
        >>> R.print(R.subtract(s1, s2))
        -1/4 - 2/15*x
        """
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            poly = s1 - s2
            if poly.degree() < self._prec:
                return poly
            return fmpq_series(poly, prec=self._prec)

        with _global_cap(self._prec):
            return s1 - s2

    def multiply(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """
        Multiply two power series.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> s1 = R.from_list([QQ(1,2), QQ(1,3)])
        >>> s2 = R.from_list([QQ(2), QQ(3)])
        >>> R.print(R.multiply(s1, s2))
        1 + 13/6*x + x**2
        """
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            poly = s1 * s2
            if poly.degree() < self._prec:
                return poly
            return fmpq_series(poly, prec=self._prec)

        with _global_cap(self._prec):
            return s1 * s2

    def multiply_ground(self, s: QQSeries, n: MPQ) -> QQSeries:
        """
        Multiply a power series by a rational number.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> x = R.gen
        >>> R.print(R.multiply_ground(x, QQ(3,2)))
        3/2*x
        """
        if isinstance(s, fmpq_poly):
            poly = s * n
            if poly.degree() < self._prec:
                return poly
            return fmpq_series(poly, prec=self._prec)

        with _global_cap(self._prec):
            return s * n

    def pow_int(self, s: QQSeries, n: int) -> QQSeries:
        """
        Raise a power series to an integer power.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> s = R.from_list([QQ(1,2), QQ(1,3)])
        >>> R.print(R.pow_int(s, 3))
        1/8 + 1/4*x + 1/6*x**2 + 1/27*x**3
        """
        if n < 0:
            raise ValueError("Power must be non-negative")

        if isinstance(s, fmpq_poly):
            poly = s**n
            if poly.degree() < self._prec:
                return poly
            return fmpq_series(poly, prec=self._prec)

        with _global_cap(self._prec):
            return s**n

    def square(self, s: QQSeries) -> QQSeries:
        """
        Return the square of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> s = R.from_list([QQ(1,2), QQ(1,3)])
        >>> R.print(R.square(s))
        1/4 + 1/3*x + 1/9*x**2
        """
        return self.pow_int(s, 2)

    def truncate(self, s: QQSeries, n: int) -> QQSeries:
        """
        Truncate a power series to the first n terms.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> s = R.from_list([QQ(1,2), QQ(2,3), QQ(3,4), QQ(4,5), QQ(5,6), QQ(1)])
        >>> R.print(s)
        1/2 + 2/3*x + 3/4*x**2 + 4/5*x**3 + 5/6*x**4 + O(x**5)
        >>> t = R.truncate(s, 3)
        >>> R.print(t)
        1/2 + 2/3*x + 3/4*x**2 + O(x**3)
        """
        if n < 0:
            raise ValueError("Truncation precision must be non-negative")

        if len(s) <= n:
            return s

        coeffs = s.coeffs()[:n]
        return fmpq_series(coeffs, prec=n)

    def compose(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """
        Compose two power series.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> s1 = R.from_list([QQ(1,2), QQ(2,3)])
        >>> s2 = R.from_list([QQ(0,1), QQ(3,4)])
        >>> R.print(R.compose(s1, s2))
        1/2 + 1/2*x
        """
        dom = self._domain
        ring_prec = self._prec
        dup1 = dup_reverse(s1.coeffs())
        dup2 = dup_reverse(s2.coeffs())

        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            comp = dup_series_compose(dup1, dup2, ring_prec, dom)
            comp = comp[::-1]

            deg1 = s1.degree()
            deg2 = s2.degree()
            if deg1 + deg2 < ring_prec:
                return fmpq_poly(comp)
            else:
                return fmpq_series(comp, prec=ring_prec)

        if dup2 and not dom.is_zero(dup2[-1]):
            raise ValueError(
                "Series composition requires the second series to have a zero constant term."
            )

        prec1 = _get_series_precision(s1) if isinstance(s1, fmpq_series) else ring_prec
        prec2 = _get_series_precision(s2) if isinstance(s2, fmpq_series) else ring_prec

        DUP1 = (dup1, prec1)
        DUP2 = (dup2, prec2)

        coeffs1, coeffs2, min_prec = _unify_prec(DUP1, DUP2, dom, ring_prec)
        comp = dup_series_compose(coeffs1, coeffs2, min_prec, dom)
        comp = comp[::-1]
        return fmpq_series(comp, prec=min_prec)

    def inverse(self, s: QQSeries) -> QQSeries:
        """
        Compute the inverse of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> s = R.from_list([QQ(1,2), QQ(1,3)])
        >>> R.print(R.inverse(s))
        2 - 4/3*x + 8/9*x**2 - 16/27*x**3 + 32/81*x**4 + O(x**5)
        """
        dom = self._domain
        ring_prec = self._prec
        coeffs = s.coeffs()

        if not coeffs or not dom.is_unit(coeffs[0]):
            raise NotReversible(
                "Series inverse requires the constant term to be a unit"
            )

        prec = _get_series_precision(s) if isinstance(s, fmpq_series) else ring_prec

        dup = dup_reverse(coeffs)
        inv = dup_revert(dup, prec, dom)
        inv = dup_truncate(inv, prec, dom)[::-1]
        return fmpq_series(inv, prec=prec)

    def reversion(self, s: QQSeries) -> QQSeries:
        """
        Compute the reversion of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> s = R.from_list([QQ(0,1), QQ(1,2)])
        >>> R.print(R.reversion(s))
        2*x + O(x**5)
        """
        dom = self._domain
        coeffs = s.coeffs()

        if not coeffs or not dom.is_zero(coeffs[0]):
            raise NotReversible(
                "Series reversion requires the constant term to be zero."
            )

        if len(coeffs) >= 2 and not dom.is_unit(coeffs[1]):
            raise NotReversible("Series reversion requires the linear term to be unit.")

        if isinstance(s, fmpq_poly):
            s = fmpq_series(coeffs, prec=self._prec)

        with _global_cap(self._prec):
            return s.reversion()

    def differentiate(self, s: QQSeries) -> QQSeries:
        """
        Compute the derivative of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> s = R.from_list([QQ(1,2), QQ(1,3)])
        >>> R.print(R.differentiate(s))
        1/3
        """
        if isinstance(s, fmpq_poly):
            poly = s.derivative()
            if poly.degree() < self._prec:
                return poly
            return fmpq_series(poly, prec=self._prec)

        poly = fmpq_poly(s.coeffs())
        derivative = poly.derivative()
        prec = _get_series_precision(s)
        return fmpq_series(derivative, prec=prec - 1)

    def integrate(self, s: QQSeries) -> QQSeries:
        """
        Compute the integral of a power series.

        Examples
        ========

        >>> from sympy.polys.domains import QQ
        >>> from sympy.polys.series import power_series_ring
        >>> R = power_series_ring(QQ, 5)
        >>> s = R.from_list([QQ(1,2), QQ(1,3)])
        >>> R.print(R.integrate(s))
        1/2*x + 1/6*x**2
        """
        if isinstance(s, fmpq_poly):
            poly = s.integral()
            if poly.degree() < self._prec:
                return poly
            return fmpq_series(poly, prec=self._prec)
        with _global_cap(self._prec):
            return s.integral()
