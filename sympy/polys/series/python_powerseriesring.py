from __future__ import annotations

from typing import Any, Union
from sympy.polys.densearith import (
    dup_add,
    dup_mul,
    dup_mul_ground,
    dup_neg,
    dup_sub,
    dup_series_mul,
    dup_series_pow,
)
from sympy.polys.densebasic import dup_degree, dup_reverse, dup_truncate, dup_from_list
from sympy.polys.densetools import (
    dup_diff,
    dup_integrate,
    dup_revert,
    dup_series_compose,
    dup_series_reversion,
)
from sympy.polys.polyerrors import NotReversible
from sympy.polys.domains import Domain, QQ, ZZ
from sympy.polys.domains.domain import Er
from sympy.polys.series.powerseriesring import _series_from_list, PowerSeriesRing
from sympy.external.gmpy import MPZ, MPQ


DUP = list[Er]
USeries = tuple[DUP[Er], Union[int, None]]


def _useries(
    coeffs: DUP[Er], series_prec: int | None, dom: Domain, ring_prec: int
) -> USeries[Er]:
    """Helper function to decide if the polynomial can be exact or it become a useries
    element."""

    if series_prec is None:
        deg = dup_degree(coeffs)
        if deg < ring_prec:
            return coeffs, None
        series_prec = ring_prec
    coeffs = dup_truncate(coeffs, series_prec, dom)
    return coeffs, series_prec


def _unify_prec(
    s1: USeries[Er], s2: USeries[Er], dom: Domain, ring_prec: int
) -> tuple[DUP[Er], DUP[Er], int]:
    """Unify the precision of two series."""
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if prec1 == prec2:
        unified_prec = prec1
    elif prec1 is None:
        unified_prec = prec2
    elif prec2 is None:
        unified_prec = prec1
    else:
        unified_prec = min(prec1, prec2)

    if unified_prec is None:
        unified_prec = ring_prec

    coeffs1 = dup_truncate(coeffs1, unified_prec, dom)
    coeffs2 = dup_truncate(coeffs2, unified_prec, dom)
    return coeffs1, coeffs2, unified_prec


def _useries_equality(s1: USeries[Er], s2: USeries[Er], dom: Domain) -> bool | None:
    """Check if two power series are equal."""
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if prec1 is None and prec2 is None:
        return s1 == s2

    coeffs1, coeffs2, _ = _unify_prec(s1, s2, dom, 0)
    if coeffs1 != coeffs2:
        return False
    return None


def _useries_equal_repr(s1: USeries[Er], s2: USeries[Er]) -> bool:
    return s1 == s2


def _useries_neg(s: USeries[Er], dom: Domain, ring_prec: int) -> USeries[Er]:
    coeffs, prec = s
    neg_coeffs = dup_neg(coeffs, dom)
    return _useries(neg_coeffs, prec, dom, ring_prec)


def _useries_add(
    s1: USeries[Er], s2: USeries[Er], dom: Domain, ring_prec: int
) -> USeries[Er]:
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2
    max_degree = max(dup_degree(coeffs1), dup_degree(coeffs2))

    if prec1 is None and prec2 is None and (max_degree < ring_prec):
        series = dup_add(coeffs1, coeffs2, dom)
        return series, None

    coeffs1, coeffs2, prec = _unify_prec(s1, s2, dom, ring_prec)
    return dup_add(coeffs1, coeffs2, dom), prec


def _useries_sub(
    s1: USeries[Er], s2: USeries[Er], dom: Domain, ring_prec: int
) -> USeries[Er]:
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2
    max_degree = max(dup_degree(coeffs1), dup_degree(coeffs2))

    if prec1 is None and prec2 is None and (max_degree < ring_prec):
        series = dup_sub(coeffs1, coeffs2, dom)
        return series, None

    coeffs1, coeffs2, prec = _unify_prec(s1, s2, dom, ring_prec)
    return dup_sub(coeffs1, coeffs2, dom), prec


def _useries_mul(
    s1: USeries[Er], s2: USeries[Er], dom: Domain, ring_prec: int
) -> USeries[Er]:
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if prec1 is None and prec2 is None:
        d = dup_degree(coeffs1) + dup_degree(coeffs2)
        if d < ring_prec:
            return dup_mul(coeffs1, coeffs2, dom), None

    coeffs1, coeffs2, min_prec = _unify_prec(s1, s2, dom, ring_prec)
    coeffs = dup_series_mul(coeffs1, coeffs2, min_prec, dom)
    return coeffs, min_prec


def _useries_mul_ground(
    s: USeries[Er], n: Any, dom: Domain, ring_prec: int
) -> USeries[Er]:
    coeffs, prec = s

    if n == 0:
        return [], prec
    if n == 1:
        return s
    if n == -1:
        return _useries_neg(s, dom, ring_prec)

    series = dup_mul_ground(coeffs, n, dom)
    return _useries(series, prec, dom, ring_prec)


def _useries_pow_int(
    s: USeries[Er], n: int, dom: Domain, ring_prec: int
) -> USeries[Er]:
    """Raise a power series to a non-negative integer power with truncation."""
    if n < 0:
        raise ValueError("Power must be a non-negative integer")

    coeffs, prec = s

    if n == 0:
        return [dom.one], prec

    if prec is None:
        deg = dup_degree(coeffs) * n
        if deg < ring_prec:
            return dup_series_pow(coeffs, n, ring_prec, dom), None
        prec = ring_prec

    series = dup_series_pow(coeffs, n, prec, dom)
    return series, prec


def _useries_truncate(s: USeries[Er], n: int, dom: Domain) -> USeries[Er]:
    """Truncate a power series to the first n terms."""
    coeffs, prec = s

    if n < 0:
        raise ValueError("Truncation precision must be non-negative")

    deg = dup_degree(coeffs)
    if deg < n:
        return coeffs, prec
    return dup_truncate(coeffs, n, dom), n


def _useries_compose(
    s1: USeries[Er], s2: USeries[Er], dom: Domain, ring_prec: int
) -> USeries[Er]:
    """Compose two power series."""
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if prec1 is None and prec2 is None:
        comp = dup_series_compose(coeffs1, coeffs2, ring_prec, dom)

        deg1 = dup_degree(coeffs1)
        deg2 = dup_degree(coeffs2)
        if deg1 + deg2 < ring_prec:
            return comp, None
        else:
            return comp, ring_prec

    if coeffs2 and not dom.is_zero(coeffs2[-1]):
        raise ValueError(
            "Series composition requires the constant term of the second series to be zero."
        )

    coeffs1, coeffs2, min_prec = _unify_prec(s1, s2, dom, ring_prec)
    comp = dup_series_compose(coeffs1, coeffs2, min_prec, dom)
    return comp, min_prec


def _useries_inversion(s: USeries[Er], dom: Domain, ring_prec: int) -> USeries[Er]:
    """Compute the series multiplicative inverse of a power series."""
    coeffs, prec = s

    if not coeffs or not dom.is_unit(coeffs[-1]):
        raise NotReversible("Series inversion requires the constant term to be a unit")

    if prec is None:
        prec = ring_prec

    inv = dup_revert(coeffs, prec, dom)
    inv = dup_truncate(inv, prec, dom)
    return inv, prec


def _useries_reversion(s: USeries[Er], dom: Domain, ring_prec: int) -> USeries[Er]:
    """Compute the composite inverse of a power series."""
    coeffs, prec = s

    if not coeffs or not dom.is_zero(coeffs[-1]):
        raise NotReversible("Series reversion requires the constant term to be zero.")

    if len(coeffs) >= 2 and not dom.is_unit(coeffs[-2]):
        raise NotReversible("Series reversion requires the linear term to be unit.")

    if prec is None:
        prec = ring_prec

    series = dup_series_reversion(coeffs, prec, dom)
    return series, prec


def _useries_derivative(s: USeries[Er], dom: Domain, ring_prec: int) -> USeries[Er]:
    """Compute the first derivative of a power series."""
    coeffs, prec = s
    series = dup_diff(coeffs, 1, dom)
    if prec:
        prec -= 1
    return _useries(series, prec, dom, ring_prec)


def _useries_integrate(s: USeries[Er], dom: Domain, ring_prec: int) -> USeries[Er]:
    """Compute the integral of a power series."""
    coeffs, prec = s
    series = dup_integrate(coeffs, 1, dom)
    if prec:
        prec += 1
    return _useries(series, prec, dom, ring_prec)


class PythonPowerSeriesRingZZ(PowerSeriesRing[USeries[MPZ]]):
    """Python implementation of power series ring over integers ring."""

    _domain = ZZ

    def __init__(self, prec: int = 6) -> None:
        if prec < 0:
            raise ValueError("Power series precision must be non-negative")
        self._prec = prec

    def __repr__(self) -> str:
        return (
            f"Python Power Series Ring over {self._domain} with precision {self._prec}"
        )

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, PythonPowerSeriesRingZZ):
            return NotImplemented
        return self._prec == other.prec

    def __hash__(self) -> int:
        return hash((self._domain, self._prec))

    @property
    def domain(self) -> Domain:
        """Return the ground domain of the power series ring."""
        return self._domain

    @property
    def prec(self) -> int:
        """Return the default precision for power series operations."""
        return self._prec

    @property
    def one(self) -> USeries[MPZ]:
        if self._prec == 0:
            return ([], 0)
        return ([self._domain.one], None)

    @property
    def zero(self) -> USeries[MPZ]:
        if self._prec == 0:
            return ([], 0)
        return ([], None)

    @property
    def gen(self) -> USeries[MPZ]:
        if self._prec < 2:
            return ([], self._prec)
        return ([self._domain.one, self._domain.zero], None)

    def pretty(self, s: USeries[MPZ]) -> str:
        coeffs, prec = s
        return _series_from_list(coeffs[::-1], prec)

    def print(self, s: USeries[MPZ]) -> None:
        print(self.pretty(s))

    def from_list(self, coeffs: list[MPZ], prec: int | None = None) -> USeries[MPZ]:
        """
        Create a power series from a list of coefficients in ascending order of
        expononets. If `prec` is not specified, it defaults to the ring's precision.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> s = R.from_list([1, 2, 3, 4, 5])
        >>> R.print(s)
        1 + 2*x + 3*x**2 + 4*x**3 + 5*x**4
        """
        coeffs = dup_reverse(coeffs)
        coeffs = dup_from_list(coeffs, self._domain)
        if prec is None and len(coeffs) > self._prec:
            prec = self._prec
            coeffs = dup_truncate(coeffs, prec, self._domain)

        return coeffs, prec

    def to_list(self, s: USeries[MPZ]) -> list[MPZ]:
        """
        Returns the list of coefficients.

        Examples
        ========
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> x = R.gen
        >>> R.to_list(x)
        [0, 1]
        """
        coeffs, _ = s
        return coeffs[::-1]

    def equal(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> bool | None:
        """Check if two power series are equal."""
        return _useries_equality(s1, s2, self._domain)

    def equal_repr(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> bool:
        """
        Check if two power series are equal coeffs and precision.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> s1 = R.from_list([1, 2, 1])
        >>> s2 = R.square(R.add(R.one, R.gen))
        >>> R.equal_repr(s1, s2)
        True
        """
        return _useries_equal_repr(s1, s2)

    def positive(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """
        Return the positive of a power series (which is the same as the series itself).

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> x = R.gen
        >>> R.print(R.positive(x))
        x
        """
        return _useries(s[0], s[1], self._domain, self._prec)

    def negative(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """
        Negate all the coeffs of power series.
        Examples
        =========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> x = R.gen
        >>> R.print(R.negative(x))
        -x
        """
        return _useries_neg(s, self._domain, self._prec)

    def add(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """Add two power series.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(3)
        >>> s1 = R.from_list([1, 2, 3, 4])
        >>> s2 = R.from_list([3, 4, 5, 6])
        >>> R.print(R.add(s1, s2))
        4 + 6*x + 8*x**2 + O(x**3)
        """
        return _useries_add(s1, s2, self._domain, self._prec)

    def subtract(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """
        Subtract two power series.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(3)
        >>> s1 = R.from_list([1, 2, 3, 4])
        >>> s2 = R.from_list([3, 4, 5, 6])
        >>> R.print(R.subtract(s1, s2))
        -2 - 2*x - 2*x**2 + O(x**3)
        """
        return _useries_sub(s1, s2, self._domain, self._prec)

    def multiply(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """
        Multiply two power series.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> s1 = R.from_list([1, 3, 3, 1])
        >>> s2 = R.from_list([1 , 4 , 6, 4, 1])
        >>> R.print(R.multiply(s1, s2))
        1 + 7*x + 21*x**2 + 35*x**3 + 35*x**4 + O(x**5)
        """
        return _useries_mul(s1, s2, self._domain, self._prec)

    def multiply_ground(self, s: USeries[MPZ], n: MPZ) -> USeries[MPZ]:
        """
        Multiply a power series by a ground element.

        Examples
        ========
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> x = R.gen
        >>> R.print(R.multiply_ground(x, 3))
        3*x
        """
        return _useries_mul_ground(s, n, self._domain, self._prec)

    def pow_int(self, s: USeries[MPZ], n: int) -> USeries[MPZ]:
        """
        Raise a power series to an integer power.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> s = R.from_list([1, 2, 1])
        >>> R.print(R.pow_int(s, 5))
        1 + 10*x + 45*x**2 + 120*x**3 + 210*x**4 + O(x**5)
        """
        return _useries_pow_int(s, n, self._domain, self._prec)

    def square(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """
        Return the square of a power series.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> s = R.from_list([1, 2, 1])
        >>> R.print(R.square(s))
        1 + 4*x + 6*x**2 + 4*x**3 + x**4
        """
        return _useries_mul(s, s, self._domain, self._prec)

    def compose(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """
        Compose two power series.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> s1 = R.from_list([1, 2])
        >>> s2 = R.from_list([0, 3])
        >>> R.print(R.compose(s1, s2))
        1 + 6*x
        """
        return _useries_compose(s1, s2, self._domain, self._prec)

    def inversion(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """
        Compute the series multiplicative inverse of a power series.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> s = R.from_list([1, 2, 3])
        >>> R.print(R.inversion(s))
        1 - 2*x + x**2 + 4*x**3 - 11*x**4 + O(x**5)
        """
        return _useries_inversion(s, self._domain, self._prec)

    def reversion(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """
        Compute the composite inverse of a power series.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> s = R.from_list([0, 1, 1])
        >>> R.print(R.reversion(s))
        x - x**2 + 2*x**3 - 5*x**4 + O(x**5)
        """
        return _useries_reversion(s, self._domain, self._prec)

    def truncate(self, s: USeries[MPZ], n: int) -> USeries[MPZ]:
        """
        Truncate a power series to the first n terms.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> s = R.from_list([1, 2, 3, 4, 5, 6])
        >>> R.print(s)
        1 + 2*x + 3*x**2 + 4*x**3 + 5*x**4 + O(x**5)
        >>> t = R.truncate(s, 3)
        >>> R.print(t)
        1 + 2*x + 3*x**2 + O(x**3)
        """
        return _useries_truncate(s, n, self._domain)

    def differentiate(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """
        Compute the derivative of a power series.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> s = R.from_list([1, 2, 3])
        >>> R.print(R.differentiate(s))
        2 + 6*x
        """
        return _useries_derivative(s, self._domain, self._prec)


class PythonPowerSeriesRingQQ(PowerSeriesRing[USeries[MPQ]]):
    """Python implementation of power series ring over rational field."""

    _domain = QQ

    def __init__(self, prec: int = 6) -> None:
        if prec < 0:
            raise ValueError("Power series precision must be non-negative")
        self._prec = prec

    def __repr__(self) -> str:
        return (
            f"Python Power Series Ring over {self._domain} with precision {self._prec}"
        )

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, PythonPowerSeriesRingQQ):
            return NotImplemented
        return self._prec == other.prec

    def __hash__(self) -> int:
        return hash((self._domain, self._prec))

    @property
    def domain(self) -> Domain:
        """Return the ground domain of the power series ring."""
        return self._domain

    @property
    def prec(self) -> int:
        """Return the default precision for power series operations."""
        return self._prec

    @property
    def one(self) -> USeries[MPQ]:
        if self._prec == 0:
            return ([], 0)
        return ([QQ(1)], None)

    @property
    def zero(self) -> USeries[MPQ]:
        if self._prec == 0:
            return ([], 0)
        return ([], None)

    @property
    def gen(self) -> USeries[MPQ]:
        if self._prec < 2:
            return ([], self._prec)
        return ([QQ(1), QQ(0)], None)

    def pretty(self, s: USeries[MPQ]) -> str:
        coeffs, prec = s
        return _series_from_list(coeffs[::-1], prec)

    def print(self, s: USeries[MPQ]) -> None:
        print(self.pretty(s))

    def from_list(self, coeffs: list[MPQ], prec: int | None = None) -> USeries[MPQ]:
        """
        Create a power series from a list of coefficients in ascending order of
        expononets. If `prec` is not specified, it defaults to the ring's precision.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s = R.from_list([QQ(1,2), QQ(3,4)])
        >>> R.print(s)
        1/2 + 3/4*x
        """
        coeffs = dup_reverse(coeffs)
        coeffs = dup_from_list(coeffs, self._domain)
        if prec is None and len(coeffs) > self._prec:
            prec = self._prec
            coeffs = dup_truncate(coeffs, prec, self._domain)

        return coeffs, prec

    def to_list(self, s: USeries[MPQ]) -> list[MPQ]:
        """
        Returns the list of coefficients.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> x = R.gen
        >>> R.to_list(x)
        [0, 1]
        """
        coeffs, _ = s
        return coeffs[::-1]

    def equal(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> bool | None:
        """Check if two power series are equal."""
        return _useries_equality(s1, s2, self._domain)

    def equal_repr(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> bool:
        """
        Check if two power series are equal coeffs and precision.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s1 = R.from_list([QQ(1), QQ(2), QQ(1)])
        >>> s2 = R.square(R.add(R.one, R.gen))
        >>> R.equal_repr(s1, s2)
        True
        """
        return _useries_equal_repr(s1, s2)

    def positive(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """
        Return the positive of a power series (which is the same as the series itself).

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> x = R.gen
        >>> R.print(R.positive(x))
        x
        """
        return _useries(s[0], s[1], self._domain, self._prec)

    def negative(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """
        Return the negative of a power series.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> x = R.gen
        >>> R.print(R.negative(x))
        -x
        """
        return _useries_neg(s, self._domain, self._prec)

    def add(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """Add two power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(3)
        >>> s1 = R.from_list([QQ(1,2), QQ(2,3)])
        >>> s2 = R.from_list([QQ(3,4), QQ(4,5)])
        >>> R.print(R.add(s1, s2))
        5/4 + 22/15*x
        """
        return _useries_add(s1, s2, self._domain, self._prec)

    def subtract(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """
        Subtract two power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(3)
        >>> s1 = R.from_list([QQ(1,2), QQ(2,3)])
        >>> s2 = R.from_list([QQ(3,4), QQ(4,5)])
        >>> R.print(R.subtract(s1, s2))
        -1/4 - 2/15*x
        """
        return _useries_sub(s1, s2, self._domain, self._prec)

    def multiply(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """
        Multiply two power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s1 = R.from_list([QQ(1,2), QQ(1,3)])
        >>> s2 = R.from_list([QQ(2), QQ(3)])
        >>> R.print(R.multiply(s1, s2))
        1 + 13/6*x + x**2
        """
        return _useries_mul(s1, s2, self._domain, self._prec)

    def multiply_ground(self, s: USeries[MPQ], n: MPQ) -> USeries[MPQ]:
        """
        Multiply a power series by a ground element.

        Examples
        ========
        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> x = R.gen
        >>> R.print(R.multiply_ground(x, QQ(3,2)))
        3/2*x
        """
        return _useries_mul_ground(s, n, self._domain, self._prec)

    def pow_int(self, s: USeries[MPQ], n: int) -> USeries[MPQ]:
        """
        Raise a power series to an integer power.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s = R.from_list([QQ(1,2), QQ(1,3)])
        >>> R.print(R.pow_int(s, 3))
        1/8 + 1/4*x + 1/6*x**2 + 1/27*x**3
        """
        return _useries_pow_int(s, n, self._domain, self._prec)

    def square(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """
        Return the square of a power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s = R.from_list([QQ(1,2), QQ(1,3)])
        >>> R.print(R.square(s))
        1/4 + 1/3*x + 1/9*x**2
        """
        return _useries_mul(s, s, self._domain, self._prec)

    def compose(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """
        Compose two power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s1 = R.from_list([QQ(1,2), QQ(2,3)])
        >>> s2 = R.from_list([QQ(0,1), QQ(3,4)])
        >>> R.print(R.compose(s1, s2))
        1/2 + 1/2*x
        """
        return _useries_compose(s1, s2, self._domain, self._prec)

    def inversion(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """
        Compute the series multiplicative inverse of a power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s = R.from_list([QQ(1,2), QQ(1,3)])
        >>> R.print(R.inversion(s))
        2 - 4/3*x + 8/9*x**2 - 16/27*x**3 + 32/81*x**4 + O(x**5)
        """
        return _useries_inversion(s, self._domain, self._prec)

    def reversion(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """
        Compute the composite inverse of a power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s = R.from_list([QQ(0), QQ(1,2)])
        >>> R.print(R.reversion(s))
        2*x + O(x**5)
        """
        return _useries_reversion(s, self._domain, self._prec)

    def truncate(self, s: USeries[MPQ], n: int) -> USeries[MPQ]:
        """
        Truncate a power series to the first n terms.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s = R.from_list([QQ(1,2), QQ(2,3), QQ(3,4), QQ(4,5), QQ(5,6)])
        >>> R.print(s)
        1/2 + 2/3*x + 3/4*x**2 + 4/5*x**3 + 5/6*x**4
        >>> t = R.truncate(s, 3)
        >>> R.print(t)
        1/2 + 2/3*x + 3/4*x**2 + O(x**3)
        """
        return _useries_truncate(s, n, self._domain)

    def differentiate(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """
        Compute the derivative of a power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s = R.from_list([QQ(1,2), QQ(2,3), QQ(3,4)])
        >>> R.print(R.differentiate(s))
        2/3 + 3/2*x
        """
        return _useries_derivative(s, self._domain, self._prec)

    def integrate(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """
        Compute the integral of a power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s = R.from_list([QQ(1,2), QQ(2,3), QQ(3,4)])
        >>> R.print(R.integrate(s))
        1/2*x + 1/3*x**2 + 1/4*x**3
        """
        return _useries_integrate(s, self._domain, self._prec)
