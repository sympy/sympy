from __future__ import annotations

from typing import Any, TypeVar
from sympy.polys.densearith import dup_add, dup_mul, dup_mul_ground, dup_neg, dup_sub
from sympy.polys.densebasic import dup_degree, dup_reverse, dup_slice
from sympy.polys.densetools import dup_diff, dup_integrate
from sympy.polys.domains import Domain, QQ, ZZ
from sympy.polys.series.powerseriesring import _series_from_list, PowerSeriesRing


T = TypeVar("T")
DUP = list[T]
USeries = tuple[DUP[T], int | None]


def _useries(
    coeffs: DUP[T], series_prec: int | None, dom: Domain, ring_prec: int
) -> USeries[T]:
    """Helper function to decide if the polynomial can be exact or it become a useries
    element."""

    if series_prec is None:
        deg = dup_degree(coeffs)
        if deg < ring_prec:
            return coeffs, None
        series_prec = ring_prec
    coeffs = dup_slice(coeffs, 0, series_prec, dom)
    return coeffs, series_prec


def _unify_prec(
    s1: USeries[T], s2: USeries[T], dom: Domain, ring_prec: int
) -> tuple[DUP[T], DUP[T], int]:
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

    coeffs1 = dup_slice(coeffs1, 0, unified_prec, dom)
    coeffs2 = dup_slice(coeffs2, 0, unified_prec, dom)
    return coeffs1, coeffs2, unified_prec


def _useries_equality(s1: USeries[T], s2: USeries[T]) -> bool | None:
    """Check if two power series are equal."""
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if prec1 is None and prec2 is None:
        return s1 == s2

    if coeffs1 != coeffs2:
        return False
    return None


def _useries_equal_repr(s1: USeries[T], s2: USeries[T]) -> bool:
    return s1 == s2


def _useries_neg(s: USeries[T], dom: Domain) -> USeries[T]:
    coeffs, prec = s
    neg_coeffs = dup_neg(coeffs, dom)
    return neg_coeffs, prec


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
    coeffs = dup_mul(coeffs1, coeffs2, dom)
    return _useries(coeffs, min_prec, dom, ring_prec)


def _useries_mul_ground(
    s: USeries[T], n: Any, dom: Domain, ring_prec: int
) -> USeries[T]:
    coeffs, prec = s

    if n == 0:
        return [], prec
    if n == 1:
        return s
    if n == -1:
        return _useries_neg(s, dom)

    series = dup_mul_ground(coeffs, n, dom)
    return _useries(series, prec, dom, ring_prec)


def _pow_recursive(coeffs: DUP[T], n: int, prec: int, dom: Domain) -> list[T]:
    """Raise a truncated power series `coeffs` to the power `n` modulo x^prec."""
    if n == 1:
        return dup_slice(coeffs, 0, prec, dom)

    q, r = divmod(n, 2)
    half = _pow_recursive(coeffs, q, prec, dom)
    square = _useries_mul((half, prec), (half, prec), dom, prec)[0]

    if r == 0:
        return square
    return _useries_mul((square, prec), (coeffs, prec), dom, prec)[0]


def _useries_pow_int(s: USeries[T], n: int, dom: Domain, ring_prec: int) -> USeries[T]:
    """Raise a power series to a non-negative integer power with truncation."""
    if n < 0:
        raise ValueError("Power must be a non-negative integer")

    coeffs, prec = s

    if n == 0:
        return [dom.one], prec

    if prec is None:
        deg = dup_degree(coeffs) * n
        if deg < ring_prec:
            return _pow_recursive(coeffs, n, ring_prec, dom), None
        prec = ring_prec

    result = _pow_recursive(coeffs, n, prec, dom)
    return result, prec


def _useries_truncate(s: USeries[T], n: int, dom: Domain) -> USeries[T]:
    """Truncate a power series to the first n terms."""
    coeffs, prec = s

    if n < 0:
        raise ValueError("Truncation precision must be non-negative")

    deg = dup_degree(coeffs)
    if deg < n:
        return coeffs, prec
    return dup_slice(coeffs, 0, n, dom), n


def _useries_derivative(s: USeries[T], dom: Domain, ring_prec: int) -> USeries[T]:
    """Compute the first derivative of a power series."""
    coeffs, prec = s
    series = dup_diff(coeffs, 1, dom)
    return _useries(series, prec, dom, ring_prec)


def _useries_integrate(s: USeries[T], dom: Domain, ring_prec: int) -> USeries[T]:
    """Compute the integral of a power series."""
    coeffs, prec = s
    series = dup_integrate(coeffs, 1, dom)
    return _useries(series, prec, dom, ring_prec)


class PythonPowerSeriesRingZZ(PowerSeriesRing[USeries[T]]):
    """Python implementation of power series ring over integers ring."""

    _domain = ZZ

    def __init__(self, prec: int = 6) -> None:
        if prec < 0:
            raise ValueError("Power series precision must be non-negative")
        self._prec = prec

    def __repr__(self) -> str:
        return (
            f"Python Power Series Ring over {self.domain} with precision {self._prec}"
        )

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, PythonPowerSeriesRingZZ):
            return NotImplemented
        return self._prec == other.prec

    def __hash__(self) -> int:
        return hash((self.domain, self._prec))

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
        if self._prec == 0:
            return ([], 0)
        return ([self.domain.one], None)

    @property
    def zero(self) -> USeries:
        if self._prec == 0:
            return ([], 0)
        return ([], None)

    @property
    def gen(self) -> USeries:
        if self._prec < 2:
            return ([], self._prec)
        return ([self.domain.one, self.domain.zero], None)

    def pretty(self, series: USeries[T]) -> str:
        coeffs, prec = series
        return _series_from_list(coeffs[::-1], prec)

    def print(self, series: USeries[T]) -> None:
        print(self.pretty(series))

    def from_list(self, coeffs: list[T], prec: int | None = None) -> USeries[T]:
        """
        Create a power series from a list of coefficients in ascending order of
        expononets. If `prec` is not specified, it defaults to the ring's precision.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> s = R.from_list([1, 2, 3, 4, 5])
        >>> R.print(s)
        1 + 2*x + 3*x**2 + 4*x**3 + 5*x**4 + O(x**5)
        """
        coeffs = dup_reverse(coeffs)
        if prec is None and len(coeffs) > self._prec:
            prec = self._prec
            coeffs = dup_slice(coeffs, 0, prec, self.domain)

        return coeffs, prec

    def to_list(self, series: USeries[T]) -> list[T]:
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
        coeffs, _ = series
        return coeffs[::-1]

    def equal(self, s1: USeries[T], s2: USeries[T]) -> bool | None:
        """Check if two power series are equal."""
        return _useries_equality(s1, s2)

    def equal_repr(self, s1: USeries[T], s2: USeries[T]) -> bool:
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

    def negative(self, s: USeries[T]) -> USeries[T]:
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
        return _useries_neg(s, self.domain)

    def add(self, s1: USeries[T], s2: USeries[T]) -> USeries[T]:
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
        return _useries_add(s1, s2, self.domain, self._prec)

    def subtract(self, s1: USeries[T], s2: USeries[T]) -> USeries[T]:
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
        return _useries_sub(s1, s2, self.domain, self._prec)

    def multiply(self, s1: USeries[T], s2: USeries[T]) -> USeries[T]:
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
        return _useries_mul(s1, s2, self.domain, self._prec)

    def multiply_ground(self, s: USeries[T], n: Any) -> USeries[T]:
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
        return _useries_mul_ground(s, n, self.domain, self._prec)

    def pow_int(self, s: USeries[T], n: int) -> USeries[T]:
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
        return _useries_pow_int(s, n, self.domain, self._prec)

    def square(self, s: USeries[T]) -> USeries[T]:
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
        return _useries_mul(s, s, self.domain, self._prec)

    def truncate(self, s: USeries[T], n: int) -> USeries[T]:
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
        return _useries_truncate(s, n, self.domain)

    def differentiate(self, s: USeries[T]) -> USeries[T]:
        """
        Compute the derivative of a power series.

        Examples
        ========

        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingZZ
        >>> R = PythonPowerSeriesRingZZ(5)
        >>> s = R.from_list([1, 2, 3])
        >>> R.print(R.differentiate(s))
        2 + 6*x + O(x**2)
        """
        return _useries_derivative(s, self.domain, self._prec)


class PythonPowerSeriesRingQQ(PowerSeriesRing[USeries[T]]):
    """Python implementation of power series ring over rational field."""

    _domain = QQ

    def __init__(self, prec: int = 6) -> None:
        if prec < 0:
            raise ValueError("Power series precision must be non-negative")
        self._prec = prec

    def __repr__(self) -> str:
        return (
            f"Python Power Series Ring over {self.domain} with precision {self._prec}"
        )

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, PythonPowerSeriesRingQQ):
            return NotImplemented
        return self._prec == other.prec

    def __hash__(self) -> int:
        return hash((self.domain, self._prec))

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
        if self._prec == 0:
            return ([], 0)
        return ([QQ(1)], None)

    @property
    def zero(self) -> USeries[T]:
        if self._prec == 0:
            return ([], 0)
        return ([], None)

    @property
    def gen(self) -> USeries[T]:
        if self._prec < 2:
            return ([], self._prec)
        return ([QQ(1), QQ(0)], None)

    def pretty(self, series: USeries[T]) -> str:
        coeffs, prec = series
        return _series_from_list(coeffs[::-1], prec)

    def print(self, series: USeries[T]) -> None:
        print(self.pretty(series))

    def from_list(self, coeffs: list[T], prec: int | None = None) -> USeries[T]:
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
        if prec is None and len(coeffs) > self._prec:
            prec = self._prec
            coeffs = dup_slice(coeffs, 0, prec, self.domain)

        return coeffs, prec

    def to_list(self, series: USeries[T]) -> list[T]:
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
        coeffs, _ = series
        return coeffs[::-1]

    def equal(self, s1: USeries[T], s2: USeries[T]) -> bool | None:
        """Check if two power series are equal."""
        return _useries_equality(s1, s2)

    def equal_repr(self, s1: USeries[T], s2: USeries[T]) -> bool:
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

    def negative(self, s: USeries[T]) -> USeries[T]:
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
        return _useries_neg(s, self.domain)

    def add(self, s1: USeries[T], s2: USeries[T]) -> USeries[T]:
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
        return _useries_add(s1, s2, self.domain, self._prec)

    def subtract(self, s1: USeries[T], s2: USeries[T]) -> USeries[T]:
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
        return _useries_sub(s1, s2, self.domain, self._prec)

    def multiply(self, s1: USeries[T], s2: USeries[T]) -> USeries[T]:
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
        return _useries_mul(s1, s2, self.domain, self._prec)

    def multiply_ground(self, s: USeries[T], n: Any) -> USeries[T]:
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
        return _useries_mul_ground(s, n, self.domain, self._prec)

    def pow_int(self, s: USeries[T], n: int) -> USeries[T]:
        """
        Raise a power series to an integer power.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s = R.from_list([QQ(1,2), QQ(1,3)])
        >>> R.print(R.pow_int(s, 3))
        1/8 + 1/4*x + 19/72*x**2 + 1/18*x**3 + 5/324*x**4 + O(x**5)
        """
        return _useries_pow_int(s, n, self.domain, self._prec)

    def square(self, s: USeries[T]) -> USeries[T]:
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
        return _useries_mul(s, s, self.domain, self._prec)

    def truncate(self, s: USeries[T], n: int) -> USeries[T]:
        """
        Truncate a power series to the first n terms.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s = R.from_list([QQ(1,2), QQ(2,3), QQ(3,4), QQ(4,5), QQ(5,6)])
        >>> R.print(s)
        1/2 + 2/3*x + 3/4*x**2 + 4/5*x**3 + 5/6*x**4 + O(x**5)
        >>> t = R.truncate(s, 3)
        >>> R.print(t)
        1/2 + 2/3*x + 3/4*x**2 + O(x**3)
        """
        return _useries_truncate(s, n, self.domain)

    def differentiate(self, s: USeries[T]) -> USeries[T]:
        """
        Compute the derivative of a power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s = R.from_list([QQ(1,2), QQ(2,3), QQ(3,4)])
        >>> R.print(R.differentiate(s))
        2/3 + 3/2*x + O(x**2)
        """
        return _useries_derivative(s, self.domain, self._prec)

    def integrate(self, s: USeries[T]) -> USeries[T]:
        """
        Compute the integral of a power series.

        Examples
        ========

        >>> from sympy import QQ
        >>> from sympy.polys.series.python_powerseriesring import PythonPowerSeriesRingQQ
        >>> R = PythonPowerSeriesRingQQ(5)
        >>> s = R.from_list([QQ(1,2), QQ(2,3), QQ(3,4)])
        >>> R.print(R.integrate(s))
        1/2*x + 1/6*x**2 + 1/12*x**3 + O(x**4)
        """
        return _useries_integrate(s, self.domain, self._prec)
