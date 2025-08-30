from __future__ import annotations

from contextlib import contextmanager
from typing import Sequence, Union, TYPE_CHECKING

from sympy.polys.densebasic import dup, dup_reverse
from sympy.polys.domains import Domain, QQ, ZZ
from sympy.polys.series.ringpython import _useries_valuation
from sympy.polys.series.base import series_pprint
from sympy.external.gmpy import GROUND_TYPES, MPZ, MPQ
from sympy.utilities.decorator import doctest_depends_on
from sympy.polys.polyerrors import NotReversible


if TYPE_CHECKING:
    from flint import fmpq_poly, fmpq_series, fmpz_poly, fmpz_series, ctx  # type: ignore
    from flint.utils.flint_exceptions import DomainError  # type: ignore
elif GROUND_TYPES == "flint":
    from flint import fmpq_poly, fmpq_series, fmpz_poly, fmpz_series, ctx
    from flint.utils.flint_exceptions import DomainError
else:
    fmpq_poly = fmpq_series = fmpz_poly = fmpz_series = ctx = None


ZZSeries = Union[fmpz_series, fmpz_poly]
QQSeries = Union[fmpq_series, fmpq_poly]


def _get_series_precision(s: fmpz_series | fmpq_series) -> int:
    """Helper function to get the precision of a series. By using the
    representation of the ring"""

    prec = getattr(s, "prec", None)
    if prec is not None:
        return prec

    # XXX: This approach is inefficient, but for python-flint <= 0.7.1, there is
    # no alternative method to extract the precision from a series element.
    # This function could be removed when python-flint 0.8 becomes the
    # minimum supported version.
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


@doctest_depends_on(ground_types=["flint"])
class FlintPowerSeriesRingZZ:
    """
    Flint implementation of power series ring over integers :ref:`ZZ`.

    This class provides high-performance power series operations over the integer ring,
    leveraging the FLINT library for optimized arithmetic and series manipulations
    precision handling and truncation.

    Parameters
    ==========

    prec : int, optional
        The default precision for power series operations. Default is 6.

    Examples
    ========

    >>> from sympy.polys.series.ringflint import FlintPowerSeriesRingZZ
    >>> R = FlintPowerSeriesRingZZ()
    >>> s = R([1, 2, 3])  # 1 + 2*x + 3*x^2
    >>> R.print(s)
    1 + 2*x + 3*x**2

    >>> s_pow = R.pow_int(s, 2)  # Square the series
    >>> R.print(s_pow)
    1 + 4*x + 10*x**2 + 12*x**3 + 9*x**4

    >>> s_inv = R.inverse(R([1, 1]))  # Inverse of 1 + x
    >>> R.print(s_inv)
    1 - x + x**2 - x**3 + x**4 - x**5 + O(x**6)

    Note
    ====

    The recommended way to create a power series ring is using the factory function
    which returns a new instance of the higher level PowerSeriesRing class with
    the ring generator:

    >>> from sympy.polys.series import power_series_ring
    >>> from sympy.polys.domains import ZZ
    >>> R, x = power_series_ring("x", ZZ, 6)
    >>> R
    Power Series Ring in x over ZZ of size 6
    >>> type(x)
    <class 'sympy.polys.series.ring.PowerSeriesElement'>

    This function automatically uses the Flint implementation if available for better
    performance, falling back to the Python implementation otherwise.

    See Also
    ========

    sympy.polys.series.ringflint.FlintPowerSeriesRingQQ
    sympy.polys.series.ringpython.PythonPowerSeriesRingZZ
    sympy.polys.series.ring.power_series_ring
    sympy.polys.series.ring.PowerSeriesRingRing
    sympy.polys.series.ring.PowerSeriesRingField
    sympy.polys.series.ring.PowerSeriesElement
    """

    _domain = ZZ

    def __init__(self, prec: int = 6) -> None:
        if prec < 0:
            raise ValueError("Power series precision must be non-negative")
        self._prec = prec

    def __repr__(self) -> str:
        return (
            f"Flint Power Series Ring over {self._domain} with precision {self._prec}"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, FlintPowerSeriesRingZZ):
            return NotImplemented
        return self._prec == other.prec

    def __hash__(self) -> int:
        return hash((self._domain, self._prec))

    def __call__(
        self, coeffs: Sequence[MPZ | int], prec: int | None = None
    ) -> ZZSeries:
        """
        Create a power series from a list of coefficients.

        If `prec` is not specified, it defaults to the ring's precision.
        """
        s: list[MPZ] = []
        for c in coeffs:
            if isinstance(c, MPZ):
                s.append(c)
            elif isinstance(c, int):
                s.append(self._domain(c))
            else:
                raise TypeError(f"Unsupported coefficient type: {type(c)}")

        return self.from_list(s, prec)

    @property
    def domain(self) -> Domain[MPZ]:
        """Return the ground domain of the power series ring."""
        return self._domain

    @property
    def prec(self) -> int:
        """Return the ring's precision."""
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

    def pretty(self, s: ZZSeries, *, symbol: str = "x", ascending: bool = True) -> str:
        """Return a pretty-printed string representation of a power series."""
        coeffs = dup_reverse(s.coeffs(), ZZ)

        if isinstance(s, fmpz_poly):
            return series_pprint(coeffs, None)

        prec = self.series_prec(s)
        return series_pprint(coeffs, prec, sym=symbol, ascending=ascending)

    def print(self, s: ZZSeries, *, symbol: str = "x", ascending: bool = True) -> None:
        """Print a pretty-printed representation of a power series."""
        print(self.pretty(s, symbol=symbol, ascending=ascending))

    def from_list(self, coeffs: list[MPZ], prec: int | None = None) -> ZZSeries:
        """
        Create a power series from a list of ground coefficients.

        If `prec` is not specified, it defaults to the ring's precision.
        """
        if prec is None:
            if len(coeffs) <= self._prec:
                return fmpz_poly(coeffs)
            prec = self._prec
        else:
            prec = min(prec, self._prec)

        return fmpz_series(coeffs, prec=prec)

    def to_list(self, s: ZZSeries) -> list[MPZ]:
        """Returns the list of series coefficients."""
        return s.coeffs()

    def to_dense(self, s: ZZSeries) -> dup[MPZ]:
        """Return the coefficients of a power series as a dense list."""
        return s.coeffs()[::-1]

    def series_prec(self, s: ZZSeries) -> int | None:
        """Return the precision of the series."""
        if isinstance(s, fmpz_poly):
            return None
        return _get_series_precision(s)

    def equal(self, s1: ZZSeries, s2: ZZSeries) -> bool | None:
        """Check if two power series are equal up to their minimum precision."""
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            return s1 == s2
        elif isinstance(s1, fmpz_poly):
            min_prec = self.series_prec(s2)
        elif isinstance(s2, fmpz_poly):
            min_prec = self.series_prec(s1)
        else:
            min_prec = min(_get_series_precision(s1), _get_series_precision(s2))

        coeffs1 = s1.coeffs()[:min_prec]
        coeffs2 = s2.coeffs()[:min_prec]

        if coeffs1 != coeffs2:
            return False
        return None

    def equal_repr(self, s1: ZZSeries, s2: ZZSeries) -> bool:
        """Check if two power series have the same representation."""
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            return s1 == s2
        elif isinstance(s1, fmpz_series) and isinstance(s2, fmpz_series):
            return s1._equal_repr(s2)
        else:
            return False

    def is_ground(self, arg: ZZSeries) -> bool | None:
        """Check if a arg is a ground element of the power series ring."""
        if self.prec == 0:
            return None

        return len(arg) <= 1

    def constant_coefficient(self, s: ZZSeries) -> MPZ:
        """Return the constant coefficient of a power series."""
        return s[0]

    def positive(self, s: ZZSeries) -> ZZSeries:
        """Return the unary positive of a power series, adjusted to the ring's precision."""
        ring_prec = self._prec
        if isinstance(s, fmpz_poly):
            if s.degree() >= ring_prec:
                return fmpz_series(s, prec=ring_prec)
            return s

        # XXX: This shold simply be: fmpz_series(s, prec=ring_prec)
        # https://github.com/flintlib/python-flint/issues/304
        prec = min(_get_series_precision(s), ring_prec)
        return fmpz_series(s.coeffs(), prec=prec)

    def negative(self, s: ZZSeries) -> ZZSeries:
        """Return the unary negative of a power series."""
        with _global_cap(self._prec):
            return self.positive(-s)

    def add(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Add two power series."""
        ring_prec = self._prec
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            poly = s1 + s2
            if poly.degree() < ring_prec:
                return poly
            return fmpz_series(poly, prec=ring_prec)

        with _global_cap(ring_prec):
            return s1 + s2

    def add_ground(self, s: ZZSeries, n: MPZ) -> ZZSeries:
        """Add a ground element to a power series."""
        with _global_cap(self._prec):
            return s + n

    def subtract(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Subtract two power series."""
        ring_prec = self._prec
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            poly = s1 - s2
            if poly.degree() < ring_prec:
                return poly
            return fmpz_series(poly, prec=ring_prec)

        with _global_cap(ring_prec):
            return s1 - s2

    def subtract_ground(self, s: ZZSeries, n: MPZ) -> ZZSeries:
        """Subtract a ground element from a power series."""
        with _global_cap(self._prec):
            return s - n

    def rsubtract_ground(self, s: ZZSeries, n: MPZ) -> ZZSeries:
        """Subtract a power series from a ground element."""
        with _global_cap(self._prec):
            return n - s

    def multiply(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Multiply two power series."""
        ring_prec = self._prec
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            deg1 = s1.degree()
            deg2 = s2.degree()

            if deg1 + deg2 < ring_prec:
                return s1 * s2
            else:
                if deg1 >= ring_prec:
                    s1 = self.truncate(s1, ring_prec)
                if deg2 >= ring_prec:
                    s2 = self.truncate(s2, ring_prec)

                return fmpz_series(s1 * s2, prec=ring_prec)

        with _global_cap(ring_prec):
            return s1 * s2

    def multiply_ground(self, s: ZZSeries, n: MPZ) -> ZZSeries:
        """Multiply a power series by a ground element."""
        ring_prec = self._prec
        if isinstance(s, fmpz_poly):
            poly = s * n
            if poly.degree() < ring_prec:
                return poly
            return fmpz_series(poly, prec=ring_prec)

        with _global_cap(ring_prec):
            return s * n

    def divide(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Divide two power series."""
        ring_prec = self._prec
        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            try:
                return s1 / s2
            except DomainError:
                ring_prec = ring_prec + _useries_valuation(
                    (s2.coeffs()[::-1], None), self._domain
                )
                s1 = fmpz_series(s1, prec=ring_prec)
                s2 = fmpz_series(s2, prec=ring_prec)

        with _global_cap(ring_prec):
            return s1 / s2

    def pow_int(self, s: ZZSeries, n: int) -> ZZSeries:
        """Raise a power series to a integer power."""
        ring_prec = self._prec
        if n < 0:
            n = -n
            s = self.pow_int(s, n)
            try:
                inv = self.inverse(s)
                return inv
            except NotReversible:
                raise ValueError("Result would not be a power series")

        if isinstance(s, fmpz_poly):
            if s.degree() * n < ring_prec:
                return s**n

            if s.degree() > ring_prec:
                s = self.truncate(s, ring_prec)
            poly = s**n
            return fmpz_series(poly, prec=ring_prec)

        with _global_cap(ring_prec):
            return s**n

    def square(self, s: ZZSeries) -> ZZSeries:
        """Compute the square of a power series."""
        return self.pow_int(s, 2)

    def compose(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Compose two power series, `s1(s2)`."""
        dom: Domain[MPZ] = self._domain
        ring_prec: int = self._prec

        if s2 and not dom.is_zero(s2[0]):
            raise ValueError(
                "Series composition requires the second series to have a zero constant term."
            )

        if isinstance(s1, fmpz_poly) and isinstance(s2, fmpz_poly):
            deg1: int = s1.degree()
            deg2: int = s2.degree()

            if deg1 * deg2 < ring_prec:
                return s1(s2)

            if deg1 <= ring_prec and deg2 <= ring_prec:
                return self.truncate(s1(s2), ring_prec)
            else:
                s1 = fmpz_series(s1, prec=ring_prec)
                s2 = fmpz_series(s2, prec=ring_prec)

                with _global_cap(ring_prec):
                    return s1(s2)

        if isinstance(s1, fmpz_poly):
            prec2 = self.series_prec(s2)
            s1 = fmpz_series(s1, prec=prec2)

        if isinstance(s2, fmpz_poly):
            prec1 = self.series_prec(s1)
            s2 = fmpz_series(s2, prec=prec1)

        with _global_cap(ring_prec):
            return s1(s2)

    def inverse(self, s: ZZSeries) -> ZZSeries:
        """Compute the multiplicative inverse of a power series."""
        dom: Domain[MPZ] = self._domain
        ring_prec: int = self._prec

        if not s or not dom.is_unit(s[0]):
            raise NotReversible(
                "Series inverse requires the constant term to be a unit"
            )

        if isinstance(s, fmpz_poly):
            if len(s) == 1:
                return 1 / s
            s = fmpz_series(s, prec=ring_prec)

        with _global_cap(ring_prec):
            return 1 / s

    def reversion(self, s: ZZSeries) -> ZZSeries:
        """Compute the compositional inverse of a power series."""
        dom = self._domain

        if not s or not dom.is_zero(s[0]):
            raise NotReversible(
                "Series compositional inverse requires the constant term to be zero."
            )

        if len(s) >= 2 and not dom.is_unit(s[1]):
            raise NotReversible(
                "Series compositional inverse requires the linear term to be unit."
            )

        if isinstance(s, fmpz_poly):
            s = fmpz_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.reversion()

    def truncate(self, s: ZZSeries, n: int) -> ZZSeries:
        """Truncate a power series to `n` terms."""
        if n < 0:
            raise ValueError("Truncation precision must be non-negative")

        if len(s) <= n:
            return s

        # XXX: This should simply be: return fmpz_series(s, prec=n)
        # https://github.com/flintlib/python-flint/issues/304
        coeffs = s.coeffs()[:n]
        return fmpz_series(coeffs, prec=n)

    def differentiate(self, s: ZZSeries) -> ZZSeries:
        """Compute the derivative of a power series."""
        if isinstance(s, fmpz_poly):
            poly = s.derivative()
            if poly.degree() < self._prec:
                return poly
            return fmpz_series(poly, prec=self._prec)

        poly = fmpz_poly(s.coeffs())
        derivative = poly.derivative()
        prec = min(_get_series_precision(s) - 1, self._prec)
        return fmpz_series(derivative, prec=prec)


@doctest_depends_on(ground_types=["flint"])
class FlintPowerSeriesRingQQ:
    """
    Flint implementation of power series ring over rational field :ref:`QQ`.

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

    >>> from sympy.polys.series.ringflint import FlintPowerSeriesRingQQ
    >>> R = FlintPowerSeriesRingQQ()
    >>> s = R([1, (1, 2), (1, 3)])  # 1 + x/2 + x^2/3
    >>> R.print(s)
    1 + 1/2*x + 1/3*x**2

    >>> s_int = R.integrate(s)  # Integration
    >>> R.print(s_int)
    x + 1/4*x**2 + 1/9*x**3

    >>> s_inv = R.inverse(R([1, (1, 2)]))  # Inverse of 1 + x/2
    >>> R.print(s_inv)
    1 - 1/2*x + 1/4*x**2 - 1/8*x**3 + 1/16*x**4 - 1/32*x**5 + O(x**6)

    Note
    ====

    The recommended way to create a power series ring is using the factory function
    which returns a new instance of the higher level PowerSeriesRing class with
    the ring generator:

    >>> from sympy.polys.series import power_series_ring
    >>> from sympy.polys.domains import QQ
    >>> R, x = power_series_ring("x", QQ, 6)
    >>> R
    Power Series Ring in x over QQ of size 6
    >>> type(x)
    <class 'sympy.polys.series.ring.PowerSeriesElement'>

    This function automatically uses the Flint implementation if available for better
    performance, falling back to the Python implementation otherwise.

    See Also
    ========

    sympy.polys.series.ringflint.FlintPowerSeriesRingZZ
    sympy.polys.series.ringpython.PythonPowerSeriesRingQQ
    sympy.polys.series.ring.power_series_ring
    sympy.polys.series.ring.PowerSeriesRingRing
    sympy.polys.series.ring.PowerSeriesRingField
    sympy.polys.series.ring.PowerSeriesElement
    """

    _domain = QQ

    def __init__(self, prec: int = 6) -> None:
        if prec < 0:
            raise ValueError("Power series precision must be non-negative")
        self._prec = prec

    def __repr__(self) -> str:
        return (
            f"Flint Power Series Ring over {self._domain} with precision {self._prec}"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, FlintPowerSeriesRingQQ):
            return NotImplemented
        return self._prec == other.prec

    def __hash__(self) -> int:
        return hash((self._domain, self._prec))

    def __call__(
        self, coeffs: Sequence[MPQ | int | tuple[int, int]], prec: int | None = None
    ) -> QQSeries:
        """
        Create a power series from a list of coefficients.

        If `prec` is not specified, it defaults to the ring's precision.
        """
        s: list[MPQ] = []
        for c in coeffs:
            if isinstance(c, MPQ):
                s.append(c)
            elif isinstance(c, int):
                s.append(self._domain(c))
            elif isinstance(c, tuple):
                s.append(self._domain(*c))
            else:
                raise TypeError(f"Unsupported coefficient type: {type(c)}")

        return self.from_list(s, prec)

    @property
    def domain(self) -> Domain[MPQ]:
        """Return the ground domain of the power series ring."""
        return self._domain

    @property
    def prec(self) -> int:
        """Return the ring's precision."""
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

    def pretty(self, s: QQSeries, *, symbol: str = "x", ascending: bool = True) -> str:
        """Return a pretty-printed string representation of a power series."""
        coeffs = dup_reverse(s.coeffs(), QQ)

        if isinstance(s, fmpq_poly):
            return series_pprint(coeffs, None)

        prec = self.series_prec(s)
        return series_pprint(coeffs, prec, sym=symbol, ascending=ascending)

    def print(self, s: QQSeries, *, symbol: str = "x", ascending: bool = True) -> None:
        """Print a pretty-printed representation of a power series."""
        print(self.pretty(s, symbol=symbol, ascending=ascending))

    def from_list(self, coeffs: list[MPQ], prec: int | None = None) -> QQSeries:
        """
        Create a power series from a list of ground coefficients.

        If `prec` is not specified, it defaults to the ring's precision.
        """
        if prec is None:
            if len(coeffs) <= self._prec:
                return fmpq_poly(coeffs)
            prec = self._prec
        else:
            prec = min(prec, self._prec)

        return fmpq_series(coeffs, prec=prec)

    def to_list(self, s: QQSeries) -> list[MPQ]:
        """Returns the list of series coefficients."""
        return s.coeffs()

    def to_dense(self, s: QQSeries) -> dup[MPQ]:
        """Return the coefficients of a power series as a dense list."""
        return s.coeffs()[::-1]

    def series_prec(self, s: QQSeries) -> int | None:
        """Return the precision of the series."""
        if isinstance(s, fmpq_poly):
            return None
        return _get_series_precision(s)

    def equal(self, s1: QQSeries, s2: QQSeries) -> bool | None:
        """Check if two power series are equal up to their minimum precision."""
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            return s1 == s2
        elif isinstance(s1, fmpq_poly):
            min_prec = self.series_prec(s2)
        elif isinstance(s2, fmpq_poly):
            min_prec = self.series_prec(s1)
        else:
            min_prec = min(_get_series_precision(s1), _get_series_precision(s2))

        coeffs1 = s1.coeffs()[:min_prec]
        coeffs2 = s2.coeffs()[:min_prec]

        if coeffs1 != coeffs2:
            return False
        return None

    def equal_repr(self, s1: QQSeries, s2: QQSeries) -> bool:
        """Check if two power series have the same representation."""
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            return s1 == s2
        elif isinstance(s1, fmpq_series) and isinstance(s2, fmpq_series):
            return s1._equal_repr(s2)
        else:
            return False

    def is_ground(self, arg: QQSeries) -> bool | None:
        """Check if a arg is a ground element of the power series ring."""
        if self.prec == 0:
            return None

        return len(arg) <= 1

    def constant_coefficient(self, s: QQSeries) -> MPQ:
        """Return the constant coefficient of a power series."""
        return s[0]

    def positive(self, s: QQSeries) -> QQSeries:
        """Return the unary positive of a power series, adjusted to the ring's precision."""
        ring_prec = self._prec
        if isinstance(s, fmpq_poly):
            if s.degree() >= ring_prec:
                return fmpq_series(s, prec=ring_prec)
            return s

        # XXX: This should simply be: fmpq_series(s, prec=ring_prec)
        # https://github.com/flintlib/python-flint/issues/304
        prec = min(_get_series_precision(s), ring_prec)
        return fmpq_series(s.coeffs(), prec=prec)

    def negative(self, s: QQSeries) -> QQSeries:
        """Return the unary negative of a power series."""
        with _global_cap(self._prec):
            return self.positive(-s)

    def add(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Add two power series."""
        ring_prec = self._prec
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            poly = s1 + s2
            if poly.degree() < ring_prec:
                return poly
            return fmpq_series(poly, prec=ring_prec)

        with _global_cap(ring_prec):
            return s1 + s2

    def add_ground(self, s: QQSeries, n: MPQ) -> QQSeries:
        """Add a ground element to a power series."""
        with _global_cap(self._prec):
            return s + n

    def subtract(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Subtract two power series."""
        ring_prec = self._prec
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            poly = s1 - s2
            if poly.degree() < ring_prec:
                return poly
            return fmpq_series(poly, prec=ring_prec)

        with _global_cap(ring_prec):
            return s1 - s2

    def subtract_ground(self, s: QQSeries, n: MPQ) -> QQSeries:
        """Subtract a ground element from a power series."""
        with _global_cap(self._prec):
            return s - n

    def rsubtract_ground(self, s: QQSeries, n: MPQ) -> QQSeries:
        """Subtract a power series from a ground element."""
        with _global_cap(self._prec):
            return n - s

    def multiply(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Multiply two power series."""
        ring_prec = self._prec
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            deg1 = s1.degree()
            deg2 = s2.degree()

            if deg1 + deg2 < ring_prec:
                return s1 * s2
            else:
                if deg1 >= ring_prec:
                    s1 = self.truncate(s1, ring_prec)
                if deg2 >= ring_prec:
                    s2 = self.truncate(s2, ring_prec)

                return fmpq_series(s1 * s2, prec=ring_prec)

        with _global_cap(ring_prec):
            return s1 * s2

    def multiply_ground(self, s: QQSeries, n: MPQ) -> QQSeries:
        """Multiply a power series by a ground element."""
        ring_prec = self._prec
        if isinstance(s, fmpq_poly):
            poly = s * n
            if poly.degree() < ring_prec:
                return poly
            return fmpq_series(poly, prec=ring_prec)

        with _global_cap(ring_prec):
            return s * n

    def divide(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Divide two power series."""
        ring_prec = self._prec
        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            try:
                return s1 / s2
            except DomainError:
                ring_prec = ring_prec + _useries_valuation(
                    (s2.coeffs()[::-1], None), self._domain
                )
                s1 = fmpq_series(s1, prec=ring_prec)
                s2 = fmpq_series(s2, prec=ring_prec)

        with _global_cap(ring_prec):
            return s1 / s2

    def pow_int(self, s: QQSeries, n: int) -> QQSeries:
        """Raise a power series to a integer power."""
        ring_prec = self._prec
        if n < 0:
            n = -n
            s = self.pow_int(s, n)
            try:
                inv = self.inverse(s)
                return inv
            except NotReversible:
                raise ValueError("Result would not be a power series")

        if isinstance(s, fmpq_poly):
            if s.degree() * n < ring_prec:
                return s**n

            if s.degree() > ring_prec:
                s = self.truncate(s, ring_prec)
            poly = s**n
            return fmpq_series(poly, prec=ring_prec)

        with _global_cap(ring_prec):
            return s**n

    def square(self, s: QQSeries) -> QQSeries:
        """Compute the square of a power series."""
        return self.pow_int(s, 2)

    def sqrt(self, s: QQSeries) -> QQSeries:
        """Compute the square root of a power series."""
        if not s:
            return s

        if isinstance(s, fmpq_poly):
            try:
                return s.sqrt()
            except DomainError:
                s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.sqrt()

    def compose(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Compose two power series, `s1(s2)`."""
        dom: Domain[MPQ] = self._domain
        ring_prec: int = self._prec

        if s2 and not dom.is_zero(s2[0]):
            raise ValueError(
                "Series composition requires the second series to have a zero constant term."
            )

        if isinstance(s1, fmpq_poly) and isinstance(s2, fmpq_poly):
            deg1: int = s1.degree()
            deg2: int = s2.degree()

            if deg1 * deg2 < ring_prec:
                return s1(s2)

            if deg1 <= ring_prec and deg2 <= ring_prec:
                return self.truncate(s1(s2), ring_prec)
            else:
                s1 = fmpq_series(s1, prec=ring_prec)
                s2 = fmpq_series(s2, prec=ring_prec)

                with _global_cap(ring_prec):
                    return s1(s2)

        if isinstance(s1, fmpq_poly):
            prec2 = self.series_prec(s2)
            s1 = fmpq_series(s1, prec=prec2)

        if isinstance(s2, fmpq_poly):
            prec1 = self.series_prec(s1)
            s2 = fmpq_series(s2, prec=prec1)

        with _global_cap(ring_prec):
            return s1(s2)

    def inverse(self, s: QQSeries) -> QQSeries:
        """Compute the multiplicative inverse of a power series."""
        dom: Domain[MPQ] = self._domain
        ring_prec: int = self._prec

        if not s or not dom.is_unit(s[0]):
            raise NotReversible(
                "Series inverse requires the constant term to be a unit"
            )

        if isinstance(s, fmpq_poly):
            if len(s) == 1:
                return 1 / s
            s = fmpq_series(s, prec=ring_prec)

        with _global_cap(ring_prec):
            return s.inv()

    def reversion(self, s: QQSeries) -> QQSeries:
        """Compute the compositional inverse of a power series."""
        dom = self._domain

        if not s or not dom.is_zero(s[0]):
            raise NotReversible(
                "Series compositional inverse requires the constant term to be zero."
            )

        if len(s) >= 2 and not dom.is_unit(s[1]):
            raise NotReversible(
                "Series compositional inverse requires the linear term to be unit."
            )

        if isinstance(s, fmpq_poly):
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.reversion()

    def truncate(self, s: QQSeries, n: int) -> QQSeries:
        """Truncate a power series to `n` terms."""
        if n < 0:
            raise ValueError("Truncation precision must be non-negative")

        if len(s) <= n:
            return s

        # XXX: This should simply be: return fmpq_series(s, prec=n)
        # https://github.com/flintlib/python-flint/issues/304
        coeffs = s.coeffs()[:n]
        return fmpq_series(coeffs, prec=n)

    def differentiate(self, s: QQSeries) -> QQSeries:
        """Compute the derivative of a power series."""
        if isinstance(s, fmpq_poly):
            poly = s.derivative()
            if poly.degree() < self._prec:
                return poly
            return fmpq_series(poly, prec=self._prec)

        poly = fmpq_poly(s.coeffs())
        derivative = poly.derivative()
        prec = min(_get_series_precision(s) - 1, self._prec)
        return fmpq_series(derivative, prec=prec)

    def integrate(self, s: QQSeries) -> QQSeries:
        """Compute the integral of a power series."""
        if isinstance(s, fmpq_poly):
            poly = s.integral()
            if poly.degree() < self._prec:
                return poly
            return fmpq_series(poly, prec=self._prec)
        else:
            s_prec = _get_series_precision(s)
            if s_prec >= self._prec:
                s = fmpq_series(s, prec=self._prec - 1)

        with _global_cap(self._prec):
            return s.integral()

    def log(self, s: QQSeries) -> QQSeries:
        """Compute the logarithm of a power series."""
        if isinstance(s, fmpq_poly):
            if s == 1:
                return fmpq_poly([])
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.log()

    def log1p(self, s: QQSeries) -> QQSeries:
        """Compute the logarithm of (1 + s) for a power series with zero constant term."""
        if s[0] != 0:
            raise ValueError(
                "Log1p requires the constant term of the input series to be zero."
            )

        with _global_cap(self._prec):
            s = s + 1

        return self.log(s)

    def exp(self, s: QQSeries) -> QQSeries:
        """Compute the exponential of a power series."""
        if isinstance(s, fmpq_poly):
            if not s:
                return fmpq_poly([1])
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.exp()

    def expm1(self, s: QQSeries) -> QQSeries:
        """Compute the exponential of a power series minus 1."""
        if isinstance(s, fmpq_poly):
            if not s:
                return fmpq_poly([])
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.exp() - 1

    def atan(self, s: QQSeries) -> QQSeries:
        """Compute the arctangent of a power series."""
        if not s:
            return s

        if isinstance(s, fmpq_poly):
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.atan()

    def atanh(self, s: QQSeries) -> QQSeries:
        """Compute the hyperbolic arctangent of a power series."""
        if not s:
            return s

        if isinstance(s, fmpq_poly):
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.atanh()

    def asin(self, s: QQSeries) -> QQSeries:
        """Compute the arcsine of a power series."""
        if not s:
            return s

        if isinstance(s, fmpq_poly):
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.asin()

    def asinh(self, s: QQSeries) -> QQSeries:
        """Compute the hyperbolic arcsine of a power series."""
        if not s:
            return s

        if isinstance(s, fmpq_poly):
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.asinh()

    def tan(self, s: QQSeries) -> QQSeries:
        """Compute the tangent of a power series."""
        if not s:
            return s

        if isinstance(s, fmpq_poly):
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.tan()

    def tanh(self, s: QQSeries) -> QQSeries:
        """Compute the hyperbolic tangent of a power series."""
        if not s:
            return s

        if isinstance(s, fmpq_poly):
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.tanh()

    def sin(self, s: QQSeries) -> QQSeries:
        """Compute the sine of a power series."""
        if not s:
            return s

        if isinstance(s, fmpq_poly):
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.sin()

    def sinh(self, s: QQSeries) -> QQSeries:
        """Compute the hyperbolic sine of a power series."""
        if not s:
            return s

        if isinstance(s, fmpq_poly):
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.sinh()

    def cos(self, s: QQSeries) -> QQSeries:
        """Compute the cosine of a power series."""
        if isinstance(s, fmpq_poly):
            if not s:
                return fmpq_poly([1])
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.cos()

    def cosh(self, s: QQSeries) -> QQSeries:
        """Compute the hyperbolic cosine of a power series."""
        if isinstance(s, fmpq_poly):
            if not s:
                return fmpq_poly([1])
            s = fmpq_series(s, prec=self._prec)

        with _global_cap(self._prec):
            return s.cosh()

    def hypot(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Compute the hypotenuse of two power series."""

        with _global_cap(self._prec):
            sqr_s1 = s1**2
            sqr_s2 = s2**2
            add = sqr_s1 + sqr_s2

        return self.sqrt(add)
