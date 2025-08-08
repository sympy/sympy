from __future__ import annotations

from typing import Sequence, TYPE_CHECKING

if TYPE_CHECKING:
    from typing import TypeAlias, Union

from sympy.polys.densearith import (
    dup_add,
    dup_mul,
    dup_mul_ground,
    dup_neg,
    dup_sub,
    dup_exquo,
    dup_series_mul,
    dup_series_pow,
    dup_rshift,
)
from sympy.polys.densebasic import dup_degree, dup_reverse, dup_truncate
from sympy.polys.densetools import (
    dup,
    dup_diff,
    dup_integrate,
    dup_revert,
    dup_series_compose,
    dup_series_reversion,
)
from sympy.polys.polyerrors import NotReversible, ExactQuotientFailed
from sympy.polys.domains import Domain, QQ, ZZ
from sympy.polys.domains.domain import Er, Ef
from sympy.polys.domains.field import Field
from sympy.polys.series.base import series_pprint
from sympy.external.gmpy import MPZ, MPQ


USeries: TypeAlias = "tuple[dup[Er], Union[int, None]]"


def _useries(
    coeffs: dup[Er], series_prec: int | None, dom: Domain[Er], ring_prec: int
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


def _useries_valuation(s: USeries[Er], dom: Domain) -> int:
    """
    Returns the valuation of this power series.

    If there are no known nonzero coefficients, returns -1.
    """
    coeffs, _ = s

    if not coeffs:
        return -1

    i = -1
    while dom.is_zero(coeffs[i]):
        i -= 1

    return -i - 1


def _unify_prec(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
) -> tuple[dup[Er], dup[Er], int]:
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


def _useries_equality(s1: USeries[Er], s2: USeries[Er], dom: Domain[Er]) -> bool | None:
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


def _useries_pos(s: USeries[Er], dom: Domain[Er], ring_prec: int) -> USeries[Er]:
    """Return the positive of a power series (which is the same as the series itself)."""
    coeffs, prec = s
    deg = dup_degree(coeffs)

    if prec is None:
        if deg >= ring_prec:
            coeffs = dup_truncate(coeffs, ring_prec, dom)
            prec = ring_prec
    else:
        prec = min(prec, ring_prec)
        if deg >= prec:
            coeffs = dup_truncate(coeffs, prec, dom)

    return coeffs, prec


def _useries_neg(s: USeries[Er], dom: Domain[Er], ring_prec: int) -> USeries[Er]:
    coeffs, prec = s
    neg_coeffs = dup_neg(coeffs, dom)
    s = neg_coeffs, prec
    return _useries_pos(s, dom, ring_prec)


def _useries_add(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
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
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
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
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
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
    s: USeries[Er], n: Er, dom: Domain[Er], ring_prec: int
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


def _useries_div_direct(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    """Direct division when divisor has zero valuation."""
    if not dom.is_unit(s2[0][-1]):
        raise ValueError("Trailing coefficient of the divisor must be a unit")

    _, _, prec = _unify_prec(s1, s2, dom, ring_prec)
    s2_inv = _useries_inverse(s2, dom, prec)
    return _useries_mul(s1, s2_inv, dom, prec)


def _useries_div(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if not coeffs2:
        raise ZeroDivisionError("Series division by zero")

    val1 = _useries_valuation(s1, dom)
    val2 = _useries_valuation(s2, dom)

    if val1 < val2:
        raise ValueError("quotient would not be a power series")

    if prec1 is None and prec2 is None:
        try:
            q = dup_exquo(coeffs1, coeffs2, dom)
            return q, None
        except ExactQuotientFailed:
            pass

    if val2 == 0:
        return _useries_div_direct(s1, s2, dom, ring_prec)
    else:
        # Shift both series to make divisor's valuation zero
        shifted_coeffs1 = dup_rshift(coeffs1, val2, dom)
        shifted_coeffs2 = dup_rshift(coeffs2, val2, dom)

        cap = ring_prec - val2

        prec1 = prec1 - val2 if prec1 is not None else ring_prec
        prec2 = prec2 - val2 if prec2 is not None else ring_prec

        shifted_s1 = _useries(shifted_coeffs1, prec1, dom, cap)
        shifted_s2 = _useries(shifted_coeffs2, prec2, dom, cap)

        q_coeffs = _useries_div_direct(shifted_s1, shifted_s2, dom, cap)[0]

        _, _, prec = _unify_prec(shifted_s1, shifted_s2, dom, ring_prec)
        return q_coeffs, prec


def _useries_pow_int(
    s: USeries[Er], n: int, dom: Domain[Er], ring_prec: int
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


def _useries_truncate(s: USeries[Er], n: int, dom: Domain[Er]) -> USeries[Er]:
    """Truncate a power series to the first n terms."""
    coeffs, prec = s

    if n < 0:
        raise ValueError("Truncation precision must be non-negative")

    deg = dup_degree(coeffs)
    if deg < n:
        return coeffs, prec
    return dup_truncate(coeffs, n, dom), n


def _useries_compose(
    s1: USeries[Er], s2: USeries[Er], dom: Domain[Er], ring_prec: int
) -> USeries[Er]:
    """Compose two power series."""
    coeffs1, prec1 = s1
    coeffs2, prec2 = s2

    if coeffs2 and not dom.is_zero(coeffs2[-1]):
        raise ValueError(
            "Series composition requires the constant term of the second series to be zero."
        )

    if prec1 is None and prec2 is None:
        comp = dup_series_compose(coeffs1, coeffs2, ring_prec, dom)

        deg1 = dup_degree(coeffs1)
        deg2 = dup_degree(coeffs2)

        if deg1 * deg2 < ring_prec:
            return comp, None
        else:
            return comp, ring_prec

    coeffs1, coeffs2, min_prec = _unify_prec(s1, s2, dom, ring_prec)
    comp = dup_series_compose(coeffs1, coeffs2, min_prec, dom)
    return comp, min_prec


def _useries_inverse(s: USeries[Er], dom: Domain[Er], ring_prec: int) -> USeries[Er]:
    """Compute the series multiplicative inverse of a power series."""
    coeffs, prec = s

    if not coeffs or not dom.is_unit(coeffs[-1]):
        raise NotReversible("Series inverse requires the constant term to be a unit")

    if prec is None:
        if len(coeffs) == 1:
            return [dom.revert(coeffs[0])], None
        prec = ring_prec

    inv = dup_revert(coeffs, prec, dom)
    inv = dup_truncate(inv, prec, dom)
    return inv, prec


def _useries_reversion(s: USeries[Er], dom: Domain[Er], ring_prec: int) -> USeries[Er]:
    """Compute the composite inverse of a power series."""
    coeffs, prec = s

    if not coeffs or not dom.is_zero(coeffs[-1]):
        raise NotReversible(
            "Series compositional inverse requires the constant term to be zero."
        )

    if len(coeffs) >= 2 and not dom.is_unit(coeffs[-2]):
        raise NotReversible(
            "Series compositional inverse requires the linear term to be unit."
        )

    if prec is None:
        prec = ring_prec

    series = dup_series_reversion(coeffs, prec, dom)
    return series, prec


def _useries_derivative(s: USeries[Er], dom: Domain[Er], ring_prec: int) -> USeries[Er]:
    """Compute the first derivative of a power series."""
    coeffs, prec = s
    series = dup_diff(coeffs, 1, dom)
    if prec:
        prec -= 1
    return _useries(series, prec, dom, ring_prec)


def _useries_integrate(s: USeries[Ef], dom: Field[Ef], ring_prec: int) -> USeries[Ef]:
    """Compute the integral of a power series."""
    coeffs, prec = s
    series = dup_integrate(coeffs, 1, dom)
    if prec:
        prec += 1
    return _useries(series, prec, dom, ring_prec)


class PythonPowerSeriesRingZZ:
    """
    Python implementation of power series ring over integers :ref:`ZZ`.

    This class provides comprehensive power series operations over the integer ring,
    supporting both series manipulations with precision handling and truncation.

    Parameters
    ==========

    prec : int, optional
        The default precision for power series operations. Default is 6.

    Examples
    ========

    >>> from sympy.polys.series.ringpython import PythonPowerSeriesRingZZ
    >>> R = PythonPowerSeriesRingZZ()
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

    The recommended way to create a power series ring is using the factory function:

    >>> from sympy.polys.series import power_series_ring
    >>> from sympy import ZZ
    >>> R = power_series_ring(ZZ, prec=6)

    This function automatically uses the Flint implementation if available for better
    performance, falling back to the Python implementation otherwise.

    See Also
    ========

    PythonPowerSeriesRingQQ
    FlintPowerSeriesRingZZ
    power_series_ring
    """

    _domain = ZZ

    def __init__(self, prec: int = 6) -> None:
        if prec < 0:
            raise ValueError("Power series precision must be non-negative")
        self._prec = prec

    def __repr__(self) -> str:
        return (
            f"Python Power Series Ring over {self._domain} with precision {self._prec}"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PythonPowerSeriesRingZZ):
            return NotImplemented
        return self._prec == other.prec

    def __hash__(self) -> int:
        return hash((self._domain, self._prec))

    def __call__(
        self, coeffs: Sequence[MPZ | int], prec: int | None = None
    ) -> USeries[MPZ]:
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

    def pretty(
        self, s: USeries[MPZ], *, symbol: str = "x", ascending: bool = True
    ) -> str:
        """Return a pretty-printed string representation of a power series."""
        coeffs, prec = s
        return series_pprint(coeffs, prec, sym=symbol, ascending=ascending)

    def print(
        self, s: USeries[MPZ], *, symbol: str = "x", ascending: bool = True
    ) -> None:
        """Print a pretty-printed representation of a power series."""
        print(self.pretty(s, symbol=symbol, ascending=ascending))

    def from_list(self, coeffs: list[MPZ], prec: int | None = None) -> USeries[MPZ]:
        """
        Create a power series from a list of ground coefficients.

        If `prec` is not specified, it defaults to the ring's precision.
        """
        coeffs = dup_reverse(coeffs)
        if prec is None and len(coeffs) > self._prec:
            prec = self._prec
            coeffs = dup_truncate(coeffs, prec, self._domain)

        return coeffs, prec

    def to_list(self, s: USeries[MPZ]) -> list[MPZ]:
        """Returns the list of series coefficients."""
        coeffs, _ = s
        return coeffs[::-1]

    def equal(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> bool | None:
        """Check if two power series are equal up to their minimum precision."""
        return _useries_equality(s1, s2, self._domain)

    def equal_repr(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> bool:
        """Check if two power series have the same representation."""
        return _useries_equal_repr(s1, s2)

    def positive(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """Return the unary positive of a power series, adjusted to the ring's precision."""
        return _useries_pos(s, self._domain, self._prec)

    def negative(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """Return the unary negative of a power series."""
        return _useries_neg(s, self._domain, self._prec)

    def add(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """Add two power series."""
        return _useries_add(s1, s2, self._domain, self._prec)

    def subtract(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """Subtract two power series."""
        return _useries_sub(s1, s2, self._domain, self._prec)

    def multiply(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """Multiply two power series."""
        return _useries_mul(s1, s2, self._domain, self._prec)

    def multiply_ground(self, s: USeries[MPZ], n: MPZ) -> USeries[MPZ]:
        """Multiply a power series by a ground element."""
        return _useries_mul_ground(s, n, self._domain, self._prec)

    def divide(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """Divide two power series."""
        return _useries_div(s1, s2, self._domain, self._prec)

    def pow_int(self, s: USeries[MPZ], n: int) -> USeries[MPZ]:
        """Raise a power series to a non-negative integer power."""
        return _useries_pow_int(s, n, self._domain, self._prec)

    def square(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """Compute the square of a power series."""
        return _useries_mul(s, s, self._domain, self._prec)

    def compose(self, s1: USeries[MPZ], s2: USeries[MPZ]) -> USeries[MPZ]:
        """Compose two power series, `s1(s2)`."""
        return _useries_compose(s1, s2, self._domain, self._prec)

    def inverse(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """Compute the multiplicative inverse of a power series."""
        return _useries_inverse(s, self._domain, self._prec)

    def reversion(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """Compute the compositional inverse of a power series."""
        return _useries_reversion(s, self._domain, self._prec)

    def truncate(self, s: USeries[MPZ], n: int) -> USeries[MPZ]:
        """Truncate a power series to `n` terms."""
        return _useries_truncate(s, n, self._domain)

    def differentiate(self, s: USeries[MPZ]) -> USeries[MPZ]:
        """Compute the derivative of a power series."""
        return _useries_derivative(s, self._domain, self._prec)


class PythonPowerSeriesRingQQ:
    """
    Python implementation of power series ring over rational field :ref:`QQ`.

    This class provides comprehensive power series operations over the rational field,
    supporting series manipulations with precision handling and truncation.
    It extends the integer ring functionality with support for rational coefficients
    and integration.

    Parameters
    ==========

    prec : int, optional
        The default precision for power series operations. Default is 6.

    Examples
    ========

    >>> from sympy.polys.series.ringpython import PythonPowerSeriesRingQQ
    >>> R = PythonPowerSeriesRingQQ()
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

    The recommended way to create a power series ring is using the factory function:

    >>> from sympy.polys.series import power_series_ring
    >>> from sympy import QQ
    >>> R = power_series_ring(QQ, prec=6)

    This function automatically uses the Flint implementation if available for better
    performance, falling back to the Python implementation otherwise.

    See Also
    ========

    PythonPowerSeriesRingZZ
    FlintPowerSeriesRingQQ
    power_series_ring
    """

    _domain = QQ

    def __init__(self, prec: int = 6) -> None:
        if prec < 0:
            raise ValueError("Power series precision must be non-negative")
        self._prec = prec

    def __repr__(self) -> str:
        return (
            f"Python Power Series Ring over {self._domain} with precision {self._prec}"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PythonPowerSeriesRingQQ):
            return NotImplemented
        return self._prec == other.prec

    def __hash__(self) -> int:
        return hash((self._domain, self._prec))

    def __call__(
        self, coeffs: Sequence[MPQ | int | tuple[int, int]], prec: int | None = None
    ) -> USeries[MPQ]:
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
        return ([self._domain.one, self._domain.zero], None)

    def pretty(
        self, s: USeries[MPQ], *, symbol: str = "x", ascending: bool = True
    ) -> str:
        """Return a pretty-printed string representation of a power series."""
        coeffs, prec = s
        return series_pprint(coeffs, prec, sym=symbol, ascending=ascending)

    def print(
        self, s: USeries[MPQ], *, symbol: str = "x", ascending: bool = True
    ) -> None:
        """Print a pretty-printed representation of a power series."""
        print(self.pretty(s, symbol=symbol, ascending=ascending))

    def from_list(self, coeffs: list[MPQ], prec: int | None = None) -> USeries[MPQ]:
        """
        Create a power series from a list of ground coefficients.

        If `prec` is not specified, it defaults to the ring's precision.
        """
        coeffs = dup_reverse(coeffs)
        if prec is None and len(coeffs) > self._prec:
            prec = self._prec
            coeffs = dup_truncate(coeffs, prec, self._domain)

        return coeffs, prec

    def to_list(self, s: USeries[MPQ]) -> list[MPQ]:
        """Return the list of series coefficients."""
        coeffs, _ = s
        return coeffs[::-1]

    def equal(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> bool | None:
        """Check if two power series are equal up to their minimum precision."""
        return _useries_equality(s1, s2, self._domain)

    def equal_repr(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> bool:
        """Check if two power series have the same representation."""
        return _useries_equal_repr(s1, s2)

    def positive(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Return the unary positive of a power series, adjusted to the ring's precision."""
        return _useries_pos(s, self._domain, self._prec)

    def negative(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Return the unary negative of a power series."""
        return _useries_neg(s, self._domain, self._prec)

    def add(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """Add two power series."""
        return _useries_add(s1, s2, self._domain, self._prec)

    def subtract(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """Subtract two power series."""
        return _useries_sub(s1, s2, self._domain, self._prec)

    def multiply(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """Multiply two power series."""
        return _useries_mul(s1, s2, self._domain, self._prec)

    def multiply_ground(self, s: USeries[MPQ], n: MPQ) -> USeries[MPQ]:
        """Multiply a power series by a ground element."""
        return _useries_mul_ground(s, n, self._domain, self._prec)

    def divide(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """Divide two power series."""
        return _useries_div(s1, s2, self._domain, self._prec)

    def pow_int(self, s: USeries[MPQ], n: int) -> USeries[MPQ]:
        """Raise a power series to a non-negative integer power."""
        return _useries_pow_int(s, n, self._domain, self._prec)

    def square(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the square of a power series."""
        return _useries_mul(s, s, self._domain, self._prec)

    def compose(self, s1: USeries[MPQ], s2: USeries[MPQ]) -> USeries[MPQ]:
        """Compose two power series, `s1(s2)`."""
        return _useries_compose(s1, s2, self._domain, self._prec)

    def inverse(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the multiplicative inverse of a power series."""
        return _useries_inverse(s, self._domain, self._prec)

    def reversion(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the compositional inverse of a power series."""
        return _useries_reversion(s, self._domain, self._prec)

    def truncate(self, s: USeries[MPQ], n: int) -> USeries[MPQ]:
        """Truncate a power series to `n` terms."""
        return _useries_truncate(s, n, self._domain)

    def differentiate(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the derivative of a power series."""
        return _useries_derivative(s, self._domain, self._prec)

    def integrate(self, s: USeries[MPQ]) -> USeries[MPQ]:
        """Compute the integral of a power series."""
        return _useries_integrate(s, self._domain, self._prec)
