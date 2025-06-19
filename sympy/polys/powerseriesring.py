from __future__ import annotations

from typing import Any, Protocol, TypeAlias, runtime_checkable

from sympy.external.gmpy import GROUND_TYPES
from sympy.polys.domains import Domain
from sympy.utilities import public


ZZSeries: TypeAlias = Any
QQSeries: TypeAlias = Any


def series_from_list(series: Any, prec: int | None = None) -> str:
    terms = []

    for i, coeff in enumerate(series):
        if coeff == 0:
            continue
        if coeff == 1 and i != 0:
            term = "x" if i == 1 else f"x**{i}"
        elif coeff == -1 and i != 0:
            term = "-x" if i == 1 else f"-x**{i}"
        else:
            if i == 0:
                term = f"{coeff}"
            elif i == 1:
                term = f"{coeff}*x"
            else:
                term = f"{coeff}*x**{i}"
        terms.append((i, term))

    if not terms:
        return "0"

    # Sort terms based on exponent depending on `prec`
    terms.sort(reverse=(prec is None))  # Descending if prec is None
    poly = " + ".join(term for _, term in terms).replace("+ -", "- ")

    if prec is not None:
        return f"{poly} + O(x**{prec})"
    return poly

@public
def power_series_ring(
    domain: Domain, prec: int = 6
) -> PowerSeriesRingZZ | PowerSeriesRingQQ:
    """
    Create a power series ring over the given domain.

    Parameters
    ----------
    domain : Domain
        The ground domain for the power series ring. Must be ZZ or QQ.
    prec : int, optional
        The default precision for power series operations. Default is 6.

    Returns
    -------
    PowerSeriesRingZZ | PowerSeriesRingQQ
        A power series ring instance. The actual implementation (Flint or Python)
        depends on whether python-flint is available.

    Examples
    --------
    >>> from sympy.polys.domains import ZZ, QQ
    >>> from sympy.polys.powerseriesring import power_series_ring
    >>> R_ZZ = power_series_ring(ZZ, 5)
    >>> R_QQ = power_series_ring(QQ, 10)
    """
    flint: bool | None = None
    if GROUND_TYPES == 'flint':
        flint = True

    if domain.is_ZZ:
        if flint:
            from sympy.polys.flint_powerseriesring import FlintPowerSeriesRingZZ
            return FlintPowerSeriesRingZZ(prec)
        from sympy.polys.python_powerseriesring import PythonPowerSeriesRingZZ
        return PythonPowerSeriesRingZZ(prec)
    if domain.is_QQ:
        if flint:
            from sympy.polys.flint_powerseriesring import FlintPowerSeriesRingQQ
            return FlintPowerSeriesRingQQ(prec)
        from sympy.polys.python_powerseriesring import PythonPowerSeriesRingQQ
        return PythonPowerSeriesRingQQ(prec)
    raise TypeError("Ground domain must be an instance of QQ or ZZ")


@runtime_checkable
class PowerSeriesRingZZ(Protocol):
    """Protocol for power series ring over integer ring.

    Attributes
    ----------
    domain : Domain
        The ground domain (ZZ).
    prec : int
        The default precision for power series operations.
    """

    domain: Domain
    prec: int

    def __init__(self, prec: int = 6) -> None:
        """Initialize a power series ring over ZZ."""
        ...

    def __repr__(self) -> str:
        """Return string representation of the ring."""
        ...

    @property
    def one(self) -> ZZSeries:
        """Return the multiplicative identity (1) as a power series."""
        ...

    @property
    def zero(self) -> ZZSeries:
        """Return the additive identity (0) as a power series."""
        ...

    @property
    def gen(self) -> ZZSeries:
        """Return the generator (x) as a power series."""
        ...

    def pretty(self, s: ZZSeries) -> str:
        """Return a pretty string representation of a power series."""
        ...

    def print(self, s: ZZSeries) -> str:
        """Return a printable string representation of a power series."""
        ...

    def equal(self, s1: ZZSeries, s2: ZZSeries) -> bool:
        """Check if two power series are equal."""
        ...

    def negative(self, s: ZZSeries) -> ZZSeries:
        """Return the additive inverse of a power series."""
        ...

    def add(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Add two power series."""
        ...

    def subtract(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Subtract two power series."""
        ...

    def multiply(self, s1: ZZSeries, s2: ZZSeries) -> ZZSeries:
        """Multiply two power series."""
        ...

    def multiply_ground(self, s: ZZSeries, c: Any) -> ZZSeries:
        """Multiply a power series by a ground element."""
        ...

    def trunc(self, s: ZZSeries, n: int) -> ZZSeries:
        """Truncate a power series to the first n terms."""
        ...


@runtime_checkable
class PowerSeriesRingQQ(Protocol):
    """Protocol for power series ring over rational field.

    This protocol defines the interface that all implementations of power series
    rings over QQ must follow. Concrete implementations include PythonPowerSeriesRingQQ
    and FlintPowerSeriesRingQQ.

    Attributes
    ----------
    domain : Domain
        The ground domain (QQ).
    prec : int
        The default precision for power series operations.
    """

    domain: Domain
    prec: int

    def __init__(self, prec: int = 6) -> None:
        """Initialize a power series ring over QQ."""
        ...

    def __repr__(self) -> str:
        """Return string representation of the ring."""
        ...

    @property
    def one(self) -> QQSeries:
        """Return the multiplicative identity (1) as a power series."""
        ...

    @property
    def zero(self) -> QQSeries:
        """Return the additive identity (0) as a power series."""
        ...

    @property
    def gen(self) -> QQSeries:
        """Return the generator (x) as a power series."""
        ...

    def pretty(self, s: QQSeries) -> str:
        """Return a pretty string representation of a power series."""
        ...

    def print(self, s: QQSeries) -> str:
        """Return a printable string representation of a power series."""
        ...

    def equal(self, s1: QQSeries, s2: QQSeries) -> bool:
        """Check if two power series are equal."""
        ...

    def negative(self, s: QQSeries) -> QQSeries:
        """Return the additive inverse of a power series."""
        ...

    def add(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Add two power series."""
        ...

    def subtract(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Subtract two power series."""
        ...

    def multiply(self, s1: QQSeries, s2: QQSeries) -> QQSeries:
        """Multiply two power series."""
        ...

    def multiply_ground(self, s: QQSeries, c: Any) -> QQSeries:
        """Multiply a power series by a ground element."""
        ...

    def trunc(self, s: QQSeries, n: int) -> QQSeries:
        """Truncate a power series to the first n terms."""
        ...
