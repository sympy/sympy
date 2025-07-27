from __future__ import annotations

from typing import Any, TYPE_CHECKING
from sympy.external.gmpy import GROUND_TYPES, MPQ, MPZ
from sympy.polys.domains.domain import Domain, Er
from sympy.polys.series.ringpython import (
    PythonPowerSeriesRingZZ,
    PythonPowerSeriesRingQQ,
)

if TYPE_CHECKING:
    from sympy.polys.series import PowerSeriesRing

# Make use of flint implementation if flint available.
FlintPowerSeriesRingZZ: type[Any] | None = None
FlintPowerSeriesRingQQ: type[Any] | None = None

if GROUND_TYPES == "flint":
    from sympy.polys.series.ringflint import (
        FlintPowerSeriesRingZZ,
        FlintPowerSeriesRingQQ,
    )


def power_series_ring(domain: Domain[Er], prec: int = 6) -> PowerSeriesRing[Any, Er]:
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
    >>> from sympy import QQ
    >>> from sympy.polys.series import power_series_ring
    >>> R_QQ = power_series_ring(QQ, 5)
    >>> inv = R_QQ.inverse(R_QQ([1, 2, (3, 4), (5, 6)]))
    >>> R_QQ.print(inv)
    1 - 2*x + 13/4*x**2 - 35/6*x**3 + 523/48*x**4 + O(x**5)

    """
    if domain.is_ZZ:
        if FlintPowerSeriesRingZZ is not None:
            return FlintPowerSeriesRingZZ(prec)
        R_ZZ: PowerSeriesRing[Any, MPZ] = PythonPowerSeriesRingZZ(prec)
        return R_ZZ  # type: ignore

    if domain.is_QQ:
        if FlintPowerSeriesRingQQ is not None:
            return FlintPowerSeriesRingQQ(prec)
        R_QQ: PowerSeriesRing[Any, MPQ] = PythonPowerSeriesRingQQ(prec)
        return R_QQ  # type: ignore

    raise TypeError("Ground domain must be an instance of QQ or ZZ")
