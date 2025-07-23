from typing import Optional, Type
from sympy.external.gmpy import GROUND_TYPES
from sympy.polys.domains.domain import Domain, Er
from sympy.polys.series.powerseriesring import PowerSeriesRing
from sympy.polys.series.python_powerseriesring import (
    PythonPowerSeriesRingZZ,
    PythonPowerSeriesRingQQ,
)

# Make use of flint implementation if flint available.
FlintPowerSeriesRingZZ: Optional[Type[PowerSeriesRing]] = None
FlintPowerSeriesRingQQ: Optional[Type[PowerSeriesRing]] = None

if GROUND_TYPES == "flint":
    from sympy.polys.series.flint_powerseriesring import (
        FlintPowerSeriesRingZZ,
        FlintPowerSeriesRingQQ,
    )


def power_series_ring(domain: Domain[Er], prec: int = 6) -> PowerSeriesRing:
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
        return PythonPowerSeriesRingZZ(prec)
    if domain.is_QQ:
        if FlintPowerSeriesRingQQ is not None:
            return FlintPowerSeriesRingQQ(prec)
        return PythonPowerSeriesRingQQ(prec)
    raise TypeError("Ground domain must be an instance of QQ or ZZ")
