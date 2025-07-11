from typing import Optional, Type
from sympy.external.gmpy import GROUND_TYPES
from sympy.polys.domains import Domain
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


def power_series_ring(domain: Domain, prec: int = 6) -> PowerSeriesRing:
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
    >>> from sympy.polys.series.seriesring import power_series_ring
    >>> R_ZZ = power_series_ring(ZZ, 5)
    >>> R_QQ = power_series_ring(QQ, 10)
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
