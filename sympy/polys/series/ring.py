from __future__ import annotations

from typing import Callable, cast
from sympy.external.gmpy import GROUND_TYPES, MPQ, MPZ
from sympy.polys.domains.domain import Domain, Er
from sympy.polys.series.base import PowerSeriesRing


_PowerSeriesRingZZ: Callable[[int], PowerSeriesRing[MPZ]]
_PowerSeriesRingQQ: Callable[[int], PowerSeriesRing[MPQ]]

if GROUND_TYPES == "flint":
    # Make use of flint implementation if flint available.
    from sympy.polys.series.ringflint import (
        FlintPowerSeriesRingZZ as _PowerSeriesRingZZ,
        FlintPowerSeriesRingQQ as _PowerSeriesRingQQ,
    )
else:
    # Fall back to python implementation.
    from sympy.polys.series.ringpython import (
        PythonPowerSeriesRingZZ as _PowerSeriesRingZZ,
        PythonPowerSeriesRingQQ as _PowerSeriesRingQQ,
    )


def power_series_ring(domain: Domain[Er], prec: int = 6) -> PowerSeriesRing[Er]:
    if domain.is_ZZ:
        return cast(PowerSeriesRing[Er], _PowerSeriesRingZZ(prec))
    elif domain.is_QQ:
        return cast(PowerSeriesRing[Er], _PowerSeriesRingQQ(prec))
    else:
        raise TypeError("Ground domain must be an instance of QQ or ZZ")
