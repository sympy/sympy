from __future__ import annotations

from typing import TYPE_CHECKING, Generic, cast
from sympy.external.gmpy import GROUND_TYPES, MPQ, MPZ
from sympy.polys.domains.domain import Domain, Er

from sympy.polys.series.base import PowerSeriesRingProto

from sympy.polys.series.ringpython import (
    USeries,
    PythonPowerSeriesRingZZ,
    PythonPowerSeriesRingQQ,
)

if TYPE_CHECKING:
    from sympy.polys.series.ringflint import (
        ZZSeries,
        QQSeries,
        FlintPowerSeriesRingZZ as _FlintPowerSeriesRingZZ,
        FlintPowerSeriesRingQQ as _FlintPowerSeriesRingQQ,
    )

    FlintPowerSeriesRingZZ: type[_FlintPowerSeriesRingZZ] | None
    FlintPowerSeriesRingQQ: type[_FlintPowerSeriesRingQQ] | None


if GROUND_TYPES == "flint":
    from sympy.polys.series.ringflint import (
        FlintPowerSeriesRingZZ,
        FlintPowerSeriesRingQQ,
    )
else:
    FlintPowerSeriesRingZZ = None
    FlintPowerSeriesRingQQ = None


class TSeries(Generic[Er]):
    """
    Dummy type for power series elements, used only for static type checking.
    Actual series elements are implementation-dependent and should be
    used via :class:`PowerSeriesRing`.
    """

    # Dummy method to ensure that TSeries is invariant in Er.
    def __invalid__(self, other: Er) -> Er:
        raise NotImplementedError


PowerSeriesRing = PowerSeriesRingProto[TSeries[Er], Er]


def PowerSeriesRingZZ(prec: int = 6) -> PowerSeriesRing[MPZ]:
    if FlintPowerSeriesRingZZ is None:
        R_python: PowerSeriesRingProto[USeries[MPZ], MPZ] = PythonPowerSeriesRingZZ(
            prec
        )
        return cast("PowerSeriesRing[MPZ]", R_python)
    else:
        R_flint: PowerSeriesRingProto[ZZSeries, MPZ] = FlintPowerSeriesRingZZ(prec)
        return cast("PowerSeriesRing[MPZ]", R_flint)


def PowerSeriesRingQQ(prec: int = 6) -> PowerSeriesRing[MPQ]:
    if FlintPowerSeriesRingQQ is None:
        R_python: PowerSeriesRingProto[USeries[MPQ], MPQ] = PythonPowerSeriesRingQQ(
            prec
        )
        return cast("PowerSeriesRing[MPQ]", R_python)
    else:
        R_flint: PowerSeriesRingProto[QQSeries, MPQ] = FlintPowerSeriesRingQQ(prec)
        return cast("PowerSeriesRing[MPQ]", R_flint)


def power_series_ring(K: Domain[Er], prec: int = 6) -> PowerSeriesRing[Er]:
    """
    Create a power series ring over the given domain.

    Parameters
    ==========

    domain : Domain
        The ground domain for the power series ring. Must be ZZ or QQ.
    prec : int, optional
        The default precision for power series operations. Default is 6.

    Returns
    =======

    PowerSeriesRingZZ | PowerSeriesRingQQ
        A power series ring instance. The actual implementation (Flint or Python)
        depends on whether python-flint is available.

    Examples
    ========

    >>> from sympy import QQ
    >>> from sympy.polys.series import power_series_ring
    >>> R_QQ = power_series_ring(QQ, 5)
    >>> inv = R_QQ.inverse(R_QQ([1, 2, (3, 4), (5, 6)]))
    >>> R_QQ.print(inv)
    """

    if K.is_ZZ:
        return cast("PowerSeriesRing[Er]", PowerSeriesRingZZ(prec))
    elif K.is_QQ:
        return cast("PowerSeriesRing[Er]", PowerSeriesRingQQ(prec))
    else:
        raise ValueError(f"Unsupported ground domain: {K}")
