from __future__ import annotations

from sympy.external.gmpy import GROUND_TYPES, MPQ, MPZ

from sympy.polys.domains.domain import Domain, Er
from sympy.polys.domains.field import Field, Ef
from sympy.polys.series.base import PowerSeriesRingProto, PowerSeriesRingFieldProto
from sympy.polys.series.ringpython import (
    USeries,
    PythonPowerSeriesRingZZ,
    PythonPowerSeriesRingQQ,
)

from typing import Generic, TYPE_CHECKING, cast, overload


if GROUND_TYPES == "flint":
    from sympy.polys.series.ringflint import (
        ZZSeries,
        QQSeries,
        FlintPowerSeriesRingZZ,
        FlintPowerSeriesRingQQ,
    )

    flint = True
else:
    if TYPE_CHECKING:
        from sympy.polys.series.ringflint import (
            ZZSeries,
            QQSeries,
            FlintPowerSeriesRingZZ,
            FlintPowerSeriesRingQQ,
        )
    else:
        ZZSeries = QQSeries = None
        FlintPowerSeriesRingZZ = FlintPowerSeriesRingQQ = None
    flint = False


class TSeriesElement(Generic[Er]):
    """Dummy type for lower power series elements, used only for static type checking."""

    # Dummy method to ensure that TSeriesElement is invariant in Er.
    def __invalid__(self, other: Er) -> Er:
        raise NotImplementedError

    flint = True


# Types for lower ring power series ring combining different ground types.
TSeriesRing = PowerSeriesRingProto[TSeriesElement[Er], Er]
TSeriesRingField = PowerSeriesRingFieldProto[TSeriesElement[Ef], Ef]


def PowerSeriesRingZZ(prec: int = 6) -> TSeriesRing[MPZ]:
    if FlintPowerSeriesRingZZ is None:
        R_python: PowerSeriesRingProto[USeries[MPZ], MPZ] = PythonPowerSeriesRingZZ(
            prec
        )
        return cast("TSeriesRing[MPZ]", R_python)
    else:
        R_flint: PowerSeriesRingProto[ZZSeries, MPZ] = FlintPowerSeriesRingZZ(prec)
        return cast("TSeriesRing[MPZ]", R_flint)


def PowerSeriesRingQQ(prec: int = 6) -> TSeriesRingField[MPQ]:
    if FlintPowerSeriesRingQQ is None:
        R_python: PowerSeriesRingProto[USeries[MPQ], MPQ] = PythonPowerSeriesRingQQ(
            prec
        )
        return cast("TSeriesRingField[MPQ]", R_python)
    else:
        R_flint: PowerSeriesRingProto[QQSeries, MPQ] = FlintPowerSeriesRingQQ(prec)
        return cast("TSeriesRingField[MPQ]", R_flint)


@overload
def _power_series_ring(K: Field[Ef], prec: int = 6) -> TSeriesRingField[Ef]: ...


@overload
def _power_series_ring(K: Domain[Er], prec: int = 6) -> TSeriesRing[Er]: ...


def _power_series_ring(
    K: Domain[Er] | Field[Ef], prec: int = 6
) -> TSeriesRing[Er] | TSeriesRingField[Ef]:
    """
    Helper function for the Power Series Ring classes to create a base ring from lower
    power series rings.
    """
    if K.is_ZZ:
        R_ZZ: TSeriesRing[MPZ] = PowerSeriesRingZZ(prec)
        return cast("TSeriesRing[Er]", R_ZZ)
    elif K.is_QQ:
        R_QQ: TSeriesRingField[MPQ] = PowerSeriesRingQQ(prec)
        return cast("TSeriesRingField[Ef]", R_QQ)
    else:
        raise ValueError(f"Unsupported ground domain: {K}")
