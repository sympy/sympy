from __future__ import annotations
from sympy.external.mpmath import repr_dps


def test_repr_dps_is_stable():
    assert repr_dps(13) == 6
    assert repr_dps(33) == 12
    assert repr_dps(53) == 17
    assert repr_dps(66) == 22
