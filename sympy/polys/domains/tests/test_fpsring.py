from sympy.core.symbol import symbols
from sympy.polys.domains import QQ, ZZ, RR
from sympy.polys.domains.fpsring import PowerSeriesRing
from sympy.polys.fps_ring import PowerSeriesElement, PowerSeriesPolyRing
from sympy.utilities.pytest import raises


def test_fpsring():
    R = PowerSeriesRing('QQ', 'x', 5)
    assert R.dtype is PowerSeriesElement
    assert type(R.gens[0]) == R.dtype
    assert str(R) == 'QQ[[x], 5]'
    assert str(R.gens[0]) == 'x + O(x**5)'

    R = PowerSeriesRing('RR', 'x')
    assert R.dtype is PowerSeriesElement
    assert type(R.gens[0]) == R.dtype
    assert str(R) == 'RR[[x], 6]'
