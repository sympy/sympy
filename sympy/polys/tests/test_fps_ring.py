from sympy.abc import x
from sympy.polys.domains import QQ, ZZ, RR
from sympy.polys.fps_ring import PowerSeriesElement, PowerSeriesPolyRing
from sympy.polys.polyerrors import GeneratorsError
from sympy.testing.pytest import raises


def test_PowerSeriesPolyRing():
    R = PowerSeriesPolyRing('x', QQ)
    R2 =  PowerSeriesPolyRing([x], QQ)
    p1, p2 = R.gens[0], R2.gens[0]
    assert isinstance(R, PowerSeriesPolyRing)
    assert hash(R) == hash(R2)
    assert R == R2
    assert p1 == p2
    assert R != PowerSeriesPolyRing('x', ZZ)
    assert R != PowerSeriesPolyRing('x', RR)
    assert R != PowerSeriesPolyRing('y', QQ)
    assert str(R) == 'Power Series Ring in x over QQ with lex order'


    raises(GeneratorsError, lambda: PowerSeriesPolyRing('', ZZ, 3))
    raises(GeneratorsError, lambda: PowerSeriesPolyRing('x, y', QQ))

    raises(ValueError, lambda: PowerSeriesPolyRing('x', QQ, -1))
    raises(ValueError, lambda: PowerSeriesPolyRing('x', QQ, 0))


def test_PowerSeriesPolyRing_from_poly():
    from sympy.polys.rings import ring
    R, x = ring('x', QQ)

    R2 = PowerSeriesPolyRing('y', QQ)
    raises(GeneratorsError, lambda: R2.from_poly(x))


def test_PowerSeriesElement():
    R = PowerSeriesPolyRing('x', QQ)
    p1 = R.from_poly(R.ring.gens[0])
    p2 = PowerSeriesElement(R, R.ring.gens[0])
    assert isinstance(p1, PowerSeriesElement)
    assert isinstance(p2, PowerSeriesElement)
    assert p1 == p2
    assert str(p1) == 'x + O(x**6)'
    assert p1.ring == R
    assert p2.ring == R
