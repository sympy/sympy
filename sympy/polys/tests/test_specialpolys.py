from sympy.polys.specialpolys import (
    fateman_poly_F_1, zzX_fateman_poly_F_1,
    fateman_poly_F_2, zzX_fateman_poly_F_2,
    fateman_poly_F_3, zzX_fateman_poly_F_3)

from sympy.polys.integerpolys import (
    zzX_from_poly)

def test_fateman_poly_F_1():
    f,g,h = fateman_poly_F_1(1)
    F,G,H = zzX_fateman_poly_F_1(1)

    assert map(zzX_from_poly, [f,g,h]) == [F,G,H]

    f,g,h = fateman_poly_F_1(3)
    F,G,H = zzX_fateman_poly_F_1(3)

    assert map(zzX_from_poly, [f,g,h]) == [F,G,H]

def test_fateman_poly_F_2():
    f,g,h = fateman_poly_F_2(1)
    F,G,H = zzX_fateman_poly_F_2(1)

    assert map(zzX_from_poly, [f,g,h]) == [F,G,H]

    f,g,h = fateman_poly_F_2(3)
    F,G,H = zzX_fateman_poly_F_2(3)

    assert map(zzX_from_poly, [f,g,h]) == [F,G,H]

def test_fateman_poly_F_3():
    f,g,h = fateman_poly_F_3(1)
    F,G,H = zzX_fateman_poly_F_3(1)

    assert map(zzX_from_poly, [f,g,h]) == [F,G,H]

    f,g,h = fateman_poly_F_3(3)
    F,G,H = zzX_fateman_poly_F_3(3)

    assert map(zzX_from_poly, [f,g,h]) == [F,G,H]
