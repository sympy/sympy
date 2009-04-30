from sympy.polys.specialpolys import (
    fateman_poly_F_1, dmp_fateman_poly_F_1,
    fateman_poly_F_2, dmp_fateman_poly_F_2,
    fateman_poly_F_3, dmp_fateman_poly_F_3,
)

from sympy.polys.algebratools import ZZ

def test_fateman_poly_F_1():
    f,g,h = fateman_poly_F_1(1)
    F,G,H = dmp_fateman_poly_F_1(1, ZZ)

    assert [ t.rep.rep for t in [f,g,h] ] == [F,G,H]

    f,g,h = fateman_poly_F_1(3)
    F,G,H = dmp_fateman_poly_F_1(3, ZZ)

    assert [ t.rep.rep for t in [f,g,h] ] == [F,G,H]

def test_fateman_poly_F_2():
    f,g,h = fateman_poly_F_2(1)
    F,G,H = dmp_fateman_poly_F_2(1, ZZ)

    assert [ t.rep.rep for t in [f,g,h] ] == [F,G,H]

    f,g,h = fateman_poly_F_2(3)
    F,G,H = dmp_fateman_poly_F_2(3, ZZ)

    assert [ t.rep.rep for t in [f,g,h] ] == [F,G,H]

def test_fateman_poly_F_3():
    f,g,h = fateman_poly_F_3(1)
    F,G,H = dmp_fateman_poly_F_3(1, ZZ)

    assert [ t.rep.rep for t in [f,g,h] ] == [F,G,H]

    f,g,h = fateman_poly_F_3(3)
    F,G,H = dmp_fateman_poly_F_3(3, ZZ)

    assert [ t.rep.rep for t in [f,g,h] ] == [F,G,H]

