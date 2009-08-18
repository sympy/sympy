from sympy.polys.specialpolys import (
    swinnerton_dyer_poly,
    fateman_poly_F_1, dmp_fateman_poly_F_1,
    fateman_poly_F_2, dmp_fateman_poly_F_2,
    fateman_poly_F_3, dmp_fateman_poly_F_3,
)

from sympy.polys.algebratools import ZZ

from sympy import symbols, Poly

x, y = symbols('x,y')

def test_swinnerton_dyer_poly():
    assert swinnerton_dyer_poly(1, x) == \
        Poly(x**2 - 2)

    assert swinnerton_dyer_poly(2, x) == \
        Poly(x**4 - 10*x**2 + 1)

    assert swinnerton_dyer_poly(3, x) == \
        Poly(x**8 - 40*x**6 + 352*x**4 - 960*x**2 + 576)

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

