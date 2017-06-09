from sympy.vector.coordsysrect import CoordSysCartesian
from sympy import cos
from sympy import symbols


x, y, z, q = symbols('x y z q')


def test_spherical_system():
    a = CoordSysCartesian('a', curv_coord_name='spherical')
    assert a.lame_coefficients() == (1, a.x, cos(a.y)*a.x)


def test_cartesian_system():
    a = CoordSysCartesian('a')
    assert a.lame_coefficients() == (1, 1, 1)
