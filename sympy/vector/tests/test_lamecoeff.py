from sympy.vector.coordsysrect import CoordSysCartesian
from sympy import sqrt, sin, cos, simplify
from sympy import symbols
from sympy.core.function import Derivative
from sympy.vector.functions import gradient, curl, express

x, y, z, q = symbols('x y z q')


def test_spherical_system():
    a = CoordSysCartesian('a', curv_coord_name='spherical')
    assert a.lame_coefficients() == (1, a.x, cos(a.y)*a.x)
    assert a.transformation_equations() == (sin(a.y)*cos(a.z)*a.x, sin(a.y)*sin(a.z)*a.x, cos(a.y)*a.x)


def test_cartesian_system():
    a = CoordSysCartesian('a')
    assert a.lame_coefficients() == (1, 1, 1)
    assert a.transformation_equations() == (1, 1, 1)
