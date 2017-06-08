from sympy.vector.coordsysrect import CoordSysCartesian
from sympy import sqrt, sin, cos
from sympy import symbols
from sympy.core.function import Derivative


x, y, z = symbols('x y z')


def test_spherical_system():
    a = CoordSysCartesian('a', 'spherical')
    assert a.lame_coefficients() == (1, a.x, cos(a.y)*a.x)
    print(a.transformation_equations())
    assert a.transformation_equations() == (sin(a.y)*cos(a.z)*a.x, sin(a.y)*sin(a.z)*a.x, cos(a.y)*a.x)


def test_cartesian_system():
    a = CoordSysCartesian('a')
    assert a.lame_coefficients() == (1, 1, 1)
    assert a.transformation_equations() == (1, 1, 1)


def test_any_system():
    a = CoordSysCartesian('a', transformation_equations=(x * sin(y) * cos(z), x * sin(y) * sin(z), x * cos(y)))
    assert a.lame_coefficients() == (sqrt(Derivative(cos(a.y)*a.x, a.x)**2 + Derivative(sin(a.y)*sin(a.z)*a.x, a.x)**2
                                          + Derivative(sin(a.y)*cos(a.z)*a.x, a.x)**2),
                                     sqrt(Derivative(cos(a.y)*a.x, a.y)**2 + Derivative(sin(a.y)*sin(a.z)*a.x, a.y)**2
                                          + Derivative(sin(a.y)*cos(a.z)*a.x, a.y)**2),
                                     sqrt(Derivative(cos(a.y)*a.x, a.z)**2 + Derivative(sin(a.y)*sin(a.z)*a.x, a.z)**2
                                          + Derivative(sin(a.y)*cos(a.z)*a.x, a.z)**2))
    assert a.transformation_equations() == (sin(a.y)*cos(a.z)*a.x, sin(a.y)*sin(a.z)*a.x, cos(a.y)*a.x)
