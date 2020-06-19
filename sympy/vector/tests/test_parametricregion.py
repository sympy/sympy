from sympy import sin, cos, pi
from sympy.vector.coordsysrect import CoordSys3D
from sympy.vector.parametricregion import ParametricRegion
from sympy.testing.pytest import raises
from sympy.abc import a, b, r, t, x, y, z, theta, phi

C = CoordSys3D('C')

def test_parametricregion():

    point = ParametricRegion((3, 4))
    assert point.definition == (3, 4)
    assert point.parameters == ()
    assert point.limits == {}
    assert point.dimensions == 0

    # line x = y
    line_xy = ParametricRegion((y, y), (y, 1, 5))
    assert line_xy .definition == (y, y)
    assert line_xy.parameters == (y,)
    assert line_xy.dimensions == 1

    # line y = z
    line_yz = ParametricRegion((x,t,t), x, (t, 1, 2))
    assert line_yz.definition == (x,t,t)
    assert line_yz.parameters == (x, t)
    assert line_yz.limits == {t: (1, 2)}
    assert line_yz.dimensions == 1

    p1 = ParametricRegion((9*a, -16*b), (a, 0, 2), (b, -1, 5))
    assert p1.definition == (9*a, -16*b)
    assert p1.parameters == (a, b)
    assert p1.limits == {a: (0, 2), b: (-1, 5)}
    assert p1.dimensions == 2

    p2 = ParametricRegion((t, t**3), t)
    assert p2.parameters == (t,)
    assert p2.limits == {}
    assert p2.dimensions == 0

    circle = ParametricRegion((r*cos(theta), r*sin(theta)), r, (theta, 0, 2*pi))
    assert circle.definition == (r*cos(theta), r*sin(theta))
    assert circle.dimensions == 1

    halfdisc = ParametricRegion((r*cos(theta), r* sin(theta)), (r, -2, 2), (theta, 0, pi))
    assert halfdisc.definition == (r*cos(theta), r*sin(theta))
    assert halfdisc.parameters == (r, theta)
    assert halfdisc.limits == {r: (-2, 2), theta: (0, pi)}
    assert halfdisc.dimensions == 2

    ellipse = ParametricRegion((a*cos(t), b*sin(t)), (t, 0, 8))
    assert ellipse.parameters == (t,)
    assert ellipse.limits == {t: (0, 8)}
    assert ellipse.dimensions == 1

    cylinder = ParametricRegion((r*cos(theta), r*sin(theta), z), (r, 0, 1), (theta, 0, 2*pi), (z, 0, 4))
    assert cylinder.parameters == (r, theta, z)
    assert cylinder.dimensions == 3

    sphere = ParametricRegion((r*sin(phi)*cos(theta),r*sin(phi)*sin(theta), r*cos(phi)),
                                r, (theta, 0, 2*pi), (phi, 0, pi))
    assert sphere.definition == (r*sin(phi)*cos(theta),r*sin(phi)*sin(theta), r*cos(phi))
    assert sphere.parameters == (r, theta, phi)
    assert sphere.dimensions == 2

    raises(ValueError, lambda: ParametricRegion((a*t**2, 2*a*t), (a, -2)))
    raises(ValueError, lambda: ParametricRegion((a, b), (a**2, sin(b)), (a, 2, 4, 6)))
