from sympy import sin, cos, pi, symbols
from sympy.vector.coordsysrect import CoordSys3D
from sympy.vector.parametricregion import ParametricRegion

C = CoordSys3D('C')
r, theta, phi, a, b, t = symbols('r theta phi a b t')

def test_parametricregion():

    point = ParametricRegion(3, 4)
    assert point.definition == [3, 4]
    assert point.parameters == []
    assert point.bounds == []
    assert point.system == None

    l = ParametricRegion(C.y, (C.y, -3, 3), system=C)
    assert l.definition == [C.y]
    assert l.parameters == [C.y]
    assert l.system == C

    p1 = ParametricRegion(9*a, -16*b, (a, 0, 3))
    assert p1.definition == [9*a, -16*b]
    assert p1.parameters == [a]
    assert p1.bounds == [(a, 0, 3)]

    circle = ParametricRegion(r*cos(theta), r*sin(theta), (theta, 0, 2*pi))
    assert circle.definition == [r*cos(theta), r*sin(theta)]
    assert circle.system == None

    halfdisc = ParametricRegion(r*cos(theta), r* sin(theta), (r, -2, 2), (theta, 0, pi))
    assert halfdisc.definition == [r*cos(theta), r*sin(theta)]
    assert halfdisc.parameters == [r, theta]
    assert halfdisc.bounds == [(r, -2, 2), (theta, 0, pi)]

    parabola = ParametricRegion(a*t**2, 2*a*t, (t, -2, 2))
    assert parabola.definition == [a*t**2, 2*a*t]
    assert parabola.parameters == [t]

    ellipse = ParametricRegion(a*cos(t), b*sin(t), (t, 0, 8))
    assert ellipse.parameters == [t]
    assert ellipse.bounds == [(t, 0, 8)]

    sphere = ParametricRegion(r*sin(phi)*cos(theta),r*sin(phi)*sin(theta), r*cos(phi), (theta, 0, 2*pi), (phi, 0, pi))
    assert sphere.definition == [r*sin(phi)*cos(theta),r*sin(phi)*sin(theta), r*cos(phi)]
    assert sphere.parameters == [theta, phi]
    assert sphere.system == None

    cylinder = ParametricRegion(r*cos(theta), r*sin(theta), C.z, (theta, 0, 2*pi), (C.z, 0, 4), system=C)
    assert cylinder.parameters == [theta, C.z]
    assert cylinder.system == C
