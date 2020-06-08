from sympy import sin, cos, pi, symbols
from sympy.vector.coordsysrect import CoordSys3D
from sympy.vector.parametricregion import ParametricRegion
from sympy.testing.pytest import raises

C = CoordSys3D('C')
r, theta, phi, a, b, t = symbols('r theta phi a b t')

def test_parametricregion():

    point = ParametricRegion((), (3, 4), {})
    assert point.definition == (3, 4)
    assert point.parameters == ()
    assert point.limits == {}
    assert point.system == None

    # line x = y
    line_xy = ParametricRegion((C.y), (C.y), limits={C.y: (-3, 3)}, system=C)
    assert line_xy .definition == (C.y,)
    assert line_xy.parameters == (C.y,)
    assert line_xy.system == C

    # line y = z
    line_yz = ParametricRegion((t), (C.x,t,t), limits={t: (1, 2)}, system=C)
    assert line_yz.definition == (C.x,t,t)
    assert line_yz.parameters == (t,)

    p1 = ParametricRegion((a, b), (9*a, -16*b), limits={a: (0, 2), b: (-1, 5)})
    assert p1.definition == (9*a, -16*b)
    assert p1.parameters == (a, b)
    assert p1.limits == {a: (0, 2), b: (-1, 5)}

    p2 = ParametricRegion(t, (t, t**3))
    assert p2.parameters == (t,)
    assert p2.limits == {}

    circle = ParametricRegion((r, theta), (r*cos(theta), r*sin(theta)), {theta: (0, 2*pi)})
    assert circle.definition == (r*cos(theta), r*sin(theta))
    assert circle.system == None

    halfdisc = ParametricRegion((r, theta), (r*cos(theta), r* sin(theta)), {r: (-2, 2), theta: (0, pi)})
    assert halfdisc.definition == (r*cos(theta), r*sin(theta))
    assert halfdisc.parameters == (r, theta)
    assert halfdisc.limits == {r: (-2, 2), theta: (0, pi)}

    ellipse = ParametricRegion((t), (a*cos(t), b*sin(t)), {t: (0, 8)})
    assert ellipse.parameters == (t,)
    assert ellipse.limits == {t: (0, 8)}

    cylinder = ParametricRegion((r, theta, C.z), (cos(theta), r*sin(theta), C.z), {theta: (0, 2*pi), C.z: (0, 4)},
                                system=C)
    assert cylinder.system == C

    sphere = ParametricRegion((r, theta, phi), (r*sin(phi)*cos(theta),r*sin(phi)*sin(theta), r*cos(phi)),
                            {theta: ( 0, 2*pi), phi: ( 0, pi)})
    assert sphere.definition == (r*sin(phi)*cos(theta),r*sin(phi)*sin(theta), r*cos(phi))
    assert sphere.parameters == (r, theta, phi)
    assert sphere.system == None

    raises(ValueError, lambda: ParametricRegion((t), (a*t**2, 2*a*t), {a:  (-2, 2)}))
    raises(ValueError, lambda: ParametricRegion((a, b), (a**2, sin(b)), {a: (2, 4, 6)}))
