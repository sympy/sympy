from sympy.differential_geometry.Rn import R2, R2_p, R2_r
from sympy.differential_geometry import ScalarField, VectorField
from sympy import (symbols, simplify, sqrt, atan2, Matrix, sin, cos, Function,
        Derivative)
# TODO explicit and well ordered imports


def test_functional_differential_geometry_ch2():
    # From "Functional Differential Geometry" as of 2011
    # by Sussman and Wisdom.
    x0, y0, r0, theta0 = symbols('x0, y0, r0, theta0', real=True)
    x, y, r, theta = symbols('x, y, r, theta', real=True)
    f = Function('f')

    assert (R2_p.point_to_coords(R2_r.point([x0, y0])) ==
                Matrix([sqrt(x0**2+y0**2), atan2(y0, x0)]))
    assert (R2_r.point_to_coords(R2_p.point([r0, theta0])) ==
                Matrix([r0*cos(theta0), r0*sin(theta0)]))
    #TODO jacobian page 12 - 32

    field = ScalarField(R2_r, [x, y], f(x, y))
    p1_in_rect = R2_r.point([x0, y0])
    p1_in_polar = R2_p.point([sqrt(x0**2 + y0**2), atan2(y0,x0)])
    assert field(p1_in_rect) == f(x0, y0)
    # TODO better simplification for the next one
    #print simplify(field(p1_in_polar))
    #assert simplify(field(p1_in_polar)) == f(x0, y0)

    p_r = R2_r.point([x0, y0])
    p_p = R2_p.point([r0, theta0])
    assert R2.x(p_r) == x0
    assert R2.x(p_p) == r0*cos(theta0)
    assert R2.r(p_p) == r0
    assert R2.r(p_r) == sqrt(x0**2 + y0**2)
    assert R2.theta(p_r) == atan2(y0, x0)

    h = R2.x*R2.r**2 + R2.y**3
    assert h(p_r) == x0*(x0**2 + y0**2) + y0**3
    assert h(p_p) == r0**3*sin(theta0)**3 + r0**3*cos(theta0)


def test_functional_differential_geometry_ch3():
    # From "Functional Differential Geometry" as of 2011
    # by Sussman and Wisdom.
    x0, y0, r0, theta0 = symbols('x0, y0, r0, theta0', real=True)
    x, y, r, theta = symbols('x, y, r, theta', real=True)
    f = Function('f')
    b1 = Function('b1')
    b2 = Function('b2')
    p_r = R2_r.point([x0, y0])

    s_field = ScalarField(R2_r, [x, y], f(x,y))
    v_field = VectorField(R2_r, [x, y], [b1(x), b2(y)])
    assert v_field(s_field)(p_r) ==  b1(x0)*Derivative(f(x0, y0), x0) + b2(y0)*Derivative(f(x0, y0), y0)

    assert R2.d_dx(R2.r**2)(p_r) == 2*x0
    v = R2.d_dx + 2*R2.d_dy
    s = R2.r**2 + 3*R2.x
    assert v(s)(p_r) == 2*x0 + 4*y0 + 3


def test_R2():
    x0, y0, r0, theta0 = symbols('x0, y0, r0, theta0', real=True)
    point_r = R2_r.point([x0, y0])
    point_p = R2_p.point([r0, theta0])

    # r**2 = x**2 + y**2
    assert (R2.r**2 - R2.x**2 - R2.y**2)(point_r) == 0
    assert simplify( (R2.r**2 - R2.x**2 - R2.y**2)(point_p) ) == 0

    assert simplify( R2.d_dr(R2.x**2+R2.y**2)(point_p) ) == 2*r0
