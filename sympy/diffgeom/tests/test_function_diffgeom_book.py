from sympy.diffgeom.Rn import R2, R2_p, R2_r
from sympy.diffgeom import intcurve_series, intcurve_diffequ, Differential
from sympy.core import symbols, Function, Derivative
from sympy.simplify import trigsimp, simplify
from sympy.functions import sqrt, atan2, sin, cos
from sympy.matrices import Matrix

# Most of the functionality is covered in the
# test_functional_diffgeom_ch* tests which are based on the
# example from the paper of Sussman and Wisdom.
# If they do not cover something, additional tests are added in other test
# functions.

def test_functional_diffgeom_ch2():
    # From "Functional Differential Geometry" as of 2011
    # by Sussman and Wisdom.
    x0, y0, r0, theta0 = symbols('x0, y0, r0, theta0', real=True)
    x, y = symbols('x, y', real=True)
    f = Function('f')

    assert (R2_p.point_to_coords(R2_r.point([x0, y0])) ==
                Matrix([sqrt(x0**2+y0**2), atan2(y0, x0)]))
    assert (R2_r.point_to_coords(R2_p.point([r0, theta0])) ==
                Matrix([r0*cos(theta0), r0*sin(theta0)]))

    assert R2_p.jacobian(R2_r, [r0, theta0]) == Matrix([[cos(theta0), -r0*sin(theta0)], [sin(theta0),  r0*cos(theta0)]])

    field = f(R2.x, R2.y)
    p1_in_rect = R2_r.point([x0, y0])
    p1_in_polar = R2_p.point([sqrt(x0**2 + y0**2), atan2(y0,x0)])
    assert field(p1_in_rect) == f(x0, y0)
    assert field(p1_in_polar) == f(x0, y0)

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


def test_functional_diffgeom_ch3():
    # From "Functional Differential Geometry" as of 2011
    # by Sussman and Wisdom.
    x0, y0 = symbols('x0, y0', real=True)
    x, y, t = symbols('x, y, t', real=True)
    f = Function('f')
    b1 = Function('b1')
    b2 = Function('b2')
    p_r = R2_r.point([x0, y0])

    s_field = f(R2.x,R2.y)
    v_field = b1(R2.x)*R2.e_x + b2(R2.y)*R2.e_y
    assert v_field(s_field)(p_r).doit() ==  b1(x0)*Derivative(f(x0, y0), x0) + b2(y0)*Derivative(f(x0, y0), y0)

    assert R2.e_x(R2.r**2)(p_r) == 2*x0
    v = R2.e_x + 2*R2.e_y
    s = R2.r**2 + 3*R2.x
    assert v(s)(p_r).doit() == 2*x0 + 4*y0 + 3

    circ = -R2.y*R2.e_x + R2.x*R2.e_y
    series = intcurve_series(circ, t, R2_r.point([1, 0]), coeffs=True)
    series_x, series_y = zip(*series)
    assert all([term == cos(t).taylor_term(i,t) for i, term in enumerate(series_x)])
    assert all([term == sin(t).taylor_term(i,t) for i, term in enumerate(series_y)])


def test_functional_diffgeom_ch4():
    # From "Functional Differential Geometry" as of 2011
    # by Sussman and Wisdom.
    x0, y0, theta0 = symbols('x0, y0, theta0', real=True)
    x, y, r, theta = symbols('x, y, r, theta', real=True)
    r0 = symbols('r0', positive=True)
    f = Function('f')
    b1 = Function('b1')
    b2 = Function('b2')
    p_r = R2_r.point([x0, y0])
    p_p = R2_p.point([r0, theta0])

    f_field = b1(R2.x,R2.y)*R2.dx + b2(R2.x,R2.y)*R2.dy
    assert f_field(R2.e_x)(p_r) == b1(x0, y0)
    assert f_field(R2.e_y)(p_r) == b2(x0, y0)

    s_field_r = f(R2.x,R2.y)
    df = Differential(s_field_r)
    assert df(R2.e_x)(p_r).doit() == Derivative(f(x0, y0), x0)
    assert df(R2.e_y)(p_r).doit() == Derivative(f(x0, y0), y0)

    s_field_p = f(R2.r,R2.theta)
    df = Differential(s_field_p)
    assert trigsimp(df(R2.e_x)(p_p).doit()) == cos(theta0)*Derivative(f(r0, theta0), r0) - sin(theta0)*Derivative(f(r0, theta0), theta0)/r0
    assert trigsimp(df(R2.e_y)(p_p).doit()) == sin(theta0)*Derivative(f(r0, theta0), r0) + cos(theta0)*Derivative(f(r0, theta0), theta0)/r0

    assert R2.dx(R2.e_x)(p_r) == 1
    assert R2.dx(R2.e_y)(p_r) == 0

    circ = -R2.y*R2.e_x + R2.x*R2.e_y
    assert R2.dx(circ)(p_r).doit() == -y0
    assert R2.dy(circ)(p_r) == x0
    assert R2.dr(circ)(p_r) == 0
    assert simplify(R2.dtheta(circ)(p_r)) == 1

    assert (circ - R2.e_theta)(s_field_r)(p_r) == 0
