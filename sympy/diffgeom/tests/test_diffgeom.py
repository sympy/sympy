from sympy.diffgeom.Rn import R2, R2_p, R2_r, R3, R3_r, R3_c, R3_s
from sympy.diffgeom import (intcurve_series, intcurve_diffequ, Differential,
        TensorProduct, WedgeProduct)
from sympy.core import symbols, Function, Derivative
from sympy.simplify import trigsimp, simplify
from sympy.functions import sqrt, atan2, sin, cos
from sympy.matrices import Matrix, eye


def test_R2():
    x0, y0, r0, theta0 = symbols('x0, y0, r0, theta0', real=True)
    point_r = R2_r.point([x0, y0])
    point_p = R2_p.point([r0, theta0])

    # r**2 = x**2 + y**2
    assert (R2.r**2 - R2.x**2 - R2.y**2)(point_r) == 0
    assert trigsimp( (R2.r**2 - R2.x**2 - R2.y**2)(point_p) ) == 0
    assert trigsimp(R2.e_r(R2.x**2+R2.y**2)(point_p).doit()) == 2*r0

    # polar->rect->polar == Id
    a, b = symbols('a b', positive=True)
    m = Matrix([[a], [b]])
    #TODO assert m == R2_r.coord_transform_to(R2_p, R2_p.coord_transform_to(R2_r, [a, b])).applyfunc(simplify)
    assert m == R2_p.coord_transform_to(R2_r, R2_r.coord_transform_to(R2_p, m)).applyfunc(simplify)


def test_R3():
    a, b, c = symbols('a b c', positive=True)
    m = Matrix([[a], [b], [c]])
    assert m == R3_c.coord_transform_to(R3_r, R3_r.coord_transform_to(R3_c, m)).applyfunc(simplify)
    #TODO assert m == R3_r.coord_transform_to(R3_c, R3_c.coord_transform_to(R3_r, m)).applyfunc(simplify)
    assert m == R3_s.coord_transform_to(R3_r, R3_r.coord_transform_to(R3_s, m)).applyfunc(simplify)
    #TODO assert m == R3_r.coord_transform_to(R3_s, R3_s.coord_transform_to(R3_r, m)).applyfunc(simplify)
    assert m == R3_s.coord_transform_to(R3_c, R3_c.coord_transform_to(R3_s, m)).applyfunc(simplify)
    #TODO assert m == R3_c.coord_transform_to(R3_s, R3_s.coord_transform_to(R3_c, m)).applyfunc(simplify)


def test_intcurve_diffequ():
    t = symbols('t')
    start_point = R2_r.point([1, 0])
    vector_field = -R2.y*R2.e_x + R2.x*R2.e_y
    equations, init_cond = intcurve_diffequ(vector_field, t, start_point)
    assert str(equations) == '[f_1(t) + Derivative(f_0(t), t), -f_0(t) + Derivative(f_1(t), t)]'
    assert str(init_cond) == '[f_0(0) - 1, f_1(0)]'
    equations, init_cond = intcurve_diffequ(vector_field, t, start_point, R2_p)
    assert str(equations) == '[Derivative(f_0(t), t), Derivative(f_1(t), t) - 1]'
    assert str(init_cond) == '[f_0(0) - 1, f_1(0)]'


def test_products():
    assert TensorProduct(R2.dx, R2.dy)(R2.e_x, R2.e_y) == R2.dx(R2.e_x)*R2.dy(R2.e_y) == 1
    assert WedgeProduct(R2.dx, R2.dy)(R2.e_x, R2.e_y) == 1


def test_differential():
    xdy = R2.x*R2.dy
    dxdy = Differential(xdy)
    assert dxdy(R2.e_x, R2.e_y) == 1

