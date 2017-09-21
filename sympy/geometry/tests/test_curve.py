from __future__ import division

from sympy import Symbol, pi, symbols, Tuple, S, sqrt, asinh, Matrix, sin, cos, sympify
from sympy.geometry import Curve, Line, Point, Ellipse, Ray, Segment, Circle, Polygon, RegularPolygon
from sympy.matrices.common import ShapeError
from sympy.utilities.pytest import raises, slow


def test_curve():
    x = Symbol('x', real=True)
    s = Symbol('s')
    z = Symbol('z')

    # this curve is independent of the indicated parameter
    c = Curve([2 * s, s ** 2], (z, 0, 2))

    assert c.parameter == z
    assert c.functions == Matrix((2 * s, s ** 2))
    assert c.arbitrary_point() == Point(2 * s, s ** 2)
    assert c.arbitrary_point(z) == Point(2 * s, s ** 2)

    # this is how it is normally used
    c = Curve([2 * s, s ** 2], (s, 0, 2))

    assert c.parameter == s
    assert c.functions == Matrix((2 * s, s ** 2))
    t = Symbol('t')
    # the t returned as assumptions
    assert c.arbitrary_point() != Point(2 * t, t ** 2)
    t = Symbol('t', real=True)
    # now t has the same assumptions so the test passes
    assert c.arbitrary_point() == Point(2 * t, t ** 2)
    assert c.arbitrary_point(z) == Point(2 * z, z ** 2)
    assert c.arbitrary_point(c.parameter) == Point(2 * s, s ** 2)
    assert c.arbitrary_point(None) == Point(2 * s, s ** 2)
    assert c.plot_interval() == [t, 0, 2]
    assert c.plot_interval(z) == [z, 0, 2]

    assert Curve([x, x], (x, 0, 1)).rotate(pi / 2, (1, 2)).scale(2, 3).translate(
        1, 3).arbitrary_point(s) == \
           Line((0, 0), (1, 1)).rotate(pi / 2, (1, 2)).scale(2, 3).translate(
               1, 3).arbitrary_point(s) == \
           Point(-2 * s + 7, 3 * s + 6)

    raises(ValueError, lambda: Curve((s), (s, 1, 2)))
    raises(ValueError, lambda: Curve((x, x * 2), (1, x)))

    raises(ValueError, lambda: Curve((s, s + t), (s, 1, 2)).arbitrary_point())
    raises(ValueError, lambda: Curve((s, s + t), (t, 1, 2)).arbitrary_point(s))


@slow
def test_free_symbols():
    a, b, c, d, e, f, s = symbols('a:f,s')
    assert Point(a, b).free_symbols == {a, b}
    assert Line((a, b), (c, d)).free_symbols == {a, b, c, d}
    assert Ray((a, b), (c, d)).free_symbols == {a, b, c, d}
    assert Ray((a, b), angle=c).free_symbols == {a, b, c}
    assert Segment((a, b), (c, d)).free_symbols == {a, b, c, d}
    assert Line((a, b), slope=c).free_symbols == {a, b, c}
    assert Curve((a * s, b * s), (s, c, d)).free_symbols == {a, b, c, d}
    assert Ellipse((a, b), c, d).free_symbols == {a, b, c, d}
    assert Ellipse((a, b), c, eccentricity=d).free_symbols == \
           {a, b, c, d}
    assert Ellipse((a, b), vradius=c, eccentricity=d).free_symbols == \
           {a, b, c, d}
    assert Circle((a, b), c).free_symbols == {a, b, c}
    assert Circle((a, b), (c, d), (e, f)).free_symbols == \
           {e, d, c, b, f, a}
    assert Polygon((a, b), (c, d), (e, f)).free_symbols == \
           {e, b, d, f, a, c}
    assert RegularPolygon((a, b), c, d, e).free_symbols == {e, a, b, c, d}


def test_transform():
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)
    c = Curve((x, x ** 2), (x, 0, 1))
    cout = Curve((2 * x - 4, 3 * x ** 2 - 10), (x, 0, 1))
    pts = [Point(0, 0), Point(1 / 2, 1 / 4), Point(1, 1)]
    pts_out = [Point(-4, -10), Point(-3, -37 / 4), Point(-2, -7)]

    assert c.scale(2, 3, (4, 5)) == cout
    assert [c.subs(x, xi / 2) for xi in Tuple(0, 1, 2)] == pts
    assert [cout.subs(x, xi / 2) for xi in Tuple(0, 1, 2)] == pts_out
    assert Curve((x + y, 3 * x), (x, 0, 1)).subs(y, S.Half) == \
           Curve((x + 1 / 2, 3 * x), (x, 0, 1))
    assert Curve((x, 3 * x), (x, 0, 1)).translate(4, 5) == \
           Curve((x + 4, 3 * x + 5), (x, 0, 1))


def test_length():
    t = Symbol('t', real=True)

    c1 = Curve((t, 0), (t, 0, 1))
    assert c1.length == 1

    c2 = Curve((t, t), (t, 0, 1))
    assert c2.length == sqrt(2)

    c3 = Curve((t ** 2, t), (t, 2, 5))
    assert c3.length == -sqrt(17) - asinh(4) / 4 + asinh(10) / 4 + 5 * sqrt(101) / 2


def test_tangent():
    t = Symbol('t', real=True)

    c1 = Curve((t,), (t, 0, 1))
    assert c1.tangent == Matrix([[1]])

    c2 = Curve((t, 3 * sin(t), 3 * cos(t), t ** 2), (t, 0, 1))
    assert c2.tangent == Matrix([[1], [3 * cos(t)], [-3 * sin(t)], [2 * t]])


def test_normal():
    t = Symbol('t', real=True)

    c1 = Curve((t,), (t, 0, 1))
    assert c1.normal == Matrix([[0]])

    c2 = Curve((t, t ** 2, cos(t), sin(t)), (t, 0, 1))
    assert c2.normal == Matrix([[0], [2], [-cos(t)], [-sin(t)]])


def test_binormal():
    t = Symbol('t', real=True)

    c1 = Curve((t, t, t), (t, 0, 1))
    assert c1.binormal == Matrix([[0], [0], [0]])

    c2 = Curve((1, t, t ** 2), (t, 0, 1))
    assert c2.binormal == Matrix([[2], [0], [0]])

    raises(ShapeError, lambda: Curve((t,), (t, 0, 1)).binormal)


def test_curvature():
    t = Symbol('t', real=True)

    c1 = Curve((t,), (t, 0, 1))
    assert c1.curvature == 0

    c2 = Curve((t, t ** 2), (t, 0, 1))
    assert c2.curvature == 2 / (4 * t ** 2 + 1) ** (sympify(3) / 2)


def test_torsion():
    a = Symbol('a', real=True)
    b = Symbol('b', real=True)
    t = Symbol('t', real=True)

    c1 = Curve((t, t ** 2, t ** 2), (t, 0, 1))
    assert c1.torsion == 0

    c2 = Curve((a * cos(t), a * sin(t), b * t), (t, 0, 1))
    assert c2.torsion == b / sqrt(a ** 2 + b ** 2)


def test_line_integral():
    t = Symbol('t', real=True)
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)
    z = Symbol('z', real=True)

    c1 = Curve((3 * cos(t), 3 * sin(t)), (t, -pi / 2, pi / 2))
    assert c1.line_integral(x * y ** 4, (x, y)) == sympify(1458) / 5

    c2 = Curve((t, t + 1, t + 2), (t, 0, 2))
    assert c2.line_integral(x + y + z, (x, y, z)) == 12 * sqrt(3)


def test_vector_line_integral():
    t = Symbol('t', real=True)
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)
    z = Symbol('z', real=True)

    c3 = Curve((t, t ** 2, t ** 3), (t, 0, 1))
    assert c3.vector_line_integral([8 * x ** 2 * y * z, 5 * z, -4 * x * y], (x, y, z)) == 1

    c4 = Curve((4 * t - 1, 2 - 2 * t, t), (t, 0, 1))
    assert c4.vector_line_integral((x * z, 0, -y * z), (x, y, z)) == 3
