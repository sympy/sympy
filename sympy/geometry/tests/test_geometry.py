from __future__ import division

from sympy import S, Symbol, oo, pi, sqrt, symbols
from sympy.geometry import Circle, Curve, Ellipse, Line, Point, Polygon, Ray, RegularPolygon, Segment, Triangle
from sympy.geometry.entity import scale
from sympy.utilities.randtest import verify_numerically
from sympy.utilities.pytest import raises

x = Symbol('x', real=True)
y = Symbol('y', real=True)
z = Symbol('z', real=True)
t = Symbol('t', real=True)
k = Symbol('k', real=True)


def test_subs():
    p = Point(x, 2)
    q = Point(1, 1)
    r = Point(3, 4)
    for o in [p,
              Segment(p, q),
              Ray(p, q),
              Line(p, q),
              Triangle(p, q, r),
              RegularPolygon(p, 3, 6),
              Polygon(p, q, r, Point(5, 4)),
              Circle(p, 3),
              Ellipse(p, 3, 4)]:
        assert 'y' in str(o.subs(x, y))
    assert p.subs({x: 1}) == Point(1, 2)
    assert Point(1, 2).subs(Point(1, 2), Point(3, 4)) == Point(3, 4)
    assert Point(1, 2).subs((1, 2), Point(3, 4)) == Point(3, 4)
    assert Point(1, 2).subs(Point(1, 2), Point(3, 4)) == Point(3, 4)
    assert Point(1, 2).subs(set([(1, 2)])) == Point(2, 2)
    raises(ValueError, lambda: Point(1, 2).subs(1))
    raises(ValueError, lambda: Point(1, 1).subs((Point(1, 1), Point(1,
           2)), 1, 2))


def test_free_symbols():
    a, b, c, d, e, f, s = symbols('a:f,s')
    assert Point(a, b).free_symbols == set([a, b])
    assert Line((a, b), (c, d)).free_symbols == set([a, b, c, d])
    assert Ray((a, b), (c, d)).free_symbols == set([a, b, c, d])
    assert Ray((a, b), angle=c).free_symbols == set([a, b, c])
    assert Segment((a, b), (c, d)).free_symbols == set([a, b, c, d])
    assert Line((a, b), slope=c).free_symbols == set([a, b, c])
    assert Curve((a*s, b*s), (s, c, d)).free_symbols == set([a, b, c, d])
    assert Ellipse((a, b), c, d).free_symbols == set([a, b, c, d])
    assert Ellipse((a, b), c, eccentricity=d).free_symbols == \
        set([a, b, c, d])
    assert Ellipse((a, b), vradius=c, eccentricity=d).free_symbols == \
        set([a, b, c, d])
    assert Circle((a, b), c).free_symbols == set([a, b, c])
    assert Circle((a, b), (c, d), (e, f)).free_symbols == \
        set([e, d, c, b, f, a])
    assert Polygon((a, b), (c, d), (e, f)).free_symbols == \
        set([e, b, d, f, a, c])
    assert RegularPolygon((a, b), c, d, e).free_symbols == set([e, a, b, c, d])


def test_geometry_transforms():
    from sympy import Tuple
    c = Curve((x, x**2), (x, 0, 1))
    pts = [Point(0, 0), Point(1/2, 1/4), Point(1, 1)]
    cout = Curve((2*x - 4, 3*x**2 - 10), (x, 0, 1))
    pts_out = [Point(-4, -10), Point(-3, -37/4), Point(-2, -7)]
    assert c.scale(2, 3, (4, 5)) == cout
    assert [c.subs(x, xi/2) for xi in Tuple(0, 1, 2)] == pts
    assert [cout.subs(x, xi/2) for xi in Tuple(0, 1, 2)] == pts_out
    assert Triangle(*pts).scale(2, 3, (4, 5)) == Triangle(*pts_out)

    assert Ellipse((0, 0), 2, 3).scale(2, 3, (4, 5)) == \
        Ellipse(Point(-4, -10), 4, 9)
    assert Circle((0, 0), 2).scale(2, 3, (4, 5)) == \
        Ellipse(Point(-4, -10), 4, 6)
    assert Ellipse((0, 0), 2, 3).scale(3, 3, (4, 5)) == \
        Ellipse(Point(-8, -10), 6, 9)
    assert Circle((0, 0), 2).scale(3, 3, (4, 5)) == \
        Circle(Point(-8, -10), 6)
    assert Circle(Point(-8, -10), 6).scale(1/3, 1/3, (4, 5)) == \
        Circle((0, 0), 2)
    assert Curve((x + y, 3*x), (x, 0, 1)).subs(y, S.Half) == \
        Curve((x + 1/2, 3*x), (x, 0, 1))
    assert Curve((x, 3*x), (x, 0, 1)).translate(4, 5) == \
        Curve((x + 4, 3*x + 5), (x, 0, 1))
    assert Circle((0, 0), 2).translate(4, 5) == \
        Circle((4, 5), 2)
    assert Circle((0, 0), 2).scale(3, 3) == \
        Circle((0, 0), 6)
    assert Point(1, 1).scale(2, 3, (4, 5)) == \
        Point(-2, -7)
    assert Point(1, 1).translate(4, 5) == \
        Point(5, 6)
    assert scale(1, 2, (3, 4)).tolist() == \
        [[1, 0, 0], [0, 2, 0], [0, -4, 1]]
    assert RegularPolygon((0, 0), 1, 4).scale(2, 3, (4, 5)) == \
        Polygon(Point(-2, -10), Point(-4, -7), Point(-6, -10), Point(-4, -13))


def test_reflect():
    b = Symbol('b')
    m = Symbol('m')
    l = Line((0, b), slope=m)
    p = Point(x, y)
    r = p.reflect(l)
    dp = l.perpendicular_segment(p).length
    dr = l.perpendicular_segment(r).length
    assert verify_numerically(dp, dr)
    t = Triangle((0, 0), (1, 0), (2, 3))
    assert t.area == -t.reflect(l).area
    e = Ellipse((1, 0), 1, 2)
    assert e.area == -e.reflect(Line((1, 0), slope=0)).area
    assert e.area == -e.reflect(Line((1, 0), slope=oo)).area
    raises(NotImplementedError, lambda: e.reflect(Line((1, 0), slope=m)))
    assert Polygon((1, 0), (2, 0), (2, 2)).reflect(Line((3, 0), slope=oo)) \
        == Triangle(Point(5, 0), Point(4, 0), Point(4, 2))
    assert Polygon((1, 0), (2, 0), (2, 2)).reflect(Line((0, 3), slope=oo)) \
        == Triangle(Point(-1, 0), Point(-2, 0), Point(-2, 2))
    assert Polygon((1, 0), (2, 0), (2, 2)).reflect(Line((0, 3), slope=0)) \
        == Triangle(Point(1, 6), Point(2, 6), Point(2, 4))
    assert Polygon((1, 0), (2, 0), (2, 2)).reflect(Line((3, 0), slope=0)) \
        == Triangle(Point(1, 0), Point(2, 0), Point(2, -2))

    # test entity overrides
    c = Circle((x, y), 3)
    cr = c.reflect(l)
    assert cr == Circle(r, -3)
    assert c.area == -cr.area
    pent = RegularPolygon((1, 2), 1, 5)
    l = Line((0, pi), slope=sqrt(2))
    rpent = pent.reflect(l)
    poly_pent = Polygon(*pent.vertices)
    assert rpent.center == pent.center.reflect(l)
    assert str([w.n(3) for w in rpent.vertices]) == (
        '[Point2D(-0.586, 4.27), Point2D(-1.69, 4.66), '
        'Point2D(-2.41, 3.73), Point2D(-1.74, 2.76), '
        'Point2D(-0.616, 3.10)]')
    assert pent.area.equals(-rpent.area)
