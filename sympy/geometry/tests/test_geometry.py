from sympy import Symbol, Rational, sqrt, pi, cos, oo, simplify, Real, raises
from sympy.geometry import Point, Polygon, convex_hull, Segment, \
        RegularPolygon, Circle, Ellipse, GeometryError, Line, intersection, \
        Ray, Triangle, are_similar, Curve

x = Symbol('x', real=True)
y = Symbol('y', real=True)
x1 = Symbol('x1', real=True)
x2 = Symbol('x2', real=True)
y1 = Symbol('y1', real=True)
y2 = Symbol('y2', real=True)
half = Rational(1,2)

def feq(a, b):
    """Test if two floating point values are 'equal'."""
    t = Real("1.0E-10")
    return -t < a-b < t

def test_curve():
    t = Symbol('t')
    z = Symbol('z')
    C = Curve([2*t, t**2], (z, 0, 2))

    assert C.parameter == z
    assert C.functions == [2*t, t**2]


def test_point():
    p1 = Point(x1, x2)
    p2 = Point(y1, y2)
    p3 = Point(0, 0)
    p4 = Point(1, 1)

    assert len(p1) == 2
    assert p2[1] == y2
    assert (p3+p4) == p4
    assert (p2-p1) == Point(y1-x1, y2-x2)
    assert p4*5 == Point(5, 5)
    assert -p2 == Point(-y1, -y2)

    assert Point.midpoint(p3, p4) == Point(half, half)
    assert Point.midpoint(p1, p4) == Point(half + half*x1, half + half*x2)
    assert Point.midpoint(p2, p2) == p2

    assert Point.distance(p3, p4) == sqrt(2)
    assert Point.distance(p1, p1) == 0
    #assert Point.distance(p3, p2) == abs(p2)

    p1_1 = Point(x1, x1)
    p1_2 = Point(y2, y2)
    p1_3 = Point(x1 + 1, x1)
    assert Point.is_collinear(p3)
    assert Point.is_collinear(p3, p4)
    assert Point.is_collinear(p3, p4, p1_1, p1_2)
    assert Point.is_collinear(p3, p4, p1_1, p1_3) == False

    x_pos = Symbol('x', real=True, positive=True)
    p2_1 = Point(x_pos, 0)
    p2_2 = Point(0, x_pos)
    p2_3 = Point(-x_pos, 0)
    p2_4 = Point(0, -x_pos)
    p2_5 = Point(x_pos, 5)
    assert Point.is_concyclic(p2_1)
    assert Point.is_concyclic(p2_1, p2_2)
    assert Point.is_concyclic(p2_1, p2_2, p2_3, p2_4)
    assert Point.is_concyclic(p2_1, p2_2, p2_3, p2_5) == False

def test_line():
    p1 = Point(0, 0)
    p2 = Point(1, 1)
    p3 = Point(x1, x1)
    p4 = Point(y1, y1)
    p5 = Point(x1, 1 + x1)

    l1 = Line(p1, p2)
    l2 = Line(p3, p4)
    l3 = Line(p3, p5)

    # Basic stuff
    assert Line(p1, p2) == Line(p2, p1)
    assert l1 == l2
    assert l1 != l3
    assert l1.slope == 1
    assert l3.slope == oo
    assert p1 in l1 # is p1 on the line l1?
    assert p1 not in l3

    assert simplify(l1.equation()) in (x-y, y-x)
    assert simplify(l3.equation()) in (x-x1, x1-x)

    assert l2.arbitrary_point() in l2
    for ind in xrange(0, 5):
        assert l3.random_point() in l3

    # Orthogonality
    p1_1 = Point(-x1, x1)
    l1_1 = Line(p1, p1_1)
    assert l1.perpendicular_line(p1) == l1_1
    assert Line.is_perpendicular(l1, l1_1)
    assert Line.is_perpendicular(l1 , l2) == False

    # Parallelity
    p2_1 = Point(-2*x1, 0)
    l2_1 = Line(p3, p5)
    assert l2.parallel_line(p1_1) == Line(p2_1, p1_1)
    assert l2_1.parallel_line(p1) == Line(p1, Point(0, 2))
    assert Line.is_parallel(l1, l2)
    assert Line.is_parallel(l2, l3) == False
    assert Line.is_parallel(l2, l2.parallel_line(p1_1))
    assert Line.is_parallel(l2_1, l2_1.parallel_line(p1))

    # Intersection
    assert intersection(l1, p1) == [p1]
    assert intersection(l1, p5) == []
    assert intersection(l1, l2) in [[l1], [l2]]
    assert intersection(l1, l1.parallel_line(p5)) == []

    # Concurrency
    l3_1 = Line(Point(5, x1), Point(-Rational(3,5), x1))
    assert Line.is_concurrent(l1, l3)
    assert Line.is_concurrent(l1, l3, l3_1)
    assert Line.is_concurrent(l1, l1_1, l3) == False

    # Projection
    assert l2.projection(p4) == p4
    assert l1.projection(p1_1) == p1
    assert l3.projection(p2) == Point(x1, 1)

    # Finding angles
    l1_1 = Line(p1, Point(5, 0))
    assert feq(Line.angle_between(l1, l1_1).evalf(), pi.evalf()/4)

    # Testing Rays and Segments (very similar to Lines)
    r1 = Ray(p1, Point(-1, 5))
    r2 = Ray(p1, Point(-1, 1))
    r3 = Ray(p3, p5)
    assert l1.projection(r1) == Ray(p1, p2)
    assert l1.projection(r2) == p1
    assert r3 != r1

    s1 = Segment(p1, p2)
    s2 = Segment(p1, p1_1)
    assert s1.midpoint == Point(Rational(1,2), Rational(1,2))
    assert s2.length == sqrt( 2*(x1**2) )
    assert s1.perpendicular_bisector() == Line(Point(0, 1), Point(1, 0))

    # Special cases of projection and intersection
    r1 = Ray(Point(1, 1), Point(2, 2))
    r2 = Ray(Point(2, 2), Point(0, 0))
    r3 = Ray(Point(1, 1), Point(-1, -1))
    r4 = Ray(Point(0, 4), Point(-1, -5))
    assert intersection(r1, r2) == [Segment(Point(1, 1), Point(2, 2))]
    assert intersection(r1, r3) == [Point(1, 1)]
    assert r1.projection(r3) == Point(1, 1)
    assert r1.projection(r4) == Segment(Point(1, 1), Point(2, 2))

    r5 = Ray(Point(0, 0), Point(0, 1))
    r6 = Ray(Point(0, 0), Point(0, 2))
    assert r5 in r6
    assert r6 in r5

    s1 = Segment(Point(0, 0), Point(2, 2))
    s2 = Segment(Point(-1, 5), Point(-5, -10))
    s3 = Segment(Point(0, 4), Point(-2, 2))
    assert intersection(r1, s1) == [Segment(Point(1, 1), Point(2, 2))]
    assert r1.projection(s2) == Segment(Point(1, 1), Point(2, 2))
    assert s3.projection(r1) == Segment(Point(0, 4), Point(-1, 3))

    l1 = Line(Point(0, 0), Point(3, 4))
    r1 = Ray(Point(0, 0), Point(3, 4))
    s1 = Segment(Point(0, 0), Point(3, 4))
    assert intersection(l1, l1) == [l1]
    assert intersection(l1, r1) == [r1]
    assert intersection(l1, s1) == [s1]
    assert intersection(r1, l1) == [r1]
    assert intersection(s1, l1) == [s1]

    entity1 = Segment(Point(-10,10), Point(10,10))
    entity2 = Segment(Point(-5,-5), Point(-5,5))
    assert intersection(entity1, entity2) == []


def test_ellipse():
    p1 = Point(0, 0)
    p2 = Point(1, 1)
    p3 = Point(x1, x2)
    p4 = Point(0, 1)
    p5 = Point(-1, 0)

    e1 = Ellipse(p1, 1, 1)
    e2 = Ellipse(p2, half, 1)
    e3 = Ellipse(p1, y1, y1)
    c1 = Circle(p1, 1)
    c2 = Circle(p2,1)
    c3 = Circle(Point(sqrt(2),sqrt(2)),1)

    # Test creation with three points
    cen,rad = Point(3*half, 2), 5*half
    assert Circle(Point(0,0), Point(3,0), Point(0,4)) == Circle(cen, rad)
    raises(GeometryError, "Circle(Point(0,0), Point(1,1), Point(2,2))")

    # Basic Stuff
    assert e1 == c1
    assert e1 != e2
    assert p4 in e1
    assert p2 not in e2
    assert e1.area == pi
    assert e2.area == pi/2
    assert e3.area == pi*(y1**2)
    assert c1.area == e1.area
    assert c1.circumference == 2*pi

    assert e2.arbitrary_point() in e2
    for ind in xrange(0, 5):
        assert e3.random_point() in e3

    # Foci
    f1,f2 = Point(sqrt(12), 0), Point(-sqrt(12), 0)
    ef = Ellipse(Point(0, 0), 4, 2)
    assert ef.foci in [(f1, f2), (f2, f1)]

    # Tangents
    v = sqrt(2) / 2
    p1_1 = Point(v, v)
    p1_2 = p2 + Point(half, 0)
    p1_3 = p2 + Point(0, 1)
    assert e1.tangent_line(p4) == c1.tangent_line(p4)
    assert e2.tangent_line(p1_2) == Line(p1_2, p2 + Point(half, 1))
    assert e2.tangent_line(p1_3) == Line(p1_3, p2 + Point(half, 1))
    assert c1.tangent_line(p1_1) == Line(p1_1, Point(0, sqrt(2)))
    assert e2.is_tangent(Line(p1_2, p2 + Point(half, 1)))
    assert e2.is_tangent(Line(p1_3, p2 + Point(half, 1)))
    assert c1.is_tangent(Line(p1_1, Point(0, sqrt(2))))
    assert e1.is_tangent(Line(Point(0, 0), Point(1, 1))) == False

    # Intersection
    l1 = Line(Point(1, -5), Point(1, 5))
    l2 = Line(Point(-5, -1), Point(5, -1))
    l3 = Line(Point(-1, -1), Point(1, 1))
    l4 = Line(Point(-10, 0), Point(0, 10))
    pts_c1_l3 = [Point(sqrt(2)/2, sqrt(2)/2), Point(-sqrt(2)/2, -sqrt(2)/2)]

    assert intersection(e2, l4) == []
    assert intersection(c1, Point(1, 0)) == [Point(1, 0)]
    assert intersection(c1, l1) == [Point(1, 0)]
    assert intersection(c1, l2) == [Point(0, -1)]
    assert intersection(c1, l3) in [pts_c1_l3, [pts_c1_l3[1], pts_c1_l3[0]]]
    assert intersection(c1, c2) in [[(1,0), (0,1)],[(0,1),(1,0)]]
    assert intersection(c1, c3) == [(sqrt(2)/2, sqrt(2)/2)]

    v = sqrt(2)
    t1 = Triangle(Point(0, v), Point(0, -v), Point(v, 0))
    points = intersection(t1, c1)
    assert len(points) == 4
    assert Point(0, 1) in points
    assert Point(0, -1) in points
    assert Point(v/2, v/2) in points
    assert Point(v/2, -v/2) in points

    e1 = Circle(Point(0, 0), 5)
    e2 = Ellipse(Point(0, 0), 5, 20)
    assert intersection(e1, e2) in \
        [[Point(5, 0), Point(-5, 0)], [Point(-5, 0), Point(5, 0)]]

    # Combinations of above
    assert e3.is_tangent(e3.tangent_line(p1 + Point(y1, 0)))

def test_polygon():
    p1 = Polygon(
        Point(0, 0), Point(3,-1),
        Point(6, 0), Point(4, 5),
        Point(2, 3), Point(0, 3))
    p2 = Polygon(
        Point(6, 0), Point(3,-1),
        Point(0, 0), Point(0, 3),
        Point(2, 3), Point(4, 5))
    p3 = Polygon(
        Point(0, 0), Point(3, 0),
        Point(5, 2), Point(4, 4))
    p4 = Polygon(
        Point(0, 0), Point(4, 4),
        Point(5, 2), Point(3, 0))

    #
    # General polygon
    #
    assert p1 == p2
    assert len(p1) == Rational(6)
    assert len(p1.sides) == 6
    assert p1.perimeter == 5+2*sqrt(10)+sqrt(29)+sqrt(8)
    assert p1.area == 22
    assert not p1.is_convex()
    assert p3.is_convex()
    assert p4.is_convex()  # ensure convex for both CW and CCW point specification

    #
    # Regular polygon
    #
    p1 = RegularPolygon(Point(0, 0), 10, 5)
    p2 = RegularPolygon(Point(0, 0), 5, 5)

    assert p1 != p2
    assert p1.interior_angle == 3*pi/5
    assert p1.exterior_angle == 2*pi/5
    assert p2.apothem == 5*cos(pi/5)
    assert p2.circumcircle == Circle(Point(0, 0), 5)
    assert p2.incircle == Circle(Point(0, 0), p2.apothem)
    assert p1.is_convex()

    #
    # Angles
    #
    angles = p4.angles
    assert feq(angles[Point(0, 0)].evalf(), Real("0.7853981633974483"))
    assert feq(angles[Point(4, 4)].evalf(), Real("1.2490457723982544"))
    assert feq(angles[Point(5, 2)].evalf(), Real("1.8925468811915388"))
    assert feq(angles[Point(3, 0)].evalf(), Real("2.3561944901923449"))

    angles = p3.angles
    assert feq(angles[Point(0, 0)].evalf(), Real("0.7853981633974483"))
    assert feq(angles[Point(4, 4)].evalf(), Real("1.2490457723982544"))
    assert feq(angles[Point(5, 2)].evalf(), Real("1.8925468811915388"))
    assert feq(angles[Point(3, 0)].evalf(), Real("2.3561944901923449"))

    #
    # Triangle
    #
    p1 = Point(0, 0)
    p2 = Point(5, 0)
    p3 = Point(0, 5)
    t1 = Triangle(p1, p2, p3)
    t2 = Triangle(p1, p2, Point(Rational(5,2), sqrt(Rational(75,4))))
    t3 = Triangle(p1, Point(x1, 0), Point(0, x1))
    s1 = t1.sides
    s2 = t2.sides
    s3 = t3.sides

    # Basic stuff
    assert t1.area == Rational(25,2)
    assert t1.is_right()
    assert t2.is_right() == False
    assert t3.is_right()
    assert p1 in t1
    assert Point(5, 5) not in t2
    assert t1.is_convex()
    assert feq(t1.angles[p1].evalf(), pi.evalf()/2)

    assert t1.is_equilateral() == False
    assert t2.is_equilateral()
    assert t3.is_equilateral() == False
    assert are_similar(t1, t2) == False
    assert are_similar(t1, t3)
    assert are_similar(t2, t3) == False

    # Bisectors
    bisectors = t1.bisectors
    assert bisectors[p1] == Segment(p1, Point(Rational(5,2), Rational(5,2)))
    ic = (250 - 125*sqrt(2)) / 50
    assert t1.incenter == Point(ic, ic)

    # Medians + Centroid
    m = t1.medians
    assert t1.centroid == Point(Rational(5,3), Rational(5,3))
    assert m[p1] == Segment(p1, Point(Rational(5,2), Rational(5,2)))
    assert t3.medians[p1] == Segment(p1, Point(x1/2, x1/2))
    assert intersection(m[p1], m[p2], m[p3]) == [t1.centroid]

    # Perpendicular
    altitudes = t1.altitudes
    assert altitudes[p1] == Segment(p1, Point(Rational(5,2), Rational(5,2)))
    assert altitudes[p2] == s1[0]
    assert altitudes[p3] == s1[2]

    # Ensure
    assert len(intersection(*bisectors.values())) == 1
    assert len(intersection(*altitudes.values())) == 1
    assert len(intersection(*m.values())) == 1

def test_convex_hull():
    p = [Point(-5,-1), Point(-2,1), Point(-2,-1), Point(-1,-3), Point(0,0),
         Point(1,1), Point(2,2), Point(2,-1), Point(3,1), Point(4,-1), Point(6,2)]
    ch = Polygon(p[0], p[3], p[9], p[10], p[6], p[1])
    #test handling of duplicate points
    p.append(p[3])

    #more than 3 collinear points
    another_p = [Point(-45, -85), Point(-45, 85), Point(-45,26),Point(-45,-24)]
    ch2 = Segment(another_p[0],another_p[1]);

    assert convex_hull(another_p) == ch2
    assert convex_hull(p) == ch
    assert convex_hull(p[0]) == p[0]
    assert convex_hull(p[0], p[1]) == Segment(p[0], p[1])

def test_concyclic_doctest_bug():
    p1,p2 = Point(-1, 0), Point(1, 0)
    p3,p4 = Point(0, 1), Point(-1, 2)
    assert Point.is_concyclic(p1, p2, p3)
    assert not Point.is_concyclic(p1, p2, p3, p4)

if __name__ == "__main__":
    from sys import modules,stderr,exc_info,excepthook
    import hotshot, hotshot.stats
    curr_module = modules[__name__]
    for member in dir(curr_module):
        if member.startswith("test_"):
            try:
                sys.stderr.write("Testing %s..." % member)
                curr_module.__dict__[member]()
                sys.stderr.write("SUCCESS!\n")
            except AssertionError, e:
                sys.stderr.write("FAILED!\n")
                sys.stderr.write('-' * 25 + '\n')
                excepthook(*exc_info())
                sys.stderr.write('-' * 25 + '\n')
