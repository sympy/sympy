from __future__ import division
import warnings

from sympy import (Abs, I, Dummy, Rational, Float, S, Symbol, cos, oo, pi,
                   simplify, sin, sqrt, symbols, Derivative, asin, acos)
from sympy.core.compatibility import range
from sympy.functions.elementary.trigonometric import tan
from sympy.geometry import (Circle, Curve, Ellipse, GeometryError, Line, Point,
                            Polygon, Ray, RegularPolygon, Segment, Triangle,
                            are_similar, convex_hull, intersection,
                            Point3D, Line3D, Ray3D, Segment3D, Plane, centroid)
from sympy.geometry.line import Undecidable
from sympy.geometry.entity import rotate, scale, translate
from sympy.geometry.polygon import _asa as asa, rad, deg
from sympy.geometry.util import idiff, are_coplanar
from sympy.sets import (Set, FiniteSet, Union, Intersection, Complement, EmptySet)
from sympy.integrals.integrals import Integral
from sympy.matrices import Matrix
from sympy.solvers.solvers import solve
from sympy.utilities.iterables import cartes
from sympy.utilities.randtest import verify_numerically
from sympy.utilities.pytest import raises, slow

x = Symbol('x', real=True)
y = Symbol('y', real=True)
z = Symbol('z', real=True)
t = Symbol('t', real=True)
k = Symbol('k', real=True)
x1 = Symbol('x1', real=True)
x2 = Symbol('x2', real=True)
x3 = Symbol('x3', real=True)
y1 = Symbol('y1', real=True)
y2 = Symbol('y2', real=True)
y3 = Symbol('y3', real=True)
z1 = Symbol('z1', real=True)
z2 = Symbol('z2', real=True)
z3 = Symbol('z3', real=True)
half = Rational(1, 2)

p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
p5, p6, p7 = map(Point, [(3, 2), (1, -1), (0, 2)])
l1 = Line(Point(0,0), Point(1,1))
l2 = Line(Point(half, half), Point(5,5))
l3 = Line(p2, p3)
l4 = Line(p3, p4)
poly1 = Polygon(p1, p2, p3, p4)
poly2 = Polygon(p5, p6, p7)
poly3 = Polygon(p1, p2, p5)

def test_booleans():
    """ test basic unions and intersections """
    assert Union(l1, l2).equals(l1)
    assert Intersection(l1, l2).equals(l1)
    assert Intersection(l1, l4) == FiniteSet(Point(1,1))
    assert Intersection(Union(l1, l4), l3) == FiniteSet(Point(-1/3, -1/3), Point(5, 1))
    assert Intersection(l1, FiniteSet(Point(7,-7))) == EmptySet()
    assert Intersection(Circle(Point(0,0), 3), Line(p1,p2)) == FiniteSet(Point(-3,0), Point(3,0))

    fs = FiniteSet(Point(1/3, 1), Point(2/3, 0), Point(9/5, 1/5), Point(7/3, 1))
    # test the intersection of polygons
    assert Intersection(poly1, poly2) == fs
    # make sure if we union polygons with subsets, the subsets go away
    assert Union(poly1, poly2, fs) == Union(poly1, poly2)
    # make sure that if we union with a FiniteSet that isn't a subset,
    # that the points in the intersection stop being listed
    assert Union(poly1, FiniteSet(Point(0,0), Point(3,5))) == Union(poly1, FiniteSet(Point(3,5)))
    # intersect two polygons that share an edge
    assert Intersection(poly1, poly3) == Union(FiniteSet(Point(3/2, 1), Point(2, 1)), Segment(Point(0, 0), Point(1, 0)))
