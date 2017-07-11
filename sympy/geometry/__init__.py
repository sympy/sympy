"""
A geometry module for the SymPy library. This module contains all of the
entities and functions needed to construct basic geometrical data and to
perform simple informational queries.

Usage:
======

Examples
========

"""
from sympy.geometry.curve import Curve
from sympy.geometry.ellipse import Circle, Ellipse
from sympy.geometry.exceptions import GeometryError
from sympy.geometry.line import Line, Line2D, Line3D, Ray, Ray2D, Ray3D, \
    Segment, Segment2D, Segment3D
from sympy.geometry.parabola import Parabola
from sympy.geometry.plane import Plane
from sympy.geometry.point import Point, Point2D, Point3D
from sympy.geometry.polygon import Polygon, RegularPolygon, Triangle, deg, rad
from sympy.geometry.util import are_similar, centroid, closest_points, \
    convex_hull, farthest_points, idiff, intersection
