"""
A geometry module for the SymPy library. This module contains all of the
entities and functions needed to construct basic geometrical data and to
perform simple informational queries.

Usage:
======


Notes:
======
    Currently the geometry module supports 2-dimensional
    and 3 -dimensional Euclidean space.

Examples
========

"""

from __future__ import division, print_function

from sympy.geometry.point import Point, Point2D, Point3D
from sympy.geometry.line import Line, Ray, Segment
from sympy.geometry.line3d import Line3D, Segment3D, Ray3D
from sympy.geometry.plane import Plane
from sympy.geometry.ellipse import Ellipse, Circle
from sympy.geometry.polygon import Polygon, RegularPolygon, Triangle, rad, deg
from sympy.geometry.util import are_similar, centroid, convex_hull, idiff, \
    intersection
from sympy.geometry.exceptions import GeometryError
from sympy.geometry.curve import Curve
