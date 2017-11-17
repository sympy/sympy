"""
A geometry module for the SymPy library. This module contains all of the
entities and functions needed to construct basic geometrical data and to
perform simple informational queries.

Usage:
======

Examples
========

"""

__all__ = []

from sympy.geometry.point import Point, Point2D, Point3D
__all__ += ["Point", "Point2D", "Point3D"]

from sympy.geometry.line import (
    Line, Ray, Segment, Line2D, Segment2D, Ray2D,
    Line3D, Segment3D, Ray3D
)
__all__ += [
    "Line", "Ray", "Segment", "Line2D", "Segment2D", "Ray2D",
    "Line3D", "Segment3D", "Ray3D"
]

from sympy.geometry.plane import Plane
__all__ += ["Plane"]

from sympy.geometry.ellipse import Ellipse, Circle
__all__ += ["Ellipse", "Circle"]

from sympy.geometry.polygon import Polygon, RegularPolygon, Triangle, rad, deg
__all__ += ["Polygon", "RegularPolygon", "Triangle", "rad", "deg"]

from sympy.geometry.util import (
    are_similar, centroid, convex_hull, idiff,
    intersection, closest_points, farthest_points
)
__all__ += [
    "are_similar", "centroid", "convex_hull", "idiff",
    "intersection", "closest_points", "farthest_points"
]

from sympy.geometry.exceptions import GeometryError
__all__ += ["GeometryError"]

from sympy.geometry.curve import Curve
__all__ += ["Curve"]

from sympy.geometry.parabola import Parabola
__all__ += ["Parabola"]
