"""
A geometry module for the SymPy library. This module contains all of the
entities and functions needed to construct basic geometrical data and to
perform simple informational queries.

Usage:
======


Notes:
======
    Currently the geometry module is restricted to the 2-dimensional
    Euclidean space.

Examples:
=========

"""
from sympy.geometry.point import Point
from sympy.geometry.line import Line, Ray, Segment
from sympy.geometry.ellipse import Ellipse, Circle
from sympy.geometry.polygon import Polygon, RegularPolygon, Triangle
from sympy.geometry.util import *
from sympy.geometry.exceptions import *
from sympy.geometry.curve import Curve
