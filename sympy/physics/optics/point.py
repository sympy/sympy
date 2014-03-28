from __future__ import print_function, division

from sympy import Matrix, sqrt, Rational
from sympy.core import S, sympify
from sympy.simplify import simplify, nsimplify
from sympy.core.compatibility import xrange

__all__ = ['Point']


class Point(object):
    """
    This object represents a Point in a Cartesian coordinate system.
    It handles both 2D and 3D computations.
    """

    def __init__(self, *args):
        if(len(args) < 2):
            raise Exception("At least two coordinates are required")
        elif(len(args) > 3):
            raise Exception("Maximum three coordinates are allowed")
        self._coords = args

    def __eq__(self, other):
        ts, to = type(self), type(other)
        if ts is not to:
            return False
        return self._coords == other._coords

    def __sub__(self, other):
        """
        Subtract coordinates of one point from other

        """

        if isinstance(other, Point):
            if(len(self._coords) != len(other._coords)):
                raise TypeError("Dimensions of both points are not same")
            return Point(*[i - j for i, j in zip(self._coords, other._coords)])
        else:
            raise TypeError("Operation not supported for foreign objects")

    def __add__(self, other):
        """
        Add coordinates of one point to other

        """

        if isinstance(other, Point):
            if(len(self._coords) != len(other._coords)):
                raise TypeError("Dimensions of both points are not same")
            return Point(*[i + j for i, j in zip(self._coords, other._coords)])
        else:
            raise TypeError("Operation not supported for foreign objects")

    @property
    def x(self):
        """
        Returns the X coordinate of the Point.

        Examples
        ========

        >>> from sympy.physics.optics import Point
        >>> p = Point(2, 1, 3)
        >>> p.x
        2
        """

        return self._coords[0]

    @property
    def y(self):
        """
        Returns the Y coordinate of the Point.

        Examples
        ========

        >>> from sympy.physics.optics import Point
        >>> p = Point(2, 1, 3)
        >>> p.y
        1
        """

        return self._coords[1]

    @property
    def z(self):
        """
        Returns the Z coordinate of the Point.

        Examples
        ========

        >>> from sympy.physics.optics import Point
        >>> p = Point(2, 1, 3)
        >>> p.z
        3
        """

        if(len(self._coords) > 2):
            return self._coords[2]
        else:
            raise AttributeError("It is a 2-D point")

    @property
    def length(self):
        """
        Treating a Point as a Line, this returns 0 for the length of a Point.

        Examples
        ========

        >>> from sympy.physics.optics import Point
        >>> p = Point(3, 1, 7)
        >>> p.length
        0
        """
        return S.Zero

    def is_collinear(*points):
        """
        This method tests whether passed points are collinear or not.

        It uses slope method for 2-D points and matrix method for 3-D
        points.

        Parameters
        ==========

        points : sequence of Point

        Returns
        =======

        is_collinear : boolean

        Examples
        ========

        >>> from sympy.physics.optics import Point
        >>> p1 = Point(2, 2, 2)
        >>> p2 = Point(-3, -3, -3)
        >>> p3 = Point(0, 0, 0)
        >>> p4 = Point(1, 1, 1)
        >>> p5 = Point(1, 2, 3)
        >>> Point.is_collinear(p1, p2, p3, p4)
        True
        >>> Point.is_collinear(p1, p2, p3, p4, p5)
        False

        """

        points = list(set(points))
        if len(points) == 0:
            return False
        if len(points) <= 2:
            return True  # two points always form a line

        if(len(points[0]._coords) == 2):
            p1 = points[0]
            p2 = points[1]
            v1 = p2 - p1
            x1, y1 = v1._coords
            rv = True
            for p3 in points[2:]:
                x2, y2 = (p3 - p1)._coords
                test = simplify(x1*y2 - y1*x2).equals(0)
                if test is False:
                    return False
                if rv and not test:
                    rv = test
            return rv
        else:
            pv1 = [i - j for i, j in zip(points[0]._coords, points[1]._coords)]
            pv2 = [i - j for i, j in zip(points[1]._coords, points[2]._coords)]
            rank = Matrix([pv1, pv2]).rank()
            if(rank != 1):
                return False
            for i in xrange(1, len(points) - 3):
                pv1 = [i - j for i, j in zip(points[i]._coords,
                                             points[i + 1]._coords)]
                pv2 = [i - j for i, j in zip(points[i + 1]._coords,
                                             points[i + 2]._coords)]
                rank = Matrix([pv1, pv2]).rank()
                if(rank != 1):
                    return False
            return True

    def distance(self, point):
        """
        To calculate Euclidian distance between self and the passed point.

        Parameters
        ==========

        point : a Point

        Returns
        =======

        A number or an expression

        Examples
        ========

        >>> from sympy.physics.optics import Point
        >>> from sympy.abc import x, y
        >>> p1 = Point(2, 2, 2)
        >>> p2 = Point(2, 2, 3)
        >>> p3 = Point(x, 3, y)
        >>> p1.distance(p2)
        1
        >>> p2.distance(p3)
        sqrt((x - 2)**2 + (y - 3)**2 + 1)

        """

        if(len(self._coords) > 2):
            return sqrt(
                (point._coords[0] - self._coords[0])**2 +
                (point._coords[1] - self._coords[1])**2 +
                (point._coords[2] - self._coords[2])**2
            )
        else:
            return sqrt(
                (point._coords[0] - self._coords[0])**2 +
                (point._coords[1] - self._coords[1])**2
            )

    def midpoint(self, point):
        """
        To calculate midpoint between self and the passed point.

        Parameters
        ==========

        point : a Point

        Returns
        =======

        A Point

        Examples
        ========

        >>> from sympy.physics.optics import Point
        >>> from sympy.abc import x, y
        >>> p1 = Point(x, y, 4)
        >>> p1.midpoint(Point(2, 3, 4))
        Point(x/2 + 1, y/2 + 3/2, 4)

        """

        if(len(self._coords) > 2):
            return Point(
                simplify(nsimplify((
                    (self._coords[0] + point._coords[0])/2),
                    rational=True)
                ),
                simplify(nsimplify((
                    (self._coords[1] + point._coords[1])/2),
                    rational=True)
                ),
                simplify(nsimplify((
                    (self._coords[2] + point._coords[2])/2),
                    rational=True)
                )
            )
        else:
            return Point(
                simplify(nsimplify((
                    (self._coords[0] + point._coords[0])/2),
                    rational=True)
                ),
                simplify(nsimplify((
                    (self._coords[1] + point._coords[1])/2),
                    rational=True)
                )
            )

    def __repr__(self):
        if(len(self._coords) > 2):
            return "Point(" + repr(self._coords[0]) + ", " + \
                repr(self._coords[1]) + ", " + repr(self._coords[2]) + ")"
        else:
            return "Point(" + repr(self._coords[0]) + ", " + \
                repr(self._coords[1]) + ")"
