from __future__ import print_function, division

from sympy.core import S, sympify

__all__ = ['Point']


class Point(object):
    """
    This object represents a Point in a Cartesian coordinate
    system.
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

    def is_collinear(self, *points):
        points = list(set(points))
        if(len(self._coords) == 2):
            if len(points) == 0:
                return False
            if len(points) <= 2:
                return True  # two points always form a line
            points = [Point(a) for a in points]
            p1 = points[0]
            p2 = points[1]
            v1 = p2 - p1
            x1, y1 = v1.args
            rv = True
            for p3 in points[2:]:
                x2, y2 = (p3 - p1).args
                test = simplify(x1*y2 - y1*x2).equals(0)
                if test is False:
                    return False
                if rv and not test:
                    rv = test
            return rv
        #TODO: 3D collinearity

    def __repr__(self):
        if(len(self._coords) > 2):
            return "Point(" + repr(self._coords[0]) + ", " + repr(self._coords[1])\
                + ", " + repr(self._coords[2]) + ")"
        else:
            return "Point(" + repr(self._coords[0]) + ", " +\
                repr(self._coords[1]) + ")"
