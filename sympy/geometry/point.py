from sympy.core.basic import Basic, S, sympify
from sympy.simplify import simplify
from sympy.geometry.exceptions import GeometryError
from sympy.functions.elementary.miscellaneous import sqrt
from entity import GeometryEntity


class Point(GeometryEntity):
    """
    A point in Euclidean N-space defined by a sequence of values. Can be
    constructed from a sequence of points or a list of points.

    Examples:
    ======
        >>> from sympy.geometry import Point
        >>> Point(1, 2)
        Point(1, 2)
        >>> Point([1, 2])
        Point(1, 2)

    Notes:
    ======
        - Currently only 2-dimensional points are supported.
    """
    def __new__(cls, *args, **kwargs):
        if isinstance(args[0], (tuple, list, set)):
            coords = tuple([sympify(x) for x in args[0]])
        else:
            coords = tuple([sympify(x) for x in args])

        if len(coords) != 2:
            raise NotImplementedError("Only two dimensional points currently supported")

        return GeometryEntity.__new__(cls, *coords)

    def is_collinear(*points):
        """
        Test whether or not a set of points are collinear. Returns True if
        the set of points are collinear, or False otherwise.

        Examples:
        =========
            >>> from sympy import *
            >>> from sympy.abc import x
            >>> p1,p2 = Point(0, 0), Point(1, 1)
            >>> p3,p4,p5 = Point(2, 2), Point(x, x), Point(1, 2)
            >>> Point.is_collinear(p1, p2, p3, p4)
            True
            >>> Point.is_collinear(p1, p2, p3, p5)
            False

        Description of method used:
        ===========================
            Slope is preserved everywhere on a line, so the slope between
            any two points on the line should be the same. Take the first
            two points, p1 and p2, and create a translated point v1
            with p1 as the origin. Now for every other point we create
            a translated point, vi with p1 also as the origin. Note that
            these translations preserve slope since everything is
            consistently translated to a new origin of p1. Since slope
            is preserved then we have the following equality:
                    v1_slope = vi_slope
              =>    v1.y/v1.x = vi.y/vi.x (due to translation)
              =>    v1.y*vi.x = vi.y*v1.x
              =>    v1.y*vi.x - vi.y*v1.x = 0           (*)
            Hence, if we have a vi such that the equality in (*) is False
            then the points are not collinear. We do this test for every
            point in the list, and if all pass then they are collinear.
        """
        points = GeometryEntity.extract_entities(points)
        if len(points) == 0: return False
        if len(points) <= 2: return True # two points always form a line

        # XXX Cross product is used now, but that only extends to three
        #     dimensions. If the concept needs to extend to greater
        #     dimensions then another method would have to be used
        p1 = points[0]
        p2 = points[1]
        v1 = p2 - p1
        for p3 in points[2:]:
            v2 = p3 - p1
            test = simplify(v1[0]*v2[1] - v1[1]*v2[0])
            if simplify(test) != 0:
                return False
        return True

    def is_concyclic(*points):
        """
        Test whether or not a set of points are concyclic (i.e., on the same
        circle). Returns True if they are concyclic, or False otherwise.

        Example:
        ========
            >>> from sympy.geometry import Point
            >>> p1,p2 = Point(-1, 0), Point(1, 0)
            >>> p3,p4 = Point(0, 1), Point(-1, 2)
            >>> Point.is_concyclic(p1, p2, p3)
            True
            >>> Point.is_concyclic(p1, p2, p3, p4)
            False

        Description of method used:
        ===========================
            No points are not considered to be concyclic. One or two points
            are definitely concyclic and three points are conyclic iff they
            are not collinear.

            For more than three points, we pick the first three points and
            attempt to create a circle. If the circle cannot be created
            (i.e., they are collinear) then all of the points cannot be
            concyclic. If the circle is created successfully then simply
            check all of the other points for containment in the circle.
        """
        points = GeometryEntity.extract_entities(points)
        if len(points) == 0: return False
        if len(points) <= 2: return True
        if len(points) == 3: return (not Point.is_collinear(*points))

        try:
            from ellipse import Circle
            c = Circle(points[0], points[1], points[2])
            for point in points[3:]:
                if point not in c:
                    return False
            return True
        except GeometryError,e:
            # Circle could not be created, because of collinearity of the
            # three points passed in, hence they are not concyclic.
            return False
        """
        # This code is from Maple
        def f(u):
            dd = u[0]**2 + u[1]**2 + 1
            u1 = 2*u[0] / dd
            u2 = 2*u[1] / dd
            u3 = (dd - 2) / dd
            return u1,u2,u3

        u1,u2,u3 = f(points[0])
        v1,v2,v3 = f(points[1])
        w1,w2,w3 = f(points[2])
        p = [v1 - u1, v2 - u2, v3 - u3]
        q = [w1 - u1, w2 - u2, w3 - u3]
        r = [p[1]*q[2] - p[2]*q[1], p[2]*q[0] - p[0]*q[2], p[0]*q[1] - p[1]*q[0]]
        for ind in xrange(3, len(points)):
            s1,s2,s3 = f(points[ind])
            test = simplify(r[0]*(s1-u1) + r[1]*(s2-u2) + r[2]*(s3-u3))
            if test != 0:
                return False
        return True
        """

    @staticmethod
    def distance(p1, p2):
        """
        Get the Euclidean distance between two points.

        Example:
        ========
            >>> from sympy.geometry import Point
            >>> p1,p2 = Point(1, 1), Point(4, 5)
            >>> Point.distance(p1, p2)
            5

        """
        return sqrt( sum([(a-b)**2 for a,b in zip(p1,p2)]) )

    @staticmethod
    def midpoint(p1, p2):
        """
        Get the midpoint of two points.

        Example:
        ========
            >>> from sympy.geometry import Point
            >>> p1,p2 = Point(1, 1), Point(13, 5)
            >>> Point.midpoint(p1, p2)
            Point(7, 3)

        """
        return Point( [simplify((a + b)*S.Half) for a,b in zip(p1,p2)] )

    def evalf(self):
        """
        Evaluate and return a Point where every coordinate is evaluated to
        a floating point number.

        Example:
        ========
            >>> from sympy import *
            >>> Point(Rational(1,2), Rational(3,2)).evalf()
            Point(0.5, 1.5)
        """
        return Point([x.evalf() for x in self])

    def intersection(self, o):
        if isinstance(o, Point):
            if self == o:
                return [self]
            return []
        raise NotImplementedError()

    def __add__(self, other):
        """
        Create a new point where each coordinate in this point is
        increased by the corresponding coordinate in other.
        """
        if isinstance(other, Point):
            if len(other) == len(self):
                return Point( [simplify(a+b) for a,b in zip(self, other)] )
            else:
                raise TypeError("Points must have the same number of dimensions")
        else:
            other = sympify(other)
            return Point( [simplify(a+other) for a in self] )

    def __sub__(self, other):
        """
        Create a new point where each coordinate in this point is
        decreased by the corresponding coordinate in other.
        """
        return self + (-other)

    def __mul__(self, factor):
        """
        Create a new point where each coordinate in this point is
        multiplied by factor.
        """
        factor = sympify(factor)
        return Point( [x*factor for x in self] )

    def __div__(self, divisor):
        """
        Create a new point where each coordinate in this point is
        divided by factor.
        """
        divisor = sympify(divisor)
        return Point( [x/divisor for x in self] )

    def __neg__(self):
        """
        Create a new point where each oordinate in this point is negated.
        """
        return Point( [-x for x in self] )

    def __abs__(self):
        """Returns the distance between this point and the origin."""
        origin = Point([0] * len(self))
        return Point.distance(origin, self)
