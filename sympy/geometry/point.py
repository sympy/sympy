from sympy.core.basic import Basic, S
from sympy.simplify import simplify
from entity import GeometryEntity

class Point(GeometryEntity):
    """
    A point in Euclidean N-space defined by a sequence of values. Can be
    constructed from a sequence of points or a list of points.

    Examples:
    ======
        >>> Point(1, 2)
        Point(1, 2)
        >>> Point([1, 2])
        Point(1, 2)

    Notes:
    ======
        Currently only 2-dimensional points are supported.
    """

    def __new__(self, *args, **kwargs):
        if isinstance(args[0], (Basic, int, float)):
            coords = tuple([Basic.sympify(x) for x in args])
        else:
            coords = tuple([Basic.sympify(x) for x in args[0]])

        if len(coords) != 2:
            raise NotImplementedError("Only two dimensional points currently supported")

        obj = GeometryEntity.__new__(self, *coords, **kwargs)
        obj._coords = coords
        return obj

    def is_collinear(*args):
        """
        Test whether or not a set of points are collinear. Returns True if
        the set of points are collinear, or False otherwise.

        Examples:
        =========
            >>> from sympy import *
            >>> x = Symbol('x')
            >>> p1,p2 = Point(0, 0), Point(1, 1)
            >>> p3,p4,p5 = Point(2, 2), Point(x, x), Point(1, 2)
            >>> Point.is_collinear(p1, p2, p3, p4)
            True
            >>> Point.is_collinear(p1, p2, p3, p5)
            False
        """
        points = GeometryEntity._normalize_args(args)
        if len(points) == 0: return False
        if len(points) <= 2: return True

        # NOTE Cross ratio is used now, but that only extends to three dimensions.
        #      If the concept needs to extend to greater dimensions then one of
        #      the other methods (as listed below) would have to be used

        # XXX Should another method be used. Cross product is used but there is
        #     also: cross ratio of distances, area of triangle, membership on line
        p1 = points[0]
        p2 = points[1]
        v1 = p1 - p2
        for ind in xrange(2, len(points)):
            p3 = points[ind]
            v2 = p1 - p3
            test = simplify(v1[0]*v2[1] - v1[1]*v2[0])
            if test != 0:
                return False
        return True

    def is_concyclic(*args):
        """
        Test whether or not a set of points are concyclic (i.e., on the same
        circle). Returns True if they are concyclic, or False otherwise.

        Example:
        ========
            >>> p1,p2 = Point(-1, 0), Point(1, 0)
            >>> p3,p4 = Point(0, 1), Point(-1, 2)
            >>> Point.is_concyclic(p1, p2, p3)
            True
            >>> Point.is_concyclic(p1, p2, p3, p4)
            False
        """
        points = GeometryEntity._normalize_args(args)
        if len(points) == 0: return False
        if len(points) <= 2: return True
        if len(points) == 3: return (not Point.is_collinear(*points))

        def f(u):
            dd = u[0]**2 + u[1]**2 + 1
            u1 = 2*u[0] / dd
            u2 = 2*u[1] / dd
            u3 = (dd - 2) / dd
            return dd,u1,u2,u3

        d1,u1,u2,u3 = f(points[0])
        d2,v1,v2,v3 = f(points[1])
        d3,w1,w2,w3 = f(points[2])
        for ind in xrange(3, len(points)):
            d4,s1,s2,s3 = f(points[ind])
            p = [v1 - u1, v2 - u2, v3 - u3]
            q = [w1 - u1, w2 - u2, w3 - u3]
            r = [p[1]*q[2] - p[2]*q[1], p[2]*q[0] - p[0]*q[2], p[0]*q[1] - p[1]*q[0]]

            test = simplify(r[0]*(s1-u1) + r[1]*(s2-u2) + r[2]*(s3-u3))
            if test != 0:
                return False
        return True

    @staticmethod
    def distance(p1, p2):
        """
        Get the Euclidean distance between two points.

        Example:
        ========
            >>> p1,p2 = Point(1, 1), Point(4, 5)
            >>> Point.distance(p1, p2)
            5
        """
        res = 0
        for ind in xrange(0, len(p1)):
            res += (p1[ind] - p2[ind])**2
        return S.Sqrt(simplify(res))

    @staticmethod
    def midpoint(p1, p2):
        """
        Get the midpoint of two points.

        Example:
        ========
            >>> p1,p2 = Point(1, 1), Point(13, 5)
            >>> Point.midpoint(p1, p2)
            Point(7, 3)
        """
        return Point( simplify((p1[ind] + p2[ind])*S.Half) for ind in xrange(0, len(p1)) )

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
        coords = []
        for x in self._coords:
            coords.append( x.evalf() )
        return Point(coords)

    def _intersection(self, o):
        if isinstance(o, Point):
            if self == o:
                return [self]
            return []
        raise NotImplementedError()

    def __getitem__(self, ind):
        """Get a specific coordinate."""
        return self._coords[ind]

    def __len__(self):
        """Return the number of coordinates for this point."""
        return len(self._coords)

    def __eq__(self, p):
        try:
            if len(p) != len(self): return False
            for ind in xrange(0, len(self)):
                if self[ind] != p[ind]: return False
            return True
        except:
            return False

    def __add__(self, a):
        """Coordinate-based addition."""
        if isinstance(a, Point):
            if len(a) == len(self):
                return Point([simplify(self[ind] + a[ind]) for ind in xrange(0, len(a))])
            else:
                raise Exception("Points must have the same number of dimensions")
        else:
            a = Basic.sympify(a)
            return Point([simplify(self[ind] + a) for ind in xrange(0, len(a))])

    def __sub__(self, a):
        """Coordinate-based subtraction."""
        if isinstance(a, Point):
            if len(a) == len(self):
                return Point([simplify(self[ind] - a[ind]) for ind in xrange(0, len(a))])
            else:
                raise Exception("Points must have the same number of dimensions")
        else:
            a = Basic.sympify(a)
            return Point([simplify(self[ind] - a) for ind in xrange(0, len(a))])

    def __mul__(self, a):
        """Coordinate-based multiplication."""
        a = Basic.sympify(a)
        return Point( [x*a for x in self] )

    def __div__(self, a):
        """Coordinate-based division."""
        a = Basic.sympify(a)
        return Point( [x/a for x in self] )

    def __neg__(self):
        """Coordinate-based negation."""
        return Point([-x for x in self])

    def __abs__(self):
        """Returns the distance between this point and the origin."""
        origin = Point([0 for x in xrange(0, len(self))])
        return Point.distance(origin, self)

    def __str__(self):
        c_str = str.join(', ', [str(x) for x in self._coords])
        return "Point(" + c_str + ")"

    def __repr__(self):
        c_str = str.join(', ', [str(x) for x in self._coords])
        return "Point(" + c_str + ")"
