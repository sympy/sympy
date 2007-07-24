from sympy import *
from sympy.modules.simplify import simplify
from entity import GeometryEntity

class Point(GeometryEntity):
    """A point in space defined by a sequence of values."""

    def __init__(self, *args, **kwargs):
        GeometryEntity.__init__(self, *args, **kwargs)
        if isinstance(args[0], (Basic, int, float)):
            self._coords = tuple([Basic.sympify(x) for x in args])
        else:
            self._coords = tuple([Basic.sympify(x) for x in args[0]])
        if len(self._coords) != 2:
            raise NotImplementedError("Only two dimensional points currently supported")

    def is_collinear(*args):
        """
        Test whether or not a set of points are collinear. Returns True if
        the set of points are collinear, or False otherwise.
        """
        points = GeometryEntity._normalize_args(args)

        for x in points:
            assert isinstance(x, Point)
        if len(points) == 0: return False
        if len(points) <= 2: return True

        # TODO Cross ratio is used now, but that only extends to three dimensions.
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
            if bool(test != 0):
                return False
        return True

    def is_concyclic(*args):
        """
        Test whether or not a set of points are concyclic (i.e., on the same
        circle). Returns True if they are concyclic, or False otherwise.
        """
        from sympy import Add,Mul

        points = GeometryEntity._normalize_args(args)
        if len(points) == 0: return False
        if len(points) <= 2: return True
        if len(points) == 3: return (not Point.is_collinear(points))

        def f(u):
            dd = u[0]**Rational(2) + u[1]**Rational(2) + Rational(1)
            u1 = Rational(2)*u[0] / dd
            u2 = Rational(2)*u[1] / dd
            u3 = (dd - Rational(2)) / dd
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
            if bool(test != 0):
                return False
        return True

    @staticmethod
    def distance(p1, p2):
        """Get the Euclidean distance between two points."""
        res = Rational(0)
        for ind in xrange(0, len(p1)):
            res += (p1[ind] - p2[ind])**Rational(2)
        return sqrt(simplify(res))

    @staticmethod
    def midpoint(p1, p2):
        """Get the midpoint of two points."""
        return Point( simplify((p1[ind] + p2[ind])*Rational(1,2)) for ind in xrange(0, len(p1)) )

    def evalf(self, precision=18):
        """
        Evaluate and return a Point where every coordinate is evaluated to
        a floating point number.
        """
        from sympy.core.numbers import Real
        coords = []
        for x in self._coords:
            coords.append( Real(x, precision) )
        return Point(coords)

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
        coord_str = str.join(', ', [str(x) for x in self._coords])
        return "Point(" + coord_str + ")"
