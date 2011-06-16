"""Geometrical Points.

Contains
--------
Point

"""

from sympy.core import S, sympify
from sympy.core.compatibility import iterable
from sympy.simplify import simplify
from sympy.geometry.exceptions import GeometryError
from sympy.functions.elementary.miscellaneous import sqrt
from entity import GeometryEntity


class Point(GeometryEntity):
    """A point in a 2-dimensional Euclidean space.

    Parameters
    ----------
    coords : sequence of 2 coordinate values.

    Attributes
    ----------
    coordinates : 2-tuple of numbers or sympy objects.
        Stored in `self`. That is self[0] is the first coordinate value, and
        self[1] is the second coordinate value.

    Raises
    ------
    NotImplementedError
        When trying to create a point with more than two dimensions.
        When `intersection` is called with object other than a Point.
    TypeError
        When trying to add or subtract points with different dimensions.

    Notes
    -----
    Currently only 2-dimensional points are supported.

    Examples
    --------
    >>> from sympy.geometry import Point
    >>> from sympy.abc import x
    >>> Point(1, 2)
    Point(1, 2)
    >>> Point([1, 2])
    Point(1, 2)
    >>> Point(0, x)
    Point(0, x)

    """
    def __new__(cls, *args, **kwargs):
        if iterable(args[0]):
            coords = tuple([sympify(x) for x in args[0]])
        else:
            coords = tuple([sympify(x) for x in args])

        if len(coords) != 2:
            raise NotImplementedError("Only two dimensional points currently supported")

        return GeometryEntity.__new__(cls, *coords)

    @property
    def x(self):
        return self[0]

    @property
    def y(self):
        return self[1]

    @property
    def free_symbols(self):
        return self.x.free_symbols.union(self.y.free_symbols)

    def _eval_subs(self, old, new):
        return type(self)(self.x.subs(old, new), self.y.subs(old, new))

    def is_collinear(*points):
        """Is a sequence of points collinear?

        Test whether or not a set of points are collinear. Returns True if
        the set of points are collinear, or False otherwise.

        Parameters
        ----------
        points : sequence of Point

        Returns
        -------
        is_collinear : boolean

        Notes
        --------------------------
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

        Examples
        --------
        >>> from sympy import Point
        >>> from sympy.abc import x
        >>> p1, p2 = Point(0, 0), Point(1, 1)
        >>> p3, p4, p5 = Point(2, 2), Point(x, x), Point(1, 2)
        >>> Point.is_collinear(p1, p2, p3, p4)
        True
        >>> Point.is_collinear(p1, p2, p3, p5)
        False

        """
        if len(points) == 0:
            return False
        if len(points) <= 2:
            return True # two points always form a line
        points = [Point(a) for a in points]

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
        """Is a sequence of points concyclic?

        Test whether or not a sequence of points are concyclic (i.e., they lie
        on a circle).

        Parameters
        ----------
        points : sequence of Points

        Returns
        -------
        is_concyclic : boolean
            True if points are concyclic, False otherwise.

        Notes
        -----
        No points are not considered to be concyclic. One or two points
        are definitely concyclic and three points are conyclic iff they
        are not collinear.

        For more than three points, create a circle from the first three
        points. If the circle cannot be created (i.e., they are collinear)
        then all of the points cannot be concyclic. If the circle is created
        successfully then simply check the remaining points for containment
        in the circle.

        Examples
        --------
        >>> from sympy.geometry import Point
        >>> p1, p2 = Point(-1, 0), Point(1, 0)
        >>> p3, p4 = Point(0, 1), Point(-1, 2)
        >>> Point.is_concyclic(p1, p2, p3)
        True
        >>> Point.is_concyclic(p1, p2, p3, p4)
        False

        """
        if len(points) == 0:
            return False
        if len(points) <= 2:
            return True
        points = [Point(p) for p in points]
        if len(points) == 3:
            return (not Point.is_collinear(*points))

        try:
            from ellipse import Circle
            c = Circle(points[0], points[1], points[2])
            for point in points[3:]:
                if point not in c:
                    return False
            return True
        except GeometryError, e:
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

    def distance(self, p):
        """The Euclidean distance from self to point p.

        Parameters
        ----------
        p : Point

        Returns
        -------
        distance : number or symbolic expression.

        Examples
        --------
        >>> from sympy.geometry import Point
        >>> p1, p2 = Point(1, 1), Point(4, 5)
        >>> p1.distance(p2)
        5

        >>> from sympy.abc import x, y
        >>> p3 = Point(x, y)
        >>> p3.distance(Point(0, 0))
        (x**2 + y**2)**(1/2)

        """
        return sqrt(sum([(a - b)**2 for a, b in zip(self, p)]))

    def midpoint(self, p):
        """The midpoint between self and point p.

        Parameters
        ----------
        p : Point

        Returns
        -------
        midpoint : Point

        Examples
        --------
        >>> from sympy.geometry import Point
        >>> p1, p2 = Point(1, 1), Point(13, 5)
        >>> p1.midpoint(p2)
        Point(7, 3)

        """
        return Point([simplify((a + b)*S.Half) for a, b in zip(self, p)])

    def evalf(self):
        """Evaluate the coordinates of the point.

        This method will, where possible, create and return a new Point
        where the coordinates are evaluated as floating point numbers.

        Returns
        -------
        point : Point

        Examples
        --------
        >>> from sympy import Point, Rational
        >>> p1 = Point(Rational(1, 2), Rational(3, 2))
        >>> p1
        Point(1/2, 3/2)
        >>> p1.evalf()
        Point(0.5, 1.5)

        """
        return Point([x.evalf() for x in self])

    def intersection(self, o):
        """The intersection between this point and another point.

        Parameters
        ----------
        other : Point

        Returns
        -------
        intersection : list of Points

        Notes
        -----
        The return value will either be an empty list if there is no
        intersection, otherwise it will contain this point.

        Examples
        --------
        >>> from sympy import Point
        >>> p1, p2, p3 = Point(0, 0), Point(1, 1), Point(0, 0)
        >>> p1.intersection(p2)
        []
        >>> p1.intersection(p3)
        [Point(0, 0)]

        """
        if isinstance(o, Point):
            if self == o:
                return [self]
            return []

        return o.intersection(self)

    @property
    def length(self):
        return S.Zero

    def __len__(self):
        return 1

    def __add__(self, other):
        """Add two points, or add a factor to this point's coordinates."""
        if isinstance(other, Point):
            if len(other.args) == len(self.args):
                return Point( [simplify(a + b) for a, b in zip(self, other)] )
            else:
                raise TypeError("Points must have the same number of dimensions")
        else:
            raise ValueError('Cannot add non-Point, %s, to a Point' % other)
            other = sympify(other)
            return Point([simplify(a + other) for a in self])

    def __sub__(self, other):
        """Subtract two points, or subtract a factor from this point's
        coordinates."""
        return self + (-other)

    def __mul__(self, factor):
        """Multiply point's coordinates by a factor."""
        factor = sympify(factor)
        return Point([x*factor for x in self])

    def __div__(self, divisor):
        """Divide point's coordinates by a factor."""
        divisor = sympify(divisor)
        return Point([x/divisor for x in self])

    def __neg__(self):
        """Negate the point."""
        return Point([-x for x in self])

    def __abs__(self):
        """Returns the distance between this point and the origin."""
        origin = Point([0]*len(self.args))
        return Point.distance(origin, self)
