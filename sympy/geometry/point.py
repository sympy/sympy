"""Geometrical Points.

Contains
========
Point

"""

from __future__ import print_function, division

from sympy.core import S, sympify
from sympy.core.compatibility import iterable
from sympy.core.containers import Tuple
from sympy.simplify import simplify, nsimplify
from sympy.geometry.exceptions import GeometryError
from sympy.functions.elementary.miscellaneous import sqrt
from .entity import GeometryEntity
from sympy.matrices import Matrix
from sympy.core.numbers import Float


class Point(GeometryEntity):
    """A point in a 2-dimensional Euclidean space.

    Parameters
    ==========

    coords : sequence of 2 coordinate values.

    Attributes
    ==========

    x
    y
    length

    Raises
    ======

    NotImplementedError
        When trying to create a point with more than two dimensions.
        When `intersection` is called with object other than a Point.
    TypeError
        When trying to add or subtract points with different dimensions.

    Notes
    =====

    Currently only 2-dimensional points are supported.

    See Also
    ========

    sympy.geometry.line.Segment : Connects two Points

    Examples
    ========

    >>> from sympy.geometry import Point
    >>> from sympy.abc import x
    >>> Point(1, 2)
    Point(1, 2)
    >>> Point([1, 2])
    Point(1, 2)
    >>> Point(0, x)
    Point(0, x)

    Floats are automatically converted to Rational unless the
    evaluate flag is False:

    >>> Point(0.5, 0.25)
    Point(1/2, 1/4)
    >>> Point(0.5, 0.25, evaluate=False)
    Point(0.5, 0.25)

    """

    def __new__(cls, *args, **kwargs):
        if iterable(args[0]):
            args = args[0]
        elif isinstance(args[0], Point):
            args = args[0].args
        coords = Tuple(*args)

        if len(coords) != 2:
            raise NotImplementedError(
                "Only two dimensional points currently supported")
        if kwargs.get('evaluate', True):
            coords = coords.xreplace(dict(
                [(f, simplify(nsimplify(f, rational=True)))
                for f in coords.atoms(Float)]))

        return GeometryEntity.__new__(cls, *coords)

    def __hash__(self):
        return super(Point, self).__hash__()

    def __eq__(self, other):
        ts, to = type(self), type(other)
        if ts is not to:
            return False
        return self.args == other.args

    def __lt__(self, other):
        return self.args < other.args

    def __contains__(self, item):
        return item == self

    @property
    def x(self):
        """
        Returns the X coordinate of the Point.

        Examples
        ========

        >>> from sympy import Point
        >>> p = Point(0, 1)
        >>> p.x
        0
        """
        return self.args[0]

    @property
    def y(self):
        """
        Returns the Y coordinate of the Point.

        Examples
        ========

        >>> from sympy import Point
        >>> p = Point(0, 1)
        >>> p.y
        1
        """
        return self.args[1]

    @property
    def length(self):
        """
        Treating a Point as a Line, this returns 0 for the length of a Point.

        Examples
        ========

        >>> from sympy import Point
        >>> p = Point(0, 1)
        >>> p.length
        0
        """
        return S.Zero

    def is_collinear(*points):
        """Is a sequence of points collinear?

        Test whether or not a set of points are collinear. Returns True if
        the set of points are collinear, or False otherwise.

        Parameters
        ==========

        points : sequence of Point

        Returns
        =======

        is_collinear : boolean

        Notes
        =====

        Slope is preserved everywhere on a line, so the slope between
        any two points on the line should be the same. Take the first
        two points, p1 and p2, and create a translated point v1
        with p1 as the origin. Now for every other point we create
        a translated point, vi with p1 also as the origin. Note that
        these translations preserve slope since everything is
        consistently translated to a new origin of p1. Since slope
        is preserved then we have the following equality:

              * v1_slope = vi_slope
              * v1.y/v1.x = vi.y/vi.x (due to translation)
              * v1.y*vi.x = vi.y*v1.x
              * v1.y*vi.x - vi.y*v1.x = 0           (*)

        Hence, if we have a vi such that the equality in (*) is False
        then the points are not collinear. We do this test for every
        point in the list, and if all pass then they are collinear.

        See Also
        ========

        sympy.geometry.line.Line

        Examples
        ========

        >>> from sympy import Point
        >>> from sympy.abc import x
        >>> p1, p2 = Point(0, 0), Point(1, 1)
        >>> p3, p4, p5 = Point(2, 2), Point(x, x), Point(1, 2)
        >>> Point.is_collinear(p1, p2, p3, p4)
        True
        >>> Point.is_collinear(p1, p2, p3, p5)
        False

        """
        # Coincident points are irrelevant and can confuse this algorithm.
        # Use only unique points.
        points = list(set(points))
        if len(points) == 0:
            return False
        if len(points) <= 2:
            return True  # two points always form a line
        points = [Point(a) for a in points]

        # XXX Cross product is used now, but that only extends to three
        #     dimensions. If the concept needs to extend to greater
        #     dimensions then another method would have to be used
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

    def is_concyclic(*points):
        """Is a sequence of points concyclic?

        Test whether or not a sequence of points are concyclic (i.e., they lie
        on a circle).

        Parameters
        ==========

        points : sequence of Points

        Returns
        =======

        is_concyclic : boolean
            True if points are concyclic, False otherwise.

        See Also
        ========

        sympy.geometry.ellipse.Circle

        Notes
        =====

        No points are not considered to be concyclic. One or two points
        are definitely concyclic and three points are conyclic iff they
        are not collinear.

        For more than three points, create a circle from the first three
        points. If the circle cannot be created (i.e., they are collinear)
        then all of the points cannot be concyclic. If the circle is created
        successfully then simply check the remaining points for containment
        in the circle.

        Examples
        ========

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
            from .ellipse import Circle
            c = Circle(points[0], points[1], points[2])
            for point in points[3:]:
                if point not in c:
                    return False
            return True
        except GeometryError:
            # Circle could not be created, because of collinearity of the
            # three points passed in, hence they are not concyclic.
            return False

    def distance(self, p):
        """The Euclidean distance from self to point p.

        Parameters
        ==========

        p : Point

        Returns
        =======

        distance : number or symbolic expression.

        See Also
        ========

        sympy.geometry.line.Segment.length

        Examples
        ========

        >>> from sympy.geometry import Point
        >>> p1, p2 = Point(1, 1), Point(4, 5)
        >>> p1.distance(p2)
        5

        >>> from sympy.abc import x, y
        >>> p3 = Point(x, y)
        >>> p3.distance(Point(0, 0))
        sqrt(x**2 + y**2)

        """
        p = Point(p)
        return sqrt(sum([(a - b)**2 for a, b in zip(self.args, p.args)]))

    def midpoint(self, p):
        """The midpoint between self and point p.

        Parameters
        ==========

        p : Point

        Returns
        =======

        midpoint : Point

        See Also
        ========

        sympy.geometry.line.Segment.midpoint

        Examples
        ========

        >>> from sympy.geometry import Point
        >>> p1, p2 = Point(1, 1), Point(13, 5)
        >>> p1.midpoint(p2)
        Point(7, 3)

        """
        return Point([simplify((a + b)*S.Half) for a, b in zip(self.args, p.args)])

    def evalf(self, prec=None, **options):
        """Evaluate the coordinates of the point.

        This method will, where possible, create and return a new Point
        where the coordinates are evaluated as floating point numbers to
        the precision indicated (default=15).

        Returns
        =======

        point : Point

        Examples
        ========

        >>> from sympy import Point, Rational
        >>> p1 = Point(Rational(1, 2), Rational(3, 2))
        >>> p1
        Point(1/2, 3/2)
        >>> p1.evalf()
        Point(0.5, 1.5)

        """
        if prec is None:
            coords = [x.evalf(**options) for x in self.args]
        else:
            coords = [x.evalf(prec, **options) for x in self.args]
        return Point(*coords, evaluate=False)

    n = evalf

    def intersection(self, o):
        """The intersection between this point and another point.

        Parameters
        ==========

        other : Point

        Returns
        =======

        intersection : list of Points

        Notes
        =====

        The return value will either be an empty list if there is no
        intersection, otherwise it will contain this point.

        Examples
        ========

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

    def rotate(self, angle, pt=None):
        """Rotate ``angle`` radians counterclockwise about Point ``pt``.

        See Also
        ========

        rotate, scale

        Examples
        ========

        >>> from sympy import Point, pi
        >>> t = Point(1, 0)
        >>> t.rotate(pi/2)
        Point(0, 1)
        >>> t.rotate(pi/2, (2, 0))
        Point(2, -1)

        """
        from sympy import cos, sin, Point

        c = cos(angle)
        s = sin(angle)

        rv = self
        if pt is not None:
            pt = Point(pt)
            rv -= pt
        x, y = rv.args
        rv = Point(c*x - s*y, s*x + c*y)
        if pt is not None:
            rv += pt
        return rv

    def scale(self, x=1, y=1, pt=None):
        """Scale the coordinates of the Point by multiplying by
        ``x`` and ``y`` after subtracting ``pt`` -- default is (0, 0) --
        and then adding ``pt`` back again (i.e. ``pt`` is the point of
        reference for the scaling).

        See Also
        ========

        rotate, translate

        Examples
        ========

        >>> from sympy import Point
        >>> t = Point(1, 1)
        >>> t.scale(2)
        Point(2, 1)
        >>> t.scale(2, 2)
        Point(2, 2)

        """
        if pt:
            pt = Point(pt)
            return self.translate(*(-pt).args).scale(x, y).translate(*pt.args)
        return Point(self.x*x, self.y*y)

    def translate(self, x=0, y=0):
        """Shift the Point by adding x and y to the coordinates of the Point.

        See Also
        ========

        rotate, scale

        Examples
        ========

        >>> from sympy import Point
        >>> t = Point(0, 1)
        >>> t.translate(2)
        Point(2, 1)
        >>> t.translate(2, 2)
        Point(2, 3)
        >>> t + Point(2, 2)
        Point(2, 3)

        """
        return Point(self.x + x, self.y + y)

    def transform(self, matrix):
        """Return the point after applying the transformation described
        by the 3x3 Matrix, ``matrix``.

        See Also
        ========
        geometry.entity.rotate
        geometry.entity.scale
        geometry.entity.translate
        """
        x, y = self.args
        return Point(*(Matrix(1, 3, [x, y, 1])*matrix).tolist()[0][:2])

    def dot(self, p2):
        """Return dot product of self with another Point."""
        p2 = Point(p2)
        x1, y1 = self.args
        x2, y2 = p2.args
        return x1*x2 + y1*y2

    def __add__(self, other):
        """Add other to self by incrementing self's coordinates by those of other.

        See Also
        ========

        sympy.geometry.entity.translate

        """

        if isinstance(other, Point):
            if len(other.args) == len(self.args):
                return Point(*[simplify(a + b) for a, b in
                               zip(self.args, other.args)])
            else:
                raise TypeError(
                    "Points must have the same number of dimensions")
        else:
            raise ValueError('Cannot add non-Point, %s, to a Point' % other)

    def __sub__(self, other):
        """Subtract two points, or subtract a factor from this point's
        coordinates."""
        return self + (-other)

    def __mul__(self, factor):
        """Multiply point's coordinates by a factor."""
        factor = sympify(factor)
        return Point([x*factor for x in self.args])

    def __div__(self, divisor):
        """Divide point's coordinates by a factor."""
        divisor = sympify(divisor)
        return Point([x/divisor for x in self.args])

    __truediv__ = __div__

    def __neg__(self):
        """Negate the point."""
        return Point([-x for x in self.args])

    def __abs__(self):
        """Returns the distance between this point and the origin."""
        origin = Point([0]*len(self.args))
        return Point.distance(origin, self)
