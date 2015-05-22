"""Geometrical Points.

Contains
========
Point3D

"""

from __future__ import print_function, division

from sympy.core import S, sympify
from sympy.core.compatibility import iterable, range
from sympy.core.containers import Tuple
from sympy.simplify import simplify, nsimplify
from sympy.geometry.point import Point
from sympy.functions.elementary.miscellaneous import sqrt
from .entity import GeometryEntity
from sympy.matrices import Matrix
from sympy.core.numbers import Float
from sympy.core.evaluate import global_evaluate


class Point3D(GeometryEntity):
    """A point in a 3-dimensional Euclidean space.

    Parameters
    ==========

    coords : sequence of 3 coordinate values.

    Attributes
    ==========

    x
    y
    z
    length

    Raises
    ======

    NotImplementedError
        When trying to create a point other than 2 or 3 dimensions.
    TypeError
        When trying to add or subtract points with different dimensions.
        When `intersection` is called with object other than a Point.

    Notes
    =====

    Currently only 2-dimensional and 3-dimensional points are supported.

    Examples
    ========

    >>> from sympy import Point3D
    >>> from sympy.abc import x
    >>> Point3D(1, 2, 3)
    Point3D(1, 2, 3)
    >>> Point3D([1, 2, 3])
    Point3D(1, 2, 3)
    >>> Point3D(0, x, 3)
    Point3D(0, x, 3)

    Floats are automatically converted to Rational unless the
    evaluate flag is False:

    >>> Point3D(0.5, 0.25, 2)
    Point3D(1/2, 1/4, 2)
    >>> Point3D(0.5, 0.25, 3, evaluate=False)
    Point3D(0.5, 0.25, 3)

    """
    def __new__(cls, *args, **kwargs):
        eval = kwargs.get('evaluate', global_evaluate[0])
        if isinstance(args[0], (Point, Point3D)):
            if not eval:
                return args[0]
            args = args[0].args
        else:
            if iterable(args[0]):
                args = args[0]
            if len(args) not in (2, 3):
                raise TypeError(
                    "Enter a 2 or 3 dimensional point")
        coords = Tuple(*args)
        if len(coords) == 2:
            coords += (S.Zero,)
        if eval:
            coords = coords.xreplace(dict(
                [(f, simplify(nsimplify(f, rational=True)))
                for f in coords.atoms(Float)]))
        return GeometryEntity.__new__(cls, *coords)

    def __contains__(self, item):
        return item == self

    @property
    def x(self):
        """
        Returns the X coordinate of the Point.

        Examples
        ========

        >>> from sympy import Point3D
        >>> p = Point3D(0, 1, 3)
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

        >>> from sympy import Point3D
        >>> p = Point3D(0, 1, 2)
        >>> p.y
        1
        """
        return self.args[1]

    @property
    def z(self):
        """
        Returns the Z coordinate of the Point.

        Examples
        ========

        >>> from sympy import Point3D
        >>> p = Point3D(0, 1, 1)
        >>> p.z
        1
        """
        return self.args[2]

    @property
    def length(self):
        """
        Treating a Point as a Line, this returns 0 for the length of a Point.

        Examples
        ========

        >>> from sympy import Point3D
        >>> p = Point3D(0, 1, 1)
        >>> p.length
        0
        """
        return S.Zero

    def direction_ratio(self, point):
        """
        Gives the direction ratio between 2 points

        Parameters
        ==========

        p : Point3D

        Returns
        =======

        list

        Examples
        ========

        >>> from sympy import Point3D
        >>> p1 = Point3D(1, 2, 3)
        >>> p1.direction_ratio(Point3D(2, 3, 5))
        [1, 1, 2]
        """
        return [(point.x - self.x),(point.y - self.y),(point.z - self.z)]

    def direction_cosine(self, point):
        """
        Gives the direction cosine between 2 points

        Parameters
        ==========

        p : Point3D

        Returns
        =======

        list

        Examples
        ========

        >>> from sympy import Point3D
        >>> p1 = Point3D(1, 2, 3)
        >>> p1.direction_cosine(Point3D(2, 3, 5))
        [sqrt(6)/6, sqrt(6)/6, sqrt(6)/3]
        """
        a = self.direction_ratio(point)
        b = sqrt(sum(i**2 for i in a))
        return [(point.x - self.x) / b,(point.y - self.y) / b,
                (point.z - self.z) / b]

    @staticmethod
    def are_collinear(*points):
        """Is a sequence of points collinear?

        Test whether or not a set of points are collinear. Returns True if
        the set of points are collinear, or False otherwise.

        Parameters
        ==========

        points : sequence of Point

        Returns
        =======

        are_collinear : boolean

        See Also
        ========

        sympy.geometry.line3d.Line3D

        Examples
        ========

        >>> from sympy import Point3D, Matrix
        >>> from sympy.abc import x
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(1, 1, 1)
        >>> p3, p4, p5 = Point3D(2, 2, 2), Point3D(x, x, x), Point3D(1, 2, 6)
        >>> Point3D.are_collinear(p1, p2, p3, p4)
        True
        >>> Point3D.are_collinear(p1, p2, p3, p5)
        False
        """
        # Coincident points are irrelevant and can confuse this algorithm.
        # Use only unique points.
        points = list(set(points))
        if not all(isinstance(p, Point3D) for p in points):
            raise TypeError('Must pass only 3D Point objects')

        if len(points) < 2:
            return False
        if len(points) == 2:
            return True  # two points always form a line
        if len(points) == 3:
            a = (points[0].direction_cosine(points[1]))
            b = (points[0].direction_cosine(points[2]))
            c = [(-1)*i for i in b]  # opposite of b

            # a and b are either equal or opposite
            if a == b or a == c:
                return True
            else:
                return False
        # XXX Cross product is used now,
        # If the concept needs to extend to more than three
        # dimensions then another method would have to be used
        for i in range(len(points) - 2):
            pv1 = [j - k for j, k in zip(points[i].args,   \
                points[i + 1].args)]
            pv2 = [j - k for j, k in zip(points[i + 1].args,
                points[i + 2].args)]
            rank = Matrix([pv1, pv2]).rank()
            if (rank != 1):
                return False
        return True

    @staticmethod
    def are_coplanar(*points):
        """

        This function tests whether passed points are coplanar or not.
        It uses the fact that the triple scalar product of three vectors
        vanishes if the vectors are coplanar. Which means that the volume
        of the solid described by them will have to be zero for coplanarity.

        Parameters
        ==========

        A set of points 3D points

        Returns
        =======

        boolean

        Examples
        ========

        >>> from sympy import Point3D
        >>> p1 = Point3D(1, 2, 2)
        >>> p2 = Point3D(2, 7, 2)
        >>> p3 = Point3D(0, 0, 2)
        >>> p4 = Point3D(1, 1, 2)
        >>> Point3D.are_coplanar(p1, p2, p3, p4)
        True
        >>> p5 = Point3D(0, 1, 3)
        >>> Point3D.are_coplanar(p1, p2, p3, p5)
        False

        """
        from sympy.geometry.plane import Plane
        points = list(set(points))
        if len(points) < 3:
            raise ValueError('At least 3 points are needed to define a plane.')
        a, b = points[:2]
        for i, c in enumerate(points[2:]):
            try:
                p = Plane(a, b, c)
                for j in (0, 1, i):
                    points.pop(j)
                return all(p.is_coplanar(i) for i in points)
            except ValueError:
                pass
        raise ValueError('At least 3 non-collinear points needed to define plane.')

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

        >>> from sympy import Point3D
        >>> p1, p2 = Point3D(1, 1, 1), Point3D(4, 5, 0)
        >>> p1.distance(p2)
        sqrt(26)

        >>> from sympy.abc import x, y, z
        >>> p3 = Point3D(x, y, z)
        >>> p3.distance(Point3D(0, 0, 0))
        sqrt(x**2 + y**2 + z**2)

        """
        p = Point3D(p)
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

        >>> from sympy import Point3D
        >>> p1, p2 = Point3D(1, 1, 1), Point3D(13, 5, 1)
        >>> p1.midpoint(p2)
        Point3D(7, 3, 1)

        """
        p = Point3D(p)
        return Point3D([simplify((a + b)*S.Half) for a, b in
            zip(self.args, p.args)])

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

        >>> from sympy import Point3D, Rational
        >>> p1 = Point3D(Rational(1, 2), Rational(3, 2), Rational(5, 2))
        >>> p1
        Point3D(1/2, 3/2, 5/2)
        >>> p1.evalf()
        Point3D(0.5, 1.5, 2.5)

        """
        coords = [x.evalf(prec, **options) for x in self.args]
        return Point3D(*coords, evaluate=False)

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

        >>> from sympy import Point3D
        >>> p1, p2, p3 = Point3D(0, 0, 0), Point3D(1, 1, 1), Point3D(0, 0, 0)
        >>> p1.intersection(p2)
        []
        >>> p1.intersection(p3)
        [Point3D(0, 0, 0)]

        """
        if isinstance(o, Point3D):
            if self == o:
                return [self]
            return []

        return o.intersection(self)

    def scale(self, x=1, y=1, z=1, pt=None):
        """Scale the coordinates of the Point by multiplying by
        ``x`` and ``y`` after subtracting ``pt`` -- default is (0, 0) --
        and then adding ``pt`` back again (i.e. ``pt`` is the point of
        reference for the scaling).

        See Also
        ========

        translate

        Examples
        ========

        >>> from sympy import Point3D
        >>> t = Point3D(1, 1, 1)
        >>> t.scale(2)
        Point3D(2, 1, 1)
        >>> t.scale(2, 2)
        Point3D(2, 2, 1)

        """
        if pt:
            pt = Point3D(pt)
            return self.translate(*(-pt).args).scale(x, y, z).translate(*pt.args)
        return Point3D(self.x*x, self.y*y, self.z*z)

    def translate(self, x=0, y=0, z=0):
        """Shift the Point by adding x and y to the coordinates of the Point.

        See Also
        ========

        rotate, scale

        Examples
        ========

        >>> from sympy import Point3D
        >>> t = Point3D(0, 1, 1)
        >>> t.translate(2)
        Point3D(2, 1, 1)
        >>> t.translate(2, 2)
        Point3D(2, 3, 1)
        >>> t + Point3D(2, 2, 2)
        Point3D(2, 3, 3)

        """
        return Point3D(self.x + x, self.y + y, self.z + z)

    def transform(self, matrix):
        """Return the point after applying the transformation described
        by the 4x4 Matrix, ``matrix``.

        See Also
        ========
        geometry.entity.rotate
        geometry.entity.scale
        geometry.entity.translate
        """
        if isinstance(matrix, Matrix):
            from sympy.matrices.expressions import Transpose
            x, y, z = self.args
            m = Transpose(matrix)
            return Point3D(*(Matrix(1, 4, [x, y, z, 1])*m).tolist()[0][:3])
        else:
            raise ValueError("You must specify a 4x4 matrix to transform by")

    def dot(self, p2):
        """Return dot product of self with another Point."""
        p2 = Point3D(p2)
        x1, y1, z1 = self.args
        x2, y2, z2 = p2.args
        return x1*x2 + y1*y2 + z1*z2

    def equals(self, other):
        if not isinstance(other, Point3D):
            return False
        return all(a.equals(b) for a, b in zip(self.args, other.args))

    def __add__(self, other):
        """Add other to self by incrementing self's coordinates by those of
        other.

        See Also
        ========

        sympy.geometry.entity.translate

        """

        if isinstance(other, Point3D):
            if len(other.args) == len(self.args):
                return Point3D(*[simplify(a + b) for a, b in
                               zip(self.args, other.args)])
            else:
                raise TypeError(
                    "Points must have the same number of dimensions")
        else:
            raise ValueError('Cannot add non-Point, %s, to a Point' % other)

    def __sub__(self, other):
        """Subtract two points, or subtract a factor from this point's
        coordinates."""
        if isinstance(other, Point3D):
            if len(other.args) == len(self.args):
                return Point3D(*[simplify(a - b) for a, b in
                               zip(self.args, other.args)])
            else:
                raise TypeError(
                    "Points must have the same number of dimensions")
        else:
            raise ValueError('Cannot subtract non-Point, %s, to a Point'
                % other)

    def __mul__(self, factor):
        """Multiply point's coordinates by a factor."""
        factor = sympify(factor)
        return Point3D([x*factor for x in self.args])

    def __div__(self, divisor):
        """Divide point's coordinates by a factor."""
        divisor = sympify(divisor)
        return Point3D([x/divisor for x in self.args])

    __truediv__ = __div__

    def __neg__(self):
        """Negate the point."""
        return Point3D([-x for x in self.args])

    def __abs__(self):
        """Returns the distance between this point and the origin."""
        origin = Point3D([0]*len(self.args))
        return Point3D.distance(origin, self)
