"""Geometrical Points.

Contains
========
Point
Point2D
Point3D

"""

from __future__ import division, print_function

from sympy.core import S, sympify
from sympy.core.numbers import Number
from sympy.core.compatibility import iterable
from sympy.core.containers import Tuple
from sympy.simplify import nsimplify, simplify
from sympy.geometry.exceptions import GeometryError
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.complexes import im
from sympy.matrices import Matrix
from sympy.core.numbers import Float
from sympy.core.evaluate import global_evaluate
from sympy.core.add import Add

from .entity import GeometryEntity


class Point(GeometryEntity):
    """A point in a n-dimensional Euclidean space.

    Parameters
    ==========

    coords : sequence of n-coordinate values. In the special
    case where n=2 or 3, a Point2D or Point3D will be created
    as appropriate.

    Attributes
    ==========

    length
    origin: A `Point` representing the origin of the
        appropriately-dimensioned space.

    Raises
    ======

    TypeError
        When trying to add or subtract points with different dimensions.

    See Also
    ========

    sympy.geometry.line.Segment : Connects two Points

    Examples
    ========

    >>> from sympy.geometry import Point
    >>> from sympy.abc import x
    >>> Point(1, 2, 3)
    Point3D(1, 2, 3)
    >>> Point([1, 2])
    Point2D(1, 2)
    >>> Point(0, x)
    Point2D(0, x)

    Floats are automatically converted to Rational unless the
    evaluate flag is False:

    >>> Point(0.5, 0.25)
    Point2D(1/2, 1/4)
    >>> Point(0.5, 0.25, evaluate=False)
    Point2D(0.5, 0.25)

    """
    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', global_evaluate[0])

        if iterable(args[0]):
            if isinstance(args[0], Point) and not evaluate:
                return args[0]
            args = args[0]

        # unpack the arguments into a friendly Tuple
        # if we were already a Point, we're doing an excess
        # iteration, but we'll worry about efficiency later
        coords = Tuple(*args)
        if any(a.is_number and im(a) for a in coords):
            raise ValueError('Imaginary coordinates not permitted.')

        # Turn any Floats into rationals and simplify
        # any expressions before we instantiate
        if evaluate:
            coords = coords.xreplace(dict(
                [(f, simplify(nsimplify(f, rational=True)))
                for f in coords.atoms(Float)]))
        if len(coords) == 0:
            raise ValueError("A point must have at least one coordinate")
        if len(coords) == 2:
            return Point2D(coords, **kwargs)
        if len(coords) == 3:
            return Point3D(coords, **kwargs)

        return GeometryEntity.__new__(cls, *coords)

    is_Point = True

    def __contains__(self, item):
        return item in self.args

    def contains(self, other):
        return Point.pointify(other) == self

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
        if Point.pointify(points[0]).ambient_dimension == 2:
            return Point.is_cospherical(*points)
        return Point.are_coplanar(*points) and Point.is_cospherical(*points)

    def is_cospherical(*points):
        """Is a sequence of points cospherical?

        A set of poins is called cospherical if there exists
        sphere on which every point lies.

        Parameters
        ==========

        points : sequence of Points

        Returns
        =======

        is_cospherical : boolean
            True if points are cospherical, False otherwise.

        Notes
        =====

        Subsets of size at most 2 are trivially cospherical.
        Otherwise, the points (a,b,c,...) are cospherical
        if there exists a point x so that |a-x|=|b-x|=|c-x|=....

        Examples
        ========

        >>> from sympy.geometry import Point
        >>> p1, p2 = Point(-1, 0), Point(1, 0)
        >>> p3, p4 = Point(0, 1), Point(-1, 2)
        >>> Point.is_cospherical(p1, p2, p3)
        True
        >>> Point.is_cospherical(p1, p2, p3, p4)
        False

        """

        points = [Point.pointify(p) for p in points]
        # get unique points only
        points = list(set(points))

        if len(points) > 0 and any(p.ambient_dimension != points[0].ambient_dimension for p in points[1:]):
            raise ValueError("Cospherical points must be of the same dimension")

        if len(points) <= 2:
            return True

        # if we're here, we have at least 3 points.  This is a very common query
        # so we'll do a quick test for collinearity to avoid potential matrix operations
        if len(points) == 3:
            return (not Point.is_collinear(*points))

        # if x is the center of the sphere and r is the radius, then
        # r**2 == (a-x).dot(a-x) == (b-x).dot(b-x) == ...
        # If a == 0, then this simplifies to
        # 0 == b.dot(b) - 2*b.dot(x) == c.dot(c) - 2*c.dot(x) == ...
        # so let's test if there's a solution to this!
        origin = points[0]
        points = [p - origin for p in points[1:]]
        m, v = Matrix([p*2 for p in points]), Matrix([[p.dot(p)] for p in points])
        # solve the equation Mx=b by using the augmented matrix [M|b] and making
        # sure the b column is not a pivot
        m_rref, pivots = m.col_insert(m.cols, v).rref(simplify=True)
        return m.cols not in pivots

    def is_collinear(*args):
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

        A non-empty set set of points (a, b, c,...) is collinear if
        (b-a, c-a,...) are all scalar multiples of each other. In
        other words, the rank of the matrix [b-a|c-a|...] is one.

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
        args = [Point.pointify(p) for p in args]

        # Coincident points are irrelevant; use only unique points.
        # XXX: the current workaround for issue #11238 relies on the points
        # being in a specific order, so don't let `set` change the order
        #args = list(set(args))

        if len(args) == 0:
            return False
        if len(args) <= 2:
            return True

        return Point.affine_rank(*args) == 1

    def is_scalar_multiple(self, p2):
        """Returns whether `p1` and `p2` are scalar multiples
        of eachother.
        """
        p1, p2 = Point.pointify(self), Point.pointify(p2)

        # 2d points happen a lot, so optimize this function call
        if p1.ambient_dimension == 2:
            (x1,y1), (x2,y2) = p1.args, p2.args
            return simplify(x1*y2-x2*y1) == 0

        # if the vectors p1 and p2 are linearly dependent, then they must
        # be scalar multiples of eachother
        m = Matrix([p1.args, p2.args])
        # XXX: issue #9480 we need `simplify=True` otherwise the
        # rank may be computed incorrectly
        return m.rank(simplify=True) < 2

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

    @property
    def origin(self):
        """A point of all zeros of the same ambient dimension
        as the current point"""
        return Point([0]*len(self))

    @property
    def is_zero(self):
        """True if every coordinate is zero, False if any coordinate is not zero,
        and None if it cannot be determined."""
        nonzero = [x.is_nonzero for x in self.args]
        if any(nonzero):
            return False
        if any(x is None for x in nonzero):
            return None
        return True

    @property
    def is_nonzero(self):
        """True if any coordinate is nonzero, False if every coordinate is zero,
        and None if it cannot be determined."""
        is_zero = self.is_zero
        if is_zero is None:
            return None
        return not is_zero

    @property
    def ambient_dimension(self):
        """The dimension of the ambient space the point is in.
        I.e., if the point is in R^n, the ambient dimension
        will be n"""
        return len(self)

    @property
    def orthogonal_direction(self):
        """Returns a non-zero point in a direction
        orthogonal to `self`."""

        dim = self.ambient_dimension
        if dim == 1 and self[0] != S.Zero:
            raise ValueError("No orthogonal direction")

        # if a coordinate is zero, we can put a 1 there and zeros elsewhere
        if self[0] == S.Zero:
            return Point([1] + (dim - 1)*[0])
        if self[1] == S.Zero:
            return Point([0,1] + (dim - 2)*[0])
        # if the first two coordinates aren't zero, we can create a non-zero
        # orthogonal vector by swapping them, negating one, and padding with zeros
        return Point([-self[1], self[0]] + (dim - 2)*[0])

    @staticmethod
    def affine_rank(*points):
        """The affine rank of a set of points is the dimension
        of the smallest affine space containing all the points.
        For example, if the points lie on a line (and are not all
        the same) their affine rank is 1.  If the points lie on a plane
        but not a line, their affine rank is 2."""

        if len(points) == 0:
            raise TypeError("Must have at least one point")

        # make sure we're genuinely points
        # and translate to every point to the origin
        # XXX: issue #11238 we need to pre-simplify otherwise
        # the rank calculation will be incorrect
        origin = simplify(Point.pointify(points[0]))
        points = [Point.pointify(p, dimension=origin.ambient_dimension) - origin for p in points[1:]]

        m = Matrix([p.args for p in points])
        # XXX: issue #9480 we need `simplify=True` otherwise the
        # rank may be computed incorrectly
        return m.rank(simplify=True)

    @staticmethod
    def project(a, b):
        """Project the point `a` onto the non-zero point `b`."""
        a, b = Point.pointify(a), Point.pointify(b)
        if b.is_zero:
            raise ValueError("Cannot project to the zero vector.")
        if a.ambient_dimension != b.ambient_dimension:
            raise ValueError("Points must be of the same dimension")

        return b*(a.dot(b) / b.dot(b))

    @staticmethod
    def pointify(a, dimension=None, ensure_geometry_entity=True, ensure_point=True, **kwargs):
        """If `a` is a `Point` or a non-iterable, it is left alone.
        Otherwise, it is converted to a point with `kwargs` passed along
        to `Point()`. If `dimension` is set (default `None`), a dimension check is carried out.
        If `ensure_geometry_entity` is set (default `True`), a type check for GeometryEntity
        is carried out and the return value will be a GeometryEntity.
        If `ensure_point` is set (default `True`), a type check for Point is carried out
        and the return value will be a Point."""

        if isinstance(a, Point):
            if dimension and a.ambient_dimension != dimension:
                raise ValueError("Point {} must be of dimension {}.".format(a, dimension))
            return a

        if iterable(a):
            # in the special case of 2 and 3 dimensions, let the constructor
            # ensure the dimension is correct.  In particular, the convention
            # for Point3D is to allow being called with only two arguments
            # and silently add a zero as the last coordinate.
            if dimension == 2:
                return Point.pointify(Point2D(a, **kwargs))
            elif dimension == 3:
                return Point.pointify(Point3D(a, **kwargs))
            return Point.pointify(Point(a, **kwargs), dimension)

        if ensure_geometry_entity and not isinstance(a, GeometryEntity):
            raise TypeError("{} must be a GeometryEntity".format(a))
        if ensure_point and not isinstance(a, Point):
            raise TypeError("{} must be a Point".format(a))

        return a

    @staticmethod
    def are_coplanar(*points):
        """Do all points lie on a plane?

        A set of points is coplanar if there exists a plane
        containing all the points.  This is trivially true
        for 2d points.  For 3d and higher points, matrix operations
        are used to compute the affine rank of the set of points.

        Parameters
        ==========

        points: A sequence of points

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
        if len(points) < 0:
            return True
        return Point.affine_rank(*points) <= 2

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
        return sqrt(sum([(a - b)**2 for a, b in zip(
            self.args, p.args if isinstance(p, Point) else p)]))

    def taxicab_distance(self, p):
        """The Taxicab Distance from self to point p.

        Returns the sum of the horizontal and vertical distances to point p.

        Parameters
        ==========

        p : Point

        Returns
        =======

        taxicab_distance : The sum of the horizontal
        and vertical distances to point p.

        See Also
        ========

        sympy.geometry.Point.distance

        Examples
        ========

        >>> from sympy.geometry import Point
        >>> p1, p2 = Point(1, 1), Point(4, 5)
        >>> p1.taxicab_distance(p2)
        7

        """
        p = Point.pointify(p)
        return sum(abs(a - b) for a, b in zip(self.args, p.args))

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
        Point2D(7, 3)

        """
        p = Point.pointify(p)
        return Point([simplify((a + b)*S.Half) for a, b in zip(self.args, p.args)])

    def evalf(self, prec=None, **options):
        """Evaluate the coordinates of the point.

        This method will, where possible, create and return a new Point
        where the coordinates are evaluated as floating point numbers to
        the precision indicated (default=15).

        Parameters
        ==========

        prec : int

        Returns
        =======

        point : Point

        Examples
        ========

        >>> from sympy import Point, Rational
        >>> p1 = Point(Rational(1, 2), Rational(3, 2))
        >>> p1
        Point2D(1/2, 3/2)
        >>> p1.evalf()
        Point2D(0.5, 1.5)

        """
        coords = [x.evalf(prec, **options) for x in self.args]
        return Point(*coords, evaluate=False)

    n = evalf

    def intersection(self, other):
        """The intersection between this point and another GeometryEntity.

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
        [Point2D(0, 0)]

        """
        other = Point.pointify(other, ensure_point=False)
        if isinstance(other, Point):
            if len(self) != len(other):
                raise ValueError("Points must be of the same dimension to intersect")
            if self == other:
                return [self]
            return []

        return other.intersection(self)

    def dot(self, p2):
        """Return dot product of self with another Point."""
        p2 = Point.pointify(p2, dimension=self.ambient_dimension)
        # sum could return native python floats, we'd like to ensure the
        # return value is a sympy object, so wrap in S()
        return S(sum(a*b for a,b in zip(self, p2)))

    def equals(self, other):
        """Returns whether the coordinates of self and other agree."""
        # a point is equal to another point if all its components are equal
        if not isinstance(other, Point) or len(self.args) != len(other.args):
            return False
        return all(a.equals(b) for a,b in zip(self.args, other.args))

    def __len__(self):
        return len(self.args)

    def __iter__(self):
        return self.args.__iter__()

    def __eq__(self, other):
        if not isinstance(other, Point) or len(self.args) != len(other.args):
            return False
        return self.args == other.args

    def __hash__(self):
        return hash(self.args)

    def __getitem__(self, key):
        return self.args[key]

    def __add__(self, other):
        """Add other to self by incrementing self's coordinates by those of other.

        See Also
        ========

        sympy.geometry.entity.translate

        """

        if iterable(other):
            coords = [simplify(a + b) for a, b in zip(self, other)]
            if len(coords) != self.ambient_dimension:
                raise ValueError("Points must have the same number of dimensions")
            if all(isinstance(c, Number) for c in coords):
                return Point(coords, evaluate=False)
            return Point(coords)
        raise GeometryError("Don't know how to add {} and a Point object".format(other))

    def __sub__(self, other):
        """Subtract two points, or subtract a factor from this point's
        coordinates."""
        return self + (-x for x in other)

    def __mul__(self, factor):
        """Multiply point's coordinates by a factor."""
        factor = sympify(factor)
        coords = [simplify(x*factor) for x in self.args]
        if all(isinstance(c, Number) for c in coords):
            return Point(coords, evaluate=False)
        return Point(coords)

    def __div__(self, divisor):
        """Divide point's coordinates by a factor."""
        divisor = sympify(divisor)
        coords = [simplify(x/divisor) for x in self.args]
        if all(isinstance(c, Number) for c in coords):
            return Point(coords, evaluate=False)
        return Point(coords)

    __truediv__ = __div__

    def __neg__(self):
        """Negate the point."""
        coords = [-x for x in self.args]
        if all(isinstance(c, Number) for c in coords):
            return Point(coords, evaluate=False)
        return Point(coords)

    def __abs__(self):
        """Returns the distance between this point and the origin."""
        origin = Point([0]*len(self))
        return Point.distance(origin, self)

class Point2D(Point):
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

    TypeError
        When trying to add or subtract points with different dimensions.
        When trying to create a point with more than two dimensions.
        When `intersection` is called with object other than a Point.

    See Also
    ========

    sympy.geometry.line.Segment : Connects two Points

    Examples
    ========

    >>> from sympy.geometry import Point2D
    >>> from sympy.abc import x
    >>> Point2D(1, 2)
    Point2D(1, 2)
    >>> Point2D([1, 2])
    Point2D(1, 2)
    >>> Point2D(0, x)
    Point2D(0, x)

    Floats are automatically converted to Rational unless the
    evaluate flag is False:

    >>> Point2D(0.5, 0.25)
    Point2D(1/2, 1/4)
    >>> Point2D(0.5, 0.25, evaluate=False)
    Point2D(0.5, 0.25)

    """
    def __new__(cls, *args, **kwargs):
        eval = kwargs.get('evaluate', global_evaluate[0])
        check = True
        if isinstance(args[0], Point2D):
            if not eval:
                return args[0]
            args = args[0].args
            check = False
        else:
            if iterable(args[0]):
                args = args[0]
        coords = Tuple(*args)
        if len(coords) != 2:
            raise ValueError(
                "Only two dimensional points may be initialized with `Point2D`")
        if check:
            if any(a.is_number and im(a) for a in coords):
                raise ValueError('Imaginary args not permitted.')
        if eval:
            coords = coords.xreplace(dict(
                [(f, simplify(nsimplify(f, rational=True)))
                for f in coords.atoms(Float)]))
        return GeometryEntity.__new__(cls, *coords)

    @property
    def x(self):
        """
        Returns the X coordinate of the Point.

        Examples
        ========

        >>> from sympy import Point2D
        >>> p = Point2D(0, 1)
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

        >>> from sympy import Point2D
        >>> p = Point2D(0, 1)
        >>> p.y
        1
        """
        return self.args[1]

    @property
    def bounds(self):
        """Return a tuple (xmin, ymin, xmax, ymax) representing the bounding
        rectangle for the geometric figure.

        """
        return (self.x, self.y, self.x, self.y)

    def rotate(self, angle, pt=None):
        """Rotate ``angle`` radians counterclockwise about Point ``pt``.

        See Also
        ========

        rotate, scale

        Examples
        ========

        >>> from sympy import Point2D, pi
        >>> t = Point2D(1, 0)
        >>> t.rotate(pi/2)
        Point2D(0, 1)
        >>> t.rotate(pi/2, (2, 0))
        Point2D(2, -1)

        """
        from sympy import cos, sin, Point

        c = cos(angle)
        s = sin(angle)

        rv = self
        if pt is not None:
            pt = Point.pointify(pt)
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

        >>> from sympy import Point2D
        >>> t = Point2D(1, 1)
        >>> t.scale(2)
        Point2D(2, 1)
        >>> t.scale(2, 2)
        Point2D(2, 2)

        """
        if pt:
            pt = Point.pointify(pt)
            return self.translate(*(-pt).args).scale(x, y).translate(*pt.args)
        return Point(self.x*x, self.y*y)

    def translate(self, x=0, y=0):
        """Shift the Point by adding x and y to the coordinates of the Point.

        See Also
        ========

        rotate, scale

        Examples
        ========

        >>> from sympy import Point2D
        >>> t = Point2D(0, 1)
        >>> t.translate(2)
        Point2D(2, 1)
        >>> t.translate(2, 2)
        Point2D(2, 3)
        >>> t + Point2D(2, 2)
        Point2D(2, 3)

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
        try:
            col, row = matrix.shape
            valid_matrix = matrix.is_square and col == 3
        except AttributeError:
            # We hit this block if matrix argument is not actually a Matrix.
            valid_matrix = False
        if not valid_matrix:
            raise ValueError("The argument to the transform function must be " \
            + "a 3x3 matrix")
        x, y = self.args
        return Point(*(Matrix(1, 3, [x, y, 1])*matrix).tolist()[0][:2])

class Point3D(Point):
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

    TypeError
        When trying to add or subtract points with different dimensions.
        When `intersection` is called with object other than a Point.

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
        coords = Tuple(*args)
        if len(coords) not in (2, 3):
            raise TypeError(
                "`Point3D` must be initilized with 2 or 3 components")
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
        point = Point.pointify(point, dimension=3)
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
        point = Point.pointify(point, dimension=3)
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

        sympy.geometry.line.Line3D

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
        return Point.is_collinear(*points)

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
        try:
            col, row = matrix.shape
            valid_matrix = matrix.is_square and col == 4
        except AttributeError:
            # We hit this block if matrix argument is not actually a Matrix.
            valid_matrix = False
        if not valid_matrix:
            raise ValueError("The argument to the transform function must be " \
            + "a 4x4 matrix")
        from sympy.matrices.expressions import Transpose
        x, y, z = self.args
        m = Transpose(matrix)
        return Point3D(*(Matrix(1, 4, [x, y, z, 1])*m).tolist()[0][:3])
