"""Geometrical Points.

Contains
========
Point

"""

from sympy.core import S, sympify
from sympy.core.compatibility import iterable
from sympy.core.containers import Tuple
from sympy.simplify import simplify
from sympy.geometry.exceptions import GeometryError
from sympy.functions.elementary.miscellaneous import sqrt
from entity import GeometryEntity
from sympy.core.numbers import Float
from sympy.core.exprtools import factor_terms


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

    """

    def __new__(cls, *args, **kwargs):
        if iterable(args[0]):
            coords = Tuple(*args[0])
        elif isinstance(args[0], Point):
            coords = args[0].args
        else:
            coords = Tuple(*args)

        if len(coords) != 2:
            raise NotImplementedError("Only two dimensional points currently supported")

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

    def is_collinear(*points, **kwargs):
        """Is a sequence of points collinear?

        Test whether or not a set of points is collinear (i.e the points lie
        on a straight line). Return True if the points are collinear for all
        values of the free variables, False if they are not and None if it
        is not known.
        
        If the coordinates are floating point values, use `tol` as
        tolerance.
        
        If ``conditions`` is True return a list of conditions equivalent to 
        the collinearity of the points or a boolean.

        Parameters
        ==========

        points : sequence of Point
        tol : tolerance for floating point comparisons
        conditions : return conditions for free variables

        Returns
        =======

        is_collinear : boolean, or list of conditions for collinearity

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

        >>> from sympy import Point, Symbol, sqrt
        >>> from sympy.abc import x
        >>> p1, p2 = Point(0, 0), Point(1, 1)
        >>> p3, p4, p5 = Point(2, 2), Point(x, x), Point(1, 2)
        >>> Point.is_collinear(p1, p2, p3, p4)
        True
        >>> Point.is_collinear(p1, p2, p3, p5)
        False
        
        >>> p_2_x = Point(2, x)
        >>> Point.is_collinear(p1, p2, p_2_x)
        False
        >>> Point.is_collinear(p1, p2, p_2_x, conditions=True)
        [x - 2]
        >>> y = Symbol('y', real=True)
        >>> p_y_3 = Point(y, 3)
        >>> Point.is_collinear(p1, p2, p_2_x, p_y_3, conditions=True)
        [x - 2, -y + 3]
        
        >>> p3 = Point(2, 2.0001)
        >>> Point.is_collinear(p1, p2, p3, tol=0.00001)
        False
        >>> Point.is_collinear(p1, p2, p3, tol=0.001)
        True
        
        >>> p1 = Point(0, 0)
        >>> p2 = Point(-25*(-5*sqrt(2)/2 + 5)**2 + (-25*sqrt(2)/2 + 25)**2, 0) # == 0
        >>> p3 = Point(1, 0)
        >>> p4 = Point(0, 1)
        >>> Point.is_collinear(p1, p2, p3, p4)
        False
        """
        conditions = kwargs.pop("conditions", False)
        tol = kwargs.pop("tol", Float("1E-15"))
        
        if len(kwargs) > 0:
            raise ValueError("Keyword arguments not understood: " + str(kwargs))
        
        if len(points) == 0:
            return False
        if len(points) <= 2:
            return True # two points always form a line
        points = [Point(a) for a in points]
        first = points.pop(0)

        
        points = [p - first for p in points]
        first = points.pop(0)
        
        # drop first Points while they are 0
        while True:
            if len(points) == 0:
                return True
            all_zero = True
            for x in first.args:
                if x.is_Float:
                    if not x.epsilon_eq(0, tol):
                        all_zero = False
                        break
                elif x.is_Number:
                    if not x.is_zero:
                        all_zero = False
                        break
                elif not x.equals(0):
                    all_zero = False
                    break
            if all_zero:
                first = points.pop(0)
            else:
                break
        
            
        failings = []
        

        for p in points:
            # for n-dimensional Points this could be used
            #collin = first.dot(p) ** 2 - first.dot(first) * p.dot(p)
            
            # but for 2 dimensions this is better
            collin = first.x * p.y - first.y * p.x
            
            if collin.is_Float:
                if collin.epsilon_eq(0, epsilon=tol):
                    continue
                else:
                    return False
            
            is_const = collin.is_constant()
            
            if is_const is True:
                #speed issue: this calls is_constant again 
                failing_exp = collin.equals(0, failing_expression=True)
                if failing_exp is True:
                    continue
                elif failing_exp is False:
                    return False
                else:
                    failings.append(failing_exp)
            elif is_const is False:
                failing_exp = collin.equals(0, failing_expression=True)
                if failing_exp is True:
                    continue
                elif failing_exp is False:
                    if conditions:
                        failings.append(collin)
                    else:
                        return False
                else:
                    if conditions:
                        failings.append(collin)
                    else:
                        return None
            else:
                if conditions:
                    failings.append(collin)
                else:
                    return None
                

        if len(failings) == 0:
            return True
        
        if conditions:
            return failings


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
            from ellipse import Circle
            c = Circle(points[0], points[1], points[2])
            for point in points[3:]:
                if point not in c:
                    return False
            return True
        except GeometryError:
            # Circle could not be created, because of collinearity of the
            # three points passed in, hence they are not concyclic.
            return False

#       """
#       # This code is from Maple
#       def f(u):
#           dd = u[0]**2 + u[1]**2 + 1
#           u1 = 2*u[0] / dd
#           u2 = 2*u[1] / dd
#           u3 = (dd - 2) / dd
#           return u1,u2,u3

#       u1,u2,u3 = f(points[0])
#       v1,v2,v3 = f(points[1])
#       w1,w2,w3 = f(points[2])
#       p = [v1 - u1, v2 - u2, v3 - u3]
#       q = [w1 - u1, w2 - u2, w3 - u3]
#       r = [p[1]*q[2] - p[2]*q[1], p[2]*q[0] - p[0]*q[2], p[0]*q[1] - p[1]*q[0]]
#       for ind in xrange(3, len(points)):
#           s1,s2,s3 = f(points[ind])
#           test = simplify(r[0]*(s1-u1) + r[1]*(s2-u2) + r[2]*(s3-u3))
#           if test != 0:
#               return False
#       return True
#       """

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
            return Point(*[x.evalf(**options) for x in self.args])
        else:
            return Point(*[x.evalf(prec, **options) for x in self.args])

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
                raise TypeError("Points must have the same number of dimensions")
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
