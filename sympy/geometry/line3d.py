"""Line-like geometrical entities.

Contains
========
Base
Line3D
Ray3D
Segment3D

"""
from __future__ import print_function, division

from sympy.core import S, C, sympify, Dummy
from sympy.functions.elementary.trigonometric import _pi_coeff as pi_coeff, \
    sqrt
from sympy.core.logic import fuzzy_and
from sympy.core.exprtools import factor_terms
from sympy.simplify.simplify import simplify
from sympy.solvers import solve
from sympy.geometry.exceptions import GeometryError
from .entity import GeometryEntity
from .point3d import Point3D
from .util import _symbol

class Base(GeometryEntity):
    """An abstract base class for all linear entities (line, ray and segment)
    in a 3-dimensional Euclidean space.

    Attributes
    ==========

    p1
    p2
    direction_ratio
    points

    Notes
    =====

    This is an abstract class and is not meant to be instantiated.
    Subclasses should implement the following methods:

    * __eq__
    * contains
    """

    def __new__(cls, p1, p2, **kwargs):
        p1 = Point3D(p1)
        p2 = Point3D(p2)
        if p1 == p2:
            # if it makes sense to return a Point, handle in subclass
            raise ValueError(
                "%s.__new__ requires two unique Points." % cls.__name__)

        return GeometryEntity.__new__(cls, p1, p2, **kwargs)

    @property
    def p1(self):
        """The first defining point of a linear entity.

        See Also
        ========

        sympy.geometry.point3d.Point3D

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point(0, 0, 0), Point(5, 3, 1)
        >>> l = Line3D(p1, p2)
        >>> l.p1
        Point3D(0, 0, 0)

        """
        return self.args[0]

    @property
    def p2(self):
        """The second defining point of a linear entity.

        See Also
        ========

        sympy.geometry.point3d.Point3D

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(5, 3, 1)
        >>> l = Line3D(p1, p2)
        >>> l.p2
        Point3D(5, 3, 1)

        """
        return self.args[1]        
        
    @property
    def direction_ratio(self):
        """The direction ratio of a given line in 3D.

        See Also
        ========

        sympy.geometry.line.Line.equation

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(5, 3, 1)
        >>> l = Line3D(p1, p2)
        >>> l.direction_ratio
        [5, 3, 1]
        """
        p1, p2 = self.points
        return p1.direction_ratio(p2)

    @property
    def direction_cosine(self):
        """The direction ratio of a given line in 3D.

        See Also
        ========

        sympy.geometry.line.Line.equation

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(5, 3, 1)
        >>> l = Line3D(p1, p2)
        >>> l.direction_cosine
        [sqrt(35)/7, 3*sqrt(35)/35, sqrt(35)/35]
        """
        p1, p2 = self.points
        return p1.direction_cosine(p2)
            
    @property
    def length(self):
        """
        The length of the line.

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(3, 5, 1)
        >>> l1 = Line3D(p1, p2)
        >>> l1.length
        oo
        """
        return S.Infinity

    @property
    def points(self):
        """The two points used to define this linear entity.

        Returns
        =======

        points : tuple of Points

        See Also
        ========

        sympy.geometry.point3d.Point3D

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(5, 11, 1)
        >>> l1 = Line3D(p1, p2)
        >>> l1.points
        (Point3D(0, 0, 0), Point3D(5, 11, 1))

        """
        return (self.p1, self.p2)

    def is_parallel(l1, l2):
        """Are two linear entities parallel?

        Parameters
        ==========

        l1 : LinearEntity
        l2 : LinearEntity

        Returns
        =======

        True : if l1 and l2 are parallel,
        False : otherwise.

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(3, 4, 5)
        >>> p3, p4 = Point3D(2, 1, 1), Point3D(8, 9, 11)
        >>> l1, l2 = Line3D(p1, p2), Line3D(p3, p4)
        >>> Line3D.is_parallel(l1, l2)
        True

        >>> p5 = Point3D(6, 6, 6)
        >>> l3 = Line3D(p3, p5)
        >>> Line3D.is_parallel(l1, l3)
        False

        """
        if l1 == l2:
            raise('Enter two distinct lines')
        return l1.direction_cosine == l2.direction_cosine
        
    def arbitrary_point(self, parameter='t'):
        """A parameterized point on the Line.

        Parameters
        ==========

        parameter : str, optional
            The name of the parameter which will be used for the parametric
            point. The default value is 't'. When this parameter is 0, the
            first point used to define the line will be returned, and when
            it is 1 the second point will be returned.

        Returns
        =======

        point : Point

        Raises
        ======

        ValueError
            When ``parameter`` already appears in the Line's definition.

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy import Point, Line
        >>> p1, p2 = Point(1, 0), Point(5, 3)
        >>> l1 = Line(p1, p2)
        >>> l1.arbitrary_point()
        Point(4*t + 1, 3*t)

        """
        t = _symbol(parameter)
        if t.name in (f.name for f in self.free_symbols):
            raise ValueError('Symbol %s already appears in object '
            'and cannot be used as a parameter.' % t.name)
        x = simplify(self.p1.x + t*(self.p2.x - self.p1.x))
        y = simplify(self.p1.y + t*(self.p2.y - self.p1.y))
        z = simplify(self.p1.z + t*(self.p2.z - self.p1.z))
        return Point3D(x, y, z)

    def is_similar(self, other):
        """
        Return True if self and other are contained in the same line.

        Examples
        ========

        >>> from sympy import Point, Line
        >>> p1, p2, p3 = Point(0, 1), Point(3, 4), Point(2, 3)
        >>> l1 = Line(p1, p2)
        >>> l2 = Line(p1, p3)
        >>> l1.is_similar(l2)
        True
        """
        return 
        
    def __contains__(self, other):
        """Return a definitive answer or else raise an error if it cannot
        be determined that other is on the boundaries of self."""
        result = self.contains(other)

        if result is not None:
            return result
        else:
            raise Undecidable(
                "can't decide whether '%s' contains '%s'" % (self, other))

    def contains(self, other):
        """Subclasses should implement this method and should return
            True if other is on the boundaries of self;
            False if not on the boundaries of self;
            None if a determination cannot be made."""
        raise NotImplementedError()

    def __eq__(self, other):
        """Subclasses should implement this method."""
        raise NotImplementedError()

    def __hash__(self):
        return super(Base, self).__hash__()    
        
class Line3D(Base):
    """An infinite 3D line in space.

    A line is declared with two distinct points or a point and direction_ratio
    as defined using keyword `direction_ratio`.

    Parameters
    ==========

    p1 : Point3D
    pt : Point3D
    direction_ratio : list

    See Also
    ========

    sympy.geometry.point3d.Point3D

    Examples
    ========

    >>> import sympy
    >>> from sympy import Point
    >>> from sympy.abc import L
    >>> from sympy.geometry import Line, Segment
    >>> L = Line(Point(2,3), Point(3,5))
    >>> L
    Line(Point(2, 3), Point(3, 5))
    >>> L.points
    (Point(2, 3), Point(3, 5))
    >>> L.equation()
    -2*x + y + 1
    >>> L.coefficients
    (-2, 1, 1)

    Instantiate with keyword ``slope``:

    >>> Line(Point(0, 0), slope=0)
    Line(Point(0, 0), Point(1, 0))

    Instantiate with another linear object

    >>> s = Segment((0, 0), (0, 1))
    >>> Line(s).equation()
    x
    """

    def __new__(cls, p1, pt=None, direction_ratio=[], **kwargs):
        if isinstance(p1, Base):
            p1, pt = p1.args
        else:
            p1 = Point3D(p1)
        if pt is not None and len(direction_ratio) == 0:
            try:
                pt = Point3D(pt)
            except NotImplementedError:
                raise ValueError('The 2nd argument was not a valid Point. '
                'If it was the direction_ratio of the desired line, enter it' 
                'with keyword "direction_ratio".')
        elif len(direction_ratio) == 3 and pt is None:
            pt = Point3D(p1.x + direction_ratio[0], p1.y + direction_ratio[1],
                         p1.z + direction_ratio[2])
        else:
            raise ValueError('A 2nd Point or keyword "direction_ratio" must' 
            'be used.')
            
        return Base.__new__(cls, p1, pt, **kwargs)

    def equation(self, x='x', y='y', z='z', k='k'):
        """The equation of the line in 3D

        Parameters
        ==========

        x : str, optional
            The name to use for the x-axis, default value is 'x'.
        y : str, optional
            The name to use for the y-axis, default value is 'y'.
        z : str, optional
            The name to use for the x-axis, default value is 'z'.    

        Returns
        =======

        equation : sympy expression

        See Also
        ========

        LinearEntity.coefficients

        Examples
        ========

        >>> from sympy import Point, Line
        >>> p1, p2 = Point(1, 0), Point(5, 3)
        >>> l1 = Line(p1, p2)
        >>> l1.equation()
        -3*x + 4*y + 3

        """
        x, y, z, k = _symbol(x), _symbol(y), _symbol(z), _symbol(k)
        p1, p2 = self.points
        a = p1.direction_ratio(p2)
        return (((x - p1.x)/a[0]) - k, ((y - p1.y)/a[1]) - k, 
                ((z - p1.z)/a[2]) - k)

    def __eq__(self, other):
        """Return True if other is equal to this Line, or False otherwise."""
        if not isinstance(other, Line3D):
            return False
        return Point3D.is_collinear(self.p1, self.p2, other.p1, other.p2)

    def __hash__(self):
        return super(Line3D, self).__hash__()
        
class Segment3D(Base):
    """A undirected line segment in space.

    Parameters
    ==========

    p1 : Point
    p2 : Point

    Attributes
    ==========

    length : number or sympy expression
    midpoint : Point

    See Also
    ========

    sympy.geometry.point.Point, Line

    Notes
    =====

    At the moment only segments in a 2D space can be declared, because
    Points can be defined only for 2D spaces.

    Examples
    ========

    >>> import sympy
    >>> from sympy import Point
    >>> from sympy.abc import s
    >>> from sympy.geometry import Segment
    >>> Segment((1, 0), (1, 1)) # tuples are interpreted as pts
    Segment(Point(1, 0), Point(1, 1))
    >>> s = Segment(Point(4, 3), Point(1, 1))
    >>> s
    Segment(Point(1, 1), Point(4, 3))
    >>> s.points
    (Point(1, 1), Point(4, 3))
    >>> s.slope
    2/3
    >>> s.length
    sqrt(13)
    >>> s.midpoint
    Point(5/2, 2)

    """

    def __new__(cls, p1, p2, **kwargs):
        # Reorder the two points under the following ordering:
        #   if p1.x != p2.x then p1.x < p2.x
        #   if p1.x == p2.x then p1.y < p2.y
        #   The z-coordinate will not come into picture while ordering
        p1 = Point3D(p1)
        p2 = Point3D(p2)
        if p1 == p2:
            return Point3D(p1)
        if (p1.x > p2.x) == True:
            p1, p2 = p2, p1
        elif (p1.x == p2.x) == True and (p1.y > p2.y) == True:
            p1, p2 = p2, p1
        return Base.__new__(cls, p1, p2, **kwargs)
        
    @property
    def length(self):
        """The length of the line segment.

        See Also
        ========

        sympy.geometry.point.Point.distance

        Examples
        ========

        >>> from sympy import Point, Segment
        >>> p1, p2 = Point(0, 0), Point(4, 3)
        >>> s1 = Segment(p1, p2)
        >>> s1.length
        5

        """
        return Point3D.distance(self.p1, self.p2)

    @property
    def midpoint(self):
        """The midpoint of the line segment.

        See Also
        ========

        sympy.geometry.point.Point.midpoint

        Examples
        ========

        >>> from sympy import Point, Segment
        >>> p1, p2 = Point(0, 0), Point(4, 3)
        >>> s1 = Segment(p1, p2)
        >>> s1.midpoint
        Point(2, 3/2)

        """
        return Point3D.midpoint(self.p1, self.p2)    