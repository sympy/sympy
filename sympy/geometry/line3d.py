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
    direction_cosine
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
        
    def is_perpendicular(l1, l2):
        """Are two linear entities perpendicular?

        Parameters
        ==========

        l1 : LinearEntity
        l2 : LinearEntity

        Returns
        =======

        True : if l1 and l2 are perpendicular,
        False : otherwise.

        See Also
        ========

        direction_ratio

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2, p3 = Point3D(0, 0, 0), Point3D(1, 1, 1), Point3D(-1, 2, 0)
        >>> l1, l2 = Line3D(p1, p2), Line(p2, p3)
        >>> l1.is_perpendicular(l2)
        True

        >>> p4 = Point3D(5, 3, 7)
        >>> l3 = Line(p1, p4)
        >>> l1.is_perpendicular(l3)
        False

        """
        a = sum([i*j for i, j in zip(l1.direction_ratio, l2.direction_ratio)])
        if a == 0:
            return True
        else:
            return False
            
    def angle_between(l1, l2):
        """The angle formed between the two linear entities.

        Parameters
        ==========

        l1 : LinearEntity
        l2 : LinearEntity

        Returns
        =======

        angle : angle in radians

        Notes
        =====

        From the dot product of vectors v1 and v2 it is known that:

            ``dot(v1, v2) = |v1|*|v2|*cos(A)``

        where A is the angle formed between the two vectors. We can
        get the directional vectors of the two lines and readily
        find the angle between the two using the above formula.

        See Also
        ========

        is_perpendicular

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2, p3 = Point3D(0, 0, 0), Point3D(1, 1, 1), Point3D(-1, 2, 0)
        >>> l1, l2 = Line3D(p1, p2), Line(p2, p3)
        >>> l1.angle_between(l2)
        pi/2

        """
        v1 = l1.p2 - l1.p1
        v2 = l2.p2 - l2.p1
        return C.acos(v1.dot(v2)/(abs(v1)*abs(v2)))
        
    def parallel_line(self, p):
        """Create a new Line parallel to this linear entity which passes
        through the point `p`.

        Parameters
        ==========

        p : Point3D

        Returns
        =======

        line : Line3D

        See Also
        ========

        is_parallel

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2, p3 = Point3D(0, 0, 0), Point3D(2, 3, 4), Point3D(-2, 2, 0)
        >>> l1 = Line3D(p1, p2)
        >>> l2 = l1.parallel_line(p3)
        >>> p3 in l2
        True
        >>> l1.is_parallel(l2)
        True

        """
        d = self.direction_ratio
        return Line3D(p, direction_ratio=d)
        
    def perpendicular_line(self, p):
        """Create a new Line perpendicular to this linear entity which passes
        through the point `p`.

        Parameters
        ==========

        p : Point3D

        Returns
        =======

        line : Line3D

        See Also
        ========

        is_perpendicular, perpendicular_segment

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2, p3 = Point3D(0, 0, 0), Point3D(2, 3, 4), Point3D(-2, 2, 0)
        >>> l1 = Line3D(p1, p2)
        >>> l2 = l1.perpendicular_line(p3)
        >>> p3 in l2
        True
        >>> l1.is_perpendicular(l2)
        True

        """
        a = self.arbitrary_point()
        b = [i - j for i, j in zip(p.args, a.args)]
        c = sum([i*j for i, j in zip(b, self.direction_ratio)])
        d = solve(c)
        e = a.subs(a.args[0], d.pop())
        return Line3D(p, e)

    def perpendicular_segment(self, p):
        """Create a perpendicular line segment from `p` to this line.

        The enpoints of the segment are ``p`` and the closest point in
        the line containing self. (If self is not a line, the point might
        not be in self.)

        Parameters
        ==========

        p : Point3D

        Returns
        =======

        segment : Segment3D

        Notes
        =====

        Returns `p` itself if `p` is on this linear entity.

        See Also
        ========

        perpendicular_line

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2, p3 = Point3D(0, 0, 0), Point3D(1, 1, 1), Point3D(0, 2, 0)
        >>> l1 = Line3D(p1, p2)
        >>> s1 = l1.perpendicular_segment(p3)
        >>> l1.is_perpendicular(s1)
        True
        >>> p3 in s1
        True
        >>> l1.perpendicular_segment(Point3D(4, 0, 0))
        Segment3D(Point3D(4/3, 4/3, 4/3), Point3D(4, 0, 0))

        """
        a = self.arbitrary_point()
        b = [i - j for i, j in zip(p.args, a.args)]
        c = sum([i*j for i, j in zip(b, self.direction_ratio)])
        d = solve(c)
        e = a.subs(a.free_symbols.pop(), d.pop())
        return Segment3D(p, e)

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

        point : Point3D

        Raises
        ======

        ValueError
            When ``parameter`` already appears in the Line's definition.

        See Also
        ========

        sympy.geometry.point3d.Point3D

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(1, 0, 0), Point3D(5, 3, 1)
        >>> l1 = Line3D(p1, p2)
        >>> l1.arbitrary_point()
        Point3D(4*t + 1, 3*t, t)

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

        >>> from sympy import Point3D, Line3D
        >>> p1, p2, p3 = Point3D(0, 0, 0), Point3D(1, 1, 1), Point3D(2, 2, 2)
        >>> l1 = Line3D(p1, p2)
        >>> l2 = Line3D(p1, p3)
        >>> l1.is_similar(l2)
        True
        """
        if isinstance(other, Line3D):
            if self.direction_cosine == other.direction_cosine and other.p1 in self:
                return True
            else:
                return False
        raise NotImplementedError()
        
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
    >>> from sympy import Point3D
    >>> from sympy.abc import L
    >>> from sympy.geometry import Line3D, Segment3D
    >>> L = Line3D(Point3D(2, 3, 4), Point3D(3, 5, 1))
    >>> L
    Line3D(Point3D(2, 3, 4), Point3D(3, 5, 1))
    >>> L.points
    (Point3D(2, 3, 4), Point3D(3, 5, 1))
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
        
    def plot_interval(self, parameter='t'):
        """The plot interval for the default geometric plot of line. Gives
        values that will produce a line that is +/- 5 units long (where a
        unit is the distance between the two points that define the line).

        Parameters
        ==========

        parameter : str, optional
            Default value is 't'.

        Returns
        =======

        plot_interval : list (plot interval)
            [parameter, lower_bound, upper_bound]

        Examples
        ========

        >>> from sympy import Point, Line
        >>> p1, p2 = Point(0, 0), Point(5, 3)
        >>> l1 = Line(p1, p2)
        >>> l1.plot_interval()
        [t, -5, 5]

        """
        t = _symbol(parameter)
        return [t, -5, 5]    

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

        equation : tuple

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(1, 0, 0), Point3D(5, 3, 0)
        >>> l1 = Line3D(p1, p2)
        >>> l1.equation()
        (-k + x/4 - 1/4, -k + y/3, zoo*z - k)

        """
        x, y, z, k = _symbol(x), _symbol(y), _symbol(z), _symbol(k)
        p1, p2 = self.points
        a = p1.direction_ratio(p2)
        return (((x - p1.x)/a[0]) - k, ((y - p1.y)/a[1]) - k, 
                ((z - p1.z)/a[2]) - k)

    def contains(self, o):
        """Return True if o is on this Line, or False otherwise."""
        if isinstance(o, Point3D):
            o = o.func(*[simplify(i) for i in o.args])
            eq = self.equation()
            a = []
            for i in range(3):
                a.append(eq[i].subs(eq[i].args[0], o.args[i]))
            if len(set(a)) == 1:
                return True
            else:
                return False
        elif not isinstance(o, Base):
            return False
        elif isinstance(o, Line3D):
            return self.__eq__(o)
        else:
            return o.p1 in self and o.p2 in self and o.p3 in self

    def distance(self, o):
        """
        Finds the shortest distance between a line and a point.

        Raises
        ======

        NotImplementedError is raised if o is not an instance of Point3D

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(1, 1, 1)
        >>> s = Line3D(p1, p2)
        >>> s.distance(Point3D(-1, 1, 1))
        2*sqrt(6)/3
        """
        if not isinstance(o, Point3D):
            raise NotImplementedError
        if o in self:
            return 0
        a = self.perpendicular_segment(o).length
        return a

    def __eq__(self, other):
        """Return True if other is equal to this Line, or False otherwise."""
        if not isinstance(other, Line3D):
            return False
        return Point3D.is_collinear(self.p1, self.p2, other.p1, other.p2)

    def __hash__(self):
        return super(Line3D, self).__hash__()

class Ray3D(Base):
    
    """
    A Ray is a semi-line in the space with a source point and a direction.

    Parameters
    ==========

    p1 : Point3D
        The source of the Ray
    p2 : Point or a direction vector
    direction_ratio: Determines the direction in which the Ray propagates.
        

    Attributes
    ==========

    source
    xdirection
    ydirection
    zdirection

    See Also
    ========

    sympy.geometry.point3d.Point3D, Line3D


    Examples
    ========

    >>> import sympy
    >>> from sympy import Point3D, pi
    >>> from sympy.abc import r
    >>> from sympy.geometry import Ray3D
    >>> r = Ray3D(Point3D(2, 3, 4), Point3D(3, 5, 0))
    >>> r
    Ray(Point3D(2, 3, 4), Point3D(3, 5, 0))
    >>> r.points
    (Point3D(2, 3, 4), Point3D(3, 5, 0))
    >>> r.source
    Point3D(2, 3, 4)
    >>> r.xdirection
    oo
    >>> r.ydirection
    oo
    >>> r.direction_ratio
    [1, 2, -4]

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
        
    @property
    def source(self):
        """The point from which the ray emanates.

        See Also
        ========

        sympy.geometry.point3d.Point3D

        Examples
        ========

        >>> from sympy import Point3D, Ray3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(4, 1, 5)
        >>> r1 = Ray3D(p1, p2)
        >>> r1.source
        Point3D(0, 0, 0)

        """
        return self.p1
        
    @property
    def xdirection(self):
        """The x direction of the ray.

        Positive infinity if the ray points in the positive x direction,
        negative infinity if the ray points in the negative x direction,
        or 0 if the ray is vertical.

        See Also
        ========

        ydirection

        Examples
        ========

        >>> from sympy import Point3D, Ray3D
        >>> p1, p2, p3 = Point3D(0, 0, 0), Point3D(1, 1, 1), Point3D(0, -1, 0)
        >>> r1, r2 = Ray3D(p1, p2), Ray3D(p1, p3)
        >>> r1.xdirection
        oo
        >>> r2.xdirection
        0

        """
        if self.p1.x < self.p2.x:
            return S.Infinity
        elif self.p1.x == self.p2.x:
            return S.Zero
        else:
            return S.NegativeInfinity

    @property
    def ydirection(self):
        """The y direction of the ray.

        Positive infinity if the ray points in the positive y direction,
        negative infinity if the ray points in the negative y direction,
        or 0 if the ray is horizontal.

        See Also
        ========

        xdirection

        Examples
        ========

        >>> from sympy import Point, Ray
        >>> p1, p2, p3 = Point3D(0, 0, 0), Point3D(-1, -1, -1), Point3D(-1, 0, 0)
        >>> r1, r2 = Ray3D(p1, p2), Ray3D(p1, p3)
        >>> r1.ydirection
        -oo
        >>> r2.ydirection
        0

        """
        if self.p1.y < self.p2.y:
            return S.Infinity
        elif self.p1.y == self.p2.y:
            return S.Zero
        else:
            return S.NegativeInfinity
            
    @property
    def zdirection(self):
        """The y direction of the ray.

        Positive infinity if the ray points in the positive y direction,
        negative infinity if the ray points in the negative y direction,
        or 0 if the ray is horizontal.

        See Also
        ========

        xdirection

        Examples
        ========

        >>> from sympy import Point3D, Ray3D
        >>> p1, p2, p3 = Point3D(0, 0, 0), Point3D(-1, -1, -1), Point3D(-1, 0, 0)
        >>> r1, r2 = Ray3D(p1, p2), Ray3D(p1, p3)
        >>> r1.ydirection
        -oo
        >>> r2.ydirection
        0
        >>> r2.zdirection
        0

        """
        if self.p1.z < self.p2.z:
            return S.Infinity
        elif self.p1.z == self.p2.z:
            return S.Zero
        else:
            return S.NegativeInfinity
            
    def plot_interval(self, parameter='t'):
        """The plot interval for the default geometric plot of the Ray. Gives
        values that will produce a ray that is 10 units long (where a unit is
        the distance between the two points that define the ray).

        Parameters
        ==========

        parameter : str, optional
            Default value is 't'.

        Returns
        =======

        plot_interval : list
            [parameter, lower_bound, upper_bound]

        Examples
        ========

        >>> from sympy import Point3D, Ray3D, pi
        >>> r = Ray((0, 0), angle=pi/4)
        >>> r.plot_interval()
        [t, 0, 10]

        """
        t = _symbol(parameter)
        return [t, 0, 10]

    def __eq__(self, other):
        """Is the other GeometryEntity equal to this Ray?"""
        if not isinstance(other, Ray):
            return False
        return (self.source == other.source) and (other.p2 in self)

    def __hash__(self):
        return super(Ray, self).__hash__()        
    
class Segment3D(Base):
    """A undirected line segment in a 3D space.

    Parameters
    ==========

    p1 : Point3D
    p2 : Point3D

    Attributes
    ==========

    length : number or sympy expression
    midpoint : Point3D

    See Also
    ========

    sympy.geometry.point.Point3D, Line3D

    Examples
    ========

    >>> import sympy
    >>> from sympy import Point3D
    >>> from sympy.abc import s
    >>> from sympy.geometry import Segment3D
    >>> Segment3D((1, 0, 0), (1, 1, 1)) # tuples are interpreted as pts
    Segment3D(Point3D(1, 0, 0), Point3D(1, 1, 1))
    >>> s = Segment3D(Point3D(4, 3, 9), Point3D(1, 1, 7))
    >>> s
    Segment3D(Point3D(1, 1, 7), Point3D(4, 3, 9))
    >>> s.points
    (Point3D(1, 1, 7), Point3D(4, 3, 9))
    >>> s.length
    sqrt(17)
    >>> s.midpoint
    Point3D(5/2, 2, 8)

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
        
    def plot_interval(self, parameter='t'):
        """The plot interval for the default geometric plot of the Segment gives
        values that will produce the full segment in a plot.

        Parameters
        ==========

        parameter : str, optional
            Default value is 't'.

        Returns
        =======

        plot_interval : list
            [parameter, lower_bound, upper_bound]

        Examples
        ========

        >>> from sympy import Point3D, Segment3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(5, 3, 0)
        >>> s1 = Segment3D(p1, p2)
        >>> s1.plot_interval()
        [t, 0, 1]

        """
        t = _symbol(parameter)
        return [t, 0, 1]    
        
    @property
    def length(self):
        """The length of the line segment.

        See Also
        ========

        sympy.geometry.point3d.Point3D.distance

        Examples
        ========

        >>> from sympy import Point3D, Segment3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(4, 3, 3)
        >>> s1 = Segment3D(p1, p2)
        >>> s1.length
        sqrt(34)

        """
        return Point3D.distance(self.p1, self.p2)

    @property
    def midpoint(self):
        """The midpoint of the line segment.

        See Also
        ========

        sympy.geometry.point3d.Point3D.midpoint

        Examples
        ========

        >>> from sympy import Point3D, Segment3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(4, 3, 3)
        >>> s1 = Segment3D(p1, p2)
        >>> s1.midpoint
        Point3D(2, 3/2, 3/2)

        """
        return Point3D.midpoint(self.p1, self.p2)
        