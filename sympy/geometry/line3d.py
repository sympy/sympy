"""Line-like geometrical entities.

Contains
========
LinearEntity3D
Line3D
Ray3D
Segment3D

"""
from __future__ import division, print_function

from sympy.core import Dummy, S, nan
from sympy.functions.elementary.trigonometric import acos
from sympy.simplify.simplify import simplify
from sympy.solvers import solve
from sympy.solvers.solveset import solveset, linsolve
from sympy.geometry.exceptions import GeometryError
from sympy.core.compatibility import is_sequence, range

from .entity import GeometryEntity
from .point import Point3D
from .util import _symbol


class LinearEntity3D(GeometryEntity):
    """An base class for all linear entities (line, ray and segment)
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

    This is a base class and is not meant to be instantiated.
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

        sympy.geometry.point.Point3D

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(5, 3, 1)
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

        sympy.geometry.point.Point3D

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
        """The normalized direction ratio of a given line in 3D.

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
        >>> sum(i**2 for i in _)
        1
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

        sympy.geometry.point.Point3D

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(5, 11, 1)
        >>> l1 = Line3D(p1, p2)
        >>> l1.points
        (Point3D(0, 0, 0), Point3D(5, 11, 1))

        """
        return (self.p1, self.p2)

    @staticmethod
    def are_concurrent(*lines):
        """Is a sequence of linear entities concurrent?

        Two or more linear entities are concurrent if they all
        intersect at a single point.

        Parameters
        ==========

        lines : a sequence of linear entities.

        Returns
        =======

        True : if the set of linear entities are concurrent,
        False : otherwise.

        Notes
        =====

        Simply take the first two lines and find their intersection.
        If there is no intersection, then the first two lines were
        parallel and had no intersection so concurrency is impossible
        amongst the whole set. Otherwise, check to see if the
        intersection point of the first two lines is a member on
        the rest of the lines. If so, the lines are concurrent.

        See Also
        ========

        sympy.geometry.util.intersection

        Examples
        ========

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(3, 5, 2)
        >>> p3, p4 = Point3D(-2, -2, -2), Point3D(0, 2, 1)
        >>> l1, l2, l3 = Line3D(p1, p2), Line3D(p1, p3), Line3D(p1, p4)
        >>> Line3D.are_concurrent(l1, l2, l3)
        True

        >>> l4 = Line3D(p2, p3)
        >>> Line3D.are_concurrent(l2, l3, l4)
        False

        """

        # Concurrency requires intersection at a single point; One linear
        # entity cannot be concurrent.
        if len(lines) <= 1:
            return False

        try:
            # Get the intersection (if parallel)
            p = lines[0].intersection(lines[1])
            if len(p) == 0:
                return False

            # Make sure the intersection is on every linear entity
            for line in lines[2:]:
                if p[0] not in line:
                    return False
            return True
        except AttributeError:
            return False

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
            return True
        a = l1.direction_cosine
        b = l2.direction_cosine
        # lines are parallel if the direction_cosines are the same or
        # differ by a constant
        rat = set()
        for i, j in zip(a, b):
            if i and j:
                rat.add(i/j)
                if len(rat) > 1:
                    return False
            elif i or j:
                return False
        return True

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
        >>> l1, l2 = Line3D(p1, p2), Line3D(p2, p3)
        >>> l1.is_perpendicular(l2)
        False

        >>> p4 = Point3D(5, 3, 7)
        >>> l3 = Line3D(p1, p4)
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
        >>> l1, l2 = Line3D(p1, p2), Line3D(p2, p3)
        >>> l1.angle_between(l2)
        acos(-sqrt(2)/3)

        """
        v1 = l1.p2 - l1.p1
        v2 = l2.p2 - l2.p1
        return acos(v1.dot(v2)/(abs(v1)*abs(v2)))

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
        p = Point3D(p)
        if p in self:
            raise NotImplementedError("Given point should not be on the line")
        t = Dummy()
        a = self.arbitrary_point(t)
        b = [i - j for i, j in zip(p.args, a.args)]
        c = sum([i*j for i, j in zip(b, self.direction_ratio)])
        d = list(solveset(c, t))
        e = a.subs(t, d[0])
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
        p = Point3D(p)
        if p in self:
            raise NotImplementedError("Given point should not be on the line")
        t = Dummy()
        a = self.arbitrary_point(t)
        b = [i - j for i, j in zip(p.args, a.args)]
        c = sum([i*j for i, j in zip(b, self.direction_ratio)])
        d = list(solveset(c, t))
        e = a.subs(t, d[0])
        return Segment3D(p, e)

    def projection(self, o):
        """Project a point, line, ray, or segment onto this linear entity.

        Parameters
        ==========

        other : Point or LinearEntity (Line, Ray, Segment)

        Returns
        =======

        projection : Point or LinearEntity (Line, Ray, Segment)
            The return type matches the type of the parameter ``other``.

        Raises
        ======

        GeometryError
            When method is unable to perform projection.

        Notes
        =====

        A projection involves taking the two points that define
        the linear entity and projecting those points onto a
        Line and then reforming the linear entity using these
        projections.
        A point P is projected onto a line L by finding the point
        on L that is closest to P. This point is the intersection
        of L and the line perpendicular to L that passes through P.

        See Also
        ========

        sympy.geometry.point.Point3D, perpendicular_line

        Examples
        ========

        >>> from sympy import Point3D, Line3D, Segment3D, Rational
        >>> p1, p2, p3 = Point3D(0, 0, 1), Point3D(1, 1, 2), Point3D(2, 0, 1)
        >>> l1 = Line3D(p1, p2)
        >>> l1.projection(p3)
        Point3D(2/3, 2/3, 5/3)

        >>> p4, p5 = Point3D(10, 0, 1), Point3D(12, 1, 3)
        >>> s1 = Segment3D(p4, p5)
        >>> l1.projection(s1)
        [Segment3D(Point3D(10/3, 10/3, 13/3), Point3D(5, 5, 6))]

        """
        tline = Line3D(self.p1, self.p2)

        def _project(p):
            """Project a point onto the line representing self."""
            if p in tline:
                return p
            l1 = tline.perpendicular_line(p)
            return tline.intersection(l1)[0]

        projected = None
        if isinstance(o, Point3D):
            return _project(o)
        elif isinstance(o, LinearEntity3D):
            n_p1 = _project(o.p1)
            n_p2 = _project(o.p2)
            if n_p1 == n_p2:
                projected = n_p1
            else:
                projected = o.__class__(n_p1, n_p2)

        # Didn't know how to project so raise an error
        if projected is None:
            n1 = self.__class__.__name__
            n2 = o.__class__.__name__
            raise GeometryError(
                "Do not know how to project %s onto %s" % (n2, n1))

        return self.intersection(projected)

    def intersection(self, o):
        """The intersection with another geometrical entity.

        Parameters
        ==========

        o : Point or LinearEntity3D

        Returns
        =======

        intersection : list of geometrical entities

        See Also
        ========

        sympy.geometry.point.Point3D

        Examples
        ========

        >>> from sympy import Point3D, Line3D, Segment3D
        >>> p1, p2, p3 = Point3D(0, 0, 0), Point3D(1, 1, 1), Point3D(7, 7, 7)
        >>> l1 = Line3D(p1, p2)
        >>> l1.intersection(p3)
        [Point3D(7, 7, 7)]

        >>> l1 = Line3D(Point3D(4,19,12), Point3D(5,25,17))
        >>> l2 = Line3D(Point3D(-3, -15, -19), direction_ratio=[2,8,8])
        >>> l1.intersection(l2)
        [Point3D(1, 1, -3)]

        >>> p6, p7 = Point3D(0, 5, 2), Point3D(2, 6, 3)
        >>> s1 = Segment3D(p6, p7)
        >>> l1.intersection(s1)
        []

        """
        if isinstance(o, Point3D):
            if o in self:
                return [o]
            else:
                return []

        elif isinstance(o, LinearEntity3D):
            if self == o:
                return [self]
            elif self.is_parallel(o):
                if isinstance(self, Line3D):
                    if o.p1 in self:
                        return [o]
                    return []
                elif isinstance(self, Ray3D):
                    if isinstance(o, Ray3D):
                        # case 1, rays in the same direction
                        if self.xdirection == o.xdirection and \
                                self.ydirection == o.ydirection and \
                                self.zdirection == o.zdirection:
                            return [self] if (self.source in o) else [o]
                        # case 2, rays in the opposite directions
                        else:
                            if o.source in self:
                                if self.source == o.source:
                                    return [self.source]
                                return [Segment3D(o.source, self.source)]
                            return []
                    elif isinstance(o, Segment3D):
                        if o.p1 in self:
                            if o.p2 in self:
                                return [o]
                            return [Segment3D(o.p1, self.source)]
                        elif o.p2 in self:
                            return [Segment3D(o.p2, self.source)]
                        return []
                elif isinstance(self, Segment3D):
                    if isinstance(o, Segment3D):
                        # A reminder that the points of Segments are ordered
                        # in such a way that the following works. See
                        # Segment3D.__new__ for details on the ordering.
                        if self.p1 not in o:
                            if self.p2 not in o:
                                # Neither of the endpoints are in o so either
                                # o is contained in this segment or it isn't
                                if o in self:
                                    return [o]
                                return []
                            else:
                                # p1 not in o but p2 is. Either there is a
                                # segment as an intersection, or they only
                                # intersect at an endpoint
                                if self.p2 == o.p1:
                                    return [o.p1]
                                return [Segment3D(o.p1, self.p2)]
                        elif self.p2 not in o:
                            # p2 not in o but p1 is. Either there is a
                            # segment as an intersection, or they only
                            # intersect at an endpoint
                            if self.p1 == o.p2:
                                return [o.p2]
                            return [Segment3D(o.p2, self.p1)]

                        # Both points of self in o so the whole segment
                        # is in o
                        return [self]

                else:  # unrecognized LinearEntity
                    raise NotImplementedError

            else:
                # If the lines are not parallel then solve their arbitrary points
                # to obtain the point of intersection
                t = t1, t2 = Dummy(), Dummy()
                a = self.arbitrary_point(t1)
                b = o.arbitrary_point(t2)
                dx = a.x - b.x
                c = linsolve([dx, a.y - b.y], t).args[0]
                d = linsolve([dx, a.z - b.z], t).args[0]
                if len(c.free_symbols) == 1 and len(d.free_symbols) == 1:
                    return []
                e = a.subs(t1, c[0])
                if e in self and e in o:
                    return [e]
                else:
                    return []

        return o.intersection(self)

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

        sympy.geometry.point.Point3D

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

class Line3D(LinearEntity3D):
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

    sympy.geometry.point.Point3D

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
        if isinstance(p1, LinearEntity3D):
            p1, pt = p1.args
        else:
            p1 = Point3D(p1)
        if pt is not None and len(direction_ratio) == 0:
            pt = Point3D(pt)
        elif len(direction_ratio) == 3 and pt is None:
            pt = Point3D(p1.x + direction_ratio[0], p1.y + direction_ratio[1],
                         p1.z + direction_ratio[2])
        else:
            raise ValueError('A 2nd Point or keyword "direction_ratio" must '
            'be used.')

        return LinearEntity3D.__new__(cls, p1, pt, **kwargs)

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

        >>> from sympy import Point3D, Line3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(5, 3, 1)
        >>> l1 = Line3D(p1, p2)
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
        (x/4 - 1/4, y/3, zoo*z, k)

        """
        x, y, z, k = _symbol(x), _symbol(y), _symbol(z), _symbol(k)
        p1, p2 = self.points
        a = p1.direction_ratio(p2)
        return (((x - p1.x)/a[0]), ((y - p1.y)/a[1]),
                ((z - p1.z)/a[2]), k)

    def contains(self, o):
        """Return True if o is on this Line, or False otherwise.

        Examples
        ========

        >>> from sympy import Line3D
        >>> a = (0, 0, 0)
        >>> b = (1, 1, 1)
        >>> c = (2, 2, 2)
        >>> l1 = Line3D(a, b)
        >>> l2 = Line3D(b, a)
        >>> l1 == l2
        False
        >>> l1 in l2
        True
        """
        if is_sequence(o):
            o = Point3D(o)
        if isinstance(o, Point3D):
            sym = list(map(Dummy, 'xyz'))
            eq = self.equation(*sym)
            a = [eq[i].subs(sym[i], o.args[i]) for i in range(3)]
            a = [i for i in a if i != nan]
            if len(a) == 1:
                return True
            first = a.pop(0)
            for i in a:
                rv = first.equals(i)
                if not rv:
                    return rv
            return True
        elif not isinstance(o, LinearEntity3D):
            return False
        elif isinstance(o, Line3D):
            return all(i in self for i in o.points)

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
        >>> s.distance((-1, 1, 1))
        2*sqrt(6)/3
        """
        if not isinstance(o, Point3D):
            if is_sequence(o):
                o = Point3D(o)
        if o in self:
            return S.Zero
        a = self.perpendicular_segment(o).length
        return a

    def equals(self, other):
        """Returns True if self and other are the same mathematical entities"""
        if not isinstance(other, Line3D):
            return False
        return Point3D.are_collinear(self.p1, other.p1, self.p2, other.p2)

class Ray3D(LinearEntity3D):
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

    sympy.geometry.point.Point3D, Line3D


    Examples
    ========

    >>> import sympy
    >>> from sympy import Point3D, pi
    >>> from sympy.abc import r
    >>> from sympy.geometry import Ray3D
    >>> r = Ray3D(Point3D(2, 3, 4), Point3D(3, 5, 0))
    >>> r
    Ray3D(Point3D(2, 3, 4), Point3D(3, 5, 0))
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
        if isinstance(p1, LinearEntity3D):
            p1, pt = p1.args
        else:
            p1 = Point3D(p1)
        if pt is not None and len(direction_ratio) == 0:
            pt = Point3D(pt)
        elif len(direction_ratio) == 3 and pt is None:
            pt = Point3D(p1.x + direction_ratio[0], p1.y + direction_ratio[1],
                         p1.z + direction_ratio[2])
        else:
            raise ValueError('A 2nd Point or keyword "direction_ratio" must'
            'be used.')

        return LinearEntity3D.__new__(cls, p1, pt, **kwargs)

    @property
    def source(self):
        """The point from which the ray emanates.

        See Also
        ========

        sympy.geometry.point.Point3D

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

        >>> from sympy import Point3D, Ray3D
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
        """The z direction of the ray.

        Positive infinity if the ray points in the positive z direction,
        negative infinity if the ray points in the negative z direction,
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

    def distance(self, o):
        """
        Finds the shortest distance between the ray and a point.

        Raises
        ======

        NotImplementedError is raised if o is not a Point

        Examples
        ========

        >>> from sympy import Point3D, Ray3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(1, 1, 2)
        >>> s = Ray3D(p1, p2)
        >>> s.distance(Point3D(-1, -1, 2))
        sqrt(6)
        >>> s.distance((-1, -1, 2))
        sqrt(6)
        """
        if not isinstance(o, Point3D):
            if is_sequence(o):
                o = Point3D(o)
        if o in self:
            return S.Zero
        s = self.perpendicular_segment(o)
        if not isinstance(s, Point3D):
            non_o = s.p1 if s.p1 == o else s.p2
            if self.contains(non_o):
                return Line3D(self).distance(o)  # = s.length but simpler
        # the following applies when neither of the above apply
        return self.source.distance(o)

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
        >>> r = Ray3D(Point3D(0, 0, 0), Point3D(1, 1, 1))
        >>> r.plot_interval()
        [t, 0, 10]

        """
        t = _symbol(parameter)
        return [t, 0, 10]

    def contains(self, o):
        """Is other GeometryEntity contained in this Ray?"""
        if isinstance(o, Ray3D):
            return (Point3D.are_collinear(self.p1, self.p2, o.p1, o.p2) and
                    self.xdirection == o.xdirection and
                    self.ydirection == o.ydirection and
                    self.zdirection == o.zdirection)
        elif isinstance(o, Segment3D):
            return o.p1 in self and o.p2 in self
        elif is_sequence(o):
            o = Point3D(o)
        if isinstance(o, Point3D):
            if Point3D.are_collinear(self.p1, self.p2, o):
                if self.xdirection is S.Infinity:
                    rv = o.x >= self.source.x
                elif self.xdirection is S.NegativeInfinity:
                    rv = o.x <= self.source.x
                elif self.ydirection is S.Infinity:
                    rv = o.y >= self.source.y
                elif self.ydirection is S.NegativeInfinity:
                    rv = o.y <= self.source.y
                elif self.zdirection is S.Infinity:
                    rv = o.z <= self.source.z
                else:
                    rv = o.z <= self.source.z
                if rv == True or rv == False:
                    return bool(rv)
                raise Undecidable(
                    'Cannot determine if %s is in %s' % (o, self))
            else:
                # Points are not collinear, so the rays are not parallel
                # and hence it is impossible for self to contain o
                return False

        # No other known entity can be contained in a Ray
        return False

    def equals(self, other):
        """Returns True if self and other are the same mathematical entities"""
        if not isinstance(other, Ray3D):
            return False
        return self.source == other.source and other.p2 in self

class Segment3D(LinearEntity3D):
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
        return LinearEntity3D.__new__(cls, p1, p2, **kwargs)

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

        sympy.geometry.point.Point3D.distance

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

        sympy.geometry.point.Point3D.midpoint

        Examples
        ========

        >>> from sympy import Point3D, Segment3D
        >>> p1, p2 = Point3D(0, 0, 0), Point3D(4, 3, 3)
        >>> s1 = Segment3D(p1, p2)
        >>> s1.midpoint
        Point3D(2, 3/2, 3/2)

        """
        return Point3D.midpoint(self.p1, self.p2)

    def distance(self, o):
        """
        Finds the shortest distance between a line segment and a point.

        Raises
        ======

        NotImplementedError is raised if o is not a Point3D

        Examples
        ========

        >>> from sympy import Point3D, Segment3D
        >>> p1, p2 = Point3D(0, 0, 3), Point3D(1, 1, 4)
        >>> s = Segment3D(p1, p2)
        >>> s.distance(Point3D(10, 15, 12))
        sqrt(341)
        >>> s.distance((10, 15, 12))
        sqrt(341)
        """
        if is_sequence(o):
            o = Point3D(o)
        if isinstance(o, Point3D):
            seg_vector = self.p2 - self.p1
            pt_vector = o - self.p1
            t = seg_vector.dot(pt_vector)/self.length**2
            if t >= 1:
                distance = Point3D.distance(self.p2, o)
            elif t <= 0:
                distance = Point3D.distance(self.p1, o)
            else:
                distance = Point3D.distance(
                    self.p1 + Point3D(t*seg_vector.x, t*seg_vector.y,
                                      t*seg_vector.y), o)
            return distance
        raise NotImplementedError()

    def contains(self, other):
        """
        Is the other GeometryEntity contained within this Segment?

        Examples
        ========

        >>> from sympy import Point3D, Segment3D
        >>> p1, p2 = Point3D(0, 1, 1), Point3D(3, 4, 5)
        >>> s = Segment3D(p1, p2)
        >>> s2 = Segment3D(p2, p1)
        >>> s.contains(s2)
        True
        """
        if is_sequence(other):
            other = Point3D(other)
        if isinstance(other, Segment3D):
            return other.p1 in self and other.p2 in self
        elif isinstance(other, Point3D):
            if Point3D.are_collinear(self.p1, self.p2, other):
                if other.distance(self.p1) + other.distance(self.p2) == self.length:
                    return True
                else:
                    return False
        return False
