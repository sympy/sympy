from sympy.core.basic import Basic, S, C
from sympy.simplify import simplify
from sympy.geometry.exceptions import GeometryError
from entity import GeometryEntity
from point import Point

class LinearEntity(GeometryEntity):
    """
    A linear entity (line, ray, segment, etc) in space.

    This is an abstract class and is not meant to be instantiated.
    Subclasses should implement the following methods:
        __eq__
        __contains__
    """
    def __new__(cls, p1, p2, **kwargs):
        if not isinstance(p1, Point) or not isinstance(p2, Point):
            raise TypeError("%s.__new__ requires Point instances" % cls.__name__)
        if p1 == p2:
            raise RuntimeError("%s.__new__ requires two distinct points" % cls.__name__)

        return GeometryEntity.__new__(cls, p1, p2, **kwargs)

    @property
    def p1(self):
        """One of the defining points of a linear entity."""
        return self.__getitem__(0)

    @property
    def p2(self):
        """One of the defining points of a linear entity."""
        return self.__getitem__(1)

    @property
    def coefficients(self):
        """The coefficients (a,b,c) for equation ax+by+c=0"""
        return (self.p1[1]-self.p2[1],
                self.p2[0]-self.p1[0],
                self.p1[0]*self.p2[1] - self.p1[1]*self.p2[0])

    def is_concurrent(*lines):
        """
        Returns True if the set of linear entities are concurrent, False
        otherwise. Two or more linear entities are concurrent if they all
        intersect at a single point.

        Description of Method Used:
        ===========================
            Simply take the first two lines and find their intersection.
            If there is no intersection, then the first two lines were
            parallel and had no intersection so concurrency is impossible
            amongst the whole set. Otherwise, check to see if the
            intersection point of the first two lines is a member on
            the rest of the lines. If so, the lines are concurrent.
        """
        _lines = lines
        lines = GeometryEntity.extract_entities(lines)

        # Concurrency requires intersection at a single point; One linear
        # entity cannot be concurrent.
        if len(lines) <= 1:
            return False

        try:
            # Get the intersection (if parallel)
            p = GeometryEntity.do_intersection(lines[0], lines[1])
            if len(p) == 0: return False

            # Make sure the intersection is on every linear entity
            for line in lines[2:]:
                if p[0] not in line:
                    return False
            return True
        except AttributeError:
            return False

    def is_parallel(l1, l2):
        """Returns True if l1 and l2 are parallel, False otherwise"""
        try:
            a1,b1,c1 = l1.coefficients
            a2,b2,c2 = l2.coefficients
            return bool(simplify(a1*b2 - b1*a2) == 0)
        except AttributeError:
            return False

    def is_perpendicular(l1, l2):
        """Returns True if l1 and l2 are perpendicular, False otherwise"""
        try:
            a1,b1,c1 = l1.coefficients
            a2,b2,c2 = l2.coefficients
            return bool(simplify(a1*a2 + b1*b2) == 0)
        except AttributeError:
            return False

    def angle_between(l1, l2):
        """
        Returns an angle formed between the two linear entities.

        Description of Method Used:
        ===========================
            From the dot product of vectors v1 and v2 it is known that:
                dot(v1, v2) = |v1|*|v2|*cos(A)
            where A is the angle formed between the two vectors. We can
            get the directional vectors of the two lines and readily
            find the angle between the two using the above formula.
        """
        v1 = l1.p2 - l1.p1
        v2 = l2.p2 - l2.p1
        return C.acos( (v1[0]*v2[0]+v1[1]*v2[1]) / (abs(v1)*abs(v2)) )

    def parallel_line(self, p):
        """
        Returns a new Line which is parallel to this linear entity and passes
        through the specified point.
        """
        d = self.p1 - self.p2
        return Line(p, p + d)

    def perpendicular_line(self, p):
        """
        Returns a new Line which is perpendicular to this linear entity and
        passes through the specified point.
        """
        d1,d2 = self.p1 - self.p2
        if d2 == 0: # If an horizontal line
            if p[1] == self.p1[1]: # if p is on this linear entity
                p2 = Point(p[0], p[1] + 1)
                return Line(p, p2)
            else:
                p2 = Point(p[0], self.p1[1])
                return Line(p, p2)
        else:
            p2 = Point(p[0] - d2, p[1] + d1)
            return Line(p, p2)

    def perpendicular_segment(self, p):
        """
        Returns a new Segment which connects p to a point on this linear
        entity and is also perpendicular to this line. Returns p itself
        if p is on this linear entity.
        """
        if p in self:
            return p
        pl = self.perpendicular_line(p)
        p2 = GeometryEntity.do_intersection(self, pl)[0]
        return Segment(p, p2)

    @property
    def slope(self):
        """
        The slope of this linear entity, or infinity if vertical.
        """
        d1,d2 = self.p1 - self.p2
        if d1 == 0:
            return S.Infinity
        return simplify(d2/d1)

    @property
    def points(self):
        """The two points used to define this linear entity."""
        return (self.p1, self.p2)

    def projection(self, o):
        """
        Project a point, line, ray, or segment onto this linear entity.
        If projection cannot be performed then a GeometryError is raised.

        Notes:
        ======
            - A projection involves taking the two points that define
              the linear entity and projecting those points onto a
              Line and then reforming the linear entity using these
              projections.
            - A point P is projected onto a line L by finding the point
              on L that is closest to P. This is done by creating a
              perpendicular line through P and L and finding its
              intersection with L.
        """
        tline = Line(self.p1, self.p2)

        def project(p):
            """Project a point onto the line representing self."""
            if p in tline: return p
            l1 = tline.perpendicular_line(p)
            return tline.intersection(l1)[0]

        projected = None
        if isinstance(o, Point):
            return project(o)
        elif isinstance(o, LinearEntity):
            n_p1 = project(o.p1)
            n_p2 = project(o.p2)
            if n_p1 == n_p2:
                projected = n_p1
            else:
                projected = o.__class__(n_p1, n_p2)

        # Didn't know how to project so raise an error
        if projected is None:
            n1 = self.__class__.__name__
            n2 = o.__class__.__name__
            raise GeometryError("Do not know how to project %s onto %s" % (n2, n1))

        return GeometryEntity.do_intersection(self, projected)[0]

    def intersection(self, o):
        if isinstance(o, Point):
            if o in self:
                return [o]
            else:
                return []
        elif isinstance(o, LinearEntity):
            a1,b1,c1 = self.coefficients
            a2,b2,c2 = o.coefficients
            t = simplify(a1*b2 - a2*b1)
            if t == 0: # are parallel?
                if isinstance(self, Line):
                    if o.p1 in self:
                        return [o]
                    return []
                elif isinstance(o, Line):
                    if self.p1 in o:
                        return [self]
                    return []
                elif isinstance(self, Ray):
                    if isinstance(o, Ray):
                        # case 1, rays in the same direction
                        if self.xdirection == o.xdirection:
                            if self.source[0] < o.source[0]:
                                return [o]
                            return [self]
                        # case 2, rays in the opposite directions
                        else:
                            if o.source in self:
                                if self.source == o.source:
                                    return [self.source]
                                return [Segment(o.source, self.source)]
                            return []
                    elif isinstance(o, Segment):
                        if o.p1 in self:
                            if o.p2 in self:
                                return [o]
                            return [Segment(o.p1, self.source)]
                        elif o.p2 in self:
                            return [Segment(o.p2, self.source)]
                        return []
                elif isinstance(self, Segment):
                    if isinstance(o, Ray):
                        return o.intersection(self)
                    elif isinstance(o, Segment):
                        # A reminder that the points of Segments are ordered
                        # in such a way that the following works. See
                        # Segment.__new__ for details on the ordering.
                        if self.p1 not in o:
                            if self.p2 not in o:
                                # Neither of the endpoints are in o so either
                                # o is contained in this segment or it isn't
                                if o in self:
                                    return [self]
                                return []
                            else:
                                # p1 not in o but p2 is. Either there is a
                                # segment as an intersection, or they only
                                # intersect at an endpoint
                                if self.p2 == o.p1:
                                    return [o.p1]
                                return [Segment(o.p1, self.p2)]
                        elif self.p2 not in o:
                            # p2 not in o but p1 is. Either there is a
                            # segment as an intersection, or they only
                            # intersect at an endpoint
                            if self.p1 == o.p2:
                                return [o.p2]
                            return [Segment(o.p2, self.p1)]

                        # Both points of self in o so the whole segment
                        # is in o
                        return [self]

                # Unknown linear entity
                return []

            # Not parallel, so find the point of intersection
            px = simplify((b1*c2 - c1*b2) / t)
            py = simplify((a2*c1 - a1*c2) / t)
            inter = Point(px, py)
            if (inter in self) and (inter in o):
                return [inter]
            return []
        raise NotImplementedError()

    def random_point(self):
        """Returns a random point on this Ray."""
        from random import randint
        from sys import maxint

        # The lower and upper
        lower, upper = -maxint-1, maxint

        if self.slope is S.Infinity:
            if isinstance(self, Ray):
                if self.ydirection is S.Infinity:
                    lower = self.p1[1]
                else:
                    upper = self.p1[1]
            elif isinstance(self, Segment):
                lower = self.p1[1]
                upper = self.p2[1]

            x = self.p1[0]
            y = randint(lower, upper)
        else:
            if isinstance(self, Ray):
                if self.xdirection is S.Infinity:
                    lower = self.p1[0]
                else:
                    upper = self.p1[0]
            elif isinstance(self, Segment):
                lower = self.p1[0]
                upper = self.p2[0]

            a,b,c = self.coefficients
            x = randint(lower, upper)
            y = simplify( (-c - a*x) / b )
        return Point(x, y)

    def __eq__(self, other):
        raise NotImplementedError()

    def __contains__(self, other):
        raise NotImplementedError()


class Line(LinearEntity):
    """A line in space.

    A line is declared with two distinct points.

    Note:
    At the moment only lines in a 2D space can be declared, because
    Points can be defined only for 2D spaces.

    *Example*

        >>> import sympy
        >>> from sympy import Point
        >>> from sympy.abc import L
        >>> from sympy.geometry import Line
        >>> L = Line(Point(2,3), Point(3,5))
        >>> L
        Line(Point(2, 3), Point(3, 5))
        >>> L.points
        (Point(2, 3), Point(3, 5))
        >>> L.equation()
        1 + y - 2*x
        >>> L.coefficients
        (-2, 1, 1)
    """

    def arbitrary_point(self, parameter_name='t'):
        """Returns a symbolic point that is on this line."""
        t = C.Symbol(parameter_name, real=True)
        x = simplify(self.p1[0] + t*(self.p2[0] - self.p1[0]))
        y = simplify(self.p1[1] + t*(self.p2[1] - self.p1[1]))
        return Point(x, y)

    def plot_interval(self, parameter_name='t'):
        """Returns the plot interval for the default geometric plot of line"""
        t = C.Symbol(parameter_name, real=True)
        return [t, -5, 5]

    def equation(self, xaxis_name='x', yaxis_name='y'):
        """
        Returns the equation for this line. Optional parameters xaxis_name
        and yaxis_name can be used to specify the names of the symbols used
        for the equation.
        """
        x = C.Symbol(xaxis_name, real=True)
        y = C.Symbol(yaxis_name, real=True)
        a,b,c = self.coefficients
        return simplify(a*x + b*y + c)

    def __contains__(self, o):
        """Return True if o is on this Line, or False otherwise."""
        if isinstance(o, Line):
            return self.__eq__(o)
        elif isinstance(o, Point):
            x = C.Symbol('x', real=True)
            y = C.Symbol('y', real=True)
            r = self.equation().subs({x: o[0], y: o[1]})
            x = simplify(r)
            return simplify(x) == 0
        else:
            return False

    def __eq__(self, other):
        """Return True if other is equal to this Line, or False otherwise."""
        if not isinstance(other, Line): return False
        return Point.is_collinear(self.p1, self.p2, other.p1, other.p2)


class Ray(LinearEntity):
    """
    A ray is a semi-line in the space. It starts at one point and
    propagates in one unique direction.

    A ray is declared with two distinct points: the first point is the source,
    whereas the second point lies on the semi-line. Therefore, the second
    point determines the direction to which the semi-line propagates.

    Note:
    At the moment only rays in a 2D space can be declared, because
    Points can be defined only for 2D spaces.

    *Example*

        >>> import sympy
        >>> from sympy import Point
        >>> from sympy.abc import r
        >>> from sympy.geometry import Ray
        >>> r = Ray(Point(2, 3), Point(3, 5))
        >>> r = Ray(Point(2, 3), Point(3, 5))
        >>> r
        Ray(Point(2, 3), Point(3, 5))
        >>> r.points
        (Point(2, 3), Point(3, 5))
        >>> r.source
        Point(2, 3)
        >>> r.xdirection
        oo
        >>> r.ydirection
        oo
        >>> r.slope
        2
    """

    @property
    def source(self):
        """The point from which the ray emanates."""
        return self.p1

    @property
    def xdirection(self):
        """
        The x direction of the ray. Positive infinity if the ray points in
        the positive x direction, negative infinity if the ray points
        in the negative x direction, or 0 if the ray is vertical.
        """
        if self.p1[0] < self.p2[0]:
            return S.Infinity
        elif self.p1[0] == self.p2[0]:
            return S.Zero
        else:
            return S.NegativeInfinity

    @property
    def ydirection(self):
        """
        The y direction of the ray. Positive infinity if the ray points in
        the positive y direction, negative infinity if the ray points
        in the negative y direction, or 0 if the ray is horizontal.
        """
        if self.p1[1] < self.p2[1]:
            return S.Infinity
        elif self.p1[1] == self.p2[1]:
            return S.Zero
        else:
            return S.NegativeInfinity

    def __eq__(self, other):
        """Return True if other is equal to this Ray, or False otherwise."""
        if not isinstance(other, Ray):
            return False
        return ((self.source == other.source) and (other.p2 in self))

    def __contains__(self, o):
        """Return True if o is on this Ray, or False otherwise."""
        if isinstance(o, Ray):
            d = o.p2 - o.p1
            return Point.is_collinear(self.p1, self.p2, o.p1, o.p2) \
                    and (self.xdirection == o.xdirection) \
                    and (self.ydirection == o.ydirection)
        elif isinstance(o, Segment):
            return ((o.p1 in self) and (o.p2 in self))
        elif isinstance(o, Point):
            if Point.is_collinear(self.p1, self.p2, o):
                if (not self.p1[0].atoms(C.Symbol)) and (not self.p1[1].atoms(C.Symbol)) \
                        and (not self.p2[0].atoms(C.Symbol)) and (not self.p2[1].atoms(C.Symbol)):
                    if self.xdirection is S.Infinity:
                        return o[0] >= self.source[0]
                    elif self.xdirection is S.NegativeInfinity:
                        return o[0] <= self.source[0]
                    elif self.ydirection is S.Infinity:
                        return o[1] >= self.source[1]
                    return o[1] <= self.source[1]
                else:
                    # There are symbols lying around, so assume that o
                    # is contained in this ray (for now)
                    return True
            else:
                # Points are not collinear, so the rays are not parallel
                # and hence it isimpossible for self to contain o
                return False

        # No other known entity can be contained in a Ray
        return False


class Segment(LinearEntity):
    """An undirected line segment in space.

    A segment is declared with two distinct points.

    Note:
    At the moment only segments in a 2D space can be declared, because
    Points can be defined only for 2D spaces.

    *Example*

        >>> import sympy
        >>> from sympy import Point
        >>> from sympy.abc import s
        >>> from sympy.geometry import Segment
        >>> s = Segment(Point(4, 3), Point(1, 1))
        >>> s
        Segment(Point(1, 1), Point(4, 3))
        >>> s.points
        (Point(1, 1), Point(4, 3))
        >>> s.slope
        2/3
        >>> s.length
        13**(1/2)
        >>> s.midpoint
        Point(5/2, 2)
    """

    def __new__(cls, p1, p2, **kwargs):
        # Reorder the two points under the following ordering:
        #   if p1[0] != p2[0] then p1[0] < p2[0]
        #   if p1[0] == p2[0] then p1[1] < p2[1]
        if p1[0] > p2[0]:
            p1, p2 = p2, p1
        elif p1[0] == p2[0] and p1[1] > p2[0]:
            p1, p2 = p2, p1
        return LinearEntity.__new__(cls, p1, p2, **kwargs)

    def arbitrary_point(self, parameter_name='t'):
        """Returns a symbolic point that is on this line segment."""
        t = C.Symbol(parameter_name, real=True)
        x = simplify(self.p1[0] + t*(self.p2[0] - self.p1[0]))
        y = simplify(self.p1[1] + t*(self.p2[1] - self.p1[1]))
        return Point(x, y)

    def plot_interval(self, parameter_name='t'):
        t = C.Symbol(parameter_name, real=True)
        return [t, 0, 1]

    def perpendicular_bisector(self, p=None):
        """
        Returns the perpendicular bisector of this segment. If no point is
        specified or the point specified is not on the bisector then the
        bisector is returned as a Line. Otherwise a Segment is returned
        that joins the point specified and the intersection of the bisector
        and the segment.
        """
        l = LinearEntity.perpendicular_line(self, self.midpoint)
        if p is None or p not in l:
            return l
        else:
            return Segment(self.midpoint, p)

    @property
    def length(self):
        """The length of the segment."""
        return Point.distance(self.p1, self.p2)

    @property
    def midpoint(self):
        """The midpoint of the segment."""
        return Point.midpoint(self.p1, self.p2)

    def __eq__(self, other):
        """Return True if other is equal to this Line, or False otherwise."""
        if not isinstance(other, Segment):
            return False
        return ((self.p1 == other.p1) and (self.p2 == other.p2))

    def __contains__(self, o):
        """Return True if o is on this Segment, or False otherwise."""
        if isinstance(o, Segment):
            return ((o.p1 in self) and (o.p2 in self))
        elif isinstance(o, Point):
            if Point.is_collinear(self.p1, self.p2, o):
                d = self.length
                if (Point.distance(self.p1,o) <= d) and (Point.distance(self.p2,o) <= d):
                    return True

        # No other known entity can be contained in a Ray
        return False
