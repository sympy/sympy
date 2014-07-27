"""Geometrical Planes.

Contains
========
Plane

"""
from __future__ import print_function, division

from sympy.core import S, C, sympify, Dummy, nan, Eq, symbols, Symbol, Rational
from sympy.core.function import expand_mul
from sympy.functions.elementary.trigonometric import _pi_coeff as pi_coeff, \
    sqrt
from sympy.core.logic import fuzzy_and
from sympy.core.exprtools import factor_terms
from sympy.simplify.simplify import simplify
from sympy.solvers import solve
from sympy.polys.polytools import cancel
from sympy.geometry.exceptions import GeometryError
from .entity import GeometryEntity
from .point3d import Point3D
from .point import Point
from .line3d import LinearEntity3D, Line3D, Segment3D, Ray3D
from .line import Line, Segment, Ray
from sympy.matrices import Matrix

class Plane(GeometryEntity):
    """
    A plane is a flat, two-dimensional surface. A plane is the two-dimensional
    analogue of a point (zero-dimensions), a line (one-dimension) and a solid
    (three-dimensions). A plane can generally be constructed by two types of
    inputs. They are three non-collinear points and a point and the plane's
    normal vector.

    Attributes
    ==========

    p1
    normal_vector

    Examples
    ========

    >>> from sympy import Plane, Point3D
    >>> from sympy.abc import x
    >>> Plane(Point3D(1, 1, 1), Point3D(2, 3, 4), Point3D(2, 2, 2))
    Plane(Point3D(1, 1, 1), (-1, 2, -1))
    >>> Plane((1, 1, 1), (2, 3, 4), (2, 2, 2))
    Plane(Point3D(1, 1, 1), (-1, 2, -1))
    >>> Plane(Point3D(1, 1, 1), normal_vector=(1,4,7))
    Plane(Point3D(1, 1, 1), (1, 4, 7))

    """
    def __new__(cls, p1, a=None, b=None, **kwargs):
        p1 = Point3D(p1)
        if not a and not b and kwargs.get('normal_vector', None):
            a = kwargs.pop('normal_vector')

        if not b and not isinstance(a, Point3D) and \
                len(a) == 3:
            normal_vector = a
        elif a and b:
            p2 = Point3D(a)
            p3 = Point3D(b)
            if Point3D.are_collinear(p1, p2, p3):
                raise NotImplementedError('Enter three non-collinear points')
            a = p1.direction_ratio(p2)
            b = p1.direction_ratio(p3)
            normal_vector = tuple(Matrix(a).cross(Matrix(b)))
        else:
            raise ValueError('Either provide 3 3D points or a point with a '
            'normal vector')
        return GeometryEntity.__new__(cls, p1, normal_vector, **kwargs)

    @property
    def p1(self):
        """The only defining point of the plane. Others can be obtained from the
        arbitrary_point method.

        See Also
        ========

        sympy.geometry.point3d.Point3D

        Examples
        ========

        >>> from sympy import Point3D, Plane
        >>> a = Plane(Point3D(1, 1, 1), Point3D(2, 3, 4), Point3D(2, 2, 2))
        >>> a.p1
        Point3D(1, 1, 1)

        """
        return self.args[0]

    @property
    def normal_vector(self):
        """Normal vector of the given plane.

        Examples
        ========

        >>> from sympy import Point3D, Plane
        >>> a = Plane(Point3D(1, 1, 1), Point3D(2, 3, 4), Point3D(2, 2, 2))
        >>> a.normal_vector
        (-1, 2, -1)
        >>> a = Plane(Point3D(1, 1, 1), normal_vector=(1, 4, 7))
        >>> a.normal_vector
        (1, 4, 7)

        """
        return self.args[1]

    def equation(self, x=None, y=None, z=None):
        """The equation of the Plane.

        Examples
        ========

        >>> from sympy import Point3D, Plane
        >>> a = Plane(Point3D(1, 1, 2), Point3D(2, 4, 7), Point3D(3, 5, 1))
        >>> a.equation()
        -23*x + 11*y - 2*z + 16
        >>> a = Plane(Point3D(1, 4, 2), normal_vector=(6, 6, 6))
        >>> a.equation()
        6*x + 6*y + 6*z - 42

        """
        x, y, z = [i if i else Symbol(j, real=True) for i, j in zip((x, y, z), 'xyz')]
        a = Point3D(x, y, z)
        b = self.p1.direction_ratio(a)
        c = self.normal_vector
        return (sum(i*j for i, j in zip(b, c)))

    def projection(self, pt):
        """Projection of any point on the plane
           Will result in a 3D point.

        Parameters
        ==========

        Point or Point3D

        Returns
        =======

        Point3D

        Notes
        =====

        For the interaction between 2D and 3D points, you should convert the 2D
        point to 3D by using this method. For example for finding the distance
        between Point(1, 2) and Point3D(1, 2, 3) you can convert the 2D point
        to Point3D(1, 2, 0) by projecting it to plane which is parallel to xy
        plane.

        Examples
        ========

        >>> from sympy import Plane, Point, Point3D
        >>> a = Plane(Point3D(1, 1, 1), normal_vector=(1, 1, 1))
        >>> b = Point(1, 2)
        >>> a.projection(b)
        Point3D(1, 2, 0)
        >>> c = Point3D(1, 1, 1)
        >>> a.projection(c)
        Point3D(4/3, 4/3, 4/3)

        """
        x, y, z = map(Dummy, 'xyz')
        k = expand_mul(self.equation(x, y, z))
        a, b, c = [k.coeff(i) for i in (x, y, z)]
        d = k.xreplace({x: pt.args[0], y: pt.args[1], z:0})
        t = -d/(a**2 + b**2 + c**2)
        if isinstance(pt, Point):
            return Point3D(pt.x + t*a, pt.y + t*b, t*c)
        if isinstance(pt, Point3D):
            return Point3D(pt.x + t*a, pt.y + t*b, pt.z + t*c)

    def projection_line(self, l):
        """Projection of any line on the plane.
           Result will be a 3D line or Ray or segment.

        Parameters
        ==========

        LinearEntity or LinearEntity3D

        Returns
        =======

        Line3D

        Notes
        =====

        For the interaction between 2D and 3D lines(segments, rays), you should
        convert the line to 3D by using this method. For example for finding the
        intersection between a 2D and a 3D line, convert the 2D line to a 3D line
        by projecting it on a required plane and then proceed to find the
        intersection between those lines.

        Examples
        ========

        >>> from sympy import Plane, Line, Line3D, Point, Point3D
        >>> a = Plane(Point3D(1, 1, 1), normal_vector=(1, 1, 1))
        >>> b = Line(Point(1, 1), Point(2, 2))
        >>> a.projection_line(b)
        Line3D(Point3D(4/3, 4/3, 1/3), Point3D(5/3, 5/3, -1/3))
        >>> c = Line3D(Point3D(1, 1, 1), Point3D(2, 2, 2))
        >>> a.projection_line(c)
        Line3D(Point3D(4/3, 4/3, 4/3), Point3D(5/3, 5/3, 5/3))

        """
        from sympy.geometry.line import LinearEntity
        from sympy.geometry.line3d import LinearEntity3D
        if not isinstance(l, (LinearEntity, LinearEntity3D)):
            raise NotImplementedError('Enter a linear entity only')
        a, b = self.projection(l.p1), self.projection(l.p2)
        if isinstance(l, (Line, Line3D)):
            return Line3D(a, b)
        if isinstance(l, (Ray, Ray3D)):
            return Ray3D(a, b)
        if isinstance(l, (Segment, Segment3D)):
            return Segment3D(a, b)

    def is_parallel(self, l):
        """Is the given geometric entity parallel to the plane?

        Parameters
        ==========

        LinearEntity3D or Plane

        Returns
        =======

        Boolean

        Examples
        ========

        >>> from sympy import Plane, Point3D
        >>> a = Plane(Point3D(1,4,6), normal_vector=(2, 4, 6))
        >>> b = Plane(Point3D(3,1,3), normal_vector=(4, 8, 12))
        >>> a.is_parallel(b)
        True

        """
        from sympy.geometry.line3d import LinearEntity3D
        if isinstance(l, LinearEntity3D):
            a = l.direction_ratio
            b = self.normal_vector
            c = sum([i*j for i, j in zip(a, b)])
            if c == 0:
                return True
            else:
                return False
        elif isinstance(l, Plane):
            a = Matrix(l.normal_vector)
            b = Matrix(self.normal_vector)
            if a.cross(b).is_zero:
                return True
            else:
                return False

    def is_perpendicular(self, l):
        """is the given geometric entity perpendicualar to the given plane?

        Parameters
        ==========

        LinearEntity3D or Plane

        Returns
        =======

        Boolean

        Examples
        ========

        >>> from sympy import Plane, Point3D
        >>> a = Plane(Point3D(1,4,6), normal_vector=(2, 4, 6))
        >>> b = Plane(Point3D(2, 2, 2), normal_vector=(-1, 2, -1))
        >>> a.is_perpendicular(b)
        True

        """
        from sympy.geometry.line3d import LinearEntity3D
        if isinstance(l, LinearEntity3D):
            a = Matrix(l.direction_ratio)
            b = Matrix(self.normal_vector)
            if a.cross(b).is_zero:
                return True
            else:
                return False
        elif isinstance(l, Plane):
           a = Matrix(l.normal_vector)
           b = Matrix(self.normal_vector)
           if a.dot(b) == 0:
               return True
           else:
               return False
        else:
            return False

    def distance(self, o):
        """Distance beteen the plane and another geometric entity.

        Parameters
        ==========

        Point3D, LinearEntity3D, Plane.

        Returns
        =======

        distance

        Notes
        =====

        This method accepts only 3D entities as it's parameter, but if you want
        to calculate the distance between a 2D entity and a plane you should
        first convert to a 3D entity by projecting onto a desired plane and
        then proceed to calculate the distance.

        Examples
        ========

        >>> from sympy import Point, Point3D, Line, Line3D, Plane
        >>> a = Plane(Point3D(1, 1, 1), normal_vector=(1, 1, 1))
        >>> b = Point3D(1, 2, 3)
        >>> a.distance(b)
        sqrt(3)
        >>> c = Line3D(Point3D(2, 3, 1), Point3D(1, 2, 2))
        >>> a.distance(c)
        0

        """
        from sympy.geometry.line3d import LinearEntity3D
        x, y, z = map(Dummy, 'xyz')
        if self.intersection(o) != []:
            return S.Zero

        if isinstance(o, Point3D):
           x, y, z = map(Dummy, 'xyz')
           k = self.equation(x, y, z)
           a, b, c = [k.coeff(i) for i in (x, y, z)]
           d = k.xreplace({x: o.args[0], y: o.args[1], z: o.args[2]})
           t = abs(d/sqrt(a**2 + b**2 + c**2))
           return t
        if isinstance(o, LinearEntity3D):
            a, b = o.p1, self.p1
            c = Matrix(a.direction_ratio(b))
            d = Matrix(self.normal_vector)
            e = c.dot(d)
            f = sqrt(sum([i**2 for i in self.normal_vector]))
            return abs(e / f)
        if isinstance(o, Plane):
            a, b = o.p1, self.p1
            c = Matrix(a.direction_ratio(b))
            d = Matrix(self.normal_vector)
            e = c.dot(d)
            f = sqrt(sum([i**2 for i in self.normal_vector]))
            return abs(e / f)

    def angle_between(self, o):
        """Angle between the plane and other geometric entity.

        Parameters
        ==========

        LinearEntity3D, Plane.

        Returns
        =======

        angle : angle in radians

        Notes
        =====

        This method accepts only 3D entities as it's parameter, but if you want
        to calculate the angle between a 2D entity and a plane you should
        first convert to a 3D entity by projecting onto a desired plane and
        then proceed to calculate the angle.

        Examples
        ========

        >>> from sympy import Point3D, Line3D, Plane
        >>> a = Plane(Point3D(1, 2, 2), normal_vector=(1, 2, 3))
        >>> b = Line3D(Point3D(1, 3, 4), Point3D(2, 2, 2))
        >>> a.angle_between(b)
        -asin(sqrt(21)/6)

        """
        from sympy.geometry.line3d import LinearEntity3D
        if isinstance(o, LinearEntity3D):
            a = Matrix(self.normal_vector)
            b = Matrix(o.direction_ratio)
            c = a.dot(b)
            d = sqrt(sum([i**2 for i in self.normal_vector]))
            e = sqrt(sum([i**2 for i in o.direction_ratio]))
            return C.asin(c/(d*e))
        if isinstance(o, Plane):
            a = Matrix(self.normal_vector)
            b = Matrix(o.normal_vector)
            c = a.dot(b)
            d = sqrt(sum([i**2 for i in self.normal_vector]))
            e = sqrt(sum([i**2 for i in o.normal_vector]))
            return C.acos(c/(d*e))


    @staticmethod
    def are_concurrent(*planes):
        """Is a sequence of Planes concurrent?

        Two or more Planes are concurrent if their intersections
        are a common line.

        Parameters
        ==========

        planes: list

        Returns
        =======

        Boolean

        Examples
        ========

        >>> from sympy import Plane, Point3D
        >>> a = Plane(Point3D(5, 0, 0), normal_vector=(1, -1, 1))
        >>> b = Plane(Point3D(0, -2, 0), normal_vector=(3, 1, 1))
        >>> c = Plane(Point3D(0, -1, 0), normal_vector=(5, -1, 9))
        >>> Plane.are_concurrent(a, b)
        True
        >>> Plane.are_concurrent(a, b, c)
        False

        """
        if len(planes) < 2:
            return False
        for i in planes:
            if not isinstance(i, Plane):
                raise ValueError('All objects should be Planes but got %s' % i.func)
        planes = list(planes)
        first = planes.pop(0)
        sol = first.intersection(planes[0])
        if sol == []:
            return False
        else:
            line = sol[0]
            for i in planes[1:]:
                l = first.intersection(i)
                if not l or not l[0] in line:
                    return False
            return True

    def perpendicular_line(self, pt):
        """A line perpendicular to the given plane.

        Parameters
        ==========

        pt: Point3D

        Returns
        =======

        Line3D

        Examples
        ========

        >>> from sympy import Plane, Point3D, Line3D
        >>> a = Plane(Point3D(1,4,6), normal_vector=(2, 4, 6))
        >>> a.perpendicular_line(Point3D(9, 8, 7))
        Line3D(Point3D(9, 8, 7), Point3D(11, 12, 13))

        """
        a = self.normal_vector
        return Line3D(pt, direction_ratio=a)

    def parallel_plane(self, pt):
        """
        Plane parallel to the given plane and passing through the point pt.

        Parameters
        ==========

        pt: Point3D

        Returns
        =======

        Plane

        Examples
        ========

        >>> from sympy import Plane, Point3D
        >>> a = Plane(Point3D(1, 4, 6), normal_vector=(2, 4, 6))
        >>> a.parallel_plane(Point3D(2, 3, 5))
        Plane(Point3D(2, 3, 5), (2, 4, 6))

        """
        a = self.normal_vector
        return Plane(pt, normal_vector=a)

    def perpendicular_plane(self, l, pt):
        """
        Plane perpendicular to the given plane and passing through a line in the
        given plane.

        Parameters
        ==========

        pt: Point3D
        l: LinearEntity or LinearEntity3D

        Returns
        =======

        Plane

        Examples
        ========

        >>> from sympy import Plane, Point3D, Line3D
        >>> a = Plane(Point3D(1, 4, 6), normal_vector=(2, 4, 6))
        >>> b = Line3D(Point3D(-27, 27, 0), Point3D(-37, 35, -2))
        >>> a.perpendicular_plane(b, Point3D(1, 2, 3))
        Plane(Point3D(1, 2, 3), (-26, -26, 26))

        """
        if l in self:
            a = l
        else:
            a = self.projection_line(l)
        b, c, d = a.p1, a.p2, pt
        e = Plane(d, b, c)
        if e.is_perpendicular(self):
            return e
        else:
            return NotImplementedError('The plane passing through the given line'
                                        'is not perpendicular to the given plane')

    def random_point(self, seed=None):
        """ Returns a random point on the Plane.

        Returns
        =======

        Point3D

        """
        import random
        if seed is not None:
            rng = random.Random(seed)
        else:
            rng = random
        t = Dummy('t')
        return self.arbitrary_point(t).subs(t, Rational(rng.random()))

    def arbitrary_point(self, t=None):
        """ Returns an arbitrary point on the Plane; varying `t` from 0 to 2*pi
        will move the point in a circle of radius 1 about p1 of the Plane.

        Examples
        ========

        >>> from sympy.geometry.plane import Plane
        >>> from sympy.abc import t
        >>> p = Plane((0, 0, 0), (0, 0, 1), (0, 1, 0))
        >>> p.arbitrary_point(t)
        Point3D(0, cos(t), sin(t))
        >>> _.distance(p.p1).simplify()
        1

        Returns
        =======

        Point3D

        """
        from sympy import cos, sin
        t = t or Dummy('t')
        x, y, z = self.normal_vector
        a, b, c = self.p1.args
        if x == y == 0:
            return Point3D(a + cos(t), b + sin(t), c)
        elif x == z == 0:
            return Point3D(a + cos(t), b, c + sin(t))
        elif y == z == 0:
            return Point3D(a, b + cos(t), c + sin(t))
        m = Dummy()
        p = self.projection(Point3D(self.p1.x + cos(t), self.p1.y + sin(t), 0)*m)
        return p.xreplace({m: solve(p.distance(self.p1) - 1, m)[0]})

    def intersection(self, o):
        """ The intersection with other geometrical entity.

        Parameters
        ==========

        Point, Point3D, LinearEntity, LinearEntity3D, Plane

        Returns
        =======

        List

        Examples
        ========

        >>> from sympy import Point, Point3D, Line, Line3D, Plane
        >>> a = Plane(Point3D(1, 2, 3), normal_vector=(1, 1, 1))
        >>> b = Point3D(1, 2, 3)
        >>> a.intersection(b)
        [Point3D(1, 2, 3)]
        >>> c = Line3D(Point3D(1, 4, 7), Point3D(2, 2, 2))
        >>> a.intersection(c)
        [Point3D(2, 2, 2)]
        >>> d = Plane(Point3D(6, 0, 0), normal_vector=(2, -5, 3))
        >>> e = Plane(Point3D(2, 0, 0), normal_vector=(3, 4, -3))
        >>> d.intersection(e)
        [Line3D(Point3D(78/23, -24/23, 0), Point3D(147/23, 321/23, 23))]

        """
        from sympy.geometry.line3d import LinearEntity3D
        from sympy.geometry.line import LinearEntity
        if isinstance(o, Point) or isinstance(o, Point3D):
            if o in self:
                return [o]
            else:
                return []
        if isinstance(o, LinearEntity3D):
            x, y, z, t = map(Dummy, 'xyzt')
            if o in self:
                return [o]
            else:
                a = o.arbitrary_point(t)
                b = self.equation(x, y, z)
                c = solve(b.subs([(x, a.x), (y, a.y), (z, a.z)]))
                if c == []:
                    return []
                else:
                    a = a.subs(t, c[0])
                    if a in o:
                        return [a]
                    else:
                        return []
        if isinstance(o, LinearEntity):
            t = Dummy('t')
            x, y, z = map(Dummy, 'xyz')
            if o in self:
                return [o]
            else:
                a = self.equation(x, y, z).subs(z,0)
                b = o.equation(x, y)
                c = solve((a,b),[x, y])
                if c is {}:
                    return []
                else:
                    return [Point(c[x], c[y])]
        if isinstance(o, Plane):
            if o == self:
                return [self]
            if self.is_parallel(o):
                return []
            else:
                x, y, z = map(Dummy, 'xyz')
                a, b = Matrix([self.normal_vector]), Matrix([o.normal_vector])
                c = list(a.cross(b))
                d = self.equation(x, y, z)
                e = o.equation(x, y, z)
                f = solve((d.subs(z,0), e.subs(z,0)),[x, y])
                if len(f) == 2:
                    return [Line3D(Point3D(f[x], f[y], 0), direction_ratio=c)]
                g = solve((d.subs(y,0), e.subs(y,0)),[x, z])
                if len(g) == 2:
                    return [Line3D(Point3D(g[x], 0, g[z]), direction_ratio=c)]
                h = solve((d.subs(x,0), e.subs(x,0)),[y, z])
                if len(h) == 2:
                    return [Line3D(Point3D(0, h[y], h[z]), direction_ratio=c)]

    def __contains__(self, o):
        from sympy.geometry.line3d import LinearEntity3D
        from sympy.geometry.line import LinearEntity
        x, y, z = map(Dummy, 'xyz')
        k = self.equation(x, y, z)
        if isinstance(o, Point3D):
            d = k.xreplace(dict(zip((x, y, z), o.args)))
            return Eq(simplify(d), 0)
        elif isinstance(o, Point):
            if self.projection(o) == Point3D(o.x, o.y, 0):
                return True
            else:
                return False
        elif isinstance(o, LinearEntity3D):
            d = o.arbitrary_point()
            e = k.subs([(x, d.x), (y, d.y), (z, d.z)])
            return e == 0
        elif isinstance(o, LinearEntity):
            if self.projection_line(o) == Line3D(Point3D(o.p1.x, o.p1.y, 0),
                                                Point3D(o.p2.x, o.p2.y, 0)):
                return True
            else:
                return False
        else:
            return False

    def is_coplanar(self, o):
        """ Returns True if `o` is coplanar with self, else False.

        Examples
        ========

        >>> from sympy import Plane, Point3D
        >>> o = (0, 0, 0)
        >>> p = Plane(o, (1, 1, 1))
        >>> p2 = Plane(o, (2, 2, 2))
        >>> p == p2
        False
        >>> p.is_coplanar(p2)
        True
        """
        if isinstance(o, Plane):
            x, y, z = map(Dummy, 'xyz')
            return not cancel(self.equation(x, y, z)/o.equation(x, y, z)).has(x, y, z)
        if isinstance(o, Point3D):
            return o in self
        elif isinstance(o, LinearEntity3D):
            return all(i in self for i in self)
        elif isinstance(o, GeometryEntity):  # XXX should only be handling 2D objects now
            return all(i == 0 for i in self.normal_vector[:2])
