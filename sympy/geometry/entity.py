"""The definition of the base geometrical entity with attributes common to
all derived geometrical entities.

Contains
========
GeometryEntity

"""

from sympy.core.compatibility import cmp

# How entities are ordered; used by __cmp__ in GeometryEntity
ordering_of_classes = [
    "Point",
    "Segment",
    "Ray",
    "Line",
    "Triangle",
    "RegularPolygon",
    "Polygon",
    "Circle",
    "Ellipse",
    "Curve"
]

class GeometryEntity(tuple):
    """The base class for all geometrical entities.

    This class doesn't represent any particular geometric entity, it only
    provides the implementation of some methods common to all subclasses.

    """

    def __new__(cls, *args, **kwargs):
        return tuple.__new__(cls, args)

    def __getnewargs__(self):
        return tuple(self)

    @property
    def free_symbols(self):
        """
        Return all but any bound symbols that are used to define the Entity.

        Examples
        ========

        >>> from sympy import Polygon, RegularPolygon, Point
        >>> from sympy.abc import x, y
        >>> t = Polygon(*RegularPolygon(Point(x, y), 1, 3).vertices)
        >>> t.free_symbols
        set([x, y])
        """
        free = set()
        for a in self.args:
            free |= a.free_symbols
        return free

    def intersection(self, o):
        """
        Returns a list of all of the intersections of self with o.

        Notes
        =====

        An entity is not required to implement this method.

        If two different types of entities can intersect, the item with
        higher index in ordering_of_classes should implement
        intersections with anything having a lower index.

        See Also
        ========

        sympy.geometry.util.intersection

        """
        raise NotImplementedError()

    def rotate(self, angle, pt=None):
        """Rotate the object about pt by the given angle (in radians).

        The default pt is the origin, Point(0, 0)

        XXX geometry needs a modify_points method which operates
        on only the points of the object

        See Also
        ========

        scale, translate

        Examples
        ========

        >>> from sympy import Point, RegularPolygon, Polygon, pi
        >>> t = Polygon(*RegularPolygon(Point(0, 0), 1, 3).vertices)
        >>> t # vertex on x axis
        Triangle(Point(1, 0), Point(-1/2, sqrt(3)/2), Point(-1/2, -sqrt(3)/2))
        >>> t.rotate(pi/2) # vertex on y axis now
        Triangle(Point(0, 1), Point(-sqrt(3)/2, -1/2), Point(sqrt(3)/2, -1/2))

        """
        from sympy import cos, sin, Point

        c = cos(angle)
        s = sin(angle)

        if isinstance(self, Point):
            rv = self
            if pt is not None:
                rv -= pt
            x, y = rv
            rv = Point(c*x-s*y, s*x+c*y)
            if pt is not None:
                rv += pt
            return rv

        newargs = []
        for a in self.args:
            if isinstance(a, GeometryEntity):
                newargs.append(a.rotate(angle, pt))
            else:
                newargs.append(a)
        return type(self)(*newargs)

    def scale(self, x=1, y=1):
        """Scale the object by multiplying the x,y-coordinates by x and y.

        See Also
        ========

        rotate, translate

        Examples
        ========

        >>> from sympy import RegularPolygon, Point, Polygon
        >>> t = Polygon(*RegularPolygon(Point(0, 0), 1, 3).vertices)
        >>> t
        Triangle(Point(1, 0), Point(-1/2, sqrt(3)/2), Point(-1/2, -sqrt(3)/2))
        >>> t.scale(2)
        Triangle(Point(2, 0), Point(-1, sqrt(3)/2), Point(-1, -sqrt(3)/2))
        >>> t.scale(2,2)
        Triangle(Point(2, 0), Point(-1, sqrt(3)), Point(-1, -sqrt(3)))

        """
        from sympy import Point
        if isinstance(self, Point):
            return Point(self[0]*x, self[1]*y)
        newargs = []
        for a in self.args:
            if isinstance(a, GeometryEntity):
                newargs.append(a.scale(x, y))
            else:
                newargs.append(a)
        return type(self)(*newargs)

    def translate(self, x=0, y=0):
        """Shift the object by adding to the x,y-coordinates the values x and y.

        See Also
        ========

        rotate, scale

        Examples
        ========

        >>> from sympy import RegularPolygon, Point, Polygon
        >>> t = Polygon(*RegularPolygon(Point(0, 0), 1, 3).vertices)
        >>> t
        Triangle(Point(1, 0), Point(-1/2, sqrt(3)/2), Point(-1/2, -sqrt(3)/2))
        >>> t.translate(2)
        Triangle(Point(3, 0), Point(3/2, sqrt(3)/2), Point(3/2, -sqrt(3)/2))
        >>> t.translate(2,2)
        Triangle(Point(3, 2), Point(3/2, sqrt(3)/2 + 2), Point(3/2, -sqrt(3)/2 + 2))

        """
        from sympy import Point
        if not isinstance(x, Point):
            pt = Point(x, y)
        else:
            pt = x
        if isinstance(self, Point):
            return self + pt
        newargs = []
        for a in self.args:
            if isinstance(a, GeometryEntity):
                newargs.append(a.translate(pt))
            else:
                newargs.append(a)
        return type(self)(*newargs)

    def encloses(self, o):
        """
        Return True if o is inside (not on or outside) the boundaries of self.

        The object will be decomposed into Points and individual Entities need
        only define an encloses_point method for their class.

        See Also
        ========

        sympy.geometry.ellipse.Ellipse.encloses_point
        sympy.geometry.polygon.Polygon.encloses_point

        Examples
        ========

        >>> from sympy import RegularPolygon, Point, Polygon
        >>> t  = Polygon(*RegularPolygon(Point(0, 0), 1, 3).vertices)
        >>> t2 = Polygon(*RegularPolygon(Point(0, 0), 2, 3).vertices)
        >>> t2.encloses(t)
        True
        >>> t.encloses(t2)
        False
        """
        from sympy.geometry.point import Point
        from sympy.geometry.line import Segment, Ray, Line
        from sympy.geometry.ellipse import Ellipse
        from sympy.geometry.polygon import Polygon, RegularPolygon

        if isinstance(o, Point):
            return self.encloses_point(o)
        elif isinstance(o, Segment):
            return all(self.encloses_point(x) for x in o.points)
        elif isinstance(o, Ray) or isinstance(o, Line):
            return False
        elif isinstance(o, Ellipse):
            return self.encloses_point(o.center) and not self.intersection(o)
        elif isinstance(o, Polygon):
            if isinstance(o, RegularPolygon):
                if not self.encloses_point(o.center):
                    return False
            return all(self.encloses_point(v) for v in o.vertices)
        raise NotImplementedError()

    def is_similar(self, other):
        """Is this geometrical entity similar to another geometrical entity?

        Two entities are similar if a uniform scaling (enlarging or
        shrinking) of one of the entities will allow one to obtain the other.

        Notes
        =====

        This method is not intended to be used directly but rather
        through the `are_similar` function found in util.py.
        An entity is not required to implement this method.
        If two different types of entities can be similar, it is only
        required that one of them be able to determine this.

        See Also
        ========

        scale

        """
        raise NotImplementedError()

    def subs(self, *args):
        """
        Substitues new for old in self.

        See Also
        ========

        sympy.core.basic.subs

        Examples
        ========

        >>> from sympy import Point, Circle
        >>> from sympy.abc import x, y, z
        >>> c = Circle(Point(0, 0), 3)
        >>> c.subs(0, 5)
        Circle(Point(5, 5), 3)
        >>> c.subs(3, z)
        Circle(Point(0, 0), z)

        """
        return type(self)(*[a.subs(*args) for a in self.args])

    def _eval_subs(self, old, new):
        return type(self)(*[a.subs(old, new) for a in self.args])

    @property
    def args(self):
        """Return whatever is contained in the object's tuple.

        The contents will not necessarily be Points. This is also
        what will be returned when one does "for x in self".

        Examples
        ========

        >>> from sympy import RegularPolygon, Point, Polygon
        >>> t = Polygon(*RegularPolygon(Point(0, 0), 1, 3).vertices)
        >>> t.args
        (Point(1, 0), Point(-1/2, sqrt(3)/2), Point(-1/2, -sqrt(3)/2))
        """

        return tuple(self)

    def __ne__(self, o):
        """Test inequality of two geometrical entities."""
        return not self.__eq__(o)

    def __radd__(self, a):
        return a.__add__(self)

    def __rsub__(self, a):
        return a.__sub__(self)

    def __rmul__(self, a):
        return a.__mul__(self)

    def __rdiv__(self, a):
        return a.__div__(self)

    def __str__(self):
        """String representation of a GeometryEntity."""
        from sympy.printing import sstr
        return type(self).__name__ + sstr(tuple(self))

    def __repr__(self):
        """String representation of a GeometryEntity that can be evaluated
        by sympy."""
        return type(self).__name__ + repr(tuple(self))

    def __cmp__(self, other):
        """Comparison of two GeometryEntities."""
        n1 = self.__class__.__name__
        n2 = other.__class__.__name__
        c = cmp(n1, n2)
        if not c:
            return 0

        i1 = -1
        for cls in self.__class__.__mro__:
            try:
                i1 = ordering_of_classes.index(cls.__name__)
                break
            except ValueError:
                i1 = -1
        if i1 == -1:
            return c

        i2 = -1
        for cls in other.__class__.__mro__:
            try:
                i2 = ordering_of_classes.index(cls.__name__)
                break
            except ValueError:
                i2 = -1
        if i2 == -1:
            return c

        return cmp(i1, i2)

    def __contains__(self, other):
        """Subclasses should implement this method for anything more complex than equality."""
        if type(self) == type(other):
            return self == other
        raise NotImplementedError()

from sympy.core.sympify import converter
converter[GeometryEntity] = lambda x: x
