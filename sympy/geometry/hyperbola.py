"""Hyperbolic geometrical entities.

Contains
* Hyperbola

"""

from sympy.core import S, sympify
from sympy.core.evalf import N
from sympy.core.logic import fuzzy_bool
from sympy.core.numbers import Rational, oo
from sympy.core.sorting import ordered
from sympy.core.symbol import Dummy, uniquely_named_symbol, _symbol, symbols
from sympy.simplify import simplify, trigsimp
from sympy.functions.elementary.miscellaneous import sqrt, Max
from sympy.functions.elementary.trigonometric import sec, tan
from .entity import GeometryEntity, GeometrySet
from .exceptions import GeometryError
from .ellipse import Circle, Ellipse
from .line import Line, Segment, Ray2D, Segment2D, Line2D, LinearEntity3D
from .point import Point, Point2D, Point3D
from .polygon import Polygon
from sympy.solvers import solve
from sympy.utilities.misc import filldedent, func_name


class Hyperbola(GeometrySet):
    """A hyperbolic GeometryEntity.

    Parameters
    ==========

    center : Point, optional
        Default value is Point(0, 0)
    hradius : number or SymPy expression, optional
    vradius : number or SymPy expression, optional
    eccentricity : number or SymPy expression, optional
        Two of `hradius`, `vradius` and `eccentricity` must be supplied to
        create an Hyperbola. The third is derived from the two supplied.

    Attributes
    ==========

    center
    hradius
    vradius
    eccentricity
    focus_distance
    foci

    Raises
    ======

    GeometryError
        When `hradius`, `vradius` and `eccentricity` are incorrectly supplied
        as parameters.
    TypeError
        When `center` is not a Point.

    See Also
    ========

    Ellipse

    Notes
    -----
    Constructed from a center and two radii, the first being the horizontal
    radius (along the x-axis) and the second being the vertical radius (along
    the y-axis).

    When symbolic value for hradius and vradius are used, any calculation that
    refers to the foci or the major or minor axis will assume that the hyperbola
    has its major radius on the x-axis. If this is not true then a manual
    rotation is necessary.

    Examples
    ========

    >>> from sympy import Hyperbola, Point
    >>> h1 = Hyperbola(Point(0, 0), 5, 3)
    >>> h1.hradius, h1.vradius
    (5, 3)
    >>> h2 = Hyperbola(Point(3, 1), hradius=3, eccentricity= 2)
    >>> h2
    Hyperbola(Point2D(3, 1), 3, 3*sqrt(3))

    """

    def __new__(
        cls, center=None, hradius=None, vradius=None, eccentricity=None, **kwargs):
        hradius = sympify(hradius)
        vradius = sympify(vradius)

        eccentricity = sympify(eccentricity)

        if center is None:
            center = Point(0, 0)
        else:
            center = Point(center, dim=2)

        if len(center) != 2:
            raise ValueError('The center of "{}" must be a two dimensional point'.format(cls))

        if len(list(filter(lambda x: x is not None, (hradius, vradius, eccentricity)))) != 2:
            raise ValueError(filldedent('''
                Exactly two arguments of "hradius", "vradius", and
                "eccentricity" must not be None.'''))

        if eccentricity is not None:
            if eccentricity <= 1:
                raise GeometryError("Eccentricity of hyperbola is always greater than 1")
            elif hradius is None:
                hradius = vradius / sqrt(eccentricity**2 - 1)
            elif vradius is None:
                vradius = hradius * sqrt(eccentricity**2 - 1)

        #if hradius == vradius:
        #    TODO: Implement Rectangular Hyperbola class for representing standard Rectangular Hyperbola Eq(x*y - c**2, 0)
        #    return Rectangular_Hyperbola(center, hradius, **kwargs)

        if S.Zero in (hradius, vradius):
            return Segment(Point(center[0] - hradius, center[1] - vradius), Point(center[0] + hradius, center[1] + vradius))

        if hradius.is_real is False or vradius.is_real is False:
            raise GeometryError("Invalid value encountered when computing hradius / vradius.")

        return GeometryEntity.__new__(cls, center, hradius, vradius, **kwargs)

    def __contains__(self, o):
        if isinstance(o, Point):
            x = Dummy('x', real=True)
            y = Dummy('y', real=True)

            res = self.equation(x, y).subs({x: o.x, y: o.y})
            return trigsimp(simplify(res)).is_zero
        elif isinstance(o, Hyperbola):
            return self == o
        return False

    def __eq__(self, o):
        """Is the other GeometryEntity the same as this hyperbola?"""
        return isinstance(o, Hyperbola) and (self.center == o.center and
                                           self.hradius == o.hradius and
                                           self.vradius == o.vradius and
                                           self.foci == o.foci)

    def __hash__(self):
        return super().__hash__()

    def _svg(self, scale_factor=1., fill_color="#66cc99"):
        """Returns SVG hyperbola element for the Hyperbola.

        Parameters
        ==========

        scale_factor : float
            Multiplication factor for the SVG stroke-width.  Default is 1.
        fill_color : str, optional
            Hex string for fill color. Default is "#66cc99".
        """

        c = N(self.center)
        h, v = N(self.hradius), N(self.vradius)
        return (
            '<hyperbola fill="{1}" stroke="#555555" '
            'stroke-width="{0}" opacity="0.6" cx="{2}" cy="{3}" rx="{4}" ry="{5}"/>'
        ).format(2. * scale_factor, fill_color, c.x, c.y, h, v)

    @property
    def ambient_dimension(self):
        return 2

    @property
    def center(self):
        """The center of the hyperbola.

        Returns
        =======

        center : number

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(Point(0, 0), 5, 3)
        >>> h1.center
        Point2D(0, 0)

        """
        return self.args[0]

    @property
    def eccentricity(self):
        """The eccentricity of the hyperbola.

        Returns
        =======

        eccentricity : number

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(Point(0, 0), 5, 3)
        >>> h1.eccentricity
        sqrt(34)/5

        """
        return self.focus_distance / self.major

    def equation(self, x='x', y='y', _slope=None, branch=None):
        """
        Returns the equation of a hyperbola aligned with the x and y axes;
        when slope is given, the equation returned corresponds to a hyperbola
        with a major axis having that slope.

        Parameters
        ==========

        x : str, optional
            Label for the x-axis. Default value is 'x'.
        y : str, optional
            Label for the y-axis. Default value is 'y'.
        _slope : Expr, optional
                The slope of the major axis. Ignored when 'None'.

        Returns
        =======

        equation : SymPy expression

        See Also
        ========

        arbitrary_point : Returns parameterized point on hyperbola

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> from sympy.abc import x, y
        >>> h1 = Hyperbola(Point(1, 0), 3, 2)
        >>> eq1 = h1.equation(x, y); eq1
        -y**2/4 + (x/3 - 1/3)**2 - 1
        >>> eq2 = h1.equation(x, y, _slope=1); eq2
        -(-x + y + 1)**2/8 + (x + y - 1)**2/18 - 1

        A point on h1 satisfies eq1. Let's use one on the x-axis:

        >>> p1 = h1.center + Point(h1.major, 0)
        >>> assert eq1.subs(x, p1.x).subs(y, p1.y) == 0
        >>> h1.equation(x, y, branch = 'right')
        x - 3*sqrt(y**2 + 4)/2 - 1
        >>> h1.equation(x, y, branch = 'left')
        x + 3*sqrt(y**2 + 4)/2 - 1
        """

        x = _symbol(x, real=True)
        y = _symbol(y, real=True)

        dx = x - self.center.x
        dy = y - self.center.y

        if branch is None:
            if _slope is not None:
                L = (dy - _slope * dx) ** 2
                l = (_slope * dy + dx) ** 2
                h = 1 + _slope ** 2
                b = h * self.major ** 2
                a = h * self.minor ** 2
                expr = l / b - L / a - 1
            else:
                t1 = (dx/self.hradius)**2
                t2 = (dy/self.vradius)**2
                expr = t1 - t2 - 1
            return expr
        else:
            if _slope is None:
                if str(branch) == 'left':
                    expr = dx + (self.hradius * sqrt((dy)**2 + (self.vradius)**2) / self.vradius)
                if str(branch) == 'right':
                    expr = dx - (self.hradius * sqrt((dy)**2 + (self.vradius)**2) / self.vradius)
                return expr
            else:
                #TODO :Implement branches when hyperbola is constructed using _slope
                pass

    def encloses_point(self, p):
        """
        Return True if p is enclosed by (is inside of) self.

        Notes
        -----
        Being on the border of self is considered False.

        Parameters
        ==========

        p : Point

        Returns
        =======

        encloses_point : True, False or None

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy import Hyperbola, S
        >>> from sympy.abc import t
        >>> h = Hyperbola((0, 0), 3, 2)
        >>> h.encloses_point((0, 0))
        True
        >>> h.encloses_point(h.arbitrary_point(t).subs(t, S.Half))
        False
        >>> h.encloses_point((4, 0))
        False

        """
        p = Point(p, dim=2)
        if p in self:
            return False

        if len(self.foci) == 2:
            x, y = symbols('x y')
            expr = self.equation(x, y) + 1
            syms = list(self.equation().free_symbols)
            expr = expr.xreplace({syms[0]: x, syms[1]: y})
            test = expr.subs([(x, p.x), (y, p.y)])
            if test >= S.One:
                return False
            else:
                return True
        else:
            test = self.radius - self.center.distance(p)

        return fuzzy_bool(test.is_negative)

    @property
    def foci(self):
        """The foci of the hyperbola.

        Notes
        -----
        The foci can only be calculated if the major/minor axes are known.

        Raises
        ======

        ValueError
            When the major and minor axis cannot be determined.

        See Also
        ========

        sympy.geometry.point.Point
        focus_distance : Returns the distance between focus and center

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(p1, 3, 1)
        >>> h1.foci
        (Point2D(-sqrt(10), 0), Point2D(sqrt(10), 0))

        """
        c = self.center
        hr, vr = self.hradius, self.vradius

        fd = sqrt(hr**2 + vr**2)
        return (c + Point(-fd, 0), c + Point(fd, 0))

    @property
    def focus_distance(self):
        """The focal distance of the hyperbola.

        The distance between the center and one focus.

        Returns
        =======

        focus_distance : number

        See Also
        ========

        foci

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(p1, 3, 1)
        >>> h1.focus_distance
        sqrt(10)

        """
        return Point.distance(self.center, self.foci[0])

    @property
    def hradius(self):
        """The horizontal radius of the hyperbola.

        Returns
        =======

        hradius : number

        See Also
        ========

        vradius, major, minor

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(p1, 3, 1)
        >>> h1.hradius
        3

        """
        return self.args[1]

    @property
    def vradius(self):
        """The vertical radius of the hyperbola.

        Returns
        =======

        vradius : number

        See Also
        ========

        hradius, major, minor

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(p1, 3, 1)
        >>> h1.vradius
        1

        """
        return self.args[2]

    @property
    def major(self):
        """Longer axis of the hyperbola (if it can be determined) else hradius.

        Returns
        =======

        major : number or expression

        See Also
        ========

        hradius, vradius, minor

        Examples
        ========

        >>> from sympy import Point, Hyperbola, Symbol
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(p1, 3, 1)
        >>> h1.major
        3

        >>> a = Symbol('a')
        >>> b = Symbol('b')
        >>> Hyperbola(p1, a, b).major
        a
        >>> Hyperbola(p1, b, a).major
        b

        >>> m = Symbol('m')
        >>> M = m + 1
        >>> Hyperbola(p1, m, M).major
        m + 1

        """
        ab = self.args[1:3]
        if len(ab) == 1:
            return ab[0]
        a, b = ab
        o = b - a < 0
        if o == True:
            return a
        elif o == False:
            return b
        return self.hradius

    @property
    def minor(self):
        """Shorter axis of the hyperbola (if it can be determined) else vradius.

        Returns
        =======

        minor : number or expression

        See Also
        ========

        hradius, vradius, major

        Examples
        ========

        >>> from sympy import Point, Hyperbola, Symbol
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(p1, 3, 1)
        >>> h1.minor
        1

        >>> a = Symbol('a')
        >>> b = Symbol('b')
        >>> Hyperbola(p1, a, b).minor
        b
        >>> Hyperbola(p1, b, a).minor
        a

        >>> m = Symbol('m')
        >>> M = m + 1
        >>> Hyperbola(p1, m, M).minor
        m

        """
        ab = self.args[1:3]
        if len(ab) == 1:
            return ab[0]
        a, b = ab
        o = a - b < 0
        if o == True:
            return a
        elif o == False:
            return b
        return self.vradius

    @property
    def semilatus_rectum(self):
        """
        Calculates the semi-latus rectum of the hyperbola.

        Semi-latus rectum is defined as one half of the chord through a
        focus parallel to the conic section directrix of a conic section.

        Returns
        =======

        semilatus_rectum : number

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(p1, 3, 1)
        >>> h1.semilatus_rectum
        1/3

        References
        ==========

        .. [1] http://mathworld.wolfram.com/SemilatusRectum.html
        .. [2] https://en.wikipedia.org/wiki/Hyperbola#Semi-latus_rectum

        """
        return self.major * (self.eccentricity ** 2 - 1)

    def asymptote(self, slope="+"):
        """
        Returns asymptote of the hyperbola. By deafult returns the asymptote with positive slope

        An asymptote to a curve is defined as a line such that the distance between
        the curve and the line approaches zero as one or both of the x or y coordinates tends
        to infinity.

        Returns
        =======

        asymptote : expression

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(p1, 3, 1)
        >>> h1.asymptote()
        -x/3 + y
        >>> h1.asymptote(slope = '-')
        x/3 + y
        >>> h2 = Hyperbola(Point(1, 0), 3, 2)
        >>> h2.asymptote()
        -2*x/3 + y + 2/3
        >>> h2.asymptote(slope = '-')
        2*x/3 + y - 2/3

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Hyperbola#Asymptotes

        """

        x = _symbol('x', real=True)
        y = _symbol('y', real=True)

        dx = x - self.center.x
        dy = y - self.center.y

        if str(slope) == "-":
            asymptote = dy + (self.vradius * dx / self.hradius)
        else:
            asymptote = dy - (self.vradius * dx / self.hradius)
        return asymptote

    def arbitrary_point(self, parameter='t'):
        """A parameterized point on the hyperbola.

        Parameters
        ==========

        parameter : str, optional
            Default value is 't'.

        Returns
        =======

        arbitrary_point : Point

        Raises
        ======

        ValueError
            When `parameter` already appears in the functions.

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> h1 = Hyperbola(Point(0, 0), 5, 3)
        >>> h1.arbitrary_point()
        Point2D(5*sec(t), 3*tan(t))

        """
        t = _symbol(parameter, real=True)
        if t.name in (f.name for f in self.free_symbols):
            raise ValueError(filldedent('Symbol %s already appears in object '
                                        'and cannot be used as a parameter.' % t.name))
        return Point(self.center.x + self.hradius*sec(t),
                     self.center.y + self.vradius*tan(t))

    def auxiliary_circle(self):
        """Returns a Circle whose diameter is the major axis of the hyperbola.

        Examples
        ========

        >>> from sympy import Hyperbola, Point, symbols
        >>> c = Point(1, 2)
        >>> Hyperbola(c, 8, 7).auxiliary_circle()
        Circle(Point2D(1, 2), 8)
        >>> a, b = symbols('a b')
        >>> Hyperbola(c, a, b).auxiliary_circle()
        Circle(Point2D(1, 2), Max(a, b))
        """
        return Circle(self.center, Max(self.hradius, self.vradius))

    def director_circle(self):
        """
        Returns a Circle consisting of all points where two perpendicular
        tangent lines to the hyperbola cross each other.

        Returns
        =======

        Circle
            A director circle returned as a geometric object.

        Examples
        ========

        >>> from sympy import Hyperbola, Point, symbols
        >>> c = Point(3,8)
        >>> Hyperbola(c, 5, 3).director_circle()
        Circle(Point2D(3, 8), 4)
        >>> a, b = symbols('a b')
        >>> Hyperbola(c, a, b).director_circle()
        Circle(Point2D(3, 8), sqrt(a**2 - b**2))

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Director_circle

        """
        return Circle(self.center, sqrt(self.hradius**2 - self.vradius**2))

    def evolute(self, x='x', y='y'):
        """The equation of evolute of the hyperbola.

        Parameters
        ==========

        x : str, optional
            Label for the x-axis. Default value is 'x'.
        y : str, optional
            Label for the y-axis. Default value is 'y'.

        Returns
        =======

        equation : SymPy expression

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> h1 = Hyperbola(Point(1, 0), 3, 2)
        >>> h1.evolute()
        -2**(2/3)*y**(2/3) + (3*x - 3)**(2/3) - 13**(2/3)
        """
        if len(self.args) != 3:
            raise NotImplementedError('Evolute of arbitrary Hyperbola is not supported.')
        x = _symbol(x, real=True)
        y = _symbol(y, real=True)
        t1 = (self.hradius*(x - self.center.x))**Rational(2, 3)
        t2 = (self.vradius*(y - self.center.y))**Rational(2, 3)
        return t1 - t2 - (self.hradius**2 + self.vradius**2)**Rational(2, 3)

    def intersection(self, o):
        """The intersection of this hyperbola and another geometrical entity
        `o`.
        Parameters
        ==========
        o : GeometryEntity
        Returns
        =======
        intersection : list of GeometryEntity objects
        Notes
        -----
        Currently supports intersections with Point, Line, Segment, Ray
        Circle, Ellipse and Hyperbola types.
        See Also
        ========
        sympy.geometry.entity.GeometryEntity
        Examples
        ========
        >>> from sympy import Hyperbola, Point, Line, Point2D, Ellipse
        >>> h = Hyperbola(Point(0, 0), 5, 7)
        >>> h.intersection(Point(0, 0))
        []
        >>> h.intersection(Point(5, 0))
        [Point2D(5, 0)]
        >>> h.intersection(Line(Point(0,0), Point(1, 0)))
        [Point2D(-5, 0), Point2D(5, 0)]
        >>> h.intersection(Line(Point(0,0), Point(0, 1)))
        []
        >>> # Checking for intersection with Asymptotes
        >>> h.intersection(Line(Point(0,0), Point(5, 7)))
        []
        >>> h1 = Hyperbola(Point(-1, 0), 4, 3)
        >>> h1.intersection(Hyperbola(Point(1, 0), 4, 3))
        []
        >>> h1.intersection(h)
        [Point2D(-5, 0), Point2D(3245/559, -84*sqrt(755)/559), Point2D(3245/559, 84*sqrt(755)/559)]
        >>> h1.intersection(Ellipse(Point(-1, 0), 4, 3))
        [Point2D(-5, 0), Point2D(3, 0)]
        >>> h1.intersection(Hyperbola(Point2D(0, -1), 4, 3))
        [Point2D(-1/2 + 2*sqrt(4081)/21, -3*sqrt(4081)/56 - 1/2), Point2D(-2*sqrt(4081)/21 - 1/2, -1/2 + 3*sqrt(4081)/56)]
        """
        # TODO: Replace solve with nonlinsolve, when nonlinsolve will be able to solve in real domain
        x = Dummy('x', real=True)
        y = Dummy('y', real=True)

        if isinstance(o, Point):
            if o in self:
                return [o]
            else:
                return []

        elif isinstance(o, (Segment2D, Ray2D)):
            hyperbola_equation = self.equation(x, y)
            result = solve([hyperbola_equation, Line(o.points[0], o.points[1]).equation(x, y)], [x, y])
            return list(ordered([Point(i) for i in result if i in o]))

        elif isinstance(o, Polygon):
            return o.intersection(self)

        elif isinstance(o, (Hyperbola, Line2D, Ellipse)):
            if o == self:
                return self
            else:
                hyperbola_equation = self.equation(x, y)
                return list(ordered([Point(i) for i in solve([hyperbola_equation, o.equation(x, y)], [x, y])]))
        elif isinstance(o, LinearEntity3D):
            raise TypeError('Entity must be two dimensional, not three dimensional')
        else:
            raise TypeError('Intersection not handled for %s' % func_name(o))

    def is_tangent(self, o):
        """Is `o` tangent to the hyperbola?

        Parameters
        ==========

        o : GeometryEntity
            A LinearEntity or Polygon

        Raises
        ======

        NotImplementedError
            When the wrong type of argument is supplied.

        Returns
        =======

        is_tangent: boolean
            True if o is tangent to the hyperbola, False otherwise.

        Examples
        ========

        >>> from sympy import Point, Hyperbola, Line
        >>> p0, p1, p2 = Point(0, 0), Point(3, 0), Point(3, 3)
        >>> h1 = Hyperbola(p0, 3, 2)
        >>> l1 = Line(p1, p2)
        >>> h1.is_tangent(l1)
        True

        """
        if isinstance(o, Point2D):
            return False
        elif isinstance(o, Line2D):
            hit = self.intersection(o)
            if not hit:
                return False
            if len(hit) == 1:
                return True
            # might return None if it can't decide
            return hit[0].equals(hit[1])
        elif isinstance(o, Ray2D):
            intersect = self.intersection(o)
            if len(intersect) == 1:
                return intersect[0] != o.source and not self.encloses_point(o.source)
            else:
                return False
        elif isinstance(o, (Segment2D, Polygon)):
            all_tangents = False
            segments = o.sides if isinstance(o, Polygon) else [o]
            for segment in segments:
                intersect = self.intersection(segment)
                if len(intersect) == 1:
                    if not any(intersect[0] in i for i in segment.points) \
                        and not any(self.encloses_point(i) for i in segment.points):
                        all_tangents = True
                        continue
                    else:
                        return False
                else:
                    return all_tangents
            return all_tangents
        elif isinstance(o, (LinearEntity3D, Point3D)):
            raise TypeError('Entity must be two dimensional, not three dimensional')
        else:
            raise TypeError('Is_tangent not handled for %s' % func_name(o))

    def reflect(self, line):
        """Override GeometryEntity.reflect since the radius
        is not a GeometryEntity.

        Examples
        ========

        >>> from sympy import Hyperbola, Line, Point
        >>> Hyperbola(Point(3, 4), 1, 3).reflect(Line(Point(0, -4), Point(5, 0)))
        Traceback (most recent call last):
        ...
        NotImplementedError:
        General Hyperbola is not supported but the equation of the reflected
        Hyperbola is given by the zeros of: f(x, y) = (9*x/41 + 40*y/41 +
        37/41)**2 - (40*x/123 - 3*y/41 - 364/123)**2 - 1

        Notes
        =====

        Until the general hyperbola (with no axis parallel to the x-axis) is
        supported a NotImplemented error is raised and the equation whose
        zeros define the rotated hyperbola is given.

        """

        if line.slope in (0, oo):
            c = self.center
            c = c.reflect(line)
            return self.func(c, -self.hradius, self.vradius)
        else:
            x, y = [uniquely_named_symbol(
                name, (self, line), modify=lambda s: '_' + s, real=True)
                for name in 'xy']
            expr = self.equation(x, y)
            p = Point(x, y).reflect(line)
            result = expr.subs(zip((x, y), p.args
                                   ), simultaneous=True)
            raise NotImplementedError(filldedent(
                'General Hyperbola is not supported but the equation '
                'of the reflected Hyperbola is given by the zeros of: ' +
                "f(%s, %s) = %s" % (str(x), str(y), str(result))))

    def rotate(self, angle=0, pt=None):
        """Rotate ``angle`` radians counterclockwise about Point ``pt``.

        Note: since the general hyperbola is not supported, only rotations that
        are integer multiples of pi/2 are allowed.

        Examples
        ========

        >>> from sympy import Hyperbola, pi
        >>> Hyperbola((1, 0), 2, 1).rotate(pi/2)
        Hyperbola(Point2D(0, 1), 1, 2)
        >>> Hyperbola((1, 0), 2, 1).rotate(pi)
        Hyperbola(Point2D(-1, 0), 2, 1)
        """
        if self.hradius == self.vradius:
            # TODO : would require Rectangular Hyperbola Class or a Conjugate Hyperbola Class
            pass
        if (angle/S.Pi).is_integer:
            return super().rotate(angle, pt)
        if (2*angle/S.Pi).is_integer:
            return self.func(self.center.rotate(angle, pt), self.vradius, self.hradius)
        raise NotImplementedError('Only rotations of pi/2 are currently supported for Hyperbola.')

    def scale(self, x=1, y=1, pt=None):
        """Override GeometryEntity.scale since it is the major and minor
        axes which must be scaled and they are not GeometryEntities.

        Examples
        ========

        >>> from sympy import Hyperbola, Point
        >>> Hyperbola((0, 0), 2, 1).scale(2, 4)
        Hyperbola(Point2D(0, 0), 4, 4)
        >>> Hyperbola((0, 0), 2, 1).scale(2)
        Hyperbola(Point2D(0, 0), 4, 1)
        >>> Hyperbola((0, 0), 2, 1).scale(2, 3, Point(1, 1))
        Hyperbola(Point2D(-1, -2), 4, 3)
        """
        c = self.center
        if pt:
            pt = Point(pt, dim=2)
            return self.translate(*(-pt).args).scale(x, y).translate(*pt.args)
        h = self.hradius
        v = self.vradius
        return self.func(c.scale(x, y), hradius=h*x, vradius=v*y)
