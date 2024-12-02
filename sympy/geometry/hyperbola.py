"""Hyperbolic geometrical entities.

Contains
* Hyperbola

"""

from __future__ import division, print_function

from sympy import Expr, Eq
from sympy.core import S, pi, sympify
from sympy.core.evaluate import global_evaluate
from sympy.core.logic import fuzzy_bool
from sympy.core.numbers import Rational, oo
from sympy.core.compatibility import ordered
from sympy.core.symbol import Dummy, _uniquely_named_symbol, _symbol
from sympy.simplify import simplify, trigsimp
from sympy.functions.elementary.miscellaneous import sqrt, Max, Min
from sympy.functions.elementary.trigonometric import tan, sec
from sympy.functions.special.elliptic_integrals import elliptic_e
from sympy.geometry.exceptions import GeometryError
from sympy.geometry.line import Ray2D, Segment2D, Line2D, LinearEntity3D
from sympy.polys import DomainError, Poly, PolynomialError
from sympy.polys.polyutils import _not_a_coeff, _nsort
from sympy.solvers import solve
from sympy.solvers.solveset import linear_coeffs
from sympy.utilities.misc import filldedent, func_name

from .entity import GeometryEntity, GeometrySet
from .point import Point, Point2D, Point3D
from .line import Line, Segment
from .util import idiff

import random
from sympy.geometry import Circle, Ellipse
from .polygon import Polygon


class Hyperbola(GeometrySet):
    """
    A Hyperbolic GeometryEntity.

    Parameters
    ==========

    center: Point optional
        Default value is Point(0, 0)
    hradius : number or SymPy expression, optional
    vradius : number or SymPy expression, optional
    eccentricity : number or SymPy expression, optional
        Two of 'hradius', 'vradius' and 'eccentricity' must be supplied to
        create a Hyperbola. The third is derived from the two supplied

    Attributes
    ==========

    center
    hradius
    vradius
    eccentricity
    focus distance
    foci

    TODO
    ====
    Rectangular Hyperbola
    encloses_point
    apoapsis
    periapsis
    is_tangent
    normal_lines
    tangent_lines
    reflect
    rotate

    GeometryError
        When `hradius`, `vradius` and `eccentricity` are incorrectly supplied
        as parameters.
    TypeError
        When `center` is not a Point.

    Examples
    ========

    >>> from sympy import Point, Rational
    >>> from sympy.geometry.hyperbola import Hyperbola
    >>> h1 = Hyperbola(Point(0, 0), 5, 1)
    >>> h1.hradius, h1.vradius
    (5, 1)
    >>> h2 = Hyperbola(Point(3, 1), hradius=3, eccentricity=Rational(5, 4))
    >>> h2
    Hyperbola(Point2D(3, 1), 3, 9/4)

    """

    def __new__(cls, center=None, hradius=None, vradius=None, eccentricity=None, *args, **kwargs):
        hradius = sympify(hradius)
        vradius = sympify(vradius)

        eccentricity = sympify(eccentricity)

        if center is None:
            center = Point(0,0)
        else:
            center = Point(center, dim=2)

        if len(center) != 2:
            raise ValueError('The center of "{0}" must be a two dimensional point'.format(cls))

        if len(list(filter(lambda x: x is not None, (hradius, vradius, eccentricity)))) != 2:
            raise ValueError(filldedent('''
                Exactly two arguments of "hradius", "vradius", and
                "eccentricity" must not be None.'''))

        if eccentricity is not None:
            if hradius is None:
                hradius = vradius / sqrt(eccentricity**2 - 1)
            elif vradius is None:
                vradius = hradius * sqrt(eccentricity**2 - 1)

        return GeometryEntity.__new__(cls, center, hradius, vradius, **kwargs)

    def __eq__(self, o):

        return isinstance(o, Hyperbola) and (self.center == o.center and
                                             self.hradius == o.hradius and
                                             self.vradius == o.vradius)

    def __hash__(self):
        return super(Hyperbola, self).__hash__()

    def __contains__(self, o):
        if isinstance(o, Point):
            x = Dummy('x', real=True)
            y = Dummy('y', real=True)

            res = self.equation(x, y).subs({x: o.x, y: o.y})
            return trigsimp(simplify(res)) is S.Zero
        elif isinstance(o, Hyperbola):
            return self == o
        return False

    def _svg(self, scale_factor=1., fill_color="#66cc99"):
        """Returns SVG ellipse element for the Ellipse.

        Parameters
        ==========

        scale_factor : float
            Multiplication factor for the SVG stroke-width.  Default is 1.
        fill_color : str, optional
            Hex string for fill color. Default is "#66cc99".
        """

        from sympy.core.evalf import N

        c = N(self.center)
        h, v = N(self.hradius), N(self.vradius)
        return (
            'hyperbola fill="{1}" stroke="#555555" '
            'stroke-width="{0}" opacity="0.6" cx="{2}" cy="{3}" rx="{4}" ry="5"/>'
        ).format(2. * scale_factor, fill_color, c.x, c.y, h, v)

    @property
    def ambient_dimension(self):
        return 2

    @property
    def center(self):
        """The center of the hyperbola

        Returns
        =======

        center : number

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy import Point
        >>> from sympy.geometry.hyperbola import Hyperbola
        >>> p1 = Point(0,0)
        >>> h1 = Hyperbola(p1, 3, 1)
        >>> h1.center
        Point2D(0, 0)

        """
        return self.args[0]

    @property
    def hradius(self):

        return self.args[1]

    @property
    def vradius(self):

        return self.args[2]

    @property
    def eccentricity(self):
        """
        The eccentricity of the ellipse.

        Returns
        =======

        eccentricity : number

        Examples
        ========

        >>> from sympy import Point, sqrt
        >>> from sympy.geometry.hyperbola import Hyperbola
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(p1, 3, sqrt(2))
        >>> h1.eccentricity
        sqrt(11)/3

        """
        return self.focus_distance / self.hradius

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

        >>> from sympy import Point
        >>> from sympy.geometry.hyperbola import Hyperbola
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

        >>> from sympy import Point, Symbol
        >>> from sympy.geometry.hyperbola import Hyperbola
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

        rv = Max(*self.args[1:3])
        if isinstance(rv, Max):
            return self.hradius
        return rv

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

        >>> from sympy import Point, Symbol
        >>> from sympy.geometry.hyperbola import Hyperbola
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

        rv = Min(*self.args[1:3])
        if isinstance(rv, Min):
            return self.vradius
        return rv

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

        >>> from sympy import Point
        >>> from sympy.geometry.hyperbola import Hyperbola
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(p1, 3, 1)
        >>> h1.focus_distance
        sqrt(10)

        """

        return Point.distance(self.center, self.foci[0])

    @property
    def semilatus_rectum(self):
        """
        Calculates the semi-latus rectum of the Hyperbola.

        Semi-latus rectum is defined as one half of the the chord through a
        focus parallel to the conic section directrix of a conic section.

        Returns
        =======

        semilatus_rectum : number

        Examples
        ========

        >>> from sympy import Point
        >>> from sympy.geometry.hyperbola import Hyperbola
        >>> p1 = Point(0, 0)
        >>> h1 = Hyperbola(p1, 3, 1)
        >>> h1.semilatus_rectum
        1/3

        """

        return self.major * (self.eccentricity**2 - 1)

    def auxilary_circle(self):
        """Returns a Circle whose diameter is the major axis of the hyperbola.

        Examples
        ========

        >>> from sympy import Circle, Point, symbols
        >>> from sympy.geometry.hyperbola import Hyperbola
        >>> c = Point(1, 2)
        >>> h1 = Hyperbola(c, 8, 7)
        >>> h1.auxilary_circle()
        Circle(Point2D(1, 2), 8)
        >>> a, b = symbols('a b')
        >>> h2 = Hyperbola(c, a, b)
        >>> h2.auxilary_circle()
        Circle(Point2D(1, 2), a)
        """

        return Circle(self.center, self.hradius)

    @property
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

        >>> from sympy import Circle, Point, symbols
        >>> from sympy.geometry.hyperbola import Hyperbola
        >>> c = Point(3,8)
        >>> Hyperbola(c, 9, 7).director_circle
        Circle(Point2D(3, 8), 4*sqrt(2))
        >>> Hyperbola(c, 7, 7).director_circle
        Point2D(3, 8)

        """

        if self.hradius > self.vradius:
            return Circle(self.center, sqrt(self.hradius**2 - self.vradius**2))
        elif self.hradius == self.vradius:
            return self.center
        else:
            return None

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

        >>> from sympy import Point
        >>> from sympy.geometry.hyperbola import Hyperbola
        >>> h1 = Hyperbola(Point(0, 0), 3, 2)
        >>> h1.arbitrary_point()
        Point2D(3*sec(t), 2*tan(t))

        """

        t = _symbol(parameter, real=True)
        if t.name in (f.name for f in self.free_symbols):
            raise ValueError(filldedent('Symbol %s already appears in object '
                                        'and cannot be used as a parameter.' % t.name))
        return Point(self.center.x + self.hradius*sec(t),
                     self.center.y + self.vradius*tan(t))

    def equation(self, x='x', y='y', _slope=None):
        """
        Returns the equation of an hyperbola aligned with the x and y axes;
        when slope is given, the equation returned corresponds to an hyperbola
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

        equation : sympy expression

        See Also
        ========

        arbitrary_point : Returns parameterized point on hyperbola

        Examples
        ========

        >>> from sympy import Point, pi
        >>> from sympy.geometry.hyperbola import Hyperbola
        >>> from sympy.abc import x, y
        >>> h1 = Hyperbola(Point(1, 0), 3, 2)
        >>> eq1 = h1.equation(x, y); eq1
        -y**2/4 + (x/3 - 1/3)**2 - 1
        >>> eq2 = h1.equation(x, y, _slope=1); eq2
        -(-x + y + 1)**2/8 + (x + y - 1)**2/18 - 1

        A point on h1 satisfies eq1. Let's use one on the x-axis:

        >>> p1 = h1.center + Point(h1.major, 0)
        >>> assert eq1.subs(x, p1.x).subs(y, p1.y) == 0

        When rotated the same as the rotated ellipse, about the center
        point of the ellipse, it will satisfy the rotated ellipse's
        equation, too:

        >>> r1 = p1.rotate(pi/4, h1.center)
        >>> assert eq2.subs(x, r1.x).subs(y, r1.y) == 0

        """

        x = _symbol(x, real=True)
        y = _symbol(y, real=True)

        dx = x - self.center.x
        dy = y - self.center.y

        if _slope is not None:
            L = (dy - _slope * dx) ** 2
            l = (_slope * dy + dx) ** 2
            h = 1 + _slope ** 2
            b = h * self.major ** 2
            a = h * self.minor ** 2
            return l / b - L / a - 1
        else:
            t1 = (dx/self.hradius)**2
            t2 = (dy/self.vradius)**2
            return t1 - t2 - 1

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

        equation : sympy expression

        Examples
        ========

        >>> from sympy import Point
        >>> from sympy.geometry.hyperbola import Hyperbola
        >>> h1 = Hyperbola(Point(1, 0), 3, 2)
        >>> h1.evolute()
        -2**(2/3)*y**(2/3) + (3*x - 3)**(2/3) - 5**(2/3)
        """

        if len(self.args) != 3:
            raise NotImplementedError('Evolute of arbitrary Hyperbola is not supported.')
        x = _symbol(x, real=True)
        y = _symbol(y, real=True)
        t1 = (self.hradius*(x - self.center.x))**Rational(2, 3)
        t2 = (self.vradius*(y - self.center.y))**Rational(2, 3)
        return t1 - t2 - (self.hradius**2 - self.vradius**2)**Rational(2, 3)

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
        Currently supports intersections with Point, Line, Segment, Ray,
        Circle and Hyperbola types.

        See Also
        ========

        sympy.geometry.entity.GeometryEntity

        Examples
        ========

        >>> from sympy import Point, Line, sqrt
        >>> from sympy.geometry.hyperbola import Hyperbola
        >>> e = Hyperbola(Point(0, 0), 5, 7)
        >>> e.intersection(Point(0, 0))
        []
        >>> e.intersection(Point(5, 0))
        [Point2D(5, 0)]
        >>> e.intersection(Line(Point(0,0), Point(0, 1)))
        []
        >>> e.intersection(Line(Point(5,0), Point(5, 1)))
        [Point2D(5, 0)]
        >>> e.intersection(Line(Point(6,0), Point(6, 1)))
        [Point2D(6, -7*sqrt(11)/5), Point2D(6, 7*sqrt(11)/5)]
        >>> e = Hyperbola(Point(-1, 0), 4, 3)
        >>> e.intersection(Hyperbola(Point(1, 0), 4, 3))
        []
        >>> e.intersection(Hyperbola(Point(5, 0), 4, 3))
        []
        >>> e.intersection(Hyperbola(Point(100500, 0), 4, 3))
        [Point2D(100499/2, -3*sqrt(10100450937)/8), Point2D(100499/2, 3*sqrt(10100450937)/8)]
        >>> e.intersection(Hyperbola(Point(0, 0), 3, 4))
        [Point2D(3, 0)]
        >>> e.intersection(Hyperbola(Point(-1, 0), 3, 4))
        []
        """

        x = Dummy('x', real=True)
        y = Dummy('y', real=True)

        if isinstance(o, Point):
            if o in self:
                return [o]
            else:
                return []
        elif isinstance(o, (Segment2D, Ray2D)):
            hyperbola_eqaution = self.equation(x, y)
            result = solve([hyperbola_eqaution, Line(o.points[0], o.points[1]).equation(x, y)], [x, y])
            return list(ordered([Point(i) for i in result if i in o]))
        elif isinstance(o, Polygon):
            return o.intersection(self)
        elif isinstance(o, (Hyperbola, Line2D)):
            if o == self:
                return self
            else:
                hyperbola_equation = self.equation(x, y)
                return list(ordered([Point(i) for i in solve([hyperbola_equation, o.equation(x, y)], [x, y])]))
        elif isinstance(o, LinearEntity3D):
            raise TypeError('Entity must be two dimensional, not three dimensional')
        else:
            raise TypeError('Intersection not handled for %s' % func_name(o))
