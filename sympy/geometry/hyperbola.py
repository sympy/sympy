"""Hyperbolical geometrical entities.

Contains
* Hyperbola

"""

from sympy.core import S, C, sympify, pi, Dummy
from sympy.core.logic import fuzzy_bool
from sympy.core.numbers import oo, zoo
from sympy.simplify import simplify, trigsimp
from sympy.functions.elementary.miscellaneous import sqrt, Max, Min
from sympy.functions.elementary.complexes import im
from sympy.geometry.exceptions import GeometryError
from sympy.polys import Poly, PolynomialError
from sympy.solvers import solve
from sympy.utilities.lambdify import lambdify
from sympy.utilities.iterables import uniq
from sympy.utilities.misc import filldedent
from .entity import GeometryEntity
from .point import Point
from .line import LinearEntity, Line
from .ellipse import Ellipse, Circle
from .util import _symbol, idiff
from sympy.mpmath import findroot as nroot

class Hyperbola(GeometryEntity):
    """A Hyperbolical GeometryEntity.

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

    """
    def __new__(
        cls, center=None, hradius=None, vradius=None, eccentricity=None,
            **kwargs):
        hradius = sympify(hradius)
        vradius = sympify(vradius)

        eccentricity = sympify(eccentricity)

        if center is None:
            center = Point(0, 0)
        else:
            center = Point(center)

        if len(list(filter(None, (hradius, vradius, eccentricity)))) != 2:
            raise ValueError('Exactly two arguments of "hradius", '
                '"vradius", and "eccentricity" must not be None."')

        if eccentricity is not None:
            if hradius is None:
                hradius = vradius / sqrt(1 + eccentricity**2)
            elif vradius is None:
                vradius = hradius * sqrt(1 + eccentricity**2)

        return GeometryEntity.__new__(cls, center, hradius, vradius, **kwargs)

    @property
    def center(self):
        """The center of the Hyperbola.

        Returns
        =======

        center : Point

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.center
        Point(0, 0)

        """
        return self.args[0]

    @property
    def hradius(self):
        """The horizontal radius of the Hyperbola.

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
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.hradius
        3

        """
        return self.args[1]

    @property
    def vradius(self):
        """The vertical radius of the Hyperbola.

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
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.vradius
        1

        """
        return self.args[2]

    @property
    def minor(self):
        """Shorter axis of the Hyperbola (if it can be determined) else vradius.

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
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.minor
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
        if rv.func is Min:
            return self.vradius
        return rv

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
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.major
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
        if rv.func is Max:
            return self.hradius
        return rv

    @property
    def eccentricity(self):
        """The eccentricity of the hyperbola.

        Returns
        =======

        eccentricity : number

        Examples
        ========

        >>> from sympy import Point, Hyperbola, sqrt
        >>> p1 = Point(0, 0)
        >>> e1 = Hyperbola(p1, 3, sqrt(2))
        >>> e1.eccentricity
        sqrt(22)/2

        """
        return self.focus_distance / self.minor

    @property
    def focus_distance(self):
        """The focale distance of the Hyperbola.

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
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.focus_distance
        sqrt(10)

        """
        return Point.distance(self.center, self.foci[0])

    @property
    def foci(self):
        """The foci of the Hyperbola.

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
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.foci
        (Point(-sqrt(10), 0), Point(sqrt(10), 0))

        """
        c = self.center
        hr, vr = self.hradius, self.vradius
        if hr == vr:
            return (c, c)

        fd = sqrt(self.major**2 + self.minor**2)
        if hr == self.minor:
            # foci on the y-axis
            return (c + Point(0, -fd), c + Point(0, fd))
        elif hr == self.major:
            # foci on the x-axis
            return (c + Point(-fd, 0), c + Point(fd, 0))

    def rotate(self, angle=0, pt=None):
        """Rotate ``angle`` radians counterclockwise about Point ``pt``.

        Note: since the general hyperbola is not supported, the axes of
        the hyperbola will not be rotated. Only the center is rotated to
        a new position.

        Examples
        ========

        >>> from sympy import Hyperbola, pi
        >>> Hyperbola((1, 0), 2, 1).rotate(pi/2)
        Hyperbola(Point(0, 1), 2, 1)
        """
        return super(Hyperbola, self).rotate(angle, pt)

    def asymptotes(self):
        """ Returns a pair of asymtotes to the hyperbola.

        Notes
        =====

        Asymptotes of a hyperbola are the lines that pass through center of the
        hyperbola and touch the hyperbola at infinity.Basically they act as
        tangents to the two branches of hyperbola at infinity.

        Examples
        ========

        >>> from sympy import Hyperbola, Line, Point
        >>> a = Hyperbola(Point(0, 0), 3, 1)
        >>> a.asymptotes()
        [Line(Point(0, 0), Point(1, 1/3)), Line(Point(0, 0), Point(1, -1/3))]

        """
        slope = self.vradius / self.hradius
        p = self.center
        return [Line(p, slope = slope), Line(p, slope = -slope)]

    def equation(self, x='x', y='y'):
        """The equation of the hyperbola.

        Parameters
        ==========

        x : str, optional
            Label for the x-axis. Default value is 'x'.
        y : str, optional
            Label for the y-axis. Default value is 'y'.

        Returns
        =======

        equation : sympy expression

        See Also
        ========

        arbitrary_point : Returns parameterized point on hyperbola

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> e1 = Hyperbola(Point(1, 0), 3, 2)
        >>> e1.equation()
        -y**2/4 + (x/3 - 1/3)**2 - 1

        """
        x = _symbol(x)
        y = _symbol(y)
        t1 = ((x - self.center.x) / self.hradius)**2
        t2 = ((y - self.center.y) / self.vradius)**2
        return t1 - t2 - 1

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
        >>> e = Hyperbola((0, 0), 3, 2)
        >>> e.encloses_point((0, 0))
        True
        >>> e.encloses_point(e.arbitrary_point(t).subs(t, S.Half))
        False
        >>> e.encloses_point((4, 0))
        False

        """
        p = Point(p)
        if p in self:
            return False

        x = C.Dummy('x', real=True)
        y = C.Dummy('y', real=True)

        res = self.equation(x, y).subs({x: p.x, y: p.y})
        if res < 0:
            return True
        else:
            return False

    def arbitrary_point(self, parameter='t'):
        """A parameterized point on the Hyperbola.

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
        >>> e1 = Hyperbola(Point(0, 0), 3, 2)
        >>> e1.arbitrary_point()
        Point(3*sec(t), 2*tan(t))

        """
        t = _symbol(parameter)
        if t.name in (f.name for f in self.free_symbols):
            raise ValueError(filldedent('Symbol %s already appears in object '
                'and cannot be used as a parameter.' % t.name))
        return Point(self.center.x + self.hradius*C.sec(t),
                     self.center.y + self.vradius*C.tan(t))

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
        Circle and Ellipse types.

        See Also
        ========

        sympy.geometry.entity.GeometryEntity

        Examples
        ========

        >>> from sympy import Hyperbola, Point, Line, sqrt, Ellipse
        >>> e = Hyperbola(Point(0, 0), 5, 7)
        >>> e.intersection(Point(0, 0))
        []
        >>> e.intersection(Point(5, 0))
        [Point(5, 0)]
        >>> e.intersection(Line(Point(0,0), Point(0, 1)))
        []
        >>> e.intersection(Line(Point(5,0), Point(5, 1)))
        [Point(5, 0)]
        >>> e.intersection(Line(Point(6,0), Point(6, 1)))
        [Point(6, -7*sqrt(11)/5), Point(6, 7*sqrt(11)/5)]
        >>> e = Hyperbola(Point(-1, 0), 4, 3)
        >>> e.intersection(Ellipse(Point(1, 0), 4, 3))
        [Point(sqrt(15), -3*15**(1/4)*sqrt(2)/4), Point(sqrt(15), 3*15**(1/4)*sqrt(2)/4)]
        >>> e.intersection(Ellipse(Point(5, 0), 4, 3))
        [Point(2 + sqrt(7), -3*sqrt(6)*7**(1/4)/4), Point(2 + sqrt(7), 3*sqrt(6)*7**(1/4)/4)]
        >>> e.intersection(Hyperbola(Point(100500, 0), 4, 3))
        [Point(100499/2, -3*sqrt(10100450937)/8), Point(100499/2, 3*sqrt(10100450937)/8)]
        >>> e.intersection(Hyperbola(Point(0, 0), 3, 4))
        [Point(3, 0)]

        >>> e.intersection(Hyperbola(Point(-1, 0), 3, 4))
        []
        """
        if isinstance(o, Point):
            if o in self:
                return [o]
            else:
                return []

        elif isinstance(o, (LinearEntity, Circle, Ellipse, Hyperbola)):
            x = Dummy('x', real=True)
            y = Dummy('y', real=True)
            seq = self.equation(x, y)
            oeq = o.equation(x, y)
            result = solve([seq, oeq], [x, y])
            return [Point(*r) for r in list(uniq(result))]

        return self.intersection(o)

    def tangent_lines(self, p):
        """Tangent lines between `p` and the hyperbola.

        If `p` is on the hyperbola, returns the tangent line through point `p`.
        Otherwise, returns the tangent line(s) from `p` to the hyperbola, or
        None if no tangent line is possible (e.g., `p` is outside the hyperbola).

        Parameters
        ==========

        p : Point

        Returns
        =======

        tangent_lines : list with 1 or 2 Lines

        Raises
        ======

        NotImplementedError.

        See Also
        ========

        sympy.geometry.point.Point, sympy.geometry.line.Line

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> e1 = Hyperbola(Point(0, 0), 3, 2)
        >>> e1.tangent_lines(Point(3, 0))
        [Line(Point(3, 0), Point(3, -12))]

        """
        p = Point(p)
        if self.encloses_point(p) is False:
            if p in self:
                delta = self.center - p
                rise = (self.vradius ** 2)*delta.x
                run = -(self.hradius ** 2)*delta.y
                p2 = Point(simplify(p.x + run),
                       simplify(p.y + rise))
                return [Line(p, p2)]
            else:
                return []
        else:
            x, y = Dummy('x', real=True), Dummy('y', real=True)
            eq = self.equation(x, y)
            dydx = idiff(eq, y, x)
            slope = Line(p, Point(x, y)).slope
            tangent_points = solve([slope - dydx, eq], [x, y])

            # handle horizontal and vertical tangent lines
            if len(tangent_points) == 1:
                assert tangent_points[0][
                    0] == p.x or tangent_points[0][1] == p.y
                return [Line(p, p + Point(1, 0)), Line(p, p + Point(0, 1))]

            return [Line(p, tangent_points[0]), Line(p, tangent_points[1])]

    def is_tangent(self, o):
        """Is `o` tangent to the hyperbola?

        Parameters
        ==========

        o : GeometryEntity
            An Ellipse, LinearEntity or Polygon

        Raises
        ======

        NotImplementedError
            When the wrong type of argument is supplied.

        Returns
        =======

        is_tangent: boolean
            True if o is tangent to the hyperbola, False otherwise.

        See Also
        ========

        tangent_lines

        Examples
        ========

        >>> from sympy import Point, Hyperbola, Line
        >>> p0, p1, p2 = Point(0, 0), Point(3, 0), Point(3, 3)
        >>> e1 = Hyperbola(p0, 3, 2)
        >>> l1 = Line(p1, p2)
        >>> e1.is_tangent(l1)
        True

        """
        inter = None
        if isinstance(o, Hyperbola):
            inter = self.intersection(o)
            if isinstance(inter, Hyperbola):
                return False
            return (inter is not None and isinstance(inter[0], Point)
                    and len(inter) == 1)
        elif isinstance(o, LinearEntity):
            inter = self.intersection(o)
            if inter is not None and len(inter) == 1:
                return inter[0] in o
            else:
                return False
        elif isinstance(o, Polygon):
            c = 0
            for seg in o.sides:
                inter = self.intersection(seg)
                c += len([True for point in inter if point in seg])
            return c == 1
        else:
            raise NotImplementedError("Unknown argument type")

    def normal_lines(self, p, prec=None):
        """Normal lines between `p` and the hyperbola.

        Parameters
        ==========

        p : Point

        Returns
        =======

        normal_lines : list with 1, 2 or 4 Lines

        """
        p = Point(p)

        if True:
            rv = []
            if p.x == self.center.x:
                rv.append(Line(self.center, slope=oo))
            if p.y == self.center.y:
                rv.append(Line(self.center, slope=0))
            if rv:
                return rv

        x, y = Dummy('x', real=True), Dummy('y', real=True)
        eq = self.equation(x, y)
        dydx = idiff(eq, y, x)
        norm = -1/dydx
        slope = Line(p, (x, y)).slope
        seq = slope - norm
        points = []
        if prec is not None:
            yis = solve(seq, y)[0]
            xeq = eq.subs(y, yis).as_numer_denom()[0].expand()
            try:
                iv = list(zip(*Poly(xeq).intervals()))[0]
                xsol = [S(nroot(lambdify(x, xeq), i, solver="anderson"))
                    for i in iv]
                points = [Point(i, solve(eq.subs(x, i), y)[0]).n(prec)
                    for i in xsol]
            except PolynomialError:
                pass
        if not points:
            points = solve((seq, eq), (x, y))
            points = [Point(i).n(prec) if prec is not None else Point(i)
                      for i in points if all(j.n(2).is_real for j in i)]
        slopes = [norm.subs(zip((x, y), pt.args)) for pt in points]
        if prec is not None:
            slopes = [i.n(prec) if i not in (-oo, oo, zoo) else i
                for i in slopes]
        return [Line(pt, slope=s) for pt,s in zip(points, slopes)]

    def plot_interval(self, parameter='t'):
        """The plot interval for the default geometric plot of the Hyperbola.

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

        >>> from sympy import Point, Hyperbola
        >>> e1 = Hyperbola(Point(0, 0), 3, 2)
        >>> e1.plot_interval()
        [t, -pi, pi]

        """
        t = _symbol(parameter)
        return [t, -S.Pi, S.Pi]

    def random_point(self, seed=None):
        """A random point on the hyperbola.

        Returns
        =======

        point : Point

        See Also
        ========

        sympy.geometry.point.Point
        arbitrary_point : Returns parameterized point on hyperbola

        Notes
        -----

        A random point may not appear to be on the hyperbola, ie, `p in e` may
        return False. This is because the coordinates of the point will be
        floating point values, and when these values are substituted into the
        equation for the hyperbola the result may not be zero because of floating
        point rounding error.

        Examples
        ========

        >>> from sympy import Point, Hyperbola, Segment
        >>> e1 = Hyperbola(Point(0, 0), 3, 2)
        >>> e1.random_point() # gives some random point
        Point(...)

        """
        import random
        from sympy import sec, tan, Rational
        t = _symbol('t')
        x, y = self.arbitrary_point(t).args
        if seed is not None:
            rng = random.Random(seed)
        else:
            rng = random
        for i in range(10):  # should be enough?
            # simplify this now or else the Float will turn s into a Float
            c = 2*Rational(rng.random()) - 1
            s = sqrt(c**2 + 1)
            p1 = Point(x.subs(sec(t), s), y.subs(tan(t), c))
            if p1 in self:
                return p1
        raise GeometryError(
            'Having problems generating a point in the hyperbola.')

    def __contains__(self, o):
        if isinstance(o, Point):
            x = C.Dummy('x', real=True)
            y = C.Dummy('y', real=True)

            res = self.equation(x, y).subs({x: o.x, y: o.y})
            return trigsimp(simplify(res)) is S.Zero
        elif isinstance(o, Hyperbola):
            return self == o
        return False
