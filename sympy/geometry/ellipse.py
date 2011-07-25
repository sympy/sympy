"""Elliptical geometrical entities.

Contains
--------
Ellipse
Circle

"""

from sympy.core import S, C, sympify, symbol, Dummy
from sympy.core.logic import fuzzy_bool
from sympy.simplify import simplify, trigsimp
from sympy.functions.elementary.miscellaneous import sqrt, Max, Min
from sympy.functions.elementary.complexes import im
from sympy.geometry.exceptions import GeometryError
from sympy.solvers import solve_poly_system, solve
from entity import GeometryEntity
from point import Point
from line import LinearEntity, Line
from util import _symbol, idiff

class Ellipse(GeometryEntity):
    """An elliptical GeometryEntity.

    Parameters
    ----------
    center : Point, optional
        Default value is Point(0, 0)
    hradius : number or sympy expression, optional
    vradius : number or sympy expression, optional
    eccentricity : number or sympy expression, optional
        Two of `hradius`, `vradius` and `eccentricity` must be supplied to
        create an Ellipse. The third is derived from the two supplied.

    Attributes
    ----------
    center
    hradius
    vradius
    area
    circumference
    eccentricity
    periapsis
    apoapsis
    focus_distance
    foci

    Raises
    ------
    GeometryError
        When `hradius`, `vradius` and `eccentricity` are incorrectly supplied as
        parameters.
    TypeError
        When `center` is not a Point.

    See Also
    --------
    Point

    Notes
    -----
    Constructed from a center and two radii, the first being the horizontal
    radius (along the x-axis) and the second being the vertical radius (along
    the y-axis).

    When symbolic value for hradius and vradius are used, any calculation that
    refers to the foci or the major or minor axis will assume that the ellipse
    has its major radius on the x-axis. If this is not true then a manual
    rotation is necessary.

    Examples
    --------
    >>> from sympy import Ellipse, Point, Rational
    >>> e1 = Ellipse(Point(0, 0), 5, 1)
    >>> e1.hradius, e1.vradius
    (5, 1)
    >>> e2 = Ellipse(Point(3, 1), hradius=3, eccentricity=Rational(4, 5))
    >>> e2
    Ellipse(Point(3, 1), 3, 9/5)

    Plotting example
    ----------------
    >>> from sympy import Circle, Plot, Segment
    >>> c1 = Circle(Point(0,0), 1)
    >>> Plot(c1)                                # doctest: +SKIP
    [0]: cos(t), sin(t), 'mode=parametric'
    >>> p = Plot()                              # doctest: +SKIP
    >>> p[0] = c1                               # doctest: +SKIP
    >>> radius = Segment(c1.center, c1.random_point())  # doctest: +SKIP
    >>> p[1] = radius                           # doctest: +SKIP
    >>> p                                       # doctest: +SKIP
    [0]: cos(t), sin(t), 'mode=parametric'
    [1]: t*cos(1.546086215036205357975518382),
    t*sin(1.546086215036205357975518382), 'mode=parametric'

    """

    def __new__(cls, center=None, hradius=None, vradius=None, eccentricity=None,
                **kwargs):
        hradius = sympify(hradius)
        vradius = sympify(vradius)

        eccentricity = sympify(eccentricity)

        if center is None:
            center = Point(0, 0)
        else:
            center = Point(center)

        if len(filter(None, (hradius, vradius, eccentricity))) != 2:
            raise ValueError('Exactly two arguments of "hradius", '\
                '"vradius", and "eccentricity" must not be None."')

        if eccentricity is not None:
            if hradius is None:
                hradius = vradius / sqrt(1 - eccentricity**2)
            elif vradius is None:
                vradius = hradius * sqrt(1 - eccentricity**2)

        if hradius == vradius:
            return Circle(center, hradius, **kwargs)

        return GeometryEntity.__new__(cls, center, hradius, vradius, **kwargs)

    @property
    def center(self):
        """The center of the ellipse.

        Returns
        -------
        center : number

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.center
        Point(0, 0)

        """
        return self.__getitem__(0)

    @property
    def hradius(self):
        """The horizontal radius of the ellipse.

        Returns
        -------
        hradius : number

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.hradius
        3

        """
        return self.__getitem__(1)

    @property
    def vradius(self):
        """The vertical radius of the ellipse.

        Returns
        -------
        vradius : number

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.vradius
        1

        """
        return self.__getitem__(2)

    @property
    def minor(self):
        """Shorter axis of the ellipse (if it can be determined) else vradius.

        Returns
        -------
        minor : number or expression

        Examples
        --------
        >>> from sympy import Point, Ellipse, Symbol
        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.minor
        1

        >>> a = Symbol('a')
        >>> b = Symbol('b')
        >>> Ellipse(p1, a, b).minor
        b
        >>> Ellipse(p1, b, a).minor
        a

        >>> m = Symbol('m')
        >>> M = m + 1
        >>> Ellipse(p1, m, M).minor
        m

        """
        rv = Min(*self[1:3])
        if rv.func is Min:
            return self.vradius
        return rv

    @property
    def major(self):
        """Longer axis of the ellipse (if it can be determined) else hradius.

        Returns
        -------
        major : number or expression

        Examples
        --------
        >>> from sympy import Point, Ellipse, Symbol
        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.major
        3

        >>> a = Symbol('a')
        >>> b = Symbol('b')
        >>> Ellipse(p1, a, b).major
        a
        >>> Ellipse(p1, b, a).major
        b

        >>> m = Symbol('m')
        >>> M = m + 1
        >>> Ellipse(p1, m, M).major
        m + 1

        """
        rv = Max(*self[1:3])
        if rv.func is Max:
            return self.hradius
        return rv


    @property
    def area(self):
        """The area of the ellipse.

        Returns
        -------
        area : number

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.area
        3*pi

        """
        return simplify(S.Pi * self.hradius * self.vradius)

    @property
    def circumference(self):
        """The circumference of the ellipse.

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.circumference
        12*Integral(((-8*_x**2/9 + 1)/(-_x**2 + 1))**(1/2), (_x, 0, 1))

        """
        if self.eccentricity == 1:
            return 2*pi*self.hradius
        else:
            x = C.Dummy('x', real=True)
            return 4*self.major*\
                   C.Integral(sqrt((1 - (self.eccentricity*x)**2)/(1 - x**2)),
                              (x, 0, 1))

    @property
    def eccentricity(self):
        """The eccentricity of the ellipse.

        Returns
        -------
        eccentricity : number

        Examples
        --------
        >>> from sympy import Point, Ellipse, sqrt
        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, sqrt(2))
        >>> e1.eccentricity
        7**(1/2)/3

        """
        return self.focus_distance / self.major

    @property
    def periapsis(self):
        """The periapsis of the ellipse.

        The shortest distance between the focus and the contour.

        Returns
        -------
        periapsis : number

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.periapsis
        -2*2**(1/2) + 3

        """
        return self.major * (1 - self.eccentricity)

    @property
    def apoapsis(self):
        """The periapsis of the ellipse.

        The greatest distance between the focus and the contour.

        Returns
        -------
        apoapsis : number

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.apoapsis
        2*2**(1/2) + 3

        """
        return self.major * (1 + self.eccentricity)

    @property
    def focus_distance(self):
        """The focale distance of the ellipse.

        The distance between the center and one focus.

        Returns
        -------
        focus_distance : number

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.focus_distance
        2*2**(1/2)

        """
        return Point.distance(self.center, self.foci[0])

    @property
    def foci(self):
        """The foci of the ellipse.

        Notes
        -----
        The foci can only be calculated if the major/minor axes are known.

        Raises
        ------
        ValueError
            When the major and minor axis cannot be determined.

        See Also
        --------
        Point

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> p1 = Point(0, 0)
        >>> e1 = Ellipse(p1, 3, 1)
        >>> e1.foci
        (Point(-2*2**(1/2), 0), Point(2*2**(1/2), 0))

        """
        c = self.center
        hr, vr = self.hradius, self.vradius
        if hr == vr:
            return (c, c)


        # calculate focus distance manually, since focus_distance calls this routine
        fd = sqrt(self.major**2 - self.minor**2)
        if hr == self.minor:
            # foci on the y-axis
            return (c + Point(0, -fd), c + Point(0, fd))
        elif hr == self.major:
            # foci on the x-axis
            return (c + Point(-fd, 0), c + Point(fd, 0))

    def encloses_point(self, p):
        """
        Return True if p is enclosed by (is inside of) self.

        Notes
        -----
        Being on the border of self is considered False.

        Parameters
        ----------
        p : Point

        Returns
        -------
        encloses_point : True, False or None

        Examples
        --------
        >>> from sympy import Ellipse, S
        >>> from sympy.abc import t
        >>> e = Ellipse((0, 0), 3, 2)
        >>> e.encloses_point((0, 0))
        True
        >>> e.encloses_point(e.arbitrary_point(t).subs(t, S.Half))
        False
        >>> e.encloses_point((4, 0))
        False

        """
        if p in self:
            return False

        if len(self.foci) == 2:
            # if the combined distance from the foci to p (h1 + h2) is less
            # than the combined distance from the foci to the minor axis
            # (which is the same as the major axis length) then p is inside
            # the ellipse
            h1, h2 = [f.distance(p) for f in self.foci]
            test = 2*self.major - (h1 + h2)
        else:
            test = self.radius - self.center.distance(p)

        return fuzzy_bool(test.is_positive)

    def tangent_lines(self, p):
        """Tangent lines between `p` and the ellipse.

        If `p` is on the ellipse, returns the tangent line through point `p`.
        Otherwise, returns the tangent line(s) from `p` to the ellipse, or
        None if no tangent line is possible (e.g., `p` inside ellipse).

        Parameters
        ----------
        p : Point

        Returns
        -------
        tangent_lines : list with 1 or 2 Lines

        Raises
        ------
        NotImplementedError
            Can only find tangent lines for a point, `p`, on the ellipse.

        See Also
        --------
        Point
        Line

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> e1 = Ellipse(Point(0, 0), 3, 2)
        >>> e1.tangent_lines(Point(3, 0))
        [Line(Point(3, 0), Point(3, -12))]

        >>> # This will plot an ellipse together with a tangent line.
        >>> from sympy import Point, Ellipse, Plot
        >>> e = Ellipse(Point(0,0), 3, 2)
        >>> t = e.tangent_lines(e.random_point()) # doctest: +SKIP
        >>> p = Plot() # doctest: +SKIP
        >>> p[0] = e # doctest: +SKIP
        >>> p[1] = t # doctest: +SKIP

        """
        from sympy import solve

        if self.encloses_point(p):
            return []

        if p in self:
            rise = (self.vradius ** 2)*(self.center[0] - p[0])
            run = (self.hradius ** 2)*(p[1] - self.center[1])
            p2 = Point(simplify(p[0] + run),
                       simplify(p[1] + rise))
            return [Line(p, p2)]
        else:
            if len(self.foci) == 2:
                f1, f2 = self.foci
                maj = self.hradius
                test = (2*maj -
                        Point.distance(f1, p) -
                        Point.distance(f2, p))
            else:
                test = self.radius - Point.distance(self.center, p)
            if test.is_number and test.is_positive:
                return []
            # else p is outside the ellipse or we can't tell. In case of the
            # latter, the solutions returned will only be valid if
            # the point is not inside the ellipse; if it is, nan will result.
            x, y = Dummy('x'), Dummy('y')
            eq = self.equation(x, y)
            dydx = idiff(eq, y, x)
            slope = Line(p, Point(x, y)).slope
            tangent_points = solve([w.as_numer_denom()[0] for w in [slope - dydx, eq]], [x, y])

            # handle horizontal and vertical tangent lines
            if len(tangent_points) == 1:
                assert tangent_points[0][0] == p[0] or tangent_points[0][1] == p[1]
                return [Line(p, Point(p[0]+1, p[1])), Line(p, Point(p[0], p[1]+1))]

            # others
            return [Line(p, tangent_points[0]), Line(p, tangent_points[1])]

    def is_tangent(self, o):
        """Is `o` tangent to the ellipse?

        Parameters
        ----------
        o : GeometryEntity
            An Ellipse, LinearEntity or Polygon

        Raises
        ------
        NotImplementedError
            When the wrong type of argument is supplied.

        Returns
        -------
        is_tangent: boolean
            True if o is tangent to the ellipse, False otherwise.

        Examples
        --------
        >>> from sympy import Point, Ellipse, Line
        >>> p0, p1, p2 = Point(0, 0), Point(3, 0), Point(3, 3)
        >>> e1 = Ellipse(p0, 3, 2)
        >>> l1 = Line(p1, p2)
        >>> e1.is_tangent(l1)
        True

        """
        inter = None
        if isinstance(o, Ellipse):
            inter = self.intersection(o)
            return (inter is not None and isinstance(inter[0], Point)
                    and len(inter) == 1)
        elif isinstance(o, LinearEntity):
            inter = self._do_line_intersection(o)
            if inter is not None and len(inter) == 1:
                return inter[0] in o
            else:
                return False
        elif isinstance(o, Polygon):
            c = 0
            for seg in o.sides:
                inter = self._do_line_intersection(seg)
                c += len([True for point in inter if point in seg])
            return c == 1
        else:
            raise NotImplementedError("Unknown argument type")

    def arbitrary_point(self, parameter='t'):
        """A parameterized point on the ellipse.

        Parameters
        ----------
        parameter : str, optional
            Default value is 't'.

        Returns
        -------
        arbitrary_point : Point

        Raises
        ------
        ValueError
            When `parameter` already appears in the functions.

        See Also
        --------
        Point

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> e1 = Ellipse(Point(0, 0), 3, 2)
        >>> e1.arbitrary_point()
        Point(3*cos(t), 2*sin(t))

        """
        t = _symbol(parameter)
        if t.name in (f.name for f in self.free_symbols):
            raise ValueError('Symbol %s already appears in object and cannot be used as a parameter.' % t.name)
        return Point(self.center[0] + self.hradius*C.cos(t),
                self.center[1] + self.vradius*C.sin(t))

    def plot_interval(self, parameter='t'):
        """The plot interval for the default geometric plot of the Ellipse.

        Parameters
        ----------
        parameter : str, optional
            Default value is 't'.

        Returns
        -------
        plot_interval : list
            [parameter, lower_bound, upper_bound]

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> e1 = Ellipse(Point(0, 0), 3, 2)
        >>> e1.plot_interval()
        [t, -pi, pi]

        """
        t = _symbol(parameter)
        return [t, -S.Pi, S.Pi]

    def random_point(self):
        """A random point on the ellipse.

        Returns
        -------
        point : Point

        See Also
        --------
        Point

        Note
        ----
        A random point may not appear to be on the ellipse, ie, `p in e` may
        return False. This is because the coordinates of the point will be
        floating point values, and when these values are substituted into the
        equation for the ellipse the result may not be zero because of floating
        point rounding error.

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> e1 = Ellipse(Point(0, 0), 3, 2)
        >>> p1 = e1.random_point()
        >>> # a random point may not appear to be on the ellipse because of
        >>> # floating point rounding error
        >>> p1 in e1 # doctest: +SKIP
        True
        >>> p1 # doctest +ELLIPSIS
        Point(...)

        """
        from random import random
        t = _symbol('t')
        p = self.arbitrary_point(t)
        # get a random value in [-pi, pi)
        subs_val = float(S.Pi)*(2*random() - 1)
        return Point(p[0].subs(t, subs_val), p[1].subs(t, subs_val))

    def equation(self, x='x', y='y'):
        """The equation of the ellipse.

        Parameters
        ----------
        x : str, optional
            Label for the x-axis. Default value is 'x'.
        y : str, optional
            Label for the y-axis. Default value is 'y'.

        Returns
        -------
        equation : sympy expression

        Examples
        --------
        >>> from sympy import Point, Ellipse
        >>> e1 = Ellipse(Point(1, 0), 3, 2)
        >>> e1.equation()
        y**2/4 + (x/3 - 1/3)**2 - 1

        """
        x = _symbol(x)
        y = _symbol(y)
        t1 = ((x - self.center[0]) / self.hradius)**2
        t2 = ((y - self.center[1]) / self.vradius)**2
        return t1 + t2 - 1

    def _do_line_intersection(self, o):
        """
        Find the intersection of a LinearEntity and the ellipse.

        All LinearEntities are treated as a line and filtered at
        the end to see that they lie in o.
        """
        def dot(p1, p2):
            sum = 0
            for ind in xrange(len(p1.args)):
                sum += p1[ind] * p2[ind]
            return simplify(sum)

        hr_sq = self.hradius ** 2
        vr_sq = self.vradius ** 2
        lp = o.points

        ldir = lp[1] - lp[0]
        diff = lp[0] - self.center
        mdir = (ldir[0] / hr_sq, ldir[1] / vr_sq)
        mdiff = (diff[0] / hr_sq, diff[1] / vr_sq)

        a = dot(ldir, mdir)
        b = dot(ldir, mdiff)
        c = dot(diff, mdiff) - 1
        det = simplify(b*b - a*c);

        result = []
        if det == 0:
            t = -b / a
            result.append(lp[0] + (lp[1] - lp[0]) * t)
        else:
            is_good = True
            try:
                is_good = (det > 0)
            except NotImplementedError: #symbolic, allow
                is_good = True

            if is_good:
                root = sqrt(det)
                t_a = (-b - root) / a
                t_b = (-b + root) / a
                result.append( lp[0] + (lp[1] - lp[0]) * t_a )
                result.append( lp[0] + (lp[1] - lp[0]) * t_b )

        return [r for r in result if r in o]

    def _do_circle_intersection(self, o):
        """The intersection of an Ellipse and a Circle.

        Private helper method for `intersection`.

        """
        variables = self.equation().atoms(C.Symbol)
        if len(variables) > 2:
            return None
        if self.center == o.center:
            a, b, r = o.major, o.minor, self.radius
            x = a*sqrt(simplify((r**2 - b**2)/(a**2 - b**2)))
            y = b*sqrt(simplify((a**2 - r**2)/(a**2 - b**2)))
            return list(set([Point(x, y), Point(x, -y), Point(-x, y),
                             Point(-x, -y)]))
        else:
            x, y = variables
            xx = solve(self.equation(), x)
            intersect = []
            for xi in xx:
                yy = solve(o.equation().subs(x, xi), y)
                for yi in yy:
                    intersect.append(Point(xi, yi))
            return list(set(intersect))

    def _do_ellipse_intersection(self, o):
        """The intersection of two ellipses.

        Private helper method for `intersection`.

        """
        x = Dummy('x')
        y = Dummy('y')
        seq = self.equation(x, y)
        oeq = o.equation(x, y)
        result = solve([seq, oeq], [x, y])
        return [Point(*r) for r in result if im(r[0]).is_zero and im(r[1]).is_zero]

    def intersection(self, o):
        """The intersection of this ellipse and another geometrical entity
        `o`.

        Parameters
        ----------
        o : GeometryEntity

        Returns
        -------
        intersection : list of GeometryEntities

        Notes
        -----
        Currently supports intersections with Point, Line, Segment, Ray,
        Circle and Ellipse types.

        Example
        -------
        >>> from sympy import Ellipse, Point, Line, sqrt
        >>> e = Ellipse(Point(0, 0), 5, 7)
        >>> e.intersection(Point(0, 0))
        []
        >>> e.intersection(Point(5, 0))
        [Point(5, 0)]
        >>> e.intersection(Line(Point(0,0), Point(0, 1)))
        [Point(0, -7), Point(0, 7)]
        >>> e.intersection(Line(Point(5,0), Point(5, 1)))
        [Point(5, 0)]
        >>> e.intersection(Line(Point(6,0), Point(6, 1)))
        []
        >>> e = Ellipse(Point(-1, 0), 4, 3)
        >>> e.intersection(Ellipse(Point(1, 0), 4, 3))
        [Point(0, -3*15**(1/2)/4), Point(0, 3*15**(1/2)/4)]
        >>> e.intersection(Ellipse(Point(5, 0), 4, 3))
        [Point(2, -3*7**(1/2)/4), Point(2, 3*7**(1/2)/4)]
        >>> e.intersection(Ellipse(Point(100500, 0), 4, 3))
        []
        >>> e.intersection(Ellipse(Point(0, 0), 3, 4))
        [Point(-363/175, -48*111**(1/2)/175), Point(-363/175, 48*111**(1/2)/175),
        Point(3, 0)]
        >>> e.intersection(Ellipse(Point(-1, 0), 3, 4))
        [Point(-17/5, -12/5), Point(-17/5, 12/5), Point(7/5, -12/5),
        Point(7/5, 12/5)]
        """
        if isinstance(o, Point):
            if o in self:
                return [o]
            else:
                return []

        elif isinstance(o, LinearEntity):
            # LinearEntity may be a ray/segment, so check the points
            # of intersection for coincidence first
            return self._do_line_intersection(o)

        elif isinstance(o, Circle):
            return self._do_circle_intersection(o)

        elif isinstance(o, Ellipse):
            if o == self:
                return self
            else:
                return self._do_ellipse_intersection(o)

        return o.intersection(self)

    def __eq__(self, o):
        """Is the other GeometryEntity the same as this ellipse?"""
        return (self.center == o.center and self.hradius == o.hradius
                and self.vradius == o.vradius)

    def __hash__(self):
        return super(Ellipse, self).__hash__()

    def __contains__(self, o):
        if isinstance(o, Point):
            x = C.Dummy('x', real=True)
            y = C.Dummy('y', real=True)

            res = self.equation(x, y).subs({x: o[0], y: o[1]})
            return trigsimp(simplify(res)) is S.Zero
        elif isinstance(o, Ellipse):
            return self == o
        return False


class Circle(Ellipse):
    """A circle in space.

    Constructed simply from a center and a radius, or from three
    non-collinear points.

    Parameters
    ----------
    center : Point
    radius : number or sympy expression
    points : sequence of three Points

    Attributes
    ----------
    radius (synonymous with hradius, vradius, major and minor)
    circumference
    equation

    Raises
    ------
    GeometryError
        When trying to construct circle from three collinear points.
        When trying to construct circle from incorrect parameters.

    See Also
    --------
    Point

    Examples
    --------
    >>> from sympy.geometry import Point, Circle
    >>> # a circle constructed from a center and radius
    >>> c1 = Circle(Point(0, 0), 5)
    >>> c1.hradius, c1.vradius, c1.radius
    (5, 5, 5)

    >>> # a circle costructed from three points
    >>> c2 = Circle(Point(0, 0), Point(1, 1), Point(1, 0))
    >>> c2.hradius, c2.vradius, c2.radius, c2.center
    (2**(1/2)/2, 2**(1/2)/2, 2**(1/2)/2, Point(1/2, 1/2))

    """
    def __new__(cls, *args, **kwargs):
        c, r = None, None
        if len(args) == 3:
            args = [Point(a) for a in args]
            if Point.is_collinear(*args):
                raise GeometryError("Cannot construct a circle from three collinear points")
            from polygon import Triangle
            t = Triangle(*args)
            c = t.circumcenter
            r = t.circumradius
        elif len(args) == 2:
            # Assume (center, radius) pair
            c = Point(args[0])
            r = sympify(args[1])

        if not (c is None or r is None):
            return GeometryEntity.__new__(cls, c, r, **kwargs)

        raise GeometryError("Circle.__new__ received unknown arguments")

    @property
    def radius(self):
        """The radius of the circle.

        Returns
        -------
        radius : number or sympy expression

        Examples
        --------
        >>> from sympy import Point, Circle
        >>> c1 = Circle(Point(3, 4), 6)
        >>> c1.radius
        6

        """
        return self.__getitem__(1)

    @property
    def major(self):
        return self.radius

    @property
    def minor(self):
        return self.radius

    @property
    def hradius(self):
        return self.radius

    @property
    def vradius(self):
        return self.radius

    @property
    def circumference(self):
        """The circumference of the circle.

        Returns
        -------
        circumference : number or sympy expression

        Examples
        --------
        >>> from sympy import Point, Circle
        >>> c1 = Circle(Point(3, 4), 6)
        >>> c1.circumference
        12*pi

        """
        return 2 * S.Pi * self.radius

    def equation(self, x='x', y='y'):
        """The equation of the circle.

        Parameters
        ----------
        x : str or Symbol, optional
            Default value is 'x'.
        y : str or Symbol, optional
            Default value is 'y'.

        Returns
        -------
        equation : sympy expression

        Examples
        --------
        >>> from sympy import Point, Circle
        >>> c1 = Circle(Point(0, 0), 5)
        >>> c1.equation()
        x**2 + y**2 - 25

        """
        x = _symbol(x)
        y = _symbol(y)
        t1 = (x - self.center[0])**2
        t2 = (y - self.center[1])**2
        return t1 + t2 - self.major**2

    def intersection(self, o):
        """The intersection of this circle with another geometrical entity.

        Parameters
        ----------
        o : GeometryEntity

        Returns
        -------
        intersection : list of GeometryEntities

        Examples
        --------
        >>> from sympy import Point, Circle, Line, Ray
        >>> p1, p2, p3 = Point(0, 0), Point(5, 5), Point(6, 0)
        >>> p4 = Point(5, 0)
        >>> c1 = Circle(p1, 5)
        >>> c1.intersection(p2)
        []
        >>> c1.intersection(p4)
        [Point(5, 0)]
        >>> c1.intersection(Ray(p1, p2))
        [Point(5*2**(1/2)/2, 5*2**(1/2)/2)]
        >>> c1.intersection(Line(p2, p3))
        []

        """
        if isinstance(o, Circle):
            if o.center == self.center:
                if o.radius == self.radius:
                    return o
                return []
            dx,dy = o.center - self.center
            d = sqrt(simplify(dy**2 + dx**2))
            R = o.radius + self.radius
            if d > R or d < abs(self.radius - o.radius):
                return []

            a = simplify((self.radius**2 - o.radius**2 + d**2) / (2*d))

            x2 = self.center[0] + (dx * a/d)
            y2 = self.center[1] + (dy * a/d)

            h = sqrt(simplify(self.radius**2 - a**2))
            rx = -dy * (h/d)
            ry = dx * (h/d)

            xi_1 = simplify(x2 + rx)
            xi_2 = simplify(x2 - rx)
            yi_1 = simplify(y2 + ry)
            yi_2 = simplify(y2 - ry)

            ret = [Point(xi_1, yi_1)]
            if xi_1 != xi_2 or yi_1 != yi_2:
                ret.append(Point(xi_2, yi_2))
            return ret

        return Ellipse.intersection(self, o)

from polygon import Polygon
