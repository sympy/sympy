"""Elliptical geometrical entities.

Contains
--------
Ellipse
Circle

"""

from sympy.core import S, C, sympify, symbol
from sympy.simplify import simplify, trigsimp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.complexes import im
from sympy.geometry.exceptions import GeometryError
from sympy.solvers import solve_poly_system, solve
from entity import GeometryEntity
from point import Point
from line import LinearEntity, Line

from sympy.abc import x

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
    Rotation is currently not supported since an ellipse is defined
    on horizontal/vertical radii

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

        if len(filter(None, (hradius, vradius, eccentricity))) != 2:
            raise ValueError, 'Exactly two arguments between "hradius", '\
                '"vradius", and "eccentricity" must be not None."'

        if eccentricity is not None:
            if hradius is None:
                hradius = vradius / sqrt(1 - eccentricity**2)
            elif vradius is None:
                vradius = hradius * sqrt(1 - eccentricity**2)
        else:
            if hradius is None and vradius is None:
                raise ValueError("At least two arguments between hradius, "
                    "vradius and eccentricity must not be none.")

        if center is None:
            center = Point(0, 0)

        if not isinstance(center, Point):
            raise TypeError("center must be a Point")

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
        12*Integral(((1 - 8*x**2/9)/(1 - x**2))**(1/2), (x, 0, 1))

        """
        if self.eccentricity == 1:
            return 2*pi*self.hradius
        else:
            return 4 * self.hradius * \
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
        return self.focus_distance / self.hradius

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
        3 - 2*2**(1/2)

        """
        return self.hradius * (1 - self.eccentricity)

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
        3 + 2*2**(1/2)

        """
        return self.hradius * (1 + self.eccentricity)

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
        The foci can only be calculated if the radii are numerical.

        Raises
        ------
        ValueError
            When the radii aren't numerical

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
        if self.hradius == self.vradius:
            return (c, c)

        hr, vr = self.hradius, self.vradius

        # calculate focus distance manually, since focus_distance calls this routine
        h = sqrt(abs(vr**2 - hr**2))
        if hr < vr:
            return (c + Point(0, -h), c + Point(0, h))
        else:
            return (c + Point(-h, 0), c + Point(h, 0))

    def tangent_line(self, p):
        """Tangent lines between `p` and the ellipse.

        If `p` is on the ellipse, returns the tangent line through point `p`.
        Otherwise, returns the tangent line(s) from `p` to the ellipse, or
        None if no tangent line is possible (e.g., `p` inside ellipse).

        Parameters
        ----------
        p : Point

        Returns
        -------
        tangent_line : Line

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
        >>> e1.tangent_line(Point(3, 0))
        Line(Point(3, 0), Point(3, -12))

        >>> # This will plot an ellipse together with a tangent line.
        >>> from sympy import Point, Ellipse, Plot
        >>> e = Ellipse(Point(0,0), 3, 2)
        >>> t = e.tangent_line(e.random_point()) # doctest: +SKIP
        >>> p = Plot() # doctest: +SKIP
        >>> p[0] = e # doctest: +SKIP
        >>> p[1] = t # doctest: +SKIP

        """
        if p in self:
            rise = (self.vradius ** 2)*(self.center[0] - p[0])
            run = (self.hradius ** 2)*(p[1] - self.center[1])
            p2 = Point(simplify(p[0] + run),
                       simplify(p[1] + rise))
            return Line(p, p2)
        else:
            # TODO If p is not on the ellipse, attempt to create the
            #      tangent(s) from point p to the ellipse..?
            raise NotImplementedError("Cannot find tangent lines when p is not on the ellipse")

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

    def arbitrary_point(self, parameter_name='t'):
        """A parametric point on the ellipse.

        Parameters
        ----------
        parameter_name : str, optional
            Default value is 't'.

        Returns
        -------
        arbitrary_point : Point

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
        t = C.Symbol(parameter_name, real=True)
        return Point(self.center[0] + self.hradius*C.cos(t),
                self.center[1] + self.vradius*C.sin(t))

    def plot_interval(self, parameter_name='t'):
        """The plot interval for the default geometric plot of the Ellipse.

        Parameters
        ----------
        parameter_name : str, optional
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
        t = C.Symbol(parameter_name, real=True)
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
        t = C.Symbol('t', real=True)
        p = self.arbitrary_point('t')
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
        -1 + (1/3 - x/3)**2 + y**2/4

        """
        if isinstance(x, basestring):   x = C.Symbol(x, real=True)
        if isinstance(y, basestring):   y = C.Symbol(y, real=True)
        t1 = ((x - self.center[0]) / self.hradius)**2
        t2 = ((y - self.center[1]) / self.vradius)**2
        return t1 + t2 - 1

    def _do_line_intersection(self, o):
        """The intersection of a LinearEntity and the ellipse.

        Private helper method for `intersection`.

        Notes
        -----
        Makes no regards to what the LinearEntity is because it assumes a
        Line. To ensure correct intersection results one must invoke
        `intersection` to remove bad results.

        """
        def dot(p1, p2):
            sum = 0
            for ind in xrange(0, len(p1)):
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
                result.append(lp[0] + (lp[1] - lp[0]) * t_a)
                result.append(lp[0] + (lp[1] - lp[0]) * t_b)
        return result

    def _do_circle_intersection(self, o):
        """The intersection of an Ellipse and a Circle.

        Private helper methode for `intersection`.

        """
        variables = self.equation().atoms(C.Symbol)
        if len(variables) > 2:
            return None
        if self.center == o.center:
            a, b, r = o.hradius, o.vradius, self.radius
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
        seq = self.equation()
        variables = self.equation().atoms(C.Symbol)
        if len(variables) > 2:
            return None
        x, y = variables
        oeq = o.equation(x=x, y=y)
        # until the following line works...
        # result = solve([seq, oeq], [x, y])
        # return [Point(*r) for r in result if im(r[0]).is_zero and im(r[1]).is_zero]
        # we do this:
        if self.center[0] == o.center[0] or self.center[1] == o.center[1]:
            result = solve_poly_system([seq, oeq], x, y)
            return [Point(*r) for r in result if im(r[0]).is_zero and im(r[1]).is_zero]

        raise NotImplementedError("Off-axis Ellipse intersection not supported.")

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
            result = self._do_line_intersection(o)
            if result is not None:
                for ind in xrange(len(result) - 1, -1, -1):
                    if result[ind] not in o:
                        del result[ind]
            return result
        elif isinstance(o, Circle):
            return self._do_circle_intersection(o)
        elif isinstance(o, Ellipse):
            if o == self:
                return self
            else:
                return self._do_ellipse_intersection(o)

        raise NotImplementedError()

    def distance_to_center(self, t):
        """The distance from an point on the ellipse to the center.

        Parameters
        ----------
        t : number
            Parameter in range [0, 2*pi)

        Examples
        --------
        >>> from sympy import Point, Ellipse, S
        >>> e1 = Ellipse(Point(0, 0), 3, 2)
        >>> # horizontal distance to center
        >>> e1.distance_to_center(0)
        3
        >>> # vertical distance to center
        >>> e1.distance_to_center(S.Pi/2)
        2

        """
        seg = Point.distance(self.center, self.arbitrary_point())

        return seg.subs('t', t)

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
    hradius
    vradius
    radius
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
        if len(args) == 3 and isinstance(args[0], Point):
            from polygon import Triangle
            t = Triangle(args[0], args[1], args[2])
            if t.area == 0:
                raise GeometryError("Cannot construct a circle from three collinear points")
            c = t.circumcenter
            r = t.circumradius
        elif len(args) == 2:
            # Assume (center, radius) pair
            c = args[0]
            r = sympify(args[1])

        if not (c is None or r is None):
            return GeometryEntity.__new__(cls, c, r, **kwargs)

        raise GeometryError("Circle.__new__ received unknown arguments")

    @property
    def hradius(self):
        """The horizontal radius of the circle.

        Returns
        -------
        hradius : number or sympy expression

        Examples
        --------
        >>> from sympy import Point, Circle
        >>> c1 = Circle(Point(3, 4), 6)
        >>> c1.hradius
        6

        """
        return self.__getitem__(1)

    @property
    def vradius(self):
        """The vertical radius of the circle.

        Returns
        -------
        vradius : number or sympy expression

        Examples
        --------
        >>> from sympy import Point, Circle
        >>> c1 = Circle(Point(3, 4), 6)
        >>> c1.vradius
        6

        """
        return self.__getitem__(1)

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
        -25 + x**2 + y**2

        """
        if isinstance(x, basestring):
            x = C.Symbol(x, real=True)
        if isinstance(y, basestring):
            y = C.Symbol(y, real=True)
        t1 = (x - self.center[0])**2
        t2 = (y - self.center[1])**2
        return t1 + t2 - self.hradius**2

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
