"""Hyperbolical geometrical entities.

Contains
* Hyperbola
* Rectangular Hyperbola

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
from .util import _symbol, idiff
from sympy.mpmath import findroot as nroot

class Hyperbola(GeometryEntity):
    """A Hyoerbolical GeometryEntity.

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
    periapsis
    apoapsis
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

        center : number

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
        sqrt(11)/3

        """
        return self.focus_distance / self.major

        @property
    def periapsis(self):
        """The periapsis of the Hyperbola.

        The shortest distance between the focus and the contour.

        Returns
        =======

        periapsis : number

        See Also
        ========

        apoapsis : Returns greatest distance between focus and contour

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.periapsis
        -2*sqrt(2) + 3

        """
        return self.major * (1 - self.eccentricity)

    @property
    def apoapsis(self):
        """The apoapsis of the Hyperbola.

        The greatest distance between the focus and the contour.

        Returns
        =======

        apoapsis : number

        See Also
        ========

        periapsis : Returns shortest distance between foci and contour

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.apoapsis
        2*sqrt(2) + 3

        """
        return self.major * (1 + self.eccentricity)

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
        2*sqrt(2)

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
        (Point(-2*sqrt(2), 0), Point(2*sqrt(2), 0))

        """
        c = self.center
        hr, vr = self.hradius, self.vradius
        if hr == vr:
            return (c, c)

        fd = sqrt(self.major**2 - self.minor**2)
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

        >>> from sympy import hyperbola, pi
        >>> Hyperbola((1, 0), 2, 1).rotate(pi/2)
        Hyperbola(Point(0, 1), 2, 1)
        """
        return super(Hyperbola, self).rotate(angle, pt)

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

        res = self.equation(x, y).subs({x: o.x, y: o.y})
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

    def __contains__(self, o):
        if isinstance(o, Point):
            x = C.Dummy('x', real=True)
            y = C.Dummy('y', real=True)

            res = self.equation(x, y).subs({x: o.x, y: o.y})
            return trigsimp(simplify(res)) is S.Zero
        elif isinstance(o, Hyperbola):
            return self == o
        return False

class RecHyperbola(Hyperbola):
    #TODO
    pass
