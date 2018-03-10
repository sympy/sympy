""" 2D geometrical entity constructed by boolean operations
    of various geometric elements.

"""

from __future__ import division, print_function
from sympy.core import S, pi, sympify
from sympy.core.compatibility import ordered
from sympy.core.symbol import _symbol
from sympy.core.numbers import Rational, oo
from sympy import symbols, simplify, solve
from sympy.geometry.entity import GeometryEntity, GeometrySet
from sympy.geometry.point import Point, Point2D
from sympy.geometry.line import Line, Line2D, Ray2D, Segment2D, LinearEntity3D
from sympy.geometry.ellipse import Ellipse, Circle
from sympy.geometry.parabola import Parabola
from sympy.functions.elementary.trigonometric import asin

from sympy import Abs, sqrt

class Shape(GeometrySet):
    """
    Parameters
    ==========

    conic : Ellipse, Circle, Parabola
    line : Line

    Attributes
    ==========

    area
    centroid
    second moment of area
    product moment of area

    Raises
    ======

    TypeError
        Wrong type of argument were put
        When The generated shape is not a closed figure
    NotImplementedError
        When `line` is neither horizontal nor vertical.

    Notes
    =====

    If the conic figure is Circle or Ellipse and if the line is vertical,
    then the shape will be segment that belongs to right side of the line.
    And if the line is horizontal, then the shape will be above the line.

    Examples
    ========

    >>> from sympy import Point, Line, Parabola, Ellipse, Circle, Shape
    >>> e = Ellipse((0, 0), 4, 2)
    >>> p = Parabola((2, 0), Line((-2, 0), (-2, 2)))
    >>> c = Circle((0, 0), 4)
    >>> l = Line((2, 0), (2, 2))
    >>> Shape(e, l)
    Shape(Ellipse(Point2D(0, 0), 4, 2), Line2D(Point2D(2, 0), Point2D(2, 2)))
    >>> Shape(p, l)
    Shape(Parabola(Point2D(2, 0), Line2D(Point2D(-2, 0), Point2D(-2, 2))), Line2D(Point2D(2, 0), Point2D(2, 2)))
    >>> Shape(c, l)
    Shape(Circle(Point2D(0, 0), 4), Line2D(Point2D(2, 0), Point2D(2, 2)))

    """

    def __new__(cls, conic, line):
        if (isinstance(conic, Circle) == 0 and isinstance(conic, Parabola) == 0 and isinstance(conic, Ellipse) == 0):
            raise TypeError('Wrong type of argument were put')

        if (line.slope != 0 and line.slope != S.Infinity):
            raise NotImplementedError('The line must be a horizontal'
                                      ' or vertical line')

        intersection = conic.intersection(line)

        if(len(intersection) < 2):
            if(isinstance(conic, Parabola)):
                raise TypeError('The shape is not a closed figure')
            else:
                raise TypeError('The conic-section is not cut by the line')

        return GeometryEntity.__new__(cls, conic, line)

    @property
    def name(self):
        if isinstance(self.conic, Circle):
            return 'Circle'
        if isinstance(self.conic, Parabola):
            return 'Parabola'
        if isinstance(self.conic, Ellipse):
            return 'Ellipse'

    @property
    def conic(self):
        return self.args[0]

    @property
    def line(self):
        return self.args[1]

    @property
    def area(self):
        """The area of the generated shape.

        See Also
        ========

        sympy.geometry.ellipse.Ellipse.area, sympy.geometry.ellipse.Circle.area

        Examples
        ========

        >>> from sympy import Point, Line, Parabola, Ellipse, Circle, Shape
        >>> e = Ellipse((0, 0), 4, 2)
        >>> p = Parabola((1, 0), Line((-3, 0), (-3, 2)))
        >>> c = Circle((0, 0), 4)
        >>> l = Line((0, 0), (0, 2))
        >>> Shape(e, l).area
        4*pi
        >>> Shape(p, l).area
        8*sqrt(2)/3
        >>> Shape(c, l).area
        8*pi

        """
        if(self.name == 'Parabola'):
            f_l = self.conic.focal_length

            if(self.line.slope == 0):
                l = Abs(self.line.p1[1] - self.conic.vertex[1])
                return ((8*sqrt(f_l)*sqrt(l)**3))/3

            if(self.line.slope == S.Infinity):
                l = Abs(self.line.p1[0] - self.conic.vertex[0])
                return ((8*sqrt(f_l)*sqrt(l)**3))/3
        else:
            a = self.conic.hradius
            b = self.conic.vradius

            if(self.line.slope == 0):
                l = (self.line.p1[1] - self.conic.center[1])
                return a*b*((S.Pi/2) - ((l/b)*sqrt(1 - (l/b)**2) + asin(l/b)))

            if(self.line.slope == S.Infinity):
                l = (self.line.p1[0] - self.conic.center[0])
                return a*b*((S.Pi/2) - ((l/a)*sqrt(1 - (l/a)**2) + asin(l/a)))

    @property
    def centroid(self):
        """The centroid of generated shape.

        Returns
        =======

        centroid : Point

        See Also
        ========

        sympy.geometry.point.Point, sympy.geometry.util.centroid

        Examples
        ========

        >>> from sympy import Point, Line, Parabola, Ellipse, Circle, Shape
        >>> e = Ellipse((0, 0), 4, 2)
        >>> p = Parabola((1, 0), Line((-3, 0), (-3, 2)))
        >>> c = Circle((0, 0), 4)
        >>> l = Line((0, 0), (0, 2))
        >>> Shape(e, l).centroid
        Point2D(16/(3*pi), 0)
        >>> Shape(p, l).centroid
        Point2D(-2/5, 0)
        >>> Shape(c, l).centroid
        Point2D(16/(3*pi), 0)

        """
        if(self.name == 'Parabola'):
            if(self.conic.directrix.slope == 0):
                return Point(self.conic.vertex[0], (2*self.conic.vertex[1] + 3*self.line.p1[1])/5)

            if(self.conic.directrix.slope == S.Infinity):
                return Point((2*self.conic.vertex[0] + 3*self.line.p1[0])/5, self.conic.vertex[1])

        else:
            a = self.conic.hradius
            b = self.conic.vradius
            if(self.line.slope == 0):
                l = (self.line.p1[1] - self.conic.center[1])
                y = (2*a*(b**2)/3)*((sqrt(1 - (l/b)**2))**3)
                return Point(self.conic.center[0], y/self.area)

            if(self.line.slope == S.Infinity):
                l = (self.line.p1[0] - self.conic.center[0])
                x = (2*b*(a**2)/3)*((sqrt(1 - (l/a)**2))**3)
                return Point( x/self.area, self.conic.center[1])

    def second_moment_of_area(self, point=None):
        """Returns the second moment and product moment of area of generated shape.

        Parameters
        ==========

        point : Point, two-tuple of sympifiable objects, or None(default=None)
            point is the point about which second moment of area is to be found.
            If "point=None" it will be calculated about the axis passing through the
            centroid of the generated shape.

        Returns
        =======

        I_xx, I_yy, I_xy : number or sympy expression
            I_xx, I_yy are second moment of area of generated shape.
            I_xy is product moment of area of generated shape.

        Examples
        ========

        >>> from sympy import Point, Line, Parabola, Ellipse, Circle, Shape
        >>> e = Ellipse((0, 0), 4, 2)
        >>> p = Parabola((1, 0), Line((-3, 0), (-3, 2)))
        >>> c = Circle((0, 0), 4)
        >>> l = Line((0, 0), (0, 2))
        >>> Shape(e, l).second_moment_of_area()
        (4*pi, -1024/(9*pi) + 16*pi, 0)
        >>> Shape(p, l).second_moment_of_area()
        (64*sqrt(2)/15, 32*sqrt(2)/175, 0)
        >>> Shape(c, l).second_moment_of_area()
        (32*pi, -2048/(9*pi) + 32*pi, 0)

        """
        if(self.name == 'Parabola'):
            a = self.conic.focal_length
            I_xy_v = 0;
            c = self.centroid
            Ar = self.area
            if(self.conic.directrix.slope == 0):
                I_xx_v = 32*sqrt((a**3))*sqrt(Abs((self.conic.vertex[1] - self.line.p1[1])**5))/15
                I_yy_v = 8*sqrt(a*(Abs((self.conic.vertex[1] - self.line.p1[1]))**7))/7

                I_xx_c = I_yy_v + Ar*((c[1] - self.conic.vertex[1])**2)
                I_yy_c = I_xx_v
                I_xy_c = I_xy_v
                if point is None:
                    return I_xx_c, I_yy_c, I_xy_c

                I_xx = I_xx_c + Ar*((c[1] - point[1])**2)
                I_yy = I_yy_c + Ar*((c[0] - point[0])**2)
                I_xy = I_xy_c + Ar*((point[0]-c[0])*(point[1]-c[1]))

                return I_xx, I_yy, I_xy

            if(self.conic.directrix.slope == S.Infinity):
                I_xx_v = 32*sqrt((a**3))*sqrt((Abs((self.conic.vertex[0] - self.line.p1[0]))**5))/15
                I_yy_v = 8*sqrt(a*(Abs((self.conic.vertex[0] - self.line.p1[0]))**7))/7

                I_xx_c = I_xx_v
                I_yy_c = I_yy_v - Ar*((c[0] - self.conic.vertex[0])**2)
                I_xy_c = I_xy_v
                if point is None:
                    return I_xx_c, I_yy_c, I_xy_c

                I_xx = I_xx_c + Ar*((c[1] - point[1])**2)
                I_yy = I_yy_c + Ar*((c[0] - point[0])**2)
                I_xy = I_xy_c + Ar*((point[0]-c[0])*(point[1]-c[1]))

                return I_xx, I_yy, I_xy

        else:
            a = self.conic.hradius
            b = self.conic.vradius
            I_xy_v = 0
            c = self.centroid
            Ar = self.area
            if(self.line.slope == 0):
                l = (self.line.p1[1] - self.conic.center[1])

                z = (l/b)*sqrt(1 - (l/b)**2)*(1 - 2*(l/b)**2)

                I_yy_v = b*(a**3)*(asin(sqrt(1 - (l/b)**2) + z))/4 - 2*l*(a*(sqrt(1 - (l/b)**2)))**3/3
                I_xx_v = a*(b**3)*((S.Pi/2) - asin(l/b) + z)/4

                I_xx_c = I_yy_v + Ar*((c[1] - self.conic.center[1])**2)
                I_yy_c = I_xx_v
                I_xy_c = I_xy_v
                if point is None:
                    return I_xx_c, I_yy_c, I_xy_c

                I_xx = I_xx_c + Ar*((c[1] - point[1])**2)
                I_yy = I_yy_c + Ar*((c[0] - point[0])**2)
                I_xy = I_xy_c + Ar*((point[0]-c[0])*(point[1]-c[1]))

                return I_xx, I_yy, I_xy

            if(self.line.slope == S.Infinity):
                l = (self.line.p1[0] - self.conic.center[0])

                z = (l/a)*sqrt(1 - (l/a)**2)*(1 - 2*(l/a)**2)

                I_yy_v = b*(a**3)*((S.Pi/2) - asin(l/a) + z)/4
                I_xx_v = a*(b**3)*(asin(sqrt(1 - (l/a)**2)) + z)/4 - 2*l*(b*(sqrt(1 - (l/a)**2)))**3/3

                I_xx_c = I_xx_v
                I_yy_c = I_yy_v - Ar*((c[0] - self.conic.center[0])**2)
                I_xy_c = I_xy_v

                if point is None:
                    return I_xx_c, I_yy_c, I_xy_c

                I_xx = I_xx_c + Ar*((c[1] - point[1])**2)
                I_yy = I_yy_c + Ar*((c[0] - point[0])**2)
                I_xy = I_xy_c + Ar*((point[0]-c[0])*(point[1]-c[1]))

                return I_xx, I_yy, I_xy
