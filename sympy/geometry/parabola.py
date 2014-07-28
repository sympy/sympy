"""Parabolic geometrical entities.

Contains
* Parabola

"""

from __future__ import print_function, division

from sympy.core import S, C, sympify, pi, Dummy, Symbol, symbols
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

class Parabola(GeometryEntity):
    """A Parabolic GeometryEntity.

    Parameters
    ==========

    center : Point, optional
        Default value is Point(0, 0)
    a : number or SymPy expression, optional
    parallel_axis: string or SymPy expression, optional

    Notes
    =====

    The general form of the parabola which is parallel to y axis is y**2 = 4*a*x
    In the above case center is (0,0), a is a constant and parallel_axis is 'y'

    Examples
    ========

    >>> from sympy import Parabola, Point
    >>> a = Parabola(Point(0, 0), 2)
    Parabola(Point(0, 0), 2, y) # Default axis is 'y'
    >>> a = Parabola(Point(2, 3), axis='x')
    Parabola(Point(2, 3), 1, y) # Default value of a is 1

    """
    def __new__(cls, center=None, a=None, parallel_axis=None, **kwargs):
        if center is None:
            center = Point(0, 0)
        else:
            center = Point(center)
        if a is None:
            a = S.One
        else:
            a = sympify(a)
        if parallel_axis is None:
            parallel_axis = Symbol('y')
        else:
            if sympify(parallel_axis) in (Symbol('x'), Symbol('y')):
                parallel_axis = Symbol(parallel_axis)
            else:
                raise NotImplementedError('The axis should be parallel to x or y axis')
        return GeometryEntity.__new__(cls, center, a, parallel_axis, **kwargs)

    @property
    def center(self):
        """Returns the center of the given parabola.

        Examples
        ========

        >>> from sympy import Parabola, Point
        >>> a = Parabola(Point(1, 2), 2, 'x')
        >>> a.center
        Point(1, 2)
        >>> a = Parabola(a=2, parallel_axis='x')
        Parabola(Point(0, 0), 2, x)

        """
        return self.args[0]

    @property
    def constant(self):
        """In the general equation of the parabola y**2 = 4*a*x, a is the constant
        By varying a we get a family of parabolas passing through the given center
        and the given axis.

        Examples
        ========

        >>> from sympy import Parabola, Point
        >>> a = Parabola(Point(0, 0), 2)
        >>> a.constant
        2

        """
        return self.args[1]

    @property
    def axis(self):
        """The axis of the parabola passing through the center and parallel to
        the given 'x' or 'y' axis

        Examples
        ========

        >>> from sympy import Parabola, Point
        >>> a = Parabola(Point(2, 3), 4)
        >>> a.axis
        Line(Point(2, 3), Point(2, 4))
        >>> a = Parabola(Point(2, 3), 4, 'x')
        >>> a.axis
        Line(Point(2, 3), Point(3, 3))

        """
        x = Symbol('x')
        if self.args[2] == x:
            return Line(self.center, self.center + Point(1, 0))
        else:
            return Line(self.center, self.center + Point(0, 1))

    @property
    def focus(self):
        """The focus of the given parabola.

        Examples
        ========

        >>> from sympy import Parabola, Point
        >>> a = Parabola(Point(1, 4), 3)
        >>> a.focus
        Point(1, 7)

        """
        x = Symbol('x')
        if self.args[2] == x:
            return self.center + Point(self.constant, 0)
        else:
            return self.center + Point(0, self.constant)

    @property
    def eccentricity(self):
        """Eccentricity of any given parabola is one.
        """
        return S.One

    @property
    def focus_distance(self):
        """ The distance between the center and the focus.

        Examples
        ========

        >>> from sympy import Parabola, Point
        >>> a = Parabola(Point(3, 3), 2)
        >>> a.focus_distance
        2

        """
        return self.constant
