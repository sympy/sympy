"""Curves in n-dimensional Euclidean space.

Contains
========
Curve

"""

from __future__ import division, print_function

from sympy import sqrt
from sympy.core import sympify, diff
from sympy.core.compatibility import is_sequence
from sympy.core.containers import Tuple
from sympy.core.symbol import _symbol
from sympy.geometry.entity import GeometryEntity, GeometrySet
from sympy.geometry.point import Point
from sympy.integrals import integrate


class Curve(GeometrySet):
    """A curve in space.

    A curve is defined by parametric functions for the coordinates, a
    parameter and the lower and upper bounds for the parameter value.

    Parameters
    ==========

    function : list of functions
    limits : 3-tuple
        Function parameter and lower and upper bounds.

    Attributes
    ==========

    functions
    parameter
    limits

    Raises
    ======

    ValueError
        When `functions` are specified incorrectly.
        When `limits` are specified incorrectly.

    See Also
    ========

    sympy.core.function.Function
    sympy.polys.polyfuncs.interpolate

    Examples
    ========

    2-Dimensional Curves:
    >>> from sympy import sin, cos, Symbol, interpolate
    >>> from sympy.abc import t, a
    >>> from sympy.geometry import Curve
    >>> C = Curve((sin(t), cos(t)), (t, 0, 2))
    >>> C.functions
    (sin(t), cos(t))
    >>> C.limits
    (t, 0, 2)
    >>> C.parameter
    t
    >>> C = Curve((t, interpolate([1, 4, 9, 16], t)), (t, 0, 1)); C
    Curve((t, t**2), (t, 0, 1))
    >>> C.subs(t, 4)
    Point2D(4, 16)
    >>> C.arbitrary_point(a)
    Point2D(a, a**2)

    Higher Dimmension Curves:
    >>> from sympy import sin, cos, Symbol
    >>> from sympy.abc import t, a
    >>> C = Curve((sin(t), cos(t), t), (t, 0, 2))
    >>> C
    Curve((sin(t), cos(t), t), (t, 0, 2))
    >>> C.functions
    (sin(t), cos(t), t)
    >>> C.limits
    (t, 0, 2)
    >>> C.parameter
    t
    >>> C.subs(t,4)
    Point3D(sin(4), cos(4), 4)
    
    """

    def __new__(cls, function, limits):
        fun = sympify(function)
        if not is_sequence(limits) or len(limits) != 3:
            raise ValueError("Limit argument should be (t, tmin, tmax) "
                "but got %s" % str(limits))

        return GeometryEntity.__new__(cls, Tuple(*fun), Tuple(*limits))

    def _eval_subs(self, old, new):
        """The function specifying the point after evaluating at the new value.

        Examples
        ========

        2-Dimensional Curves:
        >>> from sympy import sin, cos, Symbol, interpolate
        >>> from sympy.abc import t, a
        >>> from sympy.geometry import Curve
        >>> C = Curve((sin(t), cos(t)), (t, 0, 2))
        >>> C.subs(t, 4)
        Point2D(4, 16)

        Higher Dimmension Curves:
        >>> from sympy.abc import t
        >>> C = Curve((t, t**2, t**3), (t, 0, 2))
        >>> C.subs(t,4)
        Point3D(4, 16, 64)

        """
        if old == self.parameter:
            return Point(*[f.subs(old, new) for f in self.functions])

    def arbitrary_point(self, parameter='t'):
        """
        A parameterized point on the curve.

        Parameters
        ==========

        parameter : str or Symbol, optional
            Default value is 't';
            the Curve's parameter is selected with None or self.parameter
            otherwise the provided symbol is used.

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

        >>> from sympy import Symbol
        >>> from sympy.abc import s
        >>> from sympy.geometry import Curve
        >>> C = Curve([2*s, s**2], (s, 0, 2))
        >>> C.arbitrary_point()
        Point2D(2*t, t**2)
        >>> C.arbitrary_point(C.parameter)
        Point2D(2*s, s**2)
        >>> C.arbitrary_point(None)
        Point2D(2*s, s**2)
        >>> C.arbitrary_point(Symbol('a'))
        Point2D(2*a, a**2)

        """
        if parameter is None:
            return Point(*self.functions)

        tnew = _symbol(parameter, self.parameter, real=True)
        t = self.parameter
        if (tnew.name != t.name and
                tnew.name in (f.name for f in self.free_symbols)):
            raise ValueError('Symbol %s already appears in object '
                'and cannot be used as a parameter.' % tnew.name)
        return Point(*[w.subs(t, tnew) for w in self.functions])

    @property
    def free_symbols(self):
        """
        Return a set of symbols other than the bound symbols used to
        parametrically define the Curve.

        Examples
        ========

        2-Dimensional Curves:
        >>> from sympy.abc import t, a
        >>> from sympy.geometry import Curve
        >>> Curve((t, t**2), (t, 0, 2)).free_symbols
        set()
        >>> Curve((t, t**2), (t, a, 2)).free_symbols
        {a}

        Higher Dimmension Curves:
        >>> from sympy.abc import t, a
        >>> C = Curve((t, t**2, t**3), (t, 0, 2))
        >>> C.free_symbols
        set([])
        >>> Curve((t, t**2,t**3), (t, a, 2)).free_symbols
        set([a])
        
        """
        free = set()
        for a in self.functions + self.limits[1:]:
            free |= a.free_symbols
        free = free.difference({self.parameter})
        return free

    @property
    def ambient_dimension(self):
        """The dimmension of the curve

        Examples
        ========

        2-Dimensional Curves:
        >>> from sympy.abc import t
        >>> from sympy.geometry import Curve
        >>> C = Curve((t, t**2), (t, 0, 2))
        >>> C.ambient_dimension
        2

        Higher Dimmension Curves:
        >>> from sympy.abc import t
        >>> C = Curve((t, t**2, t**3), (t, 0, 2))
        >>> C.ambient_dimension
        3

        """
        return len(self.args[0])

    @property
    def functions(self):
        """The functions specifying the curve.

        Returns
        =======

        functions : list of parameterized coordinate functions.

        See Also
        ========

        parameter

        Examples
        ========

        2-Dimensional Curves:
        >>> from sympy.abc import t
        >>> from sympy.geometry import Curve
        >>> C = Curve((t, t**2), (t, 0, 2))
        >>> C.functions
        (t, t**2)

        Higher Dimmension Curves:
        >>> from sympy.abc import t
        >>> C = Curve((t, t**2, t**3), (t, 0, 2))
        >>> C.functions
        (t, t**2, t**3)

        """
        return self.args[0]

    @property
    def limits(self):
        """The limits for the curve.

        Returns
        =======

        limits : tuple
            Contains parameter and lower and upper limits.

        See Also
        ========

        plot_interval

        Examples
        ========

        2-Dimensional Curves:
        >>> from sympy.abc import t
        >>> from sympy.geometry import Curve
        >>> C = Curve([t, t**3], (t, -2, 2))
        >>> C.limits
        (t, -2, 2)

        Higher Dimmension Curves:
        >>> from sympy.abc import t
        >>> C = Curve((t, t**2, t**3), (t, 0, 2))
        >>> C.limits
        (t, 0, 2)

        """
        return self.args[1]

    @property
    def parameter(self):
        """The curve function variable.

        Returns
        =======

        parameter : SymPy symbol

        See Also
        ========

        functions

        Examples
        ========

        2-Dimensional Curves:
        >>> from sympy.abc import t
        >>> from sympy.geometry import Curve
        >>> C = Curve([t, t**2], (t, 0, 2))
        >>> C.parameter
        t

        Higher Dimmension Curves:
        >>> from sympy.abc import t
        C = Curve((t, t**2, t**3), (t, 0, 2))
        >>> C.parameter
        t

        """
        return self.args[1][0]

    @property
    def length(self):
        """The curve length.

        Examples
        ========

        2-Dimensional Curves:
        >>> from sympy.geometry.curve import Curve
        >>> from sympy import cos, sin
        >>> from sympy.abc import t
        >>> Curve((t, t), (t, 0, 1)).length
        sqrt(2)

        Higher Dimmension Curves:
        >> from sympy.abc import t
        >>> C = Curve((t, t, t), (t, 0, 2))
        >>> C.length
        2*sqrt(3)
        
        """
        integrand = sqrt(sum(diff(func, self.limits[0])**2 for func in self.functions))
        return integrate(integrand, self.limits)

    def plot_interval(self, parameter='t'):
        """The plot interval for the default geometric plot of the curve.

        Parameters
        ==========

        parameter : str or Symbol, optional
            Default value is 't';
            otherwise the provided symbol is used.

        Returns
        =======

        plot_interval : list (plot interval)
            [parameter, lower_bound, upper_bound]

        See Also
        ========

        limits : Returns limits of the parameter interval

        Examples
        ========

        2-Dimensional Curves:
        >>> from sympy import Curve, sin
        >>> from sympy.abc import x, t, s
        >>> Curve((x, sin(x)), (x, 1, 2)).plot_interval()
        [t, 1, 2]
        >>> Curve((x, sin(x)), (x, 1, 2)).plot_interval(s)
        [s, 1, 2]

        Higher Dimmension Curves:

        """
        t = _symbol(parameter, self.parameter, real=True)
        return [t] + list(self.limits[1:])

    def rotate(self, angle=0, pt=None):
        """Rotate ``angle`` radians counterclockwise about Point ``pt``.

        The default pt is the origin, Point(0, 0).

        Examples
        ========

        >>> from sympy.geometry.curve import Curve
        >>> from sympy.abc import x
        >>> from sympy import pi
        >>> Curve((x, x), (x, 0, 1)).rotate(pi/2)
        Curve((-x, x), (x, 0, 1))
        """
        from sympy.matrices import Matrix, rot_axis3
        if pt:
            pt = -Point(pt, dim=2)
        else:
            pt = Point(0,0)
        rv = self.translate(*pt.args)
        f = list(rv.functions)
        f.append(0)
        f = Matrix(1, 3, f)
        f *= rot_axis3(angle)
        rv = self.func(f[0, :2].tolist()[0], self.limits)
        if pt is not None:
            pt = -pt
            return rv.translate(*pt.args)
        return rv

    def scale(self, scaling, pt=None):
        """Override GeometryEntity.scale since Curve is not made up of Points.

        Examples
        ========
        
        2-Dimensional Curves:
        >>> from sympy.geometry.curve import Curve
        >>> from sympy import pi
        >>> from sympy.abc import x
        >>> Curve((x, x), (x, 0, 1)).scale([2,1])
        Curve((2*x, x), (x, 0, 1))

        Higher Dimmension Curves:
        >>> from sympy.abc import x, t, s
        >>> C = Curve((t, t, t), (t, 0, 2))
        >>> C.scale([2,0,4])
        Curve((2*t, 0, 4*t), (t, 0, 2))
        
        """
        if pt:
            pt = Point(pt, dim=2)
            return self.translate(*(-pt).args).scale(x, y).translate(*pt.args)
        results = []
        for i in range(len(self.args[0])):
            results.append(self.args[0][i] * scaling[i])
        return self.func(results, self.limits)

    def translate(self, translations):
        """ Now supports nDimmensional Curves
        Translate the Curve by [translations].

        Examples
        ========
        
        2-Dimensional Curves:
        >>> from sympy.geometry.curve import Curve
        >>> from sympy import pi
        >>> from sympy.abc import x
        >>> Curve((x, x), (x, 0, 1)).translate(1, 2)
        Curve((x + 1, x + 2), (x, 0, 1))

        Higher Dimmension Curves:
        >>> from sympy.abc import t
        >>> C = Curve((t, t, t), (t, 0, 2))
        >>> C.translate([1,2,3])
        Curve((t + 1, t + 2, t + 3), (t, 0, 2))
        
        """
        results = []
        for i in range(len(self.args[0])):
            results.append(self.args[0][i] + translations[i])
        return self.func(results, self.limits)
