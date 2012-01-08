"""Curves in 2-dimensional Euclidean space.

Contains
========
Curve

"""

from sympy.core import sympify
from sympy.core.compatibility import is_sequence
from sympy.geometry.entity import GeometryEntity
from sympy.geometry.point import Point
from util import _symbol


class Curve(GeometryEntity):
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

    Examples
    ========

    >>> from sympy import sin, cos, Symbol
    >>> from sympy.abc import t
    >>> from sympy.geometry import Curve
    >>> C = Curve((sin(t), cos(t)), (t, 0, 2))
    >>> C.functions
    (sin(t), cos(t))
    >>> C.limits
    (t, 0, 2)
    >>> C.parameter
    t

    """

    def __new__(cls, function, limits):
        fun = sympify(function)
        if not is_sequence(fun) or len(fun) != 2:
            raise ValueError("Function argument should be (x(t), y(t)) but got %s" % str(function))
        if not is_sequence(limits) or len(limits) != 3:
            raise ValueError("Limit argument should be (t, tmin, tmax) but got %s" % str(limits))
        return GeometryEntity.__new__(cls, tuple(sympify(fun)), tuple(sympify(limits)))

    @property
    def free_symbols(self):
        """
        Return a set of symbols other than the bound symbols used to parametrically define the Curve.

        Examples
        ========

        >>> from sympy.abc import t, a
        >>> from sympy.geometry import Curve
        >>> Curve((t, t**2), (t, 0, 2)).free_symbols
        set()
        >>> Curve((t, t**2), (t, a, 2)).free_symbols
        set([a])
        """
        free = set()
        for a in self.functions + self.limits[1:]:
            free |= a.free_symbols
        free = free.difference(set([self.parameter]))
        return free

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

        >>> from sympy.abc import t
        >>> from sympy.geometry import Curve
        >>> C = Curve((t, t**2), (t, 0, 2))
        >>> C.functions
        (t, t**2)

        """
        return self.__getitem__(0)

    @property
    def parameter(self):
        """The curve function variable.

        Returns
        =======

        parameter : sympy symbol

        See Also
        ========

        functions

        Examples
        ========

        >>> from sympy.abc import t
        >>> from sympy.geometry import Curve
        >>> C = Curve([t, t**2], (t, 0, 2))
        >>> C.parameter
        t

        """
        return self.__getitem__(1)[0]

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

        >>> from sympy.abc import t
        >>> from sympy.geometry import Curve
        >>> C = Curve([t, t**3], (t, -2, 2))
        >>> C.limits
        (t, -2, 2)

        """
        return self.__getitem__(1)

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
        Point(2*t, t**2)
        >>> C.arbitrary_point(C.parameter)
        Point(2*s, s**2)
        >>> C.arbitrary_point(None)
        Point(2*s, s**2)
        >>> C.arbitrary_point(Symbol('a'))
        Point(2*a, a**2)

        """
        if parameter is None:
            return Point(*self.functions)

        tnew = _symbol(parameter, self.parameter)
        t = self.parameter
        if tnew.name != t.name and tnew.name in (f.name for f in self.free_symbols):
            raise ValueError('Symbol %s already appears in object and cannot be used as a parameter.' % tnew.name)
        return Point(*[w.subs(t, tnew) for w in self.functions])

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

        >>> from sympy import Curve, sin
        >>> from sympy.abc import x, t, s
        >>> Curve((x, sin(x)), (x, 1, 2)).plot_interval()
        [t, 1, 2]
        >>> Curve((x, sin(x)), (x, 1, 2)).plot_interval(s)
        [s, 1, 2]

        """
        t = _symbol(parameter, self.parameter)
        return [t] + list(self.limits[1:])
