"""Curves in 2-dimensional Euclidean space.

Contains
--------
Curve

"""

from sympy.core import sympify
from sympy.geometry.exceptions import GeometryError
from entity import GeometryEntity

class Curve(GeometryEntity):
    """A curve in space.

    A curve is defined by parametric functions for the coordinates, a
    parameter and the lower and upper bounds for the parametr value.

    Parameters
    ----------
    function : list of functions
    limits : 3-tuple
        Function parameter and lower and upper bounds.

    Attributes
    ----------
    functions
    parameter
    limits

    Raises
    ------
    GeometryError
        When `functions`
    ValueError
        When `limits` are specified incorrectly.


    Examples
    --------
    >>> from sympy import sin, cos, Symbol
    >>> from sympy.abc import t
    >>> from sympy.geometry import Curve
    >>> C = Curve([sin(t), cos(t)], (t, 0, 2))
    >>> C.functions
    [sin(t), cos(t)]
    >>> C.limits
    (t, 0, 2)
    >>> C.parameter
    t

    """

    def __new__(cls, function, limits):
        fun = sympify(function)
        if not fun:
            raise GeometryError("%s.__new__ don't know how to handle" % cls.__name__)
        if not isinstance(limits, (list, tuple)) or len(limits) != 3:
            raise ValueError("Limits argument has wrong syntax")
        return GeometryEntity.__new__(cls, fun, limits)

    @property
    def functions(self):
        """The functions specifying the curve.

        Returns
        -------
        functions : list of parametrised coordinate functions.

        Examples
        --------
        >>> from sympy.abc import t
        >>> from sympy.geometry import Curve
        >>> C = Curve([t, t**2], (t, 0, 2))
        >>> C.functions
        [t, t**2]

        """
        return self.__getitem__(0)

    @property
    def parameter(self):
        """The curve function variable.

        Returns
        -------
        parameter : sympy symbol

        Examples
        --------
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
        -------
        limits : tuple
            Contains parameter and lower and upper limits.

        Examples
        --------
        >>> from sympy.abc import t
        >>> from sympy.geometry import Curve
        >>> C = Curve([t, t**3], (t, -2, 2))
        >>> C.limits
        (t, -2, 2)

        """
        return self.__getitem__(1)
