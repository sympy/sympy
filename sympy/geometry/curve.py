from sympy.core import sympify
from sympy.geometry.exceptions import GeometryError
from entity import GeometryEntity

class Curve(GeometryEntity):
    """
    A curve in space.

    Example:
    ========
        >>> from sympy import sin, cos, Symbol
        >>> t = Symbol("t")
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
            raise GeometryError("%s.__new__ don't know how to handle" % cls.__name__);
        if not isinstance(limits, (list, tuple)) or len(limits) != 3:
            raise ValueError("Limits argument has wrong syntax");
        return GeometryEntity.__new__(cls, fun, limits)

    @property
    def functions(self):
        """The functions specifying the curve."""
        return self.__getitem__(0)

    @property
    def parameter(self):
        """The curve function variable."""
        return self.__getitem__(1)[0]

    @property
    def limits(self):
        """The limits for the curve."""
        return self.__getitem__(1)
