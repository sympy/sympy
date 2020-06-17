from sympy.core.basic import Basic
from sympy.core.containers import Tuple

class ParametricRegion(Basic):
    """
    Represents a parametric region in space.

    Examples
    ========

    >>> from sympy import cos, sin, pi
    >>> from sympy.abc import r, theta, t, a, b, x, y
    >>> from sympy.vector import ParametricRegion

    >>> ParametricRegion((t, t**2), (t, -1, 2))
    ParametricRegion((t, t**2), (t, -1, 2))
    >>> ParametricRegion((x, y), (x, 3, 4), (y, 5, 6))
    ParametricRegion((x, y), (x, 3, 4), (y, 5, 6))
    >>> ParametricRegion((r*cos(theta), r*sin(theta)), (r, -2, 2), (theta, 0, pi))
    ParametricRegion((r*cos(theta), r*sin(theta)), (r, -2, 2), (theta, 0, pi))
    >>> ParametricRegion((a*cos(t), b*sin(t)), t)
    ParametricRegion((a*cos(t), b*sin(t)), t)

    >>> circle = ParametricRegion((r*cos(theta), r*sin(theta)), r, (theta, 0, pi))
    >>> circle.parameters
    (r, theta)
    >>> circle.definition
    (r*cos(theta), r*sin(theta))
    >>> circle.limits
    {theta: (0, pi)}

    Dimension of a parametric region determines whether a region is a curve, surface
    or volume region. It does not represent its dimensions in space.
    >>> circle.dimensions
    1

    Parameters
    ==========

    definition : tuple to define base scalars in terms of parameters.

    bounds : Parameter or a tuple of length 3 to define parameter and
            corresponding lower and upper bound
    """
    def __new__(cls, definition, *bounds):
        parameters = ()
        limits = {}

        if not isinstance(bounds, Tuple):
            bounds = Tuple(*bounds)

        for bound in bounds:
            if  isinstance(bound, tuple) or isinstance(bound, Tuple):
                if len(bound) != 3:
                    raise ValueError("Tuple should be in the form (parameter, lowerbound, upperbound)")
                parameters += (bound[0],)
                limits[bound[0]] = (bound[1], bound[2])
            else:
                parameters += (bound,)

        if not (isinstance(definition, tuple) or isinstance(definition, Tuple)):
            definition = (definition,)

        obj = super().__new__(cls, Tuple(*definition), *bounds)
        obj._parameters = parameters
        obj._limits = limits

        return obj

    @property
    def definition(self):
        return self.args[0]

    @property
    def limits(self):
        return self._limits

    @property
    def parameters(self):
        return self._parameters

    @property
    def dimensions(self):
        return len(self.limits)
