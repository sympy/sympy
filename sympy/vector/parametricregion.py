from sympy.core.basic import Basic
from sympy.core.containers import Dict, Tuple
from sympy.vector.coordsysrect import CoordSys3D

class ParametricRegion(Basic):
    """
    Represents a parametric region in space.

    Examples
    ========

    >>> from sympy import cos, sin, pi
    >>> from sympy.abc import r, theta, t, a, b
    >>> from sympy.vector import ParametricRegion

    >>> ParametricRegion(t, (t, t**2), limits={t: (-1, 2)})
    ParametricRegion((t,), (t, t**2), {t: (-1, 2)})
    >>> ParametricRegion((r, theta), (r*cos(theta), r*sin(theta)), {r: (-2, 2), theta: (0, pi)})
    ParametricRegion((r, theta), (r*cos(theta), r*sin(theta)), {r: (-2, 2), theta: (0, pi)})
    >>> ParametricRegion(t, (a*cos(t), b*sin(t)))
    ParametricRegion((t,), (a*cos(t), b*sin(t)), {})

    >>> circle = ParametricRegion((r, theta), (r*cos(theta), r*sin(theta)), {theta: (0, pi)})
    >>> circle.parameters
    (r, theta)
    >>> circle.definition
    (r*cos(theta), r*sin(theta))
    >>> circle.limits
    {theta: (0, pi)}

    Parameters
    ==========

    parameters_or_coordsys : parameter or a tuple of parameters or a CoordSys3d object.
                        When a CoordSys3d object is passed, its base scalars are used as parameters.

    definition : tuple to define base scalars in terms of parameters.

    limits : dict to define bounds of each parameter.
            Each Key of dictionary should be parameter and value
            is a tuple to represent corresponding lower and upper bound.`

    """
    def __new__(cls, parameters_or_coordsys, definition, limits=None):

        if isinstance(parameters_or_coordsys, CoordSys3D):
            parameters = parameters_or_coordsys.base_scalars()
        elif not (isinstance(parameters_or_coordsys, tuple) or isinstance(parameters_or_coordsys, Tuple)):
            parameters = (parameters_or_coordsys,)
        else:
            parameters = parameters_or_coordsys

        if not (isinstance(definition, tuple) or isinstance(definition, Tuple)):
            definition = (definition,)

        if limits is None:
            limits = {}

        for parameter, bounds in limits.items():
            if parameter not in parameters:
                raise ValueError("%s is not listed in parameter tuple" % parameter)
            if len(bounds) != 2:
                raise ValueError("Bounds should be in the form (lower_bound, upper_bound)")

        obj = super().__new__(cls, Tuple(*parameters), Tuple(*definition), Dict(limits))
        return obj

    @property
    def limits(self):
        return self.args[2]

    @property
    def definition(self):
        return self.args[1]

    @property
    def parameters(self):
        return self.args[0]
