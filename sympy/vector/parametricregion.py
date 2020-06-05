from sympy.core.basic import Basic
from sympy.core.containers import Dict, Tuple

class ParametricRegion(Basic):
    """
    Represents a parametric region in space.

    Examples
    ========

    >>> from sympy import symbols, cos, sin, pi
    >>> from sympy.vector import ParametricRegion

    >>> r, theta = symbols("r theta")
    >>> circle = ParametricRegion((r, theta), r*cos(theta), r*sin(theta), limits={theta: (0, pi)})
    """
    def __new__(cls, parameters, definition, limits, system=None):
        if not isinstance(parameters, tuple):
            parameters = (parameters,)
        if not isinstance(definition, tuple):
            definition = (definition,)

        for parameter, bounds in limits.items():
            if parameter not in parameters:
                raise ValueError("%s is not listed in parameter tuple" % parameter)
            if len(bounds) != 2:
                raise ValueError("Bounds should be in the form (lower_bound, upper_bound)")

        obj = super().__new__(cls, Tuple(*parameters), Tuple(*definition), Dict(limits))
        obj._system = system
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

    @property
    def system(self):
        return self._system
