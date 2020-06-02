from sympy.core.basic import Basic

class ParametricRegion(Basic):
    """
    Represents a parametric region in space.

    Examples
    ========

    >>> from sympy import symbols, cos, sin, pi
    >>> from sympy.vector import ParametricRegion
    >>> r, theta = symbols("r theta")
    >>> circle = ParametricRegion(r*cos(theta), r*sin(theta), (theta, 0, 2*pi))
    """
    def __new__(cls, *args, system=None):
        definition = list()
        bounds = list()
        parameters = list()

        for arg in args:
            if isinstance(arg, tuple):
                bounds.append(arg)
            else:
                definition.append(arg)

        for elem in bounds:
            parameters.append(elem[0])

        obj = super().__new__(cls, definition, bounds, parameters, system)
        obj._system = system

        return obj

    @property
    def bounds(self):
        return self.args[1]

    @property
    def definition(self):
        return self.args[0]

    @property
    def parameters(self):
        return self.args[2]

    @property
    def system(self):
        return self._system
