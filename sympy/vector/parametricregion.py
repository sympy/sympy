from sympy.core.basic import Basic

class ParametricRegion(Basic):
    """
    Represents a parametric region in space.
    """
    def __init__(self, system, *args):
        """
        Examples
        ========

        >>> from sympy.abc import r, theta
        >>> from sympy.vector import CoordSys3D, ParametricRegion
        >>> from sympy import cos, sin
        >>> C = CoordSys3D('C')
        >>> circle = ParametricRegion(C, r*cos(theta), r*sin(theta), (theta, 0, 2*pi))
        """
        self._definition = list()
        self._bounds = list()
        self._parameters = list()
        self.system = system

        for arg in args:
            if isinstance(arg, tuple):
                self.bounds.append(arg)
            else:
                self.definition.append(arg)

        if len(self.definition) > 3:
            raise ValueError("Only 3-D regions are supported")

        if len(self.bounds) > 3:
            raise ValueError("Only 3-D regions are supported")

        for elem in self.bounds:
            self.parameters.append(elem[0])

    @property
    def definition(self):
        return self._definition

    @property
    def bounds(self):
        return self._bounds

    @property
    def parameters(self):
        return self._parameters
