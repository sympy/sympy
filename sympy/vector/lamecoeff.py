from sympy.core import Expr
from sympy import sqrt, simplify, sin, cos, symbols, diff
from sympy.core.function import Derivative

x, y, z = symbols('x y z')


class CoeffProvider(Expr):
    """
    It provides specific information of curvilinear coordinate
    system to CoordSysCartesian.

    Ideally, users should not instantiate this class.

    """

    def __new__(cls, system):
        obj = super(CoeffProvider, cls).__new__(cls, system)
        obj.sys = system
        obj._x, obj._y, obj._z = system.x, system.y, system.z
        obj.coordinates_mapping = {
            'cartesian': CartesianCoeff,
            'spherical': SphericalCoeff
        }
        return obj

    def get_coefficients(self, *args, **kwargs):
        if len(args) == 1:
            return self.coordinates_mapping[args[0]](self.sys)
        elif len(kwargs) == 3:
            self.eq1, self.eq2, self.eq3 = kwargs['eq1'], kwargs['eq2'], kwargs['eq3']
            return CurvilinearCoeff(self.sys, **kwargs)
        else:
            raise ValueError('Wrong set of parameters')


class LameCoeff:
    """
    Base class for every predefined class with Lame coefficient.

    Ideally, users should not instantiate this class.

    """
    def __init__(self, system,  *args, **kwargs):
        pass

    def transformation_equations(self):
        pass

    def get_lame_coefficient(self):
        return self.h1, self.h2, self.h3

    def h1(self):
        pass

    def h2(self):
        pass

    def h3(self):
        pass


class CartesianCoeff(LameCoeff):
    """
    Class for Cartesian coordinate system.

    Ideally, users should not instantiate this class.
    """
    def transformation_equations(self):
        return 1, 1, 1

    def h1(self):
        return 1

    def h2(self):
        return 1

    def h3(self):
        return 1


class SphericalCoeff(LameCoeff):
    """
    Class for Spherical coordinate system.

    Ideally, users should not instantiate this class.
    """
    def __init__(self, system,  *args, **kwargs):
        LameCoeff.__init__(self, system, *args, **kwargs)
        self.x, self.y, self.z = system.x, system.y, system.z

    def transformation_equations(self):
        return self.x * sin(self.y) * cos(self.z), \
               self.x * sin(self.y) * sin(self.z), \
               self.x * cos(self.y)

    def h1(self):
        return 1

    def h2(self):
        return self.x

    def h3(self):
        return self.x * cos(self.y)


class CurvilinearCoeff(LameCoeff):
    """
    General implementation of curvilinear coordinate system.

    Ideally, users should not instantiate this class.

    """
    def __init__(self, system,  *args, **kwargs):
        LameCoeff.__init__(self, system, *args, **kwargs)
        self.x, self.y, self.z = system.x, system.y, system.z
        try:
            self.eq1, self.eq2, self.eq3 = kwargs['eq1'], \
                                           kwargs['eq2'], \
                                           kwargs['eq3']
        except KeyError:
            raise ValueError('Wrong set of parameters')

    def transformation_equations(self):
        return self.eq1, self.eq2, self.eq3

    def h1(self):
        return sqrt(Derivative(self.eq1, self.x) ** 2 +
                    Derivative(self.eq2, self.x) ** 2 +
                    Derivative(self.eq3, self.x) ** 2)

    def h2(self):
        return sqrt(Derivative(self.eq1, self.y) ** 2 +
                    Derivative(self.eq2, self.y) ** 2 +
                    Derivative(self.eq3, self.y) ** 2)

    def h3(self):
        return sqrt(Derivative(self.eq1, self.z) ** 2 +
                    Derivative(self.eq2, self.z) ** 2 +
                    Derivative(self.eq3, self.z) ** 2)
