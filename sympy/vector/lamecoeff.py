from sympy.core.basic import Basic
from sympy import sqrt, sin, cos
from sympy.core.function import Derivative


class CoeffProvider(Basic):
    """
    It provides specific information of curvilinear coordinate
    system to CoordSysCartesian.

    Ideally, users should not instantiate this class.

    """

    def __new__(cls, system, curv_coord_name=None, transformation_equations=None):
        obj = super(CoeffProvider, cls).__new__(cls, system)
        obj.sys = system
        obj._x, obj._y, obj._z = system.x, system.y, system.z
        if curv_coord_name is not None:
            class DefinedCoefficient:
                coefficient_mapping = {
                    'cartesian': (1, 1, 1),
                    'spherical': (1, obj._x, obj._x * cos(obj._y))
                }
                equation_mapping = {
                    'cartesian': (1, 1, 1),
                    'spherical': (obj._x * sin(obj._y) * cos(obj._z),
                                  obj._x * sin(obj._y) * sin(obj._z),
                                  obj._x * cos(obj._y))
                }

                def __init__(self, name):
                    self.trans_equations = self.equation_mapping[name]
                    self.lame_coefficient = self.coefficient_mapping[name]

            obj._curvilinear_parameters = DefinedCoefficient(curv_coord_name)

        elif transformation_equations is not None:
            class CurvilinearCoefficient:
                def __init__(self, equations):
                    self.eq1, self.eq2, self.eq3 = equations
                    self.lame_coefficient = (
                        sqrt(Derivative(self.eq1, obj._x) ** 2 +
                             Derivative(self.eq2, obj._x) ** 2 +
                             Derivative(self.eq3, obj._x) ** 2),
                        sqrt(Derivative(self.eq1, obj._y) ** 2 +
                             Derivative(self.eq2, obj._y) ** 2 +
                             Derivative(self.eq3, obj._y) ** 2),
                        sqrt(Derivative(self.eq1, obj._z) ** 2 +
                             Derivative(self.eq2, obj._z) ** 2 +
                             Derivative(self.eq3, obj._z) ** 2)
                    )
                    self.trans_equations = equations

            obj._curvilinear_parameters = CurvilinearCoefficient(transformation_equations)

        else:
            raise ValueError('Wrong set of parameters')

        obj.lame_coefficient = obj._curvilinear_parameters.lame_coefficient
        obj.trans_equations = obj._curvilinear_parameters.trans_equations
        return obj

    def transformation_equations(self):
        return self.trans_equations

    def get_lame_coefficient(self):
        return self.lame_coefficient

    def h1(self):
        return self.lame_coefficient[0]

    def h2(self):
        return self.lame_coefficient[1]

    def h3(self):
        return self.lame_coefficient[2]
