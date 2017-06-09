from sympy.core.basic import Basic
from sympy import sqrt, sin, cos
from sympy.core.function import Derivative


class CoeffProvider(Basic):
    """
    It provides specific information of curvilinear coordinate
    system to CoordSysCartesian.

    Ideally, users should not instantiate this class.

    """

    def __new__(cls, system, curv_coord_name):
        obj = super(CoeffProvider, cls).__new__(cls, system)
        obj.sys = system
        obj._x, obj._y, obj._z = system.x, system.y, system.z
        obj.lame_coefficient = obj._coefficient_mapping(curv_coord_name)
        obj.equations = obj._equation_mapping(curv_coord_name)
        return obj

    def _coefficient_mapping(self, curv_coord_name):
        try:
            coefficient_mapping = {
                'cartesian': (1, 1, 1),
                'spherical': (1, self._x, self._x * cos(self._y))
            }
        except:
            raise ValueError('Wrong set of parameters')
        return coefficient_mapping[curv_coord_name]

    def get_lame_coefficient(self):
        return self.lame_coefficient

    def h1(self):
        return self.lame_coefficient[0]

    def h2(self):
        return self.lame_coefficient[1]

    def h3(self):
        return self.lame_coefficient[2]
