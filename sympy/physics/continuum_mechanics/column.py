from sympy.solvers import linsolve, solve, solveset
from sympy.core import Symbol, diff, symbols
from sympy import dsolve, Function, Derivative, Eq

class Column(object):
    def __init__(self, height, elastic_modulus, second_moment, load, eccentricity=None, top="pinned", bottom="pinned", boundary_conditions=None):
        self._height = height
        self._elastic_modulus = elastic_modulus
        self._second_moment = second_moment
        self._load = load
        self._eccentricity = eccentricity
        self._moment = 0
        self._end_conditions = {'top':top, 'bottom': bottom}
        self._boundary_conditions = {'deflection': [], 'slope': []}
        self._variable = Symbol('x')
        self._deflection = None
        self._slope = None

    @property
    def height(self):
        return self._height

    @property
    def elastic_modulus(self):
        return self._elastic_modulus

    @property
    def second_moment(self):
        return self._second_moment

    @property
    def load(self):
        return self._load

    @property
    def eccentricity(self):
        return self._eccentricity

    @property
    def end_conditions(self):
        return self._end_conditions

    @property
    def boundary_conditions(self):
        return self._boundary_conditions


    def _apply_load_conditions(self):
        y = Function('y')
        x = self._variable
        P = Symbol('P')

        self._moment += P*y(x)
        if self.eccentricity:
            self._moment += P*eccentricity

        if self._end_conditions['top'] == "fixed" and self._end_conditions['bottom'] == "fixed":
            # `M` is the reaction moment
            M = Symbol('M')
            self._boundary_conditions['deflection'].append((self._height, 0))
            # moment  = P*y - M
            self._moment -= M

        if self._end_conditions['top'] == "pinned" and self._end_conditions['bottom'] == "fixed":
            # `F` is the horizontal force at the pinned end to counter the rection moment at fixed end
            F = Symbol('F')
            self._boundary_conditions['deflection'].append((self._height, 0))
            # moment = P*y - F(l - x)
            self._moment -= F*(self.height - x)

        if self._end_conditions['top'] == "free" and self._end_conditions['bottom'] == "fixed":
            # `d` is the deflection at the free end
            d = Symbol('d')
            self._boundary_conditions['deflection'].append((self._height, d))
            # moment = P*y - P*d
            self._moment -= self.load*d


    def solve_slope_deflection(self):
        y = Function('y')
        x = self._variable

        self._apply_load_conditions()
        C1, C2 = symbols('C1, C2')
        E = self._elastic_modulus
        I = self._second_moment

        # differnetial equation of Column buckling
        eq = E*I*y(x).diff(x, 2) + self._moment

        defl_sol = dsolve(Eq(eq, 0), y(x)).args[1]
        slope_sol = defl_sol.diff(x)

        constants = list(linsolve([defl_sol.subs(x, 0), slope_sol.subs(x, 0)], C1, C2).args[0])

        self._deflection = defl_sol.subs({C1: constants[0], C2:constants[1]})
        self._slope = slope_sol.subs({C1: constants[0], C2:constants[1]})


    def moment(self):
        self._apply_load_conditions()
        return self._moment

    def slope(self):
        return self._slope

    def deflection(self):
        return self._deflection
