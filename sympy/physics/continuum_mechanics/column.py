from sympy.solvers import linsolve, solve
from sympy.core import Symbol, diff, symbols
from sympy import dsolve, Function, Derivative, Eq, cos, sin, sqrt, tan
from sympy.core.symbol import Dummy
from sympy.printing import sstr

class Column(object):
    """
    A column is a structural member designed to undertake axial
    compressive loads. A column is characterized by its
    cross-sectional profile(second moment of area), its length and
    its material.

    Examples
    ========

    There is a solid round bar 3 m long with second-moment I is used as a
    column with both the ends pinned. Young's modulus of the Column is E.
    The buckling load applied is 78KN

    >>> from sympy.physics.continuum_mechanics.column import Column
    >>> from sympy import Symbol, symbols
    >>> E, I, P = symbols('E, I, P', positive=True)
    >>> c = Column(3, E, I, 78000, top="pinned", bottom="pinned")
    >>> c.end_conditions
    {'bottom': 'pinned', 'top': 'pinned'}
    >>> c.boundary_conditions
    {'deflection': [(0, 0), (3, 0)], 'slope': [(0, 0)]}
    >>> c.moment()
    78000*y(x)
    >>> c.solve_slope_deflection()
    >>> c.deflection()
    C1*sin(20*sqrt(195)*x/(sqrt(E)*sqrt(I)))
    >>> c.slope()
    20*sqrt(195)*C1*cos(20*sqrt(195)*x/(sqrt(E)*sqrt(I)))/(sqrt(E)*sqrt(I))
    >>> c.critical_load()
    pi**2*E*I/9
    """
    def __init__(self, height, elastic_modulus, second_moment, load, eccentricity=None, top="pinned", bottom="pinned", bc_slope=None, bc_deflection=None):
        """
        Parameters
        ==========

        height: Sympifyable
            A symbol or a value representing column's height

        elastic_modulus: Sympifyable
            A symbol or a value representing the Column's modulus of
            elasticity. It is a measure of the stiffness of the Column
            material.

        second_moment: Sympifyable
            A symbol or a value representing Column's second-moment of area
            It is a geometrical property of an area which reflects how its
            points are distributed with respect to its neutral axis.

        load: Sympifyable
            A symbol or a value representing the load applied on the Column.

        eccentricity: Sympifyable (default=None)
            A symbol or a value representing the eccentricity of the load
            applied. Eccentricity is the distance of the point of application
            of load from the neutral axis.

        top: string (default="pinned")
            A string representing the top-end condition of the column.
            It can be: pinned
                       fixed
                       free

        bottom: string (default="pinned")
            A string representing the bottom-end condition of the column.
            It can be: pinned
                       fixed

        bc_slope: list of tuples
            A list of tuples representing the boundary conditions of slope.
            The tuple takes two elements `location` and `value`.

        bc_deflection: list of tuples
            A list of tuples representing the boundary conditions of deflection
            The tuple consists of two elements `location` and `value`.
        """
        self._height = height
        self._elastic_modulus = elastic_modulus
        self._second_moment = second_moment
        self._load = load
        self._eccentricity = eccentricity
        self._moment = 0
        self._end_conditions = {'top':top, 'bottom': bottom}
        self._boundary_conditions = {'deflection': [], 'slope': []}
        if bc_deflection:
            self._boundary_conditions['deflection'] = bc_deflection
        if bc_slope:
            self._boundary_conditions['slope'] = bc_slope
        self._variable = Symbol('x')
        self._deflection = None
        self._slope = None
        self._critical_load = None
        self._apply_load_conditions()

    def __str__(self):
        str_sol = 'Column({}, {}, {})'.format(sstr(self._height), sstr(self._elastic_modulus), sstr(self._second_moment))
        return str_sol

    @property
    def height(self):
        """Height of the column"""
        return self._height

    @property
    def elastic_modulus(self):
        """Elastic modulus of the column"""
        return self._elastic_modulus

    @property
    def second_moment(self):
        """Second moment of the column"""
        return self._second_moment

    @property
    def load(self):
        """Load applied on the column"""
        return self._load

    @property
    def eccentricity(self):
        """Eccentricity of the load applied on the column"""
        return self._eccentricity

    @property
    def end_conditions(self):
        """End-conditions in the form of a dictionary"""
        return self._end_conditions

    @property
    def boundary_conditions(self):
        """Boundary conditions in the form of a dictionary"""
        return self._boundary_conditions


    def _apply_load_conditions(self):
        y = Function('y')
        x = self._variable
        P = Symbol('P', positive=True)

        self._moment += P*y(x)
        if self.eccentricity:
            self._moment += P*eccentricity

        # Initial boundary conditions, considering slope and deflection
        # the bottom always zero
        self._boundary_conditions['deflection'].append((0, 0))
        self._boundary_conditions['slope'].append((0, 0))

        if self._end_conditions['top'] == "pinned" and self._end_conditions['bottom'] == "pinned":
            self._boundary_conditions['deflection'].append((self._height, 0))

        elif self._end_conditions['top'] == "fixed" and self._end_conditions['bottom'] == "fixed":
            # `M` is the reaction moment
            M = Symbol('M')
            self._boundary_conditions['deflection'].append((self._height, 0))
            # moment  = P*y - M
            self._moment -= M

        elif self._end_conditions['top'] == "pinned" and self._end_conditions['bottom'] == "fixed":
            # `F` is the horizontal force at the pinned end to counter the rection moment at fixed end
            F = Symbol('F')
            self._boundary_conditions['deflection'].append((self._height, 0))
            # moment = P*y - F(l - x)
            self._moment -= F*(self.height - x)

        elif self._end_conditions['top'] == "free" and self._end_conditions['bottom'] == "fixed":
            # `d` is the deflection at the free end
            d = Symbol('d')
            self._boundary_conditions['deflection'].append((self._height, d))
            # moment = P*y - P*d
            self._moment -= P*d

        else:
            raise ValueError("{} {} end-condition is not supported".format(sstr(self._end_conditions['top']), sstr(self._end_conditions['bottom'])))


    def solve_slope_deflection(self):
        """
        Solves the differnetial equation of buckling to determine the
        deflection and slope equations.

        Examples
        ========

        A column of timber section is 8m long both ends being fixed.
        Young's modulus of the timber is E and second-moment I.
        The applied load on the column is 360KN.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy import Symbol, symbols
        >>> E, I, P = symbols('E, I, P', positive=True)
        >>> c = Column(8, E, I, 360000, top="fixed", bottom="fixed")
        >>> c.end_conditions
        {'bottom': 'fixed', 'top': 'fixed'}
        >>> c.boundary_conditions
        {'deflection': [(0, 0), (8, 0)], 'slope': [(0, 0)]}
        >>> c.moment()
        -M + 360000*y(x)
        >>> c.solve_slope_deflection()
        >>> c.deflection()
        -M*cos(600*x/(sqrt(E)*sqrt(I)))/360000 + M/360000
        >>> c.slope()
        M*sin(600*x/(sqrt(E)*sqrt(I)))/(600*sqrt(E)*sqrt(I))
        """
        y = Function('y')
        x = self._variable
        P = Symbol('P', positive=True)

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

        # if deflection is zero, no buckling occurs, which is not the case,
        # so trying to solve for the constants differently
        if self._deflection == 0:
            self._deflection = defl_sol
            self._slope = slope_sol

            defl_eqs = []
            # taking last two bounndary conditions which are actually
            # the initial boundary conditions.
            for point, value in self._boundary_conditions['deflection'][-2:]:
                defl_eqs.append(self._deflection.subs(x, point) - value)

            # solve for C1, C2 along with P
            solns = solve(defl_eqs, (P, C1, C2), dict=True)
            for sol in solns:
                if self._deflection.subs(sol) == 0:
                    # removing trivial solutions
                    solns.remove(sol)

            # checking if the constants are solved, and subtituting them in
            # the deflection and slope equation
            if C1 in solns[0].keys():
                self._deflection = self._deflection.subs(C1, solns[0][C1])
                self._slope = self._slope.subs(C1, solns[0][C1])
            if C2 in solns[0].keys():
                self._deflection = self._deflection.subs(C2, solns[0][C2])
                self._slope = self._slope.subs(C2, solns[0][C2])
            if P in solns[0].keys():
                self._critical_load = solns[0][P]


    def critical_load(self):
        """
        Detrmines the critical load (for single bow buckling condition) of
        the given column under the given conditions.

        Examples
        ========

        A solid round bar 10 long is used as a column. The young's modulus
        is E and the second moment is I. One end of the column is fixed
        and other end is free. A load of 15KN is applied on it.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy import Symbol, symbols
        >>> E, I, P = symbols('E, I, P', positive=True)
        >>> c = Column(10, E, I, 15000, top="free", bottom="fixed")
        >>> c.end_conditions
        {'bottom': 'fixed', 'top': 'free'}
        >>> c.boundary_conditions
        {'deflection': [(0, 0), (10, d)], 'slope': [(0, 0)]}
        >>> c.moment()
        -15000*d + 15000*y(x)
        >>> c.solve_slope_deflection()
        >>> c.deflection()
        -d*cos(50*sqrt(6)*x/(sqrt(E)*sqrt(I))) + d
        >>> c.slope()
        50*sqrt(6)*d*sin(50*sqrt(6)*x/(sqrt(E)*sqrt(I)))/(sqrt(E)*sqrt(I))
        >>> c.critical_load()
        pi**2*E*I/400
        """
        y = Function('y')
        x = self._variable
        P = Symbol('P', positive=True)

        if self._critical_load is None:
            defl_eqs = []
            # taking last two bounndary conditions which are actually
            # the initial boundary conditions.
            for point, value in self._boundary_conditions['deflection'][-2:]:
                defl_eqs.append(self._deflection.subs(x, point) - value)

            # C1, C2 already solved, solve for P
            self._critical_load = solve(defl_eqs, P, dict=True)[0][P]

        return self._critical_load


    def moment(self):
        """Returns the moment equation in terms of any arbitrary point ``x``
        on the column and the deflection at that point.
        """
        P = Symbol('P', positive=True)
        return self._moment.subs(P, self._load)


    def slope(self):
        """Returns the slope equation in terms of any arbitrary point ``x``
        on the column.
        """
        P = Symbol('P', positive=True)
        return self._slope.subs(P, self._load)


    def deflection(self):
        """Returns the deflection equation in terms of any arbitrary point ``x``
        on the column.
        """
        P = Symbol('P', positive=True)
        return self._deflection.subs(P, self._load)
