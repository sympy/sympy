"""
This module can be used to solve column displacement problems
using singularity functions in mechanics.
"""

from sympy import nsimplify, simplify
from sympy.core import Symbol, symbols
from sympy.core.relational import Eq
from sympy.core.sympify import sympify
from sympy.functions import SingularityFunction, factorial
from sympy.integrals import integrate
from sympy.plotting import plot
from sympy.printing import sstr
from sympy.series import limit
from sympy.solvers import linsolve

class Column:
    """
    A column is a structural element that withstands axial
    loading primarily through compression. Columns are
    characterized by their length, material and cress-sectional
    area.

    .. note::
        A consistent sign convention must be used when solving
        these problems. Applied forces are positive when aimed in
        the positive x-direction. Normal forces are positive when
        they lead to extension, and negative when they lead to
        compression.

    .. note::
        The columns are set up horizontally, from left to right.
        This is due to it then having better compatability with
        the 2-Dimensional module, where all objects are projected
        horizontally.

    Examples
    ========
    The is a column with a length of 20 meters. It has an area of
    0.75 m^2 and an elastic modulus of 20.000 kN/m^2. The column is
    fixed at both ends. From x = 0 to x = 10 there is a constant
    distributed load of 5 kN/m. At x = 12 and x = 16 there are point
    loads of 20 kN. Loads aiming to the right are positive, and positive
    axial forces lead to extension.

    >>> from sympy.physics.continuum_mechanics.column import Column
    >>> c = Column(20, 20000, 0.75)
    >>> c.apply_support(0)
    >>> c.apply_support(20)
    >>> c.apply_load(5, 0, 0, end=10)
    >>> c.apply_load(20, 12, -1)
    >>> c.apply_load(20, 16, -1)
    >>> c.applied_loads
    [(5, 0, 0, 10), (20, 12, -1, None), (20, 16, -1, None)]
    >>> c.load
    R_0*SingularityFunction(x, 0, -1) + R_20*SingularityFunction(x, 20, -1)
    + 5*SingularityFunction(x, 0, 0) - 5*SingularityFunction(x, 10, 0)
    + 20*SingularityFunction(x, 12, -1) + 20*SingularityFunction(x, 16, -1)
    >>> c.solve_for_reaction_loads()
    >>> c.reaction_loads
    {R_0: -99/2, R_20: -81/2}
    >>> c.axial_force()
    99*SingularityFunction(x, 0, 0)/2 - 5*SingularityFunction(x, 0, 1)
    + 5*SingularityFunction(x, 10, 1) - 20*SingularityFunction(x, 12, 0)
    - 20*SingularityFunction(x, 16, 0) + 81*SingularityFunction(x, 20, 0)/2
    >>> c.extension()
    0.0033*SingularityFunction(x, 0, 1) - 0.000166666666666667*SingularityFunction(x, 0, 2)
    + 0.000166666666666667*SingularityFunction(x, 10, 2) - 0.00133333333333333*SingularityFunction(x, 12, 1)
    - 0.00133333333333333*SingularityFunction(x, 16, 1) + 0.0027*SingularityFunction(x, 20, 1)
    """

    def __init__(self, length, elastic_modulus, area, variable=Symbol('x'), base_char='C'):
        """
        Initializes the Column class.

        Parameters
        ==========

        length: Sympifyable
            A Symbol or value representing the column's length.

        elastic_modulus: Sympifyable
            A Symbol or value representing the column's modulus of
            elasticity. It is a measure of the stiffness of the material.

        area: Sympifyable
            A symbol or value representing the column's cross-sectional
            area.

        variable : Symbol, optional
            A Symbol object that will be used as the variable along the column
            while representing the load, axial force, or axial deformation.
            By default, it is set to ``Symbol('x')``.

        base_char : String, optional
            A String that will be used as base character to generate sequential
            symbols for integration constants in cases where boundary conditions
            are not sufficient to solve them.
        """
        self.length = length
        self.elastic_modulus = elastic_modulus
        self.area = area
        self.variable = variable
        self._base_char = base_char
        self._bc_extension = []
        self._bc_hinge = []
        self._applied_supports = []
        self._applied_hinges = []
        self._applied_loads = []

        self._reaction_loads = {}
        self._hinge_extensions = {}
        self._integration_constants = {}

        self._load = 0
        self._axial_force = 0
        self._extension = 0

        self._is_solved = False

    def __str__(self):
        str_sol = 'Column({}, {}, {})'.format(sstr(self._length), sstr(self._elastic_modulus), sstr(self._area))
        return str_sol

    @property
    def length(self):
        """Returns the length of the column."""
        return self._length

    @length.setter
    def length(self, l):
        self._length = sympify(l)

    @property
    def elastic_modulus(self):
        """Returns the elastic modulus of the column"""
        return self._elastic_modulus

    @elastic_modulus.setter
    def elastic_modulus(self, E):
        self._elastic_modulus = sympify(E)

    @property
    def area(self):
        """Returns the cross-sectional area of the column."""
        return self._area

    @area.setter
    def area(self, A):
        self._area = sympify(A)

    @property
    def variable(self):
        return self._variable

    @variable.setter
    def variable(self, v):
        if isinstance(v, Symbol):
            self._variable = v
        else:
            raise TypeError("""The variable should be a Symbol object.""")

    @property
    def reaction_loads(self):
        """Returns the reactions loads as dictionary."""
        return self._reaction_loads

    @property
    def hinge_extensions(self):
        """Returns the hinge extensions as dictionary."""
        return self._hinge_extensions

    def apply_support(self, loc):
        """
        Add supports to the system.

        Supports can be defined in two ways in the module:

        1. Using ``apply_support()``
        Automatically applies all required boundary conditions internally and
        generates a ``Symbol(R_loc)`` representing the reaction load at the
        specified support location.

        2. Manual method
        Add a reaction symbol as a load using ``apply_load()`` and manually
        specify the corresponding boundary conditions.


        This applies a support that is fixed in the
        horizontal direction. using the apply_support() method

        Parameters
        ==========
        loc: Sympifyable
            Location at which the fixed support is applied.

        Examples
        ========
        There is a column of 10 meters. It has an area A and
        elastic modulus E. The column has fixed supports at
        both ends.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy.core.symbol import symbols
        >>> E, A = symbols('E A ')
        >>> c = Column(10, E, A)
        >>> c.apply_support(0)
        >>> c.apply_support(10)
        >>> c.apply_load(-10, 10, -1)
        >>> print(c.load)
        R_0*SingularityFunction(x, 0, -1) + R_10*SingularityFunction(x, 10, -1)
            - 10*SingularityFunction(x, 10, -1)


        This applies a support that is fixed in the
        horizontal direction. using the apply_load() method

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy.core.symbol import symbols
        >>> E, A = symbols('E A ')
        >>> R_0, R_10 = symbols('R_0 R_10')
        >>> c = Column(10, E, A)
        >>> c.apply_load(R_0,0,-1)   # Reaction applied as a point load
        >>> c.apply_load(R_10,10,-1) # Reaction applied as a poinnt load
        >>> c._bc_extension.append(0) # boundary conditions are added manually
        >>> c._bc_extension.append(10) # boundary conditions are added manually
        >>> c.apply_load(-10, 10, -1)
        >>> print(c.load)
        R_0*SingularityFunction(x, 0, -1) + R_10*SingularityFunction(x, 10, -1)
            - 10*SingularityFunction(x, 10, -1)
        """
        loc = sympify(loc)
        if loc in self._applied_supports:
            raise ValueError("There is already a support at this location.")

        if loc.is_number:
            reaction_load = Symbol(f'R_{float(loc):g}')
        else:
            reaction_load = Symbol(f'R_{str(loc)}')

        self._applied_supports.append(reaction_load)
        self._load += reaction_load * SingularityFunction(self.variable, loc, -1)
        self._bc_extension.append(loc)

    def apply_load(self, value, start, order, end=None):
        """
        This method applies a load to the column object. Only
        axial loads can be applied, no moments.

        Parameters
        ==========
        value: Sympifyable
            The value inserted should have the units [Force/(Distance**(n+1)]
            where n is the order of applied load.
            Units for applied loads:

                - For point loads: kN*m
                - For constant distributed load: kN/m
                - For ramp loads: kN/m**2
                - For parabolic ramp loads: kN/m**3
                - And so on.

        loc: Sympifyable
            The starting point of the applied load. For
            point loads this is simply the location.
        order: Integer
            The order of the singularity function of the applied load.

                - For point loads, order =-1
                - For constant distributed load, order = 0
                - For ramp loads, order = 1
                - For parabolic ramp loads, order = 2
                - ... so on.

        Examples
        ========
        There is a column of 10 meters, area A and elastic modulus E. It
        is supported only at the left side. There is a negative point load
        of 10 kN applied to the column at the right end. A positive point
        load of 5 kN is applied to the left end.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy.core.symbol import symbols
        >>> E, A = symbols('E A')
        >>> c = Column(10, E, A)
        >>> c.apply_support(0)
        >>> c.apply_load(-10, 10, -1)
        >>> c.applied_loads
        [(-10, 10, -1, None)]
        >>> c.apply_load(5, 0, -1)
        >>> c.applied_loads
        [(-10, 10, -1, None), (5, 0, -1, None)]
        """
        value = sympify(value)
        start = sympify(start)
        order = sympify(order)

        if end is not None:
            end = sympify(end)
            self._applied_loads.append((value, start, order, end))
            self._load += value * SingularityFunction(self.variable, start, order)

            if simplify(start - end).is_positive:
                self._load += self._taylor_helper(self.variable, -value, start, order, end)

            else:
                self._load -= self._taylor_helper(self.variable, value, start, order, end)

        else:
            self._applied_loads.append((value, start, order, None))
            self._load += value * SingularityFunction(self.variable, start, order)

    def remove_load(self, value, start, order, end=None):
        """
        Removes an applied load.

        Examples
        ========
        There is a column with a length of 4 meters, area A and elastic
        modulus E. A point load of -2 kN is applied at x = 2 and a
        constant distributed load of 2 kN/m is applied along the whole
        column. Later, the point load is removed.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy.core.symbol import symbols
        >>> E, A = symbols('E A')
        >>> c = Column(4, E, A)
        >>> c.apply_load(-2, 2, -1)
        >>> c.apply_load(2, 0, 0, end=4)
        >>> print(c.load)
        2*SingularityFunction(x, 0, 0) - 2*SingularityFunction(x, 2, -1) - 2*SingularityFunction(x, 4, 0)
        >>> c.remove_load(-2, 2, -1)
        >>> print(c.load)
        2*SingularityFunction(x, 0, 0) - 2*SingularityFunction(x, 4, 0)
        """
        value = sympify(value)
        start = sympify(start)
        order = sympify(order)

        if (value, start, order, end) in self._applied_loads:
            self._load -= value*SingularityFunction(self.variable, start, order)
            self._applied_loads.remove((value, start, order, end))
        else:
            msg = "No such load distribution exists on the column object."
            raise ValueError(msg)

        if end is not None:
            if simplify(start - end).is_positive:
                self._load -= self._taylor_helper(self.variable, -value, start, order, end)
            else:
                self._load += self._taylor_helper(self.variable, value, start, order, end)

    def _taylor_helper(self, x, value, start, order, end):
        """
        This helper function provides a Taylor series that
        is used to define the summation of singularity
        that is subtracted from the load equation at 'end'.
        This helper is only called for loads with a provided
        'end' value, their order has to be above 0.
        """
        if order.is_negative:
            msg = ("If 'end' is provided the 'order' of the load cannot "
                    "be negative.")
            raise ValueError(msg)

        f = value * (x - start)**order # Symbolic function

        end_load = 0
        for i in range(0, order + 1):
            end_load += (f.diff(x, i).subs(x, end) *
                            SingularityFunction(x, end, i)/factorial(i))
        return end_load

    @property
    def applied_loads(self):
        """
        Returns a list of all loads applied to the column. Each
        load in the list is a tuple of form (value, start, order).

        Examples
        ========
        There is a column of length L, area A and elastic modulus E. It
        is supported at both ends of the column and in the middle. There
        is a positive force F applied to the column at 1/4 * L and
        3/4 * L.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy.core.symbol import symbols
        >>> L, E, A, F = symbols('L E A F')
        >>> c = Column(L, E, A)
        >>> c.apply_support(0)
        >>> c.apply_support(L / 2)
        >>> c.apply_support(L)
        >>> c.apply_load(F, L / 4, -1)
        >>> c.apply_load(F, 3*L / 4, -1)
        >>> c.applied_loads
        [(F, L/4, -1, None), (F, 3*L/4, -1, None)]
        """
        return self._applied_loads

    def apply_telescope_hinge(self, loc):
        """
        Applies a telescope hinge at the given location. At this location,
        the column is free to move in the axial direction.

        Parameters
        ==========
        loc: Sympifyable
            Location at which the telescope hinge is applied.

        Examples
        ========
        There is a column with a length of 10 meters, elastic
        modulus E and area A. At x = 5 the column is loaded axially
        by a point load of 10 kN. At x = 8, there is a telescope hinge.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> c = Column(10, 20000, 0.5)
        >>> c.apply_support(0)
        >>> c.apply_support(10)
        >>> c.apply_telescope_hinge(8)
        >>> c.apply_load(10, 5, -1)
        >>> c.solve_for_reaction_loads()
        >>> c.reaction_loads
        {R_0: -10, R_10: 0}
        >>> c.hinge_extensions
        {u_8: 1/200}
        """
        E = self.elastic_modulus
        A = self.area

        loc = sympify(loc)
        self._bc_hinge.append(loc)

        if loc.is_number:
            extension = Symbol(f'u_{float(loc):g}')
        else:
            extension = Symbol(f'u_{str(loc)}')

        self._applied_hinges.append(extension)
        self._load += (E * A * extension) * SingularityFunction(self.variable, loc, -2)

    @property
    def load(self):
        """
        Returns the load equation qx of the applied loads and reaction
        forces using singularity functions.

        Examples
        ========
        There is a column of length L, area A and elastic modulus E. It
        is supported at both ends of the column and in the middle. There
        is a positive force F applied to the column at 1/4 * L and
        3/5 * L.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy.core.symbol import symbols
        >>> L, E, A, F = symbols('L E A F')
        >>> c = Column(L, E, A)
        >>> c.apply_support(0)
        >>> c.apply_support(L / 2)
        >>> c.apply_support(L)
        >>> c.apply_load(F, L / 4, -1)
        >>> c.apply_load(F, 3*L / 4, -1)
        >>> c.load
        F*SingularityFunction(x, L/4, -1) + F*SingularityFunction(x, 3*L/4, -1)
            + R_0*SingularityFunction(x, 0, -1) + R_L*SingularityFunction(x, L, -1)
            + R_L/2*SingularityFunction(x, L/2, -1)
        """
        return self._load

    def solve_for_reaction_loads(self,*reactions):
        """
        This method solves the horizontal reaction loads and unknown
        displacement jumps due to telescope hinges.

        Parameters
        ==========
        reactions: If supports of the system are applied manually using apply_load()
        method The reaction load symbols to be passed as reactions.

        Examples
        ========
        A column of 10 meters long, with area A and elastic modulus E
        is supported at both ends. A distributed load of -5 kN/m is
        applied to the whole column. A point load of -10 kN is applied
        at x = 4.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy.core.symbol import symbols
        >>> E, A = symbols('E A')
        >>> c = Column(10, A, E)
        >>> c.apply_support(0)
        >>> c.apply_support(10)
        >>> c.apply_load(-5, 0, 0, end=10)
        >>> c.apply_load(-10, 4, -1)
        >>> c.solve_for_reaction_loads()
        >>> c.reaction_loads
        {R_0: 31, R_10: 29}

        The same example when applied applied supports using manual method
        apply_load()

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy.core.symbol import symbols
        >>> E, A = symbols('E A')
        >>> R, G = symbols('R G')
        >>> c = Column(10, A, E)
        >>> c.apply_load(R,0,-1) # applied reaction load as a point load
        >>> c.apply_load(G,10,-1) # applied reaction load as a point load
        >>> c._bc_extension.append(0) # boundary conditions are added manually
        >>> c._bc_extension.append(10) # boundary conditions are added manually
        >>> c.apply_load(-5, 0, 0, end=10)
        >>> c.apply_load(-10, 4, -1)
        >>> c.solve_for_reaction_loads(R,G) # Reaction loads need to be passed.
        >>> c.reaction_loads
        {G: 29, R: 31}
        """
        x = self.variable
        qx = self._load
        L = self.length
        A = self.area
        E = self.elastic_modulus
        applied_supports = list(reactions) + self._applied_supports

        C_N, C_u = symbols('C_N, C_u')

        axial_force = -integrate(qx, x) + C_N
        extension = -(integrate(integrate(qx, x), x)) / (E * A) + C_N * x + C_u

        eq_axial_force = [
            limit(axial_force, x, 0, dir='-'),
            limit(axial_force, x, L, dir='+')
        ]

        eq_bc_displacement = [Eq(extension.subs(x, loc), 0) for loc in self._bc_extension]

        eq_bc_hinge = [Eq(limit(axial_force, x, loc, dir='+'), 0) for loc in self._bc_hinge] # Just right to avoid infinity

        total_eq = eq_axial_force + eq_bc_displacement + eq_bc_hinge
        unknowns = applied_supports + self._applied_hinges + [C_N, C_u]

        solution = list((linsolve(total_eq, unknowns).args)[0])
        solution = [nsimplify(s, rational=True) for s in solution] # To get rid of tiny residuals

        num_supports = len(applied_supports)
        reaction_solutions = solution[:num_supports]
        self._reaction_loads = dict(zip(applied_supports, reaction_solutions))

        displacement_solutions = solution[num_supports:-2]
        self._hinge_extensions = dict(zip(self._applied_hinges, displacement_solutions))

        integration_constants = solution[-2:]
        self._integration_constants = dict(zip([C_N, C_u], integration_constants))

        self._is_solved = True

    def _solved_load(self):
        """
        Helper function that fills in the solved unknowns into
        the load equation.
        """
        solved_load = self._load
        solved_load = solved_load.subs(self._reaction_loads)
        solved_load = solved_load.subs(self._hinge_extensions)
        return solved_load

    def axial_force(self):
        """
        Returns a singularity function expression that
        represents the axial forces in the column.

        Examples
        ========
        A column with a length of 10 meters, elastic modulus E and
        area A is supported at both ends and in the middle. At x = 3
        the column is loaded by a point load of magnitude -1 kN. Between
        x = 5 and x = 10 the column is loaded by a uniform distributed
        load of -1 kN/m.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy.core.symbol import symbols
        >>> E, A = symbols('E A')
        >>> c = Column(10, E, A)
        >>> c.apply_support(0)
        >>> c.apply_support(5)
        >>> c.apply_support(10)
        >>> c.apply_load(-5, 3, -1)
        >>> c.apply_load(-1, 5, 0, end=10)
        >>> c.axial_force()
        C_N - R_0*SingularityFunction(x, 0, 0) - R_10*SingularityFunction(x, 10, 0)
            - R_5*SingularityFunction(x, 5, 0) + 5*SingularityFunction(x, 3, 0)
            + SingularityFunction(x, 5, 1) - SingularityFunction(x, 10, 1)
        >>> c.solve_for_reaction_loads()
        >>> c.axial_force()
        -2*SingularityFunction(x, 0, 0) + 5*SingularityFunction(x, 3, 0)
            - 11*SingularityFunction(x, 5, 0)/2 + SingularityFunction(x, 5, 1)
            - 5*SingularityFunction(x, 10, 0)/2 - SingularityFunction(x, 10, 1)
        """
        load_equation = self._solved_load()
        x = self.variable
        C_N = self._integration_constants[Symbol('C_N')] if self._is_solved else Symbol('C_N')
        return -integrate(load_equation, x) + C_N

    def extension(self):
        """
        Returns a singularity function expression that
        represents the extension of the column.

        Examples
        ========
        A column with a length of 10 meters has an elastic
        modulus of 210000 kN/m^2 and an area of 1 m^2. It
        is supported at both ends and is loaded by a point
        load of 10 kN at x = 5.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> c = Column(10, 210000, 1)
        >>> c.apply_support(0)
        >>> c.apply_support(10)
        >>> c.apply_load(10, 5, -1)
        >>> c.extension()
        C_N*x + C_u - R_0*SingularityFunction(x, 0, 1)/210000
        - R_10*SingularityFunction(x, 10, 1)/210000 - SingularityFunction(x, 5, 1)/21000
        >>> c.solve_for_reaction_loads()
        >>> c.extension()
        SingularityFunction(x, 0, 1)/42000 - SingularityFunction(x, 5, 1)/21000
        + SingularityFunction(x, 10, 1)/42000
        """
        load_equation = self._solved_load()
        x = self.variable
        E = self.elastic_modulus
        A = self.area
        C_N = self._integration_constants[Symbol('C_N')] if self._is_solved else Symbol('C_N')
        C_u = self._integration_constants[Symbol('C_u')] if self._is_solved else Symbol('C_u')
        return -(integrate(integrate(load_equation, x), x)) / (E * A) + C_N * x + C_u

    def plot_axial_force(self):
        """
        Returns a plot for the axial forces in the column.

        Examples
        ========
        There is a column with a length of 10 meters, an elastic
        modulus of 210000 kN/m^2 and an area of 1 m^2. It is supported
        at x = 0, x = 8 and x = 10. There is a uniform distributed load
        of 5 kN applied to the whole column. Furthermore, there are
        point loads of -10 kN and 5 kN at x = 5 and x = 8, respectively.

        Negative axial force loads to compression, and positive axial
        force leads to extension.

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.column import Column
            >>> c = Column(10, 210000, 1)
            >>> c.apply_support(0)
            >>> c.apply_support(8)
            >>> c.apply_support(10)
            >>> c.apply_load(5, 0, 0, end=10)
            >>> c.apply_load(-10, 5, -1)
            >>> c.apply_load(5, 8, -1)
            >>> c.solve_for_reaction_loads()
            >>> c.plot_axial_force()
            Plot object containing:
            [0]: cartesian line: 65*SingularityFunction(x, 0, 0)/4
            - 5*SingularityFunction(x, 0, 1) + 10*SingularityFunction(x, 5, 0)
            + 75*SingularityFunction(x, 8, 0)/4 + 5*SingularityFunction(x, 10, 0)
            + 5*SingularityFunction(x, 10, 1) for x over (0.0, 10.0)
        """
        return plot(self.axial_force(), (self.variable, 0, self.length), title='Axial Force',
                xlabel=r'$\mathrm{x}$', ylabel=r'$\mathrm{N(x)}$', line_color='r')

    def plot_extension(self):
        """
        Returns a plot for the extensions in the column. To plot
        the extension, numeric values for elastic modulus E and
        are A should be provided.

        Examples
        ========
        There is a column with a length of 10 meters, an elastic
        modulus of 210000 kN/m^2 and an area of 1 m^2. It is supported
        at x = 0, x = 8 and x = 10. There is a uniform distributed load
        of 5 kN applied to the whole column. Furthermore, there are
        point loads of -10 kN and 5 kN at x = 5 and x = 8, respectively.

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.column import Column
            >>> c = Column(10, 210000, 1)
            >>> c.apply_support(0)
            >>> c.apply_support(8)
            >>> c.apply_support(10)
            >>> c.apply_load(5, 0, 0, end=10)
            >>> c.apply_load(-10, 5, -1)
            >>> c.apply_load(5, 8, -1)
            >>> c.solve_for_reaction_loads()
            >>> c.plot_extension()
            Plot object containing:
            [0]: cartesian line: 13*SingularityFunction(x, 0, 1)/168000
            - SingularityFunction(x, 0, 2)/84000 + SingularityFunction(x, 5, 1)/21000
            + SingularityFunction(x, 8, 1)/11200 + SingularityFunction(x, 10, 1)/42000
            + SingularityFunction(x, 10, 2)/84000 for x over (0.0, 10.0)
        """
        return plot(self.extension(), (self.variable, 0, self.length), title='Extension',
                xlabel=r'$\mathrm{x}$', ylabel=r'$\mathrm{u(x)}$', line_color='b')
