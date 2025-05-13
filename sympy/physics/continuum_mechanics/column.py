"""
This module can be used to solve column displacement problems
using singularity functions in mechanics.
"""

from sympy.core import Symbol
from sympy.core.sympify import sympify
from sympy.printing import sstr
from sympy.functions import SingularityFunction

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
        the 2-Dimensional module, where all beams are projected
        horizontally.

    Examples
    ========
    """

    def __init__(self, length, elastic_modulus, area, variable=Symbol('x'), base_char='C'):
        """
        Initializes the Column class.

        Parameters
        ==========

        length: Sympifyable
            A Symbol or value representing thecColumn's length.
        
        elastic_modulus: Sympifyable
            A Symbol or value representing the column's modulus of
            elasticity. It is a measure of the stiffness of the material.
            (Make function of the length??)
        
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
        self._bc_deflection = []
        self._applied_supports = []
        self._applied_loads = []
        self._load = 0
    
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
    
    def apply_support(self, loc):
        """
        This method applies a support that is fixed in the
        horizontal direction. It returns the name of the
        unknown reaction load.

        Parameters
        ==========
        loc: Sympifyable
            Location at which the fixed support is applied.
        
        Returns
        =======
        Symbol
            The unknown reaction load as a symbol.

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
        """
        loc = sympify(loc)
        if loc in self._applied_supports:
            return ValueError("There is already a support at this location.")
        
        self._applied_supports.append(loc)
        reaction_load = Symbol('R_' + str(loc))
        self._load += reaction_load * SingularityFunction(self.variable, loc, -1)
        self._bc_deflection.append((loc, 0))

        return reaction_load
    
    def apply_load(self, value, start, order):
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
        >>> E, A = symbols('E A ')
        >>> c = Column(10, E, A)
        >>> c.apply_support(0)
        >>> c.apply_load(-10, 10, -1)
        >>> print(c.applied_loads)
        [(-10, 10, -1)]
        >>> c.apply_load(5, 0, -1)
        >>> print(c.applied_loads)
        [(-10, 10, -1), (5, 0, -1)]
        """
        value = sympify(value)
        start = sympify(start)
        order = sympify(order)

        self._applied_loads.append((value, start, order))
        self._load += value*SingularityFunction(self.variable, start, order)
    
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
        >>> print(c.applied_loads)
        [(F, L/4, -1), (F, 3*L/4, -1)]
        """
        return self._applied_loads
    
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
        >>> print(c.load)
        F*SingularityFunction(x, L/4, -1) + F*SingularityFunction(x, 3*L/4, -1)
            + R_0*SingularityFunction(x, 0, -1) + R_L*SingularityFunction(x, L, -1)
            + R_L/2*SingularityFunction(x, L/2, -1)
        """
        return self._load

    


    

