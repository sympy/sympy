"""
This module can be used to solve 2D beam bending problems with
singularity functions in mechanics.

"""

from __future__ import print_function, division

from sympy.core import S, Symbol, diff
from sympy.solvers import linsolve
from sympy.printing import sstr
from sympy.functions import SingularityFunction
from sympy.core import sympify
from sympy.integrals import integrate


class Beam(object):
    """
    A Beam is a structural element that is capable of withstanding load
    primarily by resisting against bending. Beams are characterized by
    their cross sectional profile(Second moment of area), their length
    and their material.


    Examples
    ========
    There is a beam of length 4 meters. A constant distributed load of 6 Nm/m
    is applied from half of the beam till the end. There are two supports below
    the beam, one at the starting point and another at the ending point of the beam.
    The deflection of the beam at the end is resticted.

    >>> from sympy.physics.continuum_mechanics.beam import Beam
    >>> from sympy import Symbol, Piecewise
    >>> x = Symbol('x')
    >>> E = Symbol('E')
    >>> I = Symbol('I')
    >>> b = Beam(4, E, I)
    >>> b.apply_load(value=-9, order=-1, start=4)
    >>> b.apply_load(value=-3, order=-1, start=0)
    >>> b.apply_load(order=0, value=6, start=2)
    >>> b.bc_deflection = [(4, 0)]
    >>> b.boundary_conditions
    {'deflection': [(4, 0)], 'moment': [], 'slope': []}
    >>> b.load
    -3*SingularityFunction(x, 0, -1) + 6*SingularityFunction(x, 2, 0) - 9*SingularityFunction(x, 4, -1)
    >>> b.shear_force()
    -3*SingularityFunction(x, 0, 0) + 6*SingularityFunction(x, 2, 1) - 9*SingularityFunction(x, 4, 0)
    >>> b.bending_moment()
    3*SingularityFunction(x, 0, 1) - 3*SingularityFunction(x, 2, 2) + 9*SingularityFunction(x, 4, 1)
    >>> b.slope()
    (3*SingularityFunction(x, 0, 2)/2 - SingularityFunction(x, 2, 3) + 9*SingularityFunction(x, 4, 2)/2 - 7)/(E*I)
    >>> b.deflection()
    (-7*x + SingularityFunction(x, 0, 3)/2 - SingularityFunction(x, 2, 4)/4 + 3*SingularityFunction(x, 4, 3)/2)/(E*I)
    >>> b.deflection().rewrite(Piecewise)
    (-7*x + Piecewise((x**3, x > 0), (0, True))/2
          + 3*Piecewise(((x - 4)**3, x - 4 > 0), (0, True))/2
          - Piecewise(((x - 2)**4, x - 2 > 0), (0, True))/4)/(E*I)
    """

    def __init__(self, length, elastic_modulus, second_moment, variable=Symbol('x')):
        """Initializes the class.

        Parameters
        ==========
        length :
            A Symbol or value representing the Beam's length.
        elastic_modulus : Sympifyable
            A SymPy expression representing the Beam's Modulus of Elasticity.
            It is a measure of the stiffness of the Beam material.
        second_moment : Sympifyable
            A SymPy expression representing the Beam's Second moment of area.
            It is a geometrical property of an area which reflects how its
            points are distributed with respect to its neutral axis.
        variable : Symbol
            A Symbol object that will be used as the variable along the beam
            while representing the load, shear, moment, slope and deflection
            curve. By default, it is set to ``Symbol('x')``.
        """
        self._length = sympify(length)
        self._elastic_modulus = sympify(elastic_modulus)
        self._second_moment = sympify(second_moment)
        if isinstance(variable, Symbol):
            self._variable = variable
        else:
            raise TypeError("""The variable should be a Symbol object.""")
        self._boundary_conditions = {'deflection': [], 'moment': [], 'slope': []}
        self._load = 0

    def __str__(self):
        str_sol = 'Beam(%s, %s, %s)' % (sstr(self._length), sstr(self._elastic_modulus), sstr(self._second_moment))
        return str_sol

    __repr__ = __str__

    @property
    def length(self):
        """Length of the Beam."""
        return self._length

    @length.setter
    def length(self, l):
        self._length = sympify(l)

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
    def elastic_modulus(self):
        """Young's Modulus of the Beam. """
        return self._elastic_modulus

    @elastic_modulus.setter
    def elastic_modulus(self, e):
        self._elastic_modulus = sympify(e)

    @property
    def second_moment(self):
        """Second moment of area of the Beam. """
        return self._second_moment

    @second_moment.setter
    def second_moment(self, i):
        self._second_moment = sympify(i)

    @property
    def boundary_conditions(self):
        """
        Returns a dictionary of boundary conditions applied on the beam.
        The dictionary has three kewwords namely moment, slope and deflection.
        The value of each keyword is a list of tuple, where each tuple
        contains loaction and value of a boundary condition in the format
        (location, value).

        Examples
        ========
        There is a beam of length 4 meters. The bending moment at 0 should be 4
        and at 4 it should be 0. The slope of the beam should be 1 at 0. The
        deflection should be 2 at 0.

        >>> from sympy.physics.continuum_mechanics.beam import Beam
        >>> from sympy import Symbol
        >>> E = Symbol('E')
        >>> I = Symbol('I')
        >>> b = Beam(4, E, I)
        >>> b.bc_moment = [(0, 4), (4, 0)]
        >>> b.bc_deflection = [(0, 2)]
        >>> b.bc_slope = [(0, 1)]
        >>> b.boundary_conditions
        {'deflection': [(0, 2)], 'moment': [(0, 4), (4, 0)], 'slope': [(0, 1)]}

        Here the deflection of the beam should be ``2`` at ``0``.
        Similarly, the slope of the beam should be ``1`` at ``0`` and
        there are two boundary conditions for bending moment. The bending
        moment of the beam should be ``4`` at ``0`` and ``0`` at ``4``
        """
        return self._boundary_conditions

    @property
    def bc_moment(self):
        return self._boundary_conditions['moment']

    @bc_moment.setter
    def bc_moment(self, m_bcs):
        self._boundary_conditions['moment'] = m_bcs

    @property
    def bc_slope(self):
        return self._boundary_conditions['slope']

    @bc_slope.setter
    def bc_slope(self, s_bcs):
        self._boundary_conditions['slope'] = s_bcs

    @property
    def bc_deflection(self):
        return self._boundary_conditions['deflection']

    @bc_deflection.setter
    def bc_deflection(self, d_bcs):
        self._boundary_conditions['deflection'] = d_bcs

    def apply_load(self, value, start, order):
        """
        This method adds up the loads given to a particular beam object.

        Parameters
        ==========
        value :
            The magnitude of an applied load.
        start :
            The starting point of the applied load.
        end :
            An optional argument that can be used if the load have an end point
            within the length of the beam.
        order :
            The order of the applied load.

            For Moments, order = -2
            For Pointloads, order = -1
            For constant distributed load, order = 0
            For Ramp load, order = 1
            For Parabolic Ramp, order = 2
            ... so on.

        Examples
        ========
        There is a beam of length 4 meters. A moment of magnitude 3 Nm is
        applied in the clockwise direction at the starting point of the beam.
        A pointload of magnitude 4 N is applied from the top of the beam at
        2 meters from the starting point and a parabolic ramp load of magnitude
        2 Nm/m is applied below the beam starting from 3 meters away from the
        starting point of the beam.

        >>> from sympy.physics.continuum_mechanics.beam import Beam
        >>> from sympy import Symbol
        >>> E = Symbol('E')
        >>> I = Symbol('I')
        >>> b = Beam(4, E, I)
        >>> b.apply_load(value=-3, order=-2, start=0)
        >>> b.apply_load(value=4, order=-1, start=2)
        >>> b.apply_load(value=-2, order=2, start=3)
        >>> b.load
        -3*SingularityFunction(x, 0, -2) + 4*SingularityFunction(x, 2, -1) - 2*SingularityFunction(x, 3, 2)

        """
        x = self.variable
        value = sympify(value)
        start = sympify(start)
        order = sympify(order)

        self._load += value*SingularityFunction(x, start, order)

    @property
    def load(self):
        """
        Returns a Singularity Function expression which represents
        the load distribution curve of the Beam object.

        Examples
        ========
        There is a beam of length 4 meters. A moment of magnitude 3 Nm is
        applied in the clockwise direction at the starting point of the beam.
        A pointload of magnitude 4 N is applied from the top of the beam at
        2 meters from the starting point and a parabolic ramp load of magnitude
        2 Nm/m is applied below the beam starting from 3 meters away from the
        starting point of the beam.

        >>> from sympy.physics.continuum_mechanics.beam import Beam
        >>> from sympy import Symbol
        >>> E = Symbol('E')
        >>> I = Symbol('I')
        >>> b = Beam(4, E, I)
        >>> b.apply_load(value=-3, order=-2, start=0)
        >>> b.apply_load(value=4, order=-1, start=2)
        >>> b.apply_load(value=-2, order=2, start=3)
        >>> b.load
        -3*SingularityFunction(x, 0, -2) + 4*SingularityFunction(x, 2, -1) - 2*SingularityFunction(x, 3, 2)
        """
        return self._load

    def shear_force(self):
        """
        Returns a Singularity Function expression which represents
        the shear force curve of the Beam object.

        Examples
        ========
        There is a beam of length 4 meters. A moment of magnitude 3 Nm is
        applied in the clockwise direction at the starting point of the beam.
        A pointload of magnitude 4 N is applied from the top of the beam at
        2 meters from the starting point and a parabolic ramp load of magnitude
        2 Nm/m is applied below the beam starting from 3 meters away from the
        starting point of the beam. The bending moment at 0 should be 4 and
        at 4 it should be 0. The slope of the beam should be 1 at 0. The
        deflection should be 2 at 0.

        >>> from sympy.physics.continuum_mechanics.beam import Beam
        >>> from sympy import Symbol
        >>> E = Symbol('E')
        >>> I = Symbol('I')
        >>> b = Beam(4, E, I)
        >>> b.apply_load(value=-3, order=-2, start=0)
        >>> b.apply_load(value=4, order=-1, start=2)
        >>> b.apply_load(value=-2, order=2, start=3)
        >>> b.bc_moment = [(0, 4), (4, 0)]
        >>> b.bc_deflection = [(0, 2)]
        >>> b.bc_slope = [(0, 1)]
        >>> b.shear_force()
        -3*SingularityFunction(x, 0, -1) + 4*SingularityFunction(x, 2, 0) - 2*SingularityFunction(x, 3, 3)/3 - 71/24
        """
        x = self.variable
        if not self._boundary_conditions['moment']:
            return integrate(self.load, x)
        return diff(self.bending_moment(), x)

    def bending_moment(self):
        """
        Returns a Singularity Function expression which represents
        the bending moment curve of the Beam object.

        Examples
        ========
        There is a beam of length 4 meters. A moment of magnitude 3 Nm is
        applied in the clockwise direction at the starting point of the beam.
        A pointload of magnitude 4 N is applied from the top of the beam at
        2 meters from the starting point and a parabolic ramp load of magnitude
        2 Nm/m is applied below the beam starting from 3 meters away from the
        starting point of the beam. The bending moment at 0 should be 4 and
        at 4 it should be 0. The slope of the beam should be 1 at 0. The
        deflection should be 2 at 0.

        >>> from sympy.physics.continuum_mechanics.beam import Beam
        >>> from sympy import Symbol
        >>> E = Symbol('E')
        >>> I = Symbol('I')
        >>> b = Beam(4, E, I)
        >>> b.apply_load(value=-3, order=-2, start=0)
        >>> b.apply_load(value=4, order=-1, start=2)
        >>> b.apply_load(value=-2, order=2, start=3)
        >>> b.bc_moment = [(0, 4), (4, 0)]
        >>> b.bc_deflection = [(0, 2)]
        >>> b.bc_slope = [(0, 1)]
        >>> b.bending_moment()
        -71*x/24 - 3*SingularityFunction(x, 0, 0) + 4*SingularityFunction(x, 2, 1) - SingularityFunction(x, 3, 4)/6 + 7
        """
        x = self.variable
        if not self._boundary_conditions['moment']:
            return -1*integrate(self.shear_force(), x)

        C1 = Symbol('C1')
        C2 = Symbol('C2')
        load_curve = self.load
        shear_curve = integrate(load_curve, x) + C1
        moment_curve = integrate(shear_curve, x) + C2

        bc_eqs = []
        for position, value in self._boundary_conditions['moment']:
            eqs = moment_curve.subs(x, position) - value
            bc_eqs.append(eqs)

        constants = list(linsolve(bc_eqs, C1, C2))
        moment_curve = moment_curve.subs({C1: constants[0][0], C2: constants[0][1]})

        return moment_curve

    def slope(self):
        """
        Returns a Singularity Function expression which represents
        the slope the elastic curve of the Beam object.

        Examples
        ========
        There is a beam of length 4 meters. A moment of magnitude 3 Nm is
        applied in the clockwise direction at the starting point of the beam.
        A pointload of magnitude 4 N is applied from the top of the beam at
        2 meters from the starting point and a parabolic ramp load of magnitude
        2 Nm/m is applied below the beam starting from 3 meters away from the
        starting point of the beam. The bending moment at 0 should be 4 and
        at 4 it should be 0. The slope of the beam should be 1 at 0. The
        deflection should be 2 at 0.

        >>> from sympy.physics.continuum_mechanics.beam import Beam
        >>> from sympy import Symbol
        >>> E = Symbol('E')
        >>> I = Symbol('I')
        >>> b = Beam(4, E, I)
        >>> b.apply_load(value=-3, order=-2, start=0)
        >>> b.apply_load(value=4, order=-1, start=2)
        >>> b.apply_load(value=-2, order=2, start=3)
        >>> b.bc_moment = [(0, 4), (4, 0)]
        >>> b.bc_deflection = [(0, 2)]
        >>> b.bc_slope = [(0, 1)]
        >>> b.slope()
        (-71*x**2/48 + 7*x - 3*SingularityFunction(x, 0, 1) + 2*SingularityFunction(x, 2, 2) - SingularityFunction(x, 3, 5)/30 + 1)/(E*I)
        """
        x = self.variable
        E = self.elastic_modulus
        I = self.second_moment
        if not self._boundary_conditions['slope']:
            return diff(self.deflection(), x)

        C3 = Symbol('C3')
        slope_curve = integrate(self.bending_moment(), x) + C3

        bc_eqs = []
        for position, value in self._boundary_conditions['slope']:
            eqs = slope_curve.subs(x, position) - value
            bc_eqs.append(eqs)

        constants = list(linsolve(bc_eqs, C3))
        slope_curve = slope_curve.subs({C3: constants[0][0]})
        return S(1)/(E*I)*slope_curve

    def deflection(self):
        """
        Returns a Singularity Function expression which represents
        the elastic curve or deflection of the Beam object.

        Examples
        ========
        There is a beam of length 4 meters. A moment of magnitude 3 Nm is
        applied in the clockwise direction at the starting point of the beam.
        A pointload of magnitude 4 N is applied from the top of the beam at
        2 meters from the starting point and a parabolic ramp load of magnitude
        2 Nm/m is applied below the beam starting from 3 meters away from the
        starting point of the beam. The bending moment at 0 should be 4 and
        at 4 it should be 0. The slope of the beam should be 1 at 0. The
        deflection should be 2 at 0.

        >>> from sympy.physics.continuum_mechanics.beam import Beam
        >>> from sympy import Symbol
        >>> E = Symbol('E')
        >>> I = Symbol('I')
        >>> b = Beam(4, E, I)
        >>> b.apply_load(value=-3, order=-2, start=0)
        >>> b.apply_load(value=4, order=-1, start=2)
        >>> b.apply_load(value=-2, order=2, start=3)
        >>> b.bc_moment = [(0, 4), (4, 0)]
        >>> b.bc_deflection = [(0, 2)]
        >>> b.bc_slope = [(0, 1)]
        >>> b.deflection()
        (-71*x**3/144 + 7*x**2/2 + x - 3*SingularityFunction(x, 0, 2)/2 + 2*SingularityFunction(x, 2, 3)/3 - SingularityFunction(x, 3, 6)/180 + 2)/(E*I)
        """
        x = self.variable
        E = self.elastic_modulus
        I = self.second_moment
        if not self._boundary_conditions['deflection'] and not self._boundary_conditions['slope']:
            return S(1)/(E*I)*integrate(integrate(self.bending_moment(), x), x)
        elif not self._boundary_conditions['deflection']:
            return integrate(self.slope(), x)
        elif not self._boundary_conditions['slope'] and self._boundary_conditions['deflection']:
            C3 = Symbol('C3')
            slope_curve = integrate(self.bending_moment(), x) + C3
            deflection_curve = integrate(slope_curve, x)
            bc_eqs = []
            for position, value in self._boundary_conditions['deflection']:
                eqs = deflection_curve.subs(x, position) - value
                bc_eqs.append(eqs)
            constants = list(linsolve(bc_eqs, C3))
            deflection_curve = deflection_curve.subs({C3: constants[0][0]})
            return S(1)/(E*I)*deflection_curve

        C4 = Symbol('C4')
        deflection_curve = integrate((E*I)*self.slope(), x) + C4

        bc_eqs = []
        for position, value in self._boundary_conditions['deflection']:
            eqs = deflection_curve.subs(x, position) - value
            bc_eqs.append(eqs)

        constants = list(linsolve(bc_eqs, C4))
        deflection_curve = deflection_curve.subs({C4: constants[0][0]})
        return S(1)/(E*I)*deflection_curve
