"""
This module can be used to solve beam bending problems in mechanics.

"""

from __future__ import print_function, division

from sympy.core import S, Symbol
from sympy.printing import sstr
from sympy.physics.mechanics import Point
from sympy.functions import SingularityFunction
from sympy import sympify


class Beam(object):

    """
    A Beam is a structural element that is capable of withstanding load
    primarily by resisting against bending.
    Beams are characterized by their profile (Second moment of area),
    their length, and their material.

    Parameters
    ==========
    length : Sympifyable
        A SymPy expression representing the Beam's length.
    elastic_modulus : Sympifyable
        A SymPy expression representing the Beam's Modulus of Elasticity.
        It is a measure of the stiffness of the Beam material.
    second_moment : Sympifyable
        A SymPy expression representing the Beam's Second moment of area.
        It is a geometrical property of an area which reflects how its
        points are distributed with regard to an arbitrary axis.

    Examples
    ========
    >>> from sympy.physics.continuum_mechanics.beam import Beam
    >>> from sympy import Symbol
    >>> E = Symbol('E')
    >>> I = Symbol('I')
    >>> Beam(1, E, I)
    Beam(1, E, I)

    """

    def __init__(self, length, elastic_modulus, second_moment):
        self._length = length
        self._elastic_modulus = elastic_modulus
        self._second_moment = second_moment
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

    def apply_boundary_conditions(self, **bcs):
        """
        Takes the boundary conditions of the beam bending problem as input.
        The boundary conditions should be passed as keyworded arguments.
        It is suggested to use ``moment``, ``slope`` and ``deflection`` as the
        keywords while providing the boundary conditions of the moment curve,
        slope curve and the deflection curve respectively as inputs.

        Outputs a dictionary.

        Examples
        ========
        >>> from sympy.physics.continuum_mechanics.beam import Beam
        >>> from sympy import Symbol
        >>> E = Symbol('E')
        >>> I = Symbol('I')
        >>> b = Beam(4, E, I)
        >>> bcs = b.apply_boundary_conditions(moment = [(0, 4), (4, 0)], deflection = [(0, 2)], slope = [(0, 1)])
        >>> bcs
        {'deflection': [(0, 2)], 'moment': [(0, 4), (4, 0)], 'slope': [(0, 1)]}
        >>> bcs['moment']
        [(0, 4), (4, 0)]

        """
        for m_bcs in bcs['moment']:
            self._boundary_conditions['moment'].append(m_bcs)
        for s_bcs in bcs['slope']:
            self._boundary_conditions['slope'].append(s_bcs)
        for d_bcs in bcs['deflection']:
            self._boundary_conditions['deflection'].append(d_bcs)
        return self._boundary_conditions

    def apply_moment_boundary_conditions(self, *m_bcs):
        """
        Takes only the moment boundary conditions as input.

        Examples
        ========
        >>> from sympy.physics.continuum_mechanics.beam import Beam
        >>> from sympy import Symbol
        >>> E = Symbol('E')
        >>> I = Symbol('I')
        >>> b = Beam(4, E, I)
        >>> bcs = b.apply_boundary_conditions(moment = [(0, 4), (4, 0)], deflection = [(0, 2)], slope = [(0, 1)])
        >>> bcs
        {'deflection': [(0, 2)], 'moment': [(0, 4), (4, 0)], 'slope': [(0, 1)]}
        >>> b.apply_moment_boundary_conditions((4, 3), (5, 0))
        [(0, 4), (4, 0), (4, 3), (5, 0)]

        """
        for bcs in m_bcs:
            self._boundary_conditions['moment'].append(bcs)
        return self._boundary_conditions['moment']

    def apply_slope_boundary_conditions(self, *s_bcs):
        """
        Takes only the slope boundary conditions as input.

        Examples
        ========
        >>> from sympy.physics.continuum_mechanics.beam import Beam
        >>> from sympy import Symbol
        >>> E = Symbol('E')
        >>> I = Symbol('I')
        >>> b = Beam(4, E, I)
        >>> bcs = b.apply_boundary_conditions(moment = [(0, 4), (4, 0)], deflection = [(0, 2)], slope = [(0, 1)])
        >>> bcs
        {'deflection': [(0, 2)], 'moment': [(0, 4), (4, 0)], 'slope': [(0, 1)]}
        >>> b.apply_slope_boundary_conditions((4, 3), (5, 0))
        [(0, 1), (4, 3), (5, 0)]

        """
        for bcs in s_bcs:
            self._boundary_conditions['slope'].append(bcs)
        return self._boundary_conditions['slope']

    def apply_deflection_boundary_conditions(self, *d_bcs):
        """
        Takes only the slope boundary conditions as input.

        Examples
        ========
        >>> from sympy.physics.continuum_mechanics.beam import Beam
        >>> from sympy import Symbol
        >>> E = Symbol('E')
        >>> I = Symbol('I')
        >>> b = Beam(4, E, I)
        >>> bcs = b.apply_boundary_conditions(moment = [(0, 4), (4, 0)], deflection = [(0, 2)], slope = [(0, 1)])
        >>> bcs
        {'deflection': [(0, 2)], 'moment': [(0, 4), (4, 0)], 'slope': [(0, 1)]}
        >>> b.apply_deflection_boundary_conditions((4, 3), (5, 0))
        [(0, 2), (4, 3), (5, 0)]

        """
        for bcs in d_bcs:
            self._boundary_conditions['deflection'].append(bcs)
        return self._boundary_conditions['deflection']

    def boundary_conditions(self):
        """
        Returns a dictionary of boundary conditions applied on the beam.
        """
        return self._boundary_conditions

    def _load_as_SingularityFunction(self, load):
        x = Symbol('x')
        if isinstance(load, PointLoad):
            if load.moment:
                return load.value*SingularityFunction(x, load.location, S(-2))
            return load.value*SingularityFunction(x, load.location, S(-1))
        elif isinstance(load, DistributedLoad):
            return load.value*SingularityFunction(x, load.start, load.order)

    def apply_loads(self, *loads):
        """
        Takes PointLoad and DistributedLoad as input. This method would apply
        the loads, that are passed as arguments, to the beam object. Internally
        this method would represent the PointLoads and DistributedLoads using
        Singularity Functions.

        """
        for load in loads:
            if isinstance(load, PointLoad) or isinstance(load, DistributedLoad):
                self._load += self._load_as_SingularityFunction(load)
            else:
                raise TypeError("""apply_loads takes either PointLoad or DistributedLoad objects""")

    def load_distribution(self):
        """
        Returns a Singularity Function expression which represents
        the load distribution curve of the Beam object.

        Examples
        ========
        >>> from sympy.physics.continuum_mechanics.beam import Beam, DistributedLoad, PointLoad
        >>> from sympy import Symbol
        >>> from sympy.physics.mechanics import Point
        >>> E = Symbol('E')
        >>> I = Symbol('I')
        >>> b = Beam(4, E, I)
        >>> P1 = Point('0')
        >>> P2 = Point('2')
        >>> P3 = Point('3')
        >>> Load_1 = PointLoad(location = P1, value = -3, moment = True)
        >>> Load_2 = PointLoad(location = P2, value = 4)
        >>> Load_3 = DistributedLoad(start = P3, order = 2, value = -2)
        >>> b.apply_loads(Load_1, Load_2, Load_3)
        >>> b.load_distribution()
        -3*SingularityFunction(x, 0, -2) + 4*SingularityFunction(x, 2, -1) - 2*SingularityFunction(x, 3, 2)

        """
        return self._load

    def shear_force(self):
        """
        Returns a Singularity Function expression which represents
        the shear force curve of the Beam object.

        Examples
        ========

        """

    def bending_moment(self):
        """
        Returns a Singularity Function expression which represents
        the bending moment curve of the Beam object.

        Examples
        ========

        """

    def slope(self):
        """
        Returns a Singularity Function expression which represents
        the slope the elastic curve of the Beam object.

        Examples
        ========

        """

    def deflection(self):
        """
        Returns a Singularity Function expression which represents
        the elastic curve or deflection of the Beam object.

        Examples
        ========

        """


class PointLoad(object):
    """A Point Load.

    A load applied to a single, specific point. It is also known as a
    concentrated load.

    Parameters
    ==========
    location : Point
        The specific point of the applied load.
    value : Sympifyable
        A SymPy expression representing the value of the applied load.

    Examples
    ========
    >>> from sympy.physics.continuum_mechanics.beam import PointLoad
    >>> from sympy.physics.mechanics import Point
    >>> p = Point('4')
    >>> PointLoad(location = p, value = -4)
    PointLoad(4, -4, Load)

    A Moment can be defined just by passing moment=True as an argument.

    >>> PointLoad(location = p, value = -4, moment=True)
    PointLoad(4, -4, Moment)

    """

    def __init__(self, location, value, moment=False):
        self._location = location
        self._value = value
        self._moment = moment

    def __str__(self):
        if self.moment:
            str_sol = 'PointLoad(%s, %s, %s)' % (sstr(self.location), sstr(self.value), sstr('Moment'))
            return str_sol
        str_sol = 'PointLoad(%s, %s, %s)' % (sstr(self.location), sstr(self.value), sstr('Load'))
        return str_sol

    __repr__ = __str__

    @property
    def location(self):
        """Location of the applied Point Load."""
        return self._location

    @location.setter
    def location(self, l):
        if not isinstance(l, Point):
            raise TypeError("PointLoad location attribute must be a Point object.")
        self._location = l

    @property
    def value(self):
        """Value of the applied Point Load. """
        return self._value

    @value.setter
    def value(self, v):
        self._value = sympify(v)

    @property
    def moment(self):
        """Stores whether a Point Load is a load or a couple. """
        return self._moment

    @moment.setter
    def moment(self, m):
        if not isinstance(m, bool):
            raise TypeError("PointLoad moment attribute must be a bool object.")
        self._moment = m


class DistributedLoad(object):

    """A Distributed Load.

    A load applied across a length instead of at one point.

    Parameters
    ==========
    start : Point
        The starting point of the applied load.
    order : Sympifyable
        The order of the applied load.
    value : Sympifyable
        A SymPy expression representing the value of the applied load.

    Examples
    ========
    >>> from sympy.physics.continuum_mechanics.beam import DistributedLoad
    >>> from sympy.physics.mechanics import Point
    >>> from sympy import Symbol
    >>> a = Point('4')
    >>> b = 2
    >>> DistributedLoad(start = a, order = b, value = 2)
    DistributedLoad(4, 2, 2)

    """

    def __init__(self, start, order, value):

        if isinstance(start, Point):
            self._start = start
        else:
            raise TypeError("DistributedLoad start attribute must be a Point object.")
        self._order = order
        self._value = value

    def __str__(self):
        str_sol = 'DistributedLoad(%s, %s, %s)' % (sstr(self.start), sstr(self.order), sstr(self.value))
        return str_sol

    __repr__ = __str__

    @property
    def start(self):
        """The starting point of the applied load."""
        return self._start

    @start.setter
    def start(self, s):
        if not isinstance(s, Point):
            raise TypeError("DistributedLoad start attribute must be a Point object.")
        self._start = s

    @property
    def order(self):
        """The order of the applied load."""
        return self._order

    @order.setter
    def order(self, o):
        self._order = sympify(o)

    @property
    def value(self):
        """The value of the applied load."""
        return self._value

    @value.setter
    def value(self, v):
        self._value = sympify(v)
