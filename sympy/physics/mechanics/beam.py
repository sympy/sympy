"""
This module can be used to solve beam bending problems in mechanics.

"""

from __future__ import print_function, division

from sympy.printing import sstr
from sympy.physics.mechanics import Point
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
    >>> from sympy.physics.mechanics.beam import Beam
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
        It is suggested to use ``moment``, ``slope`` and ``deflection`` as the keywords
        while providing the boundary conditions of the moment curve, slope curve and the
        deflection curve respectively as inputs.

        Outputs a dictionary.

        Examples
        ========
        >>> from sympy.physics.mechanics.beam import Beam
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
        self._boundary_conditions = bcs
        return self._boundary_conditions

    def apply_moment_boundary_conditions(self, *m_bcs):
        """
        Takes only the moment boundary conditions as input.

        Examples
        ========
        >>> from sympy.physics.mechanics.beam import Beam
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


    def apply_slope_boundary_conditions(self):
        """
        Takes only the slope boundary conditions as input.

        Examples
        ========

        """

    def apply_deflection_boundary_conditions(self):
        """
        Takes only the slope boundary conditions as input.

        Examples
        ========

        """


    def apply_loads(self, *loads):
        """
        Takes PointLoad and DistributedLoad as input. This method would apply
        the loads, that are passed as arguments, to the beam object. Internally
        this method would represent the PointLoads and DistributedLoads using
        Singularity Functions.

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
    >>> from sympy.physics.mechanics.beam import PointLoad
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
        The starting point of the applied load
    end : Point
        The ending point of the applied load
    value : Sympifyable
        A SymPy expression representing the value of the applied load.

    Examples
    ========
    >>> from sympy.physics.mechanics.beam import DistributedLoad
    >>> from sympy.physics.mechanics import Point
    >>> from sympy import Symbol
    >>> a = Point('4')
    >>> b = Point('6')
    >>> x = Symbol('x')
    >>> DistributedLoad(start = a, end = b, value = 2*x)
    DistributedLoad(4, 6, 2*x)

    """

    def __init__(self, start, end, value):
        self._start = start
        self._end = end
        self._value = value

    def __str__(self):
        str_sol = 'DistributedLoad(%s, %s, %s)' % (sstr(self.start), sstr(self.end), sstr(self.value))
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
    def end(self):
        """The ending point of the applied load."""
        return self._end

    @end.setter
    def end(self, e):
        if not isinstance(e, Point):
            raise TypeError("DistributedLoad end attribute must be a Point object.")
        self._end = e

    @property
    def value(self):
        """The value of the applied load."""
        return self._value

    @value.setter
    def value(self, v):
        self._value = sympify(v)
