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
    Beams are characterized by their profile (moment of Inertia), their length,
    and their material.

    Parameters
    ==========
    length : Sympifyable
        A SymPy expression representing the Beam's length.
    E : Sympifyable
        A SymPy expression representing the Beam's Young's Modulus.
        It is a measure of the stiffness of the Beam material.
    I : Sympifyable
        A SymPy expression representing the Beam's Moment of Inertia.
        It determines the torque needed for a desired angular acceleration
        about a rotational axis of the Beam.

    Examples
    ========
    >>> from sympy.physics.mechanics.beam import Beam
    >>> from sympy import Symbol
    >>> E = Symbol('E')
    >>> I = Symbol('I')
    >>> Beam(1, E, I)
    Beam(1, E, I)

    """

    def __init__(self, length, E, I):
        self._length = length
        self._E = E
        self._I = I

    def __str__(self):
        str_sol = 'Beam(%s, %s, %s)' % (sstr(self._length), sstr(self._E), sstr(self._I))
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
    def E(self):
        """Young's Modulus of the Beam. """
        return self._E

    @E.setter
    def E(self, e):
        self._E = sympify(e)

    @property
    def I(self):
        """Moment of Inertia of the Beam. """
        return self._I

    @I.setter
    def I(self, e):
        self._I = sympify(e)


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
