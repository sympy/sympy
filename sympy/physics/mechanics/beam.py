"""
This module can be used to solve beam bending problems in mechanics.

"""

from __future__ import print_function, division

from sympy.printing import sstr


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
        self.length = length
        self.E = E
        self.I = I

    def __str__(self):
        str_sol = 'Beam(%s, %s, %s)' % (sstr(self.length), sstr(self.E), sstr(self.I))
        return str_sol

    __repr__ = __str__


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
        self.location = location
        self.value = value
        self.moment = moment

    def __str__(self):
        if self.moment:
            str_sol = 'PointLoad(%s, %s, %s)' % (sstr(self.location), sstr(self.value), sstr('Moment'))
            return str_sol
        str_sol = 'PointLoad(%s, %s, %s)' % (sstr(self.location), sstr(self.value), sstr('Load'))
        return str_sol

    __repr__ = __str__

    def location(self):
        return self._location

    def value(self):
        return self._value

    def moment(self):
        return self._moment


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
    >>> p = Point('4')
    >>> x = Symbol('x')
    >>> DistributedLoad(start = 3, end = 6, value = 2*x)
    DistributedLoad(3, 6, 2*x)

    """

    def __init__(self, start, end, value):
        self.start = start
        self.end = end
        self.value = value

    def __str__(self):
        str_sol = 'DistributedLoad(%s, %s, %s)' % (sstr(self.start), sstr(self.end), sstr(self.value))
        return str_sol

    __repr__ = __str__

    def start(self):
        return self._start

    def end(self):
        return self._end

    def value(self):
        return self._value
