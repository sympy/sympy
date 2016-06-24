"""
This module can be used to solve beam bending problems in mechanics.

"""

from __future__ import print_function, division

from sympy.printing import sstr
from sympy.physics.mechanics import Point
from sympy import sympify


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
