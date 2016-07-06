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
        The starting point of the applied load.
    order : Sympifyable
        The order of the applied load.
    value : Sympifyable
        A SymPy expression representing the value of the applied load.

    Examples
    ========
    >>> from sympy.physics.mechanics.beam import DistributedLoad
    >>> from sympy.physics.mechanics import Point
    >>> from sympy import Symbol
    >>> a = Point('4')
    >>> b = 2
    >>> DistributedLoad(start = a, order = b, value = 2)
    DistributedLoad(4, 2, 2)

    """

    def __init__(self, start, order, value):
        self._start = start
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
