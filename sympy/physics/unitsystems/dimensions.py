# -*- coding:utf-8 -*-

"""
Dimension for physical quantities.
"""

from __future__ import division
from copy import copy

from sympy.core.containers import Dict
from sympy import Rational, Number, sympify

# TODO: define the Dimension as a commutative Symbol, see paulialgebra module.
# TODO: define a dimensionless fixed dimension, which is returned instead of 1
#       in operations

class Dimension(Dict):
    """
    This class represent the dimension of a physical units.

    The dimensions should have a name and possibly a symbol. All other
    arguments are dimensional powers.

    Examples
    ========

    >>> from sympy.physics.unitsystems.mks import length, time
    >>> velocity = length / time
    >>> velocity
    {length: 1, time: -1}
    """

    def __new__(cls, *args, **kwargs):

        # before setting the dict, check if a name or a symbol are defined
        # if so, remove them from the dict
        name = kwargs.pop('name', None)
        symbol = kwargs.pop('symbol', None)

        pairs = []

        # do more advanced check to verify that the names tuple contains
        # only str and the powers tuple only integers
        for arg in args:
            if isinstance(arg, dict):
                arg = copy(arg)
                pairs.extend(arg.items())
            if isinstance(arg, (tuple, list)):
                for p in arg:
                    if len(p) != 2:
                        raise TypeError("Tuple and list argument must have "
                                        "length 2.")
                pairs.extend(arg)

        pairs.extend(kwargs.items())

        # filter dimension set to zero; this avoid the following odd result:
        # Dimension(length=1) == Dimension(length=1, mass=0) => False
        pairs = [pair for pair in pairs if pair[1] != 0]

        new = Dict.__new__(cls, *pairs)
        new.name = name
        new.symbol = symbol

        return new

    #TODO: use unit system to pretty print the dimension
    #def __str__(self):
        #if self.symbol is not None:
        #    return self.symbol
        #return dict.__str__(self)
        #s = ''.join()

    def __str__(self):
        if self.symbol is not None:
            return self.symbol
        elif self.name is not None:
            return self.name
        else:
            return repr(self)

    def __add__(self, other):
        """
        Addition of dimension has a sense only if they are the same (we don't
        add oranges and apples).
        """

        if not isinstance(other, Dimension):
            raise TypeError('Only dimension can be added.')
        if self != other:
            raise TypeError('Only dimension which are equal can be added.')

        return self

    def __sub__(self, other):
        return self + other

    def __mul__(self, other):
        if not isinstance(other, Dimension):
            #TODO: improve to not raise error
            raise TypeError('Only dimension can be multiplied.')

        d = dict(self)
        for key in other:
            try:
                d[key] += other[key]
            except KeyError:
                d[key] = other[key]
        d = Dimension(d)

        # if all dimensions are zero, then return 1 so that there is no more
        # dimensions
        if d.is_dimensionless:
            return 1
        else:
            return d

    def __div__(self, other):
        if not isinstance(other, Dimension):
            raise TypeError('Only dimension can be divided.')

        d = dict(self)
        for key in other:
            try:
                d[key] -= other[key]
            except KeyError:
                d[key] = -other[key]
        d = Dimension(d)

        if d.is_dimensionless:
            return 1
        else:
            return d

    __truediv__ = __div__

    def __rdiv__(self, other):
        if other == 1:
            return Dimension([(x, -y) for x, y in self.items()])

    __rtruediv__ = __rdiv__

    def __neg__(self):
        return self

    def __pow__(self, other):

        other = sympify(other)
        #TODO: check consistency when having rational, float...
        if isinstance(other, Number):
            return Dimension([(x, y*other) for x, y in self.items()])
        """
        elif isinstance(other, float):
            float_part = other % 1
            int_part = other / float_part
            float_power = 1 / float_part
            if int_part == 0:
                dim = self
            else:
                dim = self * int(int_part)
            return Dimension([(x, int(y / float_power)) for x, y in self.items()])
        """

    @property
    def is_dimensionless(self):
        """
        Check if the dimension object really has a dimension.
        """

        for key in self:
            if self[key] != 0:
                return False
        else:
            return True

    def dimension_is_integer(self):
        """
        Check if all the powers are integers.
        """

        for key in self:
            if int(self[key]) != self[key]:
                return False
        return True
