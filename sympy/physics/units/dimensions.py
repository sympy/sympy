# -*- coding:utf-8 -*-

"""
Dimension for physical quantities.
"""

from __future__ import division

from sympy import Rational, Number

# TODO: define the DImension as a commutative Symbol, see paulialgebra module.

class Dimension(dict):
    """
    This class represent the dimension of a physical units.

    Examples
    ========

    >>> velocity = length / time
    >>> velocity
    {length: 1, time: -1}
    """

    def __init__(self, *args, **kwargs):
        """
        Set the dimension dictionary.

        Before setting it, check if a name or a symbol are defined.
        """

        self.name = kwargs.get('name', None)
        self.symbol = kwargs.get('symbol', None)
        if self.name is not None:
            del kwargs['name']
        if self.symbol is not None:
            del kwargs['symbol']
        dict.__init__(self, *args, **kwargs)

    def __add__(self, other):
        """
        Addition of dimension has a sense only if they are the same (we don't
        add oranges and apples).
        """

        if not isinstance(other, Dimension):
            raise TypeError('Only dimension can be added.')
        if not self == other:
            raise TypeError('Only dimension which are equal can be added.')

        return self

    def __sub__(self, other):
        return self + other

    def __mul__(self, other):
        if not isinstance(other, Dimension):
            #TODO: improve to not raise error
            raise TypeError('Only dimension can be multiplied.')

        d = Dimension(self)
        for key in other:
            try:
                d[key] += other[key]
            except KeyError:
                d[key] = other[key]

        # if all dimensions are zero, then return 1 so that there is no more
        # dimensions
        if d.is_dimensionless:
            return 1
        else:
            return d

    def __div__(self, other):
        if not isinstance(other, Dimension):
            raise TypeError('Only dimension can be divided.')

        d = Dimension(self)
        for key in other:
            try:
                d[key] -= other[key]
            except KeyError:
                d[key] = -other[key]

        if d.is_dimensionless:
            return 1
        else:
            return d

        return d

    __truediv__ = __div__

    def __rdiv__(self, other):
        if other == 1:
            return Dimension([(x, -y) for x, y in self.items()])

    __rtruediv__ = __rdiv__

    def __neg__(self):
        return self

    def __pow__(self, other):

        if isinstance(other, Number) or int(other) == other:
            return Dimension([(x, y*other) for x, y in self.items()])
        elif isinstance(other, float):
            float_part = other % 1
            int_part = other / float_part
            float_power = 1 / float_part
            if int_part == 0:
                dim = self
            else:
                dim = self * int(int_part)
            return Dimension([(x, int(y / float_power)) for x, y in self.items()])

    #TODO: use unit system to pretty print the dimension
    #def __str__(self):
        #if self.symbol is not None:
        #    return self.symbol
        #return dict.__str__(self)
        #s = ''.join()

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
