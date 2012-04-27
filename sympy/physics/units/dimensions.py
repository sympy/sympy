# -*- coding:utf-8 -*-

"""
Dimension for physical quantities.
"""

from __future__ import division
from copy import copy

from sympy import Rational, Number

# TODO: define the Dimension as a commutative Symbol, see paulialgebra module.
# TODO: define a dimensionless fixed dimension, which is returned instead of 1
#       in operations

class Dimension():
    """
    This class represent the dimension of a physical units.

    Examples
    ========

    >>> velocity = length / time
    >>> velocity
    {length: 1, time: -1}
    """

    def __init__(self, *args, **kwargs):

        # before setting the dict, check if a name or a symbol are defined
        # if so, remove them from the dict
        self.name = kwargs.pop('name', None)
        self.symbol = kwargs.pop('symbol', None)

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

        self.names, self.powers = map(tuple, zip(*sorted(pairs)))

    def __getitem__(self, key):
        return self.as_dict[key]

    def __len__(self):
        return len(self.names)

    def __hash__(self):
        try:
            return self._cached_hash
        except AttributeError:
            h = self._cached_hash = hash(frozenset((self.names, self.powers)))
            return h

    #TODO: use unit system to pretty print the dimension
    #def __str__(self):
        #if self.symbol is not None:
        #    return self.symbol
        #return dict.__str__(self)
        #s = ''.join()
    def __str__(self):
        return str(self.as_dict)

    def __repr__(self):
        return "<Dimension: %s>" % self.as_dict

    def __eq__(self, other):
        if isinstance(other, Dimension):
            return self.names == other.names and self.powers == other.powers
        else:
            raise TypeError("Dimension can be compared only to dimensions.")

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

        d = self.as_dict
        for key in other.as_dict:
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

        d = self.as_dict
        for key in other.as_dict:
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

        if isinstance(other, Number) or int(other) == other:
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
    def as_dict(self):
        return dict(self.items())

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def items(self):
        return zip(self.names, self.powers)

    @property
    def is_dimensionless(self):
        """
        Check if the dimension object really has a dimension.
        """

        for val in self.powers:
            if val != 0:
                return False
        else:
            return True
