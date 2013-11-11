# -*- coding: utf-8 -*-

"""
Physical quantities.
"""

from __future__ import division
import numbers

from sympy import sympify, AtomicExpr, Number
from .units import Unit


class Quantity(AtomicExpr):

    is_commutative = True

    def __new__(cls, factor=1, unit=None, **assumptions):

        if not isinstance(factor, str):
            factor = sympify(factor)

        # if the given unit is a number (because of some operations) and
        # the factor is represented as a number, then return a number
        if ((unit is None or isinstance(unit, (Number, numbers.Real)))
                    and isinstance(factor, (Number, numbers.Real))):
            return factor * (unit or 1)

        obj = AtomicExpr.__new__(cls, **assumptions)

        #TODO: if factor is of the form "1 m", parse the factor and the unit
        if isinstance(factor, (Number, numbers.Real)):
            if isinstance(unit, Unit):
                obj.factor, obj.unit = factor, unit
            else:
                raise TypeError("'unit' should be a Unit instance; %s found"
                                % type(unit))
        else:
            raise NotImplementedError

        return obj

    def __str__(self):
        return '%g %s' % (self.factor, self.unit)

    def __repr__(self):
        return '%g %s' % (self.factor, repr(self.unit))
