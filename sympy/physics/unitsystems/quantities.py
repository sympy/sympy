# -*- coding: utf-8 -*-

"""
Physical quantities.
"""

from __future__ import division
import numbers

from sympy import sympify, AtomicExpr, Number, Mul, Pow
from .units import Unit

#TODO: in operations, interpret a Unit as a quantity with factor 1?


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
        return "%g %s" % (self.factor, self.unit)

    def __repr__(self):
        return "%g %s" % (self.factor, repr(self.unit))

    def __eq__(self, other):

        return (isinstance(other, Quantity) and self.factor == other.factor
                and self.unit == other.unit)

    def __neg__(self):
        return Quantity(-self.factor, self.unit)

    def __add__(self, other):

        if isinstance(other, Quantity):
            return Quantity(self.factor + other.convert_to(self.unit).factor,
                            self.unit)
        else:
            raise TypeError("Only quantities can be added")

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):

        if isinstance(other, Quantity):
            return Quantity(self.factor - other.convert_to(self.unit).factor,
                            self.unit)
        else:
            raise TypeError("Only quantities can be subtracted")

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):

        other = sympify(other)

        if isinstance(other, Quantity):
            return Quantity(self.factor * other.factor, self.unit * other.unit)
        elif isinstance(other, (Number, numbers.Real)):
            return Quantity(self.factor * other, self.unit)
        else:
            return Mul(self, other)

    def __rmul__(self, other):

        return self * other

    def __div__(self, other):

        other = sympify(other)
        if isinstance(other, Quantity):
            return Quantity(self.factor / other.factor, self.unit / other.unit)
        elif isinstance(other, (Number, numbers.Real)):
            return Quantity(self.factor / other, self.unit)
        else:
            return Mul(self, Pow(other, -1))

    __truediv__ = __div__

    def __rdiv__(self, other):

        other = sympify(other)
        if isinstance(other, Quantity):
            return Quantity(other.factor / self.factor, other.unit / self.unit)
        elif isinstance(other, (Number, numbers.Real)):
            return Quantity(other / self.factor, self.unit**-1)
        else:
            return Mul(self**-1, other)

    __rtruediv__ = __rdiv__

    def __pow__(self, other):

        other = sympify(other)
        if isinstance(other, (Number, numbers.Real)):
            f = self.factor**other
            # without evalf a Pow instance is returned, and it can not be
            # handled by Quantity.__new__
            return Quantity(f.evalf(), self.unit**other)
        else:
            return Pow(self, other)

    def convert_to(self, unit):
        """
        Convert the quantity to another (compatible) unit.
        """

        if self.unit.is_compatible(unit) is False:
            raise ValueError("Only compatible units can be converted; "
                                 "'%s' found" % unit.dim)

        return Quantity(self.factor * self.unit.factor / unit.factor, unit)
