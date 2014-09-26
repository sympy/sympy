# -*- coding: utf-8 -*-

"""
Physical quantities.
"""

from __future__ import division
import numbers

from sympy import sympify, Expr, Number, Mul, Pow
from .units import Unit

#TODO: in operations, interpret a Unit as a quantity with factor 1


class Quantity(Expr):
    """
    Physical quantity.

    A quantity is defined from a factor and a unit.
    """

    is_commutative = True

    def __new__(cls, factor=1, unit=None, **assumptions):

        if not isinstance(factor, str):
            factor = sympify(factor)

        # if the given unit is a number (because of some operations) and
        # the factor is represented as a number, then return a number
        if ((unit is None or isinstance(unit, (Number, numbers.Real)))
                    and isinstance(factor, (Number, numbers.Real))):
            return factor * (unit or 1)

        #TODO: if factor is of the form "1 m", parse the factor and the unit
        if isinstance(factor, (Number, numbers.Real)):
            if not isinstance(unit, Unit):
                raise TypeError("'unit' should be a Unit instance; %s found"
                                % type(unit))
        else:
            raise NotImplementedError

        obj = Expr.__new__(cls, factor, unit, **assumptions)
        obj.factor, obj.unit = factor, unit

        return obj

    def __str__(self):
        return "%g %s" % (self.factor, self.unit)

    def __repr__(self):
        return "%g %s" % (self.factor, repr(self.unit))

    def __neg__(self):
        return Quantity(-self.factor, self.unit)

    def add(self, other):
        """
        Add two quantities.

        If the other object is not a quantity, raise an error.
        Two quantities can be added only if they have the same unit: so we
        convert first the other quantity to the same unit and, if it succedded,
        then we add the factors.
        """

        if isinstance(other, Quantity):
            return Quantity(self.factor + other.convert_to(self.unit).factor,
                            self.unit)
        else:
            raise TypeError("Only quantities can be added")

    def sub(self, other):

        if isinstance(other, Quantity):
            return Quantity(self.factor - other.convert_to(self.unit).factor,
                            self.unit)
        else:
            raise TypeError("Only quantities can be subtracted")

    def mul(self, other):

        other = sympify(other)

        if isinstance(other, Quantity):
            return Quantity(self.factor * other.factor,
                            self.unit.mul(other.unit))
        elif isinstance(other, (Number, numbers.Real)):
            return Quantity(self.factor * other, self.unit)
        else:
            return Mul(self, other)

    def div(self, other):

        other = sympify(other)
        if isinstance(other, Quantity):
            return Quantity(self.factor / other.factor,
                            self.unit.div(other.unit))
        elif isinstance(other, (Number, numbers.Real)):
            return Quantity(self.factor / other, self.unit)
        else:
            return Mul(self, Pow(other, -1))

    def rdiv(self, other):

        other = sympify(other)
        if isinstance(other, Quantity):
            return Quantity(other.factor / self.factor,
                            other.unit.div(self.unit))
        elif isinstance(other, (Number, numbers.Real)):
            return Quantity(other / self.factor, self.unit.pow(-1))
        else:
            return Mul(self**-1, other)

    def pow(self, other):

        other = sympify(other)
        if isinstance(other, (Number, numbers.Real)):
            f = self.factor**other
            # without evalf a Pow instance is returned, and it can not be
            # handled by Quantity.__new__
            return Quantity(f.evalf(), self.unit.pow(other))
        else:
            return Pow(self, other)

    @property
    def as_unit(self):
        """
        Convert the quantity to a unit.
        """

        from .units import Unit
        return Unit(self.unit, factor=self.factor)

    def convert_to(self, unit):
        """
        Convert the quantity to another (compatible) unit.
        """

        if self.unit.is_compatible(unit) is False:
            raise ValueError("Only compatible units can be converted; "
                                 "'%s' found" % unit.dim)

        return Quantity(self.factor * self.unit.factor / unit.factor, unit)
