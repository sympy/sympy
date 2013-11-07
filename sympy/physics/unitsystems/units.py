# -*- coding: utf-8 -*-

"""
Unit system for physical quantities.
"""

from __future__ import division

from sympy import sympify, AtomicExpr
from .dimensions import Dimension


class Unit(AtomicExpr):
    """
    Class for the units.
    """

    is_commutative = True
    is_number = False
    # make sqrt(m**2) --> m
    is_positive = True

    def __new__(cls, dim, abbrev="", factor=1, prefix=None, **assumptions):
        """
        Create a new unit instance.

        ``dim`` can be a Dimension or Unit object. The latter allows to
        construct derived units and constants. Note that the argument prefix
        is ignored if ``dim`` is a Unit instance.
        """

        factor = sympify(factor)

        obj = AtomicExpr.__new__(cls, **assumptions)

        if isinstance(dim, Dimension):
            obj._abbrev = abbrev
            obj._factor = factor
            obj.dim = dim
            obj.prefix = prefix
        elif isinstance(dim, Unit):
            obj._abbrev = abbrev or dim._abbrev
            obj._factor = factor * dim._factor
            obj.dim = dim.dim
            obj.prefix = None
        else:
            raise TypeError("'dim' object should be Unit or Dimension "
                            "instance; found %s" % type(dim))

        return obj
