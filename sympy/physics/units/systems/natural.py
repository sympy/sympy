# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-

"""
Naturalunit system.

The natural system comes from "setting c = 1, hbar = 1". From the computer
point of view it means that we use velocity and action instead of length and
time. Moreover instead of mass we use energy.
"""

from __future__ import division

from sympy.physics.units.definitions import eV, hbar, c

from sympy.physics.units.dimensions import DimensionSystem
from sympy.physics.units.dimensions import length, mass, time, momentum,\
    force, energy, power, frequency, action, velocity
from sympy.physics.units.unitsystem import UnitSystem
from sympy.physics.units.prefixes import PREFIXES, prefix_unit

dims = (length, mass, time, momentum, force, energy, power, frequency)

# dimension system
_natural_dim = DimensionSystem(base=(action, energy, velocity), dims=dims,
                              name="Natural system")

units = prefix_unit(eV, PREFIXES)

# unit system
natural = UnitSystem(base=(hbar, eV, c), units=units, name="Natural system")
