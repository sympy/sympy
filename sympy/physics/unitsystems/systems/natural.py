# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-

"""
Naturalunit system.

The natural system comes from "setting c = 1, hbar = 1". From the computer
point of view it means that we use velocity and action instead of length and
time. Moreover instead of mass we use energy.
"""

from __future__ import division

from sympy import Symbol

from sympy.physics.unitsystems.dimensions import DimensionSystem
from sympy.physics.unitsystems.dimensions import length, mass, time, momentum,\
    force, energy, power, frequency, action, velocity
from sympy.physics.unitsystems.units import Unit, Constant, UnitSystem
from sympy.physics.unitsystems.prefixes import PREFIXES, prefix_unit

dims = (length, mass, time, momentum, force, energy, power, frequency)

# dimension system
natural_dim = DimensionSystem(base=(action, energy, velocity), dims=dims,
                              name="Natural system")

# base units
hbar = Constant(action, factor=1.05457266e-34, abbrev="hbar")
eV = Unit(energy, factor=1.60219e-19, abbrev="eV")
c = Constant(velocity, factor=299792458, abbrev="c")

units = prefix_unit(eV, PREFIXES)

# unit system
natural = UnitSystem(base=(hbar, eV, c), units=units, name="Natural system")
