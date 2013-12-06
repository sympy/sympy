# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-

"""
Naturalunit system.

The natural system comes from "setting c = 1, hbar = 1". From the computer
point of view it means that we use velocity and action instead of length and
time. Moreover instead of mass we use energy.
"""

from __future__ import division

from sympy.physics.unitsystems.dimensions import Dimension, DimensionSystem
from sympy.physics.unitsystems.units import Unit, Constant, UnitSystem
from sympy.physics.unitsystems.prefixes import PREFIXES, prefix_unit

# base dimensions
action = Dimension(name="action", symbol="A", length=2, mass=1, time=-1)
energy = Dimension(name="energy", symbol="E", length=2, mass=1, time=-2)
velocity = Dimension(name="velocity", symbol="V", length=1, time=-1)

# derived dimensions
length = Dimension(name="length", symbol="L", length=1)
mass = Dimension(name="mass", symbol="M", mass=1)
time = Dimension(name="time", symbol="T", time=1)
acceleration = Dimension(name="acceleration", length=1, time=-2)
momentum = Dimension(name="momentum", mass=1, length=1, time=-1)
force = Dimension(name="force", symbol="F", mass=1, length=1, time=-2)
power = Dimension(name="power", length=2, mass=1, time=-3)
frequency = Dimension(name="frequency", symbol="f", time=-1)

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
