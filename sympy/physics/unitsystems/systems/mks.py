# -*- coding: utf-8 -*-

"""
MKS unit system.

MKS stands for "meter, kilogram, second".
"""

from __future__ import division

from sympy.physics.unitsystems.simplifiers import qsimplify
from sympy.physics.unitsystems import (Dimension, DimensionSystem, Unit,
                                       Constant, UnitSystem)
from sympy.physics.unitsystems.prefixes import PREFIXES, prefix_unit

# base dimensions
length = Dimension(name="length", symbol="L", length=1)
mass = Dimension(name="mass", symbol="M", mass=1)
time = Dimension(name="time", symbol="T", time=1)

# derived dimensions
velocity = Dimension(name="velocity", length=1, time=-1)
acceleration = Dimension(name="acceleration", length=1, time=-2)
momentum = Dimension(name="momentum", mass=1, length=1, time=-1)
force = Dimension(name="force", symbol="F", mass=1, length=1, time=-2)
energy = Dimension(name="energy", symbol="E", mass=1, length=2, time=-2)
power = Dimension(name="power", length=2, mass=1, time=-3)
pressure = Dimension(name="pressure", mass=1, length=-1, time=-2)
frequency = Dimension(name="frequency", symbol="f", time=-1)
action = Dimension(name="action", symbol="A", length=2, mass=1, time=-1)

dims = (velocity, acceleration, momentum, force, energy, power, pressure,
        frequency, action)

# dimension system
mks_dim = DimensionSystem(base=(length, mass, time), dims=dims, name="MKS")

# base units
m = Unit(length, abbrev="m")
kg = Unit(mass, abbrev="g", prefix=PREFIXES["k"])
s = Unit(time, abbrev="s")

# gram; used to define its prefixed units
g = Unit(mass, abbrev="g")

# derived units
v = Unit(velocity)
a = Unit(acceleration)
p = Unit(momentum)
J = Unit(energy, factor=10**3, abbrev="J")
N = Unit(force, factor=10**3, abbrev="N")
W = Unit(power, factor=10**3, abbrev="W")
Pa = Unit(pressure, factor=10**3, abbrev="Pa")
Hz = Unit(frequency, abbrev="Hz")

# constants
# Newton constant
G = Constant(qsimplify(m**3*kg**-1*s**-2).as_unit, factor=6.67384e-11, abbrev="G")
# speed of light
c = Constant(velocity, factor=299792458, abbrev="c")

units = [m, g, s, J, N, W, Pa, Hz]
all_units = []
for u in units:
    all_units.extend(prefix_unit(u, PREFIXES))
all_units.extend([v, a, p, G, c])

# unit system
mks = UnitSystem(base=(m, kg, s), units=all_units, name="MKS")
