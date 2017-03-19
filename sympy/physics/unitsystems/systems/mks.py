# -*- coding: utf-8 -*-

"""
MKS unit system.

MKS stands for "meter, kilogram, second".
"""

from __future__ import division

from sympy.physics.unitsystems.dimensions import (velocity, acceleration, momentum, force, energy, power, pressure,
                                                  frequency, action, length, mass, time)

from sympy.physics.unitsystems.simplifiers import qsimplify
from sympy.physics.unitsystems import (DimensionSystem, Unit,
                                       Constant, UnitSystem)
from sympy.physics.unitsystems.prefixes import PREFIXES, prefix_unit

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

# Prefixes of units like g, J, N etc get added using `prefix_unit`
# in the for loop, but the actual units have to be added manually.
all_units.extend([g, J, N, W, Pa, Hz])

for u in units:
    all_units.extend(prefix_unit(u, PREFIXES))
all_units.extend([v, a, p, G, c])

# unit system
mks = UnitSystem(base=(m, kg, s), units=all_units, name="MKS")
