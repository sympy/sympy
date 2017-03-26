# -*- coding: utf-8 -*-

"""
MKS unit system.

MKS stands for "meter, kilogram, second".
"""

from __future__ import division

from sympy.physics.units.definitions import (m, kg, s, J, N, W, Pa, Hz, g, G, c)
from sympy.physics.units.dimensions import (velocity, acceleration, momentum, force, energy, power, pressure,
                                            frequency, action, length, mass, time)

from sympy.physics.units import DimensionSystem, UnitSystem
from sympy.physics.units.prefixes import PREFIXES, prefix_unit

dims = (velocity, acceleration, momentum, force, energy, power, pressure,
        frequency, action)

# dimension system
_mks_dim = DimensionSystem(base=(length, mass, time), dims=dims, name="MKS")

units = [m, g, s, J, N, W, Pa, Hz]
all_units = []

# Prefixes of units like g, J, N etc get added using `prefix_unit`
# in the for loop, but the actual units have to be added manually.
all_units.extend([g, J, N, W, Pa, Hz])

for u in units:
    all_units.extend(prefix_unit(u, PREFIXES))
all_units.extend([G, c])

# unit system
MKS = UnitSystem(base=(m, kg, s), units=all_units, name="MKS")
