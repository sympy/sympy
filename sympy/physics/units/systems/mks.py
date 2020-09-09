"""
MKS unit system.

MKS stands for "meter, kilogram, second".
"""

from __future__ import division

from sympy.physics.units import UnitSystem, DimensionSystem
from sympy.physics.units.definitions import G, Hz, J, N, Pa, W, c, g, kg, m, s
from sympy.physics.units.definitions.dimension_definitions import (
    acceleration, action, energy, force, frequency, momentum,
    power, pressure, velocity, length, mass, time)
from sympy.physics.units.prefixes import PREFIXES, prefix_unit
from sympy.physics.units.systems.length_weight_time import dimsys_length_weight_time

dims = (velocity, acceleration, momentum, force, energy, power, pressure,
        frequency, action)

units = [m, g, s, J, N, W, Pa, Hz]
all_units = []

# Prefixes of units like g, J, N etc get added using `prefix_unit`
# in the for loop, but the actual units have to be added manually.
all_units.extend([g, J, N, W, Pa, Hz])

for u in units:
    all_units.extend(prefix_unit(u, PREFIXES))
all_units.extend([G, c])

# unit system
MKS = UnitSystem(base_units=(m, kg, s), units=all_units, name="MKS", dimension_system=dimsys_length_weight_time)


__all__ = [
    'force', 'division', 'DimensionSystem', 'energy', 'Pa', 'MKS',
    'dimsys_length_weight_time', 'Hz', 'power', 's', 'UnitSystem', 'units',
    'mass', 'momentum', 'acceleration', 'G', 'J', 'N', 'pressure', 'W',
    'all_units', 'c', 'kg', 'g', 'dims', 'prefix_unit', 'm', 'PREFIXES',
    'length', 'frequency', 'u', 'time', 'action', 'velocity',
]
