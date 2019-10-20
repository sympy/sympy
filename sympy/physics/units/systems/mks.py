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

dimsys_MKS = DimensionSystem([
    # Dimensional dependencies for MKS base dimensions
    length,
    mass,
    time,
], dimensional_dependencies=dict(
    # Dimensional dependencies for derived dimensions
    velocity=dict(length=1, time=-1),
    acceleration=dict(length=1, time=-2),
    momentum=dict(mass=1, length=1, time=-1),
    force=dict(mass=1, length=1, time=-2),
    energy=dict(mass=1, length=2, time=-2),
    power=dict(length=2, mass=1, time=-3),
    pressure=dict(mass=1, length=-1, time=-2),
    frequency=dict(time=-1),
    action=dict(length=2, mass=1, time=-1),
    volume=dict(length=3),
))

# unit system
MKS = UnitSystem(base_units=(m, kg, s), units=all_units, name="MKS", dimension_system=dimsys_MKS)
