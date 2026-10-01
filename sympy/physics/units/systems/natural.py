"""
Naturalunit system.

The natural system comes from "setting c = 1, hbar = 1". From the computer
point of view it means that we use velocity and action instead of length and
time. Moreover instead of mass we use energy.
"""
from __future__ import annotations

from sympy.physics.units import DimensionSystem
from sympy.physics.units.definitions import c, eV, hbar, joule
from sympy.physics.units.definitions.dimension_definitions import (
    action, energy, velocity)
from sympy.physics.units.prefixes import PREFIXES, prefix_unit
from sympy.physics.units.systems.length_weight_time import dimsys_length_weight_time
from sympy.physics.units.unitsystem import UnitSystem


# dimension system
_natural_dim = DimensionSystem(
    base_dims=(action, energy, velocity),
    dimensional_dependencies={
        "length": {"action": 1, "velocity": 1, "energy": -1},
        "mass": {"energy": 1, "velocity": -2},
        "time": {"action": 1, "energy": -1},
        "acceleration": {"energy": 1, "velocity": 1, "action": -1},
        "momentum": {"energy": 1, "velocity": -1},
        "force": {"energy": 2, "action": -1, "velocity": -1},
        "power": {"energy": 2, "action": -1},
        "pressure": {"energy": 4, "action": -3, "velocity": -3},
        "frequency": {"energy": 1, "action": -1},
        "area": {"action": 2, "velocity": 2, "energy": -2},
        "volume": {"action": 3, "velocity": 3, "energy": -3},
    })

_natural_dim._quantity_dimension_map.update(dimsys_length_weight_time._quantity_dimension_map)
_natural_dim._quantity_scale_factors.update(dimsys_length_weight_time._quantity_scale_factors)

_natural_dim.set_quantity_dimension(eV, energy)
_natural_dim.set_quantity_scale_factor(eV, 1.602176634e-19*joule)

units = prefix_unit(eV, PREFIXES)

# unit system
natural = UnitSystem(base_units=(hbar, eV, c), units=units, name="Natural system", dimension_system=_natural_dim)
