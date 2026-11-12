"""
MKSA unit system.

MKSA stands for "meter, kilogram, second, ampere".
"""

from __future__ import annotations

from sympy.core.numbers import pi
from sympy.core.singleton import S
from sympy.physics.units.definitions import (
    Z0, ampere, coulomb, farad, henry, siemens, tesla, volt, weber, ohm,
    joule, meter, newton, second, speed_of_light, elementary_charge,
    magnetic_constant, vacuum_permittivity, vacuum_impedance,
    coulomb_constant)
from sympy.physics.units.definitions.dimension_definitions import (
    capacitance, charge, conductance, current, force, impedance, inductance,
    length, magnetic_density, magnetic_flux, voltage)
from sympy.physics.units.prefixes import PREFIXES, prefix_unit
from sympy.physics.units.systems.mks import MKS, dimsys_length_weight_time
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sympy.physics.units.quantities import Quantity

dims = (voltage, impedance, conductance, current, capacitance, inductance, charge,
        magnetic_density, magnetic_flux)

units = [ampere, volt, ohm, siemens, farad, henry, coulomb, tesla, weber]

all_units: list[Quantity] = []
for u in units:
    all_units.extend(prefix_unit(u, PREFIXES))
all_units.extend(units)

all_units.append(Z0)

dimsys_MKSA = dimsys_length_weight_time.extend([
    # Dimensional dependencies for base dimensions (MKSA not in MKS)
    current,
], new_dim_deps={
    # Dimensional dependencies for derived dimensions
    "voltage": {"mass": 1, "length": 2, "current": -1, "time": -3},
    "impedance": {"mass": 1, "length": 2, "current": -2, "time": -3},
    "conductance": {"mass": -1, "length": -2, "current": 2, "time": 3},
    "capacitance": {"mass": -1, "length": -2, "current": 2, "time": 4},
    "inductance": {"mass": 1, "length": 2, "current": -2, "time": -2},
    "charge": {"current": 1, "time": 1},
    "magnetic_density": {"mass": 1, "current": -1, "time": -2},
    "magnetic_flux": {"length": 2, "mass": 1, "current": -1, "time": -2},
})

One = S.One

dimsys_MKSA.set_quantity_scale_factor(ampere, One)

# derived units

dimsys_MKSA.set_quantity_scale_factor(coulomb, One)

dimsys_MKSA.set_quantity_scale_factor(volt, joule/coulomb)

dimsys_MKSA.set_quantity_scale_factor(ohm, volt/ampere)

dimsys_MKSA.set_quantity_scale_factor(siemens, ampere/volt)

dimsys_MKSA.set_quantity_scale_factor(farad, coulomb/volt)

dimsys_MKSA.set_quantity_scale_factor(henry, volt*second/ampere)

dimsys_MKSA.set_quantity_scale_factor(tesla, volt*second/meter**2)

dimsys_MKSA.set_quantity_scale_factor(weber, joule/ampere)

# elementary charge
# REF: NIST SP 959 (June 2019)

dimsys_MKSA.set_quantity_dimension(elementary_charge, charge)
dimsys_MKSA.set_quantity_scale_factor(elementary_charge, 1.602176634e-19*coulomb)

# magnetic constant:

dimsys_MKSA.set_quantity_dimension(magnetic_constant, force / current ** 2)
dimsys_MKSA.set_quantity_scale_factor(magnetic_constant, 4*pi/10**7 * newton/ampere**2)

# electric constant:

dimsys_MKSA.set_quantity_dimension(vacuum_permittivity, capacitance / length)
dimsys_MKSA.set_quantity_scale_factor(vacuum_permittivity, 1/(magnetic_constant * speed_of_light**2))

# vacuum impedance:

dimsys_MKSA.set_quantity_dimension(vacuum_impedance, impedance)
dimsys_MKSA.set_quantity_scale_factor(vacuum_impedance, magnetic_constant * speed_of_light)

# Coulomb's constant:

dimsys_MKSA.set_quantity_dimension(coulomb_constant, force * length ** 2 / charge ** 2)
dimsys_MKSA.set_quantity_scale_factor(coulomb_constant, 1/(4*pi*vacuum_permittivity))

MKSA = MKS.extend(base=(ampere,), units=all_units, name='MKSA', dimension_system=dimsys_MKSA, derived_units={
    magnetic_flux: weber,
    impedance: ohm,
    current: ampere,
    voltage: volt,
    inductance: henry,
    conductance: siemens,
    magnetic_density: tesla,
    charge: coulomb,
    capacitance: farad,
})
