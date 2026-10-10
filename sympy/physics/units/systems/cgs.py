"""
CGS-Gaussian unit system.

The base units are centimeter, gram and second. The Coulomb constant is equal
to one, so that the electromagnetic quantities are expressed in terms of
length, mass and time.
"""
from __future__ import annotations
from sympy.core.numbers import pi
from sympy.core.singleton import S
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.physics.units import UnitSystem, centimeter, gram, second, coulomb, charge, speed_of_light, current, mass, \
    length, voltage, magnetic_density, magnetic_flux, capacitance, force, impedance
from sympy.physics.units.definitions import coulombs_constant
from sympy.physics.units.definitions.unit_definitions import statcoulomb, statampere, statvolt, volt, tesla, gauss, \
    weber, maxwell, debye, oersted, ohm, siemens, farad, henry, erg, ampere, coulomb_constant, elementary_charge, \
    vacuum_permittivity, magnetic_constant, vacuum_impedance
from sympy.physics.units.systems.mks import dimsys_length_weight_time

One = S.One

dimsys_cgs = dimsys_length_weight_time.extend(
    [],
    new_dim_deps={
        # Dimensional dependencies for derived dimensions
        "impedance": {"time": 1, "length": -1},
        "conductance": {"time": -1, "length": 1},
        "capacitance": {"length": 1},
        "inductance": {"time": 2, "length": -1},
        "charge": {"mass": S.Half, "length": S(3)/2, "time": -1},
        "current": {"mass": One/2, "length": 3*One/2, "time": -2},
        "voltage": {"length": One/2, "mass": One/2, "time": -1},
        "magnetic_density": {"length": -One/2, "mass": One/2, "time": -1},
        "magnetic_flux": {"length": 3*One/2, "mass": One/2, "time": -1},
    }
)

cgs_gauss = UnitSystem(
    base_units=[centimeter, gram, second],
    units=[],
    name="cgs_gauss",
    dimension_system=dimsys_cgs,
    defining_constants={coulomb_constant: One})


cgs_gauss.set_quantity_scale_factor(coulombs_constant, 1)

cgs_gauss.set_quantity_dimension(statcoulomb, charge)
cgs_gauss.set_quantity_scale_factor(statcoulomb, centimeter**(S(3)/2)*gram**(S.Half)/second)

cgs_gauss.set_quantity_dimension(coulomb, charge)

cgs_gauss.set_quantity_dimension(statampere, current)
cgs_gauss.set_quantity_scale_factor(statampere, statcoulomb/second)

cgs_gauss.set_quantity_dimension(statvolt, voltage)
cgs_gauss.set_quantity_scale_factor(statvolt, erg/statcoulomb)

cgs_gauss.set_quantity_dimension(volt, voltage)

cgs_gauss.set_quantity_dimension(gauss, magnetic_density)
cgs_gauss.set_quantity_scale_factor(gauss, sqrt(gram/centimeter)/second)

cgs_gauss.set_quantity_dimension(tesla, magnetic_density)

cgs_gauss.set_quantity_dimension(maxwell, magnetic_flux)
cgs_gauss.set_quantity_scale_factor(maxwell, sqrt(centimeter**3*gram)/second)

# SI units expressed in CGS-gaussian units:
cgs_gauss.set_quantity_scale_factor(coulomb, 10*speed_of_light*statcoulomb)
cgs_gauss.set_quantity_scale_factor(ampere, 10*speed_of_light*statcoulomb/second)
cgs_gauss.set_quantity_scale_factor(volt, 10**6/speed_of_light*statvolt)
cgs_gauss.set_quantity_scale_factor(weber, 10**8*maxwell)
cgs_gauss.set_quantity_scale_factor(tesla, 10**4*gauss)
cgs_gauss.set_quantity_dimension(debye, charge*length)
cgs_gauss.set_quantity_scale_factor(debye, One/10**18*statcoulomb*centimeter)
cgs_gauss.set_quantity_dimension(oersted, magnetic_density)
cgs_gauss.set_quantity_scale_factor(oersted, sqrt(gram/centimeter)/second)
cgs_gauss.set_quantity_scale_factor(ohm, 10**5/speed_of_light**2*second/centimeter)
cgs_gauss.set_quantity_scale_factor(siemens, One/10**5*speed_of_light**2*centimeter/second)
cgs_gauss.set_quantity_scale_factor(farad, One/10**5*speed_of_light**2*centimeter)
cgs_gauss.set_quantity_scale_factor(henry, 10**5/speed_of_light**2/centimeter*second**2)

# Coulomb's constant:
cgs_gauss.set_quantity_dimension(coulomb_constant, 1)
cgs_gauss.set_quantity_scale_factor(coulomb_constant, 1)

# Other constants of the SI expressed in CGS-gaussian units:
cgs_gauss.set_quantity_dimension(vacuum_permittivity, capacitance/length)
cgs_gauss.set_quantity_scale_factor(vacuum_permittivity, 1/(4*pi))

cgs_gauss.set_quantity_dimension(magnetic_constant, force/current**2)
cgs_gauss.set_quantity_scale_factor(magnetic_constant, 4*pi/speed_of_light**2)

cgs_gauss.set_quantity_dimension(vacuum_impedance, impedance)
cgs_gauss.set_quantity_scale_factor(vacuum_impedance, 4*pi/speed_of_light)

cgs_gauss.set_quantity_dimension(elementary_charge, charge)
cgs_gauss.set_quantity_scale_factor(elementary_charge, 1.602176634e-19*coulomb)

__all__ = [
    'ohm', 'tesla', 'maxwell', 'speed_of_light', 'volt', 'second', 'voltage',
    'debye', 'dimsys_length_weight_time', 'centimeter', 'coulomb_constant',
    'farad', 'sqrt', 'UnitSystem', 'current', 'charge', 'weber', 'gram',
    'statcoulomb', 'gauss', 'S', 'statvolt', 'oersted', 'statampere',
    'dimsys_cgs', 'coulomb', 'magnetic_density', 'magnetic_flux', 'One',
    'length', 'erg', 'mass', 'coulombs_constant', 'henry', 'ampere',
    'cgs_gauss', 'siemens', 'elementary_charge', 'vacuum_permittivity',
    'magnetic_constant', 'vacuum_impedance', 'capacitance', 'force',
    'impedance', 'pi',
]
