import time

from sympy import S

from sympy.core.numbers import pi

from sympy.physics.units import DimensionSystem, hertz, kilogram
from sympy.physics.units.definitions import (
    G, Hz, J, N, Pa, W, c, g, kg, m, s, meter, gram, second, newton,
    joule, watt, pascal)
from sympy.physics.units.definitions.dimension_definitions import (
    acceleration, action, energy, force, frequency, momentum,
    power, pressure, velocity, length, mass, time)
from sympy.physics.units.prefixes import PREFIXES, prefix_unit
from sympy.physics.units.prefixes import (
    kibi, mebi, gibi, tebi, pebi, exbi
)
from sympy.physics.units.definitions import (
    kilogram, newton, second, meter, gram, cd, K, joule, watt, pascal, hertz,
    coulomb, volt, ohm, siemens, farad, henry, tesla, weber, dioptre, lux,
    katal, gray, becquerel, inch, liter, julian_year, gravitational_constant,
    speed_of_light, elementary_charge, planck, hbar, electronvolt,
    avogadro_number, avogadro_constant, boltzmann_constant,
    stefan_boltzmann_constant, atomic_mass_constant, molar_gas_constant,
    faraday_constant, josephson_constant, von_klitzing_constant,
    acceleration_due_to_gravity, magnetic_constant, vacuum_permittivity,
    vacuum_impedance, coulomb_constant, atmosphere, bar, pound, psi, mmHg,
    milli_mass_unit, quart, lightyear, astronomical_unit, planck_mass,
    planck_time, planck_temperature, planck_length, planck_charge, planck_area,
    planck_volume, planck_momentum, planck_energy, planck_force, planck_power,
    planck_density, planck_energy_density, planck_intensity,
    planck_angular_frequency, planck_pressure, planck_current, planck_voltage,
    planck_impedance, planck_acceleration, bit, byte, kibibyte, mebibyte,
    gibibyte, tebibyte, pebibyte, exbibyte, curie, rutherford, radian, degree,
    steradian, angular_mil, atomic_mass_unit, gee, kPa, ampere, u0, c, kelvin,
    mol, mole, candela, m, kg, s, electric_constant, G, boltzmann
)


dimsys_length_weight_time = DimensionSystem([
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


One = S.One


# Base units:
dimsys_length_weight_time.set_quantity_dimension(meter, length)
dimsys_length_weight_time.set_quantity_scale_factor(meter, One)

# gram; used to define its prefixed units
dimsys_length_weight_time.set_quantity_dimension(gram, mass)
dimsys_length_weight_time.set_quantity_scale_factor(gram, One)

dimsys_length_weight_time.set_quantity_dimension(second, time)
dimsys_length_weight_time.set_quantity_scale_factor(second, One)

# derived units

dimsys_length_weight_time.set_quantity_dimension(newton, force)
dimsys_length_weight_time.set_quantity_scale_factor(newton, kilogram*meter/second**2)

dimsys_length_weight_time.set_quantity_dimension(joule, energy)
dimsys_length_weight_time.set_quantity_scale_factor(joule, newton*meter)

dimsys_length_weight_time.set_quantity_dimension(watt, power)
dimsys_length_weight_time.set_quantity_scale_factor(watt, joule/second)

dimsys_length_weight_time.set_quantity_dimension(pascal, pressure)
dimsys_length_weight_time.set_quantity_scale_factor(pascal, newton/meter**2)

dimsys_length_weight_time.set_quantity_dimension(hertz, frequency)
dimsys_length_weight_time.set_quantity_scale_factor(hertz, One)

# Other derived units:

dimsys_length_weight_time.set_quantity_dimension(dioptre, 1 / length)
dimsys_length_weight_time.set_quantity_scale_factor(dioptre, 1/meter)

# Common volume and area units

dimsys_length_weight_time.set_quantity_dimension(liter, length ** 3)
dimsys_length_weight_time.set_quantity_scale_factor(liter, meter**3 / 1000)


# Newton constant
# REF: NIST SP 959 (June 2019)

dimsys_length_weight_time.set_quantity_dimension(gravitational_constant, length ** 3 * mass ** -1 * time ** -2)
dimsys_length_weight_time.set_quantity_scale_factor(gravitational_constant, 6.67430e-11*m**3/(kg*s**2))

# speed of light

dimsys_length_weight_time.set_quantity_dimension(speed_of_light, velocity)
dimsys_length_weight_time.set_quantity_scale_factor(speed_of_light, 299792458*meter/second)


# Planck constant
# REF: NIST SP 959 (June 2019)

dimsys_length_weight_time.set_quantity_dimension(planck, action)
dimsys_length_weight_time.set_quantity_scale_factor(planck, 6.62607015e-34*joule*second)

# Reduced Planck constant
# REF: NIST SP 959 (June 2019)

dimsys_length_weight_time.set_quantity_dimension(hbar, action)
dimsys_length_weight_time.set_quantity_scale_factor(hbar, planck / (2 * pi))
