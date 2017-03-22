# -*- coding: utf-8 -*-

"""
Dimensional analysis and unit systems.

This module defines dimension/unit systems and physical quantities. It is
based on a group-theoretical construction where dimensions are represented as
vectors (coefficients being the exponents), and units are defined as a dimension
to which we added a scale.

Quantities are built from a factor and a unit, and are the basic objects that
one will use when doing computations.

All objects except systems and prefixes can be used in sympy expressions.
Note that as part of a CAS, various objects do not combine automatically
under operations.

Details about the implementation can be found in the documentation, and we
will not repeat all the explanations we gave there concerning our approach.
Ideas about future developments can be found on the `Github wiki
<https://github.com/sympy/sympy/wiki/Unit-systems>`_, and you should consult
this page if you are willing to help.
"""

from sympy.physics.unitsystems.dimensions import Dimension, DimensionSystem
from sympy.physics.unitsystems.units import Unit, UnitSystem
from .simplifiers import convert_to
from sympy.physics.unitsystems.quantities import Quantity
from .definitions import A, C, F, G, H, Hz, J, N, Pa, \
    Quantity, S, S_singleton, T, Unit, V, W, Wb, Z0, action, ampere, \
    anomalistic_year, anomalistic_years, c, capacitance, centiliter, \
    centiliters, charge, cl, common_year, common_years, conductance, coulomb, \
    current, day, deciliter, deciliters, dl, draconic_year, draconic_years, \
    eV, energy, farad, feet, foot, force, frequency, ft, full_moon_cycle, \
    full_moon_cycles, g, gaussian_year, gaussian_years, gram, \
    gravitational_constant, hbar, henry, hertz, hour, hz, impedance, inch, \
    inches, inductance, joule, julian_year, julian_years, kg, kilogram, l, \
    length, liter, liters, m, magnetic_density, magnetic_flux, mass, meter, mi, \
    mile, miles, milliliter, milliliters, minute, ml, nautical_mile, \
    nautical_miles, newton, nmi, ohm, pascal, pi, power, pressure, s, second, \
    sidereal_year, sidereal_years, siemens, speed_of_light, tesla, time, \
    tropical_year, tropical_years, velocity, volt, voltage, watt, weber, yard, \
    yards, yd, year, years
from .dimensions import acceleration, action, \
    capacitance, charge, conductance, current, energy, \
    force, frequency, impedance, inductance, length, magnetic_density, \
    magnetic_flux, mass, momentum, nsimplify, power, pressure, time, \
    velocity, voltage
