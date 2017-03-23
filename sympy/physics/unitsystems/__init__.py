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

from .dimensions import Dimension, DimensionSystem
from .units import Unit, UnitSystem
from .simplifiers import convert_to
from .quantities import Quantity
from .definitions import A, C, F, G, H, Hz, J, N, Pa, \
    S, T, V, W, Wb, Z0, ampere, \
    anomalistic_year, anomalistic_years, c, \
    liter, liters, l, \
    deciliter, deciliters, dl, \
    centiliter, centiliters, cl, \
    milliliter, milliliters, ml, \
    common_year, common_years, coulomb, \
    day, days, \
    draconic_year, draconic_years, \
    eV, farad, \
    feet, foot, ft, \
    full_moon_cycle, \
    full_moon_cycles, g, gaussian_year, gaussian_years, gram, \
    gravitational_constant, hbar, henry, hertz, \
    h, hour, hours, \
    hz, inch, \
    inches, joule, julian_year, julian_years, kg, kilogram, l, \
    liter, liters, m, meter, mi, \
    ml, mile, miles, \
    minute, minutes, \
    nautical_mile, \
    nautical_miles, newton, nmi, ohm, pascal, pi, s, second, \
    sidereal_year, sidereal_years, siemens, speed_of_light, tesla,\
    tropical_year, tropical_years, volt, watt, weber, yard, \
    yards, yd, year, years,\
    km, kilometer, kilometers,\
    dm, decimeter, decimeters, \
    cm, centimeter, centimeters, \
    mm, millimeter, millimeters, \
    um, micrometer, micrometers, \
    nm, nanometer, nanometers, \
    pm, picometer, picometers, \
    ms, millisecond, milliseconds, \
    us, microsecond, microseconds, \
    ns, nanosecond, nanoseconds, \
    ps, picosecond, picoseconds

from .dimensions import acceleration, action, \
    capacitance, charge, conductance, current, energy, \
    force, frequency, impedance, inductance, length, magnetic_density, \
    magnetic_flux, mass, momentum, power, pressure, time, \
    velocity, voltage
