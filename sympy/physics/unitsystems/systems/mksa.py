# -*- coding: utf-8 -*-

"""
MKS unit system.

MKS stands for "meter, kilogram, second, ampere".
"""

from __future__ import division

from sympy import pi
from sympy.physics.unitsystems.dimensions import (voltage, impedance,
                                                  conductance, capacitance, inductance, charge,
                                                  magnetic_density, magnetic_flux, current)
from sympy.physics.unitsystems.units import Unit, Constant
from sympy.physics.unitsystems.prefixes import PREFIXES, prefix_unit
from sympy.physics.unitsystems.systems.mks import mks_dim, mks


dims = (voltage, impedance, conductance, capacitance, inductance, charge,
        magnetic_density, magnetic_flux)

# dimension system
mksa_dim = mks_dim.extend(base=(current,), dims=dims, name='MKSA')

# base units
A = Unit(current, abbrev='A')

# derived units
V = Unit(voltage, factor=10**3, abbrev='V')
ohm = Unit(impedance, factor=10**3, abbrev='ohm')
# siemens
S = Unit(conductance, factor=10**-3, abbrev='S')
# farad
F = Unit(capacitance, factor=10**-3, abbrev='F')
# henry
H = Unit(inductance, factor=10**3, abbrev='H')
# coulomb
C = Unit(charge, abbrev='C')
# tesla
T = Unit(magnetic_density, abbrev='T')
# weber
Wb = Unit(magnetic_flux, abbrev='Wb')

# constants
# Wave impedance of free space
Z0 = Constant(impedance, factor=119.9169832*pi, abbrev='Z_0')

units = [A, V, ohm, S, F, H, C, T, Wb]
all_units = []
for u in units:
    all_units.extend(prefix_unit(u, PREFIXES))

all_units.extend([Z0])

mksa = mks.extend(base=(A,), units=all_units, name='MKSA')
