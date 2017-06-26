# -*- coding: utf-8 -*-

"""
MKS unit system.

MKS stands for "meter, kilogram, second, ampere".
"""

from __future__ import division

from sympy.physics.units.definitions import A, V, C, S, ohm, F, H, Z0, Wb, T
from sympy.physics.units.dimensions import (voltage, impedance,
                                            conductance, capacitance, inductance, charge,
                                            magnetic_density, magnetic_flux, current)
from sympy.physics.units.prefixes import PREFIXES, prefix_unit
from sympy.physics.units.systems.mks import MKS, _mks_dim

dims = (voltage, impedance, conductance, capacitance, inductance, charge,
        magnetic_density, magnetic_flux)

# dimension system
_mksa_dim = _mks_dim.extend(base=(current,), dims=dims, name='MKSA')


units = [A, V, ohm, S, F, H, C, T, Wb]
all_units = []
for u in units:
    all_units.extend(prefix_unit(u, PREFIXES))

all_units.extend([Z0])

MKSA = MKS.extend(base=(A,), units=all_units, name='MKSA')
