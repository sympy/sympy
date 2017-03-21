# -*- coding: utf-8 -*-

"""
MKS unit system.

MKS stands for "meter, kilogram, second, ampere".
"""

from __future__ import division

from sympy.physics.unitsystems.definitions import A, V, C, S, ohm, F, H, Z0, Wb, T
from sympy.physics.unitsystems.dimensions import (voltage, impedance,
                                                  conductance, capacitance, inductance, charge,
                                                  magnetic_density, magnetic_flux, current)
from sympy.physics.unitsystems.prefixes import PREFIXES, prefix_unit
from sympy.physics.unitsystems.systems.mks import mks, _mks_dim

dims = (voltage, impedance, conductance, capacitance, inductance, charge,
        magnetic_density, magnetic_flux)

# dimension system
_mksa_dim = _mks_dim.extend(base=(current,), dims=dims, name='MKSA')


units = [A, V, ohm, S, F, H, C, T, Wb]
all_units = []
for u in units:
    all_units.extend(prefix_unit(u, PREFIXES))

all_units.extend([Z0])

mksa = mks.extend(base=(A,), units=all_units, name='MKSA')
