# -*- coding: utf-8 -*-

"""
SI unit system.
Based on MKSA, which stands for "meter, kilogram, second, ampere".
Added kelvin, candela and mole.

"""

from __future__ import division

from sympy.physics.units.definitions import K, cd, lux, mol
from sympy.physics.units.dimensions import (
    amount_of_substance, luminous_intensity, temperature, dimsys_SI)
from sympy.physics.units.prefixes import PREFIXES, prefix_unit
from sympy.physics.units.systems.mksa import MKSA, _mksa_dim

derived_dims = ()
base_dims = (amount_of_substance, luminous_intensity, temperature)

# dimension system
_si_dim = dimsys_SI


units = [mol, cd, K, lux]
all_units = []
for u in units:
    all_units.extend(prefix_unit(u, PREFIXES))

all_units.extend([mol, cd, K, lux])

SI = MKSA.extend(base=(mol, cd, K), units=all_units, name='SI')
