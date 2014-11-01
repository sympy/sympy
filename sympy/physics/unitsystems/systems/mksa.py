# -*- coding: utf-8 -*-

"""
MKS unit system.

MKS stands for "meter, kilogram, second, ampere".
"""

from __future__ import division

from sympy import pi
from sympy.physics.unitsystems.dimensions import Dimension, DimensionSystem
from sympy.physics.unitsystems.units import Unit, Constant, UnitSystem
from sympy.physics.unitsystems.prefixes import PREFIXES, prefix_unit
from sympy.physics.unitsystems.systems.mks import mks_dim, mks

# base dimensions
current = Dimension(name='current', symbol='I', current=1)

# derived dimensions
voltage = Dimension(name='voltage', symbol='U', mass=1, length=2, current=-1,
                    time=-3)
impedance = Dimension(name='impedance', symbol='Z', mass=1, length=2,
                      current=-2, time=-3)
conductance = Dimension(name='conductance', symbol='G', mass=-1, length=-2,
                      current=2, time=3)
capacitance = Dimension(name='capacitance', mass=-1, length=-2, current=2,
                        time=4)
inductance = Dimension(name='inductance', mass=1, length=2, current=-2, time=-2)
charge = Dimension(name='charge', symbol='Q', current=1, time=1)
magnetic_density = Dimension(name='charge', symbol='B', mass=1, current=-1,
                             time=-2)
magnetic_flux = Dimension(name='charge', length=2, mass=1, current=-1, time=-2)

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
