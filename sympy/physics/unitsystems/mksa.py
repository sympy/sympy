"""
MKSA unit system.

MKSA stands for "meter, ampere, second, kilogram."
"""

from __future__ import division

from math import pi
from sympy.physics.unitsystems.dimensions import Dimension
from sympy.physics.unitsystems.units import Unit, Constant
from sympy.physics.unitsystems.mks import mks

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

# base units
A = Unit(dimension=current, abbrev='A')

# derived units
V = Unit(dimension=voltage, factor=10**3, abbrev='V')
ohm = Unit(dimension=impedance, factor=10**3, abbrev='ohm')
# siemens
S = Unit(dimension=conductance, factor=10**-3, abbrev='S')
# farad
F = Unit(dimension=capacitance, factor=10**-3, abbrev='F')
# henry
H = Unit(dimension=inductance, factor=10**3, abbrev='H')
# coulomb
C = Unit(dimension=charge, abbrev='C')
# tesla
T = Unit(dimension=magnetic_density, abbrev='T')
# weber
Wb = Unit(dimension=magnetic_flux, abbrev='Wb')

# constants
# Wave impedance of free space
Z0 = Constant(dimension=impedance, factor=119.9169832*pi, abbrev='Z_0')

units = (V, ohm, S, F, H, C, T, Wb, Z0)
mksa = mks.extend(base=(A,), units=units, name='MKSA')
