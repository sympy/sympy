"""
MKS unit system.

MKS stands for "meter, kilogram, second".
"""

from __future__ import division

from sympy.physics.units.dimensions import Dimension
from sympy.physics.units.units import Unit, UnitSystem

length = Dimension(name='length', symbol='L', length=1)
mass = Dimension(name='mass', symbol='M', mass=1)
time = Dimension(name='time', symbol='T', time=1)
velocity = Dimension(name='velocity', length=1, time=-1)
energy = Dimension(name='energy', mass=1, length=2, time=-2)

m = Unit(abbrev='m', dimension=length)
kg = Unit(abbrev='kg', dimension=mass)
s = Unit(abbrev='s', dimension=time)

v = Unit(abbrev='v', dimension=velocity)
J = Unit(abbrev='J', dimension=energy)
#kg * m**2 * s-2

mks = UnitSystem(base=(m, kg, s), units=(J,), name='MKS')
