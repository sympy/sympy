"""
MKS unit system.

MKS stands for "meter, kilogram, second".
"""

from __future__ import division

from sympy.physics.unitsystems.dimensions import Dimension
from sympy.physics.unitsystems.units import Unit, UnitSystem, Constant, PREFIXES

# base dimensions
length = Dimension(name='length', symbol='L', length=1)
mass = Dimension(name='mass', symbol='M', mass=1)
time = Dimension(name='time', symbol='T', time=1)

# derived dimensions
velocity = Dimension(name='velocity', length=1, time=-1)
acceleration = Dimension(name='acceleration', length=1, time=-2)
momentum = Dimension(name='momentum', mass=1, length=1, time=-1)
force = Dimension(name='force', symbol='F', mass=1, length=1, time=-2)
energy = Dimension(name='energy', symbol='E', mass=1, length=2, time=-2)
power = Dimension(name='power', mass=1, length=2, time=-3)
pressure = Dimension(name='pressure', mass=1, length=-1, time=-2)
frequency = Dimension(name='frequency', symbol='f', time=-1)

# base units
m = Unit(dimension=length, abbrev='m')
kg = Unit(dimension=mass, abbrev='g', prefix=PREFIXES['k'])
s = Unit(dimension=time, abbrev='s')

# derived units
v = Unit(dimension=velocity)
a = Unit(dimension=acceleration)
p = Unit(dimension=momentum)
J = Unit(dimension=energy, factor=10**3, abbrev='J')
N = Unit(dimension=force, factor=10**3, abbrev='N')
W = Unit(dimension=power, factor=10**3, abbrev='W')
Pa = Unit(dimension=pressure, factor=10**3, abbrev='Pa')
Hz = Unit(dimension=frequency, abbrev='Hz')

# constants
# Newton constant
G = Constant(m**3*kg**-1*s**-2, factor=6.67384e-11, abbrev='G')
# speed of light
c = Constant(dimension=velocity, factor=299792458, abbrev='c')

units = (v, a, p, J, N, W, Pa, Hz, G, c)
mks = UnitSystem(base=(m, kg, s), units=units, name='MKS')
