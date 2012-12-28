"""
MASK unit system.

MASK stands for "meter, ampere, second, kilogram."
"""

from __future__ import division

from math import pi
from sympy.physics.unitsystems.dimensions import Dimension
from sympy.physics.unitsystems.units import Unit, UnitSystem, Constant, PREFIXES

# base dimensions
length = Dimension(name='length', symbol='L', length=1)
mass = Dimension(name='mass', symbol='M', mass=1)
time = Dimension(name='time', symbol='T', time=1)
current = Dimension(name='current', symbol='I', current=1)

# derived dimensions
velocity = Dimension(name='velocity', length=1, time=-1)
acceleration = Dimension(name='acceleration', length=1, time=-2)
momentum = Dimension(name='momentum', mass=1, length=1, time=-1)
force = Dimension(name='force', mass=1, length=1, time=-2)
voltage = Dimension(name='voltage', symbol='U', mass=1, length=2, current=-1, time=-3)
impedance = Dimension(name='impedance', symbol='Z', mass=1, length=2, current=-2, time=-3)
frequency = Dimension(name='frequecy', symbol='f', time=-1)
energy = Dimension(name='energy', mass=1, length=2, time=-2)
power = Dimension(name='power', mass=1, length=2, time=-3)
pressure = Dimension(name='pressure', mass=1, length=-1, time=-2)

# base units
m = Unit(dimension=length, abbrev='m')
A = Unit(dimension=current, abbrev='A')
s = Unit(dimension=time, abbrev='s')
kg = Unit(dimension=mass, abbrev='g', prefix=PREFIXES['k'])

# derived units
v = Unit(dimension=velocity)
a = Unit(dimension=acceleration)
p = Unit(dimension=momentum)
V = Unit(dimension=voltage, factor=10**3, abbrev='V')
ohm = Unit(dimension=impedance, factor=10**3, abbrev='ohm')
Hz = Unit(dimension=frequency, abbrev='Hz')
J = Unit(dimension=energy, factor=10**3, abbrev='J')
N = Unit(dimension=force, factor=10**3, abbrev='N')
W = Unit(dimension=power, factor=10**3, abbrev='W')
Pa = Unit(dimension=pressure, factor=10**3, abbrev='Pa')

# constants
# Newton constant
G = Constant(m**3*kg**-1*s**-2, factor=6.67384e-11, abbrev='G')
# speed of light
c = Constant(dimension=velocity, factor=299792458, abbrev='c')
# Wave impedance of free space
Z0 = Constant(dimension=impedance, factor=119.9169832*pi, abbrev='Z_0')

units = (v, a, p, V, ohm, Hz, J, W, G, c, Z0)
mask = UnitSystem(base=(m, A, s, kg), units=units, name='MASK')
