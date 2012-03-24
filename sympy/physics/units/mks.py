"""
MKS unit system.

MKS stands for "meter, kilogram, second".
"""

from __future__ import division

from dimensions import Dimension

length = Dimension(symbol='L', length=1)
mass = Dimension(symbol='M', mass=1)
time = Dimension(symbol='T', time=1)
velocity = Dimension(length=1, time=-1)
energy = Dimension(mass=1, length=2, time=-2)


