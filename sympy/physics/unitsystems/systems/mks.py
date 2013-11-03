# -*- coding: utf-8 -*-

"""
MKS unit system.

MKS stands for "meter, kilogram, second".
"""

from __future__ import division

from sympy.physics.unitsystems.dimensions import Dimension

# base dimensions
length = Dimension(name="length", symbol="L", length=1)
mass = Dimension(name="mass", symbol="M", mass=1)
time = Dimension(name="time", symbol="T", time=1)

# derived dimensions
velocity = Dimension(name="velocity", length=1, time=-1)
acceleration = Dimension(name="acceleration", length=1, time=-2)
momentum = Dimension(name="momentum", mass=1, length=1, time=-1)
force = Dimension(name="force", symbol="F", mass=1, length=1, time=-2)
energy = Dimension(name="energy", symbol="E", mass=1, length=2, time=-2)
power = Dimension(name="power", mass=1, length=2, time=-3)
pressure = Dimension(name="pressure", mass=1, length=-1, time=-2)
frequency = Dimension(name="frequency", symbol="f", time=-1)
