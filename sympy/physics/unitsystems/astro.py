"""
Astrophysical unit system.

Units usually used in astronomy and astrophysics.
"""

from __future__ import division

from math import pi
from sympy.physics.unitsystems.units import Constant, UnitSystem
from sympy.physics.unitsystems.mks import mks, length, time, mass

# sun-earth distance
ua = Constant(dimension=length, factor=1.495978707e11, abbrev='ua')
# mean radius
#earth_radius = Constant(dimension=length, factor=6.371e6, abbrev='earth_radius')
#earth_mass = Constant(dimension=mass, factor=5.9736e24, abbrev='earth_mass')

#day = Constant(dimension=time, factor=86400, abbrev='d')
#year = Constant(dimension=time, factor=365.25696*86400, abbrev='year')

#solar_mass = Constant(dimension=mass, factor=1.9891e30, abbrev='solar_mass')
#solar_radius = earth_radius = Constant(dimension=length, factor=6.96342e11, abbrev='solar_radius')

base_units = (ua,)
units = ()
astro = UnitSystem(base=base_units, units=units, name='Astrophysical')
