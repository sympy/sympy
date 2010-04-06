"""
Physical units and dimensions.

The base class is Unit, where all here defined units (~200) inherit from.
"""

from sympy import Rational, pi
from sympy.core import Atom, Expr

class Unit(Atom, Expr):
    """
    Base class for all physical units.

    Create own units like:
    m = Unit("meter", "m")
    """
    is_positive = True    # make (m**2)**Rational(1,2) --> m
    is_commutative = True

    __slots__ = ["name", "abbrev"]

    def __new__(cls, name, abbrev, **assumptions):
        obj = Expr.__new__(cls, **assumptions)
        assert isinstance(name, str),`type(name)`
        assert isinstance(abbrev, str),`type(abbrev)`
        obj.name = name
        obj.abbrev = abbrev
        return obj

    def __getnewargs__(self):
        return (self.name, self.abbrev)

    def __eq__(self, other):
        return isinstance(other, Unit) and self.name == other.name

    def _hashable_content(self):
        return (self.name,self.abbrev)

# Delete this so it doesn't pollute the namespace
del Atom

def defunit(value, *names):
    u = value
    g = globals()
    for name in names:
        g[name] = u


# Dimensionless

percent = percents = Rational(1,100)
permille = permille = Rational(1,1000)

ten = Rational(10)

yotta = ten**24
zetta = ten**21
exa   = ten**18
peta  = ten**15
tera  = ten**12
giga  = ten**9
mega  = ten**6
kilo  = ten**3
deca  = ten**1
deci  = ten**-1
centi = ten**-2
milli = ten**-3
micro = ten**-6
nano  = ten**-9
pico  = ten**-12
femto = ten**-15
atto  = ten**-18
zepto = ten**-21
yocto = ten**-24

rad = radian = radians = 1
deg = degree = degrees = pi/180


# Base units

defunit(Unit('meter', 'm'), 'm', 'meter', 'meters')
defunit(Unit('kilogram', 'kg'), 'kg', 'kilogram', 'kilograms')
defunit(Unit('second', 's'), 's', 'second', 'seconds')
defunit(Unit('ampere', 'A'), 'A', 'ampere', 'amperes')
defunit(Unit('kelvin', 'K'), 'K', 'kelvin', 'kelvins')
defunit(Unit('mole', 'mol'), 'mol', 'mole', 'moles')
defunit(Unit('candela', 'cd'), 'cd', 'candela', 'candelas')


# Derived units

defunit(1/s, 'Hz', 'hz', 'hertz')
defunit(m*kg/s**2, 'N', 'newton', 'newtons')
defunit(N*m, 'J', 'joule', 'joules')
defunit(J/s, 'W', 'watt', 'watts')
defunit(N/m**2, 'Pa', 'pa', 'pascal', 'pascals')
defunit(s*A, 'C', 'coulomb', 'coulombs')
defunit(W/A, 'v', 'V', 'volt', 'volts')
defunit(V/A, 'ohm', 'ohms')
defunit(A/V, 'S', 'siemens', 'mho', 'mhos')
defunit(C/V, 'F', 'farad', 'farads')
defunit(J/A, 'Wb', 'wb', 'weber', 'webers')
defunit(V*s/m**2, 'T', 'tesla', 'teslas')
defunit(V*s/A, 'H', 'henry', 'henrys')


# Common length units

defunit(kilo*m, 'km', 'kilometer', 'kilometers')
defunit(deci*m, 'dm', 'decimeter', 'decimeters')
defunit(centi*m, 'cm', 'centimeter', 'centimeters')
defunit(milli*m, 'mm', 'millimeter', 'millimeters')
defunit(micro*m, 'um', 'micrometer', 'micrometers', 'micron', 'microns')
defunit(nano*m, 'nm', 'nanometer', 'nanometers')
defunit(pico*m, 'pm', 'picometer', 'picometers')

defunit(Rational('0.3048')*m, 'ft', 'foot', 'feet')
defunit(Rational('25.4')*mm, 'inch', 'inches')
defunit(3*ft, 'yd', 'yard', 'yards')
defunit(5280*ft, 'mi', 'mile', 'miles')


# Common volume and area units

defunit(m**3 / 1000, 'l', 'liter', 'liters')
defunit(deci*l, 'dl', 'deciliter', 'deciliters')
defunit(centi*l, 'cl', 'centiliter', 'centiliters')
defunit(milli*l, 'ml', 'milliliter', 'milliliters')


# Common time units

defunit(milli*s, 'ms', 'millisecond', 'milliseconds')
defunit(micro*s, 'us', 'microsecond', 'microseconds')
defunit(nano*s, 'ns', 'nanosecond', 'nanoseconds')
defunit(pico*s, 'ps', 'picosecond', 'picoseconds')

defunit(60*s, 'minute', 'minutes')
defunit(60*minute, 'h', 'hour', 'hours')
defunit(24*hour, 'day', 'days')

defunit(Rational('31558149.540')*s, 'sidereal_year', 'sidereal_years')
defunit(Rational('365.24219')*day, 'tropical_year', 'tropical_years')
defunit(Rational('365')*day, 'common_year', 'common_years')
defunit(Rational('365.25')*day, 'julian_year', 'julian_years')

year = years = tropical_year


# Common mass units

defunit(kilogram / kilo, 'g', 'gram', 'grams')
defunit(milli * g, 'mg', 'milligram', 'milligrams')
defunit(micro * g, 'ug', 'microgram', 'micrograms')



#----------------------------------------------------------------------------
# Physical constants
#

c = speed_of_light = 299792458 * m/s
G = gravitational_constant = Rational('6.67428') * ten**-11 * m**3 / kg / s**2
u0 = magnetic_constant = 4*pi * ten**-7 * N/A**2
e0 = electric_constant = 1/(u0 * c**2)
Z0 = vacuum_impedance = u0 * c

planck = Rational('6.62606896') * ten**-34 * J*s
hbar = planck / (2*pi)

avogadro = (Rational('6.02214179') * 10**23) / mol
boltzmann = Rational('1.3806505') * ten**-23 * J / K

gee = gees = Rational('9.80665') * m/s**2
atmosphere = atmospheres = atm = 101325 * pascal

kPa = kilo*Pa
bar = bars = 100*kPa
pound = pounds = 0.45359237 * kg * gee #exact
psi = pound / inch ** 2
dHg0 = 13.5951 # approx value at 0 C
mmHg = dHg0 * 9.80665 * Pa
amu = amus = gram / avogadro
quart = quarts = 231 * inch**3
eV = 1.602176487e-19 * J

# Other convenient units and magnitudes

defunit(c*julian_year, 'ly', 'lightyear', 'lightyears')
defunit(149597870691*m, 'au', 'astronomical_unit', 'astronomical_units')

# Delete this so it doesn't pollute the namespace
del Rational, pi
