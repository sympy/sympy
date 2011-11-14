"""
Physical units and dimensions.

The base class is Unit, where all here defined units (~200) inherit from.
"""

from sympy import Rational, pi
from sympy.core import AtomicExpr

class Unit(AtomicExpr):
    """
    Base class for all physical units.

    Create own units like:
    m = Unit("meter", "m")
    """
    is_positive = True    # make sqrt(m**2) --> m
    is_commutative = True
    is_number = False

    __slots__ = ["name", "abbrev"]

    def __new__(cls, name, abbrev, **assumptions):
        obj = AtomicExpr.__new__(cls, **assumptions)
        assert isinstance(name, str),repr(type(name))
        assert isinstance(abbrev, str),repr(type(abbrev))
        obj.name = name
        obj.abbrev = abbrev
        return obj

    def __getnewargs__(self):
        return (self.name, self.abbrev)

    def __eq__(self, other):
        return isinstance(other, Unit) and self.name == other.name

    def __hash__(self):
        return super(Unit, self).__hash__()

    def _hashable_content(self):
        return (self.name,self.abbrev)

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

m = meter = meters = Unit('meter', 'm')
kg = kilogram = kilograms = Unit('kilogram', 'kg')
s = second = seconds = Unit('second', 's')
A = ampere = amperes = Unit('ampere', 'A')
K = kelvin = kelvins = Unit('kelvin', 'K')
mol = mole = moles = Unit('mole', 'mol')
cd = candela = candelas = Unit('candela', 'cd')


# Derived units

Hz = hz = hertz = 1/s
N = newton = newtons = m*kg/s**2
J = joule = joules = N*m
W = watt = watts = J/s
Pa = pa = pascal = pascals = N/m**2
C = coulomb = coulombs = s*A
v = V = volt = volts = W/A
ohm = ohms = V/A
S = siemens = mho = mhos = A/V
F = farad = farads = C/V
Wb = wb = weber = webers = J/A
T = tesla = teslas = V*s/m**2
H = henry = henrys = V*s/A


# Common length units

km = kilometer = kilometers = kilo*m
dm = decimeter = decimeters = deci*m
cm = centimeter = centimeters = centi*m
mm = millimeter = millimeters = milli*m
um = micrometer = micrometers = micron = microns = micro*m
nm = nanometer = nanometers = nano*m
pm = picometer = picometers = pico*m

ft = foot = feet = Rational('0.3048')*m
inch = inches = Rational('25.4')*mm
yd = yard = yards = 3*ft
mi = mile = miles = 5280*ft


# Common volume and area units

l = liter = liters = m**3 / 1000
dl = deciliter = deciliters = deci*l
cl = centiliter = centiliters = centi*l
ml = milliliter = milliliters = milli*l


# Common time units

ms = millisecond = milliseconds = milli*s
us = microsecond = microseconds = micro*s
ns = nanosecond = nanoseconds = nano*s
ps = picosecond = picoseconds = pico*s

minute = minutes = 60*s
h = hour = hours = 60*minute
day = days = 24*hour

sidereal_year = sidereal_years = Rational('31558149.540')*s
tropical_year = tropical_years = Rational('365.24219')*day
common_year = common_years = Rational('365')*day
julian_year = julian_years = Rational('365.25')*day

year = years = tropical_year


# Common mass units

g = gram = grams = kilogram / kilo
mg = milligram = milligrams = milli * g
ug = microgram = micrograms = micro * g



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

ly = lightyear = lightyears = c*julian_year
au = astronomical_unit = astronomical_units = 149597870691*m

# Delete this so it doesn't pollute the namespace
del Rational, pi
