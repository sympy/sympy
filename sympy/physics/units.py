"""
Physical units and dimensions.

The base class is Unit, where all here defined units (~200) inherit from.

The find_unit function can help you find units for a given quantity:

    >>> import sympy.physics.units as u
    >>> u.find_unit('coul')
    ['coulomb', 'coulombs']
    >>> u.find_unit(u.charge)
    ['C', 'charge', 'coulomb', 'coulombs']
    >>> u.coulomb
    A*s

Units are always given in terms of base units that have a name and
an abbreviation:

    >>> u.A.name
    'ampere'
    >>> u.ampere.abbrev
    'A'

The generic name for a unit (like 'length', 'mass', etc...)
can help you find units:

    >>> u.find_unit('magnet')
    ['magnetic_flux', 'magnetic_constant', 'magnetic_flux_density']
    >>> u.find_unit(u.magnetic_flux)
    ['Wb', 'wb', 'weber', 'webers', 'magnetic_flux']

If, for a given session, you wish to add a unit you may do so:

    >>> u.find_unit('gal')
    []
    >>> u.gal = 4*u.quart
    >>> u.gal/u.inch**3
    231

To see a given quantity in terms of some other unit, divide by the desired
unit:

    >>> mph = u.miles/u.hours
    >>> (u.m/u.s/mph).n(2)
    2.2

The units are defined in terms of base units, so when you divide similar
units you will obtain a pure number. This means, for example, that if you
divide a real-world mass (like grams) by the atomic mass unit (amu) you
will obtain Avogadro's number. To obtain the answer in moles you
should divide by the unit ``avogadro``:

    >>> u.grams/u.amu
    602214179000000000000000
    >>> _/u.avogadro
    mol

For chemical calculations the unit ``mmu`` (molar mass unit) has been
defined so this conversion is handled automatically. For example, the
number of moles in 1 kg of water might be calculated as:

    >>> u.kg/(18*u.mmu).n(3)
    55.5*mol

If you need the number of atoms in a mol as a pure number you can use
``avogadro_number`` but if you need it as a dimensional quantity you should use
``avogadro_constant``. (``avogadro`` is a shorthand for the dimensional
quantity.)

    >>> u.avogadro_number
    602214179000000000000000
    >>> u.avogadro_constant
    602214179000000000000000/mol
"""

from __future__ import print_function, division

from sympy import Rational, pi
from sympy.core import AtomicExpr


class Unit(AtomicExpr):
    """
    Base class for base unit of physical units.

    >>> from sympy.physics.units import Unit
    >>> Unit("meter", "m")
    m

    Other units are derived from base units:

    >>> import sympy.physics.units as u
    >>> cm = u.m/100
    >>> 100*u.cm
    m

    """
    is_positive = True    # make sqrt(m**2) --> m
    is_commutative = True
    is_number = False

    __slots__ = ["name", "abbrev"]

    def __new__(cls, name, abbrev, **assumptions):
        obj = AtomicExpr.__new__(cls, **assumptions)
        assert isinstance(name, str), repr(type(name))
        assert isinstance(abbrev, str), repr(type(abbrev))
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
        return (self.name, self.abbrev)

    @property
    def free_symbols(self):
        return set()

# Dimensionless

percent = percents = Rational(1, 100)
permille = permille = Rational(1, 1000)

ten = Rational(10)

yotta = ten**24
zetta = ten**21
exa = ten**18
peta = ten**15
tera = ten**12
giga = ten**9
mega = ten**6
kilo = ten**3
deca = ten**1
deci = ten**-1
centi = ten**-2
milli = ten**-3
micro = ten**-6
nano = ten**-9
pico = ten**-12
femto = ten**-15
atto = ten**-18
zepto = ten**-21
yocto = ten**-24

rad = radian = radians = 1
deg = degree = degrees = pi/180


# Base units

length = m = meter = meters = Unit('meter', 'm')
mass = kg = kilogram = kilograms = Unit('kilogram', 'kg')
time = s = second = seconds = Unit('second', 's')
current = A = ampere = amperes = Unit('ampere', 'A')
temperature = K = kelvin = kelvins = Unit('kelvin', 'K')
amount = mol = mole = moles = Unit('mole', 'mol')
luminosity = cd = candela = candelas = Unit('candela', 'cd')


# Derived units
volume = meter**3
frequency = Hz = hz = hertz = 1/s
force = N = newton = newtons = m*kg/s**2
energy = J = joule = joules = N*m
power = W = watt = watts = J/s
pressure = Pa = pa = pascal = pascals = N/m**2
charge = C = coulomb = coulombs = s*A
voltage = v = V = volt = volts = W/A
resistance = ohm = ohms = V/A
conductance = S = siemens = mho = mhos = A/V
capacitance = F = farad = farads = C/V
magnetic_flux = Wb = wb = weber = webers = J/A
magnetic_flux_density = T = tesla = teslas = V*s/m**2
inductance = H = henry = henrys = V*s/A
speed = m/s
acceleration = m/s**2
density = kg/m**3

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

avogadro_number = Rational('6.02214179') * 10**23
avogadro = avogadro_constant = avogadro_number / mol
boltzmann = Rational('1.3806505') * ten**-23 * J / K

gee = gees = Rational('9.80665') * m/s**2
atmosphere = atmospheres = atm = 101325 * pascal

kPa = kilo*Pa
bar = bars = 100*kPa
pound = pounds = 0.45359237 * kg * gee  # exact
psi = pound / inch ** 2
dHg0 = 13.5951  # approx value at 0 C
mmHg = dHg0 * 9.80665 * Pa
amu = amus = gram / avogadro / mol
mmu = mmus = gram / mol
quart = quarts = Rational(231, 4) * inch**3
eV = 1.602176487e-19 * J

# Other convenient units and magnitudes

ly = lightyear = lightyears = c*julian_year
au = astronomical_unit = astronomical_units = 149597870691*m


def find_unit(quantity):
    """
    Return a list of matching units names.
    if quantity is a string -- units containing the string `quantity`
    if quantity is a unit -- units having matching base units

    Examples
    ========

    >>> from sympy.physics import units as u
    >>> u.find_unit('charge')
    ['charge']
    >>> u.find_unit(u.charge)
    ['C', 'charge', 'coulomb', 'coulombs']
    >>> u.find_unit('volt')
    ['volt', 'volts', 'voltage']
    >>> u.find_unit(u.inch**3)[:5]
    ['l', 'cl', 'dl', 'ml', 'liter']
    """
    import sympy.physics.units as u
    rv = []
    if isinstance(quantity, str):
        rv = [i for i in dir(u) if quantity in i]
    else:
        units = quantity.as_coeff_Mul()[1]
        for i in dir(u):
            try:
                if units == eval('u.' + i).as_coeff_Mul()[1]:
                    rv.append(str(i))
            except:
                pass
    return sorted(rv, key=len)

# Delete this so it doesn't pollute the namespace
del Rational, pi
