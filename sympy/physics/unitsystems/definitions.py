from sympy.physics.units import temperature, acceleration

from sympy import pi
from sympy.physics.unitsystems import Quantity
from sympy.physics.unitsystems.dimensions import length, mass, force, energy, power, pressure, frequency, time, velocity, \
    impedance, voltage, conductance, capacitance, inductance, charge, magnetic_density, magnetic_flux, current, action
from sympy.physics.unitsystems.prefixes import kilo, milli, micro, nano, pico, deci, centi

#### UNITS ####

# Base units:
meter = m = Quantity("meter", length, 1, abbrev="m")

kilogram = kg = Quantity("kilogram", mass, kilo, abbrev="g")
second = s = Quantity("second", time, 1, abbrev="s")
ampere = A = Quantity("ampere", current, 1, abbrev='A')

# gram; used to define its prefixed units

g = gram = grams = Quantity("gram", mass, 1, abbrev="g")
mg = milligram = milligrams = Quantity("gram", mass, 1, "mg", "m")
ug = microgram = micrograms = Quantity("gram", mass, 1, "ug", "mu")

# derived units
newton = N = Quantity("newton", force, kilogram*meter/second**2, abbrev="N")
joule = J = Quantity("joule", energy, newton*meter, abbrev="J")
watt = W = Quantity("watt", power, joule/second, abbrev="W")
pascal = Pa = Quantity("pascal", pressure, newton/meter**2, abbrev="Pa")
hertz = hz = Hz = Quantity("hertz", frequency, 1, abbrev="Hz")

# MKSA extension:

# derived units

# coulomb
coulomb = C = Quantity("coulomb", charge, 1, abbrev='C')
# volt
volt = V = Quantity("volt", voltage, joule/coulomb, abbrev='V')
# ohm
ohm = Quantity("ohm", impedance, 10**3, abbrev='ohm')
# siemens
siemens = S = Quantity("siemens", conductance, 10**-3, abbrev='S')
# farad
farad = F = Quantity("farad", capacitance, coulomb/volt, abbrev='F')
# henry
henry = H = Quantity("henry", inductance, 10**3, abbrev='H')
# tesla
tesla = T = Quantity("tesla", magnetic_density, 1, abbrev='T')
# weber
weber = Wb = Quantity("weber", magnetic_flux, 1, abbrev='Wb')


# Common length units

km = kilometer = kilometers = Quantity("kilometer", length, kilo*meter, "km")
dm = decimeter = decimeters = Quantity("decimeter", length, deci*meter, "dm")
cm = centimeter = centimeters = Quantity("centimeter", length, centi*meter, "cm")
mm = millimeter = millimeters = Quantity("millimeter", length, milli*meter, "mm")
um = micrometer = micrometers = micron = microns = Quantity("micrometer", length, micro*meter, "um")
nm = nanometer = nanometers = Quantity("nanometer", length, nano*meter, "nn")
pm = picometer = picometers = Quantity("picometer", length, pico*meter, "pm")

ft = foot = feet = Quantity("foot", length, 0.3048*meter, "ft")
inch = inches = Quantity("inch", length, 0.0254*meter)
yd = yard = yards = Quantity("yard", length, 3*feet, "yd")
mi = mile = miles = Quantity("mile", length, 5280*feet)
nmi = nautical_mile = nautical_miles = Quantity("nautical_mile", length, 6076*feet)

# Common volume and area units

l = liter = liters = Quantity("liter", length**3, meter**3 / 1000)
dl = deciliter = deciliters = Quantity("deciliter", length**3, liter / 10)
cl = centiliter = centiliters = Quantity("centiliter", length**3, liter / 100)
ml = milliliter = milliliters = Quantity("milliliter", length**3, liter / 1000)

# Common time units

ms = millisecond = milliseconds = Quantity("millisecond", time, milli*second, "ms")
us = microsecond = microseconds = Quantity("microsecond", time, micro*second, "us")
ns = nanosecond = nanoseconds = Quantity("nanosecond", time, nano*second, "ns")
ps = picosecond = picoseconds = Quantity("picosecond", time, pico*second, "ps")

minute = minutes = Quantity("minute", time, 60*second)
h = hour = hours = Quantity("hour", time, 60*second)
day = days = Quantity("day", time, 24*hour)

anomalistic_year = anomalistic_years = Quantity("anomalistic_year", time, 365.259636*day)
sidereal_year = sidereal_years = Quantity("sidereal_year", time, 31558149.540)
tropical_year = tropical_years = Quantity("tropical_year", time, 365.24219*day)
common_year = common_years = Quantity("common_year", time, 365*day)
julian_year = julian_years = Quantity("julian_year", time, 365.25*day)
draconic_year = draconic_years = Quantity("draconic_year", time, 346.62*day)
gaussian_year = gaussian_years = Quantity("gaussian_year", time, 365.2568983*day)
full_moon_cycle = full_moon_cycles = Quantity("full_moon_cycle", time, 411.78443029*day)

year = years = tropical_year

#### CONSTANTS ####

# Newton constant
G = gravitational_constant = Quantity("gravitational_constant", length**3*mass**-1*time**-2, 6.67384e-11, abbrev="G")
# speed of light
c = speed_of_light = Quantity("speed_of_light", velocity, 299792458*meter/second, abbrev="c")
# Wave impedance of free space
Z0 = Quantity("WaveImpedence", impedance, 119.9169832*pi, abbrev='Z_0')
# Reduced Planck constant
hbar = Quantity("hbar", action, 1.05457266e-34, abbrev="hbar")
# Planck constant
planck = Quantity("planck", action, 2*pi*hbar.scale_factor, abbrev="h")
# Electronvolt
eV = Quantity("eV", energy, 1.60219e-19, abbrev="eV")
# Avogadro number
avogadro_number = Quantity("avogadro_number", 1, 6.022140857e23)
# Boltzmann constant
boltzmann = boltzmann_constant = Quantity("boltzmann_constant", energy/temperature, 1.38064852e-20)
# Atomic mass
atomic_mass_constant = Quantity("atomic_mass_constant", mass, 1.660539040e-24)
# gee
gee = gees = acceleration_due_to_gravity = Quantity("acceleration_due_to_gravity", acceleration, 9.80665)
#
