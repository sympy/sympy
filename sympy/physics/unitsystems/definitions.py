from sympy import pi
from sympy import S as S_singleton
from sympy.physics.unitsystems import Quantity
from sympy.physics.unitsystems import Unit
from sympy.physics.unitsystems.dimensions import length, mass, force, energy, power, pressure, frequency, time, velocity, \
    impedance, voltage, conductance, capacitance, inductance, charge, magnetic_density, magnetic_flux, current, action
from sympy.physics.unitsystems.prefixes import PREFIXES


#### UNITS ####

## Base units:
meter = m = Unit("meter", length, 1, abbrev="m")

kilogram = kg = Unit("kilogram", mass, 1, abbrev="g", prefix=PREFIXES["k"])
second = s = Unit("second", time, 1, abbrev="s")
ampere = A = Unit("ampere", current, 1, abbrev='A')

# gram; used to define its prefixed units

g = gram = grams = Unit("gram", mass, 1, abbrev="g")
mg = milligram = milligrams = Unit("gram", mass, 1, "mg", "m")
ug = microgram = micrograms = Unit("gram", mass, 1, "ug", "mu")

# derived units
joule = J = Unit("joule", energy, 10**3, abbrev="J")
newton = N = Unit("newton", force, 10**3, abbrev="N")
watt = W = Unit("watt", power, 10**3, abbrev="W")
pascal = Pa = Unit("pascal", pressure, 10**3, abbrev="Pa")
hertz = hz = Hz = Unit("hertz", frequency, 1, abbrev="Hz")

# MKSA extension:

# derived units
volt = V = Unit("volt", voltage, 10**3, abbrev='V')
ohm = Unit("ohm", impedance, 10**3, abbrev='ohm')
# siemens
siemens = S = Unit("siemens", conductance, 10**-3, abbrev='S')
# farad
farad = F = Unit("farad", capacitance, 10**-3, abbrev='F')
# henry
henry = H = Unit("henry", inductance, 10**3, abbrev='H')
# coulomb
coulomb = C = Unit("coulomb", charge, 1, abbrev='C')
# tesla
tesla = T = Unit("tesla", magnetic_density, 1, abbrev='T')
# weber
weber = Wb = Unit("weber", magnetic_flux, 1, abbrev='Wb')


# Common length units

km = kilometer = kilometers = Unit("kilometer", length, 1, "km", "k")
dm = decimeter = decimeters = Unit("decimeter", length, 1, "dm", "d")
cm = centimeter = centimeters = Unit("centimeter", length, 1, "cm", "c")
mm = millimeter = millimeters = Unit("millimeter", length, 1, "mm", "m")
um = micrometer = micrometers = micron = microns = Unit("micrometer", length, 1, "um", "mu")
nm = nanometer = nanometers = Unit("nanometer", length, 1, "nn", "n")
pm = picometer = picometers = Unit("picometer", length, 1, "pm", "p")

ft = foot = feet = Quantity("foot", length, 0.3048, "ft")
inch = inches = Quantity("inch", length, 0.0254)
yd = yard = yards = Quantity("yard", length, 3*0.3048, "yd")
mi = mile = miles = Quantity("mile", length, 5280*0.3048)
nmi = nautical_mile = nautical_miles = Quantity("nautical_mile", length, 6076*0.3048)

# Common volume and area units

l = liter = liters = Quantity("liter", length**3, S_singleton.One / 1000)
dl = deciliter = deciliters = Quantity("deciliter", length**3, S_singleton.One / 10**4)
cl = centiliter = centiliters = Quantity("centiliter", length**3, S_singleton.One / 10**5)
ml = milliliter = milliliters = Quantity("milliliter", length**3, S_singleton.One / 10**6)

# Common time units

ms = millisecond = milliseconds = Unit("millisecond", time, 1, "m")
us = microsecond = microseconds = Unit("microsecond", time, 1, "mu")
ns = nanosecond = nanoseconds = Unit("nanosecond", time, 1, "n")
ps = picosecond = picoseconds = Unit("picosecond", time, 1, "p")

minute = Quantity("minute", time, 60)
hour = Quantity("hour", time, 3600)
day = Quantity("day", time, 86400)

anomalistic_year = anomalistic_years = Quantity("anomalistic_year", time, 365.259636*86400)
sidereal_year = sidereal_years = Quantity("sidereal_year", time, 31558149.540)
tropical_year = tropical_years = Quantity("tropical_year", time, 365.24219*86400)
common_year = common_years = Quantity("common_year", time, 365*86400)
julian_year = julian_years = Quantity("julian_year", time, 365.25*86400)
draconic_year = draconic_years = Quantity("draconic_year", time, 346.62*86400)
gaussian_year = gaussian_years = Quantity("gaussian_year", time, 365.2568983*86400)
full_moon_cycle = full_moon_cycles = Quantity("full_moon_cycle", time, 411.78443029*86400)

year = years = tropical_year

#### CONSTANTS ####

# Newton constant
G = gravitational_constant = Quantity("gravitational_constant", length**3*mass**-1*time**-2, 6.67384e-11, abbrev="G")
# speed of light
c = speed_of_light = Quantity("speed_of_light", velocity, 299792458, abbrev="c")
# Wave impedance of free space
Z0 = Quantity("WaveImpedence", impedance, 119.9169832*pi, abbrev='Z_0')
# Reduced Planck constant
hbar = Quantity("hbar", action, 1.05457266e-34, abbrev="hbar")
# Planck constant
planck = Quantity("planck", action, 2*pi*hbar.scale_factor, abbrev="h")
# Electronvolt
eV = Quantity("eV", energy, 1.60219e-19, abbrev="eV")
