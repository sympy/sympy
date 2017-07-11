from sympy import pi, Rational, sqrt
from sympy.physics.units import Quantity
from sympy.physics.units.dimensions import length, mass, force, energy, power, pressure, frequency, time, velocity, \
    impedance, voltage, conductance, capacitance, inductance, charge, magnetic_density, magnetic_flux, current, action, \
    amount_of_substance, luminous_intensity, temperature, acceleration
from sympy.physics.units.prefixes import kilo, milli, micro, nano, pico, deci, centi

#### UNITS ####

# Dimensionless:

percent = percents = Quantity("percent", 1, Rational(1, 100))
permille = Quantity("permille", 1, Rational(1, 1000))

# Angular units (dimensionless)

rad = radian = radians = Quantity("radian", 1, 1)
deg = degree = degrees = Quantity("degree", 1, pi/180, "deg")
sr = steradian = steradians = Quantity("steradian", 1, 1, "sr")
mil = angular_mil = angular_mils = Quantity("angular_mil", 1, 2*pi/6400, "mil")

# Base units:

m = meter = meters = Quantity("meter", length, 1, abbrev="m")
kg = kilogram = kilograms = Quantity("kilogram", mass, kilo, abbrev="g")
s = second = seconds = Quantity("second", time, 1, abbrev="s")
A = ampere = amperes = Quantity("ampere", current, 1, abbrev='A')
K = kelvin = kelvins = Quantity('kelvin', temperature, 1, 'K')
mol = mole = moles = Quantity("mole", amount_of_substance, 1, "mol")
cd = candela = candelas = Quantity("candela", luminous_intensity, 1, "cd")

# gram; used to define its prefixed units

g = gram = grams = Quantity("gram", mass, 1, "g")
mg = milligram = milligrams = Quantity("milligram", mass, milli*gram, "mg")
ug = microgram = micrograms = Quantity("microgram", mass, micro*gram, "ug")

# derived units
newton = newtons = N = Quantity("newton", force, kilogram*meter/second**2, "N")
joule = joules = J = Quantity("joule", energy, newton*meter, "J")
watt = watts = W = Quantity("watt", power, joule/second, "W")
pascal = pascals = Pa = pa = Quantity("pascal", pressure, newton/meter**2, "Pa")
hertz = hz = Hz = Quantity("hertz", frequency, 1, "Hz")

# MKSA extension to MKS: derived units

coulomb = coulombs = C = Quantity("coulomb", charge, 1, abbrev='C')
volt = volts = v = V = Quantity("volt", voltage, joule/coulomb, abbrev='V')
ohm = ohms = Quantity("ohm", impedance, volt/ampere, abbrev='ohm')
siemens = S = mho = mhos = Quantity("siemens", conductance, ampere/volt, abbrev='S')
farad = farads = F = Quantity("farad", capacitance, coulomb/volt, abbrev='F')
henry = henrys = H = Quantity("henry", inductance, volt*second/ampere, abbrev='H')
tesla = teslas = T = Quantity("tesla", magnetic_density, volt*second/meter**2, abbrev='T')
weber = webers = Wb = wb = Quantity("weber", magnetic_flux, joule/ampere, abbrev='Wb')

# Other derived units:

optical_power = dioptre = D = Quantity("dioptre", 1/length, 1/meter)
lux = lx = Quantity("lux", luminous_intensity/length**2, steradian*candela/meter**2)

# Common length units

km = kilometer = kilometers = Quantity("kilometer", length, kilo*meter, "km")
dm = decimeter = decimeters = Quantity("decimeter", length, deci*meter, "dm")
cm = centimeter = centimeters = Quantity("centimeter", length, centi*meter, "cm")
mm = millimeter = millimeters = Quantity("millimeter", length, milli*meter, "mm")
um = micrometer = micrometers = micron = microns = Quantity("micrometer", length, micro*meter, "um")
nm = nanometer = nanometers = Quantity("nanometer", length, nano*meter, "nn")
pm = picometer = picometers = Quantity("picometer", length, pico*meter, "pm")

ft = foot = feet = Quantity("foot", length, Rational(3048, 10000)*meter, "ft")
inch = inches = Quantity("inch", length, foot/12)
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
h = hour = hours = Quantity("hour", time, 60*minute)
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
G = gravitational_constant = Quantity("gravitational_constant", length**3*mass**-1*time**-2, 6.67408e-11*m**3/(kg*s**2), abbrev="G")
# speed of light
c = speed_of_light = Quantity("speed_of_light", velocity, 299792458*meter/second, abbrev="c")
# Wave impedance of free space
Z0 = Quantity("WaveImpedence", impedance, 119.9169832*pi, abbrev='Z_0')
# Reduced Planck constant
hbar = Quantity("hbar", action, 1.05457266e-34*joule*second, abbrev="hbar")
# Planck constant
planck = Quantity("planck", action, 2*pi*hbar, abbrev="h")
# Electronvolt
eV = electronvolt = electronvolts = Quantity("electronvolt", energy, 1.60219e-19*joule, abbrev="eV")
# Avogadro number
avogadro_number = Quantity("avogadro_number", 1, 6.022140857e23)
# Avogadro constant
avogadro = avogadro_constant = Quantity("avogadro_constant", amount_of_substance**-1, avogadro_number / mol)
# Boltzmann constant
boltzmann = boltzmann_constant = Quantity("boltzmann_constant", energy/temperature, 1.38064852e-23*joule/kelvin)
# Atomic mass
amu = amus = atomic_mass_unit = atomic_mass_constant = Quantity("atomic_mass_constant", mass, 1.660539040e-24*gram)
# Molar gas constant
R = molar_gas_constant = Quantity("molar_gas_constant", energy/(temperature * amount_of_substance),
                                  8.3144598*joule/kelvin/mol, abbrev="R")
# Faraday constant
faraday_constant = Quantity("faraday_constant", charge/amount_of_substance, 96485.33289*C/mol)
# Josephson constant
josephson_constant = Quantity("josephson_constant", frequency/voltage, 483597.8525e9*hertz/V, abbrev="K_j")
# Von Klitzing constant
von_klitzing_constant = Quantity("von_klitzing_constant", voltage/current, 25812.8074555*ohm, abbrev="R_k")
# Acceleration due to gravity (on the Earth surface)
gee = gees = acceleration_due_to_gravity = Quantity("acceleration_due_to_gravity", acceleration, 9.80665*meter/second**2, "g")
# magnetic constant:
u0 = magnetic_constant = Quantity("magnetic_constant", force/current**2, 4*pi/10**7 * newton/ampere**2)
# electric constat:
e0 = electric_constant = vacuum_permittivity = Quantity("vacuum_permittivity", capacitance/length, 1/(u0 * c**2))
# vacuum impedance:
Z0 = vacuum_impedance = Quantity("vacuum_impedance", impedance, u0 * c)
# Coulomb's constant:
coulomb_constant = electric_force_constant = Quantity("coulomb_constant", force*length**2/charge**2, 1/(4*pi*vacuum_permittivity), "k_e")

atmosphere = atmospheres = atm = Quantity("atmosphere", pressure, 101325 * pascal, "atm")

kPa = kilopascal = Quantity("kilopascal", pressure, kilo*Pa, "kPa")
bar = bars = Quantity("bar", pressure, 100*kPa, "bar")
pound = pounds = Quantity("pound", mass, 0.45359237 * kg)  # exact
psi = Quantity("psi", pressure, pound * gee / inch ** 2)
dHg0 = 13.5951  # approx value at 0 C
mmHg = torr = Quantity("mmHg", pressure, dHg0 * acceleration_due_to_gravity * kilogram / meter**2)
mmu = mmus = milli_mass_unit = Quantity("milli_mass_unit", mass, atomic_mass_unit/1000)
quart = quarts = Quantity("quart", length**3, Rational(231, 4) * inch**3)

# Other convenient units and magnitudes

ly = lightyear = lightyears = Quantity("lightyear", length, speed_of_light*julian_year, "ly")
au = astronomical_unit = astronomical_units = Quantity("astronomical_unit", length, 149597870691*meter, "AU")

# Planck units:
planck_mass = Quantity("planck_mass", mass, sqrt(hbar*speed_of_light/G), "m_P")
planck_time = Quantity("planck_time", time, sqrt(hbar*G/speed_of_light**5), "t_P")
planck_temperature = Quantity("planck_temperature", temperature, sqrt(hbar*speed_of_light**5/G/boltzmann**2), "T_P")
planck_length = Quantity("planck_length", length, sqrt(hbar*G/speed_of_light**3), "l_P")
# TODO: add more from https://en.wikipedia.org/wiki/Planck_units
