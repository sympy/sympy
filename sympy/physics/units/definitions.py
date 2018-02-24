from sympy import Rational, pi, sqrt, sympify, S
from sympy.physics.units.quantities import Quantity
from sympy.physics.units.dimensions import (
    acceleration, action, amount_of_substance, capacitance, charge,
    conductance, current, energy, force, frequency, information, impedance, inductance,
    length, luminous_intensity, magnetic_density, magnetic_flux, mass, power,
    pressure, temperature, time, velocity, voltage)
from sympy.physics.units.dimensions import dimsys_default, Dimension
from sympy.physics.units.prefixes import (
    centi, deci, kilo, micro, milli, nano, pico,
    kibi, mebi, gibi, tebi, pebi, exbi)

One = S.One

#### UNITS ####

# Dimensionless:

percent = percents = Quantity("percent")
percent._set_dimension(One)
percent._set_scale_factor(Rational(1, 100))

permille = Quantity("permille")
permille._set_dimension(One)
permille._set_scale_factor(Rational(1, 1000))


# Angular units (dimensionless)

rad = radian = radians = Quantity("radian")
radian._set_dimension(One)
radian._set_scale_factor(One)

deg = degree = degrees = Quantity("degree", abbrev="deg")
degree._set_dimension(One)
degree._set_scale_factor(pi/180)

sr = steradian = steradians = Quantity("steradian", abbrev="sr")
steradian._set_dimension(One)
steradian._set_scale_factor(One)

mil = angular_mil = angular_mils = Quantity("angular_mil", abbrev="mil")
angular_mil._set_dimension(One)
angular_mil._set_scale_factor(2*pi/6400)


# Base units:

m = meter = meters = Quantity("meter", abbrev="m")
meter._set_dimension(length)
meter._set_scale_factor(One)

kg = kilogram = kilograms = Quantity("kilogram", abbrev="kg")
kilogram._set_dimension(mass)
kilogram._set_scale_factor(kilo)

s = second = seconds = Quantity("second", abbrev="s")
second._set_dimension(time)
second._set_scale_factor(One)

A = ampere = amperes = Quantity("ampere", abbrev='A')
ampere._set_dimension(current)
ampere._set_scale_factor(One)

K = kelvin = kelvins = Quantity("kelvin", abbrev='K')
kelvin._set_dimension(temperature)
kelvin._set_scale_factor(One)

mol = mole = moles = Quantity("mole", abbrev="mol")
mole._set_dimension(amount_of_substance)
mole._set_scale_factor(One)

cd = candela = candelas = Quantity("candela", abbrev="cd")
candela._set_dimension(luminous_intensity)
candela._set_scale_factor(One)


# gram; used to define its prefixed units

g = gram = grams = Quantity("gram", abbrev="g")
gram._set_dimension(mass)
gram._set_scale_factor(One)

mg = milligram = milligrams = Quantity("milligram", abbrev="mg")
milligram._set_dimension(mass)
milligram._set_scale_factor(milli*gram)

ug = microgram = micrograms = Quantity("microgram", abbrev="ug")
microgram._set_dimension(mass)
microgram._set_scale_factor(micro*gram)


# derived units
newton = newtons = N = Quantity("newton", abbrev="N")
newton._set_dimension(force)
newton._set_scale_factor(kilogram*meter/second**2)

joule = joules = J = Quantity("joule", abbrev="J")
joule._set_dimension(energy)
joule._set_scale_factor(newton*meter)

watt = watts = W = Quantity("watt", abbrev="W")
watt._set_dimension(power)
watt._set_scale_factor(joule/second)

pascal = pascals = Pa = pa = Quantity("pascal", abbrev="Pa")
pascal._set_dimension(pressure)
pascal._set_scale_factor(newton/meter**2)

hertz = hz = Hz = Quantity("hertz", abbrev="Hz")
hertz._set_dimension(frequency)
hertz._set_scale_factor(One)


# MKSA extension to MKS: derived units

coulomb = coulombs = C = Quantity("coulomb", abbrev='C')
coulomb._set_dimension(charge)
coulomb._set_scale_factor(One)

volt = volts = v = V = Quantity("volt", abbrev='V')
volt._set_dimension(voltage)
volt._set_scale_factor(joule/coulomb)

ohm = ohms = Quantity("ohm", abbrev='ohm')
ohm._set_dimension(impedance)
ohm._set_scale_factor(volt/ampere)

siemens = S = mho = mhos = Quantity("siemens", abbrev='S')
siemens._set_dimension(conductance)
siemens._set_scale_factor(ampere/volt)

farad = farads = F = Quantity("farad", abbrev='F')
farad._set_dimension(capacitance)
farad._set_scale_factor(coulomb/volt)

henry = henrys = H = Quantity("henry", abbrev='H')
henry._set_dimension(inductance)
henry._set_scale_factor(volt*second/ampere)

tesla = teslas = T = Quantity("tesla", abbrev='T')
tesla._set_dimension(magnetic_density)
tesla._set_scale_factor(volt*second/meter**2)

weber = webers = Wb = wb = Quantity("weber", abbrev='Wb')
weber._set_dimension(magnetic_flux)
weber._set_scale_factor(joule/ampere)


# Other derived units:

optical_power = dioptre = D = Quantity("dioptre")
dioptre._set_dimension(1/length)
dioptre._set_scale_factor(1/meter)

lux = lx = Quantity("lux")
lux._set_dimension(luminous_intensity/length**2)
lux._set_scale_factor(steradian*candela/meter**2)

# katal is the SI unit of catalytic activity
katal = kat = Quantity("katal")
katal._set_dimension(amount_of_substance/time)
katal._set_scale_factor(mol/second)

# gray is the SI unit of absorbed dose
gray = Gy = Quantity("gray")
gray._set_dimension(energy/mass)
gray._set_scale_factor(meter**2/second**2)

# becquerel is the SI unit of radioactivity
becquerel = Bq = Quantity("becquerel")
becquerel._set_dimension(1/time)
becquerel._set_scale_factor(1/second)


# Common length units

km = kilometer = kilometers = Quantity("kilometer", abbrev="km")
kilometer._set_dimension(length)
kilometer._set_scale_factor(kilo*meter)

dm = decimeter = decimeters = Quantity("decimeter", abbrev="dm")
decimeter._set_dimension(length)
decimeter._set_scale_factor(deci*meter)

cm = centimeter = centimeters = Quantity("centimeter", abbrev="cm")
centimeter._set_dimension(length)
centimeter._set_scale_factor(centi*meter)

mm = millimeter = millimeters = Quantity("millimeter", abbrev="mm")
millimeter._set_dimension(length)
millimeter._set_scale_factor(milli*meter)

um = micrometer = micrometers = micron = microns = Quantity("micrometer", abbrev="um")
micrometer._set_dimension(length)
micrometer._set_scale_factor(micro*meter)

nm = nanometer = nanometers = Quantity("nanometer", abbrev="nn")
nanometer._set_dimension(length)
nanometer._set_scale_factor(nano*meter)

pm = picometer = picometers = Quantity("picometer", abbrev="pm")
picometer._set_dimension(length)
picometer._set_scale_factor(pico*meter)


ft = foot = feet = Quantity("foot", abbrev="ft")
foot._set_dimension(length)
foot._set_scale_factor(Rational(3048, 10000)*meter)

inch = inches = Quantity("inch")
inch._set_dimension(length)
inch._set_scale_factor(foot/12)

yd = yard = yards = Quantity("yard", abbrev="yd")
yard._set_dimension(length)
yard._set_scale_factor(3*feet)

mi = mile = miles = Quantity("mile")
mile._set_dimension(length)
mile._set_scale_factor(5280*feet)

nmi = nautical_mile = nautical_miles = Quantity("nautical_mile")
nautical_mile._set_dimension(length)
nautical_mile._set_scale_factor(6076*feet)


# Common volume and area units

l = liter = liters = Quantity("liter")
liter._set_dimension(length**3)
liter._set_scale_factor(meter**3 / 1000)

dl = deciliter = deciliters = Quantity("deciliter")
deciliter._set_dimension(length**3)
deciliter._set_scale_factor(liter / 10)

cl = centiliter = centiliters = Quantity("centiliter")
centiliter._set_dimension(length**3)
centiliter._set_scale_factor(liter / 100)

ml = milliliter = milliliters = Quantity("milliliter")
milliliter._set_dimension(length**3)
milliliter._set_scale_factor(liter / 1000)


# Common time units

ms = millisecond = milliseconds = Quantity("millisecond", abbrev="ms")
millisecond._set_dimension(time)
millisecond._set_scale_factor(milli*second)

us = microsecond = microseconds = Quantity("microsecond", abbrev="us")
microsecond._set_dimension(time)
microsecond._set_scale_factor(micro*second)

ns = nanosecond = nanoseconds = Quantity("nanosecond", abbrev="ns")
nanosecond._set_dimension(time)
nanosecond._set_scale_factor(nano*second)

ps = picosecond = picoseconds = Quantity("picosecond", abbrev="ps")
picosecond._set_dimension(time)
picosecond._set_scale_factor(pico*second)


minute = minutes = Quantity("minute")
minute._set_dimension(time)
minute._set_scale_factor(60*second)

h = hour = hours = Quantity("hour")
hour._set_dimension(time)
hour._set_scale_factor(60*minute)

day = days = Quantity("day")
day._set_dimension(time)
day._set_scale_factor(24*hour)


anomalistic_year = anomalistic_years = Quantity("anomalistic_year")
anomalistic_year._set_dimension(time)
anomalistic_year._set_scale_factor(365.259636*day)

sidereal_year = sidereal_years = Quantity("sidereal_year")
sidereal_year._set_dimension(time)
sidereal_year._set_scale_factor(31558149.540)

tropical_year = tropical_years = Quantity("tropical_year")
tropical_year._set_dimension(time)
tropical_year._set_scale_factor(365.24219*day)

common_year = common_years = Quantity("common_year")
common_year._set_dimension(time)
common_year._set_scale_factor(365*day)

julian_year = julian_years = Quantity("julian_year")
julian_year._set_dimension(time)
julian_year._set_scale_factor(365.25*day)

draconic_year = draconic_years = Quantity("draconic_year")
draconic_year._set_dimension(time)
draconic_year._set_scale_factor(346.62*day)

gaussian_year = gaussian_years = Quantity("gaussian_year")
gaussian_year._set_dimension(time)
gaussian_year._set_scale_factor(365.2568983*day)

full_moon_cycle = full_moon_cycles = Quantity("full_moon_cycle")
full_moon_cycle._set_dimension(time)
full_moon_cycle._set_scale_factor(411.78443029*day)


year = years = tropical_year

#### CONSTANTS ####

# Newton constant
G = gravitational_constant = Quantity("gravitational_constant", abbrev="G")
gravitational_constant._set_dimension(length**3*mass**-1*time**-2)
gravitational_constant._set_scale_factor(6.67408e-11*m**3/(kg*s**2))

# speed of light
c = speed_of_light = Quantity("speed_of_light", abbrev="c")
speed_of_light._set_dimension(velocity)
speed_of_light._set_scale_factor(299792458*meter/second)

# Wave impedance of free space
Z0 = wave_impendence = Quantity("wave_impendence", abbrev='Z_0')
wave_impendence._set_dimension(impedance)
wave_impendence._set_scale_factor(119.9169832*pi)

# Reduced Planck constant
hbar = Quantity("hbar", abbrev="hbar")
hbar._set_dimension(action)
hbar._set_scale_factor(1.05457266e-34*joule*second)

# Planck constant
planck = Quantity("planck", abbrev="h")
planck._set_dimension(action)
planck._set_scale_factor(2*pi*hbar)

# Electronvolt
eV = electronvolt = electronvolts = Quantity("electronvolt", abbrev="eV")
electronvolt._set_dimension(energy)
electronvolt._set_scale_factor(1.60219e-19*joule)

# Avogadro number
avogadro_number = Quantity("avogadro_number")
avogadro_number._set_dimension(One)
avogadro_number._set_scale_factor(6.022140857e23)

# Avogadro constant
avogadro = avogadro_constant = Quantity("avogadro_constant")
avogadro_constant._set_dimension(amount_of_substance**-1)
avogadro_constant._set_scale_factor(avogadro_number / mol)

# Boltzmann constant
boltzmann = boltzmann_constant = Quantity("boltzmann_constant")
boltzmann_constant._set_dimension(energy/temperature)
boltzmann_constant._set_scale_factor(1.38064852e-23*joule/kelvin)

# Stefan-Boltzmann constant
stefan = stefan_boltzmann_constant = Quantity("stefan_boltzmann_constant")
stefan_boltzmann_constant._set_dimension(energy*time**-1*length**-2*temperature**-4)
stefan_boltzmann_constant._set_scale_factor(5.670367e-8*joule/(s*m**2*kelvin**4))

# Atomic mass
amu = amus = atomic_mass_unit = atomic_mass_constant = Quantity("atomic_mass_constant")
atomic_mass_constant._set_dimension(mass)
atomic_mass_constant._set_scale_factor(1.660539040e-24*gram)

# Molar gas constant
R = molar_gas_constant = Quantity("molar_gas_constant", abbrev="R")
molar_gas_constant._set_dimension(energy/(temperature * amount_of_substance))
molar_gas_constant._set_scale_factor(8.3144598*joule/kelvin/mol)

# Faraday constant
faraday_constant = Quantity("faraday_constant")
faraday_constant._set_dimension(charge/amount_of_substance)
faraday_constant._set_scale_factor(96485.33289*C/mol)

# Josephson constant
josephson_constant = Quantity("josephson_constant", abbrev="K_j")
josephson_constant._set_dimension(frequency/voltage)
josephson_constant._set_scale_factor(483597.8525e9*hertz/V)

# Von Klitzing constant
von_klitzing_constant = Quantity("von_klitzing_constant", abbrev="R_k")
von_klitzing_constant._set_dimension(voltage/current)
von_klitzing_constant._set_scale_factor(25812.8074555*ohm)

# Acceleration due to gravity (on the Earth surface)
gee = gees = acceleration_due_to_gravity = Quantity("acceleration_due_to_gravity", abbrev="g")
acceleration_due_to_gravity._set_dimension(acceleration)
acceleration_due_to_gravity._set_scale_factor(9.80665*meter/second**2)

# magnetic constant:
u0 = magnetic_constant = Quantity("magnetic_constant")
magnetic_constant._set_dimension(force/current**2)
magnetic_constant._set_scale_factor(4*pi/10**7 * newton/ampere**2)

# electric constat:
e0 = electric_constant = vacuum_permittivity = Quantity("vacuum_permittivity")
vacuum_permittivity._set_dimension(capacitance/length)
vacuum_permittivity._set_scale_factor(1/(u0 * c**2))

# vacuum impedance:
Z0 = vacuum_impedance = Quantity("vacuum_impedance")
vacuum_impedance._set_dimension(impedance)
vacuum_impedance._set_scale_factor(u0 * c)

# Coulomb's constant:
coulomb_constant = coulombs_constant = electric_force_constant = Quantity("coulomb_constant", abbrev="k_e")
coulomb_constant._set_dimension(force*length**2/charge**2)
coulomb_constant._set_scale_factor(1/(4*pi*vacuum_permittivity))


atmosphere = atmospheres = atm = Quantity("atmosphere", abbrev="atm")
atmosphere._set_dimension(pressure)
atmosphere._set_scale_factor(101325 * pascal)


kPa = kilopascal = Quantity("kilopascal", abbrev="kPa")
kilopascal._set_dimension(pressure)
kilopascal._set_scale_factor(kilo*Pa)

bar = bars = Quantity("bar", abbrev="bar")
bar._set_dimension(pressure)
bar._set_scale_factor(100*kPa)

pound = pounds = Quantity("pound")  # exact
pound._set_dimension(mass)
pound._set_scale_factor(0.45359237 * kg)

psi = Quantity("psi")
psi._set_dimension(pressure)
psi._set_scale_factor(pound * gee / inch ** 2)

dHg0 = 13.5951  # approx value at 0 C
mmHg = torr = Quantity("mmHg")
mmHg._set_dimension(pressure)
mmHg._set_scale_factor(dHg0 * acceleration_due_to_gravity * kilogram / meter**2)

mmu = mmus = milli_mass_unit = Quantity("milli_mass_unit")
milli_mass_unit._set_dimension(mass)
milli_mass_unit._set_scale_factor(atomic_mass_unit/1000)

quart = quarts = Quantity("quart")
quart._set_dimension(length**3)
quart._set_scale_factor(Rational(231, 4) * inch**3)


# Other convenient units and magnitudes

ly = lightyear = lightyears = Quantity("lightyear", abbrev="ly")
lightyear._set_dimension(length)
lightyear._set_scale_factor(speed_of_light*julian_year)

au = astronomical_unit = astronomical_units = Quantity("astronomical_unit", abbrev="AU")
astronomical_unit._set_dimension(length)
astronomical_unit._set_scale_factor(149597870691*meter)


# Fundamental Planck units:
planck_mass = Quantity("planck_mass", abbrev="m_P")
planck_mass._set_dimension(mass)
planck_mass._set_scale_factor(sqrt(hbar*speed_of_light/G))

planck_time = Quantity("planck_time", abbrev="t_P")
planck_time._set_dimension(time)
planck_time._set_scale_factor(sqrt(hbar*G/speed_of_light**5))

planck_temperature = Quantity("planck_temperature", abbrev="T_P")
planck_temperature._set_dimension(temperature)
planck_temperature._set_scale_factor(sqrt(hbar*speed_of_light**5/G/boltzmann**2))

planck_length = Quantity("planck_length", abbrev="l_P")
planck_length._set_dimension(length)
planck_length._set_scale_factor(sqrt(hbar*G/speed_of_light**3))

planck_charge = Quantity("planck_charge", abbrev="q_P")
planck_charge._set_dimension(charge)
planck_charge._set_scale_factor(sqrt(4*pi*electric_constant*hbar*speed_of_light))


# Derived Planck units:
planck_area = Quantity("planck_area")
planck_area._set_dimension(length**2)
planck_area._set_scale_factor(planck_length**2)

planck_volume = Quantity("planck_volume")
planck_volume._set_dimension(length**3)
planck_volume._set_scale_factor(planck_length**3)

planck_momentum = Quantity("planck_momentum")
planck_momentum._set_dimension(mass*velocity)
planck_momentum._set_scale_factor(planck_mass * speed_of_light)

planck_energy = Quantity("planck_energy", abbrev="E_P")
planck_energy._set_dimension(energy)
planck_energy._set_scale_factor(planck_mass * speed_of_light**2)

planck_force = Quantity("planck_force", abbrev="F_P")
planck_force._set_dimension(force)
planck_force._set_scale_factor(planck_energy / planck_length)

planck_power = Quantity("planck_power", abbrev="P_P")
planck_power._set_dimension(power)
planck_power._set_scale_factor(planck_energy / planck_time)

planck_density = Quantity("planck_density", abbrev="rho_P")
planck_density._set_dimension(mass/length**3)
planck_density._set_scale_factor(planck_mass / planck_length**3)

planck_energy_density = Quantity("planck_energy_density", abbrev="rho^E_P")
planck_energy_density._set_dimension(energy / length**3)
planck_energy_density._set_scale_factor(planck_energy / planck_length**3)

planck_intensity = Quantity("planck_intensity", abbrev="I_P")
planck_intensity._set_dimension(mass * time**(-3))
planck_intensity._set_scale_factor(planck_energy_density * speed_of_light)

planck_angular_frequency = Quantity("planck_angular_frequency", abbrev="omega_P")
planck_angular_frequency._set_dimension(1 / time)
planck_angular_frequency._set_scale_factor(1 / planck_time)

planck_pressure = Quantity("planck_pressure", abbrev="p_P")
planck_pressure._set_dimension(pressure)
planck_pressure._set_scale_factor(planck_force / planck_length**2)

planck_current = Quantity("planck_current", abbrev="I_P")
planck_current._set_dimension(current)
planck_current._set_scale_factor(planck_charge / planck_time)

planck_voltage = Quantity("planck_voltage", abbrev="V_P")
planck_voltage._set_dimension(voltage)
planck_voltage._set_scale_factor(planck_energy / planck_charge)

planck_impedance = Quantity("planck_impedance", abbrev="Z_P")
planck_impedance._set_dimension(impedance)
planck_impedance._set_scale_factor(planck_voltage / planck_current)

planck_acceleration = Quantity("planck_acceleration", abbrev="a_P")
planck_acceleration._set_dimension(acceleration)
planck_acceleration._set_scale_factor(speed_of_light / planck_time)


# Information theory units:
bit = bits = Quantity("bit")
bit._set_dimension(information)
bit._set_scale_factor(One)

byte = bytes = Quantity("byte")
byte._set_dimension(information)
byte._set_scale_factor(8*bit)

kibibyte = kibibytes = Quantity("kibibyte")
kibibyte._set_dimension(information)
kibibyte._set_scale_factor(kibi*byte)

mebibyte = mebibytes = Quantity("mebibyte")
mebibyte._set_dimension(information)
mebibyte._set_scale_factor(mebi*byte)

gibibyte = gibibytes = Quantity("gibibyte")
gibibyte._set_dimension(information)
gibibyte._set_scale_factor(gibi*byte)

tebibyte = tebibytes = Quantity("tebibyte")
tebibyte._set_dimension(information)
tebibyte._set_scale_factor(tebi*byte)

pebibyte = pebibytes = Quantity("pebibyte")
pebibyte._set_dimension(information)
pebibyte._set_scale_factor(pebi*byte)

exbibyte = exbibytes = Quantity("exbibyte")
exbibyte._set_dimension(information)
exbibyte._set_scale_factor(exbi*byte)


# check that scale factors are the right SI dimensions:
for _scale_factor, _dimension in zip(
        Quantity.SI_quantity_scale_factors.values(),
        Quantity.SI_quantity_dimension_map.values()):
    dimex = Quantity.get_dimensional_expr(_scale_factor)
    if dimex != 1:
        if not dimsys_default.equivalent_dims(_dimension, Dimension(dimex)):
            raise ValueError("quantity value and dimension mismatch")
del _scale_factor, _dimension
