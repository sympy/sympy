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

# Dimensional representations for the SI units:
SI_quantity_dimension_map = {}
# Scale factors in SI units:
SI_quantity_scale_factors = {}

#### UNITS ####

# Dimensionless:

percent = percents = Quantity("percent")
SI_quantity_dimension_map[percent] = One
SI_quantity_scale_factors[percent] = Rational(1, 100)

permille = Quantity("permille")
SI_quantity_dimension_map[permille] = One
SI_quantity_scale_factors[permille] = Rational(1, 1000)


# Angular units (dimensionless)

rad = radian = radians = Quantity("radian")
SI_quantity_dimension_map[radian] = One
SI_quantity_scale_factors[radian] = One

deg = degree = degrees = Quantity("degree", abbrev="deg")
SI_quantity_dimension_map[degree] = One
SI_quantity_scale_factors[degree] = pi/180

sr = steradian = steradians = Quantity("steradian", abbrev="sr")
SI_quantity_dimension_map[steradian] = One
SI_quantity_scale_factors[steradian] = One

mil = angular_mil = angular_mils = Quantity("angular_mil", abbrev="mil")
SI_quantity_dimension_map[angular_mil] = One
SI_quantity_scale_factors[angular_mil] = 2*pi/6400


# Base units:

m = meter = meters = Quantity("meter", abbrev="m")
SI_quantity_dimension_map[meter] = length
SI_quantity_scale_factors[meter] = One

kg = kilogram = kilograms = Quantity("kilogram", abbrev="kg")
SI_quantity_dimension_map[kilogram] = mass
SI_quantity_scale_factors[kilogram] = kilo

s = second = seconds = Quantity("second", abbrev="s")
SI_quantity_dimension_map[second] = time
SI_quantity_scale_factors[second] = One

A = ampere = amperes = Quantity("ampere", abbrev='A')
SI_quantity_dimension_map[ampere] = current
SI_quantity_scale_factors[ampere] = One

K = kelvin = kelvins = Quantity("kelvin", abbrev='K')
SI_quantity_dimension_map[kelvin] = temperature
SI_quantity_scale_factors[kelvin] = One

mol = mole = moles = Quantity("mole", abbrev="mol")
SI_quantity_dimension_map[mole] = amount_of_substance
SI_quantity_scale_factors[mole] = One

cd = candela = candelas = Quantity("candela", abbrev="cd")
SI_quantity_dimension_map[candela] = luminous_intensity
SI_quantity_scale_factors[candela] = One


# gram; used to define its prefixed units

g = gram = grams = Quantity("gram", abbrev="g")
SI_quantity_dimension_map[gram] = mass
SI_quantity_scale_factors[gram] = One

mg = milligram = milligrams = Quantity("milligram", abbrev="mg")
SI_quantity_dimension_map[milligram] = mass
SI_quantity_scale_factors[milligram] = milli*gram

ug = microgram = micrograms = Quantity("microgram", abbrev="ug")
SI_quantity_dimension_map[microgram] = mass
SI_quantity_scale_factors[microgram] = micro*gram


# derived units
newton = newtons = N = Quantity("newton", abbrev="N")
SI_quantity_dimension_map[newton] = force
SI_quantity_scale_factors[newton] = kilogram*meter/second**2

joule = joules = J = Quantity("joule", abbrev="J")
SI_quantity_dimension_map[joule] = energy
SI_quantity_scale_factors[joule] = newton*meter

watt = watts = W = Quantity("watt", abbrev="W")
SI_quantity_dimension_map[watt] = power
SI_quantity_scale_factors[watt] = joule/second

pascal = pascals = Pa = pa = Quantity("pascal", abbrev="Pa")
SI_quantity_dimension_map[pascal] = pressure
SI_quantity_scale_factors[pascal] = newton/meter**2

hertz = hz = Hz = Quantity("hertz", abbrev="Hz")
SI_quantity_dimension_map[hertz] = frequency
SI_quantity_scale_factors[hertz] = One


# MKSA extension to MKS: derived units

coulomb = coulombs = C = Quantity("coulomb", abbrev='C')
SI_quantity_dimension_map[coulomb] = charge
SI_quantity_scale_factors[coulomb] = One

volt = volts = v = V = Quantity("volt", abbrev='V')
SI_quantity_dimension_map[volt] = voltage
SI_quantity_scale_factors[volt] = joule/coulomb

ohm = ohms = Quantity("ohm", abbrev='ohm')
SI_quantity_dimension_map[ohm] = impedance
SI_quantity_scale_factors[ohm] = volt/ampere

siemens = S = mho = mhos = Quantity("siemens", abbrev='S')
SI_quantity_dimension_map[siemens] = conductance
SI_quantity_scale_factors[siemens] = ampere/volt

farad = farads = F = Quantity("farad", abbrev='F')
SI_quantity_dimension_map[farad] = capacitance
SI_quantity_scale_factors[farad] = coulomb/volt

henry = henrys = H = Quantity("henry", abbrev='H')
SI_quantity_dimension_map[henry] = inductance
SI_quantity_scale_factors[henry] = volt*second/ampere

tesla = teslas = T = Quantity("tesla", abbrev='T')
SI_quantity_dimension_map[tesla] = magnetic_density
SI_quantity_scale_factors[tesla] = volt*second/meter**2

weber = webers = Wb = wb = Quantity("weber", abbrev='Wb')
SI_quantity_dimension_map[weber] = magnetic_flux
SI_quantity_scale_factors[weber] = joule/ampere


# Other derived units:

optical_power = dioptre = D = Quantity("dioptre")
SI_quantity_dimension_map[dioptre] = 1/length
SI_quantity_scale_factors[dioptre] = 1/meter

lux = lx = Quantity("lux")
SI_quantity_dimension_map[lux] = luminous_intensity/length**2
SI_quantity_scale_factors[lux] = steradian*candela/meter**2

# katal is the SI unit of catalytic activity
katal = kat = Quantity("katal")
SI_quantity_dimension_map[katal] = amount_of_substance/time
SI_quantity_scale_factors[katal] = mol/second

# gray is the SI unit of absorbed dose
gray = Gy = Quantity("gray")
SI_quantity_dimension_map[gray] = energy/mass
SI_quantity_scale_factors[gray] = meter**2/second**2

# becquerel is the SI unit of radioactivity
becquerel = Bq = Quantity("becquerel")
SI_quantity_dimension_map[becquerel] = 1/time
SI_quantity_scale_factors[becquerel] = 1/second


# Common length units

km = kilometer = kilometers = Quantity("kilometer", abbrev="km")
SI_quantity_dimension_map[kilometer] = length
SI_quantity_scale_factors[kilometer] = kilo*meter

dm = decimeter = decimeters = Quantity("decimeter", abbrev="dm")
SI_quantity_dimension_map[decimeter] = length
SI_quantity_scale_factors[decimeter] = deci*meter

cm = centimeter = centimeters = Quantity("centimeter", abbrev="cm")
SI_quantity_dimension_map[centimeter] = length
SI_quantity_scale_factors[centimeter] = centi*meter

mm = millimeter = millimeters = Quantity("millimeter", abbrev="mm")
SI_quantity_dimension_map[millimeter] = length
SI_quantity_scale_factors[millimeter] = milli*meter

um = micrometer = micrometers = micron = microns = Quantity("micrometer", abbrev="um")
SI_quantity_dimension_map[micrometer] = length
SI_quantity_scale_factors[micrometer] = micro*meter

nm = nanometer = nanometers = Quantity("nanometer", abbrev="nn")
SI_quantity_dimension_map[nanometer] = length
SI_quantity_scale_factors[nanometer] = nano*meter

pm = picometer = picometers = Quantity("picometer", abbrev="pm")
SI_quantity_dimension_map[picometer] = length
SI_quantity_scale_factors[picometer] = pico*meter


ft = foot = feet = Quantity("foot", abbrev="ft")
SI_quantity_dimension_map[foot] = length
SI_quantity_scale_factors[foot] = Rational(3048, 10000)*meter

inch = inches = Quantity("inch")
SI_quantity_dimension_map[inch] = length
SI_quantity_scale_factors[inch] = foot/12

yd = yard = yards = Quantity("yard", abbrev="yd")
SI_quantity_dimension_map[yard] = length
SI_quantity_scale_factors[yard] = 3*feet

mi = mile = miles = Quantity("mile")
SI_quantity_dimension_map[mile] = length
SI_quantity_scale_factors[mile] = 5280*feet

nmi = nautical_mile = nautical_miles = Quantity("nautical_mile")
SI_quantity_dimension_map[nautical_mile] = length
SI_quantity_scale_factors[nautical_mile] = 6076*feet


# Common volume and area units

l = liter = liters = Quantity("liter")
SI_quantity_dimension_map[liter] = length**3
SI_quantity_scale_factors[liter] = meter**3 / 1000

dl = deciliter = deciliters = Quantity("deciliter")
SI_quantity_dimension_map[deciliter] = length**3
SI_quantity_scale_factors[deciliter] = liter / 10

cl = centiliter = centiliters = Quantity("centiliter")
SI_quantity_dimension_map[centiliter] = length**3
SI_quantity_scale_factors[centiliter] = liter / 100

ml = milliliter = milliliters = Quantity("milliliter")
SI_quantity_dimension_map[milliliter] = length**3
SI_quantity_scale_factors[milliliter] = liter / 1000


# Common time units

ms = millisecond = milliseconds = Quantity("millisecond", abbrev="ms")
SI_quantity_dimension_map[millisecond] = time
SI_quantity_scale_factors[millisecond] = milli*second

us = microsecond = microseconds = Quantity("microsecond", abbrev="us")
SI_quantity_dimension_map[microsecond] = time
SI_quantity_scale_factors[microsecond] = micro*second

ns = nanosecond = nanoseconds = Quantity("nanosecond", abbrev="ns")
SI_quantity_dimension_map[nanosecond] = time
SI_quantity_scale_factors[nanosecond] = nano*second

ps = picosecond = picoseconds = Quantity("picosecond", abbrev="ps")
SI_quantity_dimension_map[picosecond] = time
SI_quantity_scale_factors[picosecond] = pico*second


minute = minutes = Quantity("minute")
SI_quantity_dimension_map[minute] = time
SI_quantity_scale_factors[minute] = 60*second

h = hour = hours = Quantity("hour")
SI_quantity_dimension_map[hour] = time
SI_quantity_scale_factors[hour] = 60*minute

day = days = Quantity("day")
SI_quantity_dimension_map[day] = time
SI_quantity_scale_factors[day] = 24*hour


anomalistic_year = anomalistic_years = Quantity("anomalistic_year")
SI_quantity_dimension_map[anomalistic_year] = time
SI_quantity_scale_factors[anomalistic_year] = 365.259636*day

sidereal_year = sidereal_years = Quantity("sidereal_year")
SI_quantity_dimension_map[sidereal_year] = time
SI_quantity_scale_factors[sidereal_year] = 31558149.540

tropical_year = tropical_years = Quantity("tropical_year")
SI_quantity_dimension_map[tropical_year] = time
SI_quantity_scale_factors[tropical_year] = 365.24219*day

common_year = common_years = Quantity("common_year")
SI_quantity_dimension_map[common_year] = time
SI_quantity_scale_factors[common_year] = 365*day

julian_year = julian_years = Quantity("julian_year")
SI_quantity_dimension_map[julian_year] = time
SI_quantity_scale_factors[julian_year] = 365.25*day

draconic_year = draconic_years = Quantity("draconic_year")
SI_quantity_dimension_map[draconic_year] = time
SI_quantity_scale_factors[draconic_year] = 346.62*day

gaussian_year = gaussian_years = Quantity("gaussian_year")
SI_quantity_dimension_map[gaussian_year] = time
SI_quantity_scale_factors[gaussian_year] = 365.2568983*day

full_moon_cycle = full_moon_cycles = Quantity("full_moon_cycle")
SI_quantity_dimension_map[full_moon_cycle] = time
SI_quantity_scale_factors[full_moon_cycle] = 411.78443029*day


year = years = tropical_year

#### CONSTANTS ####

# Newton constant
G = gravitational_constant = Quantity("gravitational_constant", abbrev="G")
SI_quantity_dimension_map[gravitational_constant] = length**3*mass**-1*time**-2
SI_quantity_scale_factors[gravitational_constant] = 6.67408e-11*m**3/(kg*s**2)

# speed of light
c = speed_of_light = Quantity("speed_of_light", abbrev="c")
SI_quantity_dimension_map[speed_of_light] = velocity
SI_quantity_scale_factors[speed_of_light] = 299792458*meter/second

# Wave impedance of free space
Z0 = wave_impendence = Quantity("wave_impendence", abbrev='Z_0')
SI_quantity_dimension_map[wave_impendence] = impedance
SI_quantity_scale_factors[wave_impendence] = 119.9169832*pi

# Reduced Planck constant
hbar = Quantity("hbar", abbrev="hbar")
SI_quantity_dimension_map[hbar] = action
SI_quantity_scale_factors[hbar] = 1.05457266e-34*joule*second

# Planck constant
planck = Quantity("planck", abbrev="h")
SI_quantity_dimension_map[planck] = action
SI_quantity_scale_factors[planck] = 2*pi*hbar

# Electronvolt
eV = electronvolt = electronvolts = Quantity("electronvolt", abbrev="eV")
SI_quantity_dimension_map[electronvolt] = energy
SI_quantity_scale_factors[electronvolt] = 1.60219e-19*joule

# Avogadro number
avogadro_number = Quantity("avogadro_number")
SI_quantity_dimension_map[avogadro_number] = One
SI_quantity_scale_factors[avogadro_number] = 6.022140857e23

# Avogadro constant
avogadro = avogadro_constant = Quantity("avogadro_constant")
SI_quantity_dimension_map[avogadro_constant] = amount_of_substance**-1
SI_quantity_scale_factors[avogadro_constant] = avogadro_number / mol

# Boltzmann constant
boltzmann = boltzmann_constant = Quantity("boltzmann_constant")
SI_quantity_dimension_map[boltzmann_constant] = energy/temperature
SI_quantity_scale_factors[boltzmann_constant] = 1.38064852e-23*joule/kelvin

# Stefan-Boltzmann constant
stefan = stefan_boltzmann_constant = Quantity("stefan_boltzmann_constant")
SI_quantity_dimension_map[stefan_boltzmann_constant] = energy*time**-1*length**-2*temperature**-4
SI_quantity_scale_factors[stefan_boltzmann_constant] = 5.670367e-8*joule/(s*m**2*kelvin**4)

# Atomic mass
amu = amus = atomic_mass_unit = atomic_mass_constant = Quantity("atomic_mass_constant")
SI_quantity_dimension_map[atomic_mass_constant] = mass
SI_quantity_scale_factors[atomic_mass_constant] = 1.660539040e-24*gram

# Molar gas constant
R = molar_gas_constant = Quantity("molar_gas_constant", abbrev="R")
SI_quantity_dimension_map[molar_gas_constant] = energy/(temperature * amount_of_substance)
SI_quantity_scale_factors[molar_gas_constant] = 8.3144598*joule/kelvin/mol

# Faraday constant
faraday_constant = Quantity("faraday_constant")
SI_quantity_dimension_map[faraday_constant] = charge/amount_of_substance
SI_quantity_scale_factors[faraday_constant] = 96485.33289*C/mol

# Josephson constant
josephson_constant = Quantity("josephson_constant", abbrev="K_j")
SI_quantity_dimension_map[josephson_constant] = frequency/voltage
SI_quantity_scale_factors[josephson_constant] = 483597.8525e9*hertz/V

# Von Klitzing constant
von_klitzing_constant = Quantity("von_klitzing_constant", abbrev="R_k")
SI_quantity_dimension_map[von_klitzing_constant] = voltage/current
SI_quantity_scale_factors[von_klitzing_constant] = 25812.8074555*ohm

# Acceleration due to gravity (on the Earth surface)
gee = gees = acceleration_due_to_gravity = Quantity("acceleration_due_to_gravity", abbrev="g")
SI_quantity_dimension_map[acceleration_due_to_gravity] = acceleration
SI_quantity_scale_factors[acceleration_due_to_gravity] = 9.80665*meter/second**2

# magnetic constant:
u0 = magnetic_constant = Quantity("magnetic_constant")
SI_quantity_dimension_map[magnetic_constant] = force/current**2
SI_quantity_scale_factors[magnetic_constant] = 4*pi/10**7 * newton/ampere**2

# electric constat:
e0 = electric_constant = vacuum_permittivity = Quantity("vacuum_permittivity")
SI_quantity_dimension_map[vacuum_permittivity] = capacitance/length
SI_quantity_scale_factors[vacuum_permittivity] = 1/(u0 * c**2)

# vacuum impedance:
Z0 = vacuum_impedance = Quantity("vacuum_impedance")
SI_quantity_dimension_map[vacuum_impedance] = impedance
SI_quantity_scale_factors[vacuum_impedance] = u0 * c

# Coulomb's constant:
coulomb_constant = coulombs_constant = electric_force_constant = Quantity("coulomb_constant", abbrev="k_e")
SI_quantity_dimension_map[coulomb_constant] = force*length**2/charge**2
SI_quantity_scale_factors[coulomb_constant] = 1/(4*pi*vacuum_permittivity)


atmosphere = atmospheres = atm = Quantity("atmosphere", abbrev="atm")
SI_quantity_dimension_map[atmosphere] = pressure
SI_quantity_scale_factors[atmosphere] = 101325 * pascal


kPa = kilopascal = Quantity("kilopascal", abbrev="kPa")
SI_quantity_dimension_map[kilopascal] = pressure
SI_quantity_scale_factors[kilopascal] = kilo*Pa

bar = bars = Quantity("bar", abbrev="bar")
SI_quantity_dimension_map[bar] = pressure
SI_quantity_scale_factors[bar] = 100*kPa

pound = pounds = Quantity("pound")  # exact
SI_quantity_dimension_map[pound] = mass
SI_quantity_scale_factors[pound] = 0.45359237 * kg

psi = Quantity("psi")
SI_quantity_dimension_map[psi] = pressure
SI_quantity_scale_factors[psi] = pound * gee / inch ** 2

dHg0 = 13.5951  # approx value at 0 C
mmHg = torr = Quantity("mmHg")
SI_quantity_dimension_map[mmHg] = pressure
SI_quantity_scale_factors[mmHg] = dHg0 * acceleration_due_to_gravity * kilogram / meter**2

mmu = mmus = milli_mass_unit = Quantity("milli_mass_unit")
SI_quantity_dimension_map[milli_mass_unit] = mass
SI_quantity_scale_factors[milli_mass_unit] = atomic_mass_unit/1000

quart = quarts = Quantity("quart")
SI_quantity_dimension_map[quart] = length**3
SI_quantity_scale_factors[quart] = Rational(231, 4) * inch**3


# Other convenient units and magnitudes

ly = lightyear = lightyears = Quantity("lightyear", abbrev="ly")
SI_quantity_dimension_map[lightyear] = length
SI_quantity_scale_factors[lightyear] = speed_of_light*julian_year

au = astronomical_unit = astronomical_units = Quantity("astronomical_unit", abbrev="AU")
SI_quantity_dimension_map[astronomical_unit] = length
SI_quantity_scale_factors[astronomical_unit] = 149597870691*meter


# Fundamental Planck units:
planck_mass = Quantity("planck_mass", abbrev="m_P")
SI_quantity_dimension_map[planck_mass] = mass
SI_quantity_scale_factors[planck_mass] = sqrt(hbar*speed_of_light/G)

planck_time = Quantity("planck_time", abbrev="t_P")
SI_quantity_dimension_map[planck_time] = time
SI_quantity_scale_factors[planck_time] = sqrt(hbar*G/speed_of_light**5)

planck_temperature = Quantity("planck_temperature", abbrev="T_P")
SI_quantity_dimension_map[planck_temperature] = temperature
SI_quantity_scale_factors[planck_temperature] = sqrt(hbar*speed_of_light**5/G/boltzmann**2)

planck_length = Quantity("planck_length", abbrev="l_P")
SI_quantity_dimension_map[planck_length] = length
SI_quantity_scale_factors[planck_length] = sqrt(hbar*G/speed_of_light**3)

planck_charge = Quantity("planck_charge", abbrev="q_P")
SI_quantity_dimension_map[planck_charge] = charge
SI_quantity_scale_factors[planck_charge] = sqrt(4*pi*electric_constant*hbar*speed_of_light)


# Derived Planck units:
planck_area = Quantity("planck_area")
SI_quantity_dimension_map[planck_area] = length**2
SI_quantity_scale_factors[planck_area] = planck_length**2

planck_volume = Quantity("planck_volume")
SI_quantity_dimension_map[planck_volume] = length**3
SI_quantity_scale_factors[planck_volume] = planck_length**3

planck_momentum = Quantity("planck_momentum")
SI_quantity_dimension_map[planck_momentum] = mass*velocity
SI_quantity_scale_factors[planck_momentum] = planck_mass * speed_of_light

planck_energy = Quantity("planck_energy", abbrev="E_P")
SI_quantity_dimension_map[planck_energy] = energy
SI_quantity_scale_factors[planck_energy] = planck_mass * speed_of_light**2

planck_force = Quantity("planck_force", abbrev="F_P")
SI_quantity_dimension_map[planck_force] = force
SI_quantity_scale_factors[planck_force] = planck_energy / planck_length

planck_power = Quantity("planck_power", abbrev="P_P")
SI_quantity_dimension_map[planck_power] = power
SI_quantity_scale_factors[planck_power] = planck_energy / planck_time

planck_density = Quantity("planck_density", abbrev="rho_P")
SI_quantity_dimension_map[planck_density] = mass/length**3
SI_quantity_scale_factors[planck_density] = planck_mass / planck_length**3

planck_energy_density = Quantity("planck_energy_density", abbrev="rho^E_P")
SI_quantity_dimension_map[planck_energy_density] = energy / length**3
SI_quantity_scale_factors[planck_energy_density] = planck_energy / planck_length**3

planck_intensity = Quantity("planck_intensity", abbrev="I_P")
SI_quantity_dimension_map[planck_intensity] = mass * time**(-3)
SI_quantity_scale_factors[planck_intensity] = planck_energy_density * speed_of_light

planck_angular_frequency = Quantity("planck_angular_frequency", abbrev="omega_P")
SI_quantity_dimension_map[planck_angular_frequency] = 1 / time
SI_quantity_scale_factors[planck_angular_frequency] = 1 / planck_time

planck_pressure = Quantity("planck_pressure", abbrev="p_P")
SI_quantity_dimension_map[planck_pressure] = pressure
SI_quantity_scale_factors[planck_pressure] = planck_force / planck_length**2

planck_current = Quantity("planck_current", abbrev="I_P")
SI_quantity_dimension_map[planck_current] = current
SI_quantity_scale_factors[planck_current] = planck_charge / planck_time

planck_voltage = Quantity("planck_voltage", abbrev="V_P")
SI_quantity_dimension_map[planck_voltage] = voltage
SI_quantity_scale_factors[planck_voltage] = planck_energy / planck_charge

planck_impedance = Quantity("planck_impedance", abbrev="Z_P")
SI_quantity_dimension_map[planck_impedance] = impedance
SI_quantity_scale_factors[planck_impedance] = planck_voltage / planck_current

planck_acceleration = Quantity("planck_acceleration", abbrev="a_P")
SI_quantity_dimension_map[planck_acceleration] = acceleration
SI_quantity_scale_factors[planck_acceleration] = speed_of_light / planck_time


# Information theory units:
bit = bits = Quantity("bit")
SI_quantity_dimension_map[bit] = information
SI_quantity_scale_factors[bit] = One

byte = bytes = Quantity("byte")
SI_quantity_dimension_map[byte] = information
SI_quantity_scale_factors[byte] = 8*bit

kibibyte = kibibytes = Quantity("kibibyte")
SI_quantity_dimension_map[kibibyte] = information
SI_quantity_scale_factors[kibibyte] = kibi*byte

mebibyte = mebibytes = Quantity("mebibyte")
SI_quantity_dimension_map[mebibyte] = information
SI_quantity_scale_factors[mebibyte] = mebi*byte

gibibyte = gibibytes = Quantity("gibibyte")
SI_quantity_dimension_map[gibibyte] = information
SI_quantity_scale_factors[gibibyte] = gibi*byte

tebibyte = tebibytes = Quantity("tebibyte")
SI_quantity_dimension_map[tebibyte] = information
SI_quantity_scale_factors[tebibyte] = tebi*byte

pebibyte = pebibytes = Quantity("pebibyte")
SI_quantity_dimension_map[pebibyte] = information
SI_quantity_scale_factors[pebibyte] = pebi*byte

exbibyte = exbibytes = Quantity("exbibyte")
SI_quantity_dimension_map[exbibyte] = information
SI_quantity_scale_factors[exbibyte] = exbi*byte


from sympy.physics.units.quantities import process_scale_factor, process_dimension
SI_quantity_dimension_map = {q: process_dimension(dimension) for q, dimension
        in SI_quantity_dimension_map.items()}
for q, scale_factor in SI_quantity_scale_factors.items():
    SI_quantity_scale_factors[q] = process_scale_factor(scale_factor)
del process_scale_factor, process_dimension

# check that scale factors are the right SI dimensions:
for _scale_factor, _dimension in zip(SI_quantity_scale_factors.values(), SI_quantity_dimension_map.values()):
    dimex = Quantity.get_dimensional_expr(_scale_factor)
    if dimex != 1:
        if not dimsys_default.equivalent_dims(_dimension, Dimension(dimex)):
            raise ValueError("quantity value and dimension mismatch")
del _scale_factor, _dimension
