from sympy import Rational, pi, sqrt, S
from sympy.physics.units.prefixes import kilo, milli, micro, deci, centi, nano, pico
from sympy.physics.units.quantities import Quantity

One = S.One

#### UNITS ####

# Dimensionless:
percent = percents = Quantity("percent", latex_repr=r"\%")
percent.set_global_relative_scale_factor(Rational(1, 100), One)

permille = Quantity("permille")
permille.set_global_relative_scale_factor(Rational(1, 1000), One)


# Angular units (dimensionless)
rad = radian = radians = Quantity("radian", abbrev="rad")
deg = degree = degrees = Quantity("degree", abbrev="deg", latex_repr=r"^\circ")
sr = steradian = steradians = Quantity("steradian", abbrev="sr")
mil = angular_mil = angular_mils = Quantity("angular_mil", abbrev="mil")

# Base units:
m = meter = meters = Quantity("meter", abbrev="m")

# gram; used to define its prefixed units
g = gram = grams = Quantity("gram", abbrev="g")

# NOTE: the `kilogram` has scale factor 1000. In SI, kg is a base unit, but
# nonetheless we are trying to be compatible with the `kilo` prefix. In a
# similar manner, people using CGS or gaussian units could argue that the
# `centimeter` rather than `meter` is the fundamental unit for length, but the
# scale factor of `centimeter` will be kept as 1/100 to be compatible with the
# `centi` prefix.  The current state of the code assumes SI unit dimensions, in
# the future this module will be modified in order to be unit system-neutral
# (that is, support all kinds of unit systems).
kg = kilogram = kilograms = Quantity("kilogram", abbrev="kg")
kg.set_global_relative_scale_factor(kilo, gram)

s = second = seconds = Quantity("second", abbrev="s")
A = ampere = amperes = Quantity("ampere", abbrev='A')
K = kelvin = kelvins = Quantity("kelvin", abbrev='K')
mol = mole = moles = Quantity("mole", abbrev="mol")
cd = candela = candelas = Quantity("candela", abbrev="cd")

mg = milligram = milligrams = Quantity("milligram", abbrev="mg")
mg.set_global_relative_scale_factor(milli, gram)

ug = microgram = micrograms = Quantity("microgram", abbrev="ug", latex_repr=r"\mu\text{g}")
ug.set_global_relative_scale_factor(micro, gram)

# derived units
newton = newtons = N = Quantity("newton", abbrev="N")
joule = joules = J = Quantity("joule", abbrev="J")
watt = watts = W = Quantity("watt", abbrev="W")
pascal = pascals = Pa = pa = Quantity("pascal", abbrev="Pa")
hertz = hz = Hz = Quantity("hertz", abbrev="Hz")

# MKSA extension to MKS: derived units
coulomb = coulombs = C = Quantity("coulomb", abbrev='C')
volt = volts = v = V = Quantity("volt", abbrev='V')
ohm = ohms = Quantity("ohm", abbrev='ohm', latex_repr=r"\Omega")
siemens = S = mho = mhos = Quantity("siemens", abbrev='S')
farad = farads = F = Quantity("farad", abbrev='F')
henry = henrys = H = Quantity("henry", abbrev='H')
tesla = teslas = T = Quantity("tesla", abbrev='T')
weber = webers = Wb = wb = Quantity("weber", abbrev='Wb')

# Other derived units:
optical_power = dioptre = diopter = D = Quantity("dioptre")
lux = lx = Quantity("lux", abbrev="lx")

# katal is the SI unit of catalytic activity
katal = kat = Quantity("katal", abbrev="kat")

# gray is the SI unit of absorbed dose
gray = Gy = Quantity("gray")

# becquerel is the SI unit of radioactivity
becquerel = Bq = Quantity("becquerel", abbrev="Bq")


# Common length units

km = kilometer = kilometers = Quantity("kilometer", abbrev="km")
km.set_global_relative_scale_factor(kilo, meter)

dm = decimeter = decimeters = Quantity("decimeter", abbrev="dm")
dm.set_global_relative_scale_factor(deci, meter)

cm = centimeter = centimeters = Quantity("centimeter", abbrev="cm")
cm.set_global_relative_scale_factor(centi, meter)

mm = millimeter = millimeters = Quantity("millimeter", abbrev="mm")
mm.set_global_relative_scale_factor(milli, meter)

um = micrometer = micrometers = micron = microns = \
    Quantity("micrometer", abbrev="um", latex_repr=r'\mu\text{m}')
um.set_global_relative_scale_factor(micro, meter)

nm = nanometer = nanometers = Quantity("nanometer", abbrev="nm")
nm.set_global_relative_scale_factor(nano, meter)

pm = picometer = picometers = Quantity("picometer", abbrev="pm")
pm.set_global_relative_scale_factor(pico, meter)

ft = foot = feet = Quantity("foot", abbrev="ft")
ft.set_global_relative_scale_factor(Rational(3048, 10000), meter)

inch = inches = Quantity("inch")
inch.set_global_relative_scale_factor(Rational(1, 12), foot)

yd = yard = yards = Quantity("yard", abbrev="yd")
yd.set_global_relative_scale_factor(3, feet)

mi = mile = miles = Quantity("mile")
mi.set_global_relative_scale_factor(5280, feet)

nmi = nautical_mile = nautical_miles = Quantity("nautical_mile")
nmi.set_global_relative_scale_factor(6076, feet)


# Common volume and area units

l = liter = liters = Quantity("liter")

dl = deciliter = deciliters = Quantity("deciliter")
dl.set_global_relative_scale_factor(Rational(1, 10), liter)

cl = centiliter = centiliters = Quantity("centiliter")
cl.set_global_relative_scale_factor(Rational(1, 100), liter)

ml = milliliter = milliliters = Quantity("milliliter")
ml.set_global_relative_scale_factor(Rational(1, 1000), liter)


# Common time units

ms = millisecond = milliseconds = Quantity("millisecond", abbrev="ms")
millisecond.set_global_relative_scale_factor(milli, second)

us = microsecond = microseconds = Quantity("microsecond", abbrev="us", latex_repr=r'\mu\text{s}')
microsecond.set_global_relative_scale_factor(micro, second)

ns = nanosecond = nanoseconds = Quantity("nanosecond", abbrev="ns")
nanosecond.set_global_relative_scale_factor(nano, second)

ps = picosecond = picoseconds = Quantity("picosecond", abbrev="ps")
picosecond.set_global_relative_scale_factor(pico, second)

minute = minutes = Quantity("minute")
minute.set_global_relative_scale_factor(60, second)

h = hour = hours = Quantity("hour")
hour.set_global_relative_scale_factor(60, minute)

day = days = Quantity("day")
day.set_global_relative_scale_factor(24, hour)

anomalistic_year = anomalistic_years = Quantity("anomalistic_year")
anomalistic_year.set_global_relative_scale_factor(365.259636, day)

sidereal_year = sidereal_years = Quantity("sidereal_year")
sidereal_year.set_global_relative_scale_factor(31558149.540, seconds)

tropical_year = tropical_years = Quantity("tropical_year")
tropical_year.set_global_relative_scale_factor(365.24219, day)

common_year = common_years = Quantity("common_year")
common_year.set_global_relative_scale_factor(365, day)

julian_year = julian_years = Quantity("julian_year")
julian_year.set_global_relative_scale_factor((365 + One/4), day)

draconic_year = draconic_years = Quantity("draconic_year")
draconic_year.set_global_relative_scale_factor(346.62, day)

gaussian_year = gaussian_years = Quantity("gaussian_year")
gaussian_year.set_global_relative_scale_factor(365.2568983, day)

full_moon_cycle = full_moon_cycles = Quantity("full_moon_cycle")
full_moon_cycle.set_global_relative_scale_factor(411.78443029, day)

year = years = tropical_year


#### CONSTANTS ####

# Newton constant
# REF: NIST SP 959 (June 2019)
G = gravitational_constant = Quantity("gravitational_constant", abbrev="G")

# speed of light
c = speed_of_light = Quantity("speed_of_light", abbrev="c")

# elementary charge
# REF: NIST SP 959 (June 2019)
elementary_charge = Quantity("elementary_charge", abbrev="e")

# Planck constant
# REF: NIST SP 959 (June 2019)
planck = Quantity("planck", abbrev="h")

# Reduced Planck constant
# REF: NIST SP 959 (June 2019)
hbar = Quantity("hbar", abbrev="hbar")

# Electronvolt
# REF: NIST SP 959 (June 2019)
eV = electronvolt = electronvolts = Quantity("electronvolt", abbrev="eV")

# Avogadro number
# REF: NIST SP 959 (June 2019)
avogadro_number = Quantity("avogadro_number")

# Avogadro constant
avogadro = avogadro_constant = Quantity("avogadro_constant")

# Boltzmann constant
# REF: NIST SP 959 (June 2019)
boltzmann = boltzmann_constant = Quantity("boltzmann_constant")

# Stefan-Boltzmann constant
# REF: NIST SP 959 (June 2019)
stefan = stefan_boltzmann_constant = Quantity("stefan_boltzmann_constant")

# Atomic mass
# REF: NIST SP 959 (June 2019)
amu = amus = atomic_mass_unit = atomic_mass_constant = Quantity("atomic_mass_constant")

# Molar gas constant
# REF: NIST SP 959 (June 2019)
R = molar_gas_constant = Quantity("molar_gas_constant", abbrev="R")

# Faraday constant
faraday_constant = Quantity("faraday_constant")

# Josephson constant
josephson_constant = Quantity("josephson_constant", abbrev="K_j")

# Von Klitzing constant
von_klitzing_constant = Quantity("von_klitzing_constant", abbrev="R_k")

# Acceleration due to gravity (on the Earth surface)
gee = gees = acceleration_due_to_gravity = Quantity("acceleration_due_to_gravity", abbrev="g")

# magnetic constant:
u0 = magnetic_constant = vacuum_permeability = Quantity("magnetic_constant")

# electric constat:
e0 = electric_constant = vacuum_permittivity = Quantity("vacuum_permittivity")

# vacuum impedance:
Z0 = vacuum_impedance = Quantity("vacuum_impedance", abbrev='Z_0', latex_repr=r'Z_{0}')

# Coulomb's constant:
coulomb_constant = coulombs_constant = electric_force_constant = \
    Quantity("coulomb_constant", abbrev="k_e")


atmosphere = atmospheres = atm = Quantity("atmosphere", abbrev="atm")


kPa = kilopascal = Quantity("kilopascal", abbrev="kPa")
kilopascal.set_global_relative_scale_factor(kilo, Pa)

bar = bars = Quantity("bar", abbrev="bar")

pound = pounds = Quantity("pound")  # exact

psi = Quantity("psi")

dHg0 = 13.5951  # approx value at 0 C
mmHg = torr = Quantity("mmHg")

mmu = mmus = milli_mass_unit = Quantity("milli_mass_unit")

quart = quarts = Quantity("quart")


# Other convenient units and magnitudes

ly = lightyear = lightyears = Quantity("lightyear", abbrev="ly")

au = astronomical_unit = astronomical_units = Quantity("astronomical_unit", abbrev="AU")


# Fundamental Planck units:
planck_mass = Quantity("planck_mass", abbrev="m_P", latex_repr=r'm_\text{P}')

planck_time = Quantity("planck_time", abbrev="t_P", latex_repr=r't_\text{P}')

planck_temperature = Quantity("planck_temperature", abbrev="T_P",
                              latex_repr=r'T_\text{P}')

planck_length = Quantity("planck_length", abbrev="l_P", latex_repr=r'l_\text{P}')

planck_charge = Quantity("planck_charge", abbrev="q_P", latex_repr=r'q_\text{P}')


# Derived Planck units:
planck_area = Quantity("planck_area")

planck_volume = Quantity("planck_volume")

planck_momentum = Quantity("planck_momentum")

planck_energy = Quantity("planck_energy", abbrev="E_P", latex_repr=r'E_\text{P}')

planck_force = Quantity("planck_force", abbrev="F_P", latex_repr=r'F_\text{P}')

planck_power = Quantity("planck_power", abbrev="P_P", latex_repr=r'P_\text{P}')

planck_density = Quantity("planck_density", abbrev="rho_P", latex_repr=r'\rho_\text{P}')

planck_energy_density = Quantity("planck_energy_density", abbrev="rho^E_P")

planck_intensity = Quantity("planck_intensity", abbrev="I_P", latex_repr=r'I_\text{P}')

planck_angular_frequency = Quantity("planck_angular_frequency", abbrev="omega_P",
                                    latex_repr=r'\omega_\text{P}')

planck_pressure = Quantity("planck_pressure", abbrev="p_P", latex_repr=r'p_\text{P}')

planck_current = Quantity("planck_current", abbrev="I_P", latex_repr=r'I_\text{P}')

planck_voltage = Quantity("planck_voltage", abbrev="V_P", latex_repr=r'V_\text{P}')

planck_impedance = Quantity("planck_impedance", abbrev="Z_P", latex_repr=r'Z_\text{P}')

planck_acceleration = Quantity("planck_acceleration", abbrev="a_P",
                               latex_repr=r'a_\text{P}')


# Information theory units:
bit = bits = Quantity("bit")

byte = bytes = Quantity("byte")

kibibyte = kibibytes = Quantity("kibibyte")
mebibyte = mebibytes = Quantity("mebibyte")
gibibyte = gibibytes = Quantity("gibibyte")
tebibyte = tebibytes = Quantity("tebibyte")
pebibyte = pebibytes = Quantity("pebibyte")
exbibyte = exbibytes = Quantity("exbibyte")

# Older units for radioactivity
curie = Ci = Quantity("curie", abbrev="Ci")

rutherford = Rd = Quantity("rutherford", abbrev="Rd")
