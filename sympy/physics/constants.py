"""
Physical Constants

Values of constants are recommended by Committee on Data for Science and
Technology (CODATA) as of 2014 (to be updated in 2018).
See more at http://physics.nist.gov/cuu/Constants/Table/allascii.txt
"""

from sympy import Rational, pi
import sympy.physics.units

m = sympy.physics.units.m
s = sympy.physics.units.s
ten = sympy.physics.units.ten
kg = sympy.physics.units.kg
N = sympy.physics.units.N
A = sympy.physics.units.A
K = sympy.physics.units.K
J = sympy.physics.units.J

c = speed_of_light = 299792458 * m/s
G = gravitational_constant = Rational('6.67408') * ten**-11 * m**3 / kg / s**2
u0 = magnetic_constant = 4*pi * ten**-7 * N/A**2
e0 = electric_constant = 1/(u0 * c**2)
Z0 = vacuum_impedance = u0 * c
a0 = bohr_radius = Rational('5.2917721067') * ten**-11 * m
b = wien_displacement_constant = Rational('2.8977729') * ten**-3 * m * K

planck = Rational('6.62607004') * ten**-34 * J*s
hbar = planck / (2*pi)

avogadro_number = Rational('6.022140857') * 10**23
boltzmann = Rational('1.38064852') * ten**-23 * J / K
atomic_mass_constant = Rational('1.660539040') * ten**-27 * kg
dHg0 = 13.5951  # approx value at 0 C

# Delete this so it doesn't pollute the namespace
del m, s, ten, kg, N, A, K, J, Rational, pi