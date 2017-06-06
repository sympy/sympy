"""
Fundamental physical constants.
"""

from sympy import pi
import sympy.physics.units as u


# Sourced from: NIST Standard Reference Database 121:

Z_0 = 376.730313461 * u.ohms # Vacuum impedance
epsilon_0 = 8.854187817 * 10 ** -12 * ((u.m ** -3)*(u.kg ** -1)*(u.s ** 4)*(u.A ** 2)) # Electric vacuum constant
mu_0 = 4 * pi * 10 ** -7 * ((u.m)*(u.kg)*(u.s ** -2)*(u.A ** -2)) # Magnetic vacuum constant
G = 6.6738480 * 10 ** -11 * ((u.m ** 3)*(u.kg ** -1)*(u.s ** -2)) # Gravitational constant
h = 6.62606957 * 10 ** -34 * ((u.m ** 2)*(u.kg)*(u.s ** -1)) # Planck constant
l_p = 1.616199 * 10 ** -35 * (u.m) # Planck length
m_p = 2.17651 * 10 ** -8 * (u.kg) # Planck mass
t_p = 5.39106 * 10 ** -44 * (u.s ** -1) # Planck time
c = 299792458 * ((u.m)*(u.s ** -1)) # Speed of light in vacuum

# Sourced from: WMAPS
H_0 = 2 * 10 ** -18 * (u.s ** -1) # Hubble parameter
