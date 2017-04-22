# -*- coding: utf-8 -*-

from __future__ import division

from sympy import sstr, pi
from sympy.physics.units import G

from sympy import Add, Pow, Mul, sin, Tuple, sqrt, sympify
from sympy.physics.units import coulomb
from sympy.physics.units import hbar
from sympy.physics.units import joule
from sympy.physics.units import kelvin
from sympy.physics.units import mile, speed_of_light, meter, second, minute, hour, day
from sympy.physics.units import centimeter
from sympy.physics.units import inch
from sympy.physics.units import kilogram
from sympy.physics.units import kilometer
from sympy.physics.units import length
from sympy.physics.units import newton
from sympy.physics.units import planck
from sympy.physics.units import planck_length
from sympy.physics.units import planck_mass
from sympy.physics.units import planck_temperature
from sympy.physics.units import planck_time
from sympy.physics.units import radians, degree
from sympy.physics.units import steradian
from sympy.physics.units import time, gram
from sympy.physics.units.util import dim_simplify, convert_to


def NS(e, n=15, **options):
    return sstr(sympify(e).evalf(n, **options), full_prec=True)


L = length
T = time


def test_dim_simplify_add():
    assert dim_simplify(Add(L, L)) == L
    assert dim_simplify(L + L) == L


def test_dim_simplify_mul():
    assert dim_simplify(L*T) == L*T
    assert dim_simplify(L * T) == L*T


def test_dim_simplify_pow():
    assert dim_simplify(Pow(L, 2)) == L**2
    assert dim_simplify(L**2) == L**2


def test_dim_simplify_rec():
    assert dim_simplify(Mul(Add(L, L), T)) == L*T
    assert dim_simplify((L + L) * T) == L*T


def test_dim_simplify_dimless():
    # TODO: this should be somehow simplified on its own,
    # without the need of calling `dim_simplify`:
    assert dim_simplify(sin(L*L**-1)**2*L).get_dimensional_dependencies() == L.get_dimensional_dependencies()
    assert dim_simplify(sin(L * L**(-1))**2 * L).get_dimensional_dependencies() == L.get_dimensional_dependencies()


def test_convert_to_quantities():
    assert convert_to(3, meter) == 3

    assert convert_to(mile, kilometer) == 1.609344*kilometer
    assert convert_to(meter/second, speed_of_light) == speed_of_light/299792458
    assert convert_to(299792458*meter/second, speed_of_light) == speed_of_light
    assert convert_to(2*299792458*meter/second, speed_of_light) == 2*speed_of_light
    assert convert_to(speed_of_light, meter/second) == 299792458*meter/second
    assert convert_to(2*speed_of_light, meter/second) == 599584916*meter/second
    assert convert_to(day, second) == 86400*second
    assert convert_to(2*hour, minute) == 120*minute
    assert convert_to(mile, meter) == 1609.344*meter
    assert convert_to(mile/hour, kilometer/hour) == 25146*kilometer/(15625*hour)
    assert convert_to(3*newton, meter/second) == 3*newton
    assert convert_to(3*newton, kilogram*meter/second**2) == 3*meter*kilogram/second**2
    assert convert_to(kilometer + mile, meter) == 2609.344*meter
    assert convert_to(2*kilometer + 3*mile, meter) == 6828.032*meter
    assert convert_to(inch**2, meter**2) == 16129*meter**2/25000000
    assert convert_to(3*inch**2, meter) == 48387*meter**2/25000000
    assert convert_to(2*kilometer/hour + 3*mile/hour, meter/second) == 53344*meter/(28125*second)
    assert convert_to(2*kilometer/hour + 3*mile/hour, centimeter/second) == 213376*centimeter/(1125*second)
    assert convert_to(kilometer * (mile + kilometer), meter) == 2609344 * meter ** 2

    assert convert_to(steradian, coulomb) == steradian
    assert convert_to(radians, degree) == 180*degree/pi
    assert convert_to(radians, [meter, degree]) == 180*degree/pi
    assert convert_to(pi*radians, degree) == 180*degree
    assert convert_to(pi, degree) == 180*degree


def test_convert_to_tuples_of_quantities():
    assert convert_to(speed_of_light, [meter, second]) == 299792458 * meter / second
    assert convert_to(speed_of_light, (meter, second)) == 299792458 * meter / second
    assert convert_to(speed_of_light, Tuple(meter, second)) == 299792458 * meter / second
    assert convert_to(joule, [meter, kilogram, second]) == kilogram*meter**2/second**2
    assert convert_to(joule, [centimeter, gram, second]) == 10000000*centimeter**2*gram/second**2
    assert convert_to(299792458*meter/second, [speed_of_light]) == speed_of_light
    assert convert_to(speed_of_light / 2, [meter, second, kilogram]) == meter/second*299792458 / 2
    # This doesn't make physically sense, but let's keep it as a conversion test:
    assert convert_to(2 * speed_of_light, [meter, second, kilogram]) == 2 * 299792458 * meter / second
    assert convert_to(G, [G, speed_of_light, planck]) == 1.0*G

    assert NS(convert_to(meter, [G, speed_of_light, hbar]), n=7) == '6.187242e+34*gravitational_constant**0.5000000*hbar**0.5000000*speed_of_light**(-1.500000)'
    assert NS(convert_to(planck_mass, kilogram), n=7) == '2.176471e-8*kilogram'
    assert NS(convert_to(planck_length, meter), n=7) == '1.616229e-35*meter'
    assert NS(convert_to(planck_time, second), n=6) == '5.39116e-44*second'
    assert NS(convert_to(planck_temperature, kelvin), n=7) == '1.416809e+32*kelvin'
    assert NS(convert_to(convert_to(meter, [G, speed_of_light, planck]), meter), n=10) == '1.000000000*meter'
