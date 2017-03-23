# -*- coding: utf-8 -*-

from __future__ import division

from sympy import Add, Pow, Mul, sin
from sympy.physics.unitsystems import mile, speed_of_light, meter, second, minute, hour, day
from sympy.physics.unitsystems import centimeter
from sympy.physics.unitsystems import inch
from sympy.physics.unitsystems import kilogram
from sympy.physics.unitsystems import kilometer
from sympy.physics.unitsystems import length
from sympy.physics.unitsystems import newton
from sympy.physics.unitsystems import time
from sympy.physics.unitsystems.simplifiers import dim_simplify, convert_to
from sympy.physics.unitsystems.systems import _mks_dim

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
    assert convert_to(mile, kilometer) == 1.609344*kilometer
    assert convert_to(meter/second, speed_of_light) == speed_of_light/299792458
    assert convert_to(299792458*meter/second, speed_of_light) == speed_of_light
    assert convert_to(2*299792458*meter/second, speed_of_light) == 2*speed_of_light
    assert convert_to(speed_of_light, meter/second) == 299792458*meter/second
    assert convert_to(2*speed_of_light, meter/second) == 599584916*meter/second
    assert convert_to(day, second) == 86400*second
    assert convert_to(2*hour, minute) == 120*minute
    assert convert_to(mile, meter) == 1609.344*meter
    assert convert_to(mile/hour, kilometer/hour) == 1.609344*kilometer/hour
    assert convert_to(3*newton, meter/second) == 3*newton
    assert convert_to(3*newton, kilogram*meter/second**2) == 3*meter*kilogram/second**2
    assert convert_to(kilometer + mile, meter) == 2609.344*meter
    assert convert_to(2*kilometer + 3*mile, meter) == 6828.032*meter
    assert convert_to(inch**2, meter**2) == 0.00064516*meter**2
    assert convert_to(3*inch**2, meter) == 0.00193548*meter**2
    assert convert_to(2*kilometer/hour + 3*mile/hour, meter/second) == 1.89667555555556*meter/second
    assert convert_to(2*kilometer/hour + 3*mile/hour, centimeter/second) == 189.667555555556*centimeter/second
