# -*- coding: utf-8 -*-

from __future__ import division

from sympy import Add, Pow, Mul, sin
from sympy.physics.unitsystems import length
from sympy.physics.unitsystems import time
from sympy.physics.unitsystems.simplifiers import dim_simplify
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
