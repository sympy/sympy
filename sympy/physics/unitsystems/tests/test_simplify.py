# -*- coding: utf-8 -*-

from __future__ import division

from sympy import Add, Pow, Mul

from sympy.physics.unitsystems import dim_simplify
from sympy.physics.unitsystems.systems import mks_dim

L, T = mks_dim["length"], mks_dim["time"]


def test_dim_simplify_add():
    assert dim_simplify(Add(L, L)) == L


def test_dim_simplify_mul():
    assert dim_simplify(Mul(L, T)) == L * T


def test_dim_simplify_pow():
    assert dim_simplify(Pow(L, 2)) == L**2


def test_dim_simplify_rec():
    assert dim_simplify(Mul(Add(L, L), T)) == L * T
