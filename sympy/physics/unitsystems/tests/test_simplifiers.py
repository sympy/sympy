# -*- coding: utf-8 -*-

from __future__ import division

from sympy import Add, Pow, Mul

from sympy.physics.unitsystems.simplifiers import dim_simplify, qsimplify
from sympy.physics.unitsystems.quantities import Quantity as Q
from sympy.physics.unitsystems.systems import mks, mks_dim

L, T = mks_dim["length"], mks_dim["time"]


def test_dim_simplify_add():
    assert dim_simplify(Add(L, L)) == L


def test_dim_simplify_mul():
    assert dim_simplify(Mul(L, T)) == L.mul(T)


def test_dim_simplify_pow():
    assert dim_simplify(Pow(L, 2)) == L.pow(2)


def test_dim_simplify_rec():
    assert dim_simplify(Mul(Add(L, L), T)) == L.mul(T)


m, s = mks["m"], mks["s"]

q1 = Q(10, m)
q2 = Q(5, m)


def test_qsimplify_add():
    assert qsimplify(Add(q1, q2)) == q1.add(q2)


def test_qsimplify_mul():
    q3 = Q(2, s)

    assert qsimplify(Mul(q1, q2)) == q1.mul(q2)
    assert qsimplify(Mul(q1, q3)) == q1.mul(q3)


def test_qsimplify_pow():
    assert qsimplify(Pow(q1, 2)) == q1.pow(2)


def test_qsimplify_rec():
    q3 = Q(2, s)

    assert qsimplify(Mul(Add(q1, q2), q3)) == q1.add(q2).mul(q3)
