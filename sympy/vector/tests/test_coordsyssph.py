from sympy.vector.coordsysrect import CoordSysCartesian
from sympy.vector.coordsyssph import CoordSysSpherical
from sympy.vector.scalar import BaseScalar
from sympy import sin, cos, pi, ImmutableMatrix as Matrix, \
     symbols, simplify, zeros, expand
from sympy.vector.functions import express
from sympy.vector.point import Point
from sympy.vector.vector import Vector
from sympy.vector.orienters import (AxisOrienter, BodyOrienter,
                                    SpaceOrienter, QuaternionOrienter)

a, b, c, q = symbols('a b c q')
q1, q2, q3, q4 = symbols('q1 q2 q3 q4')


def test_func_args():
    A = CoordSysSpherical('A')
    assert A.r.func(*A.r.args) == A.r
    expr = 3*A.r + 4*A.r
    assert expr.func(*expr.args) == expr
    assert A.e_r.func(*A.e_r.args) == A.e_r
    v = A.r*A.e_r + A.theta*A.e_theta + A.phi*A.e_phi
    assert v.func(*v.args) == v
    assert A.origin.func(*A.origin.args) == A.origin


def test_coordsyscartesian_equivalence():
    A = CoordSysSpherical('A')
    A1 = CoordSysSpherical('A')
    assert A1 == A
    B = CoordSysSpherical('B')
    assert A != B
