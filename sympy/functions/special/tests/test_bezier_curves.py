from sympy.core import S, sympify
from sympy.core.symbol import (Dummy, symbols)
from sympy.functions.special.bezier_curves import bernstein_basis_polynomial,bezier_curve_point
from sympy import binomial
from sympy.geometry import Point

def test_bernstein_basis_10():
    assert(bernstein_basis_polynomial(1,0,1) == 0)

def test_bernstein_basis_high_degree():
    assert(bernstein_basis_polynomial(1000,500,0.5)==0.0252250181783608)

def test_bezier_point_degree2():
    control_points = [Point(0,1),Point(2,10), Point(4,0)]

    assert(bezier_curve_point(control_points,2,0)==Point(0,1))
    assert(bezier_curve_point(control_points,2,0.5)==Point(2,21/4))
    assert(bezier_curve_point(control_points,2,1)==Point(4,0))
    assert(bezier_curve_point(control_points,2,0.25)==Point(1,69/16))
    assert(bezier_curve_point(control_points,2,0.75)==Point(3, 61/16))

