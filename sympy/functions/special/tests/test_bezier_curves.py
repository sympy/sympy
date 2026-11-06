from sympy.core.symbol import (Dummy, symbols)
from sympy.geometry import Point
from sympy.functions.special import bezier_curves
from sympy.functions.combinatorial.factorials import binomial
from sympy import Piecewise

def test_bernstein_basis_10():
    assert(bezier_curves.bernstein_basis_polynomial(1,0,1) == 0)
    assert(bezier_curves.bernstein_basis_polynomial(1,1,0)==0)

def test_bernstein_basis_high_degree():
    assert(bezier_curves.bernstein_basis_polynomial(1000,500,0.5)==0.0252250181783608)

def test_bezier_point_degree2():
    control_points = [Point(0,1),Point(2,10), Point(4,0)]

    assert(bezier_curves.bezier_curve_point(control_points,2,0)==Point(0,1))
    assert(bezier_curves.bezier_curve_point(control_points,2,0.5)==Point(2,21/4))
    assert(bezier_curves.bezier_curve_point(control_points,2,1)==Point(4,0))
    assert(bezier_curves.bezier_curve_point(control_points,2,0.25)==Point(1,69/16))
    assert(bezier_curves.bezier_curve_point(control_points,2,0.75)==Point(3, 61/16))

def test_bernstein_coeffecients():

    assert(str(bezier_curves.bernstein_coeffecients(2))== "P0*(t - 1)**2 - 2*P1*t*(t - 1) + P2*t**2")
    assert(str(bezier_curves.bernstein_coeffecients(3))== "-P0*(t - 1)**3 + 3*P1*t*(t - 1)**2 - 3*P2*t**2*(t - 1) + P3*t**3 P0*(t - 1)**2 - 2*P1*t*(t - 1) + P2*t**2")

def test_subdvision():
    
    assert(str(bezier_curves.bezier_subdivision_half([(0,0),(2,6),(4,0)],2))=="Piecewise((3.0*x, (x >= 0) & (x < 2)), (nan, (x >= 2) & (x < 4)))")
