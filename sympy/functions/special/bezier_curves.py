from sympy.core import S, sympify
from sympy.core.symbol import (Dummy, symbols)
from sympy.functions import Piecewise, piecewise_fold
from sympy.logic.boolalg import And
from sympy.sets.sets import Interval
from sympy import binomial
from sympy.geometry import Point
from sympy import expand, simplify

def bernstein_basis_polynomial(n,v,x):
    """
    Returns bernstein basis polynomial: https://en.wikipedia.org/wiki/Bernstein_polynomial
    """
    return binomial(n, v) * (x**v) * ((1 - x)**(n - v))


def bezier_curve_point(points, degree, t):
    """
    points: array of Point objects that are weighed by bernstein basis polynomials
    degree: nth degree bezier curve parametrized by t:[0,1]
    returns a cartesian point for the bezier curve
    """
    bezier_point = Point(0,0)
    for i in range(degree+1):
        bernstein_sum = bernstein_basis_polynomial(degree,i,t) 
        
        bezier_point = Point(bezier_point.x + bernstein_sum*points[i].x,bezier_point.y + bernstein_sum*points[i].y)
        
    return bezier_point

def bernstein_coeffecients(degree):
    """
    points: array of Point objects that are weighed by bernstein basis polynomials
    degree: nth degree bezier curve parametrized by t:[0,1]
    returns a linear combination of the control points for a parameter t
    """
    bezier_control_points = symbols(f'P0:{degree+1}')
    t = symbols('t')
    bernstein_symbol = lambda i, degree, t: binomial(degree, i) * (t**i) * ((1 - t)**(degree - i))
    piecewise_bezier_curve = sum(bernstein_symbol(i,degree,t)*bezier_control_points[i]
                                 for i in range(degree + 1))
    symbolic_eqn =  simplify(piecewise_bezier_curve)
    return symbolic_eqn

control_points = [Point(0,1),Point(2,10), Point(4,0)]
print(bezier_curve_point(control_points,2,0.75))