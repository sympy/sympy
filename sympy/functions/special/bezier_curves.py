from sympy.core import S, sympify
from sympy.core.symbol import (Dummy, symbols)
from sympy.functions import Piecewise, piecewise_fold
from sympy.sets.sets import Interval
from sympy.functions.combinatorial.factorials import binomial
from sympy.geometry import Point
from sympy import expand, simplify
from sympy import interpolate


def bernstein_basis_polynomial(n,v,x):
    """
    Explanation
    ===========

    Bernstein polynomials are used for calculating weights of control points in bezier curves.

    Examples
    ========

    >>> from sympy import functions.special.bezier_curves
    >>> from bezier_curbes import bernstein_basis_polynomial
    >>> bernstein_basis_polynomial(1,0,0)
    0

    For the most basic 1d bernstein basis

    Parameters
    ==========

    n : integer
        degree of curve
    v: integer
        current control point
    x:
        parameter from [0,1]
    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Bernstein_polynomial
  
    """
    return binomial(n, v) * (x**v) * ((1 - x)**(n - v))


def bezier_curve_point(points, degree, t):
    """
    Explanation
    ===========

    This function returns a Point object for a parameter t of the bezier curve. 
    It sums up all bernstein basis weights for all points and then outputs
    a single co-ordinate on the x-y plane representing position of the 
    curve at that point.

    Examples
    ========

    >>> from sympy import functions.special.bezier_curves
    >>> from bezier_curbes import bernstein_basis_polynomial, bezier_curve_point
    >>> from sympy.geometry import Point
    >>> control_points = [Point(0,0), Point(3,3)]
    >>> bezier_control_point(control_points,1,0.5)
    Point2D(1.5,1.5)

    For a linear interpolation bezier curve

    >>> control_points = [Point(0,0), Point(2,10), Point(4,0)]
    >>> bezier_control_point(control_points,2,1)
    Point2D(4,0)

    For a quadratic interpolation bezier curve

    Parameters
    ==========

    control_points : list of Point objects
        represents control points
    degree: integer
        degree of curve
    t:
        parameter from [0,1]

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/B%C3%A9zier_curve
    """
    bezier_point = Point(0,0)

    for i in range(degree+1):

        bernstein_sum = bernstein_basis_polynomial(degree,i,t) 
        
        bezier_point = Point(bezier_point.x + bernstein_sum*points[i].x,bezier_point.y + bernstein_sum*points[i].y)
        
    return bezier_point

def bernstein_coeffecients(degree):
    """
    Explanation
    ===========

    This function returns an equation in terms of degree + 1
    bezier control points. This generic equation with bernstein
    coeffecient weights for every point can be used for differentiating 
    the curve by keeping each point constant.
    This is faster than the explicit calculation used for each control point
    in bezier_curve_point. Uses symbolic vars: t and P0 to Pn

    Examples
    ========

    >>> from sympy import functions.special.bezier_curves
    >>> from bezier_curbes import bernstein_basis_polynomial, bernstein_coeffecients
    >>>  bernstein_coeffecients(3, False)
    -P0*(t - 1)**3 + 3*P1*t*(t - 1)**2 - 3*P2*t**2*(t - 1) + P3*t**3

    Parameters
    ==========

    degree:
        degree of curve
   
    """
    bezier_control_points = symbols(f'P0:{degree+1}')

    t = symbols('t')

    bernstein_symbol = lambda i, degree, t: binomial(degree, i) * (t**i) * ((1 - t)**(degree - i))

    piecewise_bezier_curve = sum(bernstein_symbol(i,degree,t)*bezier_control_points[i]
                                 for i in range(degree + 1))
    
    symbolic_eqn =  simplify(piecewise_bezier_curve)
    
    return symbolic_eqn
    

def bezier_subdivision_half(control_points, degree):
    """
    Explanation
    ===========

    Each bezier curve is defined by a control polygon for 
    its control points. This control polygon
    can be brought closer to the curve by bisecting adjecent points and creating
    two new curves that are new piecewise functions.
    This implements subdivision for a quadratic curve

    Examples
    ========

    >>> from sympy import functions.special.bezier_curves
    >>> from bezier_curbes import bernstein_basis_polynomial, bezier_subdivision_half
    >>> control_points = [(0,0),(2,6),(4,0)]
    >>> bezier_subdivision_half(control_points, 2)
    Piecewise((3.0*x, (x >= 0) & (x < 2)), (nan, (x >= 2) & (x < 4)))

    Parameters
    ==========

    degree:
        degree of curve
    control_points
        control points of curve

    References
    ==========

    .. [1] https://www.clear.rice.edu/comp360/lectures/old/SubdivisionTextNew.pdf
    
        
    """
    midpoint_1 = ((control_points[0][0]+control_points[1][0])/2, (control_points[0][1]+control_points[1][1]/2))
    midpoint_2 = ((control_points[1][0]+control_points[2][0]/2),(control_points[1][1]+control_points[2][1])/2)

    sub_curve_1 = [control_points[0],midpoint_1,control_points[1]]
    sub_curve_2 = [control_points[1],midpoint_2,control_points[2]]

    x = symbols('x')

    interpolating_curve_1,interpolating_curve_2 = interpolate(sub_curve_1, x), interpolate(sub_curve_2, x)
    ivl1 = (x>=sub_curve_1[0][0]) & (x<sub_curve_1[2][0])
    ivl2 = (x>=sub_curve_2[0][0]) & (x<sub_curve_2[2][0])

    return Piecewise((interpolating_curve_1, ivl1), (interpolating_curve_2, ivl2))

