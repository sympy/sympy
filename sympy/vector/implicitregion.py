from sympy import tan, pi
from sympy.abc import theta
from sympy.core import expand_trig, Eq
from sympy.solvers import solveset, nonlinsolve
from sympy.simplify.trigsimp import trigsimp
from sympy.vector import ParametricRegion

def ImplicitRegion(Basic):
    """
    Represents an implicit region in space. 

    Example
    =======

    >>> from sympy.abc import x, y
    >>> from sympy.vector import ImplicitRegion, parametric_region_list, vector_integrate

    >>> circle1 = ImplicitRegion(Eq(x**2 + y**2 -1), x, y)
    >>> parametric_region_list(circle)
    [ParametricRegion((sin(2*theta), -cos(2*theta)), (theta, 0, pi))]

    >>> circle2 = ImplicitRegion((x - 3)**2 + (y + 2)**2 - 16, x, y)
    >>> parametric_region_list(circle2)
    [ParametricRegion((2*(-sqrt(7)*tan(theta) + 3)*cos(theta)**2, 3*sin(2*theta) + sqrt(7)*cos(2*theta) - 2), (theta, 0, pi))]
    >>> vector_integrate(1, circle2)
    8*pi

    >>> ellipse = ImplicitRegion(x**2/4 + y**2/16 - 1, x, y)
    >>> parametric_region_list(ellipse)
    [ParametricRegion((8*tan(theta)/(tan(theta)**2 + 4), -4 + 8*tan(theta)**2/(tan(theta)**2 + 4)), (theta, 0, pi))]
    >>> vector_integrate(2, ellipse)
    ### It gets stuck

    """
    def __init__(cls, equation, x, y, z):
        if isinstance(equation, Eq):
            equation = equation.lhs - equation.rhs
        return super.__new__(cls, equation, x, y, z)
    
    @property
    def equation(self):
        return self.args[0]
