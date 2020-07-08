from sympy import tan, pi
from sympy.abc import theta
from sympy.core import expand_trig, Eq
from sympy.solvers import solveset, nonlinsolve
from sympy.simplify.trigsimp import trigsimp
from sympy.vector import ParametricRegion


def ImplicitConicSection(equation, x, y):
    """
    Returns the parametric form of the conic section defined using implict equation.

    Example
    =======

    >>> from sympy.abc import x, y
    >>> from sympy.vector import ImplicitConicSection, vector_integrate

    >>> ImplicitConicSection(Eq(x**2 + y**2 -1), x, y)
    ParametricRegion((sin(2*theta), -cos(2*theta)), (theta, 0, pi))

    >>> c = ImplicitConicSection((x - 3)**2 + (y + 2)**2 - 16, x, y)
    >>> c
    ParametricRegion((2*(-sqrt(7)*tan(theta) + 3)*cos(theta)**2, 3*sin(2*theta) + sqrt(7)*cos(2*theta) - 2), (theta, 0, pi))
    >>> vector_integrate(1, c)
    8*pi

    >>> ellipse = ImplicitConicSection(x**2/4 + y**2/16 - 1, x, y)
    >>> ellipse
    ParametricRegion((8*tan(theta)/(tan(theta)**2 + 4), -4 + 8*tan(theta)**2/(tan(theta)**2 + 4)), (theta, 0, pi))
    >>> vector_integrate
    ### It gets stuck

    """

    if isinstance(equation, Eq):
        equation = equation.lhs - equation.rhs

    diff_x = diff(equation, x)
    diff_y = diff(equation, y)
    
    singular_points = nonlinsolve([equation, diff_x, diff_y], x, y)
    
    if len(singular_points) == 0:
        # No singular points
        raise ValueError()

    p_x, p_y = next(iter(singular_points))

    y_dash = tan(theta)*(x- p_x) + p_y
    eq2 = equation.subs(y, y_dash)
    x_par_list = solveset(eq2, x)

    for x_ in x_par_list:
        if x_ != p_x:
            x_par = x_
            break

    y_par = tan(theta)*(x_par - p_x) + p_y

    ans = {}
    ans[x] = trigsimp(expand_trig(x_par))
    ans[y] = trigsimp(expand_trig(y_par))

    p = ParametricRegion((ans[x], ans[y]), (theta, 0, pi))
    return p
