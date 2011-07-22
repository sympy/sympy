"""A module for solving all kinds of equations.

    Examples
    --------
    >>> from sympy.solvers import solve
    >>> from sympy.abc import x
    >>> solve(x**5+5*x**4+10*x**3+10*x**2+5*x+1,x)
    [-1]
"""
from .solvers import solve, solve_linear_system, solve_linear_system_LU, \
    solve_undetermined_coeffs, tsolve, nsolve, solve_linear, checksol

from .recurr import rsolve, rsolve_poly, rsolve_ratio, rsolve_hyper

from .ode import checkodesol, classify_ode, ode_order, dsolve, \
    homogeneous_order

from .polysys import solve_poly_system, solve_triangulated

from .pde import pde_separate, pde_separate_add, pde_separate_mul
