"""A module for solving all kinds of equations.

    Examples
    ========

    >>> from sympy.solvers import solve
    >>> from sympy.abc import x
    >>> solve(x**5+5*x**4+10*x**3+10*x**2+5*x+1,x)
    [-1]
"""
from .decompogen import decompogen
from .deutils import ode_order
from .diophantine import diophantine
from .inequalities import reduce_abs_inequalities, reduce_abs_inequality, \
    reduce_inequalities, solve_poly_inequality, solve_rational_inequalities, \
    solve_univariate_inequality
from .ode import checkodesol, classify_ode, dsolve, homogeneous_order
from .pde import checkpdesol, classify_pde, pde_separate, pde_separate_add, \
    pde_separate_mul, pdsolve
from .polysys import solve_poly_system, solve_triangulated
from .recurr import rsolve, rsolve_hyper, rsolve_poly, rsolve_ratio
from .solvers import check_assumptions, checksol, det_quick, inv_quick, \
    nsolve, solve, solve_linear, solve_linear_system, solve_linear_system_LU, \
    solve_undetermined_coeffs
from .solveset import linear_eq_to_matrix, linsolve, nonlinsolve, solveset, \
    substitution
