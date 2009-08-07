"""A module for solving all kinds of equations.

In [3]: solve(x**5+5*x**4+10*x**3+10*x**2+5*x+1,x)
Out[3]: [-1]
"""
from solvers import solve, solve_linear_system, solve_linear_system_LU, \
    solve_undetermined_coeffs, tsolve, msolve, nsolve

from recurr import rsolve, rsolve_poly, rsolve_ratio, rsolve_hyper

from ode import classify_ode, dsolve, deriv_degree, homogeneous_order

from polysys import solve_poly_system

from pde import pde_separate, pde_separate_add, pde_separate_mul