from sympy import Eq, Function, symbols, exp
from sympy.abc import t
from sympy.solvers.ode.laplace_solver import laplace_solve

def test_first_order_laplace():
    y = Function('y')(t)
    ode = Eq(y.diff(t) + y, 1)
    sol = laplace_solve(ode, y, y0=0)
    assert sol == 1 - exp(-t)
