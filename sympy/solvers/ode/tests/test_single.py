#
# The main tests for the code in single.py are currently located in
# sympy/solvers/tests/test_ode.py
#

from sympy.core import Function, Symbol
from sympy.solvers.ode.single import (NthAlgebraic, ODEMatchError,
    SingleODEProblem, SingleODESolver)

from sympy.testing.pytest import raises


x = Symbol('x')
f = Function('f')


def test_SingleODESolver():
    # Test that not implemented methods give NotImplementedError
    # Subclasses should override these methods.
    problem = SingleODEProblem(f(x).diff(x), f(x), x)
    solver = SingleODESolver(problem)
    raises(NotImplementedError, lambda: solver.matches())
    raises(NotImplementedError, lambda: solver.get_general_solution())
    raises(NotImplementedError, lambda: solver._matches())
    raises(NotImplementedError, lambda: solver._get_general_solution())

    # This ODE can not be solved by the NthAlgebraic solver. Here we test that
    # it does not match and the asking for a general solution gives
    # ODEMatchError

    problem = SingleODEProblem(f(x).diff(x) + f(x), f(x), x)

    solver = NthAlgebraic(problem)
    raises(ODEMatchError, lambda: solver.get_general_solution())

    solver = NthAlgebraic(problem)
    assert solver.matches() is False
